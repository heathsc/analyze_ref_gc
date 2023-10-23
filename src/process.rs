use std::{
    collections::{HashMap, VecDeque},
    ops::AddAssign,
};

use crossbeam_channel::{bounded, Receiver};
use crossbeam_utils::thread;
use serde::{Serialize, Serializer};

use crate::{
    cli::Config,
    reader::{self, GcBase, Seq},
};

#[derive(Copy, Clone, Eq, PartialOrd, PartialEq, Hash)]
pub struct GcHistKey(u32, u32);

impl GcHistKey {
    pub fn counts(&self) -> (u32, u32) { (self.0, self.1) }
}

impl Serialize for GcHistKey {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(&format!("{}:{}", self.0, self.1))
    }
}

#[derive(Serialize)]
struct GcHist {
    read_length: u32,
    counts: HashMap<GcHistKey, u64>,
}

impl GcHist {
    fn add(&mut self, other: &Self) {
        for (k, v) in other.counts.iter() {
            let e = self.counts.entry(*k).or_insert(0);
            *e += v
        }
    }

    fn new(read_length: u32) -> Self {
        Self {
            counts: HashMap::new(),
            read_length,
        }
    }
    pub fn hash(&self) -> &HashMap<GcHistKey, u64> {
        &self.counts
    }
}
#[derive(Serialize)]
pub struct GcRes {
    read_length_specific_counts: Box<[GcHist]>,
}

impl GcRes {
    pub fn new(rl: &[u32]) -> Self {
        let v: Vec<_> = rl.iter().map(|l| GcHist::new(*l)).collect();
        let inner = v.into_boxed_slice();
        Self {
            read_length_specific_counts: inner,
        }
    }

    fn add_count(&mut self, ix: usize, cts: (u32, u32)) {
        let e = self.read_length_specific_counts[ix]
            .counts
            .entry(GcHistKey(cts.0, cts.1))
            .or_insert(0);
        *e += 1
    }

    pub fn get_hash(&self, ix: usize) -> Option<&HashMap<GcHistKey, u64>> {
        self.read_length_specific_counts.get(ix).map(|g| g.hash())
    }

    pub fn len(&self) -> usize {
        self.read_length_specific_counts.len()
    }
}

impl AddAssign for GcRes {
    fn add_assign(&mut self, rhs: Self) {
        assert_eq!(
            self.read_length_specific_counts.len(),
            rhs.read_length_specific_counts.len()
        );
        for (p, q) in self
            .read_length_specific_counts
            .iter_mut()
            .zip(rhs.read_length_specific_counts.iter())
        {
            p.add(q)
        }
    }
}

#[derive(Copy, Clone)]
struct Counts {
    at: u32,
    gc: u32,
    threshold: u32,
}

impl Counts {
    fn new(threshold: u32) -> Self {
        assert!(threshold > 0);
        Self {
            at: 0,
            gc: 0,
            threshold,
        }
    }

    fn remove_gc_base(&mut self, base: &GcBase) {
        if let Some(x) = base.is_gc() {
            if x {
                assert!(self.gc > 0);
                self.gc -= 1
            } else {
                assert!(self.at > 0);
                self.at -= 1
            }
        }
    }

    fn add_gc_base(&mut self, base: &GcBase) {
        if let Some(x) = base.is_gc() {
            if x {
                self.gc += 1
            } else {
                self.at += 1
            }
        }
    }

    fn get_counts(&self) -> Option<(u32, u32)> {
        if self.at + self.gc >= self.threshold {
            Some((self.at, self.gc))
        } else {
            None
        }
    }
}

struct Work {
    buf: VecDeque<GcBase>,
    counts: Vec<Counts>,
}

impl Work {
    fn new(read_len: &[u32], threshold: f64) -> Self {
        let max_len = read_len
            .iter()
            .max()
            .map(|l| *l as usize)
            .expect("Empty read length vector");
        let mut buf = VecDeque::with_capacity(max_len);
        buf.resize_with(max_len, GcBase::default);
        let counts: Vec<_> = read_len
            .iter()
            .map(|l| Counts::new(((*l as f64) * threshold).ceil() as u32))
            .collect();

        Self { buf, counts }
    }

    fn clear(&mut self) {
        let l = self.buf.len();
        self.buf.clear();
        self.buf.resize_with(l, GcBase::default);
    }
}

fn process_seq(cfg: &Config, s: Seq, res: &mut GcRes, work: &mut Work) {
    let rl = cfg.read_lengths();
    let buf = &mut work.buf;
    let ct = &mut work.counts;
    let max_len = buf.len();
    let bnone = [GcBase::default()];
    let end = bnone.iter().cycle().take(max_len);

    for b in s.iter().chain(end) {
        // Decremeent counts from bases at start of reads
        for (l, c) in rl.iter().map(|l| *l as usize).zip(ct.iter_mut()) {
            assert!(l <= max_len);
            c.remove_gc_base(buf.get(max_len - l).unwrap());
        }
        // Remove base from start and add new base to end
        buf.pop_front();
        buf.push_back(*b);
        // Increment counts
        for (ix, c) in ct.iter_mut().enumerate() {
            c.add_gc_base(b);
            if let Some(cts) = c.get_counts() {
                res.add_count(ix, cts)
            }
        }
    }
}

fn process_thread(cfg: &Config, ix: usize, rx: Receiver<Seq>) -> anyhow::Result<GcRes> {
    debug!("Process task {ix} starting up");
    let mut res = GcRes::new(cfg.read_lengths());
    let mut work = Work::new(cfg.read_lengths(), cfg.threshold());
    while let Ok(mut s) = rx.recv() {
        trace!(
            "Process thread {ix} received new sequence of length {}",
            s.len()
        );
        process_seq(cfg, s, &mut res, &mut work);
        work.clear()
    }
    debug!("Process task {ix} shutting down");
    Ok(res)
}

pub fn process(cfg: &Config) -> anyhow::Result<GcRes> {
    let nt = cfg.threads();

    let mut error = false;
    let mut res = GcRes::new(cfg.read_lengths());

    thread::scope(|scope| {
        // Channel used to send sequences to process threads
        let (seq_send, seq_recv) = bounded(nt * 4);

        let mut process_tasks = Vec::with_capacity(nt);
        for ix in 0..nt {
            let rx = seq_recv.clone();
            let cfg = &cfg;
            process_tasks.push(scope.spawn(move |_| process_thread(cfg, ix, rx)));
        }
        drop(seq_recv);

        if let Err(e) = reader::reader(&cfg, seq_send) {
            error!("{:?}", e);
            error = true;
        }

        // Wait for analysis threads
        for jh in process_tasks.drain(..) {
            match jh.join().expect("Error joining analysis thread thread") {
                Err(e) => {
                    error!("{:?}", e);
                    error = true
                }
                Ok(r) => res += r,
            }
        }
    })
    .expect("Error in scope generation");

    if error {
        Err(anyhow!("Error occurred during processing"))
    } else {
        Ok(res)
    }
}
