use std::{
    collections::{BTreeMap, HashMap, VecDeque},
    ops::AddAssign,
};

use crossbeam_channel::{bounded, Receiver};
use crossbeam_utils::thread;
use serde::{Serialize, Serializer};

use crate::{
    cli::Config,
    reader::{self, Base, Seq},
};

#[derive(Copy, Clone, Eq, PartialOrd, PartialEq, Hash)]
pub struct GcHistKey(u32, u32);

impl GcHistKey {
    pub fn counts(&self) -> (u32, u32) {
        (self.0, self.1)
    }
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
pub struct GcHist {
    counts: HashMap<GcHistKey, u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    bisulfite_counts: Option<HashMap<GcHistKey, u64>>,
}

impl GcHist {
    fn add(&mut self, other: &Self) {
        for (k, v) in other.counts.iter() {
            let e = self.counts.entry(*k).or_insert(0);
            *e += v
        }
        if let Some(ct) = self.bisulfite_counts.as_mut() {
            let ct1 = other.bisulfite_counts.as_ref().unwrap();
            for (k, v) in ct1.iter() {
                let e = ct.entry(*k).or_insert(0);
                *e += v
            }
        }
    }

    fn new(bisulfite: bool) -> Self {
        let bisulfite_counts = if bisulfite {
            Some(HashMap::new())
        } else {
            None
        };
        Self {
            counts: HashMap::new(),
            bisulfite_counts,
        }
    }
    pub fn hash(&self) -> &HashMap<GcHistKey, u64> {
        &self.counts
    }

    pub fn bisulfite_hash(&self) -> Option<&HashMap<GcHistKey, u64>> {
        self.bisulfite_counts.as_ref()
    }
}
#[derive(Serialize)]
pub struct GcRes {
    read_length_specific_counts: BTreeMap<u32, GcHist>,
}

impl GcRes {
    pub fn new(rl: &[u32], bisulfite: bool) -> Self {
        let inner: BTreeMap<_, _> = rl.iter().map(|l| (*l, GcHist::new(bisulfite))).collect();
        Self {
            read_length_specific_counts: inner,
        }
    }

    fn add_count(&mut self, ix: u32, cts: (u32, u32)) {
        let e = self
            .read_length_specific_counts
            .get_mut(&ix)
            .unwrap()
            .counts
            .entry(GcHistKey(cts.0, cts.1))
            .or_insert(0);
        *e += 1
    }

    fn add_bs_count(&mut self, ix: u32, cts: (u32, u32)) {
        if let Some(c) = self
            .read_length_specific_counts
            .get_mut(&ix)
            .unwrap()
            .bisulfite_counts
            .as_mut()
        {
            let e = c.entry(GcHistKey(cts.0, cts.1)).or_insert(0);
            *e += 1
        }
    }

    pub fn get_gc_hist(&self, ix: u32) -> Option<&GcHist> {
        self.read_length_specific_counts.get(&ix)
    }
}

impl AddAssign for GcRes {
    fn add_assign(&mut self, rhs: Self) {
        assert_eq!(
            self.read_length_specific_counts.len(),
            rhs.read_length_specific_counts.len()
        );
        for ((p_key, p_val), (q_key, q_val)) in self
            .read_length_specific_counts
            .iter_mut()
            .zip(rhs.read_length_specific_counts.iter())
        {
            assert_eq!(*p_key, *q_key);
            p_val.add(q_val)
        }
    }
}

#[derive(Copy, Clone)]
struct Counts {
    counts: [u32; 4],
    threshold: u32,
}

impl Counts {
    fn new(threshold: u32) -> Self {
        assert!(threshold > 0);
        Self {
            counts: [0; 4],
            threshold,
        }
    }

    fn remove_base(&mut self, base: &Base) {
        if !base.is_gap() {
            let i = *base as usize;
            assert!(self.counts[i] > 0);
            self.counts[i] -= 1
        }
    }

    fn add_base(&mut self, base: &Base) {
        if !base.is_gap() {
            self.counts[*base as usize] += 1
        }
    }

    fn get_counts(&self) -> Option<(u32, u32)> {
        if self.counts.iter().sum::<u32>() >= self.threshold {
            Some((
                self.counts[Base::A as usize] + self.counts[Base::T as usize],
                self.counts[Base::C as usize] + self.counts[Base::G as usize],
            ))
        } else {
            None
        }
    }

    fn get_bs_counts(&self) -> Option<((u32, u32), (u32, u32))> {
        if self.counts.iter().sum::<u32>() >= self.threshold {
            Some((
                (self.counts[Base::T as usize], self.counts[Base::C as usize]),
                (self.counts[Base::A as usize], self.counts[Base::G as usize]),
            ))
        } else {
            None
        }
    }
}

struct Work {
    buf: VecDeque<Base>,
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
        buf.resize_with(max_len, Base::default);
        let counts: Vec<_> = read_len
            .iter()
            .map(|l| Counts::new(((*l as f64) * threshold).ceil() as u32))
            .collect();

        Self { buf, counts }
    }

    fn clear(&mut self) {
        let l = self.buf.len();
        self.buf.clear();
        self.buf.resize_with(l, Base::default);
        for c in self.counts.iter_mut() {
            c.counts = [0, 0, 0, 0];
        }
    }
}

fn process_seq(cfg: &Config, s: Seq, res: &mut GcRes, work: &mut Work) {
    let rl = cfg.read_lengths();
    work.clear();
    let buf = &mut work.buf;
    let ct = &mut work.counts;
    let max_len = buf.len();
    let bnone = [Base::default()];
    let end = bnone.iter().cycle().take(max_len);

    for b in s.iter().chain(end) {
        // Decrement counts from bases at start of reads
        for (l, c) in rl.iter().map(|l| *l as usize).zip(ct.iter_mut()) {
            assert!(l <= max_len);
            c.remove_base(buf.get(max_len - l).unwrap());
        }
        // Remove base from start and add new base to end
        buf.pop_front();
        buf.push_back(*b);
        // Increment counts
        for (ix, c) in ct.iter_mut().enumerate() {
            c.add_base(b);
            if cfg.bisulfite() {
                if let Some((cts1, cts2)) = c.get_bs_counts() {
                    let cts = (cts1.0 + cts2.0, cts1.1 + cts2.1);
                    res.add_count(rl[ix], cts);
                    res.add_bs_count(rl[ix], cts1);
                    res.add_bs_count(rl[ix], cts2);
                }
            } else if let Some(cts) = c.get_counts() {
                res.add_count(rl[ix], cts)
            }
        }
    }
}

fn process_thread(cfg: &Config, ix: usize, rx: Receiver<Seq>) -> anyhow::Result<GcRes> {
    debug!("Process task {ix} starting up");
    let mut res = GcRes::new(cfg.read_lengths(), cfg.bisulfite());
    let mut work = Work::new(cfg.read_lengths(), cfg.threshold());
    while let Ok(s) = rx.recv() {
        trace!(
            "Process thread {ix} received new sequence of length {}",
            s.len()
        );
        process_seq(cfg, s, &mut res, &mut work);
    }
    debug!("Process task {ix} shutting down");
    Ok(res)
}

pub fn process(cfg: &Config) -> anyhow::Result<GcRes> {
    let nt = cfg.threads();

    let mut error = false;
    let mut res = GcRes::new(cfg.read_lengths(), cfg.bisulfite());

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
            match jh.join().expect("Error joining analysis thread") {
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
