use anyhow::Context;
use compress_io::compress::CompressIo;
use crossbeam_channel::Sender;
use std::{io::BufRead, num::NonZeroU32, ops::Deref};

use crate::{
    cli::Config,
    kmcv,
    kmers::{KmerBuilder, KmerWork},
    regions::{Region, Regions},
};

#[derive(Default, Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd)]
#[repr(u8)]
pub enum Base {
    A = 0,
    C,
    T,
    G,
    N,
    #[default]
    Other,
}

impl Base {
    pub fn from_u8(c: u8) -> Self {
        match c {
            b'A' | b'a' => Self::A,
            b'C' | b'c' => Self::C,
            b'G' | b'g' => Self::G,
            b'T' | b't' => Self::T,
            b'N' | b'n' => Self::N,
            _ => Self::Other,
        }
    }

    pub fn is_gap(&self) -> bool {
        ((*self as usize) & 4) == 4
    }
}

struct RegionState<'a> {
    regions: &'a Regions,
    region_slice: Option<&'a [Region]>,
}

impl<'a> RegionState<'a> {
    fn new_contig(&mut self, ctg: &str) {
        debug!("Getting target regions for {ctg}");
        self.region_slice = self
            .regions
            .get(ctg)
            .and_then(|cr| {
                let v = cr.regions();
                debug!("{} regions found", v.len());
                if v.is_empty() {
                    None
                } else {
                    Some(v)
                }
            })
            .or_else(|| {
                debug!("No target regions found");
                None
            })
    }
    fn advance(&mut self) {
        if let Some(s) = self.region_slice.take() {
            if s.len() > 1 {
                self.region_slice = Some(&s[1..])
            }
        }
    }

    /// Returns true if position matches a target region
    fn check_pos(&mut self, pos: u32) -> Option<NonZeroU32> {
        while let Some(v) = self.region_slice {
            // As long as we use the API, v should always be non-empty
            let r = v[0];
            if pos > r.end() {
                self.advance()
            } else {
                return if pos >= r.start() {
                    Some(r.idx())
                } else {
                    None
                };
            }
        }
        None
    }
}

#[derive(Debug)]
pub struct Seq(Vec<Base>);

impl Deref for Seq {
    type Target = [Base];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
enum RdrState {
    Start,
    StartSeqId,
    InSeqId,
    StartSeq,
    StartSeqAfterNewLine,
    StartGap,
    InGap,
    InGapAfterNewLine,
    InSeq,
    InSeqAfterNewLine,
    InLongGap,
    InLongGapAfterNewLine,
    EndSeq,
    EndSeqAfterLongGap,
    StartSeqAfterInitialGap,
    NewContig,
}

struct Rdr<'a, R: BufRead> {
    r: R,
    state: RdrState,
    seq_id: String,
    max_read_length: u32,
    pos: u32,
    target_state: Option<RegionState<'a>>,
    k_work: KmerWork,
    kmer_build: KmerBuilder,
}

struct SeqWork<'a> {
    v: Vec<Base>,
    k_work: &'a mut KmerWork,
    k_build: &'a mut KmerBuilder,
}

impl<'a, R: BufRead> Rdr<'a, R> {
    fn new(r: R, max_read_length: u32, target_regions: Option<&'a Regions>) -> Self {
        let state = RdrState::Start;
        let seq_id = String::new();

        let target_state = target_regions.map(|r| RegionState {
            regions: r,
            region_slice: None,
        });

        let k_work = KmerWork::new();

        Self {
            r,
            state,
            seq_id,
            max_read_length,
            pos: 0,
            target_state,
            k_work,
            kmer_build: KmerBuilder::new(),
        }
    }

    fn get_seq(&mut self) -> anyhow::Result<Option<Seq>> {
        let v = Vec::new();
        let mut gap = 0;
        let mut ts = self.target_state.take();
        let mut seq_work = SeqWork {
            v,
            k_work: &mut self.k_work,
            k_build: &mut self.kmer_build,
        };

        loop {
            let buf = self.r.fill_buf()?;
            if buf.is_empty() {
                break;
            }
            let mut used = 0;
            let mut seq_ready = false;
            for (ix, c) in buf.iter().enumerate() {
                let idx = if let Some(t) = ts.as_mut() {
                    t.check_pos(self.pos)
                } else {
                    // If not targets are set, everything is on target!
                    None
                };
                trace!(
                    "pos = {}, idx = {:?}, state = {:?}",
                    self.pos,
                    idx,
                    self.state
                );
                let (new_state, inc_pos) = match self.state {
                    RdrState::Start => (proc_start(*c)?, false),
                    RdrState::StartSeqId => (proc_start_seq_id(*c, &mut self.seq_id)?, false),
                    RdrState::StartSeqAfterNewLine => proc_start_seq_after_new_line(*c)?,
                    RdrState::InSeqId => (proc_in_seq_id(*c, &mut self.seq_id)?, false),
                    RdrState::NewContig => {
                        debug!("Starting reading contig {}", self.seq_id);
                        if let Some(regs) = ts.as_mut() {
                            regs.new_contig(&self.seq_id)
                        }
                        seq_work.k_build.clear();
                        self.pos = 0;
                        proc_start_seq(*c)?
                    }
                    RdrState::StartSeq => proc_start_seq(*c)?,
                    RdrState::InSeq => {
                        gap = 0;
                        proc_in_seq(*c, Some(&mut seq_work), idx)?
                    }
                    RdrState::InSeqAfterNewLine => {
                        proc_after_new_line(*c, Some(&mut seq_work), proc_in_seq, idx)?
                    }
                    RdrState::InGapAfterNewLine => {
                        proc_after_new_line(*c, Some(&mut seq_work), proc_in_gap, idx)?
                    }
                    RdrState::InLongGapAfterNewLine => {
                        proc_after_new_line(*c, None, proc_in_long_gap, idx)?
                    }
                    RdrState::StartGap => {
                        gap = 1;
                        proc_in_gap(*c, Some(&mut seq_work), idx)?
                    }
                    RdrState::InGap => {
                        gap += 1;
                        if gap >= self.max_read_length {
                            assert!(seq_work.v.len() > gap as usize);
                            seq_work.v.truncate(seq_work.v.len() - gap as usize);
                            gap = 0;
                            proc_in_long_gap(*c, None, idx)?
                        } else {
                            proc_in_gap(*c, Some(&mut seq_work), idx)?
                        }
                    }
                    RdrState::InLongGap => proc_in_long_gap(*c, None, idx)?,
                    RdrState::EndSeq => {
                        used = ix;
                        seq_ready = true;
                        (RdrState::StartSeqId, false)
                    }
                    RdrState::EndSeqAfterLongGap => {
                        used = if ix > 0 { ix - 1 } else { ix };
                        if self.pos > 0 {
                            self.pos -= 1;
                        }
                        seq_ready = true;
                        (RdrState::StartSeq, false)
                    }
                    RdrState::StartSeqAfterInitialGap => {
                        used = if ix > 0 { ix - 1 } else { ix };
                        if self.pos > 0 {
                            self.pos -= 1;
                        }
                        seq_ready = true;
                        (RdrState::InSeq, false)
                    }
                };
                self.state = new_state;
                if inc_pos {
                    self.pos += 1
                }
                if seq_ready {
                    break;
                }
            }
            let used = if seq_ready {
                used
            } else if self.state == RdrState::EndSeqAfterLongGap {
                buf.len() - 1
            } else {
                buf.len()
            };
            self.r.consume(used);
            if seq_ready && !seq_work.v.is_empty() {
                break;
            }
        }

        self.target_state = ts;
        let SeqWork {
            mut v,
            k_work: _,
            k_build: _,
        } = seq_work;

        if gap > 0 {
            assert!(v.len() >= gap as usize);
            v.truncate(v.len() - gap as usize);
        }

        Ok(if v.is_empty() { None } else { Some(Seq(v)) })
    }
}

fn proc_in_gen(
    c: u8,
    sw: Option<&mut SeqWork>,
    s1: RdrState,
    s2: RdrState,
    s3: RdrState,
    target_idx: Option<NonZeroU32>,
) -> anyhow::Result<(RdrState, bool)> {
    if c == b'\n' {
        Ok((s1, false))
    } else if c.is_ascii_graphic() {
        let gc = Base::from_u8(c);
        if let Some(s) = sw {
            s.v.push(if target_idx.is_some() { gc } else { Base::N });
            s.k_build.add_base(gc, target_idx);
            trace!(
                "base: {:?}, kmers: {:?}, idx: {:?}",
                gc,
                s.k_build.kmers(),
                s.k_build.target_idx()
            );
            if let Some(k) = s.k_build.kmers() {
                let idx = s.k_build.target_idx();
                s.k_work.add_kmer(k[0], idx);
                s.k_work.add_kmer(k[1], idx);
            }
        } else {
            trace!("No SeqWork. Base: {:?}", gc);
        }
        Ok(if gc.is_gap() { (s2, true) } else { (s3, true) })
    } else {
        Err(anyhow!("Illegal character in sequence"))
    }
}

fn proc_in_gap(
    c: u8,
    sw: Option<&mut SeqWork>,
    target_idx: Option<NonZeroU32>,
) -> anyhow::Result<(RdrState, bool)> {
    proc_in_gen(
        c,
        sw,
        RdrState::InGapAfterNewLine,
        RdrState::InGap,
        RdrState::InSeq,
        target_idx,
    )
}

fn proc_in_long_gap(
    c: u8,
    sw: Option<&mut SeqWork>,
    target_idx: Option<NonZeroU32>,
) -> anyhow::Result<(RdrState, bool)> {
    proc_in_gen(
        c,
        sw,
        RdrState::InLongGapAfterNewLine,
        RdrState::InLongGap,
        RdrState::EndSeqAfterLongGap,
        target_idx,
    )
}

fn proc_after_new_line(
    c: u8,
    sw: Option<&mut SeqWork>,
    f: fn(
        c: u8,
        v: Option<&mut SeqWork>,
        target_idx: Option<NonZeroU32>,
    ) -> anyhow::Result<(RdrState, bool)>,
    target_idx: Option<NonZeroU32>,
) -> anyhow::Result<(RdrState, bool)> {
    if c == b'>' {
        Ok((RdrState::EndSeq, false))
    } else {
        f(c, sw, target_idx)
    }
}

fn proc_in_seq(
    c: u8,
    sw: Option<&mut SeqWork>,
    target_idx: Option<NonZeroU32>,
) -> anyhow::Result<(RdrState, bool)> {
    proc_in_gen(
        c,
        sw,
        RdrState::InSeqAfterNewLine,
        RdrState::StartGap,
        RdrState::InSeq,
        target_idx,
    )
}

fn proc_start_seq_after_new_line(c: u8) -> anyhow::Result<(RdrState, bool)> {
    if c == b'>' {
        Ok((RdrState::StartSeqId, false))
    } else {
        proc_start_seq(c)
    }
}

fn proc_start_seq(c: u8) -> anyhow::Result<(RdrState, bool)> {
    if c == b'\n' {
        Ok((RdrState::StartSeqAfterNewLine, false))
    } else if c.is_ascii_graphic() {
        let gc = Base::from_u8(c);
        Ok((
            if gc.is_gap() {
                RdrState::StartSeq
            } else {
                RdrState::StartSeqAfterInitialGap
            },
            true,
        ))
    } else {
        Err(anyhow!("Illegal character in sequence"))
    }
}

fn proc_in_seq_id(c: u8, s: &mut String) -> anyhow::Result<RdrState> {
    if c == b'\n' {
        if let Some(i) = s.find(char::is_whitespace) {
            s.truncate(i)
        }
        Ok(RdrState::NewContig)
    } else if c.is_ascii() && (c == b'\t' || !c.is_ascii_control()) {
        s.push(c as char);
        Ok(RdrState::InSeqId)
    } else {
        Err(anyhow!("Illegal character in sequence name"))
    }
}
fn proc_start_seq_id(c: u8, s: &mut String) -> anyhow::Result<RdrState> {
    s.clear();
    proc_in_seq_id(c, s)
}
fn proc_start(c: u8) -> anyhow::Result<RdrState> {
    if c == b'>' {
        Ok(RdrState::StartSeqId)
    } else {
        Err(anyhow!("Bad FASTA format: expecting '>'"))
    }
}

pub fn reader(cfg: &Config, snd: Sender<Seq>) -> anyhow::Result<()> {
    debug!(
        "Opening {} for input",
        cfg.input().and_then(|s| s.to_str()).unwrap_or("<stdin>")
    );
    let brdr = CompressIo::new()
        .opt_path(cfg.input())
        .bufreader()
        .with_context(|| "Could not open input file/stream")?;

    let max_rl = cfg.read_lengths().iter().max().unwrap();
    let mut rdr = Rdr::new(brdr, *max_rl, cfg.target_regions());

    info!("Starting to read input");
    while let Some(s) = rdr
        .get_seq()
        .with_context(|| "Error reading input sequence")?
    {
        snd.send(s)
            .with_context(|| "Error sending sequence for processing")?;
    }
    info!("Finished reading input");
    let k_work = rdr.k_work;
    info!("{k_work}");
    if let Some(reg) = cfg.target_regions() {
        info!("Outputting information on kmers");

        let output = format!("{}_kmers.km", cfg.prefix());

        kmcv::output_kmers(&output, reg, &k_work)
            .with_context(|| format!("Could not generate output kmer file {output}"))?;
    }
    Ok(())
}

mod test {
    #[allow(unused_imports)]
    use super::*;
    #[allow(unused_imports)]
    use std::io::BufReader;

    #[test]
    fn test1() {
        let s = ">seq1\nACTNNCCGT\nNACCAGTNNNNC\n>seq2\nNNN\n>seq3\nNNNNNNNNN\nNNNACTCNNN\n";
        let b = BufReader::new(s.as_bytes());
        let mut rdr = Rdr::new(b, 4, None);
        let exp_len = [16, 1, 4];
        for l in exp_len {
            let a = rdr.get_seq().unwrap().unwrap();
            println!("{:?}", a);
            assert_eq!(a.len(), l);
        }
        let a = rdr.get_seq().unwrap();
        assert!(a.is_none());
    }

    #[test]
    fn test2() {
        let s = ">seq1\nACTNNCCGT\nNACCAGTNNNNC\n>seq2\nNNN\n>seq3\nNNNNNNNNN\nNNNACTCNNN\n";
        let b = BufReader::with_capacity(16, s.as_bytes());
        let mut rdr = Rdr::new(b, 4, None);
        let exp_len = [16, 1, 4];
        for l in exp_len {
            let a = rdr.get_seq().unwrap().unwrap();
            println!("{:?}", a);
            assert_eq!(a.len(), l);
        }
        let a = rdr.get_seq().unwrap();
        assert!(a.is_none());
    }

    #[test]
    fn test3() {
        let s = ">seq1\nACTNNCCGT\nNACCAGTNNNNC\n>seq2\nNNN\n>seq3\nNNNNNNNNN\nNNNACTCNNN\n";
        let b = BufReader::with_capacity(30, s.as_bytes());
        let mut rdr = Rdr::new(b, 4, None);
        let exp_len = [16, 1, 4];
        for l in exp_len {
            let a = rdr.get_seq().unwrap().unwrap();
            println!("{:?}", a);
            assert_eq!(a.len(), l);
        }
        let a = rdr.get_seq().unwrap();
        assert!(a.is_none());
    }
}
