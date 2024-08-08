use std::{io::BufRead, ops::Deref};

use anyhow::Context;
use compress_io::compress::CompressIo;
use crossbeam_channel::Sender;

use crate::{
    cli::Config,
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
        self >= &Self::N
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
    fn check_pos(&mut self, pos: u32) -> bool {
        while let Some(v) = self.region_slice {
            // As long as we use the API, v should always be non-empty
            let r = v[0];
            if pos > r.end() {
                self.advance()
            } else {
                return pos >= r.start();
            }
        }
        false
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

#[derive(Copy, Clone, Eq, PartialEq)]
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
}

impl<'a, R: BufRead> Rdr<'a, R> {
    fn new(r: R, max_read_length: u32, target_regions: Option<&'a Regions>) -> Self {
        let state = RdrState::Start;
        let seq_id = String::new();
        let target_state = target_regions.map(|r| RegionState {
            regions: r,
            region_slice: None,
        });

        Self {
            r,
            state,
            seq_id,
            max_read_length,
            pos: 0,
            target_state,
        }
    }

    fn get_seq(&mut self) -> anyhow::Result<Option<Seq>> {
        let mut v = Vec::new();
        let mut gap = 0;
        let mut ts = self.target_state.take();
        loop {
            let buf = self.r.fill_buf()?;
            if buf.is_empty() {
                break;
            }
            let mut used = 0;
            let mut seq_ready = false;
            for (ix, c) in buf.iter().enumerate() {
                let ot = if let Some(t) = ts.as_mut() {
                    t.check_pos(self.pos)
                } else {
                    // If not targets are set, everything is on target!
                    true
                };
                let (new_state, inc_pos) = match self.state {
                    RdrState::Start => (proc_start(*c)?, false),
                    RdrState::StartSeqId => (proc_start_seq_id(*c, &mut self.seq_id)?, false),
                    RdrState::StartSeqAfterNewLine => proc_start_seq_after_new_line(*c, ot)?,
                    RdrState::InSeqId => (proc_in_seq_id(*c, &mut self.seq_id)?, false),
                    RdrState::NewContig => {
                        debug!("Starting reading contig {}", self.seq_id);
                        if let Some(regs) = ts.as_mut() {
                            regs.new_contig(&self.seq_id)
                        }
                        self.pos = 0;
                        proc_start_seq(*c, ot)?
                    }
                    RdrState::StartSeq => proc_start_seq(*c, ot)?,
                    RdrState::InSeq => {
                        gap = 0;
                        proc_in_seq(*c, Some(&mut v), ot)?
                    }
                    RdrState::InSeqAfterNewLine => {
                        proc_after_new_line(*c, Some(&mut v), proc_in_seq, ot)?
                    }
                    RdrState::InGapAfterNewLine => {
                        proc_after_new_line(*c, Some(&mut v), proc_in_gap, ot)?
                    }
                    RdrState::InLongGapAfterNewLine => {
                        proc_after_new_line(*c, None, proc_in_long_gap, ot)?
                    }
                    RdrState::StartGap => {
                        gap = 1;
                        proc_in_gap(*c, Some(&mut v), ot)?
                    }
                    RdrState::InGap => {
                        gap += 1;
                        if gap >= self.max_read_length {
                            assert!(v.len() > gap as usize);
                            v.truncate(v.len() - gap as usize);
                            gap = 0;
                            proc_in_long_gap(*c, None, ot)?
                        } else {
                            proc_in_gap(*c, Some(&mut v), ot)?
                        }
                    }
                    RdrState::InLongGap => proc_in_long_gap(*c, None, ot)?,
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
            if seq_ready && !v.is_empty() {
                break;
            }
        }
        if gap > 0 {
            assert!(v.len() >= gap as usize);
            v.truncate(v.len() - gap as usize);
        }
        self.target_state = ts;
        Ok(if v.is_empty() { None } else { Some(Seq(v)) })
    }
}

fn proc_in_gen(
    c: u8,
    v: Option<&mut Vec<Base>>,
    s1: RdrState,
    s2: RdrState,
    s3: RdrState,
    on_target: bool,
) -> anyhow::Result<(RdrState, bool)> {
    if c == b'\n' {
        Ok((s1, false))
    } else if c.is_ascii_graphic() {
        let gc = if on_target { Base::from_u8(c) } else { Base::N };
        if let Some(v) = v {
            v.push(gc)
        }
        Ok(if gc.is_gap() { (s2, true) } else { (s3, true) })
    } else {
        Err(anyhow!("Illegal character in sequence"))
    }
}

fn proc_in_gap(
    c: u8,
    v: Option<&mut Vec<Base>>,
    on_target: bool,
) -> anyhow::Result<(RdrState, bool)> {
    proc_in_gen(
        c,
        v,
        RdrState::InGapAfterNewLine,
        RdrState::InGap,
        RdrState::InSeq,
        on_target,
    )
}

fn proc_in_long_gap(
    c: u8,
    v: Option<&mut Vec<Base>>,
    on_target: bool,
) -> anyhow::Result<(RdrState, bool)> {
    proc_in_gen(
        c,
        v,
        RdrState::InLongGapAfterNewLine,
        RdrState::InLongGap,
        RdrState::EndSeqAfterLongGap,
        on_target,
    )
}

fn proc_after_new_line(
    c: u8,
    v: Option<&mut Vec<Base>>,
    f: fn(c: u8, v: Option<&mut Vec<Base>>, on_target: bool) -> anyhow::Result<(RdrState, bool)>,
    on_target: bool,
) -> anyhow::Result<(RdrState, bool)> {
    if c == b'>' {
        Ok((RdrState::EndSeq, false))
    } else {
        f(c, v, on_target)
    }
}

fn proc_in_seq(
    c: u8,
    v: Option<&mut Vec<Base>>,
    on_target: bool,
) -> anyhow::Result<(RdrState, bool)> {
    proc_in_gen(
        c,
        v,
        RdrState::InSeqAfterNewLine,
        RdrState::StartGap,
        RdrState::InSeq,
        on_target,
    )
}

fn proc_start_seq_after_new_line(c: u8, on_target: bool) -> anyhow::Result<(RdrState, bool)> {
    if c == b'>' {
        Ok((RdrState::StartSeqId, false))
    } else {
        proc_start_seq(c, on_target)
    }
}

fn proc_start_seq(c: u8, on_target: bool) -> anyhow::Result<(RdrState, bool)> {
    if c == b'\n' {
        Ok((RdrState::StartSeqAfterNewLine, false))
    } else if c.is_ascii_graphic() {
        let gc = if on_target { Base::from_u8(c) } else { Base::N };
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
