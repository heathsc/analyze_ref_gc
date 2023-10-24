use std::{io::BufRead, ops::Deref};

use anyhow::Context;
use compress_io::compress::CompressIo;
use crossbeam_channel::Sender;

use crate::cli::Config;

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
}

struct Rdr<R: BufRead> {
    r: R,
    state: RdrState,
    seq_id: String,
    max_read_length: u32,
}

impl<R: BufRead> Rdr<R> {
    fn new(r: R, max_read_length: u32) -> Self {
        let state = RdrState::Start;
        let seq_id = String::new();
        Self {
            r,
            state,
            seq_id,
            max_read_length,
        }
    }

    fn get_seq(&mut self) -> anyhow::Result<Option<Seq>> {
        let mut v = Vec::new();
        let mut gap = 0;
        loop {
            let buf = self.r.fill_buf()?;
            if buf.is_empty() {
                break;
            }
            let mut used = 0;
            let mut seq_ready = false;
            for (ix, c) in buf.iter().enumerate() {
                self.state = match self.state {
                    RdrState::Start => proc_start(*c)?,
                    RdrState::StartSeqId => proc_start_seq_id(*c, &mut self.seq_id)?,
                    RdrState::StartSeqAfterNewLine => proc_start_seq_after_new_line(*c)?,
                    RdrState::InSeqId => proc_in_seq_id(*c, &mut self.seq_id)?,
                    RdrState::StartSeq => match proc_start_seq(*c)? {
                        RdrState::EndSeq => {
                            if v.is_empty() {
                                proc_in_seq(*c, Some(&mut v))?
                            } else {
                                used = ix;
                                seq_ready = true;
                                RdrState::InSeq
                            }
                        }
                        s => s,
                    },
                    RdrState::InSeq => {
                        gap = 0;
                        proc_in_seq(*c, Some(&mut v))?
                    }
                    RdrState::InSeqAfterNewLine => {
                        proc_after_new_line(*c, Some(&mut v), proc_in_seq)?
                    }
                    RdrState::InGapAfterNewLine => {
                        proc_after_new_line(*c, Some(&mut v), proc_in_gap)?
                    }
                    RdrState::InLongGapAfterNewLine => {
                        proc_after_new_line(*c, None, proc_in_long_gap)?
                    }
                    RdrState::StartGap => {
                        gap = 1;
                        proc_in_gap(*c, Some(&mut v))?
                    }
                    RdrState::InGap => {
                        gap += 1;
                        if gap >= self.max_read_length {
                            assert!(v.len() > gap as usize);
                            v.truncate(v.len() - gap as usize);
                            gap = 0;
                            proc_in_long_gap(*c, None)?
                        } else {
                            proc_in_gap(*c, Some(&mut v))?
                        }
                    }
                    RdrState::InLongGap => proc_in_long_gap(*c, None)?,
                    RdrState::EndSeq => {
                        used = ix;
                        seq_ready = true;
                        RdrState::StartSeqId
                    }
                    RdrState::EndSeqAfterLongGap => {
                        used = if ix > 0 { ix - 1 } else { ix };
                        seq_ready = true;
                        RdrState::StartSeq
                    }
                };
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
        Ok(if v.is_empty() { None } else { Some(Seq(v)) })
    }
}

fn proc_in_gen(
    c: u8,
    v: Option<&mut Vec<Base>>,
    s1: RdrState,
    s2: RdrState,
    s3: RdrState,
) -> anyhow::Result<RdrState> {
    if c == b'\n' {
        Ok(s1)
    } else if c.is_ascii_graphic() {
        let gc = Base::from_u8(c);
        if let Some(v) = v {
            v.push(gc)
        }
        Ok(if gc.is_gap() { s2 } else { s3 })
    } else {
        Err(anyhow!("Illegal character in sequence"))
    }
}

fn proc_in_gap(c: u8, v: Option<&mut Vec<Base>>) -> anyhow::Result<RdrState> {
    proc_in_gen(
        c,
        v,
        RdrState::InGapAfterNewLine,
        RdrState::InGap,
        RdrState::InSeq,
    )
}

fn proc_in_long_gap(c: u8, v: Option<&mut Vec<Base>>) -> anyhow::Result<RdrState> {
    proc_in_gen(
        c,
        v,
        RdrState::InLongGapAfterNewLine,
        RdrState::InLongGap,
        RdrState::EndSeqAfterLongGap,
    )
}

fn proc_after_new_line(
    c: u8,
    v: Option<&mut Vec<Base>>,
    f: fn(c: u8, v: Option<&mut Vec<Base>>) -> anyhow::Result<RdrState>,
) -> anyhow::Result<RdrState> {
    if c == b'>' {
        Ok(RdrState::EndSeq)
    } else {
        f(c, v)
    }
}

fn proc_in_seq(c: u8, v: Option<&mut Vec<Base>>) -> anyhow::Result<RdrState> {
    proc_in_gen(
        c,
        v,
        RdrState::InSeqAfterNewLine,
        RdrState::StartGap,
        RdrState::InSeq,
    )
}

fn proc_start_seq_after_new_line(c: u8) -> anyhow::Result<RdrState> {
    if c == b'>' {
        Ok(RdrState::StartSeqId)
    } else {
        proc_start_seq(c)
    }
}

fn proc_start_seq(c: u8) -> anyhow::Result<RdrState> {
    if c == b'\n' {
        Ok(RdrState::StartSeqAfterNewLine)
    } else if c.is_ascii_graphic() {
        let gc = Base::from_u8(c);
        Ok(if gc.is_gap() {
            RdrState::StartSeq
        } else {
            RdrState::EndSeq
        })
    } else {
        Err(anyhow!("Illegal character in sequence"))
    }
}

fn proc_in_seq_id(c: u8, s: &mut String) -> anyhow::Result<RdrState> {
    if c == b'\n' {
        debug!("Starting reading sequence {}", s);
        Ok(RdrState::StartSeq)
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
    let mut rdr = Rdr::new(brdr, *max_rl);

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
        let mut rdr = Rdr::new(b, 4);
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
        let mut rdr = Rdr::new(b, 4);
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
        let mut rdr = Rdr::new(b, 4);
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
