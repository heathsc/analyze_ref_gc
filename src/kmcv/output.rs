use std::{io::Write, path::Path};

use anyhow::Context;
use compress_io::{
    compress::CompressIo,
    compress_type::{CompressThreads, CompressType},
};
use rand::random;

use crate::{
    kmers::{KmerVec, KmerWork, KMER_LENGTH, MAX_HITS},
    regions::Regions,
};

const MAJOR_VERSION: u8 = 2;
const MINOR_VERSION: u8 = 0;

#[inline]
fn u32_to_buf(b: &mut [u8], x: u32) {
    b.copy_from_slice(&x.to_le_bytes())
}

#[inline]
fn u64_to_buf(b: &mut [u8], x: u64) {
    b.copy_from_slice(&x.to_le_bytes())
}

struct KmcvHeader {
    buf: [u8; 52],
}

impl KmcvHeader {
    fn new(reg: &Regions, k_work: &KmerWork, rnd_id: u32) -> Self {
        let n_contigs = reg.n_contigs() as u32;
        let n_targets = reg.n_regions() as u32;
        let mapped = k_work.mapped_kmers();
        let redundant = k_work.highly_redundant_kmers();
        let on_target = k_work.on_target_kmers();
        let total_hits = k_work.total_hits();

        let mut buf = [0; 52];

        buf[0..4].copy_from_slice(&[b'K', b'M', b'C', b'V']);
        buf[4] = MAJOR_VERSION;
        buf[5] = MINOR_VERSION;
        buf[6] = KMER_LENGTH as u8;
        buf[7] = MAX_HITS as u8;
        u32_to_buf(&mut buf[8..12], rnd_id);
        u32_to_buf(&mut buf[12..16], n_contigs);
        u32_to_buf(&mut buf[16..20], n_targets);
        u64_to_buf(&mut buf[20..28], mapped);
        u64_to_buf(&mut buf[28..36], on_target);
        u64_to_buf(&mut buf[36..44], redundant);
        u64_to_buf(&mut buf[44..], total_hits);

        Self { buf }
    }

    fn write<W: Write>(&self, w: &mut W) -> anyhow::Result<()> {
        w.write_all(&self.buf)
            .with_context(|| "Error writing header to kmer file")
    }
}

fn write_target_blocks<W: Write>(w: &mut W, reg: &Regions) -> anyhow::Result<()> {
    for (ctg_ix, (_, ctg_regs)) in reg.iter().enumerate() {
        let ix = ctg_ix as u32;
        for r in ctg_regs.regions() {
            w.write_all(&ix.to_le_bytes())
                .with_context(|| "Error writing contig id for target")?;
            w.write_all(&r.start().to_le_bytes())
                .with_context(|| "Error writing target start")?;
            w.write_all(&r.end().to_le_bytes())
                .with_context(|| "Error writing target end")?;
        }
    }
    Ok(())
}

fn write_contig_blocks<W: Write>(w: &mut W, reg: &Regions) -> anyhow::Result<()> {
    for (ctg, _) in reg.iter() {
        let l = ctg.len();
        if l > u16::MAX as usize {
            return Err(anyhow!("Contig name is too long (size is {l}, max is {}", u16::MAX))
        }
        w.write_all(&(l as u16).to_le_bytes())
            .with_context(|| "Error writing contig name length")?;
        w.write_all(ctg.as_bytes())
            .with_context(|| "Error writing contig name")?;
    }
    Ok(())
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
enum KmerType {
    Unmapped,
    UniqueOnTarget,
    UniqueOffTarget,
    LowMultiMap(u8),
    HighMultiMap,
}

impl KmerType {
    fn from_kmer_vec(v: &KmerVec) -> Self {
        if v[0] == 0 {
            Self::Unmapped
        } else if (v[0] & 0x80000000) != 0 {
            Self::HighMultiMap
        } else if v[1] == 0 {
            if v[0] == 1 {
                Self::UniqueOffTarget
            } else {
                Self::UniqueOnTarget
            }
        } else {
            let mut n_hits = None;
            for (i, x) in v[2..].iter().enumerate() {
                if *x == 0 {
                    n_hits = Some(i + 2);
                    break;
                }
            }
            let n_hits = n_hits.unwrap_or(v.len()) as u8;
            Self::LowMultiMap(n_hits)
        }
    }

    fn type_code(&self) -> u8 {
        match self {
            Self::Unmapped => 15,
            Self::UniqueOnTarget => 1,
            Self::LowMultiMap(x) => *x - 1,
            Self::UniqueOffTarget => 9,
            Self::HighMultiMap => 8,
        }
    }
}

fn write_type_skip_nhits<W: Write>(w: &mut W, skip: u32, ktype: KmerType) -> std::io::Result<()> {
    let mut buf = [0u8; 8];

    let mut s = skip;
    if s < 0x0f {
        buf[0] = ((s as u8) << 4) | ktype.type_code();
        w.write_all(&buf[0..1])
    } else {
        buf[0] = 0xf0 | ktype.type_code();
        s -= 0x0f;
        if s < 0xff {
            buf[1] = s as u8;
            w.write_all(&buf[0..2])
        } else {
            buf[1] = 0xff;
            s -= 0xff;
            if s < 0xffff {
                buf[2..4].copy_from_slice(&(s as u16).to_le_bytes());
                w.write_all(&buf[0..4])
            } else {
                s -= 0xffff;
                buf[2] = 0xff;
                buf[3] = 0xff;
                buf[4..].copy_from_slice(&s.to_le_bytes());
                w.write_all(&buf)
            }
        }
    }
}

fn write_kmer_block<W: Write>(
    w: &mut W,
    v: &KmerVec,
    skip: u32,
    ktype: KmerType,
) -> anyhow::Result<()> {
    write_type_skip_nhits(w, skip, ktype)
        .with_context(|| "Error writing type, skip and nhits for kmer")?;

    if matches!(ktype, KmerType::UniqueOnTarget | KmerType::LowMultiMap(_)) {
        for x in v {
            if *x == 0 {
                break;
            }
            assert_eq!(*x & 0xf0000000, 0);
            let ix = *x - 1;
            w.write_all(&ix.to_le_bytes())
                .with_context(|| "Failed to write out kmer hit")?;
        }
    }
    Ok(())
}

fn write_kmer_blocks<W: Write>(w: &mut W, kmers: &[KmerVec]) -> anyhow::Result<()> {
    let mut prev = 0;
    for (k, v) in kmers.iter().enumerate() {
        let kmer = k as u32;
        let ktype = KmerType::from_kmer_vec(v);
        if ktype != KmerType::Unmapped {
            write_kmer_block(w, v, kmer - prev, ktype)?;
            prev = kmer
        }
    }
    Ok(())
}

fn write_close<W: Write>(w: &mut W, rnd_id: u32) -> anyhow::Result<()> {
    let mut buf: [u8; 8] = [0, 0, 0, 0, b'V', b'C', b'M', b'K'];

    u32_to_buf(&mut buf[0..4], rnd_id);
    w.write_all(&buf)
        .with_context(|| "Error writing closing block to kmer file")
}
pub fn output_kmers<P: AsRef<Path>>(
    path: P,
    reg: &Regions,
    k_work: &KmerWork,
) -> anyhow::Result<()> {
    let mut w = CompressIo::new()
        .path(path)
        .fix_path()
        .ctype(CompressType::Zstd)
        .cthreads(CompressThreads::NPhysCores)
        .bufwriter()
        .with_context(|| "Could not open kmer file for output")?;

    let rnd_id: u32 = random();
    let hdr = KmcvHeader::new(reg, k_work, rnd_id);
    hdr.write(&mut w)?;

    // Write contig blocks
    write_contig_blocks(&mut w, reg)?;

    // Write target blocks
    write_target_blocks(&mut w, reg)?;

    // write kmer blocks
    write_kmer_blocks(&mut w, k_work.kmers())?;

    write_close(&mut w, rnd_id)?;
    w.flush()
        .with_context(|| "Error flushing data to kmer file")
}
