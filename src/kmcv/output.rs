use std::{io::Write, path::Path};

use anyhow::Context;
use compress_io::{compress::CompressIo, compress_type::CompressType};
use rand::random;

use crate::kmers::{KmerType, KMER_LENGTH};

const MAJOR_VERSION: u8 = 1;
const MINOR_VERSION: u8 = 0;

// Minimum unique kmers required for a target
const MIN_KMERS: usize = 10;

#[inline]
fn u32_to_buf(b: &mut [u8], x: u32) {
    b.copy_from_slice(&x.to_le_bytes())
}

#[inline]
fn u64_to_buf(b: &mut [u8], x: u64) {
    b.copy_from_slice(&x.to_le_bytes())
}
fn write_header<W: Write>(
    w: &mut W,
    n_targets: u32,
    n_kmers: u64,
    rnd_id: u32,
) -> anyhow::Result<()> {
    let mut buf: [u8; 24] = [
        b'K',
        b'M',
        b'C',
        b'V',
        MAJOR_VERSION,
        MINOR_VERSION,
        KMER_LENGTH as u8,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ];

    u32_to_buf(&mut buf[8..12], rnd_id);
    u32_to_buf(&mut buf[12..16], n_targets);
    u64_to_buf(&mut buf[16..24], n_kmers);
    w.write_all(&buf)
        .with_context(|| "Error writing header to kmer file")
}

fn write_target_block<W: Write>(w: &mut W, v: &[KmerType]) -> anyhow::Result<()> {
    let n = v.len() as u32;
    w.write_all(&n.to_le_bytes())
        .with_context(|| "Error writing target length")?;
    for k in v.iter() {
        w.write_all(&k.to_le_bytes())
            .with_context(|| "Error writing kmer")?;
    }
    Ok(())
}

fn write_close<W: Write>(w: &mut W, rnd_id: u32) -> anyhow::Result<()> {
    let mut buf: [u8; 8] = [0, 0, 0, 0, b'V', b'C', b'M', b'K'];

    u32_to_buf(&mut buf[0..4], rnd_id);
    w.write_all(&buf)
        .with_context(|| "Error writing closing block to kmer file")
}
pub fn output_kmers<P: AsRef<Path>>(path: P, kmers: &[Vec<KmerType>]) -> anyhow::Result<()> {
    let mut w = CompressIo::new()
        .path(path)
        .fix_path()
        .ctype(CompressType::Gzip)
        .bufwriter()
        .with_context(|| "Could not open kmer file for output")?;

    let rnd_id: u32 = random();

    let (n_targets, n_kmers) = kmers
        .iter()
        .filter(|v| v.len() >= MIN_KMERS)
        .fold((0u32, 0u64), |(nt, nk), v| (nt + 1, nk + v.len() as u64));

    write_header(&mut w, n_targets, n_kmers, rnd_id)?;

    for v in kmers.iter().filter(|v| v.len() >= MIN_KMERS) {
        write_target_block(&mut w, v)?;
    }

    write_close(&mut w, rnd_id)?;
    w.flush()
        .with_context(|| "Error flushing data to kmer file")
}
