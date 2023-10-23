use std::{collections::HashMap, io::Write};

use libm::lgamma;

use crate::process::GcRes;
pub fn lbeta(a: f64, b: f64) -> f64 {
    lgamma(a) + lgamma(b) - lgamma(a + b)
}

const BINS: usize = 1000;

pub fn write_hist<W: Write>(wrt: &mut W, read_len: &[u32], res: &GcRes) -> anyhow::Result<()> {
    let l = read_len.len();

    let mut hist: Vec<_> = (0..l).map(|_| vec![0.0; BINS].into_boxed_slice()).collect();
    let mut lnp = Vec::with_capacity(BINS);
    let mut tmp = Vec::with_capacity(BINS);
    let mut t = vec![0.0; l];
    let inc = 1.0 / (BINS as f64);
    for i in 0..BINS {
        let x = inc * (0.5 + (i as f64));
        lnp.push((x, x.ln(), (1.0 - x).ln()))
    }
    for (ix, h) in hist.iter_mut().enumerate() {
        let cts = res.get_hash(ix).unwrap();
        for (b, a, x) in cts.iter().map(|(k, x)| {
            let (r, s) = k.counts();
            (r as f64, s as f64, *x as f64)
        }) {
            t[ix] += x;
            let konst = lbeta(a + 1.0, b + 1.0);
            tmp.clear();
            let mut z = 0.0;
            for (_, lnp, lnp1) in lnp.iter() {
                let p = (lnp * a + lnp1 * b - konst).exp();
                z += p;
                tmp.push(p);
            }
            for (p, q) in tmp.iter().zip(h.iter_mut()) {
                *q += x * p
            }
        }
    }
    let scale = BINS as f64;
    write!(wrt, "gc")?;
    for l in read_len {
        write!(wrt, "\tread_len:{}bp", l)?
    }
    writeln!(wrt)?;
    for i in 0..BINS {
        write!(wrt, "{}", lnp[i].0)?;
        for (h, t) in hist.iter().zip(t.iter()) {
            write!(wrt, "\t{}", h[i] * scale / t)?;
        }
        writeln!(wrt)?
    }
    Ok(())
}
