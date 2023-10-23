use std::io::Write;

use libm::lgamma;

use crate::process::GcRes;
pub fn lbeta(a: f64, b: f64) -> f64 {
    lgamma(a) + lgamma(b) - lgamma(a + b)
}

const BINS: usize = 1000;

pub fn write_hist<W: Write>(
    wrt: &mut W,
    read_len: &[u32],
    res: &GcRes,
    bisulfite: bool,
) -> anyhow::Result<()> {
    let l = read_len.len();

    let l2 = if bisulfite { l * 2 } else { l };

    let mut hist: Vec<_> = (0..l2)
        .map(|_| vec![0.0; BINS].into_boxed_slice())
        .collect();
    let mut lnp = Vec::with_capacity(BINS);
    let mut tmp = Vec::with_capacity(BINS);
    let mut t = vec![0.0; l2];
    let inc = 1.0 / (BINS as f64);
    for i in 0..BINS {
        let x = inc * (0.5 + (i as f64));
        lnp.push((x, x.ln(), (1.0 - x).ln()))
    }
    for (ix, h) in hist.iter_mut().enumerate() {
        let gc_hist = res.get_gc_hist(read_len[ix % l]).unwrap();
        let hash = if ix < l {
            gc_hist.hash()
        } else {
            gc_hist.bisulfite_hash().unwrap()
        };
        for (b, a, x) in hash.iter().map(|(ct, x)| {
            let (r, s) = ct.counts();
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
                *q += x * p / z
            }
        }
    }
    let scale = BINS as f64;
    write!(wrt, "gc")?;
    for l in read_len {
        write!(wrt, "\tread_len:{}bp", l)?;
        if bisulfite {
            write!(wrt, "\tbisulfite_read_len:{}bp", l)?
        }
    }
    writeln!(wrt)?;
    for i in 0..BINS {
        write!(wrt, "{}", lnp[i].0)?;
        for j in 0..l {
            let h = &hist[j];
            write!(wrt, "\t{}", h[i] * scale / t[j])?;
            if bisulfite {
                let h = &hist[j + l];
                write!(wrt, "\t{}", h[i] * scale / t[j + l])?;
            }
        }
        writeln!(wrt)?
    }
    Ok(())
}
