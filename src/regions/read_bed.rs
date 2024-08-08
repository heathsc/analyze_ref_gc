use std::{io::BufRead, path::Path};

use compress_io::compress::CompressIo;
use anyhow::Context;

use super::{Region, Regions};

pub fn read_bed<P: AsRef<Path>>(path: P, extend: u32) -> anyhow::Result<Regions> {
    let mut rdr = CompressIo::new().path(path).bufreader().with_context(|| "Error reading in bed files with regions")?;
    debug!("Reading regions from bed file");    
    
    let mut buf = String::new();
    let mut regs = Regions::default();
    
    let mut line = 0;
    while rdr.read_line(&mut buf).with_context(|| format!("Error reading line {} from bed file", line + 1))? > 0 {
        let mut itr = buf.trim_end().split('\t');
        
        let ctg = itr.next().ok_or_else(|| anyhow!("Missing contig information at line {}", line + 1))?;    
        let start = itr.next().ok_or_else(|| anyhow!("Missing start information at line {}", line + 1))?
            .parse::<u32>().with_context(|| format!("Bad start value at line {}", line + 1))?;
        let end = itr.next().ok_or_else(|| anyhow!("Missing end information at line {}", line + 1))?
            .parse::<u32>().with_context(|| format!("Bad end value at line {}", line + 1))?;
        
        if end <= start {
            return Err(anyhow!("End values should be larger than start value at line {}", line + 1))
        }
        
        // Avoid underflow
        let start1 = start.saturating_sub(extend);
        
        regs.get_or_insert_contig_regions(ctg)
            .add_region(Region::new(start1, end + extend - start1));
        
        buf.clear();
        line += 1;
    }
    
    debug!("Read in {line} regions. Normalizing regions");
    regs.normalize();
    
    debug!("Normalizing complete with {} non-overlapping regions retained", regs.len());
    
    Ok(regs)
}