use std::{cmp::Ordering, collections::HashMap, num::NonZeroU32};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct Region {
    start: u32, // zero offset from start of contig
    size: u32,
    idx: NonZeroU32,
}

impl PartialOrd for Region {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Region {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.start.cmp(&other.start) {
            Ordering::Equal => self.size.cmp(&other.size),
            s => s,
        }
    }
}

impl Region {
    pub fn new(start: u32, size: u32, idx: NonZeroU32) -> Self {
        Self { start, size, idx }
    }

    #[inline]
    pub fn start(&self) -> u32 {
        self.start
    }

    #[inline]
    pub fn size(&self) -> u32 {
        self.size
    }

    #[inline]
    pub fn idx(&self) -> NonZeroU32 {
        self.idx
    }

    #[inline]
    pub fn end(&self) -> u32 {
        self.start + self.size
    }
}

#[derive(Default)]
pub struct ContigRegions {
    regions: Vec<Region>,
}

impl ContigRegions {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_region(&mut self, r: Region) {
        self.regions.push(r)
    }

    pub fn regions(&self) -> &[Region] {
        &self.regions
    }

    pub(super) fn sort_and_merge(&mut self, mut ix: u32) -> u32 {
        if !self.regions.is_empty() {
            let mut r = Vec::new();
            self.regions.sort_unstable();

            let mut pending: Option<Region> = None;
            for reg in self.regions.drain(..) {
                if let Some(mut p) = pending.take() {
                    if p.end() >= reg.start() {
                        if p.end() < reg.end() {
                            p.size = reg.end() - p.start
                        }
                    } else {
                        ix += 1;
                        p.idx = ix.try_into().unwrap();
                        r.push(p);
                        pending = Some(reg)
                    }
                } else {
                    pending = Some(reg)
                }
            }
            if let Some(mut p) = pending.take() {
                ix += 1;
                p.idx = NonZeroU32::try_from(ix).unwrap();
                r.push(p)
            }
            self.regions = r
        }
        ix
    }
}

#[derive(Default)]
pub struct Regions {
    hash: HashMap<Box<str>, ContigRegions>,
}

impl Regions {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn get(&self, contig: &str) -> Option<&ContigRegions> {
        self.hash.get(contig)
    }

    pub fn get_or_insert_contig_regions(&mut self, contig: &str) -> &mut ContigRegions {
        if self.hash.contains_key(contig) {
            return self.hash.get_mut(contig).unwrap();
        }

        let ctg = contig.to_owned().into_boxed_str();
        self.hash.entry(ctg).or_default()
    }

    pub fn normalize(&mut self) -> u32 {
        let mut ix = 0;
        for r in self.hash.values_mut() {
            ix = r.sort_and_merge(ix)
        }
        ix
    }

    pub fn len(&self) -> usize {
        self.hash.values().map(|r| r.regions().len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.hash.is_empty()
    }
}
