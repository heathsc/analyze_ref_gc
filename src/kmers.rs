use std::{collections::VecDeque, fmt, num::NonZeroU32};

use crate::reader::Base;

pub type KType = u32;
pub const KMER_LENGTH: usize = 15;
pub const MAX_HITS: usize = 8;
pub type KmerVec = [u32; MAX_HITS];

pub struct KmerWork {
    kmers: Vec<KmerVec>,
    max_region: usize,
}

impl fmt::Display for KmerWork {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Number of kmers: {}", self.kmers.len())
    }
}

impl KmerWork {
    #[inline]
    pub fn new() -> Self {
        Self {
            kmers: vec![[0; MAX_HITS]; 1 << (KMER_LENGTH << 1)],
            // Get maximum region (note regions start from 1)
            max_region: (1 << 31) - 2,
        }
    }
    pub fn add_kmer(&mut self, kmer: KType, region: Option<NonZeroU32>) {
        let r: u32 = region.map(|x| x.into()).unwrap_or(0);
        assert!(r as usize <= self.max_region, "Region id too large!");

        let km = kmer as usize;
        let mut set_mm = true;
        for x in self.kmers[km].iter_mut() {
            if *x == 0 {
                *x = r + 1;
                set_mm = false;
                break;
            } else if (*x == r + 1) || (*x & 0x80000000) != 0 {
                set_mm = false;
                break;
            }
        }
        if set_mm {
            self.kmers[km] = [0x80000000, 0, 0, 0, 0, 0, 0, 0]
        }
    }

    pub fn kmers(&self) -> &[KmerVec] {
        &self.kmers
    }
}

/// Returns (x, valid)
/// Where x is 0, 1, 2, 3 for A, C, T, G and 0 otherwise (with valid being false)
fn decode_base(b: Base) -> (u8, u8) {
    let b = b as u8;
    (b & 3, ((b & 4) >> 2) ^ 1)
}

pub struct KmerBuilder {
    target_vec: VecDeque<Option<NonZeroU32>>,
    kmer: KType,
    rev_kmer: KType,
    valid: KType,
    mask: KType,
    valid_mask: KType,
}

impl KmerBuilder {
    const REV_SHIFT: usize = (KMER_LENGTH - 1) << 1;
    pub fn new() -> Self {
        let k = KMER_LENGTH;
        let nb = KType::BITS as usize;

        assert!(k + k <= nb, "Kmer length too large for KmerType");
        let mut target_vec = VecDeque::with_capacity(k);
        for _ in 0..k {
            target_vec.push_back(None)
        }

        const ZERO: KType = 0;
        Self {
            target_vec,
            kmer: 0,
            rev_kmer: 0,
            valid: 0,
            mask: (!ZERO) >> (nb - k - k),
            valid_mask: (!ZERO) >> (nb - k),
        }
    }

    pub fn clear(&mut self) {
        for p in self.target_vec.iter_mut() {
            *p = None
        }
        self.valid = 0;
        self.kmer = 0;
        self.rev_kmer = 0;
    }

    pub fn add_base(&mut self, base: Base, region_idx: Option<NonZeroU32>) {
        let _ = self.target_vec.pop_front().unwrap();
        let (x, valid) = decode_base(base);
        let rev_x = (x + 2) & 3;
        self.target_vec.push_back(region_idx);
        self.kmer = ((self.kmer << 2) & self.mask) | (x as KType);
        self.rev_kmer = (self.rev_kmer >> 2) | ((rev_x as KType) << Self::REV_SHIFT);
        self.valid = ((self.valid << 1) & self.valid_mask) | (valid as KType);
    }

    pub fn target_idx(&self) -> Option<NonZeroU32> {
        let mut itr = self.target_vec.iter().copied();
        let first = itr.next()?;
        if first.is_some() {
            if itr.all(|i| i == first) {
                first
            } else {
                None
            }
        } else {
            None
        }
    }

    /// Check if kmer is composed entirely of valid (i.e., A, C, G, T) bases
    #[inline]
    pub fn kmers(&self) -> Option<[KType; 2]> {
        if self.valid == self.valid_mask {
            Some([self.kmer, self.rev_kmer])
        } else {
            None
        }
    }
}
