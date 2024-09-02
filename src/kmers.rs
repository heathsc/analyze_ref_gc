use std::{
    cmp,
    collections::{HashMap, VecDeque},
    fmt,
    mem::size_of,
    num::NonZeroU32,
    ops,
};

use crate::reader::Base;

pub type KmerType = u32;
pub const KMER_LENGTH: usize = 16;

pub struct KmerWork<T> {
    kmers: HashMap<KmerType, T>,
    target_unique_count: Vec<u32>,
    max_region: usize,
    on_target_unique: usize,
    off_target_unique: usize,
    non_unique: usize,
}

impl<T> fmt::Display for KmerWork<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "On Target Unique: {}, Off Target Unique: {}, Non Unique: {}",
            self.on_target_unique, self.off_target_unique, self.non_unique
        )
    }
}

impl<T> KmerWork<T> {
    #[inline]
    pub fn new(n_targets: usize) -> Self {
        let sz = size_of::<T>();
        assert!(sz > 0, "ZST cannot be used for Kmer Work");
        let target_unique_count = vec![0; n_targets];
        Self {
            kmers: HashMap::new(),
            target_unique_count,
            // Get number of bits available (knowing that we need 1 for internal use)
            max_region: (1 << ((sz << 3) - 1)) - 1,
            on_target_unique: 0,
            off_target_unique: 0,
            non_unique: 0,
        }
    }
}

impl<T> KmerWork<T>
where
    T: Copy
        + fmt::Debug
        + TryFrom<u32>
        + Into<u64>
        + Eq
        + ops::Shl<Output = T>
        + ops::Shr<Output = T>
        + ops::BitAnd<Output = T>
        + ops::BitOrAssign,
    u64: From<T>,
{
    pub fn add_kmer(&mut self, kmer: KmerType, region: Option<NonZeroU32>) {
        let r: u32 = region.map(|x| x.into()).unwrap_or(0);
        assert!(r as usize <= self.max_region, "Region id oo large!");

        let cnv = |x| {
            T::try_from(x)
                .ok()
                .expect("Conversion to T should be infallible")
        };

        let one = cnv(1);
        if if let Some(p) = self.kmers.get_mut(&kmer) {
            let zero = cnv(0);
            if (*p & one) == zero {
                if *p != zero {
                    assert!(self.on_target_unique > 0);
                    self.on_target_unique -= 1;
                    let ix: u64 = (*p >> one).into() - 1;
                    assert!(self.target_unique_count[ix as usize] > 0);
                    self.target_unique_count[ix as usize] -= 1;
                } else {
                    assert!(self.off_target_unique > 0);
                    self.off_target_unique -= 1;
                }
                self.non_unique += 1;
                *p |= one;
                true
            } else {
                false
            }
        } else {
            self.kmers.insert(kmer, cnv(r) << one);
            if region.is_none() {
                self.off_target_unique += 1;
            } else {
                self.on_target_unique += 1;
                self.target_unique_count[(r - 1) as usize] += 1;
            }
            true
        } {
            trace!("{self}");
        }
    }
    /*
       pub fn unique(&self) -> impl Iterator<Item = (KmerType, u64)> + '_ {
           self.kmers
               .iter()
               .map(|(kmer, x)| (*kmer, u64::from(*x)))
               .filter(|(_, z)| (*z & 1) == 0)
               .map(|(kmer, x)| (kmer, x >> 1))
       }

       pub fn unique_on_target(&self) -> impl Iterator<Item = (KmerType, u64)> + '_ {
           self.kmers
               .iter()
               .map(|(kmer, x)| (*kmer, u64::from(*x)))
               .filter(|(_, z)| (*z & 1) == 0 && *z > 0)
               .map(|(kmer, x)| (kmer, x >> 1))
       }

       pub fn unique_off_target(&self) -> impl Iterator<Item = KmerType> + '_ {
           self.kmers
               .iter()
               .map(|(kmer, x)| (*kmer, u64::from(*x)))
               .filter(|(km, z)| *z == 0)
               .map(|(kmer, _)| kmer)
       }

       pub fn target_unique_count(&self) -> &[u32] {
           &self.target_unique_count
       }

    */

    pub fn on_target_kmers(mut self) -> Vec<Vec<KmerType>> {
        let n = self.target_unique_count.len();

        let mut v = Vec::with_capacity(n);
        for i in self.target_unique_count {
            v.push(Vec::with_capacity(i as usize))
        }

        for (kmer, val) in self.kmers.drain() {
            let x = u64::from(val);
            if (x & 1) == 0 && x > 0 {
                let ix = ((x >> 1) - 1) as usize;
                v[ix].push(kmer)
            }
        }

        v
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
    kmer: KmerType,
    rev_kmer: KmerType,
    valid: KmerType,
    mask: KmerType,
    valid_mask: KmerType,
}

impl KmerBuilder {
    const REV_SHIFT: usize = (KMER_LENGTH - 1) << 1;
    pub fn new() -> Self {
        let k = KMER_LENGTH;
        let nb = KmerType::BITS as usize;

        assert!(k + k <= nb, "Kmer length too large for KmerType");
        let mut target_vec = VecDeque::with_capacity(k);
        for _ in 0..k {
            target_vec.push_back(None)
        }

        const ZERO: KmerType = 0;
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
        self.kmer = ((self.kmer << 2) & self.mask) | (x as KmerType);
        self.rev_kmer = (self.rev_kmer >> 2) | ((rev_x as KmerType) << Self::REV_SHIFT);
        self.valid = ((self.valid << 1) & self.valid_mask) | (valid as KmerType);
    }

    pub fn target_idx(&self) -> Option<NonZeroU32> {
        let mut itr = self.target_vec.iter().copied();
        let first = itr.next().unwrap();
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
    pub fn kmers(&self) -> Option<[KmerType; 2]> {
        if self.valid == self.valid_mask {
            Some([self.kmer, self.rev_kmer])
        } else {
            None
        }
    }

    /// Returns Some(idx) if all bases are on target to the same region idx,
    /// otherwise returns None
    pub fn get_region_idx(&self) -> Option<NonZeroU32> {
        let mut itr = self.target_vec.iter();
        match itr.next().unwrap() {
            None => None,
            i => {
                if itr.all(|j| i == j) {
                    *i
                } else {
                    None
                }
            }
        }
    }
}
