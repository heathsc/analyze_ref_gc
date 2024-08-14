/// kmcv is a compact ninary output format with selected hashes for coverage estimation.  All multibyte constants stores in network order
///
/// HEADER
///
/// magic: [u8; 4]  - "KMCV"
/// major: u8  - Version
/// minor: u8
/// kmer length: u8
/// unused: u8 (should be set to zero)
/// rnd_id: u32
/// n_targets: u32
/// n_kmers: u64
///
/// TARGET BLOCK (repeated n_targets times)
///
/// n_kmers: u32
/// kmers: [u8;...]  KMERS are encoded as 2 bits per base, first bases in MSB. Unused bits should be set to zero. Bytes stored in network order.
///
/// EOF
/// rnd_id: u32 (should be the same as in the header)
/// magic: [u8; 4] - "VCMK"
pub mod output;
pub use output::output_kmers;
