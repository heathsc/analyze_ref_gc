/// kmcv is a compact ninary output format with selected hashes for coverage estimation.  All multibyte constants stores in little endian order
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
/// n_contigs: u32
///
/// CONTIG BLOCK (Repeated n_contig times)
///   name_length: u32
///   name: name_length * u8
///
/// TARGET BLOCK (repeated n_targets times)
///   contig_id: u32 (starts from 1)
///   start: u32
///   end: u32
///
/// KMER BLOCK
///   type_skip_nhits: u8 (see below)
///   [optional extension byte]: u8 (see below)
///   targets mapping: nmap * u32 with the ids of the targets, except in the case where the kmer maps
///     uniquely off target, in which case nothing is written. An off target hit is marked using a
///     contig id of 0.
///
/// kmers: [u8;...]  KMERS are encoded as 2 bits per base, first bases in MSB. Unused bits should be set to zero. Bytes stored in network order.
///
/// EOF block
/// rnd_id: u32 (should be the same as in the header)
/// magic: [u8; 4] - "VCMK"
///
/// NOTES for KMER Block
///
/// type_skip_nmap is given by the first 1 or 2 bytes of the block
/// Bits 0-3 of the first byte are interpreted as follows:
///   0: kmer maps once (on target)
///   1: kmer maps twice (can be on and off target)
///   ...
///   7: kmer maps 8 times
///   8: kmer maps >8 times
///   9: kmer maps once off target
///   10-14: undefined
///   15: kmer unmapped
///
/// Bits 4-7 + optional extension byte give the number of kmers to skip (a skipped kmer is
///   implicitly assumed to be unmapped)
///
/// Meaning of bits 4-7
///   0-14: This is the skip number. The extension byte does not exist
///   15: Read the extension (next) byte and use this as the skip number 15 + (0-255). If the skip
///     is more than 255 + 15 but less than 15 + 255 + 65535 then put 255 and use the next 2 bytes
///     to store the remainder. Otherwise, put 255 in the first extension byte, 65535 in the next 2
///     extension bytes and the remainder in the next 4 bytes.
///
/// i.e. if we have a skip of 4 and then a kmer that maps 3 times to targets 20, 58914 and 1001 as
///   well as mapping off target would be encoded as:
///   0x43, 0x00000014, 0x0000e622, 0x000003e9, 0x00000000
///
/// while if we have a skip of 20 and a kmer which matches uniquely off target we would have
///   0x9f, 0x05
///
/// while if we have a skip of 302 and a kmer that matches uniquely to target 426437:
///   0xff, 0xff,
///   0xf0, 0x11, 0x000681c5
///
pub mod output;
pub use output::output_kmers;
