## analyze_ref_gc
Analyze distribution GC content of reference genome in FASTA format for a range of read lengths.  Intended to allow comparison of GC content of short read, fixed length,  NGS data with expectations under a model of random fragment selection.

- [Introduction](#intro)
- [Installation](#install)
- [General usage](#usage)
  - [Command line options](#cli)

## <a name="intro"></a>Introduction

Given a reference genome and a set of read lengths, analyze_ref_gc will generate all sub-sequences with the given lengths.  For each subset which has
more than a given (configurable) percentage of ACTG bases (i.e., not Ns)
the number of gc and non-gc bases is stored in a hash table. and a count is made of how often each
combination occurs in the genome.  A separate hash table is used for each desired 
read length.  Given this information, the exact expected distribution og GC content in reads of a given length
can be calculated.

For bisulfite (or enzymatically) converted datasets, the GC content of the reads will be
distorted by the conversion.  Most protocols result in two types of reads, those with the majority
of Cs converted to Ts (C depleted), and those with Gs converted to As (G depleted).
With non-converted reads, we assess GC content by looking at C+G vs A+T, whereas for converted reads
we need to look at G vs A in C depleted reads, and C vs T in G depleted reads.
To generate the expected distribution for converted datasets we therefore compare G:A and C:T counts, rather than
(C+G):(A+T) as in the normal case.  If the bisulfite option is active (which is the default) then
the expected distributions for both normal and converted reads are generated.

## <a name="install"></a>Installation

To compile you will need an up-to-date copy of rust.  This can be
installed locally following the instructions [here](https://www.rust-lang.org/learn/get-started).  
Note that if you have rust already installed you should update it
using ``rustup update`` before trying to compile analyze_ref_gc.

Clone the repository and then from the baldur directory
use cargo to compile the application:
```
    git clone https://github.com/heathsc/analyze_ref_gc.git
    cd baldur
    cargo build --release
```
After successful the executable will be found in target/release/.  It
should be copied somewhere where it can be found by the shell.  
Once installed, basic help can be found by invoking analyze_ref_gc with
the -h flag.

## <a name="usage"></a>General usage

analyze_ref_gc is invoked with a reference fasta file, which can be compressed or uncompressed.  Handled compression formats
depend on what compression software is installed, but the program will recognize and handle (if the 
appropriate software is installed) compress, gzip, bgzip, bzip, and xz).

### <a name="cli"></a>Command line options

analyze_ref_gc has several command line options for controlling the operation process.

| Short | Long         | Description                                           | Default                   |
|-------|--------------|-------------------------------------------------------|---------------------------|
|       |              |                                                       |                           |
| T     | threshold    | Minimum proportion of valid bases                     | 0.8                       |
| r     | read-lengths | Set read lengths to analyze                           | 50 75 100 150 200 250 300 |
|       | no-bisulfite | Do not analyze bisulfite converted genome             | false                     |
| p     | prefix       | Set prefix for output names                           | analyze_gc                |
| i     | identifier   | Set identifier for reference                          |                           |
| t     | threads      | Set number of threads to use                          | No of cores               |
| l     | loglevel     | Set log level (none, error, warn, info, debug, trace) | info                      |
| V     | version      | Display version number and exit                       |                           |
| h     | help         | Display help text and exit                            |                           |
|       | quiet        | Silence all output to stderr                          | false                     |

#
# Changes
0.3.0 - Slight tweaks to JSON output format  
0.2.0 - Add bisulfite reference capability  
0.1.0 - Initial commit