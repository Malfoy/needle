<!--
SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

# Needle [![build status][1]][2] [![codecov][3]][4] [![install with bioconda][5]][6]

[1]: https://img.shields.io/github/actions/workflow/status/seqan/needle/ci_linux.yml?branch=main&style=flat&logo=github&label=CI "Open GitHub actions page"
[2]: https://github.com/seqan/needle/actions?query=branch%3Amain
[3]: https://codecov.io/gh/seqan/needle/graph/badge.svg?token=W109QS58E0 "Open Codecov page"
[4]: https://codecov.io/gh/seqan/needle
[5]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
[6]: #install-with-bioconda-linux

### A fast and space-efficient pre-filter for estimating the quantification of very large collections of nucleotide sequences

Needle is a tool for semi-quantitative analysis of very large collections of nucleotide sequences.

Needle stores its data in multiple Interleaved Bloom Filter (IBF), a fast and space efficient probabilistic data structure and uses a windowing scheme (also called minimisers) to reduce the amount of data to store. How many Interleaved Bloom Filter are used is defined by the user.

Each IBF has a so-called expression threshold and stores minimisers with an occurrence greater than or equal to its own expression threshold and smaller than the next biggest expression threshold (if there is no bigger expression threshold, all greater than or equal to the threshold are stored). These expression thresholds are then used during the query (called estimate) to approximate the expression values of given transcripts.

## Citation

In your academic works (also comparisons and pipelines) please cite:
  * *Needle: a fast and space-efficient prefilter for estimating the quantification of very large collections of expression experiments*;
    Mitra Darvish, Enrico Seiler, Svenja Mehringer, René Rahn, and Knut Reinert; Bioinformatics, Volume 38, Issue 17, 1 September 2022, Pages 4100–4108.
    doi: https://doi.org/10.1093/bioinformatics/btac492

## Download and Installation

### From Source

**Prerequisites**:
* CMake >= 3.25
* GCC >= 12
* LLVM Clang >= 17
* Intel oneAPI DPC++/C++ Compiler >= 2025.0 (older versions may work, but are not tested)
* git

Refer to the [Seqan3 Setup Tutorial](https://docs.seqan.de/seqan3/main_user/setup.html) for more in depth
information.

Needle can be built by following these commands:

```bash
git clone https://github.com/seqan/needle.git
mkdir build-needle && cd build-needle
cmake ../needle -DCMAKE_BUILD_TYPE=Release
make
```

Run tests to check if Needle is working as intended. All tests should pass.

```bash
make check
```

### Install with [bioconda](https://bioconda.github.io/recipes/needle/README.html) (Linux)

```bash
conda install -c bioconda -c conda-forge needle
```

## Usage

### Build a Needle index
To build a Needle index, several sequence files have to be given. All sequence file formats supported by SeqAn3 are accepted as an input (fasta, fastq, embl,... and their compressed forms).

The flag `--paired` in the example below indicates that the given sequence files are paired-end experiments. Furthermore, the false positive rate has to be specified with the parameter `f`.

Use -h/--help for more information and to see further parameters. The flag `-c` can be used to build a compressed Needle index.

The following example creates a compressed Needle index for two paired-end experiments for the expression thresholds 4 and 32.

```bash
./bin/needle ibf ../needle/test/data/exp_*.fasta --paired -e 16 -e 32 -f 0.3 -c -o example
```

Even though this works, it is recommended to calculate the minimisers beforehand by using the option `minimisers`. It calculates the minimisers of given experiments and stores their hash values and their occurrences in a binary file named ".minimiser".

The following command calculates the minimisers in the two experiments.
```bash
./bin/needle minimiser ../needle/test/data/exp_*.fasta --paired
```

A minimiser file is a binary file containing the following data:
- number of minimisers (uint64_t)
- kmer-size (uint8_t)
- window-size (uint32_t)
- seed (uint64_t)
- flag which is true, if shape is ungapped (bool)
- shape (uint64_t), if flag is false
- all minimiser hashes (uint64_t) with their occurrences (uint16_t)

Based on the minimiser files, the Needle index can be computed by using the following command:
```bash
./bin/needle ibfmin exp*.minimiser -e 16 -e 32  -f 0.3 -c -o example
```

### Build from unitigs (Logan, GGCAT, or Cuttlefish)

`needle minimiser --unitigs` imports precomputed mean k-mer abundances from unitig
FASTA headers, bypassing occurrence counting in Needle. Each header must contain
`ka:f:<abundance>` or `km:f:<abundance>`, for example:

```text
>unitig_0 ka:f:12.5
ACGTTGCAACGTTGCAACGTTGCAACGTTGCA
```

Use one unitig file per sample; each file becomes one index bin. Assemble samples
separately when sample-specific abundances are needed. A pooled graph's single
abundance per unitig cannot represent separate sample abundances.

The examples below use k=31. Set Needle's `-k` to the assembly k and `-w` to the same
value to retain every k-mer. Hashing parameters otherwise keep their usual defaults.
The generated `.minimiser` files use Needle's existing format and can be passed to
`ibfmin` or `insertmin`.

#### Logan

[Logan unitigs](https://github.com/IndexThePlanet/Logan/blob/main/Unitigs.md) are
assembled with k=31 and already carry mean abundance in `ka:f:` (or, for some
accessions, `km:f:`). Download an accession and decompress its `.zst` file first:

```sh
wget https://s3.amazonaws.com/logan-pub/u/SRR11905265/SRR11905265.unitigs.fa.zst
zstd -d SRR11905265.unitigs.fa.zst
needle minimiser SRR11905265.unitigs.fa --unitigs -k 31 -w 31
needle ibfmin SRR11905265.unitigs.minimiser -e 2 -e 5 -e 10 -f 0.05 -o logan_index
```

#### GGCAT

Build [GGCAT with its optional `kmer-counters` feature](https://github.com/algbio/ggcat#additional-opt-in-features)
to include mean abundance in the `km:f:` FASTA tag. For example, with Rust and Cargo
installed:

```sh
git clone https://github.com/algbio/ggcat.git
cargo install --path ggcat/crates/cmdline --locked --features kmer-counters
```

Assemble one sample's reads into an uncompressed FASTA file, then import it:

```sh
ggcat build -k 31 -j 8 sample_R1.fastq.gz sample_R2.fastq.gz -o sample.ggcat.unitigs.fa
needle minimiser sample.ggcat.unitigs.fa --unitigs -k 31 -w 31
needle ibfmin sample.ggcat.unitigs.minimiser -e 2 -e 5 -e 10 -f 0.05 -o ggcat_index
```

Needle reads the mean `km:f:` value, not the total `KC:i:` count. GGCAT output
without the abundance feature is insufficient for this mode. Use FASTA output;
decompress existing `.lz4` output with `lz4 -d` before importing it.

#### Cuttlefish

Use an abundance-enabled Cuttlefish2, such as the
[modified version used by Logan](https://github.com/rchikhi/cuttlefish), which writes
`ka:f:` tags. Standard Cuttlefish output without abundance annotations cannot be
used directly with `--unitigs`: the original read abundance cannot be recovered
from unitig sequences alone.

With the abundance-enabled `cuttlefish` executable on your PATH:

```sh
cuttlefish build --read -s sample.fastq.gz -k 31 -t 8 -o sample.cuttlefish.unitigs
needle minimiser sample.cuttlefish.unitigs.fa --unitigs -k 31 -w 31
needle ibfmin sample.cuttlefish.unitigs.minimiser -e 2 -e 5 -e 10 -f 0.05 -o cuttlefish_index
```

#### Multiple samples and abundance handling

For several samples assembled with the same k, import their unitigs together and
build one index. Every input file remains a separate sample:

```sh
needle minimiser sample1.unitigs.fa sample2.unitigs.fa --unitigs -k 31 -w 31 -t 2
needle ibfmin sample1.unitigs.minimiser sample2.unitigs.minimiser -e 2 -e 5 -e 10 -f 0.05 -o unitig_index
```

Choose expression thresholds (`-e`) appropriate to your samples, or use `-l` for
automatic threshold selection. New samples can be added with `insertmin`; import
them using the same k, window, shape, and seed as the existing uncompressed index.

Abundances are unitig averages, not exact per-k-mer counts. They are rounded to the
nearest integer (half up) and capped at 65534 to match Needle's count storage. A
repeated hash keeps the maximum abundance across unitigs and grouped files;
repeated occurrences do not increase abundance. The default cutoff is 0 in this
mode, independent of file size; `--cutoff` discards rounded abundances less than or
equal to the supplied value. Include/exclude filters and sample grouping remain
available. Missing or invalid abundance tags produce an error.

### Estimate
To estimate the expression value of one transcript, a sequence file has to be given

Use the parameter "-i" to define where the Needle index can be found (should be equal with "-o" in the previous commands).

Use -h/--help for more information and to see further parameters.

The following example searches for one gene, which is expressed in the first experiment with expression 6 and in the second with expression 37. Therefore, it should be found only in the second experiment but not the first when using expression levels of 16 and 32.

```bash
./bin/needle estimate ../needle/test/data/gene.fasta -i example
```

The created file "expressions.out" (if you prefer a different name, use "-o") should contain the following:
```text
GeneA   0      32
```

### Insert into an existing Needle index
It is possible to insert new sequence files into an uncompressed Needle index.

Similar to the build step, this can be done by either using the sequence files as input directly or the minimiser files outputted by `needle minimiser`.

Most options are the same as the ones from the build step, however as the Needle index already exists, neither the false positive rate nor the number of hash functions can be changed.

It is necessary to specify `i` to the directory, where the existing Needle index can be found.

The following example inserts into the Needle index build above for two paired-end experiments.

```bash
# Create Index
./bin/needle ibf ../needle/test/data/exp_0*.fasta --paired -e 16 -e 32 -f 0.3 -c -o example
# Insert into created index
./bin/needle insert ../needle/test/data/exp_1*.fasta --paired -i example
```

Based on minimiser files, an insertion to the Needle index can be achieved by using the following command:
```bash
# Create Index
./bin/needle ibf ../needle/test/data/exp_0*.fasta --paired -e 16 -e 32 -f 0.3 -c -o example
# Insert into created index
./bin/needle insertmin exp*.minimiser -i example
```

The insert methods based on minimiser or on sequence files is independent of the way the index was created.

### Delete experiments from an existing Needle index
It is possible to delete sequence files from an uncompressed Needle index by specifying the position of the experiment, which should be deleted.

These deleted experiments won't change the size of the index, as the space is kept for later insertions.
```bash
# Create Index
./bin/needle ibf ../needle/test/data/exp_*.fasta --paired -e 16 -e 32 -f 0.3 -c -o example
# Delete first experiment exp_0 (with position 0) from index
./bin/needle delete  -i example 0
```

## Note
This app was created with the [SeqAn app-template](https://github.com/seqan/app-template).
