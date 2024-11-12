# Fastaq
Manipulate FASTA and FASTQ files

[![Build Status](https://travis-ci.org/sanger-pathogens/Fastaq.svg?branch=master)](https://travis-ci.org/sanger-pathogens/Fastaq)   
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/Fastaq/blob/master/LICENSE)   

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
    * [Using pip3](#using-pip3)
    * [From source](#from-source)
    * [Running the tests](#running-the-tests)
  * [Usage](#usage)
    * [Examples](#examples)
    * [Available commands](#available-commands)
    * [For developers](#for-developers)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)

## Introduction
Python3 script to manipulate FASTA and FASTQ (and other format) files, plus API for developers

## Installation
There are a number of ways to install Fastaq and details are provided below. If you encounter an issue when installing Fastaq please contact your local system administrator. If you encounter a bug please log it [here](https://github.com/sanger-pathogens/Fastaq/issues) or email us at path-help@sanger.ac.uk.

### Pip install

Install from PyPi

```bash
pip3 install pyfastaq
```

Or pip install the latest development version directly from this repo.

```bash
pip3 install git+https://github.com/sanger-pathogens/Fastaq.git
```

### From source

If you want to edit the codebase, clone this repo and install in editable mode.

```bash
# Clone and install from this repository:
git clone https://github.com/sanger-pathogens/Fastaq.git && cd Fastaq && pip install -e ".[tests]"
```

### Running the tests

The test can be run from the top level directory:  

`pytest tests`

### Runtime dependencies

These must be available in your path at run time:
  * samtools 0.1.19
  * gzip
  * gunzip

## Usage

The installation will put a single script called `fastaq` in your path.
The usage is:

`fastaq <command> [options]`

Key points:
 * To list the available commands and brief descriptions, just run `fastaq`
 * Use `fastaq command -h` or `fastaq command --help` to get a longer description and the usage of that command.
 * The type of input file is automatically detected. Currently supported:
   `FASTA`, `FASTQ`, `GFF3`, `EMBL`, `GBK`, `Phylip`.
 * `fastaq` only manipulates sequences (and
   quality scores if present), so annotation is ignored where present in the input.
 * Input and output files can be gzipped. An input file is assumed to be gzipped if its name ends with .gz. To gzip an output file, just name it with .gz at the end.
 * You can use a minus sign for a filename to use stdin or stdout, so commands can be piped together. See the example below.

### Examples

Reverse complement all sequences in a file:

`fastaq reverse_complement in.fastq out.fastq`

Reverse complement all sequences in a gzipped file, then translate each sequence:

`fastaq reverse_complement in.fastq.gz - | fastaq translate - out.fasta`


### Available commands

| Command               | Description                                                          |
|-----------------------|----------------------------------------------------------------------|
| acgtn_only            | Replace every non acgtnACGTN with an N                               |
| add_indels            | Deletes or inserts bases at given position(s)                        |
| caf_to_fastq          | Converts a CAF file to FASTQ format                                  |
| capillary_to_pairs    | Converts file of capillary reads to paired and unpaired files        |
| chunker               | Splits sequences into equal sized chunks                             |
| count_sequences       | Counts the sequences in input file                                   |
| deinterleave          | Splits interleaved paired file into two separate files               |
| enumerate_names       | Renames sequences in a file, calling them 1,2,3... etc               |
| expand_nucleotides    | Makes every combination of degenerate nucleotides                    |
| fasta_to_fastq        | Convert FASTA and .qual to FASTQ                                     |
| filter                | Filter sequences to get a subset of them                             |
| get_ids               | Get the ID of each sequence                                          |
| get_seq_flanking_gaps | Gets the sequences flanking gaps                                     |
| interleave            | Interleaves two files, output is alternating between fwd/rev reads   |
| make_random_contigs   | Make contigs of random sequence                                      |
| merge                 | Converts multi sequence file to a single sequence                    |
| replace_bases         | Replaces all occurrences of one letter with another                  |
| reverse_complement    | Reverse complement all sequences                                     |
| scaffolds_to_contigs  | Creates a file of contigs from a file of scaffolds                   |
| search_for_seq        | Find all exact matches to a string (and its reverse complement)      |
| sequence_trim         | Trim exact matches to a given string off the start of every sequence |
| sort_by_name          | Sorts sequences in lexographical (name) order                        |
| sort_by_size          | Sorts sequences in length order                                      |
| split_by_base_count   | Split multi sequence file into separate files                        |
| strip_illumina_suffix | Strips /1 or /2 off the end of every read name                       |
| to_fake_qual          | Make fake quality scores file                                        |
| to_fasta              | Converts a variety of input formats to nicely formatted FASTA format |
| to_mira_xml           | Create an xml file from a file of reads, for use with Mira assembler |
| to_orfs_gff           | Writes a GFF file of open reading frames                             |
| to_perfect_reads      | Make perfect paired reads from reference                             |
| to_random_subset      | Make a random sample of sequences (and optionally mates as well)     |
| to_tiling_bam         | Make a BAM file of reads uniformly spread across the input reference |
| to_unique_by_id       | Remove duplicate sequences, based on their names. Keep longest seqs  |
| translate             | Translate all sequences in input nucleotide sequences                |
| trim_Ns_at_end        | Trims all Ns at the start/end of all sequences                       |
| trim_contigs          | Trims a set number of bases off the end of every contig              |
| trim_ends             | Trim fixed number of bases of start and/or end of every sequence     |
| version               | Print version number and exit                                        |


### For developers

Here is a template for counting the sequences in a FASTA or FASTQ file:

```python
from pyfastaq import sequences
seq_reader = sequences.file_reader(infile)
count = 0
for seq in seq_reader:
    count += 1
print(count)
```

Hopefully you get the idea and there are plenty of examples in tasks.py. Detection of the input file type and whether gzipped or not is automatic. See help(sequences) for the various methods already defined in the classes Fasta and Fastq.

## License
Fastaq is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/Fastaq/blob/master/LICENSE).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/sanger-pathogens/Fastaq/issues) or email path-help@sanger.ac.uk.
