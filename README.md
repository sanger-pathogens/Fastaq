Fastaq
======

Python3 script to manipulate FASTA and FASTQ files, plus API for developers

Installation
------------

Run the tests:

    python3 setup.py test

Install:

    python3 setup.py install

Notes:
 * A few scripts assume that samtools is installed and in your path. This is NOT tested in the tests, because most scripts don't need it.
 * The installation will put all scripts in your path and are named fastaq_*.

Usage
-----

General usage: `fastaq <command> [options]`


Key points:
 * To list the available commands and brief descriptions, just run `fastaq`
 * Use `fastaq command -h` or `fastaq command --help` to get a longer description and the usage of that command.
 * The type of input file is automatically detected. Currently supported:
   FASTA, FASTQ, GFF3, EMBL, GBK, Phylip.
 * `fastaq` only manipulates sequences (and
   quality scores if present), so annotation is ignored where present in the input.
 * Input and output files can be gzipped. An input file is assumed to be gzipped if its name ends with .gz. To gzip an output file, just name it with .gz at the end.
 * You can use a minus sign for a filename to use stdin or stdout, so scripts can be piped together. See the example below.


Examples
--------

Reverse complement all sequences in a file:

    fastaq reverse_complement in.fastq out.fastq

Reverse complement all sequences in a gzipped file, then translate each sequence

    fastaq reverse_complement in.fastq.gz - | fastaq translate - out.fasta


Available commands
------------------

| Command               | Description                                                          |
|-----------------------|----------------------------------------------------------------------|
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
| long_read_simulate    | Simulates long reads from reference                                  |
| make_random_contigs   | Make contigs of random sequence                                      |
| merge                 | Converts multi sequence file to a single sequence                    |
| replace_bases         | Replaces all occurences of one letter with another                   |
| reverse_complement    | Reverse complement all sequences                                     |
| scaffolds_to_contigs  | Creates a file of contigs from a file of scaffolds                   |
| search_for_seq        | Find all exact matches to a string (and its reverse complement)      |
| sequence_trim         | Trim exact matches to a given string off the start of every sequence |
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


For developers
--------------

Here is a template for counting the sequences in a FASTA or FASTQ file:

    from fastaq import sequences
    seq_reader = sequences.file_reader(infile)
    count = 0
    for seq in seq_reader:
        count += 1
    print(count)

Hopefully you get the idea and there are plenty of examples in tasks.py. Detection of the input file type and whether gzipped or not is automatic. See help(sequences) for the various methods already defined in the classes Fasta and Fastq.
