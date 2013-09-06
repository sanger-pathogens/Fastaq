Fastaq
======

Python3 scripts to manipulate FASTA and FASTQ files, plus API for developers

Installation
------------

Run the tests:

    python3 setup.py test

Install:

    python3 setup.py install

Notes:
 * A few scripts assume that samtools is installed and in your path. This is NOT tested in the tests, because most scripts don't need it.
 * The installation will put all scripts in your path and are named fastaq_*.

Scripts
-------

Key points:
 * Use -h or --help with a script to get its usage.
 * All scripts automatically detect whether the input is a FASTA or FASTQ file.
 * Input and output files can be gzipped. An input file is assumed to be gzipped if its name ends with .gz. To gzip an output file, just name it with .gz at the end.
 * You can use a minus sign for a filename to use stdin or stdout, so scripts can be piped together. See the following examples.

Reverse complement all sequences in a file:

    fastaq_reverse_complement in.fastq out.fastq

Reverse complement all sequences in a gzipped file, then translate each sequence

    fastaq_reverse_complement in.fastq.gz - | fastaq_translate - out.fasta

For developers
--------------

Here is a template for counting the sequences in a FASTA or FASTQ file:

    from fastaq import sequences
    seq_reader = sequences.file_reader(infile)
    count = 0
    for seq in seq_reader:
        count += 1
    print(count)

Hopefully you get the idea and there are plenty of examples in tasks.py. Detection of FASTA or FASTQ and gzipped or not input file 'infile' is automatic. See help(sequences) for the various methods already defined in the classes Fasta and Fastq.
