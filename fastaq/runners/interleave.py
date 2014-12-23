import argparse
from fastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = description,
        usage = 'fastaq interleave <fasta/q 1> <fasta/q 2> <outfile>')
    parser.add_argument('infile_1', help='Name of first input fasta/q file')
    parser.add_argument('infile_2', help='Name of second input fasta/q file')
    parser.add_argument('outfile', help='Name of output fasta/q file of interleaved reads')
    options = parser.parse_args()
    tasks.interleave(options.infile_1, options.infile_2, options.outfile)
