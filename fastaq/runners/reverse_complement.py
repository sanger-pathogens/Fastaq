import argparse
from fastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = description,
        usage = 'fastaq reverse_complement <fasta/q in> <fasta/q out>')
    parser.add_argument('infile', help='Name of input fasta/q file')
    parser.add_argument('outfile', help='Name of output fasta/q file')
    options = parser.parse_args()
    tasks.reverse_complement(options.infile, options.outfile)
