import argparse
from fastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = description,
        usage = 'fastaq strip_illumina_suffix <fasta/q in> <fasta/q out>')
    parser.add_argument('infile', help='Name of input fasta/q file')
    parser.add_argument('outfile', help='Name of output fasta/q file')
    options = parser.parse_args()
    tasks.strip_illumina_suffix(options.infile, options.outfile)
