import argparse
from fastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Removes duplicate sequences from a fasta/q file, based on their names. If the same name is found more than once, then the longest sequence is kept. Order of sequences is preserved in output',
        usage = 'fastaq to_unique_by_id <infile> <outfile>')
    parser.add_argument('infile', help='Name of input fasta/q file')
    parser.add_argument('outfile', help='Name of output fasta/q file')
    options = parser.parse_args()
    tasks.to_unique_by_id(options.infile, options.outfile)
