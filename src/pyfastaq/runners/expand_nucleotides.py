import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Makes all combinations of sequences in input file by using all possibilities of redundant bases. e.g. ART could be AAT or AGT. Assumes input is nucleotides, not amino acids',
        usage = 'fastaq expand_nucleotides <infile> <outfile>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output file')
    options = parser.parse_args()
    tasks.expand_nucleotides(
        options.infile,
        options.outfile,
    )
