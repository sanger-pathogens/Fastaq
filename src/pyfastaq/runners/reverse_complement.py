import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = description,
        usage = 'fastaq reverse_complement <infile> <outfile>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output file')
    options = parser.parse_args()
    tasks.reverse_complement(options.infile, options.outfile)
