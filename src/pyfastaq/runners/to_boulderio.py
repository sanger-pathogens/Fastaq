import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Converts input sequence file into "Boulder-IO" format, which is used by primer3',
        usage = 'fastaq to_boulderio <infile> <outfile>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output files')
    options = parser.parse_args()
    tasks.to_boulderio(options.infile, options.outfile)

