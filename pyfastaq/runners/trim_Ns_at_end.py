import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Trims any Ns off each sequence in input file. Does nothing to gaps in the middle, just trims the ends',
        usage = 'fastaq trim_Ns_at_end <infile> <outfile>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output file')
    options = parser.parse_args()
    tasks.trim_Ns_at_end(options.infile, options.outfile)
