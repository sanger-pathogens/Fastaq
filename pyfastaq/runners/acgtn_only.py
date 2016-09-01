import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
    description = 'Replaces any character that is not one of acgtACGTnN with an N',
    usage = 'fastaq acgtn_only [options] <infile> <outfile>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output file')
    options = parser.parse_args()
    tasks.acgtn_only(options.infile, options.outfile)

