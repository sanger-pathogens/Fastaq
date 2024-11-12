import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Gets IDs from each sequence in input file',
        usage = 'fastaq get_ids <infile> <outfile>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output file')
    options = parser.parse_args()
    tasks.get_ids(options.infile, options.outfile)
