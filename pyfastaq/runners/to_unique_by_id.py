import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Removes duplicate sequences from input file, based on their names. If the same name is found more than once, then the longest sequence is kept. Order of sequences is preserved in output',
        usage = 'fastaq to_unique_by_id <infile> <outfile>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output file')
    options = parser.parse_args()
    tasks.to_unique_by_id(options.infile, options.outfile)
