import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = description,
        usage = 'fastaq sort_by_size [options] <infile> <outfile>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output file')
    parser.add_argument('-r', '--reverse', action='store_true', help='Sort by shortest first instead of the default of longest first')
    options = parser.parse_args()
    tasks.sort_by_size(
        options.infile,
        options.outfile,
        smallest_first=options.reverse
    )
