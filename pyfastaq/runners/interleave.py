import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = description,
        usage = 'fastaq interleave <infile_1> <infile_2> <outfile>')
    parser.add_argument('infile_1', help='Name of first input file')
    parser.add_argument('infile_2', help='Name of second input file')
    parser.add_argument('outfile', help='Name of output file of interleaved reads')
    options = parser.parse_args()
    tasks.interleave(options.infile_1, options.infile_2, options.outfile)
