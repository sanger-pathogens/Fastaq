import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = description,
        usage = 'fastaq replace_bases <infile> <outfile> <old> <new>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output file')
    parser.add_argument('old', help='Base to be replaced')
    parser.add_argument('new', help='Replace with this letter')
    options = parser.parse_args()
    tasks.replace_bases(options.infile, options.outfile, options.old, options.new)
