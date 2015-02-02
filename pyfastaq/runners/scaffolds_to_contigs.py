import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Creates a file of contigs from a file of scaffolds - i.e. breaks at every gap in the input',
        usage = 'fastaq scaffolds_to_contigs [options] <infile> <outfile>')
    parser.add_argument('--number_contigs', action='store_true', help='Use this to enumerate contig names 1,2,3,... within each scaffold')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output contigs file')
    options = parser.parse_args()
    tasks.scaffolds_to_contigs(options.infile, options.outfile, number_contigs=options.number_contigs)
