import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = description,
        usage = 'fastaq trim_ends <infile> <bases off start> <bases off end> <outfile>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('start_trim', type=int, help='Number of bases to trim off start')
    parser.add_argument('end_trim', type=int, help='Number of bases to trim off end')
    parser.add_argument('outfile', help='Name of output file')
    options = parser.parse_args()
    tasks.trim(options.infile, options.outfile, options.start_trim, options.end_trim)
