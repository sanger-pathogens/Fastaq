import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = description,
        usage = 'fastaq merge [options] <infile> <outfile>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output file')
    parser.add_argument('-n', '--name', help='Name of sequence in output file [%(default)s]', default='union')
    options = parser.parse_args()
    tasks.merge_to_one_seq(
        options.infile,
        options.outfile,
        seqname=options.name
    )
