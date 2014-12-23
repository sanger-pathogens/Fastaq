import argparse
from fastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Makes fake quality scores file from a fasta/q file',
        usage = 'fastaq to_fake_qual [options] <infile> <outfile>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output file')
    parser.add_argument('-q', '--qual', type=int, help='Quality score to assign to all bases [%(default)s]', default=40)
    options = parser.parse_args()
    tasks.fastaq_to_fake_qual(
        options.infile,
        options.outfile,
        q=options.qual
    )

