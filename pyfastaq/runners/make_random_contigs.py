import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Makes a multi-FASTA file of random sequences, all of the same length. Each base has equal chance of being A,C,G or T',
        usage = 'fastaq make_random_contigs [options] <contigs> <length> <outfile>')
    parser.add_argument('--first_number', type=int, help='If numbering the sequences, the first sequence gets this number [%(default)s]', default=1)
    parser.add_argument('--name_by_letters', action='store_true', help='Name the contigs A,B,C,... will start at A again if you get to Z')
    parser.add_argument('--prefix', help='Prefix to add to start of every sequence name', default='')
    parser.add_argument('--seed', type=int, help='Seed for random number generator. Default is to use python\'s default', default=None)
    parser.add_argument('contigs', type=int, help='Number of contigs to make')
    parser.add_argument('length', type=int, help='Length of each contig')
    parser.add_argument('outfile', help='Name of output file')
    options = parser.parse_args()
    tasks.make_random_contigs(
        options.contigs,
        options.length,
        options.outfile,
        name_by_letters=options.name_by_letters,
        prefix=options.prefix,
        seed=options.seed,
        first_number=options.first_number
    )
