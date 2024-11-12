import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Given a file of capillary reads, makes an interleaved file of read pairs (where more than read from same ligation, takes the longest read) and a file of unpaired reads. Replaces the .p1k/.q1k part of read names to denote fwd/rev reads with /1 and /2',
        usage = 'fastaq capillary_to_pairs <infile> <outfiles prefix>')
    parser.add_argument('infile', help='Name of input fasta/q file')
    parser.add_argument('outprefix', help='Prefix of output files', metavar='outfiles prefix')
    options = parser.parse_args()
    tasks.capillary_to_pairs(options.infile, options.outprefix)

