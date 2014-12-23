import argparse
from fastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = description,
        usage = 'fastaq get_seq_flanking_gaps [options] <fasta/q in> <fasta/q out>')
    parser.add_argument('--left', type=int, help='Number of bases to get to left of gap [%(default)s]', default=25, metavar='INT')
    parser.add_argument('--right', type=int, help='Number of bases to get to right of gap [%(default)s]', default=25, metavar='INT')
    parser.add_argument('infile', help='Name of input fasta/q file')
    parser.add_argument('outfile', help='Name of output fasta/q file')
    options = parser.parse_args()
    tasks.get_seqs_flanking_gaps(options.infile, options.outfile, options.left, options.right)
