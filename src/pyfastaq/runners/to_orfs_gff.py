import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Writes a GFF file of open reading frames from a sequence file',
        usage = 'fastaq to_orfs_gff [options] <infile> <outfile>')
    parser.add_argument('--min_length', type=int, help='Minimum length of ORF, in nucleotides [%(default)s]', default=300, metavar='INT')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output GFF file')
    options = parser.parse_args()
    tasks.fastaq_to_orfs_gff(options.infile, options.outfile, min_length=options.min_length)
