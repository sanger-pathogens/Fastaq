import argparse
from fastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Writes a GFF file of open reading frames from a fasta/q file',
        usage = 'fastaq to_orfs_gff [options] <fasta/q in> <gff_out>')
    parser.add_argument('--min_length', type=int, help='Minimum length of ORF, in nucleotides [%(default)s]', default=300, metavar='INT')
    parser.add_argument('infile', help='Name of input fasta/q file')
    parser.add_argument('gff_out', help='Name of output gff file')
    options = parser.parse_args()
    tasks.fastaq_to_orfs_gff(options.infile, options.gff_out, min_length=options.min_length)
