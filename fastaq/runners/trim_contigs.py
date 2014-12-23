import argparse
from fastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Trims a set number of bases off the end of every contig, so gaps get bigger and contig ends are removed. Bases are replaced with Ns. Any sequence that ends up as all Ns is lost',
        usage = 'fastaq trim_contigs [options] <fasta/q in> <fasta/q out>')
    parser.add_argument('--trim_number', type=int, help='Number of bases to trim around each gap, and off ends of each sequence [%(default)s]', default=100)
    parser.add_argument('infile', help='Name of input fasta/q file')
    parser.add_argument('outfile', help='Name of output fasta/q file')
    options = parser.parse_args()
    tasks.trim_contigs(options.infile, options.outfile, options.trim_number)
