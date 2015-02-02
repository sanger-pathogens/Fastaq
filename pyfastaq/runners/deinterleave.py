import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
    description = 'Deinterleaves sequence file, so that reads are written alternately between two output files',
    usage = 'fastaq deinterleave [options] <infile> <out_fwd> <out_rev>')
    parser.add_argument('--fasta_out', action='store_true', help='Use this to write output as fasta (default is same as input)', default=False)
    parser.add_argument('infile', help='Name of fasta/q file to be deinterleaved')
    parser.add_argument('out_fwd', help='Name of output fasta/q file of forwards reads')
    parser.add_argument('out_rev', help='Name of output fasta/q file of reverse reads')
    options = parser.parse_args()
    tasks.deinterleave(options.infile, options.out_fwd, options.out_rev, fasta_out=options.fasta_out)
