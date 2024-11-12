import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = description,
        usage = 'fastaq fasta_to_fastq <fasta in> <qual in> <fastq out>')
    parser.add_argument('fasta', help='Name of input FASTA file', metavar='fasta in')
    parser.add_argument('qual', help='Name of input quality scores file', metavar='qual in')
    parser.add_argument('outfile', help='Name of output FASTQ file', metavar='fastq out')
    options = parser.parse_args()
    tasks.fasta_to_fastq(options.fasta, options.qual, options.outfile)
