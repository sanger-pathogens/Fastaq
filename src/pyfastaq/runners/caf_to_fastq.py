import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Converts CAF file to FASTQ format',
        usage = 'fastaq caf_to_fastq [options] <infile> <outfile>')
    parser.add_argument('infile', help='Name of input CAF file.')
    parser.add_argument('outfile', help='Name of output FASTQ file')
    parser.add_argument('-c', '--clip', action='store_true', help='Use clipping info to clip reads, if present in the input CAF file (as lines of the form "Clipping QUAL start end"). Default is to not clip')
    parser.add_argument('-l', '--min_length', type=int, help='Minimum length of sequence to output [%(default)s]', default=1, metavar='INT')
    options = parser.parse_args()

    tasks.caf_to_fastq(
        options.infile,
        options.outfile,
        trim=options.clip,
        min_length=options.min_length
    )
