import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Splits a multi sequence file into separate files. Does not split sequences. Puts up to max_bases into each split file. The exception is that any sequence longer than max_bases is put into its own file.',
        usage = 'fastaq split_by_base_count [options] <infile> <outprefix> <max_bases>')
    parser.add_argument('infile', help='Name of input file to be split')
    parser.add_argument('outprefix', help='Name of output file')
    parser.add_argument('max_bases', type=int, help='Max bases in each output split file', metavar='max_bases')
    parser.add_argument('--max_seqs', type=int, help='Max number of sequences in each output split file [no limit]', metavar='INT')

    options = parser.parse_args()
    tasks.split_by_base_count(options.infile, options.outprefix, options.max_bases, options.max_seqs)
