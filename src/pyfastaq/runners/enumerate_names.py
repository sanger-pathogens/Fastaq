import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
    description = 'Renames sequences in a file, calling them 1,2,3... etc',
    usage = 'fastaq enumerate_names [options] <infile> <outfile>')
    parser.add_argument('--start_index', type=int, help='Starting number [%(default)s]', default=1)
    parser.add_argument('--rename_file', help='If used, will write a file of old name to new name')
    parser.add_argument('--keep_suffix', action='store_true', help='Use this to keep a /1 or /2 suffix at the end of each name')
    parser.add_argument('--suffix', help='Add the given string to the end of every name', default=None)
    parser.add_argument('infile', help='Name of fasta/q file to be read')
    parser.add_argument('outfile', help='Name of output fasta/q file')
    options = parser.parse_args()
    tasks.enumerate_names(options.infile,
        options.outfile,
        start_index=options.start_index,
        keep_illumina_suffix=options.keep_suffix,
        rename_file=options.rename_file,
        suffix=options.suffix)
