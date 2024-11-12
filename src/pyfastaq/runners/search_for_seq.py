import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Searches for an exact match on a given string and its reverse complement, in every sequence of input sequence file. Case insensitive. Guaranteed to find all hits',
        usage = 'fastaq search_for_seq [options] <infile> <outfile> <search_string>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of outputfile. Tab-delimited output: sequence name, position, strand')
    parser.add_argument('search_string', help='String to search for in the sequences')
    options = parser.parse_args()
    tasks.search_for_seq(options.infile, options.outfile, options.search_string)
