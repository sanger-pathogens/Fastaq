import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Prints the number of sequences in input file to stdout',
        usage = 'fastaq count_sequences <infile>')
    parser.add_argument('infile', help='Name of input file')
    options = parser.parse_args()
    print(tasks.count_sequences(options.infile))
