import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Translates all sequences in input file. Output is always FASTA format',
        usage = 'fastaq translate [options] <infile> <outfile>')
    parser.add_argument('--frame', type=int, choices=[0,1,2], help='Frame to translate [%(default)s]', default=0)
    parser.add_argument('infile', help='Name of file to be translated')
    parser.add_argument('outfile', help='Name of output FASTA file')
    options = parser.parse_args()
    tasks.translate(options.infile, options.outfile, frame=options.frame)
