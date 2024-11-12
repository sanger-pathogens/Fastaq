import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Trims sequences off the start of all sequences in a pair of sequence files, whenever there is a perfect match. Only keeps a read pair if both reads of the pair are at least a minimum length after any trimming',
        usage = 'fastaq sequence_trim [options] <infile_1> <infile_2> <outfile_1> <outfile_2> <trim_seqs>')
    parser.add_argument('--min_length', type=int, help='Minimum length of output sequences [%(default)s]', default=50, metavar='INT')
    parser.add_argument('--revcomp', action='store_true', help='Trim the end of each sequence if it matches the reverse complement. This option is intended for PCR primer trimming') 
    parser.add_argument('infile_1', help='Name of forward fasta/q file to be trimmed')
    parser.add_argument('infile_2', help='Name of reverse fasta/q file to be trimmed')
    parser.add_argument('outfile_1', help='Name of output forward fasta/q file')
    parser.add_argument('outfile_2', help='Name of output reverse fasta/q file')
    parser.add_argument('trim_seqs', help='Name of file of sequences to search for at the start of each input sequence')
    options = parser.parse_args()
    tasks.sequence_trim(
        options.infile_1,
        options.infile_2,
        options.outfile_1,
        options.outfile_2,
        options.trim_seqs,
        min_length=options.min_length,
        check_revcomp=options.revcomp
    )
