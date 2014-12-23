import argparse
from fastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Filters a fasta/q file by sequence length and/or by name matching a regular expression',
        usage = 'fastaq filter [options] <infile> <outfile>')
    parser.add_argument('--min_length', type=int, help='Minimum length of sequence to keep [%(default)s]', default=0, metavar='INT')
    parser.add_argument('--max_length', type=float, help='Maximum length of sequence to keep [%(default)s]', default=float('inf'), metavar='INT')
    parser.add_argument('--regex', help='If given, only reads with a name matching the regular expression will be kept')
    parser.add_argument('--ids_file', help='If given, only reads whose ID is in th given file will be used. One ID per line of file.')
    parser.add_argument('-v', '--invert', action='store_true', help='Keep sequences that do not match the filters')
    parser.add_argument('infile', help='Name of fasta/q file to be filtered')
    parser.add_argument('outfile', help='Name of output fasta/q file')
    options = parser.parse_args()
    tasks.filter(options.infile,
                 options.outfile,
                 minlength=options.min_length,
                 maxlength=options.max_length,
                 regex=options.regex,
                 ids_file=options.ids_file,
                 invert=options.invert
    )
