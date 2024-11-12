import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = description,
        usage = 'fastaq to_fasta [options] <infile> <outfile>')
    parser.add_argument('infile', help='Name of input file. Can be any of FASTA, FASTQ, GFF3, EMBL, GBK, Phylip')
    parser.add_argument('outfile', help='Name of output file')
    parser.add_argument('-l', '--line_length', type=int, help='Number of bases on each sequence line of output file. Set to zero for no linebreaks in sequences [%(default)s]', default=60)
    parser.add_argument('-s', '--strip_after_whitespace', action='store_true', help='Remove everything after first whitespace in every sequence name')
    parser.add_argument('-u', '--check_unique', action='store_true', help='Die if any of the output sequence names are not unique')
    options = parser.parse_args()

    tasks.to_fasta(
        options.infile,
        options.outfile,
        line_length=options.line_length,
        strip_after_first_whitespace=options.strip_after_whitespace,
        check_unique=options.check_unique
    )

