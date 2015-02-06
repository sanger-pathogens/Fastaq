import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Splits a multi sequence file into separate files. Splits sequences into chunks of a fixed size. Aims for chunk_size chunks in each file, but allows a little extra, so chunk can be up to (chunk_size + tolerance), to prevent tiny chunks made from the ends of sequences',
        usage = 'fastaq chunker [options] <infile> <out> <chunk size> <tolerance>')
    parser.add_argument('infile', help='Name of input file to be split')
    parser.add_argument('out', help='Prefix of output file. If --onefile used, then name of single output file')
    parser.add_argument('chunk_size', type=int, help='Size of each chunk')
    parser.add_argument('tolerance', type=int, help='Tolerance allowed in chunk size')
    parser.add_argument('--onefile', action='store_true', help='Output all the sequences in one file')
    parser.add_argument('--skip_all_Ns', action='store_true', help='Do not output any sequence that consists of all Ns')
    options = parser.parse_args()
    if options.onefile:
        tasks.split_by_fixed_size_onefile(
            options.infile,
            options.out,
            options.chunk_size,
            options.tolerance,
            skip_if_all_Ns=options.skip_all_Ns
        )
    else:
        tasks.split_by_fixed_size(
            options.infile,
            options.out,
            options.chunk_size,
            options.tolerance,
            skip_if_all_Ns=options.skip_all_Ns
        )
