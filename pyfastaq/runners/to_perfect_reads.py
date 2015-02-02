import argparse
import random
from math import floor, ceil
import sys
from pyfastaq import sequences, utils

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Makes perfect paired end fastq reads from a sequence file, with insert sizes sampled from a normal distribution. Read orientation is innies. Output is an interleaved FASTQ file.',
        usage = 'fastaq to_perfect_reads [options] <infile> <outfile> <mean insert size> <insert std deviation> <mean coverage> <read length>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output file')
    parser.add_argument('mean_insert', type=int, help='Mean insert size of read pairs', metavar='mean insert size')
    parser.add_argument('insert_std', type=float, help='Standard devation of insert size', metavar='insert std deviation')
    parser.add_argument('coverage', type=float, help='Mean coverage of the reads', metavar='mean coverage')
    parser.add_argument('readlength', type=int, help='Length of each read', metavar='read length')
    parser.add_argument('--fragments', help='Write FASTA sequences of fragments (i.e. read pairs plus sequences in between them) to the given filename', metavar='FILENAME')
    parser.add_argument('--no_n', action='store_true', help='Don\'t allow any N or n characters in the reads')
    parser.add_argument('--seed', type=int, help='Seed for random number generator. Default is to use python\'s default', default=None, metavar='INT')
    options = parser.parse_args()

    random.seed(a=options.seed)

    seq_reader = sequences.file_reader(options.infile)
    fout = utils.open_file_write(options.outfile)
    pair_counter = 1

    if options.fragments:
        fout_frags = utils.open_file_write(options.fragments)

    for ref in seq_reader:
        # check if current seq is long enough
        if len(ref) < options.mean_insert + 4 * options.insert_std:
            print('Warning, sequence ', ref.id, ' too short.  Skipping it...', file=sys.stderr)
            continue

        # work out how many reads to simulate
        read_pairs = int(0.5 * options.coverage * len(ref) / options.readlength)

        # it's possible that we pick the same fragment twice, in which case the
        # reads would get the same name. So remember the frag coords
        used_fragments = {}  # (middle_position, length) => count

        # do the simulation:  pick insert size from normal distribution, and
        # position in genome from uniform distribution
        x = 0
        while x < read_pairs:
            isize = int(random.normalvariate(options.mean_insert, options.insert_std))
            while isize > len(ref) or isize < options.readlength:
                isize = int(random.normalvariate(options.mean_insert, options.insert_std))
            middle_pos = random.randint(ceil(0.5 *isize), floor(len(ref) - 0.5 * isize))
            read_start1 = int(middle_pos - ceil(0.5 * isize))
            read_start2 = read_start1 + isize - options.readlength

            readname = ':'.join([ref.id, str(pair_counter), str(read_start1+1), str(read_start2+1)])

            fragment = (middle_pos, isize)
            if fragment in used_fragments:
                used_fragments[fragment] += 1
                readname += '.dup.' + str(used_fragments[fragment])
            else:
                used_fragments[fragment] = 1

            read1 = sequences.Fastq(readname + '/1', ref.seq[read_start1:read_start1 + options.readlength], 'I' * options.readlength)
            read2 = sequences.Fastq(readname + '/2', ref.seq[read_start2:read_start2 + options.readlength], 'I' * options.readlength)


            if options.no_n and ('n' in read1.seq or 'N' in read1.seq or 'n' in read2.seq or 'N' in read2.seq):
                continue

            read2.revcomp()

            print(read1, file=fout)
            print(read2, file=fout)

            if options.fragments:
                frag = sequences.Fasta(readname, ref.seq[read_start1:read_start2 + options.readlength])
                print(frag, file=fout_frags)

            pair_counter += 1
            x += 1

    utils.close(fout)
    if options.fragments:
        utils.close(fout_frags)
