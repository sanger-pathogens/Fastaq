import argparse
import sys
import os
from pyfastaq import sequences, utils

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Takes a sequence file. Makes a BAM file containing perfect (unpaired) reads tiling the whole genome',
        usage = 'fastaq to_tiling_bam [options] <infile> <read_length> <read_step> <read_prefix> <outfile>',
        epilog = 'Important: assumes that samtools is in your path')
    parser.add_argument('infile', help='Name of input fasta/q file')
    parser.add_argument('read_length', type=int, help='Length of reads')
    parser.add_argument('read_step', type=int, help='Distance between start of each read')
    parser.add_argument('read_prefix', help='Prefix of read names')
    parser.add_argument('outfile', help='Name of output BAM file')
    parser.add_argument('--qual_char', help='Character to use for quality score [%(default)s]', default='I')
    parser.add_argument('--read_group', help='Add the given read group ID to all reads [%(default)s]' ,default='42')
    options = parser.parse_args()

    # make a header first  - we need to add the @RG line to the default header made by samtools
    tmp_empty_file = options.outfile + '.tmp.empty'
    f = utils.open_file_write(tmp_empty_file)
    utils.close(f)
    try:
        f = os.popen('samtools view -H -T ' + options.infile + ' ' + tmp_empty_file)
    except IOError:
        print('Error making tmp header file', file=sys.stderr)
        sys.exit(1)

    header_lines = f.readlines()
    header_lines.append('@RG\tID:' + options.read_group + '\tSM:FAKE')
    f.close()
    os.unlink(tmp_empty_file)

    seq_reader = sequences.file_reader(options.infile)
    try:
        f = os.popen('samtools view -hbS - > ' + options.outfile, 'w')
    except IOError:
        print("Error opening for writing BAM file '" + options.outfile + "'", file=sys.stderr)
        sys.exit(1)

    print(''.join(header_lines), file=f)

    for seq in seq_reader:
        end_range = len(seq)
        if len(seq) < options.read_length:
            end_range = 1
        for i in range(0, end_range, options.read_step):
            if len(seq) <= options.read_length:
                start = 0
                end = len(seq) - 1
            else:
                start = i
                end = start + options.read_length - 1

                if end > len(seq) - 1:
                    end  = len(seq) - 1
                    start = end - options.read_length + 1

            read = sequences.Fastq(options.read_prefix + ':' + seq.id + ':' + str(start + 1) + ':' + str(end + 1), seq[start:end+1], options.qual_char * (end - start + 1))

            print ('\t'.join([read.id,
                             '0',
                             seq.id,
                             str(start + 1),
                             '60',
                             str(len(read)) + 'M',
                             '*',
                             '*',
                             '*',
                             read.seq,
                             read.qual,
                             'RG:Z:' + options.read_group]), file=f)

            if end == len(seq) - 1:
                break

    f.close()
