import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = 'Simulates long reads from a sequence file. Can optionally make insertions into the reads, like pacbio does. If insertions made, coverage calculation is done before the insertions (so total read length may appear longer then expected).',
        usage = 'fastaq long_read_simulate [options] <infile> <outfile>')

    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output FASTA file')

    parser.add_argument('--method', help='How to sample the read positions and lengths. Choose from 1) "tiling", where reads of fixed length are taken at equal intervals from the reference. 2) "unfiform", where reads of fixed length taken at positions sampled uniformly. 3) "gamma", where reads lengths are taken from a gamma distribution, and positions sampled uniformly. [%(default)s]', default='tiling', choices=['tiling', 'uniform', 'gamma'], metavar='tiling|uniform|gamma')
    parser.add_argument('--seed', type=int, help='Seed for random number generator [default: use python\'s default]', metavar='INT')
    parser.add_argument('--qual', help='Write a file of fake quality scores called outfile.qual, all bases same quality [%(default)s]', metavar='INT')
    parser.add_argument('--fixed_read_length', type=int, help='Length of each read. Only applies if method is tile or uniform. [%(default)s]', default=20000, metavar='INT')
    parser.add_argument('--coverage', type=float, help='Read coverage. Only applies if method is gamma or uniform. [%(default)s]', default=2, metavar='FLOAT')


    tiling_group = parser.add_argument_group('tiling options')
    tiling_group.add_argument('--tile_step', type=int, help='Distance between start of each read [%(default)s]', default=10000, metavar='INT')

    gamma_group = parser.add_argument_group('gamma options')
    gamma_group.add_argument('--gamma_shape', type=float, help='Shape parameter of gamma distribution [%(default)s]', default=1.2, metavar='FLOAT')
    gamma_group.add_argument('--gamma_scale', type=float, help='Scale parameter of gamma distribution [%(default)s]', default=6000, metavar='FLOAT')
    gamma_group.add_argument('--gamma_min_length', type=int, help='Minimum read length [%(default)s]', default=20000, metavar='INT')

    ins_group = parser.add_argument_group('options to add insertions to reads')
    ins_group.add_argument('--ins_skip', type=int, help='Insert a random base every --skip bases plus or minus --ins_window. If this option is used, must also use --ins_window.', metavar='INT')
    ins_group.add_argument('--ins_window', type=int, help='See --ins_skip. If this option is used, must also use --ins_skip.', metavar='INT')


    options = parser.parse_args()
    tasks.make_long_reads(
        options.infile,
        options.outfile,
        method=options.method,
        fixed_read_length=options.fixed_read_length,
        coverage=options.coverage,
        tile_step=options.tile_step,
        gamma_shape=options.gamma_shape,
        gamma_scale=options.gamma_scale,
        gamma_min_length=options.gamma_min_length,
        seed=options.seed,
        ins_skip=options.ins_skip,
        ins_window=options.ins_window
    )

    if options.qual:
        tasks.fastaq_to_fake_qual(options.outfile, options.outfile + '.qual', q=options.qual)
