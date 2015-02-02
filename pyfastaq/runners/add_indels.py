import argparse
import sys
import random
from pyfastaq import sequences, utils, intervals

def run(description):
    parser = argparse.ArgumentParser(
        description = description,
        usage = 'fastaq add_indels [options] <infile> <outfile>')
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output file')
    parser.add_argument('-d','--delete', action='append', help='Delete the given bases from the given sequence. Format same as samtools view: name:start-end. This option can be used multiple times (once for each region to delete). Overlapping coords will be merged before deleting', metavar='Name:start:bases')
    parser.add_argument('--delete_range', help='Deletes bases starting at position P in each sequence of the input file. Deletes start + (n-1)*step bases from sequence n.', metavar='P,start,step')
    parser.add_argument('-i','--insert', action='append', help='Insert a random string of bases at the given position. Format is name:position:number_to_add. Bases are added after the position. This option can be used multiple times', metavar='Name:start:bases')
    parser.add_argument('--insert_range', help='Inserts random bases starting after position P in each sequence of the input file. Inserts start + (n-1)*step bases into sequence n.', metavar='P,start,step')
    options = parser.parse_args()

    test_ops = [int(x is not None) for x in [options.delete, options.insert, options.delete_range, options.insert_range]]

    if sum(test_ops) != 1:
        print('Must use one of --delete, --insert, --delete_range, --insert_range. Cannot continue', file=sys.stderr)
        sys.exit(1)


    def range2dic(range_in):
        if range_in is None:
            return {}
        (pos, start, step) = range_in.split(',')
        d = {}
        d['pos'] = int(pos) - 1
        d['bases'] = int(start)
        d['step'] = int(step)
        return d

    delete_range = range2dic(options.delete_range)
    insert_range = range2dic(options.insert_range)


    # convert the -d regions into sequence name, start and end coords
    to_delete = {}
    if options.delete:
        for s in options.delete:
            id, coords = s.rsplit(':')
            start, end = [int(x)-1 for x in coords.split('-')]
            if id not in to_delete:
                to_delete[id] = []
            to_delete[id].append(intervals.Interval(start, end))


    to_insert = {}
    if options.insert:
        for s in options.insert:
            id, pos, bases = s.rsplit(':',2)
            pos = int(pos) - 1
            bases = int(bases)
            if id not in to_insert:
                to_insert[id] = []
            to_insert[id].append((pos, bases))


    assert len(to_delete) * len(to_insert) == 0

    # merge overlapping regions to be deleted
    for l in to_delete.values():
        intervals.merge_overlapping_in_list(l)

    # sort positions to be inserted
    for l in to_insert.values():
        l.sort()

    # read in the fasta/q file and print outfile with deleted sequences
    seq_reader = sequences.file_reader(options.infile)
    f = utils.open_file_write(options.outfile)

    for seq in seq_reader:
        if seq.id in to_delete:
            # delete regions for this sequence, but start at the end so the
            # coords don't get messed up after the first deletion
            for inter in reversed(to_delete[seq.id]):
                seq.seq = seq.seq[:inter.start] + seq.seq[inter.end + 1:]
        elif options.delete_range:
            seq.seq = seq.seq[:delete_range['pos']] + seq.seq[delete_range['pos'] + delete_range['bases']:]
            delete_range['bases'] += delete_range['step']
        elif seq.id in to_insert:
            for pos, bases in reversed(to_insert[seq.id]):
                seq.seq = seq.seq[:pos + 1] + ''.join([random.choice('ACGT') for x in range(bases)]) + seq.seq[pos + 1:]
        elif options.insert_range:
            seq.seq = seq.seq[:insert_range['pos'] + 1] + ''.join([random.choice('ACGT') for x in range(insert_range['bases'])]) +  seq.seq[insert_range['pos'] + 1:]
            insert_range['bases'] += insert_range['step']

        print(seq, file=f)

    utils.close(f)
