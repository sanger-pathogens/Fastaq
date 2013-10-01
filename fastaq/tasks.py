import re
import copy
import random
from fastaq import sequences, utils

class Error (Exception): pass

def capillary_to_pairs(infile, outprefix):
    # hash the sequences, only taking longest where an end has been sequenced more than once
    seq_reader = sequences.file_reader(infile)
    fwd_seqs = {}
    rev_seqs = {}
    unpaired_seqs = {}

    for seq in seq_reader:
        id_info = seq.split_capillary_id()
        if id_info['dir'] == 'fwd':
            seq.id = id_info['prefix'] + '/1'
            h = fwd_seqs
        elif id_info['dir'] == 'rev':
            seq.id = id_info['prefix'] + '/2'
            h = rev_seqs
        else:
            seq.id = id_info['prefix']
            h = unpaired_seqs

        key = id_info['prefix']

        if key not in h or len(h[key]) < len(seq):
            h[key] = copy.copy(seq)

    # write the output files
    f_pe = utils.open_file_write(outprefix + '.paired.gz')
    f_up = utils.open_file_write(outprefix + '.unpaired.gz')

    for id in fwd_seqs:
        if id in rev_seqs:
            print(fwd_seqs[id], file=f_pe)
            print(rev_seqs[id], file=f_pe)
            del rev_seqs[id]
        else:
            print(fwd_seqs[id], file=f_up)

    for seq in rev_seqs.values():
        print(seq, file=f_up)

    for seq in unpaired_seqs.values():
        print(seq, file=f_up)

    utils.close(f_pe)
    utils.close(f_up)


def count_sequences(infile):
    '''Returns the number of sequences in a file'''
    seq_reader = sequences.file_reader(infile)
    n = 0
    for seq in seq_reader:
        n += 1
    return n


def deinterleave(infile, outfile_1, outfile_2, fasta_out=False):
    seq_reader = sequences.file_reader(infile)
    f_1 = utils.open_file_write(outfile_1)
    f_2 = utils.open_file_write(outfile_2)
    for seq in seq_reader:
        if fasta_out:
            print(sequences.Fasta(seq.id, seq.seq), file=f_1)
        else:
            print(seq, file=f_1)
        try:
            next(seq_reader)
        except StopIteration:
            utils.close(f_1)
            utils.close(f_2)
            raise Error('Error getting mate for sequence. Cannot continue')
        if fasta_out:
            print(sequences.Fasta(seq.id, seq.seq), file=f_2)
        else:
            print(seq, file=f_2)

    utils.close(f_1)
    utils.close(f_2)


def enumerate_names(infile, outfile, start_index=1, keep_illumina_suffix=False, rename_file=None):
    seq_reader = sequences.file_reader(infile)
    fout_seqs = utils.open_file_write(outfile)
    counter = start_index

    if keep_illumina_suffix:
        sequence_suffixes = ['/1', '/2']
    else:
        sequence_suffixes = []


    if rename_file is not None:
        fout_rename = utils.open_file_write(rename_file)
        print('#old\tnew', file=fout_rename)

    for seq in seq_reader:
        old_id = seq.id
        seq.id = str(counter)

        for suff in sequence_suffixes:
            if old_id.endswith(suff):
                seq.id += suff
                break

        if rename_file is not None:
            print(old_id, seq.id, sep='\t', file=fout_rename)

        print(seq, file=fout_seqs)
        counter += 1

    utils.close(fout_seqs)

    if rename_file is not None:
        utils.close(fout_rename)


def extend_gaps(infile, outfile, trim):
    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)

    for seq in seq_reader:
        if len(seq) < 2 * trim:
            continue

        gaps = seq.gaps()
        bases = list(seq.seq)

        # extend the length of each gap
        for gap in gaps:
            left_start = max(gap.start - trim, 0)
            right_end = min(gap.end + trim + 1, len(seq))

            for i in range(left_start, gap.start):
                bases[i] = 'N'

            for i in range(gap.end, right_end):
                bases[i] = 'N'

        seq.seq = ''.join(bases)

        # trim start/end bases and tidy up any resulting Ns at either end of the trimmed seq
        seq.trim(trim, trim)
        seq.trim_Ns()

        # check that there is some non-N sequence left over
        regex = re.compile('[^nN]')
        if regex.search(seq.seq) is not None:
            print(seq, file=fout)

    utils.close(fout)


def fasta_to_fastq(fasta_in, qual_in, outfile):
    fa_reader = sequences.file_reader(fasta_in)
    qual_reader = sequences.file_reader(qual_in, read_quals=True)
    f_out = utils.open_file_write(outfile)

    for seq in fa_reader:
        qual = next(qual_reader)
        if seq.id != qual.id:
            utils.close(f_out)
            raise Error('Mismatch in names from fasta and qual file', seq.id, qual.id)

        qual.seq = [int(x) for x in qual.seq.split()]
        print(seq.to_Fastq(qual.seq), file=f_out)

    utils.close(f_out)


def fastaq_to_mira_xml(infile, outfile):
    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)
    print('<?xml version="1.0"?>', '<trace_volume>', sep='\n', file=fout)

    for seq in seq_reader:
        print('    <trace>',
              '        <trace_name>' + seq.id + '</trace_name>',
              '        <clip_quality_right>' + str(len(seq)) + '</clip_quality_right>',
              '        <clip_vector_left>1</clip_vector_left>',
              '    </trace>', sep='\n', file=fout)


    print('</trace_volume>', file=fout)
    utils.close(fout)


def file_to_dict(infile, d):
    seq_reader = sequences.file_reader(infile)
    for seq in seq_reader:
        d[seq.id] = copy.copy(seq)


def filter(infile, outfile, minlength=0, maxlength=float('inf'), regex=None):
    seq_reader = sequences.file_reader(infile)
    f_out = utils.open_file_write(outfile)
    if regex is not None:
        r = re.compile(regex)

    for seq in seq_reader:
        if minlength <= len(seq) <= maxlength and (regex is None or r.search(seq.id)):
            print(seq, file=f_out)
    utils.close(f_out)
    

def get_ids(infile, outfile):
    seq_reader = sequences.file_reader(infile)
    f_out = utils.open_file_write(outfile)
    for seq in seq_reader:
        print(seq.id, file=f_out)
    utils.close(f_out)
    

def get_seqs_flanking_gaps(infile, outfile, left, right):
    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)

    print('#id', 'gap_start', 'gap_end', 'left_bases', 'right_bases', sep='\t', file=fout)

    for seq in seq_reader:
        gaps = seq.gaps()

        for gap in gaps:
            left_start = max(gap.start - left, 0)
            right_end = min(gap.end + right + 1, len(seq))
            print(seq.id,
                  gap.start + 1,
                  gap.end + 1,
                  seq.seq[left_start:gap.start],
                  seq.seq[gap.end + 1:right_end],
                  sep='\t', file=fout)

    utils.close(fout)


def interleave(infile_1, infile_2, outfile):
    seq_reader_1 = sequences.file_reader(infile_1)
    seq_reader_2 = sequences.file_reader(infile_2)
    f_out = utils.open_file_write(outfile)

    for seq_1 in seq_reader_1:
        try:
            seq_2 = next(seq_reader_2)
        except:
            utils.close(f_out)
            raise Error('Error getting mate for sequence', seq_1.id, ' ... cannot continue')

        print(seq_1, file=f_out)
        print(seq_2, file=f_out)

    try:
        seq_2 = next(seq_reader_2)
    except:
        seq_2 = None

    if seq_2 is not None:
        utils.close(f_out)
        raise Error('Error getting mate for sequence', seq_2.id, ' ... cannot continue')

    utils.close(f_out)


def make_random_contigs(contigs, length, outfile, name_by_letters=False, prefix='', seed=None, first_number=1):
    '''Makes a multi fasta file of random sequences, all the same length'''
    random.seed(a=seed)
    fout = utils.open_file_write(outfile)
    letters = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    letters_index = 0

    for i in range(contigs):
        if name_by_letters:
            name = letters[letters_index]
            letters_index += 1
            if letters_index == len(letters):
                letters_index = 0
        else:
            name = str(i + first_number)

        fa = sequences.Fasta(prefix + name, ''.join([random.choice('ACGT') for x in range(length)]))
        print(fa, file=fout)

    utils.close(fout)


def reverse_complement(infile, outfile):
    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)

    for seq in seq_reader:
        seq.revcomp()
        print(seq, file=fout)

    utils.close(fout)


def search_for_seq(infile, outfile, search_string):
    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)

    for seq in seq_reader:
        hits = seq.search(search_string)
        for hit in hits:
            print(seq.id, hit[0]+1, hit[1], sep='\t', file=fout)

    utils.close(fout)


def translate(infile, outfile, frame=0):
    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)

    for seq in seq_reader:
        print(seq.translate(frame=frame), file=fout)

    utils.close(fout)
    

def trim(infile, outfile, start, end):
    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)

    for seq in seq_reader:
        seq.trim(start, end)
        if len(seq):
            print(seq, file=fout)

    utils.close(fout)


def trim_Ns_at_end(infile, outfile):
    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)

    for seq in seq_reader:
        seq.trim_Ns()
        if len(seq):
            print(seq, file=fout)

    utils.close(fout)


def lengths_from_fai(fai_file, d):
    f = utils.open_file_read(fai_file)
    for line in f:
        (id, length) = line.rstrip().split()[:2]
        d[id] = int(length)
    utils.close(f)


def split_by_base_count(infile, outfiles_prefix, max_bases, max_seqs=None):
    '''Splits a fasta/q file into separate files, file size determined by number of bases.

    Puts <= max_bases in each split file The exception is a single sequence >=max_bases
    is put in its own file.  This does not split sequences.
    '''
    seq_reader = sequences.file_reader(infile)
    base_count = 0
    file_count = 1
    seq_count = 0
    fout = None
    if max_seqs is None:
        max_seqs = float('inf')

    for seq in seq_reader:
        if base_count == 0:
            fout = utils.open_file_write(outfiles_prefix + '.' + str(file_count))
            file_count += 1

        if base_count + len(seq) > max_bases or seq_count >= max_seqs:
            if base_count == 0:
                print(seq, file=fout)
                utils.close(fout)
            else:
                utils.close(fout)
                fout = utils.open_file_write(outfiles_prefix + '.' + str(file_count))
                print(seq, file=fout)
                base_count = len(seq)
                file_count += 1
                seq_count = 1
        else:
            base_count += len(seq)
            seq_count += 1
            print(seq, file=fout)

    utils.close(fout)


def split_by_fixed_size(infile, outfiles_prefix, chunk_size, tolerance, skip_if_all_Ns=False):
    '''Splits  fasta/q file into separate files, with up to (chunk_size + tolerance) bases in each file'''
    file_count = 1
    coords = []
    small_sequences = []  # sequences shorter than chunk_size
    seq_reader = sequences.file_reader(infile)
    f_coords = utils.open_file_write(outfiles_prefix + '.coords')

    for seq in seq_reader:
        if skip_if_all_Ns and seq.is_all_Ns():
             continue
        if len(seq) < chunk_size:
            small_sequences.append(copy.copy(seq))
        elif len(seq) <= chunk_size + tolerance:
            f = utils.open_file_write(outfiles_prefix + '.' + str(file_count))
            print(seq, file=f)
            utils.close(f)
            file_count += 1
        else:
            # make list of chunk coords
            chunks = [(x,x+chunk_size) for x in range(0, len(seq), chunk_size)]
            if chunks[-1][1] - 1 > len(seq):
                chunks[-1] = (chunks[-1][0], len(seq))
            if len(chunks) > 1 and (chunks[-1][1] - chunks[-1][0]) <= tolerance:
                chunks[-2] = (chunks[-2][0], chunks[-1][1])
                chunks.pop()

            # write one output file per chunk
            offset = 0
            for chunk in chunks:
                if not(skip_if_all_Ns and seq.is_all_Ns(start=chunk[0], end=chunk[1]-1)):
                    f = utils.open_file_write(outfiles_prefix + '.' + str(file_count))
                    chunk_id = seq.id + ':' + str(chunk[0]+1) + '-' + str(chunk[1])
                    print(sequences.Fasta(chunk_id, seq[chunk[0]:chunk[1]]), file=f)
                    print(chunk_id, seq.id, offset, sep='\t', file=f_coords)
                    utils.close(f)
                    file_count += 1

                offset += chunk[1] - chunk[0]

    # write files of small sequences
    if len(small_sequences):
        f = utils.open_file_write(outfiles_prefix + '.' + str(file_count))
        file_count += 1
        base_count = 0
        for seq in small_sequences:
            if base_count > 0 and base_count + len(seq) > chunk_size + tolerance:
                utils.close(f)
                f = utils.open_file_write(outfiles_prefix + '.' + str(file_count))
                file_count += 1
                base_count = 0
              
            print(seq, file=f)
            base_count += len(seq)

        utils.close(f)


def replace_bases(infile, outfile, old, new):
    seq_reader = sequences.file_reader(infile)
    f_out = utils.open_file_write(outfile)

    for seq in seq_reader:
        seq.replace_bases(old, new)
        print(seq, file=f_out)

    utils.close(f_out)


def strip_illumina_suffix(infile, outfile):
    seq_reader = sequences.file_reader(infile)
    f_out = utils.open_file_write(outfile)

    for seq in seq_reader:
        seq.strip_illumina_suffix()
        print(seq, file=f_out)

    utils.close(f_out)


def to_fasta(infile, outfile, line_length=60, strip_after_first_whitespace=False):
    seq_reader = sequences.file_reader(infile)
    f_out = utils.open_file_write(outfile)
    original_line_length = sequences.Fasta.line_length
    sequences.Fasta.line_length = line_length

    for seq in seq_reader:
        if strip_after_first_whitespace:
            seq.strip_after_first_whitespace()

        if type(seq) == sequences.Fastq:
            print(sequences.Fasta(seq.id, seq.seq), file=f_out)
        else:
            print(seq, file=f_out)

    utils.close(f_out)
    sequences.Fasta.line_length = original_line_length


def to_quasr_primers(infile, outfile):
    seq_reader = sequences.file_reader(infile)
    f_out = utils.open_file_write(outfile)

    for seq in seq_reader:
        seq2 = copy.copy(seq)
        seq2.revcomp()
        print(seq.seq, seq2.seq, sep='\t', file=f_out)

    utils.close(f_out)

    
def to_unique_by_id(infile, outfile):
    seq_reader = sequences.file_reader(infile)
    seqs = {}
    ids_in_order = []

    # has the reads, keeping the longest one when we get the same
    # name more than once
    for seq in seq_reader:
        if len(seq) == 0:
           continue
        if seq.id not in seqs:
            seqs[seq.id] = copy.copy(seq)
            ids_in_order.append(seq.id)
        elif len(seqs[seq.id]) < len(seq):
            seqs[seq.id] = copy.copy(seq)

    # write the output
    f_out = utils.open_file_write(outfile)
    for id in ids_in_order:
        print(seqs[id], file=f_out)
    utils.close(f_out)
