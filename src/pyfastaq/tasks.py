import re
import sys
import copy
import random
from pyfastaq import sequences, utils, caf

class Error (Exception): pass


class IncompatibleParametersError(Exception):
    pass


def acgtn_only(infile, outfile):
    '''Replace every non-acgtn (case insensitve) character with an N'''
    f = utils.open_file_write(outfile)
    for seq in sequences.file_reader(infile):
        seq.replace_non_acgt()
        print(seq, file=f)
    utils.close(f)


def caf_to_fastq(infile, outfile, min_length=0, trim=False):
    '''Convert a CAF file to fastq. Reads shorter than min_length are not output. If clipping information is in the CAF file (with a line Clipping QUAL ...) and trim=True, then trim the reads'''
    caf_reader = caf.file_reader(infile)
    fout = utils.open_file_write(outfile)

    for c in caf_reader:
        if trim:
            if c.clip_start is not None and c.clip_end is not None:
                c.seq.seq = c.seq.seq[c.clip_start:c.clip_end + 1]
                c.seq.qual = c.seq.qual[c.clip_start:c.clip_end + 1]
            else:
                print('Warning: no clipping info for sequence', c.id, file=sys.stderr)


        if len(c.seq) >= min_length:
            print(c.seq, file=fout)

    utils.close(fout)


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


def enumerate_names(infile, outfile, start_index=1, keep_illumina_suffix=False, rename_file=None, suffix=None):
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

        if suffix is not None:
            seq.id += suffix

        print(seq, file=fout_seqs)
        counter += 1

    utils.close(fout_seqs)

    if rename_file is not None:
        utils.close(fout_rename)


def expand_nucleotides(infile, outfile):
    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)

    for seq in seq_reader:
        seqs = seq.expand_nucleotides()
        if len(seqs) > 1:
            for s in seqs:
                print(s, file=fout)
        else:
            print(seq, file=fout)


def trim_contigs(infile, outfile, trim):
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


def fastaq_to_fake_qual(infile, outfile, q=40):
    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)

    for seq in seq_reader:
        print('>' + seq.id, file=fout)
        if sequences.Fasta.line_length == 0:
            print(' '.join([str(q)] * len(seq)), file=fout)
        else:
            for i in range(0, len(seq), sequences.Fasta.line_length):
                print(' '.join([str(q)] * min(sequences.Fasta.line_length, len(seq) - i)), file=fout)

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


def fastaq_to_orfs_gff(infile, outfile, min_length=300, tool_name='fastaq'):
    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)
    for seq in seq_reader:
        orfs = seq.all_orfs(min_length=min_length)
        for coords, revcomp in orfs:
            if revcomp:
                strand = '-'
            else:
                strand = '+'

            print(seq.id, tool_name, 'CDS', coords.start+1, coords.end+1, '.', strand, '.', sep='\t', file=fout)

    utils.close(fout)


def file_to_dict(infile, d):
    seq_reader = sequences.file_reader(infile)
    for seq in seq_reader:
        d[seq.id] = copy.copy(seq)


def filter(
      infile,
      outfile,
      minlength=0,
      maxlength=float('inf'),
      regex=None,
      ids_file=None,
      invert=False,
      mate_in=None,
      mate_out=None,
      both_mates_pass=True,
      check_comments=False
    ):
    if check_comments and not regex:
        raise IncompatibleParametersError(
            "--check_comments can only be passed with --regex"
        )

    ids_from_file = set()
    if ids_file is not None:
        f = utils.open_file_read(ids_file)
        for line in f:
            ids_from_file.add(line.rstrip())
        utils.close(f)

    if mate_in:
        if mate_out is None:
            raise Error('Error in filter! mate_in provided. Must also provide mate_out')

        seq_reader_mate = sequences.file_reader(mate_in)
        f_out_mate = utils.open_file_write(mate_out)

    seq_reader = sequences.file_reader(infile)
    f_out = utils.open_file_write(outfile)
    if regex is not None:
        r = re.compile(regex)


    def passes(seq, name_regex):
        # remove trailing comments from FASTQ readname lines
        matches = name_regex.match(seq.id)
        if matches is not None and not check_comments:
            clean_seq_id = matches.group(1)
        else:
            clean_seq_id = seq.id
        
        return minlength <= len(seq) <= maxlength \
              and (regex is None or r.search(clean_seq_id) is not None) \
              and (ids_file is None or clean_seq_id in ids_from_file)
        
    name_regex = re.compile(r'^([^\s]+).*?$')
	
    for seq in seq_reader:
        seq_passes = passes(seq, name_regex)
        if mate_in:
            try:
                seq_mate = next(seq_reader_mate)
            except:
                utils.close(f_out)
                raise Error('Error getting mate for sequence', seq.id, ' ... cannot continue')

            mate_passes = passes(seq_mate, name_regex)
            want_the_pair = (seq_passes and mate_passes) \
                            or (( seq_passes or mate_passes) and not both_mates_pass)
            if want_the_pair != invert:
                print(seq, file=f_out)
                print(seq_mate, file=f_out_mate)
        elif seq_passes != invert:
            print(seq, file=f_out)
    utils.close(f_out)
    if mate_in:
        utils.close(f_out_mate)


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


def interleave(infile_1, infile_2, outfile, suffix1=None, suffix2=None):
    '''Makes interleaved file from two sequence files. If used, will append suffix1 onto end
    of every sequence name in infile_1, unless it already ends with suffix1. Similar for sufffix2.'''
    seq_reader_1 = sequences.file_reader(infile_1)
    seq_reader_2 = sequences.file_reader(infile_2)
    f_out = utils.open_file_write(outfile)

    for seq_1 in seq_reader_1:
        try:
            seq_2 = next(seq_reader_2)
        except:
            utils.close(f_out)
            raise Error('Error getting mate for sequence', seq_1.id, ' ... cannot continue')

        if suffix1 is not None and not seq_1.id.endswith(suffix1):
            seq_1.id += suffix1
        if suffix2 is not None and not seq_2.id.endswith(suffix2):
            seq_2.id += suffix2

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


def mean_length(infile, limit=None):
    '''Returns the mean length of the sequences in the input file. By default uses all sequences. To limit to the first N sequences, use limit=N'''
    total = 0
    count = 0
    seq_reader = sequences.file_reader(infile)
    for seq in seq_reader:
        total += len(seq)
        count += 1
        if limit is not None and count >= limit:
            break

    assert count > 0
    return total / count


def merge_to_one_seq(infile, outfile, seqname='union'):
    '''Takes a multi fasta or fastq file and writes a new file that contains just one sequence, with the original sequences catted together, preserving their order'''
    seq_reader = sequences.file_reader(infile)
    seqs = []

    for seq in seq_reader:
        seqs.append(copy.copy(seq))

    new_seq = ''.join([seq.seq for seq in seqs])

    if type(seqs[0]) == sequences.Fastq:
        new_qual = ''.join([seq.qual for seq in seqs])
        seqs[:] = []
        merged = sequences.Fastq(seqname, new_seq, new_qual)
    else:
        merged = sequences.Fasta(seqname, new_seq)
        seqs[:] = []

    f = utils.open_file_write(outfile)
    print(merged, file=f)
    utils.close(f)


def reverse_complement(infile, outfile):
    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)

    for seq in seq_reader:
        seq.revcomp()
        print(seq, file=fout)

    utils.close(fout)


def scaffolds_to_contigs(infile, outfile, number_contigs=False):
    '''Makes a file of contigs from scaffolds by splitting at every N.
       Use number_contigs=True to add .1, .2, etc onto end of each
       contig, instead of default to append coordinates.'''
    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)

    for seq in seq_reader:
        contigs = seq.contig_coords()
        counter = 1
        for contig in contigs:
            if number_contigs:
                name = seq.id + '.' + str(counter)
                counter += 1
            else:
                name = '.'.join([seq.id, str(contig.start + 1), str(contig.end + 1)])
            print(sequences.Fasta(name, seq[contig.start:contig.end+1]), file=fout)

    utils.close(fout)


def search_for_seq(infile, outfile, search_string):
    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)

    for seq in seq_reader:
        hits = seq.search(search_string)
        for hit in hits:
            print(seq.id, hit[0]+1, hit[1], sep='\t', file=fout)

    utils.close(fout)


def sequence_trim(infile_1, infile_2, outfile_1, outfile_2, to_trim_file, min_length=50, check_revcomp=False):
    to_trim_seqs = {}
    file_to_dict(to_trim_file, to_trim_seqs)
    trim_seqs = [x.seq for x in to_trim_seqs.values()]
    if check_revcomp:
        for seq in to_trim_seqs.values():
            seq.revcomp()
        trim_seqs_revcomp = [x.seq for x in to_trim_seqs.values()]
    else:
        trim_seqs_revcomp = []

    seq_reader_1 = sequences.file_reader(infile_1)
    seq_reader_2 = sequences.file_reader(infile_2)
    f_out_1 = utils.open_file_write(outfile_1)
    f_out_2 = utils.open_file_write(outfile_2)

    for seq_1 in seq_reader_1:
        try:
            seq_2 = next(seq_reader_2)
        except:
            utils.close(f_out)
            raise Error('Error getting mate for sequence', seq_1.id, ' ... cannot continue')

        for seq in seq_1, seq_2:
            for trim_seq in trim_seqs:
                if seq.seq.startswith(trim_seq):
                    seq.trim(len(trim_seq),0)
                    break

            for trim_seq in trim_seqs_revcomp:
                if seq.seq.endswith(trim_seq):
                    seq.trim(0,len(trim_seq))
                    break

        if len(seq_1) >= min_length and len(seq_2) >= min_length:
            print(seq_1, file=f_out_1)
            print(seq_2, file=f_out_2)


    utils.close(f_out_1)
    utils.close(f_out_2)


def sort_by_size(infile, outfile, smallest_first=False):
    '''Sorts input sequence file by biggest sequence first, writes sorted output file. Set smallest_first=True to have smallest first'''
    seqs = {}
    file_to_dict(infile, seqs)
    seqs = list(seqs.values())
    seqs.sort(key=lambda x: len(x), reverse=not smallest_first)
    fout = utils.open_file_write(outfile)
    for seq in seqs:
        print(seq, file=fout)
    utils.close(fout)


def sort_by_name(infile, outfile):
    '''Sorts input sequence file by sort -d -k1,1, writes sorted output file.'''
    seqs = {}
    file_to_dict(infile, seqs)
    #seqs = list(seqs.values())
    #seqs.sort()
    fout = utils.open_file_write(outfile)
    for name in sorted(seqs):
        print(seqs[name], file=fout)
    utils.close(fout)


def to_fastg(infile, outfile, circular=None):
    '''Writes a FASTG file in SPAdes format from input file. Currently only whether or not a sequence is circular is supported. Put circular=set of ids, or circular=filename to make those sequences circular in the output. Puts coverage=1 on all contigs'''
    if circular is None:
        to_circularise = set()
    elif type(circular) is not set:
        f = utils.open_file_read(circular)
        to_circularise = set([x.rstrip() for x in f.readlines()])
        utils.close(f)
    else:
        to_circularise = circular

    seq_reader = sequences.file_reader(infile)
    fout = utils.open_file_write(outfile)
    nodes = 1

    for seq in seq_reader:
        new_id = '_'.join([
            'NODE', str(nodes),
            'length', str(len(seq)),
            'cov', '1',
            'ID', seq.id
        ])

        if seq.id in to_circularise:
            seq.id = new_id + ':' + new_id + ';'
            print(seq, file=fout)
            seq.revcomp()
            seq.id = new_id + "':" + new_id + "';"
            print(seq, file=fout)
        else:
            seq.id = new_id + ';'
            print(seq, file=fout)
            seq.revcomp()
            seq.id = new_id + "';"
            print(seq, file=fout)

        nodes += 1

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


def length_offsets_from_fai(fai_file):
    '''Returns a dictionary of positions of the start of each sequence, as
       if all the sequences were catted into one sequence.
       eg if file has three sequences, seq1 10bp, seq2 30bp, seq3 20bp, then
       the output would be: {'seq1': 0, 'seq2': 10, 'seq3': 40}'''
    positions = {}
    total_length = 0
    f = utils.open_file_read(fai_file)

    for line in f:
        try:
            (name, length) = line.rstrip().split()[:2]
            length = int(length)
        except:
            raise Error('Error reading the following line of fai file ' + fai_file + '\n' + line)

        positions[name] = total_length
        total_length += length

    utils.close(f)
    return positions


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


def split_by_fixed_size_onefile(infile, outfile, chunk_size, tolerance, skip_if_all_Ns=False):
    '''Splits each sequence in infile into chunks of fixed size, last chunk can be up to
       (chunk_size + tolerance) in length'''
    seq_reader = sequences.file_reader(infile)
    f_out = utils.open_file_write(outfile)
    for seq in seq_reader:
        for i in range(0, len(seq), chunk_size):
            if i + chunk_size + tolerance >= len(seq):
                end = len(seq)
            else:
                end = i + chunk_size

            subseq = seq.subseq(i, end)
            if not (skip_if_all_Ns and subseq.is_all_Ns()):
                subseq.id += '.' + str(i+1) + '_' + str(end)
                print(subseq, file=f_out)

            if end == len(seq):
                break

    utils.close(f_out)


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


def stats_from_fai(infile):
    '''Returns dictionary of length stats from an fai file. Keys are: longest, shortest, mean, total_length, N50, number'''
    f = utils.open_file_read(infile)
    try:
        lengths = sorted([int(line.split('\t')[1]) for line in f], reverse=True)
    except:
        raise Error('Error getting lengths from fai file ' + infile)
    utils.close(f)

    stats = {}
    if len(lengths) > 0:
        stats['longest'] = max(lengths)
        stats['shortest'] = min(lengths)
        stats['total_length'] = sum(lengths)
        stats['mean'] = stats['total_length'] / len(lengths)
        stats['number'] = len(lengths)

        cumulative_length = 0
        for length in lengths:
            cumulative_length += length
            if cumulative_length >= 0.5 * stats['total_length']:
                stats['N50'] = length
                break
    else:
        stats = {x: 0 for x in ('longest', 'shortest', 'mean', 'N50', 'total_length', 'number')}

    return stats


def to_boulderio(infile, outfile):
    '''Converts input sequence file into a "Boulder-IO format", as used by primer3'''
    seq_reader = sequences.file_reader(infile)
    f_out = utils.open_file_write(outfile)

    for sequence in seq_reader:
        print("SEQUENCE_ID=" + sequence.id, file=f_out)
        print("SEQUENCE_TEMPLATE=" + sequence.seq, file=f_out)
        print("=", file=f_out)

    utils.close(f_out)


def to_fasta(infile, outfile, line_length=60, strip_after_first_whitespace=False, check_unique=False):
    seq_reader = sequences.file_reader(infile)
    f_out = utils.open_file_write(outfile)
    original_line_length = sequences.Fasta.line_length
    sequences.Fasta.line_length = line_length
    if check_unique:
        used_names = {}

    for seq in seq_reader:
        if strip_after_first_whitespace:
            seq.strip_after_first_whitespace()

        if check_unique:
            used_names[seq.id] = used_names.get(seq.id, 0) + 1

        if type(seq) == sequences.Fastq:
            print(sequences.Fasta(seq.id, seq.seq), file=f_out)
        else:
            print(seq, file=f_out)

    utils.close(f_out)
    sequences.Fasta.line_length = original_line_length

    if check_unique:
        all_unique = True

        for name, count in used_names.items():
            if count > 1:
                print('Sequence name "' + name + '" not unique. Found', count, 'times', file=sys.stderr)
                all_unique = False

        if not all_unique:
            raise Error('Not all sequence names unique. Cannot continue')



def to_fasta_union(infile, outfile, seqname='union'):
    seq_reader = sequences.file_reader(infile)
    new_seq = []

    for seq in seq_reader:
        new_seq.append(seq.seq)

    f_out = utils.open_file_write(outfile)
    print(sequences.Fasta(seqname, ''.join(new_seq)), file=f_out)
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
