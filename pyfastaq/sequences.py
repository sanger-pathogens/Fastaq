import copy
import re
import string
import random
import itertools

from pyfastaq import utils, intervals, genetic_codes

class Error (Exception): pass


# python 3's seek is glacially slow. When we read a fasta file, we know
# we've reached the end of a sequence when we get a new line starting with
# '>'. Instead of using seek and tell, we just remember the previous line
# of the file, for any given filehandle
previous_lines = {}

genetic_code = 1

redundant_nts = {
    'R': ('A', 'G'),
    'Y': ('C', 'T'),
    'S': ('C', 'G'),
    'W': ('A', 'T'),
    'K': ('G', 'T'),
    'M': ('A', 'C'),
    'B': ('C', 'G', 'T'),
    'D': ('A', 'G', 'T'),
    'H': ('A', 'C', 'T'),
    'V': ('A', 'C', 'G'),
    'N': ('A', 'C', 'G', 'T')
}

def file_reader(fname, read_quals=False):
    '''Iterates over a FASTA or FASTQ file, yielding the next sequence in the file until there are no more sequences'''
    f = utils.open_file_read(fname)
    line = f.readline()
    phylip_regex = re.compile('^\s*[0-9]+\s+[0-9]+$')
    gbk_regex = re.compile('^LOCUS\s+\S')

    if line.startswith('>'):
        seq = Fasta()
        previous_lines[f] = line
    elif line.startswith('##gff-version 3'):
        seq = Fasta()
        # if a GFF file, need to skip past all the annotation
        # and get to the fasta sequences at the end of the file
        while not line.startswith('>'):
            line = f.readline()
            if not line:
                utils.close(f)
                raise Error('No sequences found in GFF file "' + fname + '"')

        seq = Fasta()
        previous_lines[f] = line
    elif line.startswith('ID   ') and line[5] != ' ':
        seq = Embl()
        previous_lines[f] = line
    elif gbk_regex.search(line):
        seq = Embl()
        previous_lines[f] = line
    elif line.startswith('@'):
        seq = Fastq()
        previous_lines[f] = line
    elif phylip_regex.search(line):
        # phylip format could be interleaved or not, need to look at next
        # couple of lines to figure that out. Don't expect these files to
        # be too huge, so just store all the sequences in memory
        number_of_seqs, bases_per_seq = line.strip().split()
        number_of_seqs = int(number_of_seqs)
        bases_per_seq = int(bases_per_seq)
        got_blank_line = False

        first_line = line
        seq_lines = []
        while 1:
            line = f.readline()
            if line == '':
                break
            elif line == '\n':
                got_blank_line = True
            else:
                seq_lines.append(line.rstrip())
        utils.close(f)

        if len(seq_lines) == 1 or len(seq_lines) == number_of_seqs:
            sequential = True
        elif seq_lines[0][10] != ' ' and seq_lines[1][10] == ' ':
            sequential = True
        else:
            sequential = False

        # if the 11th char of second sequence line is a space,  then the file is sequential, e.g.:
        # GAGCCCGGGC AATACAGGGT AT
        # as opposed to:
        # Salmo gairAAGCCTTGGC AGTGCAGGGT
        if sequential:
            current_id = None
            current_seq = ''
            for line in seq_lines:
                if len(current_seq) == bases_per_seq or len(current_seq) == 0:
                    if current_id is not None:
                        yield Fasta(current_id, current_seq.replace('-', ''))
                    current_seq = ''
                    current_id, new_bases = line[0:10].rstrip(), line.rstrip()[10:]
                else:
                    new_bases = line.rstrip()

                current_seq += new_bases.replace(' ','')

            yield Fasta(current_id, current_seq.replace('-', ''))
        else:
            # seaview files start all seqs at pos >=12. Other files start
            # their sequence at the start of the line
            if seq_lines[number_of_seqs + 1][0] == ' ':
                first_gap_pos = seq_lines[0].find(' ')
                end_of_gap = first_gap_pos
                while seq_lines[0][end_of_gap] == ' ':
                    end_of_gap += 1
                first_seq_base = end_of_gap
            else:
                first_seq_base = 10

            seqs = []
            for i in range(number_of_seqs):
                name, bases = seq_lines[i][0:first_seq_base].rstrip(), seq_lines[i][first_seq_base:]
                seqs.append(Fasta(name, bases))

            for i in range(number_of_seqs, len(seq_lines)):
                seqs[i%number_of_seqs].seq += seq_lines[i]

            for fa in seqs:
                fa.seq = fa.seq.replace(' ','').replace('-','')
                yield fa

        return
    elif line == '':
        utils.close(f)
        return
    else:
        utils.close(f)
        raise Error('Error determining file type from file "' + fname + '". First line is:\n' + line.rstrip())

    try:
        while seq.get_next_from_file(f, read_quals):
            yield seq
    finally:
        utils.close(f)


class Fasta:
    '''Class to store and manipulate FASTA sequences. They have two things: a name and a sequence'''
    # this defines the line length when printing sequences
    line_length = 60

    def _get_id_from_header_line(self, line):
        if line.startswith('>'):
            return line.rstrip()[1:]
        else:
            raise Error('Error! expected line starting with ">", but got this:\n', line)


    def __init__(self, id_in=None, seq_in=None):
        self.id = id_in
        self.seq = seq_in

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return len(self.seq)

    def subseq(self, start, end):
        '''Returns Fasta object with the same name, of the bases from start to end, but not including end'''
        return Fasta(self.id, self.seq[start:end])

    def split_capillary_id(self):
        '''Gets the prefix and suffix of an name of a capillary read, e.g. xxxxx.p1k or xxxx.q1k. Returns a tuple (prefix, suffx)'''
        try:
            a = self.id.rsplit('.', 1)
            if a[1].startswith('p'):
                dir = 'fwd'
            elif a[1].startswith('q'):
                dir = 'rev'
            else:
                dir = 'unk'

            return {'prefix': a[0], 'dir': dir, 'suffix':a[1]}
        except:
            raise Error('Error in split_capillary_id() on ID', self.id)

    def expand_nucleotides(self):
        '''Assumes sequence is nucleotides. Returns list of all combinations of redundant nucleotides. e.g. R is A or G, so CRT would have combinations CAT and CGT'''
        s = list(self.seq)
        for i in range(len(s)):
            if s[i] in redundant_nts:
                s[i] = ''.join(redundant_nts[s[i]])

        seqs = []
        for x in itertools.product(*s):
            seqs.append(Fasta(self.id + '.' + str(len(seqs) + 1), ''.join(x)))
        return seqs

    def strip_after_first_whitespace(self):
        '''Removes everything in the name after the first whitespace character'''
        self.id = self.id.split()[0]

    def strip_illumina_suffix(self):
        '''Removes any trailing /1 or /2 off the end of the name'''
        if self.id.endswith('/1') or self.id.endswith('/2'):
            self.id = self.id[:-2]

    def revcomp(self):
        '''Reverse complements the sequence'''
        self.seq = self.seq.translate(str.maketrans("ATCGatcg", "TAGCtagc"))[::-1]

    def is_all_Ns(self, start=0, end=None):
        '''Returns true if the sequence is all Ns (upper or lower case)'''
        if end is not None:
            if start > end:
                raise Error('Error in is_all_Ns. Start coord must be <= end coord')
            end += 1
        else:
            end = len(self)

        if len(self) == 0:
            return False
        else:
            return re.search('[^Nn]', self.seq[start:end]) is None

    def trim_Ns(self):
        '''Removes any leading or trailing N or n characters from the sequence'''
        self.seq = self.seq.strip('Nn')

    def add_insertions(self, skip=10, window=1, test=False):
        '''Adds a random base within window bases around every skip bases. e.g. skip=10, window=1 means a random base added somwhere in theintervals [9,11], [19,21] ... '''
        assert 2 * window < skip
        new_seq = list(self.seq)
        for i in range(len(self) - skip, 0, -skip):
            pos = random.randrange(i - window, i + window + 1)
            base = random.choice(['A', 'C', 'G', 'T'])
            if test:
                base = 'N'
            new_seq.insert(pos, base)

        self.seq = ''.join(new_seq)

    def replace_bases(self, old, new):
        '''Replaces all occurrences of 'old' with 'new' '''
        self.seq = self.seq.replace(old, new)


    def replace_non_acgt(self):
        '''Replace all non acgt characters with an N (case insensitive)'''
        self.seq = re.sub(r'''[^acgtACGTnN]''', 'N', self.seq)


    def replace_interval(self, start, end, new):
        '''Replaces the sequence from start to end with the sequence "new"'''
        if start > end or start > len(self) - 1 or end > len(self) - 1:
            raise Error('Error replacing bases ' + str(start) + '-' + str(end) + ' in sequence ' + self.id)

        self.seq = self.seq[0:start] + new + self.seq[end + 1:]

    def gaps(self, min_length = 1):
        '''Finds the positions of all gaps in the sequence that are at least min_length long. Returns a list of Intervals. Coords are zero-based'''
        gaps = []
        regex = re.compile('N+', re.IGNORECASE)
        for m in regex.finditer(self.seq):
             if m.span()[1] - m.span()[0] + 1 >= min_length:
                 gaps.append(intervals.Interval(m.span()[0], m.span()[1] - 1))
        return gaps

    def contig_coords(self):
        '''Finds coords of contigs, i.e. everything that's not a gap (N or n). Returns a list of Intervals. Coords are zero-based'''
        # contigs are the opposite of gaps, so work out the coords from the gap coords
        gaps = self.gaps()

        if len(gaps) == 0:
            return [intervals.Interval(0, len(self) - 1)]

        coords = [0]
        for g in gaps:
            if g.start == 0:
                coords = [g.end + 1]
            else:
                coords += [g.start - 1, g.end + 1]

        if coords[-1] < len(self):
            coords.append(len(self) - 1)

        return [intervals.Interval(coords[i], coords[i+1]) for i in range(0, len(coords)-1,2)]


    def orfs(self, frame=0, revcomp=False):
        '''Returns a list of ORFs that the sequence has, starting on the given
           frame. Each returned ORF is an interval.Interval object.
           If revomp=True, then finds the ORFs of the reverse complement
           of the sequence.'''
        assert frame in [0,1,2]
        if revcomp:
            self.revcomp()

        aa_seq = self.translate(frame=frame).seq.rstrip('X')
        if revcomp:
            self.revcomp()

        orfs = _orfs_from_aa_seq(aa_seq)
        for i in range(len(orfs)):
            if revcomp:
                start = len(self) - (orfs[i].end * 3 + 3) - frame
                end = len(self) - (orfs[i].start * 3) - 1 - frame
            else:
                start = orfs[i].start * 3 + frame
                end = orfs[i].end * 3 + 2 + frame

            orfs[i] = intervals.Interval(start, end)

        return orfs


    def all_orfs(self, min_length=300):
        '''Finds all open reading frames in the sequence, that are at least as
           long as min_length. Includes ORFs on the reverse strand.
           Returns a list of ORFs, where each element is a tuple:
           (interval.Interval, bool)
           where bool=True means on the reverse strand'''
        orfs = []
        for frame in [0,1,2]:
            for revcomp in [False, True]:
                orfs.extend([(t, revcomp) for t in self.orfs(frame=frame, revcomp=revcomp) if len(t)>=min_length])

        return sorted(orfs, key=lambda t:t[0])


    def is_complete_orf(self):
        '''Returns true iff length is >= 6, is a multiple of 3, and there is exactly one stop codon in the sequence and it is at the end'''
        if len(self) %3 != 0 or len(self) < 6:
            return False

        orfs = self.orfs()
        complete_orf = intervals.Interval(0, len(self) - 1)
        for orf in orfs:
            if orf == complete_orf:
                return True
        return False


    def looks_like_gene(self):
        '''Returns true iff: length >=6, length is a multiple of 3, first codon is start, last codon is a stop and has no other stop codons'''
        return self.is_complete_orf() \
          and len(self) >= 6 \
          and len(self) %3 == 0 \
          and self.seq[0:3].upper() in genetic_codes.starts[genetic_code]


    def make_into_gene(self):
        '''Tries to make into a gene sequence. Tries all three reading frames and both strands. Returns a tuple (new sequence, strand, frame) if it was successful. Otherwise returns None.'''
        for reverse in [True, False]:
            for frame in range(3):
                new_seq = copy.copy(self)
                if reverse:
                    new_seq.revcomp()
                new_seq.seq = new_seq[frame:]
                if len(new_seq) % 3:
                    new_seq.seq = new_seq.seq[:-(len(new_seq) % 3)]

                new_aa_seq = new_seq.translate()
                if len(new_aa_seq) >= 2 and new_seq[0:3] in genetic_codes.starts[genetic_code] and new_aa_seq[-1] == '*' and '*' not in new_aa_seq[:-1]:
                    strand = '-' if reverse else '+'
                    return new_seq, strand, frame

        return None


    # Fills the object with the next sequence in the file. Returns
    # True if this was successful, False if no more sequences in the file.
    # If reading a file of quality scores, set read_quals = True
    def get_next_from_file(self, f, read_quals=False):
        if f in previous_lines:
            if previous_lines[f] == None:
                self.id = self.seq = None
                return False
            else:
                self.id = self._get_id_from_header_line(previous_lines[f])
        else:
            line = '\n'
            while line == '\n':
                line = f.readline()
            self.id = self._get_id_from_header_line(line)

        self.seq = ''
        seq_lines = [] # much faster to store the seq lines in an array,
                       # then join at the end

        while 1:
            line = f.readline()

            if line.startswith('>'):
                previous_lines[f] = line.rstrip()
                break
            elif line == '':
                previous_lines[f] = None
                break
            else:
                 seq_lines.append(line.rstrip())

        if read_quals:
            self.seq = ' '.join(seq_lines)
        else:
            self.seq = ''.join(seq_lines)
        return True

    def __str__(self):
        if Fasta.line_length == 0:
            return '>' + self.id + '\n' + self.seq
        else:
            return '>' + self.id + '\n' + '\n'.join(self.seq[i:i+Fasta.line_length] for i in range(0, len(self), Fasta.line_length))

    def __getitem__(self, index):
        return self.seq[index]

    def trim(self, start, end):
        '''Removes first 'start'/'end' bases off the start/end of the sequence'''
        self.seq = self.seq[start:len(self.seq) - end]

    # qual_scores should be a list of quality scores
    def to_Fastq(self, qual_scores):
        '''Returns a Fastq object. qual_scores expected to be a list of numbers, like you would get in a .qual file'''
        if len(self) != len(qual_scores):
            raise Error('Error making Fastq from Fasta, lengths differ.', self.id)
        return Fastq(self.id, self.seq, ''.join([chr(max(0, min(x, 93)) + 33) for x in qual_scores]))

    def search(self, search_string):
        '''Finds every occurrence (including overlapping ones) of the search_string, including on the reverse strand. Returns a list where each element is a tuple (position, strand) where strand is in ['-', '+']. Positions are zero-based'''
        seq = self.seq.upper()
        search_string = search_string.upper()
        pos = 0
        found = seq.find(search_string, pos)
        hits = []

        while found != -1:
            hits.append((found, '+'))
            pos = found + 1
            found = seq.find(search_string, pos)


        pos = 0
        search_string = Fasta('x', search_string)
        search_string.revcomp()
        search_string = search_string.seq
        found = seq.find(search_string, pos)

        while found != -1:
            hits.append((found, '-'))
            pos = found + 1
            found = seq.find(search_string, pos)

        return hits

    def translate(self, frame=0):
        '''Returns a Fasta sequence, translated into amino acids. Starts translating from 'frame', where frame expected to be 0,1 or 2'''
        return Fasta(self.id, ''.join([genetic_codes.codes[genetic_code].get(self.seq[x:x+3].upper(), 'X') for x in range(frame, len(self)-1-frame, 3)]))


class Embl(Fasta):
    '''Exactly the same as Fasta, but reading seqs from a file works differently'''
    def __eq__(self, other):
        return type(other) in [Fasta, Embl] and  type(self) in [Fasta, Embl] and self.__dict__ == other.__dict__

    def _get_id_from_header_line(self, line):
        if line.startswith('ID   ') and line[5] != ' ':
            return line.split()[1].rstrip(';')
        elif line.startswith('LOCUS'):
            return line.split()[1]
        else:
            raise Error('Error! expected line starting with "ID" or "LOCUS", but got this:\n', line)

    def get_next_from_file(self, f, read_quals=False):
        if f in previous_lines:
            line = ''
            if previous_lines[f] == None:
                self.id = self.seq = None
                return False
            else:
                self.id = self._get_id_from_header_line(previous_lines[f])
        else:
            line = '\n'
            while line == '\n':
                line = f.readline()
            self.id = self._get_id_from_header_line(line)

        self.seq = ''
        seq_lines = []

        while not (line.startswith('SQ') or line.rstrip() == 'ORIGIN'):
            line = f.readline()
            if line == '':
                raise Error('Error! No SQ or ORIGIN line found for sequence ' + self.id)

        line = f.readline()

        while not line.startswith('//'):
            if line == '' or line[0] != ' ':
                raise Error('Error! Did not find end of sequence ' + self.id)
            seq_lines.append(''.join(line.rstrip().strip(' 0123456789').split()))
            line = f.readline()


        while 1:
            if line.startswith('ID') or line.startswith('LOCUS'):
                previous_lines[f] = line.rstrip()
                break
            elif line == '':
                previous_lines[f] = None
                break

            line = f.readline()

        self.seq = ''.join(seq_lines)
        return True

class Fastq(Fasta):
    '''Class to store and manipulate FASTQ sequences. They have three things: a name, sequence and string of quality scores'''
    def __init__(self, id_in=None, seq_in=None, qual_in=None):
        super().__init__(id_in, seq_in)
        self.qual = qual_in
        if (not self.seq == self.qual == None) and len(self.qual) != len(self.seq):
            raise Error('Error constructing Fastq. Mismatch in sequence and quality length\n' + str(self))

    def __str__(self):
        return '@' + self.id + '\n' + self.seq + '\n+\n' + self.qual

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def subseq(self, start, end):
        '''Returns Fastq object with the same name, of the bases from start to end, but not including end'''
        return Fastq(self.id, self.seq[start:end], self.qual[start:end])

    def get_next_from_file(self, f, read_quals=False):
        if f in previous_lines:
            line = previous_lines[f]
            del previous_lines[f]
        else:
            line = f.readline()

        while line == '\n':
            line = f.readline()

        if not line:
            self = Fastq('', '', '')
            return False

        if not line.startswith('@'):
            raise Error('Error getting next sequence from fastq file. Got line:\n' + line)

        self.id = line.rstrip()[1:]
        line = f.readline()
        if not line:
            raise Error('Error getting next sequence from fastq file, sequence has ID ' + self.id)

        self.seq = line.strip()

        line = f.readline()
        if not (line and line.startswith('+')):
            raise Error('Error getting next sequence from fastq file, no line starting with +,  sequence has ID ' + self.id)

        line = f.readline()
        if not line:
            raise Error('Error getting next sequence from fastq file, sequence has ID ' + self.id)

        self.qual = line.rstrip()
        return True

    def revcomp(self):
        '''Reverse complements the sequence'''
        super().revcomp()
        self.qual = self.qual[::-1]

    def trim(self, start, end):
        '''Removes first 'start'/'end' bases off the start/end of the sequence'''
        super().trim(start, end)
        self.qual = self.qual[start:len(self.qual) - end]

    def to_Fasta_and_qual(self):
        quals = [ord(x) - 33 for x in self.qual]
        return (Fasta(self.id, self.seq), quals)

    def expand_nucleotides(self):
        return [Fastq(x.id, x.seq, self.qual) for x in super().expand_nucleotides()]

    def trim_Ns(self):
        '''Removes any leading or trailing N or n characters from the sequence'''
        # get index of first base that is not an N
        i = 0
        while i < len(self) and self.seq[i] in 'nN':
            i += 1

        # strip off start of sequence and quality
        self.seq = self.seq[i:]
        self.qual = self.qual[i:]

        # strip the ends
        self.seq = self.seq.rstrip('Nn')
        self.qual = self.qual[:len(self.seq)]

    def replace_interval(self, start, end, new, qual_string):
        '''Replaces the sequence from start to end with the sequence "new"'''
        if len(new) != len(qual_string):
            raise Error('Length of new seq and qual string in replace_interval() must be equal. Cannot continue')
        super().replace_interval(start, end, new)
        self.qual = self.qual[0:start] + qual_string + self.qual[end + 1:]

    def translate(self):
        '''Returns a Fasta sequence, translated into amino acids. Starts translating from 'frame', where frame expected to be 0,1 or 2'''
        fa = super().translate()
        return Fastq(fa.id, fa.seq, 'I'*len(fa.seq))

    def make_into_gene(self):
        got = super().make_into_gene()
        if got is None:
            return None
        seq, strand, frame = got
        new_seq = copy.copy(self)

        if strand == '-':
            new_seq.revcomp()

        new_seq.seq = new_seq.seq[frame:frame + len(seq)]
        new_seq.qual = new_seq.qual[frame:frame + len(seq)]
        return new_seq, strand, frame


def _orfs_from_aa_seq(seq):
    orfs = []
    pos = 0
    while pos < len(seq):
        next_stop = seq.find('*', pos)
        if next_stop == -1:
            orfs.append(intervals.Interval(pos, len(seq)-1))
            break
        elif next_stop > pos:
            orfs.append(intervals.Interval(pos, next_stop))
        pos = next_stop + 1
    return orfs
