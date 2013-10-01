#!/usr/bin/env python3

import sys
import filecmp
import os
import unittest
from fastaq import sequences, utils, intervals

modules_dir = os.path.dirname(os.path.abspath(sequences.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class Error (Exception): pass

expected_embl = [
    'aaacaaaccaaatatggattttattgtagccatatttgctctgtttgttattagctcattcacaattacttccacaaatgcagttgaagcttctactcttcttgacataggtaacctgagtcggagcagttttcctcgtggcttcatctttggtgctggatcttcagcataccaatttgaaggtgcagtaaacgaaggcggtagaggaccaagtatttgggataccttcacccataaatatccagaaaaaataagggatggaagcaatgcagacatcacggttgaccaatatcaccgctacaaggaagatgttgggattatgaaggatcaaaatatggattcgtatagattctcaatctcttggccaagaatactcccaaagggaaagttgagcggaggcataaatcacgaaggaatcaaatattacaacaaccttatcaacgaactattggctaacggtatacaaccatttgtaactctttttcattgggatcttccccaagtcttagaagatgagtatggtggtttcttaaactccggtgtaataaatgattttcgagactatacggatctttgcttcaaggaatttggagatagagtgaggtattggagtactctaaatgagccatgggtgtttagcaattctggatatgcactaggaacaaatgcaccaggtcgatgttcggcctccaacgtggccaagcctggtgattctggaacaggaccttatatagttacacacaatcaaattcttgctcatgcagaagctgtacatgtgtataagactaaataccaggcatatcaaaagggaaagataggcataacgttggtatctaactggttaatgccacttgatgataatagcataccagatataaaggctgccgagagatcacttgacttccaatttggattgtttatggaacaattaacaacaggagattattctaagagcatgcggcgtatagttaaaaaccgattacctaagttctcaaaattcgaatcaagcctagtgaatggttcatttgattttattggtataaactattactcttctagttatattagcaatgccccttcacatggcaatgccaaacccagttactcaacaaatcctatgaccaatatttcatttgaaaaacatgggatacccttaggtccaagggctgcttcaatttggatatatgtttatccatatatgtttatccaagaggacttcgagatcttttgttacatattaaaaataaatataacaatcctgcaattttcaatcactgaaaatggtatgaatgaattcaacgatgcaacacttccagtagaagaagctcttttgaatacttacagaattgattactattaccgtcacttatactacattcgttctgcaatcagggctggctcaaatgtgaagggtttttacgcatggtcatttttggactgtaatgaatggtttgcaggctttactgttcgttttggattaaactttgtagattagaaagatggattaaaaaggtaccctaagctttctgcccaatggtacaagaactttctcaaaagaaactagctagtattattaaaagaactttgtagtagattacagtacatcgtttgaagttgagttggtgcacctaattaaataaaagaggttactcttaacatatttttaggccattcgttgtgaagttgttaggctgttatttctattatactatgttgtagtaataagtgcattgttgtaccagaagctatgatcataactataggttgatccttcatgtatcagtttgatgttgagaatactttgaattaaaagtctttttttatttttttaaaaaaaaaaaaaaaaaaaaaaaaaaaaa',
    'aaacaaaccaaatatggattttattgtagccatatttgctctgtttgttattagctcattcacaattacttccacaaatgcagttgaagcttctactcttcttgacataggtaacctgagtcggagcagttttcctcgtggcttcatctttggtgctggatcttcagcataccaatttgaaggtgcagtaaacgaaggcggtagaggaccaagtatttgggataccttcacccataaatatccagaaaaaataagggatggaagcaatgcagacatcacggttgaccaatatcaccgctacaaggaagatgttgggattatgaaggatcaaaatatggattcgtatagattctcaatctcttggccaagaatactcccaaagggaaagttgagcggaggcataaatcacgaaggaatcaaatattacaacaaccttatcaacgaactattggctaacggtatacaaccatttgtaactctttttcattgggatcttccccaagtcttagaagatgagtatggtggtttcttaaactccggtgtaataaatgattttcgagactatacggatctttgcttcaaggaatttggagatagagtgaggtattggagtactctaaatgagccatgggtgtttagcaattctggatatgcactaggaacaaatgcaccaggtcgatgttcggcctccaacgtggccaagcctggtgattctggaacaggaccttatatagttacacacaatcaaattcttgctcatgcagaagctgtacatgtgtataagactaaataccaggcatatcaaaagggaaagataggcataacgttggtatctaactggttaatgccacttgatgataatagcataccagatataaaggctgccgagagatcacttgacttccaatttggattgtttatggaacaattaacaacaggagattattctaagagcatgcggcgtatagttaaaaaccgattacctaagttctcaaaattcgaatcaagcctagtgaatggttcatttgattttattggtataaactattactcttctagttatattagcaatgccccttcacatggcaatgccaaacccagttactcaacaaatcctatgaccaatatttcatttgaaaaacatgggatacccttaggtccaagggctgcttcaatttggatatatgtttatccatatatgtttatccaagaggacttcgagatcttttgttacatattaaaaataaatataacaatcctgcaattttcaatcactgaaaatggtatgaatgaattcaacgatgcaacacttccagtagaagaagctcttttgaatacttacagaattgattactattaccgtcacttatactacattcgttctgcaatcagggctggctcaaatgtgaagggtttttacgcatggtcatttttggactgtaatgaatggtttgcaggctttactgttcgttttggattaaactttgtagattagaaagatggattaaaaaggtaccctaagctttctgcccaatggtacaagaactttctcaaaagaaactagctagtattattaaaagaactttgtagtagattacagtacatcgtttgaagttgagttggtgcacctaattaaataaaagaggttactcttaacatatttttaggccattcgttgtgaagttgttaggctgttatttctattatactatgttgtagtaataagtgcattgttgtaccagaagctatgatcataactataggttgatccttcatgtatcagtttgatgttgagaatactttgaattaaaagtctttttttatttttttaaaaaaaaaaaaaaaaaaaaccccccccc',
]
class TestFasta(unittest.TestCase):
    def setUp(self):
        self.fasta = sequences.Fasta('ID', 'ACGTA')

    def test_equality(self):
        self.assertTrue(self.fasta == sequences.Fasta('ID', 'ACGTA'))
        self.assertFalse(self.fasta == sequences.Fasta('I', 'ACGTA'))
        self.assertFalse(self.fasta == sequences.Fasta('ID', 'ACGT'))
        self.assertFalse(self.fasta != sequences.Fasta('ID', 'ACGTA'))
        self.assertTrue(self.fasta != sequences.Fasta('I', 'ACGTA'))
        self.assertTrue(self.fasta != sequences.Fasta('ID', 'ACGT'))

    def test_init(self):
        '''__init__ should get the ID and sequence correctly'''
        self.assertEqual(self.fasta.id, 'ID')
        self.assertEqual(self.fasta.seq, 'ACGTA')

    def test_get_next_from_file(self):
        '''get_next_from_file() should read seqs from OK, including weirdness in file'''
        f_in = utils.open_file_read(os.path.join(data_dir, 'sequences_test.fa'))
        fa = sequences.Fasta()
        counter = 1

        while fa.get_next_from_file(f_in):
            self.assertEqual(fa, sequences.Fasta(str(counter), 'ACGTA'))
            counter += 1

        utils.close(f_in)

    def test_get_id_from_header_line(self):
        '''Check that can get ID from header line or die properly'''
        self.assertEqual(sequences.Fasta._get_id_from_header_line(self.fasta, '>X'), 'X')
        with self.assertRaises(sequences.Error):
            self.assertEqual(sequences.Fasta._get_id_from_header_line(self.fasta, 'X'), 'X')

    def test_getitem(self):
        '''getitem() should return the right subsequence'''
        seq = 'AACGTGTCA'
        fa = sequences.Fasta('x', seq)
        self.assertEqual(seq[1], fa[1])
        self.assertEqual(seq[0:2], fa[0:2])
        self.assertEqual(seq[1:], fa[1:])

    def test_len(self):
        '''len() should return the length of the sequence'''
        self.assertEqual(5, len(self.fasta))

    def test_print_line_length(self):
        '''__str__ should be formatted correctly with the right number of chars per line of sequence'''
        line_lengths = [0, 3]
        correct_files = [os.path.join(data_dir, x) for x in ['sequences_test_one-per-line.fa', 'sequences_test_3-per-line.fa']]

        for i in range(len(line_lengths)):
            seq_reader = sequences.file_reader(os.path.join(data_dir, 'sequences_test_one-per-line.fa'))
            sequences.Fasta.line_length = line_lengths[i]
            tmp_out = 'tmp.line_length_test.fa'
            f = utils.open_file_write(tmp_out)
            for s in seq_reader:
                print(s, file=f)
            utils.close(f)
            self.assertTrue(filecmp.cmp(correct_files[i], tmp_out))
            os.unlink(tmp_out)

        sequences.Fasta.line_length = 60

    def test_strip_after_first_whitespace(self):
        '''Test strip_after_first_whitespace()'''
        seqs = [
            sequences.Fasta('name', 'A'),
            sequences.Fasta('name foo', 'A'),
            sequences.Fasta('name foo bar', 'A'),
            sequences.Fasta('name\tfoo', 'A'),
        ]

        for seq in seqs:
            seq.strip_after_first_whitespace()

        for seq in seqs:
            self.assertEqual(seq.id, 'name')

    def test_strip_illumina_suffix(self):
        '''Check that /1 and /2 removed correctly from IDs'''
        seqs = [sequences.Fasta('name/1', 'A'),
                sequences.Fasta('name/2', 'A'),
                sequences.Fasta('name', 'A'),
                sequences.Fasta('name/1/2', 'A'),
                sequences.Fasta('name/2/1', 'A'),
                sequences.Fasta('name/3', 'A')]

        correct_names = ['name', 'name', 'name', 'name/1', 'name/2', 'name/3']

        for seq in seqs:
            seq.strip_illumina_suffix()

        for i in range(len(seqs)):
            self.assertEqual(seqs[i].id, correct_names[i])

    def test_revcomp(self):
        '''revcomp() should correctly reverse complement a sequence'''
        fa = sequences.Fasta('ID', 'ACGTNacgtn')
        fa.revcomp()
        self.assertEqual(fa, sequences.Fasta('ID', 'nacgtNACGT'))

    def test_gaps(self):
        '''gaps() should find the gaps in a sequence correctly'''
        test_seqs = [sequences.Fasta('ID', 'ACGT'),
                     sequences.Fasta('ID', 'NACGT'),
                     sequences.Fasta('ID', 'NACGTN'),
                     sequences.Fasta('ID', 'ANNCGT'),
                     sequences.Fasta('ID', 'NANNCGTNN')]

        correct_gaps = [[],
                        [intervals.Interval(0, 0)],
                        [intervals.Interval(0, 0), intervals.Interval(5, 5)],
                        [intervals.Interval(1, 2)],
                        [intervals.Interval(0, 0), intervals.Interval(2, 3), intervals.Interval(7, 8)]]

        for i in range(len(test_seqs)):
            gaps = test_seqs[i].gaps()
            self.assertListEqual(correct_gaps[i], gaps)

    def test_contig_coords(self):
        '''contig_coords() should get the coords of all contigs in a sequence correctly'''
        test_seqs = [sequences.Fasta('ID', 'ACGT'),
                     sequences.Fasta('ID', 'NACGT'),
                     sequences.Fasta('ID', 'NNACGT'),
                     sequences.Fasta('ID', 'ACGTN'),
                     sequences.Fasta('ID', 'ACGTNN'),
                     sequences.Fasta('ID', 'NANNCGT'),
                     sequences.Fasta('ID', 'ANNCGTNNAAAAA')]

        correct_coords = [[intervals.Interval(0,3)],
                         [intervals.Interval(1, 4)],
                         [intervals.Interval(2, 5)],
                         [intervals.Interval(0, 3)],
                         [intervals.Interval(0, 3)],
                         [intervals.Interval(1, 1), intervals.Interval(4,6)],
                         [intervals.Interval(0, 0), intervals.Interval(3, 5), intervals.Interval(8, 12)]]

        for i in range(len(test_seqs)):
            gaps = test_seqs[i].contig_coords()
            self.assertListEqual(correct_coords[i], gaps)

    def test_is_all_Ns(self):
        '''Test is_all_Ns()'''
        self.assertTrue(sequences.Fasta('ID', 'n').is_all_Ns())
        self.assertTrue(sequences.Fasta('ID', 'N').is_all_Ns())
        self.assertTrue(sequences.Fasta('ID', 'nNn').is_all_Ns())
        self.assertFalse(sequences.Fasta('ID', 'a').is_all_Ns())
        self.assertFalse(sequences.Fasta('ID', '').is_all_Ns())
        self.assertFalse(sequences.Fasta('ID', 'anNg').is_all_Ns())
        self.assertFalse(sequences.Fasta('ID', 'naN').is_all_Ns())
        self.assertFalse(sequences.Fasta('ID', 'anNg').is_all_Ns(start=0, end=0))
        self.assertFalse(sequences.Fasta('ID', 'anNg').is_all_Ns(start=0, end=1))
        self.assertTrue(sequences.Fasta('ID', 'anNg').is_all_Ns(start=1, end=1))
        self.assertTrue(sequences.Fasta('ID', 'anNg').is_all_Ns(start=1, end=2))
        self.assertFalse(sequences.Fasta('ID', 'anNg').is_all_Ns(start=1))
        self.assertTrue(sequences.Fasta('ID', 'anN').is_all_Ns(start=1))
        self.assertFalse(sequences.Fasta('ID', 'anNg').is_all_Ns(end=1))
        self.assertTrue(sequences.Fasta('ID', 'nNA').is_all_Ns(end=1))

        with self.assertRaises(sequences.Error):
            sequences.Fasta('ID', 'anNg').is_all_Ns(start=1, end=0)

    def test_trim_Ns(self):
        '''trim_Ns() should do the right trimming of a sequence'''
        fa = sequences.Fasta('ID', 'ANNANA')
        test_seqs = [sequences.Fasta('ID', 'ANNANA'),
                     sequences.Fasta('ID', 'NANNANA'),
                     sequences.Fasta('ID', 'NANNANAN'),
                     sequences.Fasta('ID', 'ANNANAN'),
                     sequences.Fasta('ID', 'NNNNNNANNANAN'),
                     sequences.Fasta('ID', 'NNANNANANn')]

        for s in test_seqs:
            s.trim_Ns()
            self.assertEqual(fa, s)

    def test_replace_bases(self):
        '''Check that bases get replaced correctly'''
        fa = sequences.Fasta('X', 'AUCGTUUACT')
        fa.replace_bases('U', 'T')
        self.assertEqual(fa, sequences.Fasta('X', 'ATCGTTTACT'))

    def test_replace_interval(self):
        '''Test replace_interval()'''
        fa = sequences.Fasta('ID', 'ACGTA')
        fa.replace_interval(0, 0, 'NEW')
        self.assertEqual(fa, sequences.Fasta('ID', 'NEWCGTA'))

        fa = sequences.Fasta('ID', 'ACGTA')
        fa.replace_interval(4, 4, 'NEW')
        self.assertEqual(fa, sequences.Fasta('ID', 'ACGTNEW'))

        fa = sequences.Fasta('ID', 'ACGTA')
        fa.replace_interval(2, 3, 'NEW')
        self.assertEqual(fa, sequences.Fasta('ID', 'ACNEWA'))

        fa = sequences.Fasta('ID', 'ACGTA')
        with self.assertRaises(sequences.Error):
            fa.replace_interval(3,2,'x')
        with self.assertRaises(sequences.Error):
            fa.replace_interval(1,5,'x')
        with self.assertRaises(sequences.Error):
            fa.replace_interval(5,10,'x')

        fq = sequences.Fastq('ID', 'ACGTA', 'ABCDE')
        fq.replace_interval(0, 0, 'NEW', 'III')
        self.assertEqual(fq, sequences.Fastq('ID', 'NEWCGTA', 'IIIBCDE'))

        fq = sequences.Fastq('ID', 'ACGTA', 'ABCDE')
        fq.replace_interval(4, 4, 'NEW', 'III')
        self.assertEqual(fq, sequences.Fastq('ID', 'ACGTNEW', 'ABCDIII'))

        fq = sequences.Fastq('ID', 'ACGTA', 'ABCDE')
        fq.replace_interval(2, 3, 'NEW', 'III')
        self.assertEqual(fq, sequences.Fastq('ID', 'ACNEWA', 'ABIIIE'))

        with self.assertRaises(sequences.Error):
            fq.replace_interval(1,1,'x', 'xx')

    def test_search_string(self):
        '''Check that search_string() finds all the hits'''
        fa = sequences.Fasta('X', 'AAA')
        hits = fa.search('G')
        self.assertTrue(len(hits) == 0)
        hits = fa.search('AAA')
        self.assertListEqual(hits, [(0, '+')])
        hits = fa.search('AA')
        self.assertListEqual(hits, [(0, '+'), (1, '+')])
        hits = fa.search('TTT')
        self.assertListEqual(hits, [(0, '-')])

    def test_to_Fastq(self):
        '''Check to_Fastq converts OK, including out of range quality scores'''
        fa = sequences.Fasta('X', 'AAAAA')
        quals = [-1, 0, 40, 93, 94]
        self.assertEqual(sequences.Fastq('X', 'AAAAA', '!!I~~'), fa.to_Fastq(quals))
        with self.assertRaises(sequences.Error):
            fa.to_Fastq('AAAAAAAAAAAAA')


    def test_translate(self):
        '''Test nucleotide -> amino acid conversion works on Fasta'''
        fa = sequences.Fasta('ID', 'GCAGCCGCGGCTAGAAGGCGACGCCGGCGTAACAATGACGATTGCTGTGAAGAGCAACAGGGAGGCGGGGGTCACCATATAATCATTTTATTGCTACTCCTGCTTAAAAAGATGTTCTTTCCACCCCCGCCTAGCAGTTCATCCTCGTCTACAACCACGACTTGGTACTATGTAGTCGTGGTTTAATAGTGA')
        self.assertEqual(sequences.Fasta('ID', 'AAAARRRRRRNNDDCCEEQQGGGGHHIIILLLLLLKKMFFPPPPSSSSSSTTTTWYYVVVV***'), fa.translate())
        self.assertEqual(sequences.Fasta('ID', 'QPRLEGDAGVTMTIAVKSNREAGVTI*SFYCYSCLKRCSFHPRLAVHPRLQPRLGTM*SWFNS'), fa.translate(frame=1))
        print(fa.translate(frame=1))
        self.assertEqual(sequences.Fasta('ID', 'SRG*KATPA*Q*RLL*RATGRRGSPYNHFIATPA*KDVLSTPA*QFILVYNHDLVLCSRGLIV'), fa.translate(frame=2))


    def test_split_capillary_id(self):
        '''Tests that we get information from a sanger capillary read name OK'''
        ids = ['abcde.p1k', 'abcde.x.p1k', 'abcde.p1ka', 'abcde.q1k', 'abcde.w2k']
        expected = [{'prefix': 'abcde', 'dir': 'fwd', 'suffix': 'p1k'},
                    {'prefix': 'abcde.x', 'dir': 'fwd', 'suffix': 'p1k'},
                    {'prefix': 'abcde', 'dir': 'fwd', 'suffix': 'p1ka'},
                    {'prefix': 'abcde', 'dir': 'rev', 'suffix': 'q1k'},
                    {'prefix': 'abcde', 'dir': 'unk', 'suffix': 'w2k'}]

        for i in range(len(ids)):
            fa = sequences.Fasta(ids[i], 'A')
            self.assertEqual(fa.split_capillary_id(), expected[i])

        with self.assertRaises(sequences.Error):
            fa = sequences.Fasta('name', 'A')
            fa.split_capillary_id()


class TestEmbl(unittest.TestCase):
    def test_get_id_from_header_line(self):
        '''Test get id from header line of EMBL'''
        embl = sequences.Embl('ID', 'ACGT')
        self.assertEqual(embl._get_id_from_header_line('ID   X; blah'), 'X')
        with self.assertRaises(sequences.Error):
            self.assertEqual(embl._get_id_from_header_line('ID X;'), 'X')
        with self.assertRaises(sequences.Error):
            self.assertEqual(embl._get_id_from_header_line('XX   X;'), 'X')

    def test_get_next_from_file(self):
        f_in = utils.open_file_read(os.path.join(data_dir, 'sequences_test.embl'))
        embl = sequences.Embl()
        counter = 1

        while embl.get_next_from_file(f_in):
            self.assertEqual(embl, sequences.Fasta('seq' + str(counter), expected_embl[counter-1]))
            counter += 1

        utils.close(f_in)


class TestFastq(unittest.TestCase):
    def setUp(self):
        self.fastq = sequences.Fastq('ID', 'ACGTA', 'IIIII')

    def test_init(self):
        '''__init__ should get the ID, sequence and quality correctly'''
        self.assertEqual(self.fastq.id, 'ID')
        self.assertEqual(self.fastq.seq, 'ACGTA')
        self.assertEqual(self.fastq.qual, 'IIIII')

    def test_init_length_mismatch(self):
        '''__init__ should raise an error when length of seq and quality not the same'''
        with self.assertRaises(sequences.Error):
            sequences.Fastq('X', 'A', 'II')

    def test_get_next_from_file(self):
        '''get_next_from_file() should read seqs from OK, and raise error at badly formatted file'''
        bad_files = ['sequences_test_fail_no_AT.fq',
                     'sequences_test_fail_no_seq.fq',
                     'sequences_test_fail_no_plus.fq',
                     'sequences_test_fail_no_qual.fq']

        bad_files = [os.path.join(data_dir, x) for x in bad_files]

        for fname in bad_files:
            f_in = utils.open_file_read(fname)
            fq = sequences.Fastq()
            with self.assertRaises(sequences.Error):
                while fq.get_next_from_file(f_in):
                    pass

            utils.close(f_in)

        fname = os.path.join(data_dir, 'sequences_test_good_file.fq')
        try:
            f_in = open(fname)
        except IOError:
            print("Error opening '" + fname + "'", file=sys.stderr)
            sys.exit(1)

        fq = sequences.Fastq()
        while fq.get_next_from_file(f_in):
            self.assertEqual(fq, sequences.Fastq('ID', 'ACGTA', 'IIIII'))
        utils.close(f_in)

    def test_revcomp(self):
        '''revcomp() should correctly reverse complement a sequence'''
        fq = sequences.Fastq('ID', 'ACGTNacgtn', '1234567890')
        fq.revcomp()
        self.assertEqual(fq, sequences.Fastq('ID', 'nacgtNACGT', '0987654321'))

    def test_trim_Ns(self):
        '''trim_Ns() should do the right trimming of a fastq sequence'''
        fq = sequences.Fastq('ID', 'ANNANA', '111111')
        test_seqs = [sequences.Fastq('ID', 'ANNANA', '111111'),
                     sequences.Fastq('ID', 'NANNANA', '1111111'),
                     sequences.Fastq('ID', 'NANNANAN', '11111111'),
                     sequences.Fastq('ID', 'ANNANAN', '1111111'),
                     sequences.Fastq('ID', 'NNNNNNANNANAN', '1111111111111'),
                     sequences.Fastq('ID', 'NNANNANANn', '1111111111')]

        for s in test_seqs:
            s.trim_Ns()
            self.assertEqual(fq, s)

    def test_trim(self):
        '''trim() should trim the right number of bases off start and end'''
        fq = sequences.Fastq('ID', '1234567890', '1234567890')
        fq.trim(0, 0)
        self.assertEqual(fq, sequences.Fastq('ID', '1234567890', '1234567890'))

        fq = sequences.Fastq('ID', '1234567890', '1234567890')
        fq.trim(1, 0)
        self.assertEqual(fq, sequences.Fastq('ID', '234567890', '234567890'))

        fq = sequences.Fastq('ID', '1234567890', '1234567890')
        fq.trim(0, 1)
        self.assertEqual(fq, sequences.Fastq('ID', '123456789', '123456789'))

        fq = sequences.Fastq('ID', '1234567890', '1234567890')
        fq.trim(2, 2)
        self.assertEqual(fq, sequences.Fastq('ID', '345678', '345678'))

    def test_to_Fasta_and_qual(self):
        '''Check to_Fasta_and_qual converts quality scores correctly'''
        fq = sequences.Fastq('ID', 'ACGT', '>ADI')
        (fa, qual) = fq.to_Fasta_and_qual()
        self.assertEqual(fa, sequences.Fasta('ID', 'ACGT'))
        self.assertListEqual(qual, [29, 32, 35, 40])


    def test_translate(self):
        '''Test nucleatide -> amino acid conversion works on Fasta'''
        fq = sequences.Fastq('ID', 'GCAGCCGCGGCTAGAAGGCGACGCCGGCGTAACAATGACGATTGCTGTGAAGAGCAACAGGGAGGCGGGGGTCACCATATAATCATTTTATTGCTACTCCTGCTTAAAAAGATGTTCTTTCCACCCCCGCCTAGCAGTTCATCCTCGTCTACAACCACGACTTGGTACTATGTAGTCGTGGTTTAATAGTGA', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII')

        self.assertEqual(sequences.Fastq('ID', 'AAAARRRRRRNNDDCCEEQQGGGGHHIIILLLLLLKKMFFPPPPSSSSSSTTTTWYYVVVV***', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'), fq.translate())

class TestFileReader(unittest.TestCase):
    def test_file_reader_fasta(self):
        '''file_reader should iterate through a fasta file correctly'''
        reader = sequences.file_reader(os.path.join(data_dir, 'sequences_test.fa'))
        counter = 1
        for seq in reader:
            self.assertEqual(seq, sequences.Fasta(str(counter), 'ACGTA'))
            counter += 1

    def test_file_reader_fastq(self):
        '''file_reader should iterate through a fastq file correctly'''
        reader = sequences.file_reader(os.path.join(data_dir, 'sequences_test_good_file.fq'))
        for seq in reader:
            self.assertEqual(seq, sequences.Fastq('ID', 'ACGTA', 'IIIII'))

    def test_file_reader_bad_format(self):
        '''file_reader should die properly when not given fasta or fastq file'''
        with self.assertRaises(sequences.Error):
            reader = sequences.file_reader(os.path.join(data_dir, 'sequences_test_not_a_fastaq_file'))
            for seq in reader:
                pass

    def test_file_reader_gff(self):
        '''Test read gff file'''
        reader = sequences.file_reader(os.path.join(data_dir, 'sequences_test_gffv3.gff'))
        counter = 1
        for seq in reader:
            self.assertEqual(seq, sequences.Fasta('seq' + str(counter), 'ACGTACGTAC'))
            counter += 1
        
        bad_files = [
            'sequences_test_gffv3.no_seq.gff',
            'sequences_test_gffv3.no_seq.2.gff'
        ]
        bad_files = [os.path.join(data_dir, x) for x in bad_files]

        for filename in bad_files:
            with self.assertRaises(sequences.Error):
                reader = sequences.file_reader(filename)
                for seq in reader:
                    pass

    def test_file_reader_embl(self):
        '''Test read embl file'''
        reader = sequences.file_reader(os.path.join(data_dir, 'sequences_test.embl'))

        counter = 1
        for seq in reader:
            self.assertEqual(seq, sequences.Fasta('seq' + str(counter), expected_embl[counter-1]))
            counter += 1
        
        bad_files = [
            'sequences_test.embl.bad',
            'sequences_test.embl.bad2',
        ]
        bad_files = [os.path.join(data_dir, x) for x in bad_files]

        for filename in bad_files:
            with self.assertRaises(sequences.Error):
                reader = sequences.file_reader(filename)
                for seq in reader:
                    pass

    def test_file_reader_phylip(self):
        '''Test read phylip file'''
        test_files = [
            'sequences_test_phylip.interleaved',
            'sequences_test_phylip.interleaved2',
            'sequences_test_phylip.sequential'
        ]

        test_files = [os.path.join(data_dir, f) for f in test_files]

        expected_seqs = [
            sequences.Fasta('Turkey', 'AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT'),
            sequences.Fasta('Salmo gair', 'AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT'),
            sequences.Fasta('H. Sapiens', 'ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA')
        ]

        for fname in test_files:
            reader = sequences.file_reader(fname)
            i = 0
            for seq in reader:
                self.assertEqual(expected_seqs[i].seq, seq.seq)
                self.assertEqual(expected_seqs[i].id, seq.id)
                #self.assertEqual(expected_seqs[i], seq)
                i += 1
        

if __name__ == '__main__':
    unittest.main()

