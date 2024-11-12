import filecmp
import os
import sys
import unittest

from pyfastaq import sequences, utils, intervals, tasks

data_dir = "tests/data"

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

    def test_subseq(self):
        '''Test subseq'''
        fa = sequences.Fasta('name', 'ACGTA')
        self.assertEqual(fa.subseq(1,4), sequences.Fasta('name', 'CGT'))
        self.assertEqual(fa.subseq(None,4), sequences.Fasta('name', 'ACGT'))
        self.assertEqual(fa.subseq(1,None), sequences.Fasta('name', 'CGTA'))

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
                     sequences.Fasta('ID', 'ACNNNGTNA'),
                     sequences.Fasta('ID', 'ANNCGTNNAAAAA')]

        correct_coords = [[intervals.Interval(0,3)],
                         [intervals.Interval(1, 4)],
                         [intervals.Interval(2, 5)],
                         [intervals.Interval(0, 3)],
                         [intervals.Interval(0, 3)],
                         [intervals.Interval(1, 1), intervals.Interval(4,6)],
                         [intervals.Interval(0, 1), intervals.Interval(5, 6), intervals.Interval(8, 8)],
                         [intervals.Interval(0, 0), intervals.Interval(3, 5), intervals.Interval(8, 12)]]

        for i in range(len(test_seqs)):
            gaps = test_seqs[i].contig_coords()
            self.assertListEqual(correct_coords[i], gaps)




    def test_orfs(self):
        '''Test orfs()'''
        test_seqs = [(sequences.Fasta('ID', 'AAACCCGG'), 0, False, [intervals.Interval(0,5)]),
                     (sequences.Fasta('ID', 'AAAACCCGG'), 1, False, [intervals.Interval(1,6)]),
                     (sequences.Fasta('ID', 'AAAAACCCGG'), 2, False, [intervals.Interval(2,7)]),
                     (sequences.Fasta('ID', 'CCGGGTTT'), 0, True, [intervals.Interval(2,7)]),
                     (sequences.Fasta('ID', 'CCGGGTTTT'), 1, True, [intervals.Interval(2,7)]),
                     (sequences.Fasta('ID', 'CCGGGTTTTT'), 2, True, [intervals.Interval(2,7)]),
                     (sequences.Fasta('ID', 'AAACCCTGA'), 0, False, [intervals.Interval(0,8)]),
                     (sequences.Fasta('ID', 'AAACCCTGATAG'), 0, False, [intervals.Interval(0,8)]),
                     (sequences.Fasta('ID', 'AAACCCTGA'), 1, False, [intervals.Interval(1,6)]),
                     (sequences.Fasta('ID', ''), 0, False, []),
                     (sequences.Fasta('ID', 'A'), 0, False, []),
                     (sequences.Fasta('ID', 'AA'), 0, False, []),
                     (sequences.Fasta('ID', 'AAA'), 0, False, [intervals.Interval(0,2)]),
                     (sequences.Fasta('ID', 'AAAAAA'), 0, False, [intervals.Interval(0,5)]),
                     (sequences.Fasta('ID', 'AAA'), 1, False, []),
                     (sequences.Fasta('ID', 'AAA'), 2, False, []),
                     (sequences.Fasta('ID', 'AAA'), 0, True, [intervals.Interval(0,2)]),
                     (sequences.Fasta('ID', 'AAA'), 1, True, []),
                     (sequences.Fasta('ID', 'AAA'), 2, True, []),
                     (sequences.Fasta('ID', 'TAA'), 0, False, []),
                     (sequences.Fasta('ID', 'CTA'), 0, True, [])]


        for t in test_seqs:
            orfs = t[0].orfs(frame=t[1], revcomp=t[2])
            self.assertListEqual(orfs, t[3])

    def test_all_orfs(self):
        '''Test all_orfs()'''
        d = {}
        tasks.file_to_dict(os.path.join(data_dir, 'sequences_test_orfs.fa'), d)
        seq = d['1']
        orfs = seq.all_orfs(min_length=120)
        expected = [
            (intervals.Interval(27, 221), False),
            (intervals.Interval(44, 226), False),
            (intervals.Interval(48, 170), True),
            (intervals.Interval(109, 240), False),
            (intervals.Interval(143, 265), True),
            (intervals.Interval(227, 421), False),
            (intervals.Interval(277, 432), True),
            (intervals.Interval(286, 477), False),
            (intervals.Interval(288, 518), True),
            (intervals.Interval(562, 702), False),
            (intervals.Interval(600, 758), False),
            (intervals.Interval(605, 817), False),
            (intervals.Interval(818, 937), False),
            (intervals.Interval(835, 987), False),
            (intervals.Interval(864, 998), False)
        ]

        self.assertEqual(len(orfs), len(expected))

        for i in range(len(orfs)):
            print(orfs[i][0], expected[i][0])
            self.assertEqual(orfs[i][0], expected[i][0])
            self.assertEqual(orfs[i][1], expected[i][1])


    def test_is_complete_orf(self):
        '''Test is_complete_orf'''
        tests = [
            (sequences.Fasta('ID', 'TTT'), False),
            (sequences.Fasta('ID', 'TTTTAA'), True),
            (sequences.Fasta('ID', 'TTTTAATAA'), False),
            (sequences.Fasta('ID', 'TTGTAA'), True),
            (sequences.Fasta('ID', 'TTTAAC'), True),
            (sequences.Fasta('ID', 'TGA'), False),
            (sequences.Fasta('ID', 'TGAA'), False),
        ]

        for t in tests:
            self.assertEqual(t[0].is_complete_orf(), t[1])


    def test_looks_like_gene(self):
        '''Test looks_like_gene'''
        tests = [
            (sequences.Fasta('ID', 'TTT'), False),
            (sequences.Fasta('ID', 'TTGTAA'), True),
            (sequences.Fasta('ID', 'ttgTAA'), True),
            (sequences.Fasta('ID', 'TTGTTTTAA'), True),
            (sequences.Fasta('ID', 'TTGTAATTTTAA'), False),
            (sequences.Fasta('ID', 'TTGTTTTGAA'), False),
        ]

        for t in tests:
            self.assertEqual(t[0].looks_like_gene(), t[1])

        sequences.genetic_code = 1
        self.assertFalse(sequences.Fasta('ID', 'ATTCAGTAA').looks_like_gene())
        sequences.genetic_code = 11
        self.assertTrue(sequences.Fasta('ID', 'ATTCAGTAA').looks_like_gene())
        sequences.genetic_code = 1


    def test_make_into_gene_fasta(self):
        '''Test make_into_gene fasta'''
        print('sequences.genetic_code', sequences.genetic_code)
        tests = [
            (sequences.Fasta('ID', 'T'), None),
            (sequences.Fasta('ID', 'TT'), None),
            (sequences.Fasta('ID', 'TTT'), None),
            (sequences.Fasta('ID', 'TTG'), None),
            (sequences.Fasta('ID', 'TAA'), None),
            (sequences.Fasta('ID', 'TTGAAATAA'), (sequences.Fasta('ID', 'TTGAAATAA'), '+', 0)),
            (sequences.Fasta('ID', 'TTGAAATAT'), None),
            (sequences.Fasta('ID', 'TTGTAA'), (sequences.Fasta('ID', 'TTGTAA'), '+', 0)),
            (sequences.Fasta('ID', 'TTGTAAA'), (sequences.Fasta('ID', 'TTGTAA'), '+', 0)),
            (sequences.Fasta('ID', 'TTGTAAAA'), (sequences.Fasta('ID', 'TTGTAA'), '+', 0)),
            (sequences.Fasta('ID', 'TTGTAAAAA'), None),
            (sequences.Fasta('ID', 'ATTGTAA'), (sequences.Fasta('ID', 'TTGTAA'), '+', 1)),
            (sequences.Fasta('ID', 'ATTGTAAA'), (sequences.Fasta('ID', 'TTGTAA'), '+', 1)),
            (sequences.Fasta('ID', 'ATTGTAAAA'), (sequences.Fasta('ID', 'TTGTAA'), '+', 1)),
            (sequences.Fasta('ID', 'ATTGTAAAAA'), None),
            (sequences.Fasta('ID', 'AATTGTAA'), (sequences.Fasta('ID', 'TTGTAA'), '+', 2)),
            (sequences.Fasta('ID', 'AATTGTAAA'), (sequences.Fasta('ID', 'TTGTAA'), '+', 2)),
            (sequences.Fasta('ID', 'AATTGTAAAA'), (sequences.Fasta('ID', 'TTGTAA'), '+', 2)),
            (sequences.Fasta('ID', 'AATTGTAAAAA'), None),
            (sequences.Fasta('ID', 'TTACAA'), (sequences.Fasta('ID', 'TTGTAA'), '-', 0)),
            (sequences.Fasta('ID', 'ATTACAA'), (sequences.Fasta('ID', 'TTGTAA'), '-', 0)),
            (sequences.Fasta('ID', 'AATTACAA'), (sequences.Fasta('ID', 'TTGTAA'), '-', 0)),
            (sequences.Fasta('ID', 'AAATTACAA'), None),
            (sequences.Fasta('ID', 'TTACAAA'), (sequences.Fasta('ID', 'TTGTAA'), '-', 1)),
            (sequences.Fasta('ID', 'ATTACAAA'), (sequences.Fasta('ID', 'TTGTAA'), '-', 1)),
            (sequences.Fasta('ID', 'AATTACAAA'), (sequences.Fasta('ID', 'TTGTAA'), '-', 1)),
            (sequences.Fasta('ID', 'AAATTACAAA'), None),
            (sequences.Fasta('ID', 'TTACAAAA'), (sequences.Fasta('ID', 'TTGTAA'), '-', 2)),
            (sequences.Fasta('ID', 'ATTACAAAA'), (sequences.Fasta('ID', 'TTGTAA'), '-', 2)),
            (sequences.Fasta('ID', 'AATTACAAAA'), (sequences.Fasta('ID', 'TTGTAA'), '-', 2)),
            (sequences.Fasta('ID', 'AAATTACAAAA'), None),
        ]

        for seq, expected in tests:
            self.assertEqual(seq.make_into_gene(), expected)


    def test_make_into_gene_fastq(self):
        '''Test make_into_gene fastq'''
        print('sequences.genetic_code', sequences.genetic_code)
        tests = [
            (sequences.Fastq('ID', 'T', '1'), None),
            (sequences.Fastq('ID', 'TT', '12'), None),
            (sequences.Fastq('ID', 'TTT', '123'), None),
            (sequences.Fastq('ID', 'TTG', '123'), None),
            (sequences.Fastq('ID', 'TAA', '123'), None),
            (sequences.Fastq('ID', 'TTGAAATAA', '123456789'), (sequences.Fastq('ID', 'TTGAAATAA', '123456789'), '+', 0)),
            (sequences.Fastq('ID', 'TTGAAATAT', '123456789'), None),
            (sequences.Fastq('ID', 'TTGTAA', '123456'), (sequences.Fastq('ID', 'TTGTAA', '123456'), '+', 0)),
            (sequences.Fastq('ID', 'TTGTAAA', '1234567'), (sequences.Fastq('ID', 'TTGTAA', '123456'), '+', 0)),
            (sequences.Fastq('ID', 'TTGTAAAA', '12345678'), (sequences.Fastq('ID', 'TTGTAA', '123456'), '+', 0)),
            (sequences.Fastq('ID', 'TTGTAAAAA', '123456789'), None),
            (sequences.Fastq('ID', 'ATTGTAA', '1234567'), (sequences.Fastq('ID', 'TTGTAA', '234567'), '+', 1)),
            (sequences.Fastq('ID', 'ATTGTAAA', '12345678'), (sequences.Fastq('ID', 'TTGTAA', '234567'), '+', 1)),
            (sequences.Fastq('ID', 'ATTGTAAAA', '123456789'), (sequences.Fastq('ID', 'TTGTAA', '234567'), '+', 1)),
            (sequences.Fastq('ID', 'ATTGTAAAAA', '123456789A'), None),
            (sequences.Fastq('ID', 'AATTGTAA', '12345678'), (sequences.Fastq('ID', 'TTGTAA', '345678'), '+', 2)),
            (sequences.Fastq('ID', 'AATTGTAAA', '123456789'), (sequences.Fastq('ID', 'TTGTAA', '345678'), '+', 2)),
            (sequences.Fastq('ID', 'AATTGTAAAA', '123456789A'), (sequences.Fastq('ID', 'TTGTAA', '345678'), '+', 2)),
            (sequences.Fastq('ID', 'AATTGTAAAAA', '123456789AB'), None),
            (sequences.Fastq('ID', 'TTACAA', '123456'), (sequences.Fastq('ID', 'TTGTAA', '654321'), '-', 0)),
            (sequences.Fastq('ID', 'ATTACAA', '1234567'), (sequences.Fastq('ID', 'TTGTAA', '765432'), '-', 0)),
            (sequences.Fastq('ID', 'AATTACAA', '12345678'), (sequences.Fastq('ID', 'TTGTAA', '876543'), '-', 0)),
            (sequences.Fastq('ID', 'AAATTACAA', '123456789'), None),
            (sequences.Fastq('ID', 'TTACAAA', '1234567'), (sequences.Fastq('ID', 'TTGTAA', '654321'), '-', 1)),
            (sequences.Fastq('ID', 'ATTACAAA', '12345678'), (sequences.Fastq('ID', 'TTGTAA', '765432'), '-', 1)),
            (sequences.Fastq('ID', 'AATTACAAA', '123456789'), (sequences.Fastq('ID', 'TTGTAA', '876543'), '-', 1)),
            (sequences.Fastq('ID', 'AAATTACAAA', '123456789A'), None),
            (sequences.Fastq('ID', 'TTACAAAA', '12345678'), (sequences.Fastq('ID', 'TTGTAA', '654321'), '-', 2)),
            (sequences.Fastq('ID', 'ATTACAAAA', '123456789'), (sequences.Fastq('ID', 'TTGTAA', '765432'), '-', 2)),
            (sequences.Fastq('ID', 'AATTACAAAA', '123456789A'), (sequences.Fastq('ID', 'TTGTAA', '876543'), '-', 2)),
            (sequences.Fastq('ID', 'AAATTACAAAA', '123456789AB'), None),
        ]

        for seq, expected in tests:
            self.assertEqual(seq.make_into_gene(), expected)

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

    def test_add_insertions(self):
        '''Test add_insertions'''
        fa = sequences.Fasta('X', 'acgtacgtacgt')
        fa.add_insertions(skip=4, window=0, test=True)
        self.assertEqual(fa, sequences.Fasta('X', 'acgtNacgtNacgt'))

    def test_replace_bases(self):
        '''Check that bases get replaced correctly'''
        fa = sequences.Fasta('X', 'AUCGTUUACT')
        fa.replace_bases('U', 'T')
        self.assertEqual(fa, sequences.Fasta('X', 'ATCGTTTACT'))


    def test_replace_non_acgt(self):
        '''test replace_non_acgt'''
        tests = [
            ('acgtACGTnN', 'acgtACGTnN'),
            ('abc.g-T?aRC1T', 'aNcNgNTNaNCNT')
        ]

        for seq, expected in tests:
            fa = sequences.Fasta('id', seq)
            fa.replace_non_acgt()
            self.assertEqual(expected, fa.seq)


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

    def test_expand_nucleotides(self):
        '''Test expand_nucleotides'''
        tests = [
            (sequences.Fasta('1', 'A'), [sequences.Fasta('1.1', 'A')]),
            (sequences.Fasta('2', 'C'), [sequences.Fasta('2.1', 'C')]),
            (sequences.Fasta('3', 'G'), [sequences.Fasta('3.1', 'G')]),
            (sequences.Fasta('4', 'T'), [sequences.Fasta('4.1', 'T')]),
            (sequences.Fasta('6', 'R'), [sequences.Fasta('6.1', 'A'), sequences.Fasta('6.2', 'G')]),
            (sequences.Fasta('7', 'Y'), [sequences.Fasta('7.1', 'C'), sequences.Fasta('7.2', 'T')]),
            (sequences.Fasta('8', 'S'), [sequences.Fasta('8.1', 'C'), sequences.Fasta('8.2', 'G')]),
            (sequences.Fasta('9', 'W'), [sequences.Fasta('9.1', 'A'), sequences.Fasta('9.2', 'T')]),
            (sequences.Fasta('10', 'K'), [sequences.Fasta('10.1', 'G'), sequences.Fasta('10.2', 'T')]),
            (sequences.Fasta('11', 'M'), [sequences.Fasta('11.1', 'A'), sequences.Fasta('11.2', 'C')]),
            (sequences.Fasta('12', 'B'), [sequences.Fasta('12.1', 'C'), sequences.Fasta('12.2', 'G'), sequences.Fasta('12.3', 'T')]),
            (sequences.Fasta('13', 'D'), [sequences.Fasta('13.1', 'A'), sequences.Fasta('13.2', 'G'), sequences.Fasta('13.3', 'T')]),
            (sequences.Fasta('14', 'H'), [sequences.Fasta('14.1', 'A'), sequences.Fasta('14.2', 'C'), sequences.Fasta('14.3', 'T')]),
            (sequences.Fasta('15', 'V'), [sequences.Fasta('15.1', 'A'), sequences.Fasta('15.2', 'C'), sequences.Fasta('15.3', 'G')]),
            (sequences.Fasta('16', 'N'), [sequences.Fasta('16.1', 'A'), sequences.Fasta('16.2', 'C'), sequences.Fasta('16.3', 'G'), sequences.Fasta('16.4', 'T')]),
            (sequences.Fasta('17', 'ART'), [sequences.Fasta('17.1', 'AAT'), sequences.Fasta('17.2', 'AGT')]),
            (sequences.Fasta('18', 'ARRT'), [sequences.Fasta('18.1', 'AAAT'), sequences.Fasta('18.2', 'AAGT'), sequences.Fasta('18.3', 'AGAT'), sequences.Fasta('18.4', 'AGGT')]),
            (sequences.Fasta('19', 'ARTR'), [sequences.Fasta('19.1', 'AATA'), sequences.Fasta('19.2', 'AATG'), sequences.Fasta('19.3', 'AGTA'), sequences.Fasta('19.4', 'AGTG')]),
            (sequences.Fastq('20', 'ART', 'GHI'), [sequences.Fastq('20.1', 'AAT', 'GHI'), sequences.Fastq('20.2', 'AGT', 'GHI')]),
        ]

        for t in tests:
            self.assertListEqual(t[0].expand_nucleotides(), t[1])

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

    def test_gc_content(self):
        """Test GC content calculation works as expected"""
        tests = [
            (sequences.Fasta('ID', 'cgCG'), 1.0),
            (sequences.Fasta('ID', 'tTaA'), 0.0),
            (sequences.Fasta('ID', 'GCAT'), 0.5),
            (sequences.Fasta('ID', 'GCATNN'), 0.5),
            (sequences.Fasta('ID', 'GCATNNS'), 0.6),
            (sequences.Fasta('ID', 'GCATNNSK'), 0.5)
        ]
        for test, answer in tests:
            self.assertAlmostEqual(test.gc_content(), answer)
            self.assertAlmostEqual(test.gc_content(as_decimal=False), answer * 100)

class TestEmbl(unittest.TestCase):
    def test_get_id_from_header_line(self):
        '''Test get id from header line of EMBL'''
        embl = sequences.Embl('ID', 'ACGT')
        self.assertEqual(embl._get_id_from_header_line('ID   X; blah'), 'X')
        self.assertEqual(embl._get_id_from_header_line('LOCUS   X foo'), 'X')
        with self.assertRaises(sequences.Error):
            self.assertEqual(embl._get_id_from_header_line('ID X;'), 'X')
        with self.assertRaises(sequences.Error):
            self.assertEqual(embl._get_id_from_header_line('XX   X;'), 'X')


    def test_get_next_from_embl_file(self):
        f_in = utils.open_file_read(os.path.join(data_dir, 'sequences_test.embl'))
        embl = sequences.Embl()
        counter = 1

        while embl.get_next_from_file(f_in):
            self.assertEqual(embl, sequences.Fasta('seq' + str(counter), expected_embl[counter-1]))
            counter += 1

        utils.close(f_in)


    def test_get_next_from_gbk_file(self):
        f_in = utils.open_file_read(os.path.join(data_dir, 'sequences_test.gbk'))
        embl = sequences.Embl()
        counter = 1
        expected = [
            'gatcctccatatacaacggtatctccacctcaggtttagatctcaacaacggaaccattgccgacatgagacagttaggtatcgtcgagagttacaagctaaaacgagcagtagtcagctctgcatctgaagccgctgaagttctactaagggtggataacatcatccgtgcaagaccaatgccatgactcagattctaattttaagctattcaatttctctttgatc',
            'gatcctccatatacaacggtatctccacctcaggtttagatctcaacaacggaaccattgccgacatgagacagttaggtatcgtcgagagttacaagctaaaacgagcagtagtcagctctgcatctgaagccgctgaagttctactaagggtggataacatcatccgtgcaagaccaatgccatgactcagattctaattttaagctattcaatttctctttgaaa']

        while embl.get_next_from_file(f_in):
            self.assertEqual(embl, sequences.Fasta('NAME' + str(counter), expected[counter-1]))
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

    def test_subseq(self):
        '''Test subseq'''
        fq = sequences.Fastq('name', 'ACGTA', 'FGHIJ')
        self.assertEqual(fq.subseq(1,4), sequences.Fastq('name', 'CGT', 'GHI'))
        self.assertEqual(fq.subseq(None,4), sequences.Fastq('name', 'ACGT', 'FGHI'))
        self.assertEqual(fq.subseq(1,None), sequences.Fastq('name', 'CGTA', 'GHIJ'))

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
        good_files = [
            'sequences_test_gffv3.gff',
            'sequences_test_gffv3.no_FASTA_line.gff'
        ]
        good_files = [os.path.join(data_dir, x) for x in good_files]

        for f in good_files:
            reader = sequences.file_reader(f)
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
            sequences.Fasta('Turkey', 'AACTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT'),
            sequences.Fasta('Salmo_gair', 'AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT'),
            sequences.Fasta('H. Sapiens', 'ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA')
        ]

        for fname in test_files:
            reader = sequences.file_reader(fname)
            i = 0
            for seq in reader:
                self.assertEqual(expected_seqs[i], seq)
                i += 1

        # files made by seaview are a little different in the first line.
        # Test one of these
        expected_seqs = [
            sequences.Fasta('seq1', 96 * 'G' + 'T'),
            sequences.Fasta('seq2', 94 * 'A' + 'G')
        ]

        reader = sequences.file_reader(os.path.join(data_dir, 'sequences_test_phylip.made_by_seaview'))
        i = 0
        for seq in reader:
            print(seq)
            self.assertEqual(expected_seqs[i], seq)
            i += 1


class TestOther(unittest.TestCase):
    def test_orfs_from_aa_seq(self):
        '''Test _orfs_from_aa_seq()'''
        test_seqs = ['',
                     '*',
                     '**',
                     'A',
                     'A*A*A',
                     'AB**CDE*AB',
                     '*ABCDE*',
                     '**ABCDE**']

        correct_coords = [[],
                          [],
                          [],
                          [intervals.Interval(0, 0)],
                          [intervals.Interval(0, 1), intervals.Interval(2, 3),intervals.Interval(4, 4)],
                          [intervals.Interval(0, 2), intervals.Interval(4, 7), intervals.Interval(8, 9)],
                          [intervals.Interval(1, 6)],
                          [intervals.Interval(2, 7)]]

        for i in range(len(test_seqs)):
            orfs = sequences._orfs_from_aa_seq(test_seqs[i])
            self.assertListEqual(correct_coords[i], orfs)


if __name__ == '__main__':
    unittest.main()
