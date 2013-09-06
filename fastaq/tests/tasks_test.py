#!/usr/bin/env python3

import sys
import filecmp
import os
import unittest
from fastaq import tasks, sequences

modules_dir = os.path.dirname(os.path.abspath(sequences.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class Error (Exception): pass


class TestCapillaryToPairs(unittest.TestCase):
    def test_capillary_to_pairs(self):
        '''Check that capillary reads file converted to paired and unpaired'''
        tmp_prefix = 'tmp.cap_to_pairs'
        tasks.capillary_to_pairs(os.path.join(data_dir, 'sequences_test_cap_to_read_pairs.fa'), tmp_prefix)
        # sequences have been hashed, so could be in any order in
        # output files. So need to check contents of files are OK
        d_correct_paired = {}
        d_correct_unpaired = {}
        tasks.file_to_dict(os.path.join(data_dir, 'sequences_test_cap_to_read_pairs.fa.paired.gz'), d_correct_paired)
        tasks.file_to_dict(os.path.join(data_dir, 'sequences_test_cap_to_read_pairs.fa.unpaired.gz'), d_correct_unpaired)
        d_test_paired = {}
        d_test_unpaired = {}
        tasks.file_to_dict(tmp_prefix + '.paired.gz', d_test_paired)
        tasks.file_to_dict(tmp_prefix + '.unpaired.gz', d_test_unpaired)
        self.assertDictEqual(d_test_paired, d_correct_paired)
        self.assertDictEqual(d_test_unpaired, d_correct_unpaired)
        os.unlink(tmp_prefix + '.paired.gz')
        os.unlink(tmp_prefix + '.unpaired.gz')


class TestDeinterleave(unittest.TestCase):
    def test_deinterleave(self):
        '''deinterleave should deal with an interleaved file correctly'''
        tmp_1 = 'tmp.deinterleaved_1.fa'
        tmp_2 = 'tmp.deinterleaved_2.fa'
        tasks.deinterleave(os.path.join(data_dir, 'sequences_test_interleaved.fa'), tmp_1, tmp_2)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_deinterleaved_1.fa'), tmp_1))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_deinterleaved_2.fa'), tmp_2))

        tasks.deinterleave(os.path.join(data_dir, 'sequences_test_interleaved.fq'), tmp_1, tmp_2, fasta_out=True)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_deinterleaved_1.fa'), tmp_1))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_deinterleaved_2.fa'), tmp_2))

        with self.assertRaises(tasks.Error):
            tasks.deinterleave(os.path.join(data_dir, 'sequences_test_interleaved_bad.fa'), tmp_1, tmp_2)
        os.unlink(tmp_1)
        os.unlink(tmp_2)


class TestEnumerateNames(unittest.TestCase):
    def test_enumerate_names(self):
        '''Test enomereate_names works with all options'''
        outfile = 'tmp.enumerate_seqs.fa'
        rename_out = outfile + '.rename'
        tasks.enumerate_names(os.path.join(data_dir, 'sequences_test_enumerate_names.fa'), outfile)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_enumerate_names.fa.out.start.1'), outfile))
        tasks.enumerate_names(os.path.join(data_dir, 'sequences_test_enumerate_names.fa'), outfile, rename_file=rename_out)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_enumerate_names.fa.out.start.1'), outfile))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_enumerate_names.fa.out.start.1.rename_file'), rename_out))
        tasks.enumerate_names(os.path.join(data_dir, 'sequences_test_enumerate_names.fa'), outfile, start_index=2)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_enumerate_names.fa.out.start.2'), outfile))
        tasks.enumerate_names(os.path.join(data_dir, 'sequences_test_enumerate_names.fa'), outfile, keep_illumina_suffix=True)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_enumerate_names.fa.out.keep_suffix'), outfile))
        os.unlink(outfile)
        os.unlink(rename_out)


class TestExtendGaps(unittest.TestCase):
    def test_extend_gaps(self):
        '''Test that gap extension works'''
        outfile = 'tmp.gap_extend.fa'
        tasks.extend_gaps(os.path.join(data_dir, 'sequences_test_extend_gaps.fa'), outfile, trim=2)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_extend_gaps.fa.out'), outfile))
        os.unlink(outfile)


class TestFastqToMiraXml(unittest.TestCase):
    def test_fastaq_to_mira_xml(self):
        '''check that fastaq_to_mira_xml makes the correct xml file from a fastq file'''
        tmp = 'tmp.mira.xml'
        tasks.fastaq_to_mira_xml(os.path.join(data_dir, 'sequences_test_good_file.fq'), tmp)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_good_file_mira.xml'), tmp))
        os.unlink(tmp)


class TestFilter(unittest.TestCase):
    def test_length_filter(self):
        '''Check that filtering by length works as expected'''
        infile = os.path.join(data_dir, 'sequences_test_length_filter.fa')
        correct_files = [os.path.join(data_dir, 'sequences_test_length_filter.min-0.max-1.fa'),
                         os.path.join(data_dir, 'sequences_test_length_filter.min-0.max-inf.fa'),
                         os.path.join(data_dir, 'sequences_test_length_filter.min-4.max-4.fa')]
        cutoffs = [(0, 1), (0, float('inf')), (4, 4)]

        for i in range(len(cutoffs)):
            outfile = 'tmp.length_filter.fa'
            tasks.filter(infile, outfile, minlength=cutoffs[i][0], maxlength=cutoffs[i][1])
            self.assertTrue(filecmp.cmp(correct_files[i], outfile))
            os.unlink(outfile)

    def test_regex_filter(self):
        '''Check that filtering by name regex works as expected'''
        infile = os.path.join(data_dir, 'sequences_test_filter_by_regex.fa')
        correct_files = [os.path.join(data_dir, 'sequences_test_filter_by_regex.numeric.fa'),
                         os.path.join(data_dir, 'sequences_test_filter_by_regex.first-of-pair.fa'),
                         os.path.join(data_dir, 'sequences_test_filter_by_regex.first-char-a.fa')]
        regexes = ['^[0-9]+$', '/1$', '^a']

        for i in range(len(regexes)):
            outfile = 'tmp.regex_filter.fa'
            tasks.filter(infile, outfile, regex=regexes[i])
            self.assertTrue(filecmp.cmp(correct_files[i], outfile))
            os.unlink(outfile)


class TestGetSeqsFlankingGaps(unittest.TestCase):
    def test_get_seqs_flanking_gaps(self):
        outfile = 'tmp.seqs_flanking_gaps'
        tasks.get_seqs_flanking_gaps(os.path.join(data_dir, 'sequences_test_get_seqs_flanking_gaps.fa'), outfile, 3, 3)
        self.assertTrue(filecmp.cmp(outfile, os.path.join(data_dir, 'sequences_test_get_seqs_flanking_gaps.fa.out')))
        os.unlink(outfile)

class TestInterleave(unittest.TestCase):
    def test_interleave(self):
        '''Check that interleave works as expected'''
        tmp = 'tmp.interleaved.fa'
        tasks.interleave(os.path.join(data_dir, 'sequences_test_deinterleaved_1.fa'),
                         os.path.join(data_dir, 'sequences_test_deinterleaved_2.fa'),
                         tmp)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_interleaved.fa'), tmp))

        with self.assertRaises(tasks.Error):
            tasks.interleave(os.path.join(data_dir, 'sequences_test_deinterleaved_bad_1.fa'),
                             os.path.join(data_dir, 'sequences_test_deinterleaved_bad_2.fa'),
                             tmp)

        with self.assertRaises(tasks.Error):
            tasks.interleave(os.path.join(data_dir, 'sequences_test_deinterleaved_bad2_1.fa'),
                             os.path.join(data_dir, 'sequences_test_deinterleaved_bad2_2.fa'),
                             tmp)
        os.unlink(tmp)


class TestReverseComplement(unittest.TestCase):
    def test_reverse_complement(self):
        '''reverse_complement should correctly reverse complement each seq in a file'''
        tmp = 'tmp.revcomp.fa'
        tasks.reverse_complement(os.path.join(data_dir, 'sequences_test.fa'), tmp)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_revcomp.fa'), tmp))
        os.unlink(tmp)


class TestSearchForSeq(unittest.TestCase):
    def test_search_for_seq(self):
        '''Test that sequence search finds all hits'''
        tmp = 'tmp.search.fa'
        tasks.search_for_seq(os.path.join(data_dir, 'sequences_test_search_string.fa'), tmp, 'AGA')
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_search_string.fa.hits'), tmp))
        os.unlink(tmp)


class TestTranslate(unittest.TestCase):
    def test_translate(self):
        '''Test translate works in each frame'''
        tmp = 'tmp.translated.fa'
        for i in range(3):
            tasks.translate(os.path.join(data_dir, 'sequences_test_translate.fa'), tmp, frame=i)
            self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_translate.fa.frame' + str(i)), tmp))

        os.unlink(tmp)


class TestTrim(unittest.TestCase):
    def test_trim(self):
        '''trim should correctly trim each seq in a file'''
        tmp = 'tmp.trim.fq'
        tasks.trim(os.path.join(data_dir, 'sequences_test_untrimmed.fq'), tmp, 2, 1)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_trimmed.fq'), tmp))
        os.unlink(tmp)


    def test_trim_Ns_at_end(self):
        '''Test Ns at ends of sequences trimmed OK'''
        tmp = 'tmp.trim.fa'
        tasks.trim_Ns_at_end(os.path.join(data_dir, 'sequences_test_trim_Ns_at_end.fa'), tmp)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_trim_Ns_at_end.fa.trimmed'), tmp))
        os.unlink(tmp)


class TestFileToDict(unittest.TestCase):
    def test_file_to_dict(self):
        '''check file_to_dict fills dictionary correctly'''
        d_test = {}
        d = {}
        tasks.file_to_dict(os.path.join(data_dir, 'sequences_test.fa'), d_test)
        for i in range(1,5):
            d[str(i)] = sequences.Fasta(str(i),'ACGTA')

        self.assertSequenceEqual(d_test.keys(),d.keys())
        for i in range(1,5):
            key = str(i)
            self.assertEqual(d_test[key].id, d[key].id)
            self.assertEqual(d_test[key].seq, d[key].seq)


class TestLengthsFromFai(unittest.TestCase):
    def test_lengths_from_fai(self):
        '''Check lengths_from_fai gets the length of each seq OK'''
        d = {}
        lengths = {str(x):x for x in range(1,5)}
        tasks.lengths_from_fai(os.path.join(data_dir, 'sequences_test_fai_test.fa.fai'), d)
        self.assertSequenceEqual(d.keys(), lengths.keys())
        for i in d:
            self.assertEqual(int(i), d[i])


class TestSplit(unittest.TestCase):
    def test_split_by_base_count(self):
        '''Check that fasta/q files get split by base count correctly'''
        infile = os.path.join(data_dir, 'sequences_test_split_test.fa')
        outprefix = 'tmp.sequences_test_split_test.fa.test'
        length2files = {2: ['1','2','3','4'],
                        3: ['1','2','3'],
                        4: ['1', '2', '3'],
                        6: ['1', '2']}
        for l in length2files:
            tasks.split_by_base_count(infile, outprefix, l)
            for x in range(len(length2files[l])):
                file_index = str(length2files[l][x])
                fname = outprefix + '.' + file_index
                self.assertTrue(filecmp.cmp(fname, infile + '.' + str(l) + '.' + file_index))
                os.unlink(fname)

        # check that limiting the number of files works
        tasks.split_by_base_count(infile, outprefix, 6, 2)
        for i in range(1,4):
            test_file = outprefix + '.' + str(i)
            self.assertTrue(filecmp.cmp(test_file, os.path.join(data_dir, 'sequences_test_split_test.fa.6.limit2.') + str(i)))
            os.unlink(test_file)

        # check big sequence not broken
        tasks.split_by_base_count(os.path.join(data_dir, 'sequences_test_split_test.long.fa'), outprefix, 2)
        self.assertTrue(filecmp.cmp(outprefix + '.1', os.path.join(data_dir, 'sequences_test_split_test.long.fa.2.1')))
        self.assertTrue(filecmp.cmp(outprefix + '.2', os.path.join(data_dir, 'sequences_test_split_test.long.fa.2.2')))
        os.unlink(outprefix + '.1')
        os.unlink(outprefix + '.2')

    def test_split_by_fixed_size(self):
        '''Test fasta/q file split by fixed size'''
        infile = os.path.join(data_dir, 'sequences_test_split_fixed_size.fa')
        outprefix = 'tmp.sequences_test_split'
        tasks.split_by_fixed_size(infile, outprefix, 4, 1)
  
        for i in range(1,7,1):
            correct = os.path.join(data_dir, 'sequences_test_split_fixed_size.fa.split.' + str(i))
            test = outprefix + '.' + str(i)
            self.assertTrue(filecmp.cmp(test, correct))
            os.unlink(test)
 
        test_coords = outprefix + '.coords'
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_split_fixed_size.fa.split.coords'), test_coords))
        os.unlink(test_coords)



class TestCountSequences(unittest.TestCase):
    def test_count_sequences(self):
        '''Check that count_sequences does as expected'''
        self.assertEqual(2, tasks.count_sequences(os.path.join(data_dir, 'sequences_test_good_file.fq')))
        self.assertEqual(4, tasks.count_sequences(os.path.join(data_dir, 'sequences_test.fa')))
        self.assertEqual(0, tasks.count_sequences(os.path.join(data_dir, 'sequences_test_empty_file')))

class TestGetIds(unittest.TestCase):
    def test_get_ids(self):
        '''Check that IDs extracted correctly from fasta/q file'''
        tmpfile = 'tmp.ids'
        tasks.get_ids(os.path.join(data_dir, 'sequences_test.fa'), tmpfile)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test.fa.ids'), tmpfile))
        os.unlink(tmpfile)

class TestFastaToFastq(unittest.TestCase):
    def test_fasta_to_fastq(self):
        '''Check fasta_to_fastq converts files as expected'''
        tasks.fasta_to_fastq(os.path.join(data_dir, 'sequences_test.fa'),
                             os.path.join(data_dir, 'sequences_test.fa.qual'),
                             'tmp.fq')
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test.fasta_to_fastq.fq'), 'tmp.fq'))

        with self.assertRaises(tasks.Error):
            tasks.fasta_to_fastq(os.path.join(data_dir, 'sequences_test.fa'),
                                 os.path.join(data_dir, 'sequences_test.fa.qual.bad'),
                                 'tmp.fq')

        os.unlink('tmp.fq')


class TestReplaceBases(unittest.TestCase):
    def test_sequences_replace_bases(self):
        '''Check that fasta file gets all bases replaced OK'''
        tmpfile = 'tmp.replace_bases.fa'
        tasks.replace_bases(os.path.join(data_dir, 'sequences_test_fastaq_replace_bases.fa'), tmpfile, 'T', 'X')
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_fastaq_replace_bases.expected.fa'), tmpfile))
        os.unlink(tmpfile)


class TestStripIlluminaSuffix(unittest.TestCase):
    def test_strip_illumina_suffix(self):
        '''Check illumina suffixes stripped correctly off read names'''
        tmpfile = 'tmp.stripped.fa'
        tasks.strip_illumina_suffix(os.path.join(data_dir, 'sequences_test_strip_illumina_suffix.fq'), tmpfile)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_strip_illumina_suffix.fq.stripped'), tmpfile))
        os.unlink(tmpfile)
      

class TestToQuasrPrimers(unittest.TestCase):
    def test_to_quasr_primers(self):
        '''Check that fasta file gets converted to QUASR sequence file'''
        tmpfile = 'tmp.primers'
        tasks.to_quasr_primers(os.path.join(data_dir, 'sequences_test_fastaq_to_quasr_primers.fa'), tmpfile)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_fastaq_to_quasr_primers.expected'), tmpfile))
        os.unlink(tmpfile)


class TestToUniqueByID(unittest.TestCase):
    def test_to_unique_by_id(self):
        '''Test to_unique_by_id()'''
        tmpfile = 'tmp.unique_by_id.fa'
        tasks.to_unique_by_id(os.path.join(data_dir, 'sequences_test_to_unique_by_id.fa'), tmpfile)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_to_unique_by_id.fa.out'), tmpfile))
        os.unlink(tmpfile)


if __name__ == '__main__':
    unittest.main()

