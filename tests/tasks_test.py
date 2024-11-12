import sys
import filecmp
import os
import tempfile
import unittest

from pyfastaq import tasks, sequences

data_dir = "tests/data"

class Error (Exception): pass


class TestACGTN_only(unittest.TestCase):
    def test_acgtn_only(self):
        '''Test acgtn_only'''
        tmpfile = 'tmp.test_acgtn_only.fa'
        infile = os.path.join(data_dir, 'test_acgtn_only.in.fa')
        expected = os.path.join(data_dir, 'test_acgtn_only.expected.fa')
        tasks.acgtn_only(infile, tmpfile)
        self.assertTrue(filecmp.cmp(expected, tmpfile, shallow=False))
        os.unlink(tmpfile)


class TestCafToFastq(unittest.TestCase):
    def test_caf_to_fastq_default(self):
        '''Test caf_to_fastq with no filtering'''
        tmpfile = 'tmp.fq'
        tasks.caf_to_fastq(os.path.join(data_dir, 'caf_test.caf'), tmpfile)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'caf_test.to_fastq.no_trim.min_length_0.fq'), tmpfile, shallow=False))
        os.unlink(tmpfile)

    def test_caf_to_fastq_trim_and_min_length(self):
        '''Test caf_to_fastq with trimming and min_length'''
        tmpfile = 'tmp.fq'
        tasks.caf_to_fastq(os.path.join(data_dir, 'caf_test.caf'), tmpfile, trim=True, min_length=6)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'caf_test.to_fastq.trim.min_length_6.fq'), tmpfile, shallow=False))
        os.unlink(tmpfile)


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

        tasks.enumerate_names(os.path.join(data_dir, 'sequences_test_enumerate_names.fa'), outfile, suffix='.SUFFIX')
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_enumerate_names.fa.out.add_suffix'), outfile, shallow=False))
        os.unlink(outfile)
        os.unlink(rename_out)


class TestExpandNucleotides(unittest.TestCase):
    def test_expand_nucleoties(self):
        '''Test expand_nucleotides'''
        tmp = 'tmp.expanded'
        fq_in = os.path.join(data_dir, 'tasks_test_expend_nucleotides.in.fq')
        fa_in = os.path.join(data_dir, 'tasks_test_expend_nucleotides.in.fa')
        fq_expected = os.path.join(data_dir, 'tasks_test_expend_nucleotides.out.fq')
        fa_expected = os.path.join(data_dir, 'tasks_test_expend_nucleotides.out.fa')
        tasks.expand_nucleotides(fq_in, tmp)
        self.assertTrue(filecmp.cmp(fq_expected, tmp, shallow=False))
        os.unlink(tmp)
        tasks.expand_nucleotides(fa_in, tmp)
        self.assertTrue(filecmp.cmp(fa_expected, tmp, shallow=False))
        os.unlink(tmp)


class TestExtendGaps(unittest.TestCase):
    def test_trim_contigs(self):
        '''Test that gap extension works'''
        outfile = 'tmp.gap_extend.fa'
        tasks.trim_contigs(os.path.join(data_dir, 'sequences_test_trim_contigs.fa'), outfile, trim=2)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_trim_contigs.fa.out'), outfile))
        os.unlink(outfile)


class TestFastqToMiraXml(unittest.TestCase):
    def test_fastaq_to_mira_xml(self):
        '''check that fastaq_to_mira_xml makes the correct xml file from a fastq file'''
        tmp = 'tmp.mira.xml'
        tasks.fastaq_to_mira_xml(os.path.join(data_dir, 'sequences_test_good_file.fq'), tmp)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_good_file_mira.xml'), tmp))
        os.unlink(tmp)


class TestFastaqToOrfsGFF(unittest.TestCase):
    def test_fastaq_to_orfs_gff(self):
        '''Test fastaq_to_orfs_gff'''
        outfile = 'tmp.orfs.gff'
        tasks.fastaq_to_orfs_gff(os.path.join(data_dir, 'sequences_test_orfs.fa'), outfile, min_length=120)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_orfs.gff'), outfile, shallow=False))
        os.unlink(outfile)


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

    def test_regex_check_comments_filter(self):
        '''When check_comments is true, and the regex is in the comment'''
        infile = tempfile.NamedTemporaryFile(suffix=".fa", mode="w+")
        infile.write(
            ">read1 foo=bar\nAGCT\n>read2 bar=foo\nGGG\n>read3\nGGGG\n>read4 foo=ba\n"
            "GCA\n>read5foo=bar\nGCAT"
        )
        infile.seek(0)
        regex = '\sfoo=bar'
        outfile = tempfile.NamedTemporaryFile(suffix=".fa", mode="w+")

        tasks.filter(infile.name, outfile.name, regex=regex, check_comments=True)
        with open(outfile.name) as handle:
            actual = handle.read()

        expected = ">read1 foo=bar\nAGCT\n"

        self.assertEqual(actual, expected)

    def test_ids_from_file_filter(self):
        '''Test that can extract reads from a file of read names'''
        infile = os.path.join(data_dir, 'sequences_test_filter_by_ids_file.fa')
        outfile = 'tmp.ids_file_filter.fa'
        tasks.filter(infile, outfile, ids_file=infile + '.ids')
        self.assertTrue(filecmp.cmp(infile + '.filtered', outfile))
        os.unlink(outfile)
		
    def test_ids_with_comments_from_file_filter(self):
        '''Test that can extract reads from a file of read names where the read names have extra data after space'''
        infile = os.path.join(data_dir, 'readnames_with_comments.fastq')
        outfile = 'tmp.ids_file_filter.fastq'
        tasks.filter(infile, outfile, ids_file=infile + '.ids')
        self.assertTrue(filecmp.cmp(infile + '.filtered', outfile))
        os.unlink(outfile)

    def test_invert_filter(self):
        '''Test that inverting filtering works'''
        infile = os.path.join(data_dir, 'sequences_test_filter_by_ids_file.fa')
        outfile = 'tmp.ids_file_filter.fa'
        tasks.filter(infile, outfile, ids_file=infile + '.ids', invert=True)
        self.assertTrue(filecmp.cmp(infile + '.filtered.invert', outfile))
        os.unlink(outfile)


    def test_paired_both_pass(self):
        '''Test filter with paired file both pass'''
        infile1 = os.path.join(data_dir, 'tasks_test_filter_paired_both_pass.in_1.fa')
        infile2 = os.path.join(data_dir, 'tasks_test_filter_paired_both_pass.in_2.fa')
        outfile1 = 'tmp.filter_both_pass_1.fa'
        outfile2 = 'tmp.filter_both_pass_2.fa'
        expected1 = os.path.join(data_dir, 'tasks_test_filter_paired_both_pass.out_1.fa')
        expected2 = os.path.join(data_dir, 'tasks_test_filter_paired_both_pass.out_2.fa')
        tasks.filter(infile1, outfile1, mate_in=infile2, mate_out=outfile2, minlength=3)
        self.assertTrue(filecmp.cmp(outfile1, expected1, shallow=False))
        self.assertTrue(filecmp.cmp(outfile2, expected2, shallow=False))
        os.unlink(outfile1)
        os.unlink(outfile2)


    def test_paired_one_pass(self):
        '''Test filter with paired file one pass'''
        infile1 = os.path.join(data_dir, 'tasks_test_filter_paired_one_pass.in_1.fa')
        infile2 = os.path.join(data_dir, 'tasks_test_filter_paired_one_pass.in_2.fa')
        outfile1 = 'tmp.filter_one_pass_1.fa'
        outfile2 = 'tmp.filter_one_pass_2.fa'
        expected1 = os.path.join(data_dir, 'tasks_test_filter_paired_one_pass.out_1.fa')
        expected2 = os.path.join(data_dir, 'tasks_test_filter_paired_one_pass.out_2.fa')
        tasks.filter(infile1, outfile1, mate_in=infile2, mate_out=outfile2, both_mates_pass=False, minlength=3)
        self.assertTrue(filecmp.cmp(outfile1, expected1, shallow=False))
        self.assertTrue(filecmp.cmp(outfile2, expected2, shallow=False))
        os.unlink(outfile1)
        os.unlink(outfile2)


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

        tasks.interleave(os.path.join(data_dir, 'sequences_test_deinterleaved_no_suffixes_1.fa'),
                         os.path.join(data_dir, 'sequences_test_deinterleaved_no_suffixes_2.fa'),
                         tmp, suffix1='/1', suffix2='/2')
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_interleaved_with_suffixes.fa'), tmp))

        with self.assertRaises(tasks.Error):
            tasks.interleave(os.path.join(data_dir, 'sequences_test_deinterleaved_bad_1.fa'),
                             os.path.join(data_dir, 'sequences_test_deinterleaved_bad_2.fa'),
                             tmp)

        with self.assertRaises(tasks.Error):
            tasks.interleave(os.path.join(data_dir, 'sequences_test_deinterleaved_bad2_1.fa'),
                             os.path.join(data_dir, 'sequences_test_deinterleaved_bad2_2.fa'),
                             tmp)
        os.unlink(tmp)


class TestMakeRandomContigs(unittest.TestCase):
    def test_make_random_contigs(self):
        '''Test make_random_contigs()'''
        # Can't guarantee same results from random (even using same seed), so
        # just check sequence names and lengths
        def files_are_equal(file1, file2):
            seqs1 = {}
            seqs2 = {}
            tasks.file_to_dict(file1, seqs1)
            tasks.file_to_dict(file2, seqs2)
            if len(seqs1) != len(seqs2):
                return False

            for name in seqs1:
                seq1 = seqs1[name]
                seq2 = seqs2[name]
                if seq1.id != seq2.id:
                    return False
                if len(seq1) != len(seq2):
                    return False

            return True

        tmp = 'tmp.random_contigs.fa'
        tasks.make_random_contigs(2, 3, tmp)
        self.assertTrue(files_are_equal(os.path.join(data_dir, 'sequences_test_make_random_contigs.default.fa'), tmp))
        tasks.make_random_contigs(2, 3, tmp, prefix='p')
        self.assertTrue(files_are_equal(os.path.join(data_dir, 'sequences_test_make_random_contigs.prefix-p.fa'), tmp))
        tasks.make_random_contigs(2, 3, tmp, first_number=42)
        self.assertTrue(files_are_equal(os.path.join(data_dir, 'sequences_test_make_random_contigs.first-42.fa'), tmp))
        tasks.make_random_contigs(28, 3, tmp, name_by_letters=True)
        self.assertTrue(files_are_equal(os.path.join(data_dir, 'sequences_test_make_random_contigs.name-by-letters.fa'), tmp))
        os.unlink(tmp)


class TestMeanLength(unittest.TestCase):
    def test_mean_length(self):
        '''Test mean_length'''
        expected = [3, 2, 3, 4, 4]
        limits = [1, 2, 3, 4, None]
        assert len(expected) == len(limits)
        for i in range(len(expected)):
            mean = tasks.mean_length(os.path.join(data_dir, 'tasks_test_mean_length.fa'), limit=limits[i])
            self.assertEqual(expected[i], mean)


class TestMergeToOneSeq(unittest.TestCase):
    def test_merge_to_one_seq_fa(self):
        '''Test merge_to_one_seq with fasta'''
        tmp = 'tmp.merged.fa'
        tasks.merge_to_one_seq(os.path.join(data_dir, 'sequences_test_merge_to_one_seq.fa'), tmp)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_merge_to_one_seq.merged.fa'), tmp, shallow=False))
        os.unlink(tmp)

    def test_merge_to_one_seq_fq(self):
        '''Test merge_to_one_seq with fastq'''
        tmp = 'tmp.merged.fq'
        tasks.merge_to_one_seq(os.path.join(data_dir, 'sequences_test_merge_to_one_seq.fq'), tmp)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_merge_to_one_seq.merged.fq'), tmp, shallow=False))
        os.unlink(tmp)

class TestReverseComplement(unittest.TestCase):
    def test_reverse_complement(self):
        '''reverse_complement should correctly reverse complement each seq in a file'''
        tmp = 'tmp.revcomp.fa'
        tasks.reverse_complement(os.path.join(data_dir, 'sequences_test.fa'), tmp)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_revcomp.fa'), tmp))
        os.unlink(tmp)


class TestScaffoldsToContigs(unittest.TestCase):
    def test_scaffolds_to_contigs(self):
        '''Test scaffolds_to_contigs'''
        tmp = 'tmp.contigs.fa'
        tasks.scaffolds_to_contigs(os.path.join(data_dir, 'utils_test_scaffolds.fa'), tmp)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'utils_test_scaffolds.fa.to_contigs.fa'), tmp))
        os.unlink(tmp)

    def test_scaffolds_to_contigs_number_contigs(self):
        '''Test scaffolds_to_contigs with contig numbering'''
        tmp = 'tmp.contigs.fa'
        tasks.scaffolds_to_contigs(os.path.join(data_dir, 'utils_test_scaffolds.fa'), tmp, number_contigs=True)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'utils_test_scaffolds.fa.to_contigs.number_contigs.fa'), tmp))
        os.unlink(tmp)


class TestSearchForSeq(unittest.TestCase):
    def test_search_for_seq(self):
        '''Test that sequence search finds all hits'''
        tmp = 'tmp.search.fa'
        tasks.search_for_seq(os.path.join(data_dir, 'sequences_test_search_string.fa'), tmp, 'AGA')
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_search_string.fa.hits'), tmp))
        os.unlink(tmp)


class TestSequenceTrim(unittest.TestCase):
    def test_sequence_trim(self):
        '''Test sequence_trim'''
        tmp1 = 'tmp.trimmed_1.fa'
        tmp2 = 'tmp.trimmed_2.fa'
        in1 = os.path.join(data_dir, 'tasks_test_sequence_trim_1.fa')
        in2 = os.path.join(data_dir, 'tasks_test_sequence_trim_2.fa')
        to_trim = os.path.join(data_dir, 'tasks_test_sequences_to_trim.fa')
        expected1 = os.path.join(data_dir, 'tasks_test_sequence_trim_1.trimmed.fa')
        expected2 = os.path.join(data_dir, 'tasks_test_sequence_trim_2.trimmed.fa')
        tasks.sequence_trim(in1, in2, tmp1, tmp2, to_trim, min_length=10, check_revcomp=True)
        self.assertTrue(filecmp.cmp(expected1, tmp1))
        self.assertTrue(filecmp.cmp(expected2, tmp2))
        os.unlink(tmp1)
        os.unlink(tmp2)


class ToFastg(unittest.TestCase):
    def test_to_fastg_ids_set(self):
        '''Test to_fastg when ids are a set'''
        infile = os.path.join(data_dir, 'tasks_test_to_fastg.fasta')
        tmpfile = 'tmp.to_fastg.fastg'
        expected = os.path.join(data_dir, 'tasks_test_to_fastg.fastg')
        ids = {'seq2'}
        tasks.to_fastg(infile, tmpfile, circular=ids)
        self.assertTrue(filecmp.cmp(expected, tmpfile, shallow=False))
        os.unlink(tmpfile)


    def test_to_fastg_ids_file(self):
        '''Test to_fastg when ids in a file'''
        infile = os.path.join(data_dir, 'tasks_test_to_fastg.fasta')
        tmpfile = 'tmp.to_fastg.fastg'
        expected = os.path.join(data_dir, 'tasks_test_to_fastg.fastg')
        ids_file = os.path.join(data_dir, 'tasks_test_to_fastg.ids_to_circularise')
        tasks.to_fastg(infile, tmpfile, circular=ids_file)
        self.assertTrue(filecmp.cmp(expected, tmpfile, shallow=False))
        os.unlink(tmpfile)


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


class TestLengthOffsetsFromFai(unittest.TestCase):
    def test_length_offsets_from_fai(self):
        '''Test length_offsets_from_fai'''
        got = tasks.length_offsets_from_fai(os.path.join(data_dir, 'tasks_test_length_offsets_from_fai.fa.fai'))
        expected = {'seq1': 0, 'seq2': 42, 'seq3': 43}
        self.assertEqual(expected, got)


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

    def test_split_by_fixed_size_exclude_Ns(self):
        infile = os.path.join(data_dir, 'sequences_test_split_fixed_size.fa')
        outprefix = 'tmp.sequences_test_split'
        tasks.split_by_fixed_size(infile, outprefix, 4, 1, skip_if_all_Ns=True)

        for i in range(1,5,1):
            correct = os.path.join(data_dir, 'sequences_test_split_fixed_size.fa.split.skip_if_all_Ns.' + str(i))
            test = outprefix + '.' + str(i)
            self.assertTrue(filecmp.cmp(test, correct))
            os.unlink(test)

        test_coords = outprefix + '.coords'
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_split_fixed_size.fa.split.skip_if_all_Ns.coords'), test_coords))
        os.unlink(test_coords)

    def test_split_by_fixed_size_onefile(self):
        infile = os.path.join(data_dir, 'sequences_test_split_fixed_size_onefile.fa')
        tmp_out = 'tmp.sequences_test_split_fixed_size_onefile.fa'
        expected =  os.path.join(data_dir, 'sequences_test_split_fixed_size_onefile.out.fa')
        tasks.split_by_fixed_size_onefile(infile, tmp_out, chunk_size=3, tolerance=1)
        self.assertTrue(filecmp.cmp(expected, tmp_out))
        os.unlink(tmp_out)

    def test_split_by_fixed_size_onefile_exclude_Ns(self):
        infile = os.path.join(data_dir, 'sequences_test_split_fixed_size_onefile.fa')
        tmp_out = 'tmp.sequences_test_split_fixed_size_onefile.skip_Ns.fa'
        expected =  os.path.join(data_dir, 'sequences_test_split_fixed_size_onefile.skip_Ns.out.fa')
        tasks.split_by_fixed_size_onefile(infile, tmp_out, chunk_size=3, tolerance=1, skip_if_all_Ns=True)
        self.assertTrue(filecmp.cmp(expected, tmp_out))
        os.unlink(tmp_out)

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


class TestFastaToFakeQual(unittest.TestCase):
    def test_fasta_to_fake_qual(self):
        '''Test fasta_to_fake_qual'''
        tmpfile = 'tmp.qual'
        infile = os.path.join(data_dir, 'tasks_test_fasta_to_fake_qual.in.fa')
        tasks.fastaq_to_fake_qual(infile, tmpfile)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'tasks_test_fasta_to_fake_qual.out.default.qual'), tmpfile, shallow=False))
        os.unlink(tmpfile)
        tasks.fastaq_to_fake_qual(infile, tmpfile, q=42)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'tasks_test_fasta_to_fake_qual.out.q42.qual'), tmpfile, shallow=False))
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



class TestSortBySize(unittest.TestCase):
    def test_sort_by_size(self):
        '''Test sort_by_size'''
        infile = os.path.join(data_dir, 'tasks_test_sort_by_size.in.fa')
        tmpfile = 'tmp.sorted.fa'
        tasks.sort_by_size(infile, tmpfile)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'tasks_test_sort_by_size.out.fa'), tmpfile, shallow=False))
        tasks.sort_by_size(infile, tmpfile, smallest_first=True)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'tasks_test_sort_by_size.out.rev.fa'), tmpfile, shallow=False))
        os.unlink(tmpfile)

class TestSortByName(unittest.TestCase):
    def test_sort_by_name(self):
        '''Test sort_by_name'''
        infile = os.path.join(data_dir, 'tasks_test_sort_by_name.in.fa')
        tmpfile = 'tmp.sort_by_name.fa'
        tasks.sort_by_name(infile, tmpfile)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'tasks_test_sort_by_name.out.fa'), tmpfile, shallow=False))
        os.unlink(tmpfile)


class TestStripIlluminaSuffix(unittest.TestCase):
    def test_strip_illumina_suffix(self):
        '''Check illumina suffixes stripped correctly off read names'''
        tmpfile = 'tmp.stripped.fa'
        tasks.strip_illumina_suffix(os.path.join(data_dir, 'sequences_test_strip_illumina_suffix.fq'), tmpfile)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_strip_illumina_suffix.fq.stripped'), tmpfile))
        os.unlink(tmpfile)


class TestStatsFromFai(unittest.TestCase):
    def test_stats_from_fai_nonempty(self):
        '''Test task stats_from_fai non-empty file'''
        infile = os.path.join(data_dir, 'tasks_test_stats_from_fai.in.fai')
        got = tasks.stats_from_fai(infile)
        expected = {
            'longest': 10,
            'shortest': 1,
            'N50': 4,
            'mean': 4.2,
            'number': 5,
            'total_length': 21
        }
        self.assertEqual(expected, got)


    def test_stats_from_fai_empty(self):
        '''Test task stats_from_fai empty file'''
        infile = os.path.join(data_dir, 'tasks_test_stats_from_fai.in.empty.fai')
        got = tasks.stats_from_fai(infile)
        expected = {
            'longest': 0,
            'shortest': 0,
            'N50': 0,
            'mean': 0,
            'number': 0,
            'total_length': 0
        }
        self.assertEqual(expected, got)


class TestToBoulderio(unittest.TestCase):
    def test_to_boulderio(self):
        '''Test task to_boulderio'''
        tmpfile = 'tmp.boulder'
        tasks.to_boulderio(os.path.join(data_dir, 'tasks_test_to_boulderio.in.fa'), tmpfile)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'tasks_test_to_boulderio.out.boulder'), tmpfile, shallow=False))
        os.unlink(tmpfile)


class TestToFasta(unittest.TestCase):
    def test_to_fasta(self):
        '''Test to_fasta'''
        tmpfile = 'tmp.to_fasta'
        infiles = [
            'sequences_test_good_file.fq',
            'sequences_test_gffv3.gff',
            'sequences_test_gffv3.no_FASTA_line.gff',
            'sequences_test.embl',
            'sequences_test.gbk',
            'sequences_test_phylip.interleaved',
            'sequences_test_phylip.interleaved2',
            'sequences_test_phylip.sequential'
        ]
        infiles = [os.path.join(data_dir, x) for x in infiles]
        expected_outfiles = [x + '.to_fasta' for x in infiles]

        for i in range(len(infiles)):
            tasks.to_fasta(infiles[i], tmpfile)
            self.assertTrue(filecmp.cmp(expected_outfiles[i], tmpfile))

        tasks.to_fasta(os.path.join(data_dir, 'sequences_test.fa'), tmpfile, line_length=3)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test.line_length3.fa'), tmpfile))
        tasks.to_fasta(os.path.join(data_dir, 'sequences_test_strip_after_whitespace.fa'), tmpfile, strip_after_first_whitespace=True)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_strip_after_whitespace.fa.to_fasta'), tmpfile))
        os.unlink(tmpfile)

    def test_to_fasta_strip_after_whitespace_non_unique(self):
        '''Test strip_after_whitespace with non-unique names'''
        tmpfile = 'tmp.strip_after_whitespace.fa'
        infile = os.path.join(data_dir, 'sequences_test.to_fasta.strip_after_whitespace_non_unique.in.fa')
        expected = os.path.join(data_dir, 'sequences_test.to_fasta.strip_after_whitespace_non_unique.out.fa')

        with self.assertRaises(tasks.Error):
            tasks.to_fasta(infile, tmpfile, strip_after_first_whitespace=True, check_unique=True)

        tasks.to_fasta(infile, tmpfile, strip_after_first_whitespace=True, check_unique=False)
        self.assertTrue(filecmp.cmp(tmpfile, expected, shallow=False))
        os.unlink(tmpfile)

    def test_to_fasta_strip_after_whitespace_unique(self):
        '''Test strip_after_whitespace with unique names'''
        tmpfile = 'tmp.strip_after_whitespace.fa'
        infile = os.path.join(data_dir, 'sequences_test.to_fasta.strip_after_whitespace_unique.in.fa')
        expected = os.path.join(data_dir, 'sequences_test.to_fasta.strip_after_whitespace_unique.out.fa')
        tasks.to_fasta(infile, tmpfile, strip_after_first_whitespace=True, check_unique=True)
        self.assertTrue(filecmp.cmp(tmpfile, expected, shallow=False))
        os.unlink(tmpfile)

class TestToUniqueByID(unittest.TestCase):
    def test_to_unique_by_id(self):
        '''Test to_unique_by_id()'''
        tmpfile = 'tmp.unique_by_id.fa'
        tasks.to_unique_by_id(os.path.join(data_dir, 'sequences_test_to_unique_by_id.fa'), tmpfile)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_to_unique_by_id.fa.out'), tmpfile))
        os.unlink(tmpfile)


class TestToFastaUnion(unittest.TestCase):
    def test_to_fasta_union(self):
        '''Test to_fasta_union'''
        tmpfile = 'tmp.to_fasta_union'
        tasks.to_fasta_union(os.path.join(data_dir, 'sequences_test_to_fasta_union.in.fa'), tmpfile, seqname='testname')
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sequences_test_to_fasta_union.out.fa'), tmpfile, shallow=False))
        os.unlink(tmpfile)


if __name__ == '__main__':
    unittest.main()


