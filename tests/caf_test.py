import os
import unittest

from pyfastaq import caf, utils, sequences

data_dir = "tests/data"

class TestCaf(unittest.TestCase):
    def test_get_next_from_file(self):
        '''Test get_next_from_file()'''

        f_in = utils.open_file_read(os.path.join(data_dir, 'caf_test.caf'))

        c = caf.Caf()
        c.get_next_from_file(f_in)
        read = caf.Caf()
        read.id = 'read1.p1k'
        read.seq = sequences.Fasta(read.id, 'NACGTAN')
        read.seq = read.seq.to_Fastq([4, 24, 42, 43, 40, 30, 8])
        read.insert_min = 2000
        read.insert_max = 4000
        read.ligation = '12345'
        read.clone = 'clone1'
        read.clip_start = 1
        read.clip_end = 5
        self.assertEqual(c, read)

        c.get_next_from_file(f_in)
        read = caf.Caf()
        read.id = 'read2.p1k'
        read.seq = sequences.Fasta(read.id, 'CGACGTT')
        read.seq = read.seq.to_Fastq([9, 9, 40, 41, 42, 42, 4])
        read.insert_min = 2000
        read.insert_max = 4000
        read.ligation = '23456'
        read.clone = 'clone2'
        read.clip_start = None
        read.clip_end = None
        self.assertEqual(c, read)

        utils.close(f_in)


if __name__ == '__main__':
    unittest.main()
