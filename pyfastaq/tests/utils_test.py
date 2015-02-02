#!/usr/bin/env python3

import sys
import os
import filecmp
import unittest
from pyfastaq import utils

modules_dir = os.path.dirname(os.path.abspath(utils.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestUtils(unittest.TestCase):
    def test_write_and_read(self):
        '''open_file_write() and open_file_read() should do the right thing depending gzipped or not'''
        for filename in ['utils.tmp', 'utils.tmp.gz', 'utils.tmp.bgz']:
            f = utils.open_file_write(filename)
            for i in range(3):
                print(i, file=f)
            utils.close(f)

            counter = 0

            f = utils.open_file_read(filename)
            for line in f:
                self.assertEqual(counter, int(line.strip()))
                counter += 1
            utils.close(f)

            os.unlink(filename)

        f = utils.open_file_read('-')
        self.assertEqual(sys.stdin, f)
        f = utils.open_file_write('-')
        self.assertEqual(sys.stdout, f)

    def test_raise_exception(self):
        '''open_file_write() and open_file_read() should raise an exception when can't do the opening'''
        with self.assertRaises(utils.Error):
            utils.open_file_read('this_file_is_not_here_so_throw_error')
        with self.assertRaises(utils.Error):
            utils.open_file_read('this_file_is_not_here_so_throw_error.gz')
        with self.assertRaises(utils.Error):
            utils.open_file_read(os.path.join(data_dir, 'utils_test_not_really_zipped.gz'))

        with self.assertRaises(utils.Error):
            utils.open_file_write(os.path.join('not_a_directory', 'this_file_is_not_here_so_throw_error'))
        with self.assertRaises(utils.Error):
            utils.open_file_write(os.path.join('not_a_directory', 'this_file_is_not_here_so_throw_error.gz'))

    def test_file_transpose(self):
        '''Test that file_transpose() does what it should'''
        infile = os.path.join(data_dir, 'utils_test_file_transpose.txt')
        tmp_out = 'utils_test_file_transpose.tmp'
        correct_file = os.path.join(data_dir, 'utils_test_file_transposed.txt')
        utils.file_transpose(infile, tmp_out)
        self.assertTrue(filecmp.cmp(tmp_out, correct_file))
        os.unlink(tmp_out)

    def test_system_call(self):
        '''Test that system call appears to work and die as it should'''
        test_file = os.path.join(data_dir, 'utils_test_system_call.txt')
        tmp_out = 'utils_test_syscall.tmp'
        utils.syscall('cat ' + test_file + ' > ' + tmp_out)
        self.assertTrue(filecmp.cmp(tmp_out, test_file))
        os.unlink(tmp_out)

        with self.assertRaises(utils.Error):
            utils.syscall('thisisveryunlikelytoebarealcommandandshouldthrowerror')

        utils.syscall('echo "this is not the right string" > ' + tmp_out)
        self.assertFalse(filecmp.cmp(tmp_out, test_file))
        os.unlink(tmp_out)

        s = utils.syscall_get_stdout('echo bingo')
        self.assertListEqual(["bingo"], s)
        with self.assertRaises(utils.Error):
            utils.syscall_get_stdout('thisisveryunlikelytoebarealcommandandshouldthrowerror')

if __name__ == '__main__':
    unittest.main()
