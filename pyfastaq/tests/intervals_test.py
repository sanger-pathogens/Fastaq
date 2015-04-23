#!/usr/bin/env python3

import unittest
from pyfastaq import intervals

class TestIntervals(unittest.TestCase):
    def test_init(self):
        '''Throw error if try to construct genome_interval from a non-int, or end<start'''
        with self.assertRaises(intervals.Error):
            intervals.Interval('a', 1)
        with self.assertRaises(intervals.Error):
            intervals.Interval(1, 'a')
        with self.assertRaises(intervals.Error):
            intervals.Interval('a', 'a')
        with self.assertRaises(intervals.Error):
            intervals.Interval(3, 2)

    def test_comparisons(self):
        '''<, <=, == should work as expected'''
        self.assertTrue(intervals.Interval(1,2) < intervals.Interval(2,2))
        self.assertTrue(intervals.Interval(1,2) <= intervals.Interval(2,2))
        self.assertFalse(intervals.Interval(2,2) <= intervals.Interval(1,2))
        self.assertFalse(intervals.Interval(2,2) < intervals.Interval(1,2))
        self.assertFalse(intervals.Interval(2,2) < intervals.Interval(2,2))
        self.assertTrue(intervals.Interval(1,2) == intervals.Interval(1,2))
        self.assertFalse(intervals.Interval(1,2) == intervals.Interval(1,3))
        self.assertTrue(intervals.Interval(1,2) != intervals.Interval(1,3))
        self.assertFalse(intervals.Interval(1,2) != intervals.Interval(1,2))

    def test_len(self):
        self.assertEqual(len(intervals.Interval(1,2)), 2)
        self.assertEqual(len(intervals.Interval(1,1)), 1)
        self.assertEqual(len(intervals.Interval(10,20)), 11)

    def test_distance_to_point(self):
        '''Test distance_to_point'''
        self.assertEqual(0, intervals.Interval(42, 50).distance_to_point(42))
        self.assertEqual(0, intervals.Interval(42, 50).distance_to_point(44))
        self.assertEqual(0, intervals.Interval(42, 50).distance_to_point(50))
        self.assertEqual(1, intervals.Interval(42, 50).distance_to_point(41))
        self.assertEqual(1, intervals.Interval(42, 50).distance_to_point(51))
        self.assertEqual(5, intervals.Interval(42, 50).distance_to_point(55))
        self.assertEqual(5, intervals.Interval(42, 50).distance_to_point(37))

    def test_intersects(self):
        '''Intersection of two intervals should do the right thing'''
        a = intervals.Interval(5, 10)
        no_intersect = [intervals.Interval(3, 4),
                        intervals.Interval(11,20)]
        intersect = [intervals.Interval(3,5),
                     intervals.Interval(3,6),
                     intervals.Interval(9,12),
                     intervals.Interval(10,12),
                     intervals.Interval(6,7),
                     intervals.Interval(1,20)]

        for i in no_intersect:
            self.assertFalse(a.intersects(i), 'shouldn\'t intersect: ' + str(a) + ', ' + str(i))

        for i in intersect:
            self.assertTrue(a.intersects(i), 'should intersect: ' + str(a) + ', ' + str(i))

    def test_contains(self):
        '''Check that contains() works as expected'''
        a = intervals.Interval(5, 10)
        not_contained = [intervals.Interval(1,2),
                         intervals.Interval(4,5),
                         intervals.Interval(4,10),
                         intervals.Interval(4,11),
                         intervals.Interval(5,11),
                         intervals.Interval(1,2),
                         intervals.Interval(9,11),
                         intervals.Interval(10,11),
                         intervals.Interval(11,20)]


        contained = [intervals.Interval(5,5),
                     intervals.Interval(5,10),
                     intervals.Interval(6,7),
                     intervals.Interval(6,10),
                     intervals.Interval(10,10)]

        for i in not_contained:
            self.assertFalse(a.contains(i), 'shouldn\'t contain: ' + str(a) + ', ' + str(i))

        for i in contained:
            self.assertTrue(a.contains(i), 'should contain: ' + str(a) + ', ' + str(i))

    def test_union(self):
        '''Union should either return None or the correct union'''
        a = intervals.Interval(5, 10)
        b = intervals.Interval(8, 15)
        c = intervals.Interval(12, 20)
        d = intervals.Interval(21,22)
        self.assertEqual(a.union(c), None)
        self.assertEqual(c.union(a), None)
        self.assertEqual(a.union(b), intervals.Interval(5,15))
        self.assertEqual(b.union(a), intervals.Interval(5,15))
        self.assertEqual(c.union(d), intervals.Interval(12,22))
        self.assertEqual(d.union(c), intervals.Interval(12,22))

    def test_union_flll_gap(self):
        '''union_fill_gap() should ignore intersections and return the maximum range of coords'''
        a = intervals.Interval(5, 10)
        b = intervals.Interval(8, 15)
        c = intervals.Interval(12, 20)
        d = intervals.Interval(21,22)
        self.assertEqual(a.union_fill_gap(c), intervals.Interval(5,20))
        self.assertEqual(c.union_fill_gap(a), intervals.Interval(5,20))
        self.assertEqual(a.union_fill_gap(b), intervals.Interval(5,15))
        self.assertEqual(b.union_fill_gap(a), intervals.Interval(5,15))
        self.assertEqual(c.union_fill_gap(d), intervals.Interval(12,22))
        self.assertEqual(d.union_fill_gap(c), intervals.Interval(12,22))


    def test_intersection(self):
        '''Intersection should either return None or the correct intersection'''
        a = intervals.Interval(5, 10)
        b = intervals.Interval(8, 15)
        c = intervals.Interval(12, 20)
        self.assertEqual(a.intersection(c), None)
        self.assertEqual(a.intersection(b), intervals.Interval(8,10))

class Test_intersection(unittest.TestCase):
    def test_intersection(self):
        '''intersection() should correctly intersect two lists of intervals'''
        a = [intervals.Interval(1,2),
             intervals.Interval(10,20),
             intervals.Interval(51,52),
             intervals.Interval(54,55),
             intervals.Interval(57,58)]

        b = [intervals.Interval(5,6),
             intervals.Interval(9,11),
             intervals.Interval(13,14),
             intervals.Interval(17,18),
             intervals.Interval(20,25),
             intervals.Interval(50,60)]

        c = [intervals.Interval(100,200)]

        i = [intervals.Interval(10,11),
             intervals.Interval(13,14),
             intervals.Interval(17,18),
             intervals.Interval(20,20),
             intervals.Interval(51,52),
             intervals.Interval(54,55),
             intervals.Interval(57,58)]

        self.assertSequenceEqual(intervals.intersection(a,b), i)
        self.assertSequenceEqual(intervals.intersection(b,a), i)
        self.assertSequenceEqual(intervals.intersection(c,a), [])
        self.assertEqual(intervals.intersection([],a), [])
        self.assertEqual(intervals.intersection(a,[]), [])

class Test_merge_overlapping_in_list(unittest.TestCase):
    def test_merge_overlapping_in_list(self):
        '''merge_overlapping_in_list() merges correctly'''
        a = [intervals.Interval(1,2),
             intervals.Interval(51,60),
             intervals.Interval(10,20),
             intervals.Interval(20,30),
             intervals.Interval(20,30),
             intervals.Interval(29,50),
             intervals.Interval(65,70)]

        b = [intervals.Interval(1,2),
             intervals.Interval(10,60),
             intervals.Interval(65,70)]

        intervals.merge_overlapping_in_list(a)
        self.assertSequenceEqual(a, b)

class Test_remove_contained_in_list(unittest.TestCase):
    def test_remove_contained_in_list(self):
        '''test_remove_contained_in_list removes the right elements of list'''
        a = [intervals.Interval(1,2),
             intervals.Interval(4,4),
             intervals.Interval(4,5),
             intervals.Interval(5,6),
             intervals.Interval(7,9),
             intervals.Interval(8,10),
             intervals.Interval(9,11),
             intervals.Interval(20,25),
             intervals.Interval(20,24),
             intervals.Interval(20,26),
             intervals.Interval(30,38),
             intervals.Interval(30,37),
             intervals.Interval(30,36),
             intervals.Interval(30,35),
             intervals.Interval(30,35),
             intervals.Interval(32,33),
             intervals.Interval(38,50),
             intervals.Interval(65,70),
             intervals.Interval(67,70)]

        b = [intervals.Interval(1,2),
             intervals.Interval(4,5),
             intervals.Interval(5,6),
             intervals.Interval(7,9),
             intervals.Interval(8,10),
             intervals.Interval(9,11),
             intervals.Interval(20,26),
             intervals.Interval(30,38),
             intervals.Interval(38,50),
             intervals.Interval(65,70)]

        intervals.remove_contained_in_list(a)
        self.assertSequenceEqual(a, b)

class Test_length_sum_from_list(unittest.TestCase):
    def test_length_sum_from_list(self):
        '''Test that total length of intervals is summed correctly'''
        a = [intervals.Interval(1,2),
             intervals.Interval(4,5),
             intervals.Interval(10,19)]

        self.assertEqual(14, intervals.length_sum_from_list(a))


if __name__ == '__main__':
    unittest.main()
