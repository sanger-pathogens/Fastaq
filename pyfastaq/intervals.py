class Error (Exception): pass


class Interval:
    '''A class to deal with intervals in a genome. Can do things like intersections, unions etc'''
    def __init__(self, start, end):
        try:
            self.start = int(start)
            self.end = int(end)
        except ValueError:
            raise Error('Error making interval from :"' + str(start) + '" and "' + str(end) + '"')

        if self.end < self.start:
            raise Error('Error making interval ' + str(self) + '.  end < start.')

    def __len__(self):
        return self.end - self.start + 1

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return '(' + str(self.start) + ',' + str(self.end) + ')'

    def __lt__(self, i):
        return self.start < i.start or (self.start == i.start and self.end < i.end)

    def __le__(self, i):
        return self.start < i.start or (self.start == i.start and self.end <= i.end)

    def distance_to_point(self, p):
        '''Returns the distance from the point to the interval. Zero if the point lies inside the interval.'''
        if self.start <= p <= self.end:
            return 0
        else:
            return min(abs(self.start - p), abs(self.end - p))

    def intersects(self, i):
        '''Returns true iff this interval intersects the interval i'''
        return self.start <= i.end and i.start <= self.end

    def contains(self, i):
        '''Returns true iff this interval contains the interval i'''
        return self.start <= i.start and i.end <= self.end

    def union(self, i):
        '''If intervals intersect, returns their union, otherwise returns None'''
        if self.intersects(i) or self.end + 1 == i.start or i.end + 1 == self.start:
            return Interval(min(self.start, i.start), max(self.end, i.end))
        else:
            return None

    def union_fill_gap(self, i):
        '''Like union, but ignores whether the two intervals intersect or not'''
        return Interval(min(self.start, i.start), max(self.end, i.end))

    def intersection(self, i):
        '''If intervals intersect, returns their intersection, otherwise returns None'''
        if self.intersects(i):
            return Interval(max(self.start, i.start), min(self.end, i.end))
        else:
            return None


def intersection(l1, l2):
    '''Returns intersection of two lists.  Assumes the lists are sorted by start positions'''
    if len(l1) == 0 or len(l2) == 0:
        return []

    out = []
    l2_pos = 0

    for l in l1:
        while l2_pos < len(l2) and l2[l2_pos].end < l.start:
            l2_pos += 1

        if l2_pos == len(l2):
            break

        while l2_pos < len(l2) and l.intersects(l2[l2_pos]):
            out.append(l.intersection(l2[l2_pos]))
            l2_pos += 1

        l2_pos = max(0, l2_pos - 1)

    return out


def merge_overlapping_in_list(l):
    '''Sorts list, merges any overlapping intervals, and also adjacent intervals. e.g.
       [0,1], [1,2] would be merge to [0,.2].'''
    i = 0
    l.sort()

    while i < len(l) - 1:
        u = l[i].union(l[i+1])
        if u is not None:
            l[i] = u
            l.pop(i+1)
        else:
            i += 1


def remove_contained_in_list(l):
    '''Sorts list in place, then removes any intervals that are completely
       contained inside another interval'''
    i = 0
    l.sort()

    while i < len(l) - 1:
       if l[i+1].contains(l[i]):
           l.pop(i)
       elif l[i].contains(l[i+1]):
           l.pop(i+1)
       else:
           i += 1


def length_sum_from_list(l):
    '''Returns total length of intervals from a list'''
    return sum([len(x) for x in l])
