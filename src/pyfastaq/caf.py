from pyfastaq import sequences, utils

class Error (Exception): pass

def file_reader(fname):
    f = utils.open_file_read(fname)
    c = Caf()

    while c.get_next_from_file(f):
        yield c

    utils.close(f)


class Caf:
    def __init__(self):
        self.id = None
        self.seq = None
        self.insert_min = None
        self.insert_max = None
        self.ligation = None
        self.clone = None
        self.clip_start = None
        self.clip_end = None


    def __eq__(self, other):
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False


    def get_next_from_file(self, f):
        self.__init__()
        line = f.readline()
        if not line:
            return None
        while line == '\n':
            line = f.readline()

        if not line.startswith('DNA : '):
            raise  Error("Error reading caf file. Expected line starting with 'DNA : ...'")

        self.id = line.rstrip().split()[2]

        line = f.readline()
        seq = []

        while line != '\n':
            seq.append(line.rstrip())
            line = f.readline()

        self.seq = sequences.Fasta(self.id, ''.join(seq))

        line = f.readline()
        if not line.startswith('BaseQuality : '):
            raise  Error("Error reading caf file. Expected line starting with 'BaseQuality : ...'")

        quals = [int(x) for x in f.readline().rstrip().split()]
        self.seq = self.seq.to_Fastq(quals)

        line = f.readline()
        assert line == '\n'
        line = f.readline()

        while line not in ['', '\n']:
            a = line.rstrip().split()
            if a[0] == 'Insert_size':
                self.insert_min, self.insert_max = int(a[1]), int(a[2])
            elif a[0] == 'Ligation_no':
                self.ligation = a[1]
            elif a[0] == 'Clone':
                self.clone = a[1]
            elif a[0] == 'Clipping' and a[1] == 'QUAL':
                self.clip_start, self.clip_end = int(a[2]) - 1, int(a[3]) - 1

            line = f.readline()

        return True
