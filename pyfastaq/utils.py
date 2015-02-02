import os
import sys
import subprocess
import shlex

class Error (Exception): pass

def open_file_read(filename):
    if filename == '-':
        f = sys.stdin
    elif filename.endswith('.gz'):
        # first check that the file is OK according to gunzip
        retcode = subprocess.call('gunzip -t ' + filename, shell=True)
        if retcode != 0:
            raise Error("Error opening for reading gzipped file '" + filename + "'")

        # now open the file
        f = os.popen('gunzip -c ' + filename)
    else:
        try:
            f = open(filename)
        except:
            raise Error("Error opening for reading file '" + filename + "'")

    return f


def open_file_write(filename):
    if filename == '-':
        f = sys.stdout
    elif filename.endswith('.gz'):
        if not os.path.exists(os.path.abspath(os.path.dirname(filename))):
            raise Error("Error opening for writing gzipped file '" + filename + "'")

        try:
            f = os.popen('gzip -9 -c > ' + filename, 'w')
        except:
            raise Error("Error opening for writing gzipped file '" + filename + "'")
    else:
        try:
            f = open(filename, 'w')
        except:
            raise Error("Error opening for writing file '" + filename + "'")

    return f


def close(filehandle):
    if filehandle not in [sys.stdout, sys.stderr]:
        filehandle.close()


def file_transpose(f_in, f_out, sep_in=None, sep_out='\t'):
    rows = []
    f = open_file_read(f_in)
    for line in f:
        rows.append(line.rstrip().split(sep_in))
    close(f)

    columns_out = max([len(x) for x in rows])

    for r in rows:
        r += ['.'] * (columns_out - len(r))

    f = open_file_write(f_out)
    for i in range(columns_out):
        print(sep_out.join([str(rows[x][i]) for x in range(len(rows))]), file=f)

    close(f)


def syscall(cmd):
    retcode = subprocess.call(cmd, shell=True)

    if retcode != 0:
        raise Error("Error in system call. Command was:\n" + cmd)


def syscall_get_stdout(cmd):
    try:
        out = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE).communicate()[0].decode('utf-8').rstrip()
        return out.split('\n')
    except:
        raise Error('Error in system call. I tried to run:\n' + str(cmd))


