import os
import glob
import sys
from setuptools import setup, find_packages


try:
    import numpy
except ImportError:
    print("Error! numpy for Python3 not found.\nPlease install it (e.g. apt-get install python3-numpy)", file=sys.stderr)
    sys.exit(1)

setup(
    name='pyfastaq',
    version='3.2.0',
    description='Script to manipulate FASTA and FASTQ files, plus API for developers',
    packages = find_packages(),
    author='Martin Hunt',
    author_email='mh12@sanger.ac.uk',
    url='https://github.com/sanger-pathogens/Fastaq',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    install_requires=['nose >= 1.3'],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
