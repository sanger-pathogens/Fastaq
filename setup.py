import os
import glob
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='Fastaq',
    version='1.5.0',
    description='Scripts to manipulate FASTA and FASTQ files, plus API for developers',
    long_description=read('README.md'),
    packages = find_packages(),
    author='Martin Hunt',
    author_email='mh12@sanger.ac.uk',
    url='https://github.com/martinghunt/Fastaq',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    install_requires=['nose >= 1.3'],
    license='GPLv3',
)
