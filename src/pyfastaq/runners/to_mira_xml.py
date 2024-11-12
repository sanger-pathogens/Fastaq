import argparse
from pyfastaq import tasks

def run(description):
    parser = argparse.ArgumentParser(
        description = description,
        usage = 'fastaq to_mira_xml <infile> <xml_out>')
    parser.add_argument('infile', help='Name of input fasta/q file')
    parser.add_argument('xml_out', help='Name of output xml file')
    options = parser.parse_args()
    tasks.fastaq_to_mira_xml(options.infile, options.xml_out)
