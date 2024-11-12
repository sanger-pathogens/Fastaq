import argparse
import sys

tasks = {
    'acgtn_only':             'Replace every non acgtnACGTN with an N',
    'add_indels':             'Deletes or inserts bases at given position(s)',
    'caf_to_fastq':           'Converts a CAF file to FASTQ format',
    'capillary_to_pairs':     'Converts file of capillary reads to paired and unpaired files',
    'chunker':                'Splits sequences into equal sized chunks',
    'count_sequences':        'Counts the sequences in input file',
    'deinterleave':           'Splits interleaved paired file into two separate files',
    'enumerate_names':        'Renames sequences in a file, calling them 1,2,3... etc',
    'expand_nucleotides':     'Makes every combination of degenerate nucleotides',
    'fasta_to_fastq':         'Convert FASTA and .qual to FASTQ',
    'filter':                 'Filter sequences to get a subset of them',
    'get_ids':                'Get the ID of each sequence',
    'get_seq_flanking_gaps':  'Gets the sequences flanking gaps',
    'interleave':             'Interleaves two files, output is alternating between fwd/rev reads',
    'make_random_contigs':    'Make contigs of random sequence',
    'merge':                  'Converts multi sequence file to a single sequence',
    'replace_bases':          'Replaces all occurrences of one letter with another',
    'reverse_complement':     'Reverse complement all sequences',
    'scaffolds_to_contigs':   'Creates a file of contigs from a file of scaffolds',
    'search_for_seq':         'Find all exact matches to a string (and its reverse complement)',
    'sequence_trim':          'Trim exact matches to a given string off the start of every sequence',
    'sort_by_name':           'Sorts sequences in lexographical (name) order',  
    'sort_by_size':           'Sorts sequences in length order', 
    'split_by_base_count':    'Split multi sequence file into separate files',
    'strip_illumina_suffix':  'Strips /1 or /2 off the end of every read name',
    'to_boulderio':           'Converts to Boulder-IO format, used by primer3',
    'to_fasta':               'Converts a variety of input formats to nicely formatted FASTA format',
    'to_fake_qual':           'Make fake quality scores file',
    'to_mira_xml':            'Create an xml file from a file of reads, for use with Mira assembler',
    'to_orfs_gff':            'Writes a GFF file of open reading frames',
    'to_perfect_reads':       'Make perfect paired reads from reference',
    'to_random_subset':       'Make a random sample of sequences (and optionally mates as well)',
    'to_tiling_bam':          'Make a BAM file of reads uniformly spread across the input reference',
    'to_unique_by_id':        'Remove duplicate sequences, based on their names. Keep longest seqs',
    'translate':              'Translate all sequences in input nucleotide sequences',
    'trim_contigs':           'Trims a set number of bases off the end of every contig',
    'trim_ends':              'Trim fixed number of bases of start and/or end of every sequence',
    'trim_Ns_at_end':         'Trims all Ns at the start/end of all sequences',
    'version':                'Print version number and exit',
}


def print_usage_and_exit():
    print('Usage: fastaq <command> [options]', file=sys.stderr)
    print('\nTo get minimal usage for a command use:\nfastaq command', file=sys.stderr)
    print('\nTo get full help for a command use one of:\nfastaq command -h\nfastaq command --help\n', file=sys.stderr)
    print('\nAvailable commands:\n', file=sys.stderr)
    max_task_length = max([len(x) for x in list(tasks.keys())])
    for task in sorted(tasks):
        print('{{0: <{}}}'.format(max_task_length).format(task), tasks[task], sep='  ', file=sys.stderr)
    sys.exit(1)


def main():
    if len(sys.argv) == 1 or sys.argv[1] in ['-h', '-help', '--help']:
        print_usage_and_exit()

    task = sys.argv.pop(1)

    if task not in tasks:
        print('Task "' + task + '" not recognised. Cannot continue.\n', file=sys.stderr)
        print_usage_and_exit()


    exec('import pyfastaq.runners.' + task)
    exec('pyfastaq.runners.' + task + '.run("' + tasks[task] + '")')


if __name__ == "__main__":
    main()