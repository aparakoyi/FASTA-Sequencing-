#!/usr/bin/env python3
# nt_fasta_stats.py
"""
Abigail Parakoyi
Assignment 3
February 28, 2022

Calculate the nucleotide statistics for FASTA files
"""

import sys
import argparse
import re


def main():
    """ Business Logic """
    args = get_cli_args()
    infile = args.infile

    # Hardcode the output files where protein and secondary structure are printed out
    outfile = "influenza.stats.txt"
    # Use the get_filehandle() function to open the input file and two output files
    fh_in = get_filehandle(infile, mode="r")
    fh_out = get_filehandle(outfile, mode="w")
    list_header, list_seq = get_fasta_lists(fh_in)
    output_seq_statistics(list_header, list_seq, fh_out)
    fh_in.close()
    fh_out.close()


def get_filehandle(file, mode='r+'):
    """
    Get the filehandle from inputted file with its defined mode of opening (reading or writing)
    :param file: Desired FASTA file
    :param mode: Reading('r) or Writing ('w)
    :return: filehandle of desired FASTA file
    """
    try:
        filehandle = open(file, mode)  # open the file
        return filehandle
    except ValueError as error:  # if a value error occurs accept in cases
        if mode not in {'r', 'w', 'r+', 'w+'}:  # where mode is not any of the defined
            print(f'Wrong {mode} passed in opening')
            raise error
    except OSError as error:  # if an OS error print cannot open file
        print(f'{file} cannot be opened')
        raise error


def get_fasta_lists(fh_in):
    """
    Using the file handle of FASTA file create a list of headers and sequences
    as well as verifying that the lists are equal by calling helper function
    :param fh_in: File handle of FASTA file
    :return: list of headers and list of sequences from FASTA file
    """
    list_headers = []  # initialize an empty list for headers
    list_seqs = []  # initialize an empty list for sequences
    lines = fh_in.readlines()  # read the contents of the input file
    for pos, line in enumerate(lines):
        line = line.rstrip("\n")  # remove all newline character on right side of lines
        if re.match('^>', line):  # if the line has a > then that is a header
            list_headers.append(line)  # append that line to the list of headers
            seq_pos = pos + 1  # initialize the counter for sequences so it never stops
            seqs = ''  # create an empty string to hold sequences
            # As long as the length of the lines is more than the positions of the sequence
            # and is not the header then add that line to the seqs list and add 1 to the seq_pos
            # counter
            while len(lines) > seq_pos and not re.match('^>', lines[seq_pos]):
                clean_lines = lines[seq_pos].replace('\n', '')
                seqs += clean_lines
                seq_pos += 1
            else:
                if seqs != '':  # if the seqs is not an empty list then append it to list of seqs
                    list_seqs.append(seqs)
                else:
                    # and check that the lists are the same length to ensure 1-to-1 correspondence
                    continue
    _verify_lists(list_headers, list_seqs)
    return list_headers, list_seqs


def _verify_lists(list_headers, list_seqs):
    """
    Helper function to verify that the list of header and sequences are equal
    :param list_headers: the list of headers in the FASTA file
    :param list_seqs: list of sequences within the FASTA file
    :return: Boolean True that they are equal or exit system with message when not equal
    """
    if len(list_headers) != len(list_seqs):  # if the length of headers is not = to length of seqs
        # exit the program and print the statement
        sys.exit(f'{list_headers} and {list_seqs} are different in size')
    else:
        return True  # if they are the same just return true


def output_seq_statistics(list_headers, list_seqs, fh_out):
    """
    Calculates the number in nucleotide bases, length, and GC%  per accessions
    :param list_headers: the list of headers from the FASTA file
    :param list_seqs: the list of sequences from the FASTA file
    :param fh_out: the influenza.txt outfile tp print table of results
    :return: the fh_out file with a list of the statistics
    """
    # initialize empty lists to hold statistics table headers values
    accessions_of_headers = []
    a_counts = []
    t_counts = []
    c_counts = []
    g_counts = []
    n_counts = []
    length_count = []
    gc_count = []
    # formatting for the header of table
    fh_out.write("Numbers\tAccession\tA's\tG's\tC's\tT's\tN's\tLength\tGC%\n")
    for header in list_headers:
        header_accession = _get_ncbi_accession(header)  # call this helper function to get header and splice it
        accessions_of_headers.append(header_accession)
    for num, seq in enumerate(list_seqs):
        # call the helper function _get_num_nucleotides to get the number of nucleotides
        # for all bases then append that value to the nucleotide count list
        # follow thi same pattern for the other three statistical measure (length, gc% and number)
        a_nt = _get_num_nucleotides('A', seq)
        a_counts.append(a_nt)
        t_nt = _get_num_nucleotides('T', seq)
        t_counts.append(t_nt)
        c_nt = _get_num_nucleotides('C', seq)
        c_counts.append(c_nt)
        g_nt = _get_num_nucleotides('G', seq)
        g_counts.append(g_nt)
        n_nt = _get_num_nucleotides('N', seq)
        n_counts.append(n_nt)
        length = len(seq)
        length_count.append(length)
        gc_percentage = round((((c_nt + g_nt)/length) * 100), 1)
        gc_count.append(gc_percentage)
        fh_out.write(f'{num + 1}\t{accessions_of_headers[num]}\t{a_counts[num]}\t'
                     f'{g_counts[num]}\t{c_counts[num]}\t{t_counts[num]}\t{n_counts[num]}\t'
                     f'{length_count[num]}\t{gc_count[num]}\n')


def _get_num_nucleotides(base, seq):
    """
    Counts the number of our specified nucleotides in our reference. If not there exit the program
    :param base: the nucleotide base pairs
    :param seq: the list_seq
    :return: a count of the specified nucleotide
    """
    if base not in {'A', 'G', 'C', 'T', 'N'}:  # if the base is not any of the indicated exit the program with err mes
        sys.exit(f'Did not code this condition')
    else:
        # if in specified nucleotides allocate it to a new variable
        temp_nucleotides = seq.count(base)
    return temp_nucleotides


def _get_ncbi_accession(header_string):
    """
    Just get the accession of the header not the whole header
    :param header_string: this is the list of headers
    :return: the list of headers spliced to just the accession variable
    """
    # header_string is re-assigned to new variable header_accessions
    header_accessions = header_string
    header_spliced = header_accessions[1:9]  # splices the list_headers to only get the accession value
    return header_spliced


def get_cli_args():
    """
    Access the code through the command line
    :return: Instance of argparse argument
    """
    parser = argparse.ArgumentParser(description='Provide a FASTA file to '
                                                 'generate nucleotide statistics')
    parser.add_argument('-i', '--infile',
                        dest='infile',
                        type=str,
                        help='Path to file to open',
                        required=True)
    parser.add_argument('-o', '--outfile',
                        dest='outfile',
                        type=str,
                        help='Path to file to write',
                        required=True)
    return parser.parse_args()


if __name__ == '__main__':
    main()
