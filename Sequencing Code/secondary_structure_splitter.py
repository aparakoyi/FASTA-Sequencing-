#!/usr/bin/env python3
# secondary_structure_splitter.py
"""
Abigail Parakoyi
Assignment 3
February 28, 2022

Get the number of proteins and the number of secondary structures from a FASTA file
"""
# Import system and argparse
import re
import sys
import argparse


def main():
    """ Business Logic """

    args = get_cli_args()
    infile = args.infile

    # Hardcode the output files where protein and secondary structure are printed out
    outfile1 = "pbd_protein.fasta"
    outfile2 = "pbd_ss.fasta"
    # Use the get_filehandle() function to open the input file and two output files
    fh_in = get_filehandle(infile, mode="r")
    fh_out1 = get_filehandle(file=outfile1, mode="w")
    fh_out2 = get_filehandle(file=outfile2, mode="w")
    # Use the get_fasta_lists() function to get the list of headers and list of seq
    list_headers, list_seqs = get_fasta_lists(fh_in)
    # Use the output_results_to_file() function to get the # of proteins ans ss in list
    num_proteins, num_ss = output_results_to_file(headers=list_headers, seqs=list_seqs,
                                                  fh_out1=fh_out1, fh_out2=fh_out2)
    # writing the results to CLI through STDERR
    sys.stderr.write(f"Found {num_proteins} protein sequences\n")
    sys.stderr.write(f"Found {num_ss} ss sequences\n")
    # closing of files
    fh_in.close()
    fh_out1.close()
    fh_out2.close()


# Get file handle to sequence in the fasta file
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


def output_results_to_file(headers, seqs, fh_out1, fh_out2):
    """
    Write the information from the list of headers and seqs to the 2 output files
    :param headers: list of FASTA file headers
    :param seqs: list of FASTA file sequences
    :param fh_out1: outfile 1 'pbd_protein.fasta'
    :param fh_out2: outfile 2 'pbd_ss.fasta'
    :return: the number of proteins and number of ss
    """
    num_protein = 0  # initialize a counter for the # of proteins
    num_ss = 0  # initialize a counter for the # of ss
    for base, headers in enumerate(headers):
        if re.match('^.*(sequence)$', headers):
            fh_out1.write(headers + '\n')  # write the header to outfile 1
            fh_out1.write(seqs[base] + '\n')  # write the sequence to outfile 1
            num_protein += 1  # increment the protein counter by 1
        elif re.match('^.*(secstr)$', headers):
            fh_out2.write(headers + '\n')  # write the header to outfile 2
            fh_out2.write(seqs[base] + '\n')  # write the sequence to outfile 2
            num_ss += 1  # increment the ss counter by 1
    return num_protein, num_ss


def get_cli_args():
    """
    Just to get the command line options for this program through argpase
    :return: Instance of argparse arguments
    """
    parser = argparse.ArgumentParser(description='Provide a FASTA file to perform splitting '
                                                 'on sequence and secondary structure')
    parser.add_argument('-i', '--infile',
                        dest='infile',
                        type=str,
                        help='Path to file to open',
                        required=True)
    return parser.parse_args()


if __name__ == '__main__':
    main()
