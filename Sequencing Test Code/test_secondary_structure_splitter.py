"""Test suite for secondary_structure_splitter.py"""
import pytest
import os


from secondary_structure_splitter import (get_filehandle, get_fasta_lists, _verify_lists)

# ignore all "Missing function or method docstring" since this is a unit test
# pylint: disable=C0116
# ignore all "Function name "test_get_filehandle_4_OSError
# " doesn't conform to snake_case naming style"
# pylint: disable=C0103

FILE_TO_TEST = "test_file.txt"
FILE_TO_TEST_PARSING = "test_file.fasta"
FASTA_STR_TO_TEST = """\
>101M:A:sequence
MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGA
ILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKEL
GYQG
>101M:A:secstr
    HHHHHHHHHHHHHHGGGHHHHHHHHHHHHHHH GGGGGG TTTTT  SHHHHHH HHHHHHHHHHHHHHHH
HHTTTT  HHHHHHHHHHHHHTS   HHHHHHHHHHHHHHHHHH GGG SHHHHHHHHHHHHHHHHHHHHHHHHT
T   
>102L:A:sequence
MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAA
VRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRV
ITTFRTGTWDAYKNL
>102L:A:secstr
  HHHHHHHHH  EEEEEE TTS EEEETTEEEESSS TTTHHHHHHHHHHTS  TTB  HHHHHHHHHHHHHHH
HHHHHH TTHHHHHHHS HHHHHHHHHHHHHHHHHHHHT HHHHHHHHTT HHHHHHHHHSSHHHHHSHHHHHHH
HHHHHHSSSGGG 
"""


def test_existing_get_filehandle_for_reading():
    # does it open a file for reading
    # create a test file
    _create_file_for_testing(FILE_TO_TEST)
    # test
    test = get_filehandle(FILE_TO_TEST, "r")
    assert hasattr(test, "readline") is True, "Not able to open for reading"
    test.close()
    os.remove(FILE_TO_TEST)


def _create_file_for_testing(file):
    open(file, "r").close()


def test_existing_get_filehandle_for_writing():
    _create_file_for_testing_writing(FILE_TO_TEST)
    test = get_filehandle(FILE_TO_TEST)
    assert hasattr(test, "write") is True, "Not able to open for writing"
    test.close()
    os.remove(FILE_TO_TEST)


def _create_file_for_testing_writing(file):
    open(file, "w").close()


def test_get_filehandle_4_OSError():
    # does it raise OSError
    # this should exit
    with pytest.raises(OSError):
        get_filehandle("does_not_exist.txt", "r")


def test_get_filehandle_4_ValueError():
    # does it raise OSError
    # this should exit
    with pytest.raises(ValueError):
        get_filehandle(FILE_TO_TEST, "rrrr")


def test_get_fasta_lists():
    # is it able to parse through a file
    _create_fasta_file_for_testing()
    test = get_fasta_lists(FILE_TO_TEST_PARSING)
    print(test)


def _create_fasta_file_for_testing():
    with open(FILE_TO_TEST_PARSING, "w") as fh:
        fh.write(FASTA_STR_TO_TEST)
        fh.close()


def test__verify_lists():
    _create_fasta_file_for_testing()
    test = _verify_lists(list_headers=[275628, 364821764, 3789173], list_seqs=['Adeyagedyiq'])
    print(test)
