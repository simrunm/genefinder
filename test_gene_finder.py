from collections import Counter
import pytest

from gene_finder import (
    get_complement,
    get_reverse_complement,
    rest_of_orf,
    find_all_orfs_one_frame,
    find_all_orfs,
    find_all_orfs_both_strands,
    find_longest_orf,
    encode_amino_acids,
)


# Define sets of test cases.
get_complement_cases = [
    # Check that the complement of A is T.
    ("A", "T"),
    # Check that the complement of C is G.
    ("C", "G"),
]

get_reverse_complement_cases = [
    # Check a single nucleotide, which should be the same as the complement.
    ("A", "T"),
    # Check that the strand is reversed.
    ("ATG", "CAT"),
]

rest_of_orf_cases = [
    # Check a start followed by a stop.
    ("ATGTGA", "ATG"),
    # Check a case without a stop codon.
    ("ATGAAA", "ATGAAA"),
    # Check a case without a stop codon where the length is not a multiple of 3.
    ("ATGA", "ATGA"),
    # Check to make sure it doesn't count stop codons that are
    # not in multiples of 3.
    ("ATGATGATTTGA", "ATGATGATT"),
    # Check a case with a stop codon where the length is not a multiple 3.
    ("ATGTGAA", "ATG"),
]

find_all_orfs_one_frame_cases = [
    # Check a strand with a single ORF.
    ("ATGTGA", ["ATG"]),
    # Check a strand with two ORFs.
    ("ATGTAAATGAAATAA", ["ATG", "ATGAAA"]),
    # Check a strand with a nested ORF that doesn't start with a start codon.
    ("AATATGATGAAATGA", ["ATGATGAAA"]),
    # Check a strand with no ORFs.
    ("AAAAAAAA", []),
    # Check multiple start followed by stop codons with one nested ORF.
    ("ATGATGTGAATGTGA", ["ATGATG", "ATG"]),
    # Check a strand that ends with a start codon.
    ("ATGAATTGAATG", ["ATGAAT", "ATG"]),
    # Check a strand with one ORF in frame and one ORF out of frame.
    ("ATGTGAAATGTGA", ["ATG"]),
]

find_all_orfs_cases = [
    # This case from find_all_orfs has no ORFs in other frames, so it should
    # return the same result as in the one_frame case.
    ("ATGTAAATGAAATAA", ["ATG", "ATGAAA"]),
    # Check a strand with one ORF in each frame.
    ("AATGGAATGAATTGATGTGA", ["ATGAAT", "ATGGAA", "ATG"]),
]

find_all_orfs_both_strands_cases = [
    # Test a short strand starting with a start codon whose reverse complement
    # is itself. Thus this should return two copies of the same ORF.
    ("ATGCAT", ["ATGCAT", "ATGCAT"]),
    # Test a strand with codons in multiple frames for itself
    # and its reverse complement.
    ("TCAATGAATGTGACTTGACAT", ["ATGAATGTGACT", "ATG", "ATGTCAAGTCACATTCAT"]),
]

get_longest_orf_cases = [
    # An ORF covering the whole strand is by default the longest ORF.
    ("ATGAAAAAAAAA", "ATGAAAAAAAAA"),
    # A strand with multiple ORFs.
    ("ATGATGTGAATGTGA", "ATGATG"),
    # A strand with the longest ORF in the reverse complement.
    ("ATGTGATTAAAACAT", "ATGTTT")
]

coding_strand_to_aa_cases = [
    # Check a single start codon.
    ("ATG", "M"),
    # Check a case in which the length is not a multiple of 3.
    ("ATGCCCGCTTT", "MPA"),
    # Check a case with only two nucleotides.
    ("AA", ""),
]


# Define additional testing lists and functions that check other properties of
# functions in gene_finder.py.
@pytest.mark.parametrize("nucleotide", ["A", "T", "C", "G"])
def test_double_complement(nucleotide):
    """
    Check that taking the complement of a complement of a nucleotide produces
    the original nucleotide.

    Args:
        nucleotide: A single-character string representing one of the four DNA
            nucleotides.
    """
    assert get_complement(get_complement(nucleotide)) == nucleotide


################################################################################
# Don't change anything below these lines.
################################################################################


# Define standard testing functions to check functions' outputs given certain
# inputs defined above.
@pytest.mark.parametrize("nucleotide,complement", get_complement_cases)
def test_get_complement(nucleotide, complement):
    assert get_complement(nucleotide) == complement


@pytest.mark.parametrize("strand,reverse_complement",
                         get_reverse_complement_cases)
def test_get_reverse_complement(strand, reverse_complement):
    assert get_reverse_complement(strand) == reverse_complement


@pytest.mark.parametrize("strand,rest", rest_of_orf_cases)
def test_rest_of_orf(strand, rest):
    assert rest_of_orf(strand) == rest


@pytest.mark.parametrize("strand,orfs", find_all_orfs_one_frame_cases)
def test_find_all_orfs_oneframe(strand, orfs):
    assert Counter(find_all_orfs_one_frame(strand)) == Counter(orfs)


@pytest.mark.parametrize("strand,orfs", find_all_orfs_cases)
def test_find_all_orfs(strand, orfs):
    assert Counter(find_all_orfs(strand)) == Counter(orfs)


@pytest.mark.parametrize("strand,orfs", find_all_orfs_both_strands_cases)
def test_find_all_orfs_both_strands(strand, orfs):
    assert Counter(find_all_orfs_both_strands(strand)) == Counter(orfs)


@pytest.mark.parametrize("strand,orf", get_longest_orf_cases)
def test_get_longest_orf(strand, orf):
    assert find_longest_orf(strand) == orf


@pytest.mark.parametrize("strand,protein", coding_strand_to_aa_cases)
def test_coding_strand_to_aa(strand, protein):
    assert encode_amino_acids(strand) == protein
