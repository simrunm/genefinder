"""
Library for finding potential genes in a strand of DNA.
"""
import helpers


def get_complement(nucleotide):
    """Return the complementary DNA nucleotide.

    Each DNA nucleotide A,T,C or G has a complementary nucleotide where
    A complements T and C complements G.

    Args:
        nucleotide: A string with a single character that is either A,T,C or G.

    Returns:
        A string with a single character character containing
        the complementary nucleotide.
    """

    if nucleotide == "A":
        return "T"
    if nucleotide == "T":
        return "A"
    if nucleotide == "C":
        return "G"
    if nucleotide == "G":
        return "C"


def get_reverse_complement(strand):
    """Return the reverse complementary DNA strand.

    The reverse complement of a DNA strand is when the strand is
    flipped and then the complement of each nucleotide is returned.

    Args:
        strand: A string repesenting a single strand of DNA.

    """
    # Initialize an empty string
    complementary_strand = ""
    # Loop through the strand in reverse order and add the complement
    # to the empty strand
    for nucleotide in strand[::-1]:
        complementary_strand += (get_complement(nucleotide))
    return complementary_strand


def rest_of_orf(strand):
    """Return the sequence of nucleotides representing the rest of the orf
       not including the stop codon.

    An orf is a series of necleotides within a DNA strand that are
    translated into proteins. Orfs begin and end with specific start and end
    codons. Codons are translated in sequences of 3.

    Args:
        strand: A string repesenting a single strand of DNA.

    Returns:
        A string with the sequence of nucleotides that represent the codon not
        including the stop codon. If there is no stop codon in strand, then the
        strand is returned.
    """
    from helpers import amino_acid
    # Create a list with all the codons which are sequences of
    # three nucleotides.
    for i in range(0, len(strand), 3):
        if amino_acid(strand[i:i + 3]) == "*" and len(strand[i:i + 3]) == 3:
            return strand[:i]
    # return the whole strand if no stop codon is found.
    return strand


def find_all_orfs_one_frame(strand):
    """Return a list of strings representing all in-frame ORFs
    found in that strand.

    Each returned ORF is a multiple of three nucleotides from the
    start of the strand making it in-frame. Nested orfs are
    not returned.

    Args:
        strand: A string repesenting a single strand of DNA.

    Returns:
        A list of strings representing all in-frame ORFs found in the strand.
    """
    # importing the helper function amino acid
    from helpers import amino_acid
    # creating an empty list where orfs will added
    orfs = []
    i = 0
    # looping throught the strand for orfs
    while i < len(strand):
        # checking for a start codon in the strand
        if amino_acid(strand[i:i + 3]) == "M":
            # adding the orf from the point a start codon is found to the list
            # and skippping to the point after that orf to avoid adding
            # nested orfs
            orfs.append(rest_of_orf(strand[i:]))
            i = i + len(rest_of_orf(strand[i:]))
        # if no start codon is found, running the loop from the next codon
        # which is 3 nucleotides away.
        i += 3
    return orfs


def find_all_orfs(strand):
    """Return a list of strings representing all ORFs found in that strand.

     Orfs from different frames where the frame is moved by one or two
     nucleotides are also included. Nested orfs within a single frame are
     included.

    Args:
        strand: A string repesenting a single strand of DNA.

    Returns:
        A list of strings with all the orfs found in that strand.
    """
    all_orfs = []
    # find the orfs with the index shifted by one and two nucleotides.
    for i in range(3):
        for elements in find_all_orfs_one_frame(strand[i:]):
            # don't add the element if it isn't an empty list.
            if elements != []:
                all_orfs.append(elements)
    return all_orfs


def find_all_orfs_both_strands(strand):
    """Return a list of strings representing all ORFs found in the strand
       or its reverse complement.

    Args:
        strand: A string repesenting a single strand of DNA.
    """
    return find_all_orfs(strand) + find_all_orfs(get_reverse_complement(strand))


def find_longest_orf(strand):
    """Return the longest ORF found in the DNA strand or its reverse complement.

    Args:
        strand: A string repesenting a single strand of DNA.

    Returns:
        A string that represents the longest ORF.
    """

    return max(find_all_orfs_both_strands(strand), key=len)


def noncoding_orf_threshold(strand, num_trials):
    """Return the length of the shortest ORF out of the longest ORFs.

    A single DNA strand is randomly shuffled and the longest ORF from
    each shuffled version is chosen. Out of all the longest ORFs found,
     the length of the shortest one is chosen.

    Args:
        strand: A string repesenting a single strand of DNA.
        num_trials: An integer representing the number of times
        the DNA strand is shuffled.

    Returns:
        An integer that represents length of the shortest long ORF
    """
    import random
    from helpers import shuffle
    # initializing shortest as the whole strand beacuse the any orf
    # will be shorter than that
    shortest_long_orf = strand
    # finds the longest orf for each shuffled version
    for _ in range(num_trials):
        shuffled_strand = shuffle(strand)
        long_orf = find_longest_orf(shuffled_strand)
    # compares the length of the longest orf with the current shortest orf
        if len(long_orf) < len(shortest_long_orf):
            shortest_long_orf = long_orf
    return len(shortest_long_orf)


def encode_amino_acids(orf):
    """Return a string representing the sequence of amino acids
    corresponding to the inputted ORF.

    Args:
        orf: A string repesenting an orf
    """
    from helpers import amino_acid
    # Create a list with all the codons which are sequences of
    # three nucleotides.
    codons = []
    for i in range(0, len(orf), 3):
        codons.append(orf[i:i + 3])
    # Create a string of amino acids sequences from the codons if the codons
    # are 3 characters
    amino_acid_sequence = ""
    for codon in codons:
        if len(codon) == 3:
            amino_acid_sequence += amino_acid(codon)
    return amino_acid_sequence


def find_genes(path):
    """Returns a list of amino acid sequences that represent ORFs that are
    longer than the cutoff length.

    Args:
        path: a string path representing the location of a file in FASTA
        format

    Returns:
        A list containing all the amino acid sequences for ORFs
        longer than cutoff length
    """
    from helpers import load_fasta_file
    # creating the strand of DNA
    dna_strand = load_fasta_file(path)
    # finding the cutoff length
    cutoff_length = noncoding_orf_threshold(dna_strand, 1500)
    # finding all the orfs
    all_orfs = find_all_orfs_both_strands(dna_strand)
    all_amino_acids = []
    # if the orfs are longer than the cutoff length then find the amino acid
    # sequence that corresponds to that.
    for orf in all_orfs:
        if len(orf) > cutoff_length:
            all_amino_acids.append(encode_amino_acids(orf))
    return all_amino_acids
