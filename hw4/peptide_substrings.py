#!/usr/bin/env python3


def main():
    dna = input()
    peptide = input()
    result = peptide_encoding_sequences(dna, peptide)
    for r in result:
        print(r)


AMINO_TABLE = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAT', 'AAC'],
    'D': ['GAT', 'GAC'],
    'C': ['TGT', 'TGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'K': ['AAA', 'AAG'],
    'M': ['ATG'],
    'F': ['TTT', 'TTC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
}


def _amino_to_dna(amino: str) -> list:
    """
    Convert an amino acid to a DNA sequence that encodes this amino acid
    :param amino:
    :return: A list of all DNA sequences encoding the given amino acid
    """
    return AMINO_TABLE[amino]


RC_PIPELINE_A = ['A', 'a', 'T']
RC_PIPELINE_C = ['C', 'c', 'G']
RC_PIPELINE_T = ['T', 'A']
RC_PIPELINE_G = ['G', 'C']


def _reverse_complement(genome: str) -> str:
    """
    Get reverse complement of the given genome
    :param genome: sequence to transform
    :return: RC
    """
    return genome.\
        replace(RC_PIPELINE_A[0], RC_PIPELINE_A[1]).replace(RC_PIPELINE_C[0], RC_PIPELINE_C[1]).\
        replace(RC_PIPELINE_T[0], RC_PIPELINE_T[1]).replace(RC_PIPELINE_G[0], RC_PIPELINE_G[1]).\
        replace(RC_PIPELINE_A[1], RC_PIPELINE_A[2]).replace(RC_PIPELINE_C[1], RC_PIPELINE_C[2])\
        [::-1]


def _peptide_to_dna(protein: str) -> list:
    """
    Convert a peptide into a list of DNA sequences that encode it
    """
    result = ['']
    for amino in protein:
        next_dna_variants = _amino_to_dna(amino)
        next_result = []
        for prev_dna in result:
            for dna_variant in next_dna_variants:
                next_result.append(prev_dna + dna_variant)
        result = next_result
    return result


def peptide_encoding_sequences(dna: str, peptide: str) -> list:
    """
    Find all (sub)sequences encoding the given 'peptide' in the given 'dna'
    :return: a list of sequences found
    """
    all_peptide_encoding_sequences = set(_peptide_to_dna(peptide))

    result = []
    for ss_start_i in range(len(dna) - len(peptide) * 3 + 1):
        ss = dna[ss_start_i:(ss_start_i + len(peptide) * 3)]
        if ss in all_peptide_encoding_sequences or _reverse_complement(ss) in all_peptide_encoding_sequences:
            result.append(ss)

    return result


if __name__ == '__main__':
    main()
