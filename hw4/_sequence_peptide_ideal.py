#!/usr/bin/env python3

import sys


def main():
    spectrum = list(map(int, input().split()))
    result = sequence_peptide(spectrum)
    print(' '.join(['-'.join(list(map(str, peptide))) for peptide in result]))


AMINO_MASSES = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]


def _attach_amino_mass(peptide: list) -> list:
    """
    "Attach" all possible amino acid masses to the given peptide and return the resulting peptides
    """
    result = [peptide[:] for _ in range(len(AMINO_MASSES))]
    for i in range(len(result)):
        result[i].append(AMINO_MASSES[i])
    return result


def _is_peptide_consistent_with_spectrum(peptide: list, spectrum: set) -> bool:
    """
    Check if the given peptide consistent with the given spectrum
    """
    for amino in peptide:
        if amino not in spectrum:
            return False
    return True


def _is_peptide_cyclospectrum_equal_to_spectrum(peptide: list, spectrum: set) -> bool:
    """
    Check if the given peptide consistent with the given spectrum
    """
    cyclospectrum = set()
    for cycle_len in range(1, len(peptide)):
        peptide_extended = peptide + peptide[:(cycle_len - 1)]
        for cycle_start_pos in range(len(peptide)):
            cyclospectrum.add(_peptide_mass(peptide_extended[cycle_start_pos:(cycle_start_pos + cycle_len)]))

    cyclospectrum.add(0)
    cyclospectrum.add(_peptide_mass(peptide))

    return cyclospectrum == spectrum


def _peptide_mass(peptide: list):
    """
    Calculate peptide mass
    """
    return sum(peptide)


def sequence_peptide(spectrum: list) -> list:
    """
    Sequence a peptide from its ideal spectrum
    """
    parent_mass = spectrum[-1]
    spectrum = set(spectrum)

    peptides = [[]]
    result = []

    while len(peptides) > 0:
        print(len(peptides[0]), len(peptides), file=sys.stderr)
        next_peptides = []
        for peptide in peptides:
            new_peptides = _attach_amino_mass(peptide)
            for new_peptide in new_peptides:
                new_peptide_mass = _peptide_mass(new_peptide)
                if new_peptide_mass < parent_mass:
                    if _is_peptide_consistent_with_spectrum(new_peptide, spectrum):
                        next_peptides.append(new_peptide)
                elif new_peptide_mass == parent_mass:
                    if _is_peptide_cyclospectrum_equal_to_spectrum(new_peptide, spectrum):
                        result.append(new_peptide)

        peptides = next_peptides

    return result


if __name__ == '__main__':
    main()
