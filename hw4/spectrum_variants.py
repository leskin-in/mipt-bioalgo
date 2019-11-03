#!/usr/bin/env python3


def main():
    m = int(input())
    result = spectrum_variants(m)
    print(result)


AMINO_MASSES = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]


def _attach_amino_mass(mass: int) -> list:
    """
    "Attach" all possible amino acid masses to the given one and return the masses of results
    """
    return list(map(lambda m: m + mass, AMINO_MASSES))


def spectrum_variants(mass: int) -> int:
    """
    Calculate the number of all peptides with the given mass
    :param mass: mass of the peptide to analyse
    :return: the number of peptides, each of which has the given 'mass'
    """
    result = 0

    masses = {0: 1}
    while len(masses) > 0:
        next_masses = {}
        for current_mass in masses.keys():
            new_masses = _attach_amino_mass(current_mass)
            for m in new_masses:
                if m == mass:
                    result += masses.get(current_mass)
                elif m < mass:
                    next_masses[m] = next_masses.get(m, 0) + masses.get(current_mass)
        masses = next_masses

    return result


if __name__ == '__main__':
    main()
