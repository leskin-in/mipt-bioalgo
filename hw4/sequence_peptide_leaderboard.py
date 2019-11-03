#!/usr/bin/env python3


def main():
    lb_size = int(input())
    spectrum = list(map(int, input().split()))
    result = sequence_peptide(spectrum, lb_size)
    print('-'.join(list(map(str, result))))


AMINO_MASSES = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]


def _attach_amino_mass(peptide: tuple) -> list:
    """
    "Attach" all possible amino acid masses to the given peptide and return the resulting peptides
    """
    return [peptide + (mass,) for mass in AMINO_MASSES]


def _cyclospectrum_score(peptide: tuple, spectrum: set) -> int:
    """
    Calculate the cyclospectrum score for the given peptide
    """
    cyclospectrum = set()
    for cycle_len in range(1, len(peptide)):
        peptide_extended = peptide + peptide[:(cycle_len - 1)]
        for cycle_start_pos in range(len(peptide)):
            cyclospectrum.add(_peptide_mass(peptide_extended[cycle_start_pos:(cycle_start_pos + cycle_len)]))

    cyclospectrum.add(0)
    cyclospectrum.add(_peptide_mass(peptide))

    return len(spectrum.intersection(cyclospectrum))


def _peptide_mass(peptide: tuple):
    """
    Calculate peptide mass
    """
    return sum(peptide)


def _expand_leaderboard(leaderboard: list, spectrum: set):
    """
    Expand the leaderboard
    :param leaderboard: leaderboard of peptides
    :return: a new leaderboard (noncut), ordered DESC
    """
    lb_len = len(leaderboard)
    for i in range(lb_len):
        p_pair = leaderboard[i]
        new_peptides = _attach_amino_mass(p_pair[1])
        for new_peptide in new_peptides:
            leaderboard.append([_cyclospectrum_score(new_peptide, spectrum), new_peptide])

    for i in range(lb_len):
        leaderboard.pop(0)

    leaderboard.sort(key=lambda p: p[0], reverse=True)


def sequence_peptide(spectrum: list, lb_size: int) -> list:
    """
    Sequence a peptide using a leaderboad algorithm
    :param spectrum: spectrum of a peptide (not an ideal one)
    :param lb_size: leaderboard size
    :return: the leading peptide
    """
    parent_mass = spectrum[-1]
    spectrum = set(spectrum)

    leaderboard = [[0, ()]]
    leader = [0, ()]

    while len(leaderboard) > 0:
        _expand_leaderboard(leaderboard, spectrum)

        i = 0
        while True:
            if i >= len(leaderboard) or i >= lb_size:
                break

            p_pair = leaderboard[i]
            peptide_mass = _peptide_mass(p_pair[1])

            if peptide_mass == parent_mass:
                if p_pair[0] > leader[0]:
                    leader = p_pair
                else:
                    leaderboard.pop(i)
                    continue
            elif peptide_mass > parent_mass:
                leaderboard.pop(i)
                continue

            i += 1

        leaderboard = leaderboard[:lb_size]

    return leader[1]


if __name__ == '__main__':
    main()
