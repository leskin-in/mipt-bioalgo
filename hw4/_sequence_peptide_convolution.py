#!/usr/bin/env python3


def main():
    alphabet_size = int(input())
    lb_size = int(input())
    spectrum = list(map(int, input().split()))

    amino_alphabet = calculate_amino_alphabet(spectrum, alphabet_size)
    result = sequence_peptide(spectrum, lb_size, amino_alphabet)
    print('-'.join(list(map(str, result))))


def _attach_amino_mass(peptide: tuple, amino_alphabet: list) -> list:
    """
    "Attach" all possible amino acid masses to the given peptide and return the resulting peptides
    """
    return [peptide + (mass,) for mass in amino_alphabet]


def _cyclospectrum_score(peptide: tuple, spectrum: dict) -> int:
    """
    Calculate the cyclospectrum score for the given peptide
    """
    cyclospectrum = {}
    for cycle_len in range(1, len(peptide)):
        peptide_extended = peptide + peptide[:(cycle_len - 1)]
        for cycle_start_pos in range(len(peptide)):
            current_mass = _peptide_mass(peptide_extended[cycle_start_pos:(cycle_start_pos + cycle_len)])
            cyclospectrum[current_mass] = cyclospectrum.get(current_mass, 0) + 1

    current_mass = _peptide_mass(peptide)
    cyclospectrum[current_mass] = cyclospectrum.get(current_mass, 0) + 1

    result = 0
    for mass in cyclospectrum.keys():
        result += spectrum.get(mass, 0)  # spectrum.get(mass, 0) * cyclospectrum[mass]
    return result


def _peptide_mass(peptide: tuple):
    """
    Calculate peptide mass
    """
    return sum(peptide)


def _expand_leaderboard(leaderboard: list, spectrum: dict, amino_alphabet: list):
    """
    Expand the leaderboard
    :param leaderboard: leaderboard of peptides
    :param spectrum: spectrum of the peptide to find
    :param amino_alphabet: an alphabet of amino acids to use to expand the leaderboard
    :return: a new leaderboard (noncut), ordered DESC
    """
    lb_len = len(leaderboard)
    for i in range(lb_len):
        p_pair = leaderboard[i]
        new_peptides = _attach_amino_mass(p_pair[1], amino_alphabet)
        for new_peptide in new_peptides:
            leaderboard.append([_cyclospectrum_score(new_peptide, spectrum), new_peptide])

    for i in range(lb_len):
        leaderboard.pop(0)

    leaderboard.sort(key=lambda p: p[0], reverse=True)


def _build_spectrum(spectrum: list) -> dict:
    """
    Convert a list-spectrum into a map-spectrum
    """
    result = {}
    for amino in spectrum:
        result[amino] = 1
    return result


def sequence_peptide(spectrum: list, lb_size: int, amino_alphabet: list) -> list:
    """
    Sequence a peptide using a leaderboad algorithm
    :param spectrum: spectrum of a peptide (not an ideal one)
    :param lb_size: leaderboard size
    :param amino_alphabet: alphabet of amino acids
    :return: the leading peptide
    """
    parent_mass = spectrum[-1]
    spectrum = _build_spectrum(spectrum)

    leaderboard = [[0, ()]]
    leader = [0, ()]

    while len(leaderboard) > 0:
        _expand_leaderboard(leaderboard, spectrum, amino_alphabet)

        i = 0
        while i < len(leaderboard):
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
        print(leaderboard)

        leaderboard_result = []
        j = 0
        while j < len(leaderboard):
            leaderboard_result.append(leaderboard[j])
            j += 1
            if j == len(leaderboard):
                break
            if j >= lb_size and leaderboard[j - 1][0] > leaderboard[j][0]:
                break

        leaderboard = leaderboard_result
        print(leaderboard)

    return leader[1]


def calculate_amino_alphabet(spectrum: list, alphabet_size: int) -> list:
    """
    Calculate an alphabet of amino acids using convolutions method over the given 'spectrum'
    """
    alphabet_set = {}
    spectrum = spectrum[:]
    spectrum.insert(0, 0)

    for i in range(len(spectrum)):
        c1 = spectrum[i]
        for j in range(i, len(spectrum)):
            if i == j:
                continue
            c2 = spectrum[j]
            conv = c2 - c1
            if 57 <= conv <= 200:
                alphabet_set[conv] = alphabet_set.get(conv, 0) + 1

    alphabet_list = []
    for k in alphabet_set.keys():
        v = alphabet_set.get(k)
        alphabet_list.append((v, k))
    alphabet_list.sort(key=lambda p: p[0], reverse=True)

    result = []
    i = 0
    while i < len(alphabet_list):
        result.append(alphabet_list[i][1])
        i += 1
        if i == len(alphabet_list):
            break
        if i >= alphabet_size and alphabet_list[i - 1][0] > alphabet_list[i][0]:
            break

    return result


if __name__ == '__main__':
    main()
