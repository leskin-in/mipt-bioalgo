#!/usr/bin/env python3


def main():
    genome = input().upper()
    k, mismatches = list(map(int, input().split()))

    frequents = frequent_words_with_mismatches_and_rc(genome, k, mismatches)
    for word in frequents:
        print(word, end=' ')
    print()


LIST_A = ['C', 'T', 'G']
LIST_C = ['A', 'T', 'G']
LIST_T = ['C', 'A', 'G']
LIST_G = ['C', 'T', 'A']


def _generate_immediate_neighbours(pattern: str) -> list:
    """
    Generate immediate (different by one mismatch) neighbours of the given genome pattern
    :param pattern: a pattern to examine
    :return: neighbourhood, NOT including the given pattern
    """
    generated = []
    for i in range(len(pattern)):
        if pattern[i] == 'A':
            generated.extend([pattern[:i] + c + pattern[i + 1:] for c in LIST_A])
        elif pattern[i] == 'C':
            generated.extend([pattern[:i] + c + pattern[i + 1:] for c in LIST_C])
        elif pattern[i] == 'T':
            generated.extend([pattern[:i] + c + pattern[i + 1:] for c in LIST_T])
        elif pattern[i] == 'G':
            generated.extend([pattern[:i] + c + pattern[i + 1:] for c in LIST_G])

    return generated


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


def generate_neighbours(pattern: str, mismatches: int) -> set:
    """
    Generate neighbours for the given pattern (genome string)
    :param pattern: genome pattern
    :param mismatches: number of mismatches to generate neighbours
    :return: a set of patterns in the neighbourhood, including the 'pattern' itself
    """
    neighbourhood = set()
    neighbourhood.add(pattern)

    curr_patterns = [pattern]
    next_patterns = []

    for curr_mismatches in range(mismatches):
        for curr_pattern in curr_patterns:
            for neighbour in _generate_immediate_neighbours(curr_pattern):
                if neighbour not in neighbourhood:
                    neighbourhood.add(neighbour)
                    next_patterns.append(neighbour)

        curr_patterns = next_patterns
        next_patterns = []

    return neighbourhood


def frequent_words_with_mismatches_and_rc(genome: str, k: int, mismatches: int) -> list:
    """
    Find most frequent patterns of length 'k' in 'genome' with at most 'mismatches' mismatches, counting reverse complements as original string or its mismatches
    :param genome: genome to examine
    :param k: pattern length
    :param mismatches: maximum number of mismatches
    :return: a list of patterns
    """
    frequencies = {}

    for sstart in range(len(genome) - k + 1):
        sequence = genome[sstart:sstart + k]
        neighbours = list(generate_neighbours(sequence, mismatches))
        for neighbour in neighbours:
            neighbour_rc = _reverse_complement(neighbour)
            frequencies[neighbour] = frequencies.get(neighbour, 0) + 1
            frequencies[neighbour_rc] = frequencies.get(neighbour_rc, 0) + 1

    max_frequency = max(frequencies.values())

    result = []
    for pattern, frequency in frequencies.items():
        if frequency == max_frequency:
            result.append(pattern)

    return result


if __name__ == '__main__':
    main()
