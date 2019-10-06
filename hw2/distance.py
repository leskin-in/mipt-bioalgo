#!/usr/bin/env python3


def main():
    pattern = input()
    genomes = input().split()

    distance = hamming_distance(pattern, genomes)

    print(distance)


def _true_hamming_distance(a: str, b: str) -> int:
    """
    Calculate classic hamming distance
    """
    distance = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            distance += 1
    return distance


def hamming_distance(pattern: str, genomes: list) -> int:
    """
    Get hamming distance between 'pattern' and 'genomes'
    """
    distance = 0
    for genome in genomes:
        genome_distance = None
        for i in range(len(genome) - len(pattern) + 1):
            subpattern = genome[i:i + len(pattern)]
            subpattern_distance = _true_hamming_distance(subpattern, pattern)
            if genome_distance is None or subpattern_distance < genome_distance:
                genome_distance = subpattern_distance
        distance += genome_distance
    return distance


if __name__ == '__main__':
    main()
