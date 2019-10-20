#!/usr/bin/env python3

import sys


def main():
    k = None
    d = None
    pairs = []
    for l in sys.stdin:
        if k is None:
            k, d = list(map(int, l[:-1].split()))
            continue
        pairs.append(l[:-1])

    result = reconstruct_gapped_dna(pairs, k, d)
    print(result)


def reconstruct_gapped_dna(patterns: list, k: int, d: int) -> str:
    """
    Reconstruct the DNA from its paired De Bruijn graph path
    """

    patterns = [pattern.split('|') for pattern in patterns]

    # Convert the path found to mers
    euler_mers = [patterns[0][0]]
    euler_last_mers = [None for _ in range(d + 1)]
    euler_last_mers[-1] = patterns[0][1]
    for pattern in patterns[1:]:
        euler_mers.append(pattern[0][-1])
        euler_last_mers.pop(0)
        euler_last_mers.append(pattern[1])

    for i in range(len(euler_last_mers) - 1):
        euler_mers.append(euler_last_mers[i][0])
    euler_mers.append(euler_last_mers[-1])

    # Form the result
    result = ''.join(euler_mers)

    return result


if __name__ == '__main__':
    main()
