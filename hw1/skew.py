#!/usr/bin/env python3


SKEW_DECREASER = 'C'
SKEW_INCREASER = 'G'


def main():
    genome = input().upper()

    min_skews = min_skew_positions(genome)
    for c in min_skews:
        print(c, end=' ')
    print()


def min_skew_positions(genome: str) -> list:
    """
    Find all positions where substrings (starting from 0) of genome are such that skew (defined by SKEW_INCREASER and SKEW_INCREASER occurences) is minimal over the whole string
    :param genome: string to examine
    :return: a list with all indexes where skew is minimized
    """

    min_skew = 0
    result = [0]

    curr_skew = 0
    for i in range(len(genome)):
        n = genome[i]

        if n == SKEW_INCREASER:
            curr_skew += 1
        elif n == SKEW_DECREASER:
            curr_skew -= 1

        if curr_skew < min_skew:
            min_skew = curr_skew
            result.clear()

        if curr_skew == min_skew:
            result.append(i + 1)

    return result


if __name__ == '__main__':
    main()
