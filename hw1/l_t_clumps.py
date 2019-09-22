#!/usr/bin/env python3

from collections import deque


def main():
    genome = input().upper()
    k, l, t = list(map(int, input().split()))

    clumps = l_t_clumps(genome, k, l, t)
    for c in clumps:
        print(c, end=' ')
    print()


def l_t_clumps(genome: str, k: int, l: int, t: int):
    """
    Find patterns of length 'k' occuring at leats 't' times in each (continuous) substring of length 'l'
    :param genome: string to examing
    :param k: pattern length
    :param l: window (substring) length
    :param t: number of occurrences
    :return: a list containing all the named patterns
    """
    clump_set = set()

    pqueue_size_max = l - k + 1
    pqueue = deque(maxlen=pqueue_size_max)

    pdict = {}

    for sstart in range(len(genome) - k + 1):
        seq_current = genome[sstart:sstart + k]

        if sstart >= pqueue_size_max:
            seq_previous = pqueue.popleft()
            if pdict[seq_previous] == 1:
                pdict.pop(seq_previous)
            else:
                pdict[seq_previous] -= 1

        pqueue.append(seq_current)
        pdict[seq_current] = pdict.get(seq_current, 0) + 1

        if pdict[seq_current] >= t:
            clump_set.add(seq_current)

    return list(clump_set)


if __name__ == '__main__':
    main()
