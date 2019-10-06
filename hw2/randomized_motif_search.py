#!/usr/bin/env python3

from functools import reduce
from random import randint


def main():
    k, total_genomes = list(map(int, input().split()))
    genomes = []
    for _ in range(total_genomes):
        genomes.append(input().upper())

    # best_motifs = randomized_motif_search(k, genomes)
    best_motifs = summarized_randomized_motif_search(k, genomes, 1000)

    for motif in best_motifs:
        print(motif)


ACTG_MAP = {
    'A': 0,
    'C': 1,
    'T': 2,
    'G': 3,
}

ACTG_REVERSE_MAP = {
    0: 'A',
    1: 'C',
    2: 'T',
    3: 'G'
}


def _recount_profile(profile: list, genome: str) -> None:
    """
    Include 'genome' into given 'profile' counts matrix
    """
    for i in range(len(genome)):
        profile[i][ACTG_MAP[genome[i]]] += 1


def _reform_profile(genomes: list) -> list:
    """
    Generate a profile (counts) of the given 'genomes'
    """
    result = [[1 for _ in range(4)] for _ in range(len(genomes[0]))]
    for genome in genomes:
        _recount_profile(result, genome)
    return result


def _profile_to_probs(profile: list, total_genomes: int) -> list:
    """
    Convert a profile (with counts) to a profile of probabilities
    """
    result = []
    for bases in profile:
        result.append([float(bases[i]) / (float(total_genomes) + 4.) for i in range(4)])
    return result


def _profile_consensus_mer(profile: list) -> str:
    """
    Get a consensus mer from the given profile
    """
    result = []
    for bases in profile:
        min_base = max(range(len(bases)), key=bases.__getitem__)
        result.append(ACTG_REVERSE_MAP[min_base])
    return ''.join(result)


def _score_mers(mers: list, consensus: str) -> int:
    """
    Calculate the score of a mer list
    """
    result = 0
    for mer in mers:
        for i in range(len(mer)):
            if mer[i] != consensus[i]:
                result += 1
    return result


def _profile_most_probable_mer(k: int, profile_probs: list, genome: str) -> str:
    """
    Return profile most-probable k-mer
    """
    best_mer = None
    best_probability = None
    for k_start in range(len(genome) - k + 1):
        curr_mer = genome[k_start:k_start + k]
        curr_probability = reduce(
            lambda x, y: x * y,
            [profile_probs[i][ACTG_MAP[curr_mer[i]]] for i in range(len(curr_mer))]
        )
        if best_probability is None or curr_probability > best_probability:
            best_mer = curr_mer
            best_probability = curr_probability
    return best_mer


def randomized_motif_search(k: int, genomes: list) -> list:
    """
    Find profile-best motifs of length 'k' in each of the 'genomes' using a randomized algorithm
    """
    best_motifs = []
    for genome in genomes:
        r = randint(0, len(genomes[0]) - k)
        best_motifs.append(genome[r:r + k])
    tmp_profile = _reform_profile(best_motifs)
    best_motifs_profile_probs = _profile_to_probs(tmp_profile, len(best_motifs))
    best_motifs_score = _score_mers(best_motifs, _profile_consensus_mer(_profile_to_probs(tmp_profile, len(best_motifs))))

    while True:
        curr_motifs = []
        for genome in genomes:
            most_probable_motif = _profile_most_probable_mer(k, best_motifs_profile_probs, genome)
            curr_motifs.append(most_probable_motif)

        curr_motifs_profile = _reform_profile(curr_motifs)
        curr_motifs_score = _score_mers(curr_motifs, _profile_consensus_mer(_profile_to_probs(curr_motifs_profile, len(curr_motifs))))
        if curr_motifs_score < best_motifs_score:
            best_motifs = curr_motifs
            best_motifs_profile_probs = _profile_to_probs(curr_motifs_profile, len(curr_motifs))
            best_motifs_score = curr_motifs_score
        else:
            break

    return best_motifs


def summarized_randomized_motif_search(k: int, genomes: list, n: int) -> list:
    """
    Perform 'randomized_motif_search()' 'N' times
    """
    randomized_search_results = []
    for _ in range(n):
        randomized_search_results.append(randomized_motif_search(k, genomes))

    result = []
    for i in range(len(genomes)):
        randomized_search_counts = {}
        randomized_search_max_count = 0
        randomized_search_max_motif = None

        for search_result in randomized_search_results:
            count = randomized_search_counts.get(search_result[i], 0) + 1
            if count > randomized_search_max_count:
                randomized_search_max_count = count
                randomized_search_max_motif = search_result[i]
            randomized_search_counts[search_result[i]] = count

        result.append(randomized_search_max_motif)
        print(randomized_search_max_count)

    return result


if __name__ == '__main__':
    main()
