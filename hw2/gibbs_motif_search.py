#!/usr/bin/env python3

from functools import reduce
from random import randint, choices


def main():
    k, total_genomes, n = list(map(int, input().split()))
    genomes = []
    for _ in range(total_genomes):
        genomes.append(input().upper())

    best_motifs = summarized_gibbs_motif_search(k, genomes, n, 20)

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


def _profile(genomes: list) -> list:
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


def _mer_probability(profile_probs: list, genome: str) -> float:
    """
    Return probability of the given 'genome' in the given 'profile_probs'
    """
    return reduce(
        lambda x, y: x * y,
        [profile_probs[i][ACTG_MAP[genome[i]]] for i in range(len(genome))]
    )


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


def _random_probable_motif(profile_probs: list, k: int, genome: str):
    """
    Return a random most-probable by the given profile motif from the given genome
    """
    mers = []
    probabilities = []

    for k_start in range(len(genome) - k + 1):
        mers.append(genome[k_start:k_start + k])
        probabilities.append(reduce(
            lambda x, y: x * y,
            [profile_probs[i][ACTG_MAP[mers[-1][i]]] for i in range(len(mers[-1]))]
        ))

    return choices(mers, probabilities)[0]


def gibbs_motif_search(k: int, genomes: list, n: int) -> list:
    """
    Find profile-best motifs of length 'k' in each of the 'genomes' using a randomized algorithm
    """
    best_motifs = []
    for genome in genomes:
        r = randint(0, len(genomes[0]) - k)
        best_motifs.append(genome[r:r + k])
    best_motifs_score = _score_mers(best_motifs, _profile_consensus_mer(_profile_to_probs(_profile(best_motifs), len(best_motifs))))

    for j in range(n):
        i_rotate = randint(0, len(genomes) - 1)

        rotated_best_motifs = best_motifs[0:i_rotate]
        rotated_best_motifs.extend(best_motifs[i_rotate + 1:len(best_motifs)])

        profile_probs = _profile_to_probs(_profile(rotated_best_motifs), len(rotated_best_motifs))
        probable_motif = _random_probable_motif(profile_probs, k, genomes[i_rotate])

        rotated_best_motifs.insert(i_rotate, probable_motif)
        rotated_score = _score_mers(rotated_best_motifs, _profile_consensus_mer(_profile_to_probs(_profile(rotated_best_motifs), len(rotated_best_motifs))))
        if rotated_score < best_motifs_score:
            best_motifs = rotated_best_motifs
            best_motifs_score = rotated_score

    return best_motifs


def summarized_gibbs_motif_search(k: int, genomes: list, n: int, samples: int) -> list:
    """
    Perform 'gibbs_motif_search()' 'samples' times
    """
    randomized_search_results = []
    for _ in range(samples):
        randomized_search_results.append(gibbs_motif_search(k, genomes, n))

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

    return result


if __name__ == '__main__':
    main()
