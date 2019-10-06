#!/usr/bin/env python3

from functools import reduce


def main():
    k, total_genomes = list(map(int, input().split()))
    genomes = []
    for _ in range(total_genomes):
        genomes.append(input().upper())

    best_motifs = greedy_motif_search(k, genomes)

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


def _profile_to_probs(profile: list, total_genomes: int) -> list:
    """
    Convert a profile (with counts) to a profile of probabilities
    """
    result = []
    for bases in profile:
        result.append([float(bases[i]) / float(total_genomes) for i in range(4)])
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


def greedy_motif_search(k: int, genomes: list) -> list:
    """
    Find profile-best motifs of length 'k' in each of the 'genomes'
    """
    best_motifs = [genome[:k] for genome in genomes]
    profile_temp = [[0 for _ in range(4)] for _ in range(k)]
    for best_motif_temp in best_motifs:
        _recount_profile(profile_temp, best_motif_temp)
    best_motifs_score = _score_mers(best_motifs, _profile_consensus_mer(profile_temp))

    for k_start in range(len(genomes[0]) - k + 1):
        current_motifs = [genomes[0][k_start:k_start + k]]
        profile = [[0 for _ in range(4)] for _ in range(k)]

        for total_genomes_in_profile in range(1, len(genomes)):
            _recount_profile(profile, current_motifs[-1])
            profile_probs = _profile_to_probs(profile, total_genomes_in_profile)
            current_most_probable_motif = _profile_most_probable_mer(k, profile_probs, genomes[total_genomes_in_profile])
            current_motifs.append(current_most_probable_motif)

        _recount_profile(profile, current_motifs[-1])
        consensus = _profile_consensus_mer(profile)  # We can use usual profile; for purposes of this function, it is the same as the probabilities-profile
        current_motifs_score = _score_mers(current_motifs, consensus)
        if best_motifs_score is None or current_motifs_score < best_motifs_score:
            best_motifs = current_motifs
            best_motifs_score = current_motifs_score

    return best_motifs


if __name__ == '__main__':
    main()
