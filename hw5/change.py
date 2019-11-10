#!/usr/bin/env python3


def main():
    sum_to_form = int(input())
    denominations = list(map(int, input().split(',')))
    result = calculate_change_number_dyn(sum_to_form, denominations)
    print(result)


def calculate_change_number_dyn(sum_to_form: int, denominations: list) -> int:
    """
    Calculate the minimum number of 'denominations' required to make 'sum'
    :param sum_to_form: a sum to form
    :param denominations: possible denominations ("steps")
    :return: a total amount of denominations to use to reach 'sum'
    """
    results = [-1 for _ in range(sum_to_form + 1)]
    results[0] = 0

    for s in range(sum_to_form):
        next_sum_indexes = [s + d for d in denominations]
        next_sum = results[s] + 1
        for next_sum_index in next_sum_indexes:
            if next_sum_index > sum_to_form:
                continue
            if results[next_sum_index] < 0 or results[next_sum_index] > next_sum:
                results[next_sum_index] = next_sum

    return results[-1]


if __name__ == '__main__':
    main()
