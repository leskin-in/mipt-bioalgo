#!/usr/bin/env python3


def main():
    mapping = {}
    letters = list(input().split())
    for i in range(len(letters)):
        mapping[letters[i]] = i
    print(mapping)


if __name__ == '__main__':
    main()
