#!/usr/bin/env python3
"""Exact finite checks for Supplementary Appendix G (Python 3.11+, stdlib).

Reads the actual interval/projection tables. Uses exact rational arithmetic
and exhaustive finite graph checks, without numerical solvers. This checks
finite data; the spectral theorems and the asymptotic argument remain
analytical proofs, not computer-certified theorems.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json
import re

BASE = Path(__file__).resolve().parent.parent
SOURCE = (BASE / "threshold_supplement.tex").read_text()
SQRT5_UPPER = F("2.24")


def table_rows(label):
    tables = re.findall(r"\\begin\{table\}.*?\\end\{table\}", SOURCE, re.S)
    table, = [t for t in tables if rf"\label{{{label}}}" in t]
    body = table.split(r"\midrule", 1)[1].split(r"\bottomrule", 1)[0]
    return [[c.strip().strip("$") for c in line.rstrip("\\").split("&")]
            for line in body.splitlines() if "&" in line]


def check_intervals():
    cover = table_rows("tab:threshold-cover")
    assert len(cover) == 8
    assert F("2.235") ** 2 < 5 < SQRT5_UPPER ** 2
    assert F("2.5") ** 3 < 16
    assert F("4.325") ** 3 < 81
    assert F("6.08") ** 2 < 37 < F(73, 12) ** 2
    assert F(1845, 421) < F(22, 5)
    b_values = [F(1), F(1), F(2), F(2), F(2), F(2), F("2.235"), F("2.5")]
    # Strict rational lower bounds, except for the exact capacity of K7.
    h_lower = [None, F("6.235"), F("6.235"), F("6.325"),
               F("6.47"), F("6.705"), F(7), None]
    q_H = [F(73, 12), 4 + SQRT5_UPPER, 4 + SQRT5_UPPER, F(32, 5),
           2 + 2 * SQRT5_UPPER, 3 * SQRT5_UPPER, F(7), F(85, 12)]
    h_H = [F(37, 6), F(13, 2), F(13, 2), F(13, 2), F(7), F(15, 2), F(7), F(43, 6)]
    h_J = [F(1), F(1), F(2), F(2), F(2), F(2), F(5, 2), None]
    q_J = [F(1), F(1), F(2), F(2), F(2), F(2), SQRT5_UPPER, F(8, 3)]
    partners = ["R", r"K_4\join C_5", r"K_4\join C_5", r"K_2\join Z",
                r"K_2\join2C_5", "3C_5", "K_7", r"K_1\join R"]
    auxiliaries = ["K_1", "K_1", "K_2", "K_2", "K_2", "K_2", "C_5", "Q"]
    report = []
    for index, (row, b, lower) in enumerate(zip(cover, b_values, h_lower)):
        assert row[:2] == [auxiliaries[index], partners[index]]
        theta_H, haemers, theta_B, upper, endpoint = map(F, row[2:])
        assert theta_H >= q_H[index]
        assert theta_B >= 4 + q_J[index]
        h_B_bound = F(16, 3) + h_J[index] if h_J[index] is not None else h_H[index]
        assert haemers >= min(h_B_bound, h_H[index])
        denominator = haemers - 6 + theta_H - theta_B
        assert haemers > 6 and theta_H > theta_B
        weight = (theta_H - theta_B) / denominator
        assert 0 < weight < 1
        assert upper == (haemers * theta_H - 6 * theta_B) / denominator
        assert 6 * weight + (1 - weight) * theta_H == upper
        assert haemers * weight + (1 - weight) * theta_B == upper
        assert upper < endpoint
        assert endpoint > b and endpoint * (endpoint - b) < 32
        if index in (0, 7):
            assert (endpoint - (index == 7)) ** 2 < 37
        elif index == 6:
            assert endpoint < lower
        else:
            assert endpoint <= lower
        report.append({"row": index + 1, "derived_weight": str(weight),
                       "U": str(upper), "L": str(endpoint)})
    overlaps = [F(cover[i][6]) - F(cover[i + 1][5]) for i in range(7)]
    shift_overlap = F(cover[-1][6]) - F(cover[0][5]) - 1
    assert min(overlaps) == F(1, 1480) > 0
    assert shift_overlap == F(2, 225) > 0
    assert F("6.02") - F(cover[0][5]) == F(2, 225) > 0
    return {"rows": report, "consecutive_overlaps": list(map(str, overlaps)),
            "integer_shift_overlap": str(shift_overlap)}


def clique(words, adjacent):
    assert len(words) == len(set(words))
    assert all(any(adjacent(x, y) for x, y in zip(u, v))
               for u, v in combinations(words, 2))


def check_binary():
    vertices = [x for x in range(64) if x.bit_count() % 2]
    basis = [1 << i for i in range(6)]
    assert all((x & y).bit_count() % 2 == 0 for x, y in combinations(basis, 2))
    for v in range(64):
        if v.bit_count() % 2:
            continue
        transform = lambda x: x ^ (v if (x & v).bit_count() % 2 else 0)
        assert {transform(x) for x in vertices} == set(vertices)
        assert all(transform(transform(x)) == x for x in vertices)
        assert all((transform(x) & transform(y)).bit_count() % 2 == (x & y).bit_count() % 2
                   for x, y in combinations(vertices, 2))
    # Permutations are transitive within each weight class; these two
    # transvections send e1 into each of the other odd weight classes.
    assert (1 ^ 15).bit_count() == 3
    assert (1 ^ 63).bit_count() == 5
    # U6's displayed representation is rational: q=4, q^(-1/2)=1/2.
    vectors = {
        x: [F(1, 2)] + [
            F(int(z == x) - int(z == (x ^ 63)), 2)
            + F((-1) ** ((x & z).bit_count() % 2), 8)
            for z in vertices
        ] for x in vertices
    }
    for x in vertices:
        assert sum(a * a for a in vectors[x]) == 1
        for y in vertices:
            if x != y and (x & y).bit_count() % 2:
                assert sum(a * b for a, b in zip(vectors[x], vectors[y])) == 0
    return {"transvections_checked": 32, "odd_weight_orbits_joined": [1, 3, 5],
            "independent_set_in_complement": 6, "squared_handle": "1/4"}


def check_circulant_and_cycles():
    match = re.search(r"S=\\\{([0-9,]+)\\\}", SOURCE)
    distances = set(map(int, match.group(1).split(",")))
    adjacency = lambda x, y: min((x - y) % 37, (y - x) % 37) in distances
    for x, y in combinations(range(37), 2):
        assert adjacency(x, y) != adjacency(6 * x % 37, 6 * y % 37)
    six = [0, 6, 12, 18, 24, 30]
    assert all(adjacency(x, y) for x, y in combinations(six, 2))
    independent = [6 * x % 37 for x in six]
    assert not any(adjacency(x, y) for x, y in combinations(independent, 2))
    assert Counter((x + t) % 37 for x in independent for t in range(37)) == Counter({x: 6 for x in range(37)})
    clique([(i, 6 * i % 37) for i in range(37)], adjacency)
    clique([(i, 2 * i % 5) for i in range(5)], lambda x, y: (x - y) % 5 in (1, 4))
    clique([(x, y, (2 * x + 4 * y) % 9) for x in range(9) for y in range(9)],
           lambda x, y: (x - y) % 9 not in (0, 1, 8))
    return {"R37_code_pairs": 666, "C5_code_pairs": 10, "complement_C9_code_pairs": 3240}


def check_clebsch():
    # Polynomial basis modulo z^4+z^3+z^2+z+1, with beta=z.
    def multiply(a, b):
        value = 0
        while b:
            if b & 1:
                value ^= a
            a <<= 1
            if a & 16:
                a ^= 31
            b >>= 1
        return value
    subgroup = {1, 2, 4, 8, 15}
    assert multiply(8, 2) == 15 and multiply(15, 2) == 1
    a = 3
    clique([(x, multiply(a, x), multiply(multiply(a, a), x)) for x in range(16)],
           lambda x, y: (x ^ y) in subgroup)
    gram = [[F(int(x == y)) + F(int(x != y and (x ^ y) not in subgroup), 2)
             for y in range(16)] for x in range(16)]
    for u in range(16):
        character = [(-1) ** ((u & x).bit_count() % 2) for x in range(16)]
        eigenvalue = F(6) if u == 0 else F(1 - (4 - 2 * u.bit_count() + (-1) ** u.bit_count()), 2)
        assert eigenvalue in (0, 2, 6)
        assert all(sum(gram[x][y] * character[y] for y in range(16)) == eigenvalue * character[x]
                   for x in range(16))
    assert all(sum(row) == 6 for row in gram)
    assert F(6 ** 2, 16 * 6) == F(3, 8)
    return {"code_pairs": 120, "Gram_eigenvalues": [0, 2, 6], "squared_handle": "3/8"}


if __name__ == "__main__":
    results = {"intervals": check_intervals(), "binary_U6": check_binary(),
               "circulant_and_cycles": check_circulant_and_cycles(),
               "Clebsch": check_clebsch(), "status": "All exact finite checks passed"}
    print(json.dumps(results, indent=2))
