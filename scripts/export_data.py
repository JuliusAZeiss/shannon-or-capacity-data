#!/usr/bin/env python3
"""Export exact finite data from this package; Python 3.11+, standard library.

The copied manuscript and verification modules are authoritative. These JSON
files are deterministic convenience exports, not a substitute for the proofs.
"""

import argparse
from fractions import Fraction
import hashlib
import importlib.util
from itertools import combinations
import json
import os
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parent.parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "verification" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def matrix_digest(matrix):
    return hashlib.sha256(json.dumps(matrix, separators=(",", ":")).encode("ascii")).hexdigest()


def graph(labels, adjacency):
    require(len(labels) == len(adjacency), "Wrong label count")
    require(len(set(labels)) == len(labels), "Repeated vertex labels")
    n = len(labels)
    require(all(len(row) == n for row in adjacency), "Nonsquare adjacency")
    require(all(adjacency[i][i] == 0 for i in range(n)), "Graph loop")
    require(all(adjacency[i][j] in (0, 1) and adjacency[i][j] == adjacency[j][i]
                for i in range(n) for j in range(n)), "Invalid adjacency")
    return {"vertex_labels": labels, "order": n, "size": sum(map(sum, adjacency)) // 2,
            "adjacency_rows": ["".join(map(str, row)) for row in adjacency],
            "adjacency_matrix_sha256": matrix_digest(adjacency)}


def code(adjacency, words):
    words = [list(word) for word in words]
    require(words and len({tuple(word) for word in words}) == len(words), "Repeated codeword")
    length = len(words[0])
    require(all(len(word) == length and all(0 <= v < len(adjacency) for v in word)
                for word in words), "Invalid codeword")
    require(all(any(adjacency[x][y] for x, y in zip(u, v))
                for u, v in combinations(words, 2)), "Nonadjacent codewords")
    return {"alphabet": "zero-based row indices in the associated graph",
            "block_length": length, "word_count": len(words), "words": words,
            "pairs_checked": len(words) * (len(words) - 1) // 2}


def complement(adjacency):
    return [[int(i != j and not adjacency[i][j]) for j in range(len(adjacency))]
            for i in range(len(adjacency))]


def disjoint_union(left, right):
    return [row + [0] * len(right) for row in left] + [[0] * len(left) + row for row in right]


def binary_graph(dimension):
    labels = [x for x in range(1 << dimension) if x.bit_count() % 2]
    adjacency = [[int(x != y and (x & y).bit_count() % 2 == 0) for y in labels] for x in labels]
    return labels, adjacency


def subsix():
    check = load("check_subsix")
    bd, hc = check.bd, check.hc
    checked = check.main()
    b, d = bd.adjacency()
    h, labels = hc.read_h_from_tex(ROOT / "subsix_supplement.tex")
    require(h == [[int(v) for v in row] for row in check.H_ROWS], "Printed H differs from module")
    g = disjoint_union(b, d)
    pos = {v: i for i, v in enumerate(labels)}
    h_words = [(i, 12 * i % 29) for i in range(29) if i not in (1, 17, 28)] + [(28, 12)]
    scalar = [[pow(check.dot3(u, u), -1, 3) * check.dot3(u, v) % 3
               for v in bd.VECTORS] for u in bd.VECTORS]
    assigned = [[a for a, colour in enumerate(check.COLOR_CLASSES) if i in colour]
                for i in range(len(labels))]
    basis = [[int(k == a) for k in range(11)] for pair in assigned for a in pair]
    block_fit = [[sum(a * b for a, b in zip(u, v)) % 3 for v in basis] for u in basis]
    generators = tuple(bd.permutation(cycles) for cycles in bd.CYCLES)
    certificates = {}
    for name, adj, scale, diagonal, entries in (
            ("B", b, 24, 104, bd.B_DATA), ("D", d, 100000000, 419004112, bd.D_DATA)):
        matrix = bd.reconstruct(adj, scale, diagonal, entries, generators)
        require(matrix_digest(matrix) == checked["matrix_SHA256"][name], f"{name} matrix digest")
        certificates[name] = {"matrix": matrix, "scale": scale, "diagonal": diagonal,
                              "capacity_upper_bound": str(Fraction(diagonal + scale, scale)),
                              "matrix_sha256": matrix_digest(matrix), "orbit_entries": entries,
                              "automorphism_generators": generators}
    hm = hc.compact_h_certificate(h, labels)["matrix"]
    require(matrix_digest(hm) == checked["matrix_SHA256"]["H"], "H matrix digest")
    certificates["H"] = {"matrix": hm, "scale": hc.SCALE, "diagonal": hc.DIAGONAL,
                          "capacity_upper_bound": str(hc.Q), "matrix_sha256": matrix_digest(hm),
                          "orbit_entries_native_labels": hc.H_ORBIT_ENTRIES}
    return {"schema_version": 1, "paper_anchor": "sec:subsixproof (Appendix F)",
            "graphs": {"B": graph(list(range(27)), b), "D": graph(list(range(27)), d),
                       "G": graph(list(range(54)), g), "H": graph(labels, h)},
            "G_vertex_convention": "indices 0..26 are B; 27..53 are D in the same vector order",
            "source_vectors_F3": bd.VECTORS,
            "OR_codes": {"G": code(g, [(i, i + 27) for i in range(27)]),
                         "H": code(h, [(pos[x], pos[y]) for x, y in h_words])},
            "F3_fitting_matrices": {
                "B": {"local_block_size": 1, "rank": check.rank3(scalar), "matrix": scalar},
                "H": {"local_block_size": 2, "rank": check.rank3(block_fit), "matrix": block_fit,
                      "matrix_row_order": "vertex row index, then increasing zero-based colour index",
                      "twofold_colours_per_vertex": assigned,
                      "independent_colour_classes_row_indices": check.COLOR_CLASSES}},
            "integer_Gram_certificates": certificates,
            "interpolation": {"weight": str(hc.WEIGHT), "strict_product_upper": str(hc.L),
                              "factor_capacity_lower": "sqrt(27)",
                              "threshold_interval": "[5193957/1000000, sqrt(27))"}}


def thresholds():
    check = load("check_exact")
    checks = {"intervals": check.check_intervals(), "binary_U6": check.check_binary(),
              "circulant_and_cycles": check.check_circulant_and_cycles(), "Clebsch": check.check_clebsch()}
    labels, a = binary_graph(6)
    d = complement(a)
    source = (ROOT / "threshold_supplement.tex").read_text(encoding="utf-8")
    match = re.search(r"S=\\\{([0-9,]+)\\\}", source)
    require(match is not None, "Missing circulant distances")
    distances = list(map(int, match[1].split(",")))
    r = [[int(min((x - y) % 37, (y - x) % 37) in distances) for y in range(37)] for x in range(37)]
    c5 = [[int((x - y) % 5 in (1, 4)) for y in range(5)] for x in range(5)]
    z = [[int((x - y) % 9 not in (0, 1, 8)) for y in range(9)] for x in range(9)]
    subgroup = {1, 2, 4, 8, 15}
    q = [[int((x ^ y) in subgroup) for y in range(16)] for x in range(16)]

    def multiply(x, y):
        result = 0
        while y:
            if y & 1:
                result ^= x
            x <<= 1
            if x & 16:
                x ^= 31
            y >>= 1
        return result

    table = []
    for row in check.table_rows("tab:threshold-cover"):
        table.append(dict(zip(("J", "H", "s", "a", "d", "U", "L"),
                              row[:2] + [str(Fraction(value)) for value in row[2:]])))
    graphs = {"U6": graph(labels, a), "complement_U6": graph(labels, d),
              "R": graph(list(range(37)), r), "C5": graph(list(range(5)), c5),
              "Z": graph(list(range(9)), z), "Q": graph(list(range(16)), q)}
    return {"schema_version": 1, "paper_anchor": "sec:thresholdproof (Appendix G)",
            "graphs": graphs, "R_undirected_cyclic_distances": distances,
            "U6_label_convention": "binary digits of the integer label are coordinates in F2^6",
            "Q_label_convention": "polynomial basis in F16 modulo z^4+z^3+z^2+z+1",
            "OR_codes": {"R": code(r, [(i, 6 * i % 37) for i in range(37)]),
                         "C5": code(c5, [(i, 2 * i % 5) for i in range(5)]),
                         "Z": code(z, [(x, y, (2 * x + 4 * y) % 9) for x in range(9) for y in range(9)]),
                         "Q": code(q, [(x, multiply(3, x), multiply(multiply(3, 3), x)) for x in range(16)])},
            "coverage_table": table, "finite_check_summary": checks,
            "closed_cutoff": {"value": "60101/10000", "interpolation_weight": "93/106",
                              "optimized_endpoint": "6*exp(2*d^2/(ell+3*d))",
                              "d": "log(sqrt(37)/6)", "ell": "log(6/5)"}}


def simple_example():
    labels, u = binary_graph(8)
    complement_u = complement(u)
    f9 = disjoint_union(u, complement_u)
    codewords = [(i, 128 + i) for i in range(128)]
    require(len(labels) == 128, "Wrong U8 order")
    return {"schema_version": 1, "paper_anchor": "sec:simple (Section 2)",
            "origin": "Generated directly from the analytical binary-graph definition in main.tex",
            "graphs": {"U8": graph(labels, u), "complement_U8": graph(labels, complement_u),
                       "F9": graph([f"U8:{x}" for x in labels] + [f"complement_U8:{x}" for x in labels], f9),
                       "K9": graph(list(range(9)), [[int(i != j) for j in range(9)] for i in range(9)])},
            "binary_label_convention": "binary digits of the integer label are coordinates in F2^8",
            "OR_codes": {"F9": code(f9, codewords)},
            "scope": "Checks the finite graph definition and 128-word code; capacity upper bounds use the paper's analytical proofs."}


def main():
    require(sys.version_info >= (3, 11), "Python 3.11 or later is required")
    require(__debug__ and not os.environ.get("PYTHONOPTIMIZE"),
            "Assertions must be enabled: remove -O/-OO and unset PYTHONOPTIMIZE")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="compare committed data with fresh exports without writing")
    args = parser.parse_args()
    sys.dont_write_bytecode = True
    datasets = {"subsix.json": subsix(), "thresholds.json": thresholds(), "simple-example.json": simple_example()}
    for name, data in datasets.items():
        path = ROOT / "data" / name
        content = json.dumps(data, indent=2, ensure_ascii=True) + "\n"
        if args.check:
            require(path.is_file() and path.read_text(encoding="utf-8") == content, f"Outdated export: {name}")
            print(f"PASS: {name} matches exact reconstruction")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
            print(f"Wrote data/{name}")


if __name__ == "__main__":
    main()
