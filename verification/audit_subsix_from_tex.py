#!/usr/bin/env python3
"""Independently audit the sub-six certificate data printed in the manuscript.

Only the Python standard library is used. This checker parses the current TeX
as inert text, reconstructs the certificates, and verifies them by rational
LDL^T factorization. It neither imports the embedded-data checker nor executes
any text read from the manuscript. Run from any working directory.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import re
import sys


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def digest(value):
    data = json.dumps(value, separators=(",", ":")).encode("ascii")
    return sha256(data).hexdigest()


def between(text, first, last):
    return text.split(first, 1)[1].split(last, 1)[0]


def rational_ldlt(matrix):
    """Exact LDL^T; check positivity and reconstruct every original entry."""
    n = len(matrix)
    lower = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    diagonal = []
    for j in range(n):
        pivot = Fraction(matrix[j][j]) - sum(
            lower[j][k] ** 2 * diagonal[k] for k in range(j)
        )
        require(pivot > 0, f"LDL pivot {j} is not positive")
        diagonal.append(pivot)
        for i in range(j + 1, n):
            lower[i][j] = (
                Fraction(matrix[i][j])
                - sum(lower[i][k] * lower[j][k] * diagonal[k] for k in range(j))
            ) / pivot
    for i in range(n):
        for j in range(n):
            reconstructed = sum(
                lower[i][k] * diagonal[k] * lower[j][k] for k in range(n)
            )
            require(reconstructed == matrix[i][j], f"LDL reconstruction at {i},{j}")
    principal = Fraction(1)
    principal_hashes = []
    for j, pivot in enumerate(diagonal):
        principal *= pivot
        require(principal.denominator == 1, f"Noninteger leading determinant {j + 1}")
        principal_hashes.append(sha256(str(principal.numerator).encode()).hexdigest())
    return {
        "method": "rational LDL^T with full exact matrix reconstruction",
        "positive_pivots": len(diagonal),
        "pivot_fractions": [str(x) for x in diagonal],
        "leading_principal_determinant_sha256": principal_hashes,
        "determinant": str(principal.numerator),
        "reconstructed_entries": n * n,
    }


def rank_mod_three(matrix):
    a = [[entry % 3 for entry in row] for row in matrix]
    row = 0
    for col in range(len(a[0])):
        candidates = [i for i in range(row, len(a)) if a[i][col]]
        if not candidates:
            continue
        selected = candidates[0]
        a[row], a[selected] = a[selected], a[row]
        inverse = pow(a[row][col], -1, 3)
        a[row] = [(inverse * entry) % 3 for entry in a[row]]
        for i in range(len(a)):
            if i != row and a[i][col]:
                multiplier = a[i][col]
                a[i] = [(x - multiplier * y) % 3 for x, y in zip(a[i], a[row])]
        row += 1
        if row == len(a):
            break
    return row


def graph_record(labels, adjacency):
    n = len(labels)
    require(len(set(labels)) == n, "Duplicate graph labels")
    require(all(adjacency[i][i] == 0 for i in range(n)), "Graph loop")
    require(all(adjacency[i][j] == adjacency[j][i]
                for i in range(n) for j in range(n)), "Asymmetric adjacency")
    return {
        "vertex_labels": labels,
        "order": n,
        "size": sum(map(sum, adjacency)) // 2,
        "adjacency_sha256": digest(adjacency),
    }


def main():
    source = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        Path(__file__).resolve().parent.parent / "subsix_supplement.tex"
    )
    output = Path(sys.argv[2]) if len(sys.argv) > 2 else (
        Path(__file__).resolve().parent / "audit_subsix_from_tex.json"
    )
    tex = source.read_text(encoding="utf-8")

    # Ordered F_3^5 vectors are read from the displayed table itself.
    vector_text = between(tex, r"\(i\)&\(w_i\)", r"\end{tabular}")
    vector_rows = re.findall(r"(\d+)&\(([012]),([012]),([012]),([012]),([012])\)", vector_text)
    vector_by_label = {int(row[0]): list(map(int, row[1:])) for row in vector_rows}
    require(len(vector_rows) == len(vector_by_label) == 27, "Vector row count")
    require(sorted(vector_by_label) == list(range(27)), "Vector labels")
    vectors = [vector_by_label[i] for i in range(27)]
    require(len({tuple(x) for x in vectors}) == 27, "Repeated vectors")
    dot = lambda u, v: sum(a * b for a, b in zip(u, v)) % 3
    norms = [dot(w, w) for w in vectors]
    require(all(norms), "An isotropic vector cannot be normalized")
    b_adj = [[int(i != j and dot(vectors[i], vectors[j]) == 0)
              for j in range(27)] for i in range(27)]
    d_adj = [[int(i != j and not b_adj[i][j]) for j in range(27)] for i in range(27)]
    scalar_fit = [[pow(norms[i], -1, 3) * dot(vectors[i], vectors[j]) % 3
                   for j in range(27)] for i in range(27)]
    require(all(scalar_fit[i][i] == 1 for i in range(27)), "Ternary fit diagonal")
    require(all(scalar_fit[i][j] == 0 for i in range(27) for j in range(27)
                if b_adj[i][j]), "Ternary fit edge zero")
    scalar_rank = rank_mod_three(scalar_fit)
    require(scalar_rank <= 5, "Ternary fit exceeds stated rank")
    g_adj = [[0 for _ in range(54)] for _ in range(54)]
    for i in range(27):
        for j in range(27):
            g_adj[i][j] = b_adj[i][j]
            g_adj[i + 27][j + 27] = d_adj[i][j]
    g_formula = re.search(
        r"\\mathcal K_G=\\\{\(i,(\d+)\+i\):0\\leq i\\leq(\d+)\\\}", tex
    )
    require(g_formula is not None, "Cannot parse the displayed G code formula")
    g_offset, g_last = map(int, g_formula.groups())
    g_code = [(i, i + g_offset) for i in range(g_last + 1)]
    require(len(g_code) == len(set(g_code)) == 27, "G code word count")
    require(all(0 <= x < 54 and 0 <= y < 54 for x, y in g_code), "G code outside vertex set")
    g_pairs = list(combinations(g_code, 2))
    require(all(any(g_adj[x][y] for x, y in zip(u, v)) for u, v in g_pairs), "G code pair")

    # Derive H's native residue labels and exceptional codeword from the TeX.
    h_vertex_formula = re.search(
        r"V\(H\)=\\mathbb Z_\{(\d+)\}\\setminus\\\{([\d,]+)\\\}", tex
    )
    require(h_vertex_formula is not None, "Cannot parse H vertex set")
    h_vertex_modulus = int(h_vertex_formula[1])
    h_vertex_excluded = set(map(int, h_vertex_formula[2].split(",")))
    expected_labels = [i for i in range(h_vertex_modulus) if i not in h_vertex_excluded]
    require(h_vertex_modulus == 29 and h_vertex_excluded == {1,17}, "H residue vertex set")
    h_text = between(tex, r"\begin{longtable}", r"\end{longtable}")
    h_rows = re.findall(r"^\s*(\d+)&\\\((.*?)\\\)\\\\", h_text, re.M)
    labels = [int(label) for label, _ in h_rows]
    require(labels == expected_labels and len(labels) == 27, "H native vertex order")
    position = {label: i for i, label in enumerate(labels)}
    h_adj = [[0 for _ in labels] for _ in labels]
    for label_text, neighbours_text in h_rows:
        label = int(label_text)
        neighbours = [int(x) for x in re.findall(r"\d+", neighbours_text)]
        require(neighbours == sorted(set(neighbours)), f"H row {label} order")
        for neighbour in neighbours:
            require(neighbour in position and neighbour > label, f"Invalid H edge {label},{neighbour}")
            i, j = position[label], position[neighbour]
            h_adj[i][j] = h_adj[j][i] = 1
    require(sum(map(sum, h_adj)) // 2 == 176, "H edge count")
    def h_edge(v, w):
        return h_adj[position[v]][position[w]]
    h_formula_text = re.sub(r"\s+", "", between(tex, r"\mathcal K_H=", r"\]"))
    h_formula = re.fullmatch(
        r"\\\{\(i,(\d+)i\\bmod(\d+)\):i\\in\\mathbbZ_\{(\d+)\}"
        r"\\setminus\\\{([\d,]+)\\\}\\\}\\cup\\\{\((\d+),(\d+)\)\\\}",
        h_formula_text
    )
    require(h_formula is not None, "Cannot parse the displayed H code formula")
    h_multiplier, h_modulus, h_ring = map(int, h_formula.groups()[:3])
    h_excluded = set(map(int, h_formula[4].split(",")))
    h_extra = [tuple(map(int, h_formula.groups()[4:]))]
    require(h_modulus == h_ring == h_vertex_modulus, "H code modulus and domain disagree")
    h_code = [(i, h_multiplier * i % h_modulus) for i in range(h_ring) if i not in h_excluded]
    h_code.extend(h_extra)
    require(len(h_code) == len(set(h_code)) == 27, "H code word count")
    require(all(x in position and y in position for x, y in h_code), "H code outside vertex set")
    require(len({x for x, _ in h_code}) == len(h_code), "H code first coordinates repeat")
    require([i for i in range(h_modulus) if h_multiplier * i % h_modulus == 1] == [17], "First modular inverse")
    require([i for i in range(h_modulus) if h_multiplier * i % h_modulus == 17] == [28], "Second modular inverse")
    h_pairs = list(combinations(h_code, 2))
    require(all(h_edge(u[0], v[0]) or h_edge(u[1], v[1]) for u, v in h_pairs), "H code pair")

    # Read all eleven sets, then verify both orientations of every block zero.
    colour_text = between(tex, r"\(a\)&Independent set", r"\end{tabular}")
    colour_rows = re.findall(r"^\s*(\d+)&\\\{(.*?)\\\}\\\\", colour_text, re.M)
    colour_sets = {int(a): list(map(int, re.findall(r"\d+", entries)))
                   for a, entries in colour_rows}
    require(sorted(colour_sets) == list(range(1, 12)), "Colour set labels")
    memberships = {h: [a for a in colour_sets if h in colour_sets[a]] for h in labels}
    for a, vertices in colour_sets.items():
        require(len(vertices) == len(set(vertices)), f"Colour set {a} repeats a vertex")
        require(all(v in position for v in vertices), f"Colour set {a} has invalid vertex")
        require(all(not h_edge(u, v) for u, v in combinations(vertices, 2)), f"Colour set {a} is not independent")
    require(all(len(indices) == 2 for indices in memberships.values()), "Two-fold coverage")
    block_pairs = 0
    for h in labels:
        require(memberships[h][0] != memberships[h][1], "Block diagonal not identity")
        for k in labels:
            if h_edge(h, k):
                require(all(a != b for a in memberships[h] for b in memberships[k]), "Block edge entry")
                block_pairs += 1
    # Construct the 54-by-54 block-fitting matrix from the colour incidences.
    block_columns = [a for h in labels for a in memberships[h]]
    h_fit = [[int(a == b) for b in block_columns] for a in block_columns]
    for i in range(len(labels)):
        for j in range(len(labels)):
            if i == j or h_adj[i][j]:
                for a in range(2):
                    for b in range(2):
                        require(h_fit[2*i+a][2*j+b] == int(i == j and a == b), "H block fitting entry")
    h_block_rank = rank_mod_three(h_fit)
    require(h_block_rank <= 11, "H block-fitting rank exceeds eleven")

    # Parse scalar definitions and the normalization table independently.
    parameters = {}
    for name in ("b", "q", "L"):
        match = re.search(r"\b" + name + r"=\\frac\{(\d+)\}\{(\d+)\}=(\d+\.\d+)", tex)
        require(match is not None, f"Missing parameter {name}")
        parameters[name] = Fraction(*map(int, match.groups()[:2]))
        require(parameters[name] == Fraction(match[3]), f"Printed decimal for {name}")
    table_text = between(tex, r"Graph&Order&Scale", r"\end{tabular}")
    normalization = {}
    for name, order, scale, diagonal, parameter in re.findall(
        r"\$(B|D|H)\$&(\d+)&\$(\d+)\$&\$(\d+)\$&\$(.*?)\$", table_text
    ):
        t = Fraction(parameter) if parameter not in parameters else parameters[parameter]
        normalization[name] = (int(order), int(scale), int(diagonal), t)
    require(set(normalization) == {"B", "D", "H"}, "Normalization table")

    # The source symmetries are read from the actual displayed cycle notation.
    cycle_text = between(tex, r"\begin{align*}", r"\end{align*}")
    permutations = []
    for name in ("a", "b_0", "c"):
        line = re.search(r"^\s*" + name + r"&=(.*)$", cycle_text, re.M)
        require(line is not None, f"Missing permutation {name}")
        pairs = [tuple(map(int, re.findall(r"\d+", pair)))
                 for pair in re.findall(r"\(([^()]*)\)", line[1])]
        require(len(pairs) == 9 and all(len(pair) == 2 for pair in pairs), "Cycle shape")
        perm = list(range(27)); moved = set()
        for i, j in pairs:
            require(0 <= i < 27 and 0 <= j < 27 and i != j
                    and i not in moved and j not in moved, "Overlapping permutation cycles")
            perm[i], perm[j] = j, i; moved.update((i, j))
        require(all(b_adj[i][j] == b_adj[perm[i]][perm[j]]
                    for i in range(27) for j in range(27)), "Not a graph automorphism")
        permutations.append(tuple(perm))
    group = {tuple(range(27))}; todo = list(group)
    while todo:
        g = todo.pop()
        for h in permutations:
            gh = tuple(g[h[i]] for i in range(27))
            if gh not in group:
                group.add(gh); todo.append(gh)
    require(len(group) == 16, "Wrong generated group order")

    def printed_orbits(name):
        block = between(tex, "% BEGIN COMPACT " + name, "% END COMPACT " + name)
        block = between(block, r"\begin{Verbatim}", r"\end{Verbatim}").split("\n", 1)[1]
        rows, current = [], None
        for line in block.splitlines():
            if not line.strip(): continue
            head = re.match(r"\s*(\d+):\s*(.*)", line)
            if head:
                current, values = int(head[1]), head[2]
            else:
                require(current is not None, "Orphan continuation")
                values = line.strip()
            require(re.fullmatch(r"\d+=-?\d+(?:\s+\d+=-?\d+)*", values), "Unrecognized orbit data")
            rows.extend((current, int(j), int(value)) for j, value in re.findall(r"(\d+)=(-?\d+)", values))
        return rows

    adjacency = {"B": b_adj, "D": d_adj, "H": h_adj}
    vertex_orders = {"B": list(range(27)), "D": list(range(27)), "H": labels}
    sigma_rule = re.search(
        r"On \$H\$, let \$\\sigma\$ fix \$([\d,]+)\$ and send every other vertex\s*"
        r"\$v\$ to \$(\d+)-v\$", tex
    )
    require(sigma_rule is not None, "Cannot parse the printed H symmetry rule")
    sigma_fixed = set(map(int, sigma_rule[1].split(",")))
    sigma_modulus = int(sigma_rule[2])
    sigma = {x: x if x in sigma_fixed else sigma_modulus-x for x in labels}
    require(set(sigma.values()) == set(labels), "H symmetry does not permute its vertices")
    require(all(sigma[sigma[x]] == x for x in labels), "H symmetry is not an involution")
    require(all(h_edge(x, y) == h_edge(sigma[x], sigma[y])
                for x in labels for y in labels), "H involution does not preserve adjacency")
    matrix_data, matrix_audits, orbit_data = {}, {}, {}
    for name in ("B", "D", "H"):
        n, scale, diagonal, parameter = normalization[name]
        order = vertex_orders[name]
        pos = {x: i for i, x in enumerate(order)}
        graph = adjacency[name]
        require(n == len(order), "Normalization order mismatch")
        entries = printed_orbits(name)
        require(len(entries) == {"B":23, "D":19, "H":91}[name], "Orbit entry count")
        orbit_data[name] = entries
        matrix = [[diagonal if i == j else -scale if graph[i][j] else None
                   for j in range(n)] for i in range(n)]
        seen = set()
        for x, y, value in entries:
            orbit = ({tuple(sorted((g[x], g[y]))) for g in group} if name != "H"
                     else {(x,y), tuple(sorted((sigma[x], sigma[y])))})
            require((x,y) == min(orbit) and x < y, "Noncanonical orbit representative")
            require(not seen.intersection(orbit), "Repeated free orbit")
            for a, b in orbit:
                i, j = pos[a], pos[b]
                require(i != j and graph[i][j] == 0, "Orbit on an edge")
                matrix[i][j] = matrix[j][i] = value
            seen.update(orbit)
        require(all(x is not None for row in matrix for x in row), "Incomplete matrix")
        require(Fraction(diagonal, scale) == parameter - 1, "Matrix normalization")
        require(all(matrix[i][j] == -scale for i in range(n) for j in range(n)
                    if graph[i][j]), "Wrong edge entry")
        matrix_data[name] = matrix
        matrix_audits[name] = {"order": n, "scale": scale, "diagonal": diagonal,
                              "parameter": str(parameter), "free_pair_orbits": len(entries),
                              "matrix_sha256": digest(matrix), **rational_ldlt(matrix)}
    # Independently reconstruct the displayed small B-block minor table.
    from math import gcd
    from functools import reduce
    def permute(g, vector):
        result = [0]*27
        for i in range(27): result[g[i]] = vector[i]
        return result
    def compose(g, h): return tuple(g[h[i]] for i in range(27))
    aa, bb, cc = permutations
    zz = compose(compose(aa,bb), compose(aa,bb))
    ident = tuple(range(27))
    require(all(compose(g,g) == ident for g in (aa,bb,cc,zz)), "Involution relation")
    require(compose(aa,bb) == compose(zz,compose(bb,aa)), "AB=ZBA relation")
    require(all(compose(g,h)==compose(h,g) for g in (cc,zz) for h in (aa,bb)), "Centrality")
    def independent_columns(columns):
        pivots, basis = {}, []
        for column in columns:
            w = list(map(Fraction,column))
            for i in sorted(pivots):
                factor = w[i]
                w = [x-factor*y for x,y in zip(w,pivots[i])]
            i = next((i for i,x in enumerate(w) if x),None)
            if i is not None:
                pivots[i] = [x/w[i] for x in w]
                basis.append(column)
        return basis
    block_table = between(tex, r"$\varepsilon$&$\delta$", r"\end{tabular}")
    block_rows = re.findall(r"^\$([+-])\$&\$([+-])\$&\$([+-])\$&(\$[+-]\$|--)&(\d+)&(\d+)&\$([\d,]+)\$", block_table, re.M)
    require(len(block_rows)==8,"Printed small-block row count")
    block_record=[]
    for e,d,a,b,order,multiplicity,minor_text in block_rows:
        signs=[1 if x=='+' else -1 for x in (e,d,a)]
        beta=None if b=='--' else 1 if b=='$+$' else -1
        eps,delta,alpha=signs
        require((eps==1)==(beta is not None),"Wrong projector case")
        factors=[(eps,zz),(delta,cc),(alpha,aa)] + ([] if beta is None else [(beta,bb)])
        columns=[]
        for i in range(27):
            column=[int(j==i) for j in range(27)]
            for sign,perm in factors:
                shifted=permute(perm,column)
                column=[x+sign*y for x,y in zip(column,shifted)]
            divisor=reduce(gcd,(abs(x) for x in column),0)
            if divisor:columns.append([x//divisor for x in column])
        basis=independent_columns(columns)
        require(len(basis)==int(order),"Printed block order")
        require(int(multiplicity)==(1 if eps==1 else 2),"Printed block multiplicity")
        raw=[[sum(u[i]*matrix_data['B'][i][j]*w[j] for i in range(27) for j in range(27))
              for w in basis] for u in basis]
        divisor=reduce(gcd,(abs(x) for row in raw for x in row),0)
        matrix=[[x//divisor for x in row] for row in raw]
        record=rational_ldlt(matrix);minor=Fraction(1);minors=[]
        for pivot in record['pivot_fractions']:
            minor*=Fraction(pivot);require(minor.denominator==1,"Noninteger B-block minor")
            minors.append(minor.numerator)
        require(minors==list(map(int,minor_text.split(','))),"Printed B-block minors")
        block_record.append({'signs':signs+[beta],'order':len(basis),'multiplicity':int(multiplicity),'minors':minors})
    require(sum(x['order']*x['multiplicity'] for x in block_record)==27,"Block dimension sum")

    b, q, bound = (parameters[key] for key in ("b", "q", "L"))
    weight = re.search(r"\\lambda=(\d+)/(\d+)", tex)
    require(weight is not None, "Missing rational interpolation weight")
    numerator, denominator = map(int, weight.groups())
    remainder = denominator - numerator
    require(0 < numerator < denominator, "Interpolation weight outside (0,1)")
    power_text = re.sub(r"\s+", "", between(
        tex, r"\begin{align}", r"\end{align}"
    ))
    power_one = re.search(
        r"(\d+)\^\{(\d+)\}\\left\(\\frac\{(\d+)\}\{(\d+)\}\\right\)\^\{(\d+)\}"
        r"&<\\left\(\\frac\{(\d+)\}\{(\d+)\}\\right\)\^\{(\d+)\},\\label\{subsix:eq:power1\}",
        power_text
    )
    power_two = re.search(
        r"\\left\(\\frac\{(\d+)\}\{(\d+)\}\\right\)\^\{(\d+)\}"
        r"\\left\(\\frac\{(\d+)\}\{(\d+)\}\\right\)\^\{(\d+)\}"
        r"&<\\left\(\\frac\{(\d+)\}\{(\d+)\}\\right\)\^\{(\d+)\}\.\\label\{subsix:eq:power2\}",
        power_text
    )
    require(power_one is not None and power_two is not None, "Cannot parse both displayed power inequalities")
    p1, p2 = list(map(int, power_one.groups())), list(map(int, power_two.groups()))
    require(p1[0] == scalar_rank == 5 and p1[1] == numerator
            and Fraction(p1[2], p1[3]) == q and p1[4] == remainder
            and Fraction(p1[5], p1[6]) == bound and p1[7] == denominator,
            "First displayed power inequality differs from interpolation data")
    require(Fraction(p2[0], p2[1]) == Fraction(h_block_rank, 2) == Fraction(11, 2)
            and p2[2] == numerator and Fraction(p2[3], p2[4]) == b
            and p2[5] == remainder and Fraction(p2[6], p2[7]) == bound
            and p2[8] == denominator,
            "Second displayed power inequality differs from interpolation data")
    square = re.search(r"(\d+)\^2=(\d+)<(\d+)=(\d+)\\cdot10\^\{(\d+)\}", tex)
    require(square is not None, "Cannot parse displayed integer-square comparison")
    sq_base, sq_value, sq_upper, sq_factor, sq_exp = map(int, square.groups())
    comparisons = {
        "five_profile_strictly_below_L": Fraction(5) ** numerator * q ** remainder < bound ** denominator,
        "eleven_halves_profile_strictly_below_L": Fraction(11, 2) ** numerator * b ** remainder < bound ** denominator,
        "L_squared_less_than_27": bound ** 2 < 27,
        "q_squared_greater_than_27": q ** 2 > 27,
        "b_less_than_16_over_3": b < Fraction(16, 3),
        "factor_upper_bounds_below_six": max(Fraction(16, 3), b, q) < 6,
        "printed_integer_square": (sq_base ** 2 == sq_value < sq_upper
                                   and sq_upper == sq_factor * 10 ** sq_exp
                                   and sq_base == bound.numerator
                                   and 10 ** sq_exp == bound.denominator ** 2
                                   and sq_factor == len(h_code)),
    }
    require(all(comparisons.values()), "A strict scalar comparison failed")

    report = {
        "status": "PASS",
        "input": source.name,
        "input_sha256": sha256(source.read_bytes()).hexdigest(),
        "arithmetic": "integers and fractions.Fraction only; no floating point",
        "independence": "All finite certificate data parsed directly from TeX; no checker imported",
        "graphs": {name: graph_record(vertex_orders[name], adjacency[name]) for name in ("B", "D", "H")},
        "G": graph_record(list(range(54)), g_adj),
        "vectors_sha256": digest(vectors),
        "nonisotropic_vectors": len(norms),
        "ternary_scalar_fit": {"rank": scalar_rank, "diagonal_entries": 27,
                               "oriented_edge_zeros": sum(map(sum, b_adj)), "sha256": digest(scalar_fit)},
        "codes": {"G_words": len(g_code), "G_pair_tests": len(g_pairs),
                  "H_words": len(h_code), "H_pair_tests": len(h_pairs),
                  "parsed_G_formula": {"offset": g_offset, "last_index": g_last},
                  "parsed_H_formula": {"multiplier": h_multiplier, "modulus": h_modulus,
                                       "excluded_domain_labels": sorted(h_excluded),
                                       "exceptional_words": h_extra},
                  "G_sha256": digest(g_code), "H_sha256": digest(h_code)},
        "colouring": {"sets": colour_sets, "memberships": memberships,
                      "oriented_block_edge_tests": block_pairs,
                      "scalar_edge_entry_tests": 4 * block_pairs,
                      "F3_block_fitting_rank": h_block_rank,
                      "block_fitting_matrix_sha256": digest(h_fit),
                      "sha256": digest(list(colour_sets.values()))},
        "matrices": matrix_audits,
        "printed_B_block_table": block_record,
        "symmetry": {"group_order": len(group), "source_generators": permutations, "orbit_data": orbit_data, "H_involution": sigma, "H_vertex_labels": labels},
        "parameters": {key: str(value) for key, value in parameters.items()},
        "interpolation": {"numerator": numerator, "denominator": denominator,
                          "complementary_numerator": remainder,
                          "proof_method": "existing global spectral interpolation inequality"},
        "scalar_comparisons": comparisons,
        "scope": "Verifies the printed finite data and exact scalar comparisons; the spectral interpolation and symmetry-decomposition arguments require mathematical review.",
    }
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "input_sha256": report["input_sha256"],
                      "graphs": {name: {k: v for k, v in graph.items() if k in ("order", "size")}
                                 for name, graph in report["graphs"].items()},
                      "positive_LDL_pivots": {name: item["positive_pivots"] for name, item in matrix_audits.items()},
                      "scalar_comparisons": comparisons, "report": str(output)}, indent=2))


if __name__ == "__main__":
    main()
