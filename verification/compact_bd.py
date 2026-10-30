#!/usr/bin/env python3
"""Exact small-orbit Gram certificates for the sub-six construction.

Standard library only. B uses a new denominator-24 certificate; D retains
its original bound 5.19004112. Three involutions determine all free entries.
Their invariant spaces reduce positivity to matrices of order at most six.
Run this file to verify the graph, group, edge constraints, and positive
principal minors, and print a JSON report. No floating point is used.
"""
from fractions import Fraction
from functools import reduce
from itertools import product
from math import gcd
import hashlib
import json

VECTORS = (
 (0,0,0,0,1),(0,0,0,1,1),(0,0,0,1,2),(0,0,1,2,0),
 (0,1,0,1,0),(0,1,0,2,0),(0,1,1,0,0),(0,1,1,1,1),
 (0,1,1,1,2),(0,1,2,0,0),(0,1,2,2,1),(0,1,2,2,2),
 (1,0,0,0,0),(1,0,1,0,0),(1,0,1,1,2),(1,0,2,0,0),
 (1,0,2,2,2),(1,1,1,1,0),(1,1,1,2,1),(1,1,1,2,2),
 (1,1,2,1,1),(1,1,2,2,0),(1,2,1,1,0),(1,2,1,2,1),
 (1,2,1,2,2),(1,2,2,1,1),(1,2,2,2,0),
)
CYCLES = (
 ((4,5),(6,9),(7,11),(8,10),(17,22),(18,23),(19,24),(20,25),(21,26)),
 ((0,17),(1,19),(4,23),(6,25),(7,12),(8,26),(11,16),(13,24),(14,22)),
 ((0,16),(2,15),(4,25),(5,20),(6,23),(7,22),(9,18),(11,17),(12,14)),
)
B_DATA = (
 (0,1,24),(0,2,6),(0,7,15),(0,8,22),(0,14,-9),(0,16,12),
 (0,18,18),(0,19,27),(0,20,-2),(1,3,-31),(1,4,4),(1,19,11),
 (1,21,5),(2,3,-4),(2,4,20),(2,8,27),(3,4,35),(4,6,2),
 (4,8,-20),(4,9,19),(4,23,23),(4,26,34),(8,10,6),
)
D_DATA = (
 (0,3,92395858),(0,4,145133660),(0,6,-30688165),(0,12,49988928),
 (0,13,182327057),(0,15,132541340),(0,17,84507591),(0,21,76315484),
 (1,2,3087050),(1,6,138042021),(1,8,109278319),(1,13,-87574439),
 (2,6,127958661),(2,15,-3681718),(3,8,40039826),(4,5,47798171),
 (4,10,83490663),(4,25,123407703),(8,26,-16181236),
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def permutation(cycles, n=27):
    p = list(range(n))
    used = set()
    for i, j in cycles:
        require(i != j and i not in used and j not in used, 'overlapping cycles')
        used.update((i, j))
        p[i], p[j] = j, i
    return tuple(p)


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def permutation_group(generators):
    identity = tuple(range(len(generators[0])))
    group = {identity}
    todo = [identity]
    while todo:
        p = todo.pop()
        for q in generators:
            r = compose(p, q)
            if r not in group:
                group.add(r)
                todo.append(r)
    return tuple(sorted(group))


def adjacency(vectors=VECTORS):
    b = [[int(i != j and sum(x*y for x,y in zip(u,v)) % 3 == 0)
          for j,v in enumerate(vectors)] for i,u in enumerate(vectors)]
    d = [[int(i != j and not b[i][j]) for j in range(len(b))] for i in range(len(b))]
    return b, d


def reconstruct(graph, scale, diagonal, data, generators):
    n = len(graph)
    group = permutation_group(generators)
    matrix = [[diagonal if i == j else -scale for j in range(n)] for i in range(n)]
    seen = set()
    for i, j, value in data:
        orbit = {tuple(sorted((p[i],p[j]))) for p in group}
        require(not seen.intersection(orbit), 'repeated orbit')
        for u,v in orbit:
            require(u != v and graph[u][v] == 0, 'orbit is not a free entry')
            matrix[u][v] = matrix[v][u] = value
        seen.update(orbit)
    expected = {(i,j) for i in range(n) for j in range(i+1,n) if graph[i][j] == 0}
    require(seen == expected, 'free orbits do not cover all nonedges')
    require(all(matrix[p[i]][p[j]] == matrix[i][j]
                for p in generators for i in range(n) for j in range(n)),
            'matrix is not group invariant')
    return matrix


def apply_permutation(p, v):
    out = [0] * len(v)
    for i, j in enumerate(p):
        out[j] = v[i]
    return out


def primitive(v):
    divisor = reduce(gcd, (abs(x) for x in v), 0)
    return [x // divisor for x in v] if divisor else v


def first_independent_columns(columns):
    """Keep the first linearly independent columns, using exact row reduction."""
    kept, echelon = [], {}
    for original in columns:
        v = [Fraction(x) for x in original]
        for pivot, w in sorted(echelon.items()):
            if v[pivot]:
                factor = v[pivot]
                v = [x-factor*y for x,y in zip(v,w)]
        pivot = next((i for i,x in enumerate(v) if x), None)
        if pivot is not None:
            factor = v[pivot]
            echelon[pivot] = [x/factor for x in v]
            kept.append(original)
    return kept


def projector_basis(z, c, a, b, epsilon, delta, alpha, beta=None):
    """Primitive columns of (I+beta B)(I+alpha A)(I+delta C)(I+epsilon Z).

    When epsilon=-1, omit the B factor: B interchanges the two A signs.
    Factors are unnormalised; their image is all that is required.
    """
    factors = [(epsilon,z),(delta,c),(alpha,a)]
    if beta is not None:
        factors.append((beta,b))
    columns = []
    for i in range(len(a)):
        v = [int(i == j) for j in range(len(a))]
        for sign,p in factors:
            w = apply_permutation(p,v)
            v = [x+sign*y for x,y in zip(v,w)]
        columns.append(primitive(v))
    return first_independent_columns(columns)


def gram_block(matrix, basis):
    return [[sum(u[i]*matrix[i][j]*v[j]
                 for i in range(len(matrix)) for j in range(len(matrix)))
             for v in basis] for u in basis]


def leading_minors(matrix):
    n = len(matrix)
    require(n > 0 and all(len(row) == n for row in matrix), 'bad matrix order')
    require(all(matrix[i][j] == matrix[j][i] for i in range(n) for j in range(n)),
            'not symmetric')
    work = [row[:] for row in matrix]
    previous, minors = 1, []
    for k in range(n):
        pivot = work[k][k]
        require(pivot > 0, 'nonpositive leading principal minor')
        minors.append(pivot)
        for i in range(k+1,n):
            for j in range(i,n):
                numerator = pivot*work[i][j]-work[i][k]*work[k][j]
                require(numerator % previous == 0, 'nonexact Bareiss division')
                work[i][j] = work[j][i] = numerator // previous
        previous = pivot
    return minors


def reduced_certificate(matrix, generators):
    a,b,c = generators
    z = compose(compose(a,b),compose(a,b))
    identity = tuple(range(len(a)))
    require(all(compose(g,g) == identity for g in (a,b,c,z)), 'non-involutory generator')
    require(compose(a,b) != compose(b,a), 'expected noncommuting reflections')
    require(compose(c,a) == compose(a,c) and compose(c,b) == compose(b,c), 'C not central')
    require(compose(z,a) == compose(a,z) and compose(z,b) == compose(b,z), 'Z not central')
    report, all_bases = [], []
    for epsilon,delta,alpha in product((1,-1),repeat=3):
        betas = (1,-1) if epsilon == 1 else (None,)
        for beta in betas:
            basis = projector_basis(z,c,a,b,epsilon,delta,alpha,beta)
            if not basis:
                continue
            all_bases.extend(basis)
            block = gram_block(matrix,basis)
            divisor = reduce(gcd,(abs(x) for row in block for x in row),0)
            primitive_block = [[x//divisor for x in row] for row in block]
            minors = leading_minors(primitive_block)
            if epsilon == 1 or alpha == 1:
                report.append({'signs':[epsilon,delta,alpha,beta],
                               'dimension':len(basis),'multiplicity':1 if epsilon == 1 else 2,
                               'positive_scale':divisor,'matrix':primitive_block,
                               'positive_leading_minors':minors})
    require(len(all_bases) == len(matrix), 'decomposition has wrong dimension')
    require(len(first_independent_columns(all_bases)) == len(matrix), 'decomposition is singular')
    # Positivity on all generated blocks was checked above, including the
    # partners omitted from the report. The orthogonal splitting also follows
    # algebraically from the commuting involutions on each Z eigenspace.
    require(sum(x['dimension']*x['multiplicity'] for x in report) == len(matrix),
            'reported multiplicities have wrong dimension')
    return report


def verify(vectors=VECTORS, cycles=CYCLES, b_data=B_DATA, d_data=D_DATA):
    b,d = adjacency(vectors)
    generators = tuple(permutation(c,len(vectors)) for c in cycles)
    group = permutation_group(generators)
    require(len(group) == 16, 'unexpected group order')
    require(all(b[p[i]][p[j]] == b[i][j] for p in generators
                for i in range(len(b)) for j in range(len(b))), 'not graph automorphisms')
    output = {'group_order':len(group),'certificate_type':'exact integer arithmetic',
              'B':{},'D':{}}
    for name,graph,scale,diagonal,data in (
            ('B',b,24,104,b_data),('D',d,100000000,419004112,d_data)):
        matrix = reconstruct(graph,scale,diagonal,data,generators)
        full_minors = leading_minors(matrix)
        blocks = reduced_certificate(matrix,generators)
        output[name] = {'order':len(matrix),'scale':scale,'diagonal':diagonal,
                        'capacity_bound':str(Fraction(diagonal+scale,scale)),
                        'free_pair_orbits':len(data),'full_positive_minors':len(full_minors),
                        'matrix_sha256':hashlib.sha256(json.dumps(matrix,separators=(',',':')).encode()).hexdigest(),
                        'blocks':blocks}
    return output


if __name__ == '__main__':
    print(json.dumps(verify(),indent=2))
