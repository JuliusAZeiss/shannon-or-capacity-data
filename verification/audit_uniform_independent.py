"""Independent exact audit of the uniform cover, reconstructed from the manuscript.

No import from, or reliance on the outputs of, the existing manuscript checkers.
All assertions use integers/Fraction; Decimal is used only for a labelled
diagnostic value of the optimized transcendental bound. Run with Python 3.11+.
This finite checker supplements, and does not replace, the general spectral
and integer-shift arguments proved in Supplementary Appendix G.
"""
from collections import Counter
from decimal import Decimal, localcontext
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
source = (ROOT / 'threshold_supplement.tex').read_text()
out = {'status': 'PASS', 'source_sha256': sha256(source.encode()).hexdigest()}


def pair_check(code, adjacent):
    assert len(code) == len(set(code))
    return all(any(adjacent(x, y) for x, y in zip(a, b))
               for a, b in combinations(code, 2))


def psd_rank(matrix):
    """Exact Schur elimination: zero pivots require a zero residual row."""
    a = [list(map(F, row)) for row in matrix]
    n, rank = len(a), 0
    assert all(a[i][j] == a[j][i] for i in range(n) for j in range(n))
    for k in range(n):
        p = a[k][k]
        assert p >= 0
        if not p:
            assert all(a[k][j] == 0 for j in range(k + 1, n))
            continue
        rank += 1
        for i in range(k + 1, n):
            for j in range(i, n):
                a[i][j] -= a[i][k] * a[k][j] / p
                a[j][i] = a[i][j]
    return rank


def binary_rank(matrix):
    a = [sum(int(v) << i for i, v in enumerate(row)) for row in matrix]
    rank = 0
    for col in range(len(matrix[0])):
        pivot = next((j for j in range(rank, len(a)) if a[j] >> col & 1), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        for j in range(len(a)):
            if j != rank and a[j] >> col & 1:
                a[j] ^= a[rank]
        rank += 1
    return rank


def dot(x, y):
    return (x & y).bit_count() % 2


# A=U6, D=complement: dot matrix, six-clique, transvections and theta Gram.
vertices = [x for x in range(64) if x.bit_count() % 2]
adj_A = lambda x, y: x != y and not dot(x, y)
adj_D = lambda x, y: x != y and bool(dot(x, y))
M = [[dot(x, y) for y in vertices] for x in vertices]
assert binary_rank(M) == 6
basis = [1 << i for i in range(6)]
assert all(adj_A(x, y) for x, y in combinations(basis, 2))
maps = []
for v in range(64):
    if v.bit_count() % 2:
        continue
    mapping = {x: x ^ (v if dot(x, v) else 0) for x in vertices}
    assert set(mapping.values()) == set(vertices)
    assert all(mapping[mapping[x]] == x for x in vertices)
    assert all(dot(mapping[x], mapping[y]) == dot(x, y)
               for x in vertices for y in vertices)
    maps.append(mapping)
orbit = {1}
while True:
    enlarged = orbit | {mapping[x] for x in orbit for mapping in maps}
    if enlarged == orbit:
        break
    orbit = enlarged
assert orbit == set(vertices)
gram_D = [[F(int(x == y) - int(y == (x ^ 63)))
           + F(1 + (-1)**dot(x, y), 4) for y in vertices] for x in vertices]
assert all(gram_D[i][i] == 1 for i in range(32))
assert all(not gram_D[i][j] for i, x in enumerate(vertices)
           for j, y in enumerate(vertices) if adj_D(x, y))
augmented = [[F(1)] + [F(1, 2)] * 32] + [
    [F(1, 2)] + row for row in gram_D]
rank_D = psd_rank(augmented)
clique_D = [32 ^ v for v in (0, 3, 12, 15)]
assert all(adj_D(x, y) for x, y in combinations(clique_D, 2))
out['binary_orthogonality'] = {
    'vertices': 32, 'binary_fitting_rank': 6,
    'transvection_orbit': len(orbit), 'complement_clique': 4,
    'theta_augmented_gram_rank': rank_D, 'squared_handle_overlap': '1/4',
    'fractional_coloring_of_D': '32/6 = 16/3'}

# R37: check complementing automorphism on every pair, codes and six-sets.
S = {6, 7, 8, 9, 10, 12, 13, 15, 18}
distance = lambda a, b: min((a-b) % 37, (b-a) % 37)
adj_R = lambda a, b: a != b and distance(a, b) in S
assert all(adj_R(6*x % 37, 6*y % 37) != adj_R(x, y)
           for x, y in combinations(range(37), 2))
assert pair_check([(x, 6*x % 37) for x in range(37)], adj_R)
B = [0, 6, 12, 18, 24, 30]
assert all(adj_R(x, y) for x, y in combinations(B, 2))
I = [6*x % 37 for x in B]
assert not any(adj_R(x, y) for x, y in combinations(I, 2))
assert Counter((i+j) % 37 for i in I for j in range(37)) == Counter({x: 6 for x in range(37)})
assert F('6.08')**2 < 37 < F(73, 12)**2
out['circulant37'] = {'self_complementing_multiplier': 6, 'square_code': 37,
                      'clique_size': 6, 'independent_size': 6,
                      'fractional_coloring': '37/6'}

# Clebsch cube code, constructed directly in F2[z]/(z4+z3+z2+z+1).
def mul16(a, b):
    ans = 0
    while b:
        if b & 1:
            ans ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 31
    return ans

H = {1, 2, 4, 8, 15}
assert {mul16(2, x) for x in H} == H
assert mul16(2, 15) == 1
a = 3
assert a not in H
a2 = mul16(a, a)
cosets = [H, {mul16(a, x) for x in H}, {mul16(a2, x) for x in H}]
assert set.union(*cosets) == set(range(1, 16))
assert sum(map(len, cosets)) == len(set.union(*cosets))
adj_Q = lambda x, y: x != y and (x ^ y) in H
code_Q = [(x, mul16(a, x), mul16(a2, x)) for x in range(16)]
assert pair_check(code_Q, adj_Q)
gram_Q = [[F(int(x == y)) + F(int(x != y and not adj_Q(x, y)), 2)
           for y in range(16)] for x in range(16)]
assert all(sum(row) == 6 for row in gram_Q)
rank_Q = psd_rank(gram_Q)
out['clebsch'] = {'cube_code': 16, 'theta_gram_rank': rank_Q,
                  'theta_gram_row_sum': 6, 'squared_handle_overlap': str(F(6, 16))}

# Pentagonal square and complement-of-nine-cycle cube codes.
adj_C5 = lambda x, y: (x-y) % 5 in (1, 4)
assert pair_check([(i, 2*i % 5) for i in range(5)], adj_C5)
adj_Z = lambda x, y: (x-y) % 9 not in (0, 1, 8)
code_Z = [(x, y, (2*x + 4*y) % 9) for x in range(9) for y in range(9)]
assert pair_check(code_Z, adj_Z)
assert F('2.235')**2 < 5 < F('2.24')**2
assert F('2.5')**3 < 16
assert F('4.325')**3 < 81
assert F(9)*F(205, 216)/(1 + F(205, 216)) == F(1845, 421) < F(22, 5)
out['cycle_codes'] = {'pentagon_square': 5, 'complement9_cube': 81,
                       'theta9_rational_enclosure': '1845/421 < 22/5'}

# Parse the printed table; derive projection bounds independently from graph data.
table = source.split('\\midrule', 1)[1].split('\\bottomrule', 1)[0]
rows = []
for line in table.splitlines():
    if '&' not in line:
        continue
    fields = [field.strip().strip('$') for field in line.rstrip('\\').strip().split('&')]
    rows.append((fields[0], fields[1], *map(F, fields[2:])))
assert len(rows) == 8
hD, hC5, hZ, hR = F(16, 3), F(5, 2), F(9, 2), F(37, 6)
qD, qC5, qZ, qQ, qR = F(4), F('2.24'), F(22, 5), F(8, 3), F(73, 12)
expected_projection = [
    (qR, hR, qD+1),
    (4+qC5, min(hD+1, 4+hC5), qD+1),
    (4+qC5, min(hD+2, 4+hC5), qD+2),
    (2+qZ, min(hD+2, 2+hZ), qD+2),
    (2+2*qC5, min(hD+2, 2+2*hC5), qD+2),
    (3*qC5, min(hD+2, 3*hC5), qD+2),
    (F(7), F(7), qD+qC5),
    (1+qR, 1+hR, qD+qQ),
]
b_lower = [F(1), F(1), F(2), F(2), F(2), F(2), F('2.235'), F('2.5')]
# Every bound below is strict, except the complete graph K7 in row 7.
partner_lower = [F('6.08'), 4+F('2.235'), 4+F('2.235'), 2+F('4.325'),
                 2+2*F('2.235'), 3*F('2.235'), F(7), 1+F('6.08')]
row_results = []
for idx, ((J, Hname, s, a, d, U, L), projection, b, hp) in enumerate(
        zip(rows, expected_projection, b_lower, partner_lower), 1):
    assert (s, a, d) == projection
    assert a > 6 and s > d
    lam = (s-d)/(a-6+s-d)
    assert 0 < lam < 1
    assert U == (a*s-6*d)/(a-6+s-d)
    assert 6*lam+s*(1-lam) == a*lam+d*(1-lam) == U
    assert U < L and L > b and L*(L-b) < 32
    assert L <= hp if idx != 7 else L < hp
    row_results.append({'row': idx, 'J': J, 'H': Hname,
                        'lambda': str(lam), 'U': str(U), 'L': str(L),
                        'interval_width': str(L-U),
                        'first_factor_quadratic_margin': str(32-L*(L-b))})
out['table_rows'] = row_results
overlaps = [rows[i][6] - rows[i+1][5] for i in range(7)]
assert all(x > 0 for x in overlaps)
assert all(rows[i][5] < rows[i+1][5] and rows[i][6] < rows[i+1][6] for i in range(7))
assert min(overlaps) == F(1, 1480)
c, endpoint = rows[0][5], rows[-1][6]
assert c == F(541, 90) and endpoint == F('7.02')
assert endpoint-c-1 == F(2, 225) > 0
out['coverage'] = {'successive_overlaps': list(map(str, overlaps)),
                   'smallest_overlap': str(min(overlaps)),
                   'unit_shift_overlap': str(endpoint-c-1),
                   'base_open_interval': [str(c), str(endpoint)],
                   'integer_shift_rule': 'j=ceil(t-c)-1; 0<t-c-j<=1'}

# Exact rational root enclosures for both branches at lambda=93/106.
def root_enclosure(num, den, power, digits=14):
    scale = 10**digits
    lo, hi = 6*scale, 7*scale
    assert lo**power * den < num * scale**power < hi**power * den
    while hi-lo > 1:
        mid = (hi+lo)//2
        if mid**power * den < num * scale**power:
            lo = mid
        else:
            hi = mid
    assert lo**power * den < num * scale**power < hi**power * den
    return F(lo, scale), F(hi, scale)

cutoff = F('6.0101')
assert 6**186 * 37**13 * 10000**212 < 60101**212
assert 37**93 * 5**13 * 10000**106 < 6**93 * 60101**106
branch1 = root_enclosure(6**186 * 37**13, 1, 212)
branch2 = root_enclosure(37**93 * 5**13, 6**93, 106)
assert branch1[1] < cutoff and branch2[1] < cutoff
assert cutoff < c and c*c < 37
# rho_1>sqrt37: z(z-1)<32 at z=sqrt37 iff sqrt37>5.
assert 37 > 25
with localcontext() as ctx:
    ctx.prec = 80
    d = (Decimal(37).sqrt()/6).ln()
    ell = (Decimal(6)/5).ln()
    lam = (ell+d)/(ell+3*d)
    T0 = 6*(2*d*d/(ell+3*d)).exp()
    diagnostic = {'T0': str(T0), 'lambda_star': str(lam),
                  'label': 'Decimal diagnostic only; assertions above certify the closed cutoff.'}
out['optimized_cutoff'] = {'rational_weight': '93/106',
                           'branch1_exact_root_enclosure': list(map(str, branch1)),
                           'branch2_exact_root_enclosure': list(map(str, branch2)),
                           'closed_cutoff': str(cutoff), 'decimal_diagnostic': diagnostic}
out['limits'] = [
    'Finite arithmetic and graph checks are not a standalone proof of spectral duality/interpolation.',
    'The all-real coverage follows from the exact overlaps plus the symbolic integer-shift argument.',
    'No exact noninteger factor capacity is inferred from a separating threshold.',
    'No claim is made at the optimized open endpoint T0.'
]
(HERE / 'audit_uniform_independent.json').write_text(json.dumps(out, indent=2) + '\n')
print('PASS: independent primitive graph, eight-row, overlap, and closed-cutoff checks.')
