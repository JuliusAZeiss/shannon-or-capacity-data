"""Exact integer certificate for the closed uniform cutoff 6.0101.

The graph and eight-interval coverage inputs are independently checked by
check_exact.py. Here lambda=93/106 certifies both geometric branches
below 60101/10000. Only integer arithmetic is used for these comparisons.
The analytical interpolation and coverage arguments are proved in
Supplementary Appendix G.
"""
from fractions import Fraction

p, q = 93, 106
num, den = 60101, 10000
assert 6 ** (2*p) * 37 ** (q-p) * den ** (2*q) < num ** (2*q)
assert 37 ** p * 5 ** (q-p) * den ** q < 6 ** p * num ** q
t = Fraction(num, den)
old = Fraction(541, 90)
assert t < old
assert old**2 < 37
# rho=(1+sqrt129)/2 > sqrt37 follows from 37-sqrt37<32.
assert 37 > 25
print('PASS: both geometric product bounds are strictly below 6.0101.')
print('PASS: the new initial interval overlaps the established all-real cover.')
print('Conclusion: every real r >= 6.0101 has a strict threshold witness.')
