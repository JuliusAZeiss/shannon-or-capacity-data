# Coverage of the current paper

This inventory refers to the accompanying manuscript dated 30 October 2026. The repository contains its two appendix data sources, not the full manuscript. It covers every explicitly referenced verification program and its runtime data dependencies. The appendix labels below remain useful if theorem numbering changes in a later manuscript.

| Paper location | Finite data or claim | Accompanying files |
| --- | --- | --- |
| Data and code availability statement | Verification programs distributed with the source | All seven `verification/*.py` files and `verification/README.txt` |
| Appendix F.2 (`subsix:sec:graphs`) | The 27 ternary vectors, B and its complement D, G = B disjoint union D, the 27-vertex H adjacency table, and both 27-word OR-square codes | `compact_bd.py`, `compact_h.py`, `check_subsix.py`, `audit_subsix_from_tex.py`, `subsix_supplement.tex` |
| Appendix F.3 (`subsix:sec:rank`) | Rank-five scalar fitting matrix and rank-eleven two-block fitting matrix over F3; H's twofold colouring | `check_subsix.py`, `audit_subsix_from_tex.py` |
| Appendices F.4 and F.7 (`subsix:sec:gram`, `subsix:app:matrices`) | Integer Gram matrices, symmetry permutations and orbit entries, positivity and rational bounds | `compact_bd.py`, `compact_h.py`, `check_subsix.py`, `audit_subsix_from_tex.py` |
| Appendix F.6, Equations F.5–F.6 (`subsix:eq:power1`, `subsix:eq:power2`), Proposition F.1 (`subsix:thm:main`) | Exact power comparisons giving the strict bound below 5.193957 and the local separating interval | `check_subsix.py`, `compact_h.py`, `audit_subsix_from_tex.py` |
| Appendix F.8 (`subsix:app:verification`) | Reproduction instructions and distinction between finite verification and analytical arguments | `verification/README.txt`, the two sub-six checkers |
| Appendix E, Lemmas E.1 and E.3 (`lem:binary`, `lem:clebsch`), and Appendix G | U6, circulant R on 37 vertices, C5, complement of C9, Clebsch graph, and the finite codes used by the threshold construction | `check_exact.py`, `audit_uniform_independent.py`; graph definitions also occur in the manuscript |
| Appendix G, Table G.1 (`tab:threshold-cover`) | Eight rational interpolation rows, consecutive overlaps and integer-shift overlap | `check_exact.py`, `audit_uniform_independent.py`, `threshold_supplement.tex` |
| Appendix G, Equation G.4 (`eq:uniform-integer-check`) | The two exact integer inequalities for the closed cutoff 6.0101 | `check_uniform_cutoff.py`, `audit_uniform_independent.py` |
| Appendix G, Equation G.3 (`eq:uniform-optimized`) | Decimal approximation of the analytically optimized open endpoint | Diagnostic output of `audit_uniform_independent.py`; it is not a pass/fail proof of the endpoint |
| Section 2 (`sec:simple`), Appendix E (`prop:nine`) | Analytical U8 / F9 example and its 128-word code | Definitions and proof in the manuscript; the new exporter checks the finite 128-word code and provides `data/simple-example.json` |

Paths in the last column are relative to `verification/` unless a TeX file or other directory is specified. The [data guide](data/README.md) describes the additional machine-readable exports and their indexing conventions.

## What is verified

The five original checking commands perform exact finite checks. In the sub-six construction, one route checks embedded compact data and a second independently reconstructs the printed TeX data and checks full matrices with rational LDL factorization. The uniform audit imports none of the other checking programs. Fresh reports are retained in `results/`.

The original five programs instantiate U6, not U8. The new F9 data export and code check are repository helpers; they do not replace the analytical proof of the F9 capacity bounds. Theorem 1's extension to every integer minimum at least seven also remains analytical.

The graph H in this snapshot has **27 vertices and 176 edges**. Historical 54-by-28 data in the broader research collection is not the graph used by the current manuscript.

## What remains analytical

The general graph-capacity identities, spectral duality and interpolation, symmetry-decomposition arguments, extension of the finite cover to infinitely many thresholds, complete-factor families, and asymptotic theorem are proved in the manuscript. Running the programs does not formally verify these proofs. Finite graph exports do not enumerate infinite families.

## Excluded material

- The full manuscript (`main.tex`), manuscript PDFs, and arXiv submission archives.

- Historical manuscripts, superseded 28-vertex certificates, and exploratory search outputs unrelated to the current printed data.
- Chat archives, private working notes, submission forms, correspondence, and manuscript editing/rendering history.
- Solver installations, discovery runs, temporary files, Python caches, and LaTeX build intermediates.

The current paper names no additional discovery script or solver dependency. All its explicitly referenced verification files were located. This package supplies the finite certificate verification promised by the paper; it is not an archive of every exploratory computation that preceded the results.
