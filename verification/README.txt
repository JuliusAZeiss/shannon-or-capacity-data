Exact verification programs

These programs accompany the supplementary appendices of
"Hedetniemi's conjecture for Shannon OR-capacity: counterexamples and
asymptotic separation" by Julius A. Zeiss.

REQUIREMENTS AND COMMANDS

Use Python 3.11 or later. Only the standard library is required. Run with
assertions enabled: do not pass -O or set PYTHONOPTIMIZE. No numerical
solver, external data download, or network access is needed.

From the package root (the directory containing the two appendix data sources), run:

  python3 verification/check_exact.py
  python3 verification/check_uniform_cutoff.py
  python3 verification/check_subsix.py
  python3 verification/audit_subsix_from_tex.py
  python3 verification/audit_uniform_independent.py

All five commands must complete successfully with exit status zero and a
success/PASS message. A failed assertion or exception indicates a failed
check and must not be treated as successful verification.

Keep the verification/ directory next to threshold_supplement.tex and
subsix_supplement.tex. The programs find these sources relative to their
own file locations; the five commands may otherwise be run from any
working directory by specifying the appropriate script path.

FILES AND SCOPE

check_exact.py
  Reads the interval and projection tables in threshold_supplement.tex.
  Checks their rational bounds and overlaps, as well as the finite graph
  and OR-power code inputs used in Appendix G.

check_uniform_cutoff.py
  Checks the two strict power inequalities for the closed uniform cutoff
  6.0101, using integer arithmetic and the rational interpolation weight
  93/106. Its graph and interval inputs are checked by check_exact.py.

check_subsix.py
  Imports compact_bd.py and compact_h.py to check the finite sub-six
  construction: graph data, OR-square codes, fitting matrices over F_3,
  rational Gram certificates, and strict scalar comparisons. Its output
  includes the exact bound supporting the interval [5.193957,sqrt(27)).

compact_bd.py and compact_h.py
  Contain the finite vectors, permutations, and orbit-entry data needed by
  check_subsix.py, together with their exact reconstruction and checking
  routines. The programs are self-contained; no separate matrix files
  are required.

audit_subsix_from_tex.py
  Independently parses subsix_supplement.tex as text, reconstructs its
  printed finite data, and checks the certificates using rational LDL^T
  factorization with full matrix reconstruction. It imports neither of
  the compact certificate modules nor the other checking programs.

audit_uniform_independent.py
  Independently reconstructs the finite graph bounds, checks the interval
  table in threshold_supplement.tex, verifies its overlaps and integer-shift
  overlap inequality, and checks the strict closed-cutoff inequalities.
  It imports none of the other checking programs.

ARITHMETIC AND LIMITS

All mathematical pass/fail comparisons use exact integers, finite-field
arithmetic, or fractions.Fraction. audit_uniform_independent.py also prints
a Decimal logarithm-based diagnostic for the optimized open endpoint;
that diagnostic is not used to certify any inequality. Python 3.11 or later
is specified for a common environment including this diagnostic.

These programs verify finite inputs and explicit inequalities. They are
not a formal verification of the complete paper. In particular, the
general spectral interpolation, graph-capacity identities, infinite
coverage deduction, and asymptotic theorem rely on the analytical proofs
in the supplementary appendices.

OUTPUT FILES

The first three programs print their reports to standard output.
audit_subsix_from_tex.py prints a summary and writes
verification/audit_subsix_from_tex.json. audit_uniform_independent.py prints
a success message and writes verification/audit_uniform_independent.json.
The JSON reports are regenerated outputs, not inputs. The two audit
programs overwrite their own corresponding report on a successful run.
Python may also create a __pycache__/ directory when importing the data
modules; this directory is not an input and need not be distributed.
