# Reference verification run

All **five manuscript checking programs passed** on 30 October 2026 under CPython 3.12.14 on macOS/arm64. The complete run took approximately half a second on the collection machine. Runtime is informational, not a performance guarantee.

[`summary.json`](summary.json) records the exact interpreter version, platform, input hashes, command arguments, exit statuses, timings, and agreement checks. Each program has separate `*.stdout.txt` and `*.stderr.txt` files; empty stderr files are retained to record the absence of diagnostics.

| Program | Status | Main output |
| --- | --- | --- |
| `check_exact.py` | PASS | [Finite graph and interval report](check_exact.stdout.txt) |
| `check_uniform_cutoff.py` | PASS | [Exact cutoff comparisons](check_uniform_cutoff.stdout.txt) |
| `check_subsix.py` | PASS | [Embedded-data certificate report](check_subsix.stdout.txt) |
| `audit_subsix_from_tex.py` | PASS | [Independent full source-derived audit](audit_subsix_from_tex.json) |
| `audit_uniform_independent.py` | PASS | [Independent graph, interval and cutoff audit](audit_uniform_independent.json) |

The runner additionally confirmed that the two sub-six checking routes agree on the SHA-256 hashes of the three complete integer matrices, the factor orders 54 and 27, and the numbers of codeword pairs checked.

To reproduce the run from the package root without changing this reference directory:

```sh
python3 run_checks.py
python3 scripts/export_data.py --check
```

Use Python 3.11 or later, with assertions enabled. Fresh run reports go to `reproduced-results/`. To deliberately replace a reference run, use `python3 run_checks.py --output-dir results`; refresh the integrity manifest afterward.

During packaging, the full suite and export comparison also passed from a detached temporary copy, launched from an unrelated working directory. The exported sub-six graphs, both fitting matrices, both codes (converted back to native vertex labels), and all three Gram matrices were compared with the independent source-derived audit and agreed. The nine computational input files (two appendix data sources and seven verification modules) are recorded by SHA-256 in the run summary. The full manuscript is not an input to these finite checks.

These are finite verification and portability results. They do not formally verify the analytical arguments in the manuscript. Decimal values in the uniform audit are explicitly labelled diagnostics.
