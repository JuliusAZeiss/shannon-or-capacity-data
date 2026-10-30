# Shannon OR-capacity: accompanying code and data

Exact verification programs and finite data accompanying **Hedetniemi's conjecture for Shannon OR-capacity: counterexamples and asymptotic separation**, by Julius A. Zeiss.

Repository: [JuliusAZeiss/shannon-or-capacity-data](https://github.com/JuliusAZeiss/shannon-or-capacity-data).

This repository provides the finite graphs, codes, matrix certificates, and checking programs used by the paper, with a snapshot dated **30 October 2026**. It includes the two appendix sources that the exact verifiers read as data. The full manuscript, manuscript PDF, and arXiv submission archive are not distributed here.

## Reproduce the checks

Use **Python 3.11 or later**. Only the Python standard library is needed: there are no packages to install, numerical solvers to license, or external data to download.

From the repository root, run:

```sh
python3 run_checks.py
python3 scripts/export_data.py --check
```

All five manuscript checks must pass. The runner verifies its Python environment, executes the programs in a temporary copy, and saves reports in `reproduced-results/`. It returns a nonzero exit status if a check fails. The second command verifies that the three distributed JSON exports match fresh exact reconstructions. The distributed [reference reports](results/README.md) remain unchanged.

**Keep assertions enabled.** Do not use `python -O`, `python -OO`, or set `PYTHONOPTIMIZE`. Some original checkers use assertions for mathematical checks; the runner refuses optimized execution.

For the individual commands and their precise scope, see [`verification/README.txt`](verification/README.txt). These commands may also be run directly, but the two independent audits then write their reports into `verification/`.

## Contents

| Path | Contents |
| --- | --- |
| [`verification/`](verification/) | Five checking programs and two compact certificate-data modules supplied with the paper. |
| [`data/`](data/) | Deterministic JSON exports of finite graphs, codes, fitting matrices, Gram certificates, and threshold tables; see the [data conventions](data/README.md). |
| [`scripts/export_data.py`](scripts/export_data.py) | Regenerates the JSON data from the appendix sources and verification modules. |
| [`results/`](results/) | Captured output, independent audit reports, and a machine-readable summary from a fresh verification run. |
| [`subsix_supplement.tex`](subsix_supplement.tex) | Appendix F data source: graphs, codes, finite-field and real matrix certificates for the sub-six example. |
| [`threshold_supplement.tex`](threshold_supplement.tex) | Appendix G data source: finite graph inputs, threshold tables, and uniform cutoff. |
| [`COVERAGE.md`](COVERAGE.md) | Map from paper claims to accompanying files, including the limits of the checks. |
| [`PROVENANCE.md`](PROVENANCE.md), [`provenance.json`](provenance.json) | Source locations, inclusion decisions, and hashes of the distributed inputs. |
| [`SHA256SUMS`](SHA256SUMS) | Integrity hashes of the distributed files. |

Keep the two appendix sources immediately above `verification/`: the source-derived checks read their printed data as text. No LaTeX installation is required. The exported JSON files are convenient for reuse; the original checkers continue to verify the embedded data and appendix sources.

## Results and scope

The finite checks cover the current **54-by-27** construction, including two 27-word OR-square codes, ranks over the field with three elements, exact rational Gram certificates, and the strict product-capacity bound below `5.193957`. They also cover the finite graph inputs, eight-row interval cover, overlaps, and exact integer inequalities used for the uniform cutoff `6.0101`.

Mathematical pass/fail comparisons use exact arithmetic. The independent uniform audit also reports a Decimal approximation of the optimized open endpoint; that approximation is diagnostic and is not used to certify the cutoff.

These programs verify finite inputs and explicit inequalities. The general spectral arguments, infinite threshold-coverage deduction, and asymptotic separation theorem rely on the paper's analytical proofs. The exact capacities of the sub-six factors are not computed, and the repository does not claim an optimal threshold.

## Regenerate data and check file integrity

```sh
python3 scripts/export_data.py
```

The data export is deterministic. See [`data/README.md`](data/README.md) for its conventions and validation scope.

To check an untouched checkout on macOS:

```sh
shasum -a 256 -c SHA256SUMS
```

On systems providing `sha256sum`, use `sha256sum -c SHA256SUMS`. The manifest excludes itself, Git metadata, and the ignored `reproduced-results/` directory. Hashes verify file integrity, not mathematical correctness.

## Citation

Julius A. Zeiss, *Shannon OR-capacity: accompanying code and data*, GitHub repository, 2026. <https://github.com/JuliusAZeiss/shannon-or-capacity-data>.

The accompanying manuscript is *Hedetniemi's conjecture for Shannon OR-capacity: counterexamples and asymptotic separation*. No DOI or reuse license has been assigned to this code/data snapshot.
