# Provenance

This code/data repository was prepared on **30 October 2026** from the local accompanying package for *Hedetniemi's conjecture for Shannon OR-capacity: counterexamples and asymptotic separation* by Julius A. Zeiss.

The manuscript inputs originate from the collection-relative directory:

```text
Paper and Submission/shannon_or_hedetniemi_jctb/
```

## Computational inputs

[`provenance.json`](provenance.json) records original collection-relative paths, destination paths, byte counts, and SHA-256 hashes for the distributed source inputs:

- `subsix_supplement.tex` and `threshold_supplement.tex`, retained as exact data sources for the independent auditors;
- all seven Python programs/modules in `verification/`, unchanged from the manuscript package;
- `verification/README.txt`, whose package-root description was adapted for this code/data-only layout.

The full manuscript (`main.tex`), manuscript PDFs, and arXiv submission folder/archives are excluded. None is required to run the finite checks. The appendix sources contain the printed certificate data; they do not form a standalone compiled paper in this repository.

## Repository helpers and generated data

The README, this provenance guide, coverage map, `.gitignore`, verification runner, data exporter, data guide, and integrity manifest support reproduction and reuse. The runner was adapted from the local accompanying package to copy and hash only the two appendix sources; the full manuscript is not a computational input.

Exported JSON data is derived from the appendix sources and the compact verification modules. The Section 2 example is generated from its explicit graph definition as implemented by the exporter. These exports do not replace the paper's proofs or independently checked certificate inputs.

`results/` contains a fresh verification run and its machine-readable summary. The summary records Python/platform information, input hashes, command outcomes, and runtimes. The runner executes the unchanged checking programs in a temporary copy so that generated reports do not alter the distributed inputs.

The original programs contain assertions and must run without optimization. The runner and exporter refuse optimized execution. Exact mathematical checks and approximate endpoint diagnostics retain their distinct roles in the documentation.

## Integrity

`SHA256SUMS` hashes distributed files except the manifest itself. It excludes Git metadata and local reproduced results. It is an integrity record, not a proof certificate. Default reruns write to the ignored `reproduced-results/` directory.

Repository: <https://github.com/JuliusAZeiss/shannon-or-capacity-data>. No DOI or license grant is added by this provenance record.
