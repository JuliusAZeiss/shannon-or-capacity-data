#!/usr/bin/env python3
"""Run the paper's five exact checkers using only the Python standard library.

The distributed inputs remain untouched: each checker runs in a temporary
copy, and its output is collected in reproduced-results/ (or --output-dir).
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parent
CHECKS = (
    ("check_exact.py", "All exact finite checks passed", None),
    ("check_uniform_cutoff.py", "PASS: both geometric product bounds", None),
    ("check_subsix.py", "PASS: compact exact 54-by-27 certificate", None),
    ("audit_subsix_from_tex.py", '"status": "PASS"', "audit_subsix_from_tex.json"),
    ("audit_uniform_independent.py", "PASS: independent primitive graph", "audit_uniform_independent.json"),
)


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def cross_check(module_report, source_report):
    """Bind the independently reconstructed certificates to the embedded data."""
    matrix_hashes = {
        name: module_report["matrix_SHA256"][name]
        == source_report["matrices"][name]["matrix_sha256"]
        for name in ("B", "D", "H")
    }
    orders = {
        "G": module_report["graph_orders"]["G"] == source_report["G"]["order"] == 54,
        "H": module_report["graph_orders"]["H"] == source_report["graphs"]["H"]["order"] == 27,
    }
    pairs = module_report["square_code_pairs_checked"] == [
        source_report["codes"]["G_pair_tests"],
        source_report["codes"]["H_pair_tests"],
    ]
    passed = all(matrix_hashes.values()) and all(orders.values()) and pairs
    return {"status": "PASS" if passed else "FAIL", "matrix_hash_agreement": matrix_hashes,
            "graph_order_agreement": orders, "code_pair_count_agreement": pairs}


def main():
    if sys.version_info < (3, 11):
        print("ERROR: Python 3.11 or later is required.", file=sys.stderr)
        return 2
    if not __debug__ or sys.flags.optimize or os.environ.get("PYTHONOPTIMIZE"):
        print("ERROR: assertions must be enabled; remove -O/-OO and unset PYTHONOPTIMIZE.",
              file=sys.stderr)
        return 2

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "reproduced-results",
                        help="directory for reports and captured output (default: reproduced-results/ beside this script)")
    parser.add_argument("--timeout", type=float, default=600,
                        help="maximum seconds per checker (default: 600)")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")

    inputs = [ROOT / name for name in ("subsix_supplement.tex", "threshold_supplement.tex")]
    inputs += [ROOT / "verification" / name for name, _, _ in CHECKS]
    inputs += [ROOT / "verification" / name for name in ("compact_bd.py", "compact_h.py")]
    missing = [p.relative_to(ROOT).as_posix() for p in inputs if not p.is_file()]
    if missing:
        print("ERROR: missing input files: " + ", ".join(missing), file=sys.stderr)
        return 2

    output_dir = args.output_dir.resolve()
    if any(p.resolve().is_relative_to(output_dir) for p in inputs):
        parser.error("--output-dir must not contain the distributed source inputs")
    output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "schema_version": 1,
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": {"system": platform.system(), "release": platform.release(),
                     "machine": platform.machine()},
        "assertions_enabled": True,
        "input_sha256": {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in sorted(inputs)},
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "checks": [],
        "generated_audit_reports": [],
        "scope": "Exact finite graph, code, matrix-certificate and scalar checks. The paper's analytical arguments are not formally verified by these programs.",
        "arithmetic": "Proof comparisons use integer, finite-field or rational arithmetic. Decimal logarithm diagnostics in the independent uniform audit do not certify inequalities.",
        "log_path_redaction": "Temporary package paths in captured stdout/stderr are replaced by <verification-workdir>. No mathematical values are changed.",
    }
    started = time.monotonic()
    env = os.environ.copy()
    env.pop("PYTHONOPTIMIZE", None)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"
    captured = {}
    audits = {}

    with tempfile.TemporaryDirectory(prefix="shannon-or-checks-") as temporary:
        work = Path(temporary).resolve()
        for source in inputs:
            destination = work / source.relative_to(ROOT)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
        for name, success_marker, audit_name in CHECKS:
            command = [sys.executable, "-B", "verification/" + name]
            if name == "audit_subsix_from_tex.py":
                command += ["subsix_supplement.tex", "verification/" + audit_name]
            begin = time.monotonic()
            timed_out = False
            try:
                process = subprocess.run(command, cwd=work, env=env, capture_output=True,
                                         text=True, encoding="utf-8", timeout=args.timeout)
                stdout, stderr, returncode = process.stdout, process.stderr, process.returncode
            except subprocess.TimeoutExpired as exc:
                def as_text(value):
                    return value.decode("utf-8", errors="replace") if isinstance(value, bytes) else (value or "")
                stdout, stderr, returncode = as_text(exc.stdout), as_text(exc.stderr), None
                stderr += "\nERROR: checker exceeded the configured timeout.\n"
                timed_out = True
            # Keep recorded output portable, including exception tracebacks.
            for temporary_path in (str(work), str(temporary)):
                stdout = stdout.replace(temporary_path, "<verification-workdir>")
                stderr = stderr.replace(temporary_path, "<verification-workdir>")
            stem = Path(name).stem
            stdout_name, stderr_name = stem + ".stdout.txt", stem + ".stderr.txt"
            (output_dir / stdout_name).write_text(stdout, encoding="utf-8")
            (output_dir / stderr_name).write_text(stderr, encoding="utf-8")
            passed = returncode == 0 and success_marker in stdout
            item = {"script": "verification/" + name,
                    "command": ["python3"] + command[1:],
                    "returncode": returncode, "timed_out": timed_out,
                    "success_marker_found": success_marker in stdout,
                    "status": "PASS" if passed else "FAIL",
                    "elapsed_seconds": round(time.monotonic() - begin, 6),
                    "stdout": stdout_name, "stderr": stderr_name}
            if passed and audit_name:
                artifact = work / "verification" / audit_name
                if artifact.is_file():
                    audits[name] = json.loads(artifact.read_text(encoding="utf-8"))
                    shutil.copyfile(artifact, output_dir / audit_name)
                    report["generated_audit_reports"].append(audit_name)
                else:
                    item["status"] = "FAIL"
                    item["error"] = "Expected audit JSON was not produced."
            captured[name] = stdout
            report["checks"].append(item)
            print(item["status"] + ": " + name + " (" + str(item["elapsed_seconds"]) + " s)", flush=True)

    if all(item["status"] == "PASS" for item in report["checks"]):
        try:
            report["independent_agreement"] = cross_check(
                json.loads(captured["check_subsix.py"]), audits["audit_subsix_from_tex.py"])
        except (KeyError, ValueError, TypeError) as exc:
            report["independent_agreement"] = {"status": "FAIL", "error": str(exc)}
    else:
        report["independent_agreement"] = {"status": "NOT_RUN", "reason": "At least one checker failed."}
    passed = (all(item["status"] == "PASS" for item in report["checks"])
              and report["independent_agreement"]["status"] == "PASS")
    report["status"] = "PASS" if passed else "FAIL"
    report["elapsed_seconds"] = round(time.monotonic() - started, 6)
    report["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    write_json(output_dir / "summary.json", report)
    print(report["status"] + ": five finite checkers and independent certificate agreement. "
          + "Report: " + str(output_dir / "summary.json"))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
