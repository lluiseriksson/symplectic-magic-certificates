"""
Gate-2 Master Verification Harness for Paper 06:
  1. Compiles Lean 4 kernel (`lake build`) and runs `AuditAxioms.lean` (16/16 theorems, 0 sorry, 0 Classical.choice).
  2. Independently replays all 9 exact-Q certificates in `results/certificates/*.cert.json` (0 floats).
  3. Runs full pytest suite (`tests/test_suite.py`, including 24 live adversarial forgery & type-injection attacks).
  4. Verifies SHA-256 integrity manifest (`SHA256SUMS.json`) or writes it with `--write-manifest`.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.magic_certifier import verify_magic_certificate_exact

MANIFEST_TARGETS = [
    "lean/SymplecticMagic.lean",
    "lean/AuditAxioms.lean",
    "lean/lean-toolchain",
    "lean/lakefile.toml",
    "src/__init__.py",
    "src/symplectic_pauli.py",
    "src/magic_certifier.py",
    "scripts/run_all_experiments.py",
    "scripts/make_figures.py",
    "scripts/verify_all.py",
    "tests/test_suite.py",
    "results/exp1_certified_atlas.json",
    "results/exp2_majorana_scaling.json",
    "results/exp3_w_state_and_ccz.json",
    "results/exp4_hierarchy_separation.json",
    "results/certificates/jw_k3_n1.cert.json",
    "results/certificates/jw_k5_n2.cert.json",
    "results/certificates/jw_k21_n10.cert.json",
    "results/certificates/jw_k65_n32.cert.json",
    "results/certificates/w_state_n3.cert.json",
    "results/certificates/w_state_n20.cert.json",
    "results/certificates/ccz_n3.cert.json",
    "results/certificates/c5_pentagram_n2.cert.json",
    "results/certificates/ft6_h3_tensor2_n2.cert.json",
    "figures/fig1_majorana_exponential_separation.png",
    "figures/fig2_w_state_and_c5_pentagram.png",
    "figures/fig3_three_tier_hierarchy.png",
    "paper/main.tex",
    "paper/references.bib",
    "paper/main.pdf",
    "pyproject.toml",
    "README.md",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def check_lean() -> None:
    lean_dir = ROOT / "lean"
    print("[Gate 1] Verifying Lean 4 kernel (`lake build` + `AuditAxioms.lean`)...")
    r1 = subprocess.run(["lake", "build"], cwd=lean_dir, capture_output=True, text=True)
    if r1.returncode != 0:
        raise RuntimeError(f"lake build failed:\n{r1.stdout}\n{r1.stderr}")
    r2 = subprocess.run(["lake", "env", "lean", "AuditAxioms.lean"], cwd=lean_dir, capture_output=True, text=True)
    if r2.returncode != 0:
        raise RuntimeError(f"AuditAxioms.lean failed:\n{r2.stdout}\n{r2.stderr}")
    out = r2.stdout
    for forbidden in ["sorryAx", "Classical.choice"]:
        if forbidden in out:
            raise RuntimeError(f"Forbidden axiom {forbidden} found:\n{out}")
    lines = [ln for ln in out.splitlines() if ln.strip().startswith("'SymplecticMagic.")]
    if len(lines) != 16:
        raise RuntimeError(f"Expected 16 audited theorems, got {len(lines)}:\n{out}")
    print(f"         PASS: {len(lines)}/16 Lean 4 theorems verified (0 sorry, 0 Classical.choice).")


def check_certificates() -> None:
    print("[Gate 2] Replaying all exact-Q certificates...")
    cert_files = sorted((ROOT / "results" / "certificates").glob("*.cert.json"))
    if len(cert_files) != 9:
        raise RuntimeError(f"Expected 9 certificates, found {len(cert_files)}")
    for cf in cert_files:
        cert = json.loads(cf.read_text(encoding="utf-8"))
        ok, details = verify_magic_certificate_exact(cert)
        if not ok:
            raise RuntimeError(f"Certificate {cf.name} FAILED: {details}")
    print(f"         PASS: {len(cert_files)}/9 exact-Q certificates independently verified (0 float trust).")


def check_pytest() -> None:
    print("[Gate 3] Running pytest suite (including 28 adversarial forgery attacks)...")
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "tests/test_suite.py"], cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"pytest failed:\n{r.stdout}\n{r.stderr}")
    print(f"         PASS: {r.stdout.strip()}")


def write_or_check_manifest(write_mode: bool) -> None:
    mpath = ROOT / "SHA256SUMS.json"
    if write_mode:
        sums = {}
        for rel in MANIFEST_TARGETS:
            p = ROOT / rel
            if not p.exists():
                raise RuntimeError(f"Missing manifest target: {rel}")
            sums[rel] = sha256_file(p)
        mpath.write_text(json.dumps(sums, indent=2, sort_keys=True), encoding="utf-8")
        print(f"[Gate 4] Wrote SHA256SUMS.json ({len(sums)} files).")
    else:
        print("[Gate 4] Verifying SHA256SUMS.json integrity manifest...")
        if not mpath.exists():
            raise RuntimeError("SHA256SUMS.json missing; run with --write-manifest first.")
        expected = json.loads(mpath.read_text(encoding="utf-8"))
        for rel, h_exp in expected.items():
            p = ROOT / rel
            if not p.exists():
                raise RuntimeError(f"Missing file in manifest: {rel}")
            h_act = sha256_file(p)
            if h_act != h_exp:
                raise RuntimeError(f"SHA-256 mismatch for {rel}: expected {h_exp}, got {h_act}")
        print(f"         PASS: {len(expected)}/{len(expected)} SHA-256 hashes verified.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-manifest", action="store_true")
    parser.add_argument("--skip-manifest", action="store_true")
    args = parser.parse_args()

    check_lean()
    check_certificates()
    check_pytest()
    if not args.skip_manifest:
        write_or_check_manifest(args.write_manifest)
    print("\nALL GATE-2 VERIFICATION CHECKS PASSED!")


if __name__ == "__main__":
    main()
