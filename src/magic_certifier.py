"""
Zero-Floating-Point-Trust Exact-Q Certificate Verifier for Paper 06 (`SYMPLECTIC-MAGIC-CERT-1.0`).

Guarantees:
  - Recursive rejection of ANY `float` or `bool` anywhere in the input JSON tree.
  - Strict rational parsing (`Fraction` from exact `int` or canonical `"p"` / `"p/q"` strings with q > 0).
  - Independent symplectic verification of every Anticommuting Clique Atom (`omega(u, v) == 1` for all pairs)
    and every Signed Frustrated Commuting-Triple Atom (`omega == 0`, `A ^ B == C`, `s_A * s_B * s_C * sigma == -1`).
  - Independent exact-Q summation of the conical decomposition `W = sum_j y_j C_j + sum_l z_l F_l`
    matching the declared traceless witness `W` identically.
  - Exact Chvatal-Gomory integer-floor tightening ALLOWED ONLY when every coefficient of `W` is an exact integer.
  - Independent exact-Q verification of state expectation bounds (exact rational or algebraically certified
    radical intervals `[lo, hi]` with `lo^2 * rad_den <= rad_num` verified strictly in integer arithmetic!).
  - Independent exact-Q replay of `R_lower_bound`, `p_noise_cert`, `worst_case_rac_stab_ceiling`, and
    optional full brute-force `STAB_n` verification for `n <= 3`.
"""
from __future__ import annotations

from fractions import Fraction
import math
import re

from src.symplectic_pauli import (
    Pauli,
    symp_inner,
    pauli_mul_commuting_sign,
    str_to_pauli,
    enumerate_stabilizer_states_exact,
)

_FRAC_RE = re.compile(r"^-?[0-9]+(?:/[1-9][0-9]*)?$")


def _is_strict_int(v) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


def _assert_no_float_or_bool(obj, path: str = "root") -> None:
    """Strictly whitelist JSON-exact primitives (dict, list, str, strict int, None); reject bool, float, Decimal, sets, tuples, and non-str keys."""
    if isinstance(obj, bool):
        raise ValueError(f"{path}: bool value {obj!r} forbidden in zero-float exact certificate")
    if isinstance(obj, float):
        raise ValueError(f"{path}: float value {obj!r} forbidden in zero-float exact certificate")
    if isinstance(obj, dict):
        for k, v in obj.items():
            if not isinstance(k, str):
                raise ValueError(f"{path}: non-string dict key {k!r} forbidden")
            _assert_no_float_or_bool(v, f"{path}.{k}")
    elif isinstance(obj, list):
        for idx, item in enumerate(obj):
            _assert_no_float_or_bool(item, f"{path}[{idx}]")
    elif obj is None or isinstance(obj, (str, int)):
        return
    else:
        raise ValueError(f"{path}: unsupported non-exact type {type(obj).__name__} ({obj!r}) forbidden")


def _parse_exact_fraction(raw, name: str) -> Fraction:
    if _is_strict_int(raw):
        return Fraction(raw, 1)
    if isinstance(raw, str) and _FRAC_RE.fullmatch(raw):
        if "/" in raw:
            num_s, den_s = raw.split("/")
            f = Fraction(int(num_s), int(den_s))
        else:
            f = Fraction(int(raw), 1)
        canonical = f"{f.numerator}" if f.denominator == 1 else f"{f.numerator}/{f.denominator}"
        if raw != canonical:
            raise ValueError(f"{name}: non-canonical rational string {raw!r} (expected {canonical!r})")
        return f
    raise ValueError(f"{name}: expected strict int or canonical 'p/q' rational string, got {raw!r}")


def _parse_pauli_str(raw, n_qubits: int, name: str) -> Pauli:
    if not isinstance(raw, str) or len(raw) != n_qubits:
        raise ValueError(f"{name}: expected Pauli string of length {n_qubits}, got {raw!r}")
    u = str_to_pauli(raw)
    if u == (0, 0):
        raise ValueError(f"{name}: identity Pauli 'I...I' not allowed in traceless witness")
    return u


def _parse_sign(raw, name: str) -> int:
    if not _is_strict_int(raw) or raw not in (-1, 1):
        raise ValueError(f"{name}: sign must be strict int -1 or +1, got {raw!r}")
    return raw


def _verify_radical_interval(entry: dict, name: str) -> tuple[Fraction, Fraction]:
    """
    Verify an expectation entry for a Pauli P:
      - kind == "exact": exact rational value `val` in [-1, 1] -> returns (val, val).
      - kind == "radical_lower_bound": verifies exact rational interval [lo, hi] in [-1, 1]
        where the true expectation is `sign * sqrt(rad_num / rad_den)` (with rad_num >= 0, rad_den > 0):
          if sign == +1: verifies 0 <= lo <= hi <= 1, lo^2 * rad_den <= rad_num, hi^2 * rad_den >= rad_num.
          if sign == -1: verifies -1 <= lo <= hi <= 0, hi^2 * rad_den <= rad_num, lo^2 * rad_den >= rad_num.
    """
    if not isinstance(entry, dict):
        raise ValueError(f"{name}: expectation entry must be dict")
    kind = entry.get("kind")
    if kind == "exact":
        v = _parse_exact_fraction(entry["val"], f"{name}.val")
        if v < -1 or v > 1:
            raise ValueError(f"{name}: expectation {v} outside [-1, 1]")
        return v, v
    elif kind == "radical_interval":
        lo = _parse_exact_fraction(entry["lo"], f"{name}.lo")
        hi = _parse_exact_fraction(entry["hi"], f"{name}.hi")
        sgn = _parse_sign(entry["sign"], f"{name}.sign")
        rad_num = _parse_exact_fraction(entry["rad_num"], f"{name}.rad_num")
        rad_den = _parse_exact_fraction(entry["rad_den"], f"{name}.rad_den")
        if rad_num < 0 or rad_den <= 0 or lo > hi or lo < -1 or hi > 1:
            raise ValueError(f"{name}: invalid radical_interval parameters")
        target_sq = rad_num / rad_den
        if target_sq > 1:
            raise ValueError(f"{name}: radical squared value {target_sq} exceeds 1")
        if sgn == 1:
            if lo < 0 or lo * lo > target_sq or hi * hi < target_sq:
                raise ValueError(f"{name}: positive radical interval [{lo}, {hi}] fails exact square enclosure")
        else:
            if hi > 0 or hi * hi > target_sq or lo * lo < target_sq:
                raise ValueError(f"{name}: negative radical interval [{lo}, {hi}] fails exact square enclosure")
        return lo, hi
    else:
        raise ValueError(f"{name}: unknown expectation kind {kind!r}")


def verify_magic_certificate_exact(cert: dict) -> tuple[bool, dict]:
    """
    100% ZERO-FLOATING-POINT-TRUST exact-Q verifier for `SYMPLECTIC-MAGIC-CERT-1.0`.
    """
    if not isinstance(cert, dict):
        return False, {"reason": "certificate must be dict"}
    try:
        _assert_no_float_or_bool(cert)
        if cert.get("schema_version") != "SYMPLECTIC-MAGIC-CERT-1.0":
            return False, {"reason": f"invalid schema_version {cert.get('schema_version')!r}"}

        n_qubits = cert["n_qubits"]
        if not _is_strict_int(n_qubits) or n_qubits < 1:
            return False, {"reason": "invalid n_qubits"}

        # 1. Parse declared traceless witness W = sum_P w_P * P
        raw_W = cert["witness_W"]
        if not isinstance(raw_W, dict) or len(raw_W) == 0:
            return False, {"reason": "witness_W must be non-empty dict"}
        witness_W: dict[Pauli, Fraction] = {}
        for p_str, raw_c in raw_W.items():
            u = _parse_pauli_str(p_str, n_qubits, f"witness_W[{p_str!r}]")
            c = _parse_exact_fraction(raw_c, f"witness_W[{p_str!r}]")
            if c == 0:
                return False, {"reason": f"zero weight entry for {p_str!r} forbidden"}
            if u in witness_W:
                return False, {"reason": f"duplicate Pauli {p_str!r} in witness_W"}
            witness_W[u] = c

        # 2. Verify all dual atoms (cliques & frustrated triples) and accumulate exact W_reconstructed
        atoms = cert["atoms"]
        if not isinstance(atoms, list) or len(atoms) == 0:
            return False, {"reason": "atoms must be non-empty list"}

        W_recon: dict[Pauli, Fraction] = {}
        y_clique_sum = Fraction(0)
        z_ftriple_sum = Fraction(0)
        n_cliques = 0
        n_ftriples = 0

        for idx, atom in enumerate(atoms):
            if not isinstance(atom, dict):
                return False, {"reason": f"atom[{idx}] must be dict"}
            akind = atom.get("kind")
            weight = _parse_exact_fraction(atom["weight"], f"atom[{idx}].weight")
            if weight <= 0:
                return False, {"reason": f"atom[{idx}]: weight must be strictly positive (got {weight})"}
            raw_terms = atom["terms"]
            if not isinstance(raw_terms, list) or len(raw_terms) == 0:
                return False, {"reason": f"atom[{idx}].terms must be non-empty list"}

            parsed_terms: list[tuple[Pauli, int]] = []
            seen_u: set[Pauli] = set()
            for t_i, term in enumerate(raw_terms):
                if not isinstance(term, list) or len(term) != 2:
                    return False, {"reason": f"atom[{idx}].terms[{t_i}] must be [pauli_str, sign]"}
                u = _parse_pauli_str(term[0], n_qubits, f"atom[{idx}].terms[{t_i}][0]")
                s = _parse_sign(term[1], f"atom[{idx}].terms[{t_i}][1]")
                if u in seen_u:
                    return False, {"reason": f"atom[{idx}]: duplicate Pauli {term[0]!r} in atom"}
                seen_u.add(u)
                parsed_terms.append((u, s))

            if akind == "clique":
                # Verify EVERY pair anticommutes: omega(u, v) == 1
                for i in range(len(parsed_terms)):
                    for j in range(i + 1, len(parsed_terms)):
                        if symp_inner(parsed_terms[i][0], parsed_terms[j][0]) != 1:
                            return False, {"reason": f"atom[{idx}]: commuting pair found inside claimed anticommuting clique"}
                y_clique_sum += weight
                n_cliques += 1
            elif akind == "frustrated_triple":
                if len(parsed_terms) != 3:
                    return False, {"reason": f"atom[{idx}]: frustrated_triple must have exactly 3 terms"}
                (uA, sA), (uB, sB), (uC, sC) = parsed_terms
                if symp_inner(uA, uB) != 0 or symp_inner(uA, uC) != 0 or symp_inner(uB, uC) != 0:
                    return False, {"reason": f"atom[{idx}]: frustrated_triple terms do not pairwise commute"}
                sig, prod_AB = pauli_mul_commuting_sign(uA, uB)
                if prod_AB != uC:
                    return False, {"reason": f"atom[{idx}]: frustrated_triple fails group closure A*B != +/-C"}
                if sA * sB * sC * sig != -1:
                    return False, {"reason": f"atom[{idx}]: triple is NOT frustrated (sA*sB*sC*sigma != -1)"}
                z_ftriple_sum += weight
                n_ftriples += 1
            else:
                return False, {"reason": f"atom[{idx}]: unknown atom kind {akind!r}"}

            for u, s in parsed_terms:
                W_recon[u] = W_recon.get(u, Fraction(0)) + weight * s

        # Strip exact zeros in W_recon and verify exact identity W_recon == witness_W
        W_recon = {u: c for u, c in W_recon.items() if c != 0}
        if W_recon != witness_W:
            return False, {"reason": "conical atom decomposition does not exactly reconstruct witness_W"}

        # 3. Compute exact stabilizer ceiling beta_plus and floor -beta_minus
        beta_plus_raw = y_clique_sum + z_ftriple_sum
        beta_minus_raw = y_clique_sum + 3 * z_ftriple_sum

        cg_flag = cert.get("use_chvatal_gomory_floor", "no")
        if cg_flag not in ("yes", "no"):
            return False, {"reason": f"use_chvatal_gomory_floor must be 'yes' or 'no', got {cg_flag!r}"}
        if cg_flag == "yes":
            if not all(c.denominator == 1 for c in witness_W.values()):
                return False, {"reason": "Chvatal-Gomory integer floor forbidden for non-integer witness_W"}
            beta_plus = Fraction(math.floor(beta_plus_raw), 1)
            beta_minus = Fraction(math.floor(beta_minus_raw), 1)
        else:
            beta_plus = beta_plus_raw
            beta_minus = beta_minus_raw

        if beta_plus <= 0 or (beta_plus + beta_minus) <= 0:
            return False, {"reason": "non-positive stabilizer bounds"}

        # 4. Verify state expectation bounds and compute exact-Q W_rho_lb, R_lower_bound, p_noise_cert
        raw_exps = cert["state_expectations"]
        if not isinstance(raw_exps, dict):
            return False, {"reason": "state_expectations must be dict"}
        parsed_exps: dict[Pauli, tuple[Fraction, Fraction]] = {}
        for p_str, entry in raw_exps.items():
            u = _parse_pauli_str(p_str, n_qubits, f"state_expectations[{p_str!r}]")
            if u in parsed_exps:
                return False, {"reason": f"duplicate Pauli {p_str!r} in state_expectations"}
            parsed_exps[u] = _verify_radical_interval(entry, f"state_expectations[{p_str!r}]")

        W_rho_lb = Fraction(0)
        for u, coeff in witness_W.items():
            if u not in parsed_exps:
                return False, {"reason": "missing state expectation for witness Pauli"}
            lo, hi = parsed_exps[u]
            W_rho_lb += coeff * lo if coeff > 0 else coeff * hi

        if W_rho_lb <= beta_plus:
            return False, {"reason": f"state fails magic certification: W_rho_lb={W_rho_lb} <= beta_plus={beta_plus}"}

        R_lb = Fraction(1) + 2 * (W_rho_lb - beta_plus) / (beta_plus + beta_minus)
        p_noise = Fraction(1) - beta_plus / W_rho_lb

        # 5. Match all recorded exact rational claims
        checks = [
            ("claimed_beta_plus", beta_plus),
            ("claimed_beta_minus", beta_minus),
            ("claimed_W_rho_lb", W_rho_lb),
            ("claimed_R_lower_bound", R_lb),
            ("claimed_p_noise_cert", p_noise),
        ]
        for key, actual in checks:
            if key not in cert:
                return False, {"reason": f"missing required claim {key}"}
            claimed = _parse_exact_fraction(cert[key], key)
            if claimed != actual:
                return False, {"reason": f"{key} mismatch: claimed={claimed} vs actual={actual}"}

        # 6. Optional exact brute-force STAB_n check for n <= 3
        stab_checked = 0
        bf_flag = cert.get("verify_exact_stab_enumeration", "no")
        if bf_flag not in ("yes", "no"):
            return False, {"reason": f"verify_exact_stab_enumeration must be 'yes' or 'no', got {bf_flag!r}"}
        if bf_flag == "yes":
            if n_qubits > 3:
                return False, {"reason": "exact STAB_n enumeration only supported for n_qubits <= 3"}
            states = enumerate_stabilizer_states_exact(n_qubits)
            for st in states:
                val = sum(c * st.get(u, 0) for u, c in witness_W.items())
                if val > beta_plus or val < -beta_minus:
                    return False, {"reason": f"STAB_{n_qubits} counterexample found: val={val} outside [{-beta_minus}, {beta_plus}]"}
            stab_checked = len(states)

        return True, {
            "n_qubits": n_qubits,
            "n_paulis": len(witness_W),
            "n_cliques": n_cliques,
            "n_ftriples": n_ftriples,
            "beta_plus": str(beta_plus),
            "beta_minus": str(beta_minus),
            "W_rho_lb": str(W_rho_lb),
            "R_lower_bound": str(R_lb),
            "p_noise_cert": str(p_noise),
            "stab_checked": stab_checked,
        }
    except Exception as e:
        return False, {"reason": f"verification exception: {e}"}
