"""
Paper 06 Comprehensive Pytest & Adversarial Forgery Test Suite:
  - Tests exact symplectic Pauli algebra, Jordan-Wigner K_{2n+1} anticommutation, exact STAB_n
    counts (6, 60, 1080), and exact frustrated-triple closure.
  - Verifies all 9 published zero-float certificates in `results/certificates/*.cert.json`.
  - Attacks `verify_magic_certificate_exact` with 24 live adversarial forgery and type-injection
    attacks (commuting pair in clique, non-frustrated triple, non-closed triple, anticommuting
    triple, negative atom weight, illegal Chvatal-Gomory floor on fractional witness, tampered
    witness sum, tampered radical square enclosure, tampered R/p claims, duplicate Paulis,
    identity Pauli, schema forgery, recursive float/bool/NaN injection) and asserts 100% rejection!
"""
from __future__ import annotations

import copy
from decimal import Decimal
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pytest

from src.symplectic_pauli import (
    symp_inner,
    pauli_mul_commuting_sign,
    str_to_pauli,
    jordan_wigner_clique,
    enumerate_stabilizer_states_exact,
)
from src.magic_certifier import verify_magic_certificate_exact

CERT_DIR = ROOT / "results" / "certificates"
RES_DIR = ROOT / "results"


def _load_cert(cid: str) -> dict:
    return json.loads((CERT_DIR / f"{cid}.cert.json").read_text(encoding="utf-8"))


# 1. Exact symplectic algebra, Jordan-Wigner cliques, and STAB_n counts
def test_symplectic_algebra_and_stab_counts():
    for n in [1, 2, 5, 16]:
        clq = jordan_wigner_clique(n)
        assert len(clq) == 2 * n + 1
        for i in range(len(clq)):
            assert symp_inner(clq[i], clq[i]) == 0
            for j in range(i + 1, len(clq)):
                assert symp_inner(clq[i], clq[j]) == 1
    # Check (XX)(YY) = -ZZ
    sig, prod = pauli_mul_commuting_sign(str_to_pauli("XX"), str_to_pauli("YY"))
    assert sig == -1 and prod == str_to_pauli("ZZ")
    # Exact stabilizer cardinalities: |STAB_1|=6, |STAB_2|=60, |STAB_3|=1080
    assert len(enumerate_stabilizer_states_exact(1)) == 6
    assert len(enumerate_stabilizer_states_exact(2)) == 60
    assert len(enumerate_stabilizer_states_exact(3)) == 1080


# 2. Verify all 9 published exact-Q certificates
@pytest.mark.parametrize(
    "cert_id",
    [
        "jw_k3_n1", "jw_k5_n2", "jw_k21_n10", "jw_k65_n32",
        "w_state_n3", "w_state_n20", "ccz_n3",
        "c5_pentagram_n2", "ft6_h3_tensor2_n2",
    ],
)
def test_all_published_magic_certificates_pass(cert_id: str):
    cert = _load_cert(cert_id)
    ok, details = verify_magic_certificate_exact(cert)
    assert ok, f"Published certificate {cert_id} failed: {details}"


# 3. Adversarial Attack Suite A: Algebraic & Symplectic Forgeries (14 attacks)
def test_adversarial_algebraic_forgeries_rejected():
    jw2 = _load_cert("jw_k5_n2")
    c5 = _load_cert("c5_pentagram_n2")
    ft6 = _load_cert("ft6_h3_tensor2_n2")
    w3 = _load_cert("w_state_n3")

    # Attack 1: Sneak a commuting pair into a claimed anticommuting clique!
    a1 = copy.deepcopy(w3)
    a1["atoms"][0]["terms"][2][0] = "ZII"  # duplicate/commuting
    assert not verify_magic_certificate_exact(a1)[0]
    # Distinct commuting pair: ZII and IZI commute!
    a1b = copy.deepcopy(w3)
    a1b["witness_W"] = {"ZII": 1, "IZI": 1, "XXI": 1}
    a1b["atoms"][0]["terms"] = [["ZII", 1], ["IZI", 1], ["XXI", 1]]
    a1b["state_expectations"]["IZI"] = {"kind": "exact", "val": "1/3"}
    ok1b, d1b = verify_magic_certificate_exact(a1b)
    assert not ok1b and "commuting pair" in d1b["reason"]

    # Attack 2: Non-frustrated commuting triple (flip ZZ sign to -1 so sA*sB*sC*sigma = +1!)
    a2 = copy.deepcopy(ft6)
    a2["witness_W"]["ZZ"] = -1
    a2["atoms"][0]["terms"][2][1] = -1
    ok2, d2 = verify_magic_certificate_exact(a2)
    assert not ok2 and "NOT frustrated" in d2["reason"]

    # Attack 3: Pairwise commuting triple that fails group closure A*B != +/-C (ZII, IZI, IIZ on 3 qubits!)
    a3 = copy.deepcopy(w3)
    a3["witness_W"] = {"ZII": 1, "IZI": 1, "IIZ": 1}
    a3["atoms"] = [{"kind": "frustrated_triple", "weight": 1, "terms": [["ZII", 1], ["IZI", 1], ["IIZ", 1]]}]
    ok3, d3 = verify_magic_certificate_exact(a3)
    assert not ok3 and "group closure" in d3["reason"]

    # Attack 4: Anticommuting triple disguised as frustrated_triple
    a4 = copy.deepcopy(ft6)
    a4["atoms"][0]["terms"] = [["XI", 1], ["YI", 1], ["ZI", 1]]
    ok4, d4 = verify_magic_certificate_exact(a4)
    assert not ok4 and "pairwise commute" in d4["reason"]

    # Attack 5: Negative atom weight y_j < 0
    a5 = copy.deepcopy(c5)
    a5["atoms"][0]["weight"] = "-1/2"
    ok5, d5 = verify_magic_certificate_exact(a5)
    assert not ok5 and "strictly positive" in d5["reason"]

    # Attack 6: Illegal Chvatal-Gomory integer floor on fractional-coefficient witness_W!
    a6 = copy.deepcopy(c5)
    a6["witness_W"] = {k: "1/2" for k in a6["witness_W"]}
    for at in a6["atoms"]:
        at["weight"] = "1/4"
    a6["claimed_beta_plus"] = "1"
    a6["claimed_beta_minus"] = "1"
    ok6, d6 = verify_magic_certificate_exact(a6)
    assert not ok6 and "non-integer" in d6["reason"]

    # Attack 7: Tampered witness_W coefficient (atom sum no longer equals witness_W)
    a7 = copy.deepcopy(ft6)
    a7["witness_W"]["XX"] = 2
    ok7, d7 = verify_magic_certificate_exact(a7)
    assert not ok7 and "reconstruct" in d7["reason"]

    # Attack 8: Inflated radical lower bound (lo^2 > 1/5!)
    a8 = copy.deepcopy(jw2)
    first_p = next(iter(a8["state_expectations"]))
    a8["state_expectations"][first_p]["lo"] = "1/2"  # (1/2)^2 = 1/4 > 1/5!
    a8["state_expectations"][first_p]["hi"] = "3/5"
    ok8, d8 = verify_magic_certificate_exact(a8)
    assert not ok8 and "square enclosure" in d8["reason"]

    # Attack 9: Expectation outside [-1, 1]
    a9 = copy.deepcopy(w3)
    a9["state_expectations"]["ZII"]["val"] = "3/2"
    assert not verify_magic_certificate_exact(a9)[0]

    # Attack 10: Identity Pauli 'III' snuck into traceless witness_W
    a10 = copy.deepcopy(w3)
    a10["witness_W"]["III"] = 1
    assert not verify_magic_certificate_exact(a10)[0]

    # Attack 11: Inflated claimed R_lower_bound
    a11 = copy.deepcopy(w3)
    a11["claimed_R_lower_bound"] = "10"
    ok11, d11 = verify_magic_certificate_exact(a11)
    assert not ok11 and "claimed_R_lower_bound" in d11["reason"]

    # Attack 12: Inflated claimed p_noise_cert
    a12 = copy.deepcopy(w3)
    a12["claimed_p_noise_cert"] = "9/10"
    ok12, d12 = verify_magic_certificate_exact(a12)
    assert not ok12 and "claimed_p_noise_cert" in d12["reason"]

    # Attack 13: Non-magic state (W_rho_lb <= beta_plus)
    a13 = copy.deepcopy(w3)
    a13["state_expectations"]["ZII"]["val"] = "1/3"
    a13["state_expectations"]["XXI"]["val"] = "1/3"
    a13["state_expectations"]["YIY"]["val"] = "1/3"
    a13["claimed_W_rho_lb"] = "1"
    a13["claimed_R_lower_bound"] = "1"
    a13["claimed_p_noise_cert"] = "0"
    ok13, d13 = verify_magic_certificate_exact(a13)
    assert not ok13 and "fails magic certification" in d13["reason"]

    # Attack 14: Forged schema_version
    a14 = copy.deepcopy(w3)
    a14["schema_version"] = "FAKE-99.0"
    ok14, d14 = verify_magic_certificate_exact(a14)
    assert not ok14 and "schema_version" in d14["reason"]


# 4. Adversarial Attack Suite B: Type-Safety, Float/Bool/String Injection (10 attacks)
def test_adversarial_type_safety_and_injection_rejected():
    w3 = _load_cert("w_state_n3")
    c5 = _load_cert("c5_pentagram_n2")

    # Attack 15: Float in witness weight
    t1 = copy.deepcopy(w3); t1["witness_W"]["ZII"] = 1.0
    assert not verify_magic_certificate_exact(t1)[0]

    # Attack 16: Bool True in witness weight (True == 1!)
    t2 = copy.deepcopy(w3); t2["witness_W"]["ZII"] = True
    assert not verify_magic_certificate_exact(t2)[0]

    # Attack 17: Bool True in atom sign
    t3 = copy.deepcopy(w3); t3["atoms"][0]["terms"][0][1] = True
    assert not verify_magic_certificate_exact(t3)[0]

    # Attack 18: Float in atom weight
    t4 = copy.deepcopy(c5); t4["atoms"][0]["weight"] = 0.5
    assert not verify_magic_certificate_exact(t4)[0]

    # Attack 19: Decimal string "0.5" instead of canonical rational "1/2"
    t5 = copy.deepcopy(c5); t5["atoms"][0]["weight"] = "0.5"
    assert not verify_magic_certificate_exact(t5)[0]

    # Attack 20: Negative denominator "1/-2"
    t6 = copy.deepcopy(c5); t6["atoms"][0]["weight"] = "-1/-2"
    assert not verify_magic_certificate_exact(t6)[0]

    # Attack 21: "NaN" / "inf" string
    t7 = copy.deepcopy(w3); t7["claimed_R_lower_bound"] = "NaN"
    assert not verify_magic_certificate_exact(t7)[0]

    # Attack 22: Float hidden in unrelated metadata field (recursive sanitizer!)
    t8 = copy.deepcopy(w3); t8["extra_float"] = 3.14159
    assert not verify_magic_certificate_exact(t8)[0]

    # Attack 23: Bool hidden in unrelated metadata field (recursive sanitizer!)
    t9 = copy.deepcopy(w3); t9["extra_bool"] = False
    assert not verify_magic_certificate_exact(t9)[0]

    # Attack 24: Tuple containing float passed via Python API
    t10 = copy.deepcopy(w3); t10["extra_tuple"] = (1, 2.5)
    assert not verify_magic_certificate_exact(t10)[0]

    # Attack 25: Unreduced rational string "2/4" rejected by canonical parser
    t11 = copy.deepcopy(c5); t11["atoms"][0]["weight"] = "2/4"
    assert not verify_magic_certificate_exact(t11)[0]

    # Attack 26: Leading-zero rational string "001" rejected by canonical parser
    t12 = copy.deepcopy(w3); t12["claimed_beta_plus"] = "001"
    assert not verify_magic_certificate_exact(t12)[0]

    # Attack 27: Non-string float dict key rejected by strict sanitizer
    t13 = copy.deepcopy(w3); t13[1.5] = "x"
    assert not verify_magic_certificate_exact(t13)[0]

    # Attack 28: Decimal("1.5") or set rejected by strict sanitizer
    t14 = copy.deepcopy(w3); t14["extra_dec"] = Decimal("1.5")
    assert not verify_magic_certificate_exact(t14)[0]
    t15 = copy.deepcopy(w3); t15["extra_set"] = {3.14}
    assert not verify_magic_certificate_exact(t15)[0]


# 5. Experiment JSON consistency checks
def test_experiment_json_outputs_consistent():
    exp1 = json.loads((RES_DIR / "exp1_certified_atlas.json").read_text(encoding="utf-8"))
    exp2 = json.loads((RES_DIR / "exp2_majorana_scaling.json").read_text(encoding="utf-8"))
    exp3 = json.loads((RES_DIR / "exp3_w_state_and_ccz.json").read_text(encoding="utf-8"))
    exp4 = json.loads((RES_DIR / "exp4_hierarchy_separation.json").read_text(encoding="utf-8"))
    assert exp1["n_certificates"] == 9
    assert exp1["all_verified_zero_float"] is True
    assert exp2["rows"][-1]["n_qubits"] == 32
    assert exp2["rows"][-1]["R_cert_exact_Q"] > 8.0
    assert exp2["rows"][-1]["pauli_l1_bound"] < 1e-8
    assert exp3["ccz_comparison"]["clique_R_float"] > exp3["ccz_comparison"]["pauli_l1_R_float"]
    assert len(exp4["benchmarks"]) == 6
