"""
Pure-Python Zero-Floating-Point Exact Symplectic Pauli Algebra Engine (arbitrary n qubits).

Conventions:
  - An n-qubit Hermitian Pauli string P(x, z) is represented as a pair of non-negative Python ints
    (x, z) with 0 <= x, z < 2^n, where bit k of x and bit k of z encode qubit k:
      (0, 0) -> I,  (1, 0) -> X,  (0, 1) -> Z,  (1, 1) -> Y = i X Z.
  - Symplectic inner product over F_2:
      omega((x1, z1), (x2, z2)) = popcount((x1 & z2) ^ (z1 & x2)) mod 2.
    omega = 0 iff P1, P2 commute; omega = 1 iff P1, P2 anticommute.
  - For commuting Paulis P(u), P(v), their exact Hermitian product satisfies P(u) P(v) = sigma * P(u ^ v)
    with sigma in {-1, +1} computed strictly via exact integer bitcounts.
  - Includes exact Jordan-Wigner Majorana clique generator K_{2n+1} for arbitrary n >= 1,
    and exact stabilizer-state enumerator for n = 1 (6), n = 2 (60), and n = 3 (1080).
"""
from __future__ import annotations

import itertools

Pauli = tuple[int, int]


def symp_inner(u: Pauli, v: Pauli) -> int:
    """Return exact symplectic inner product omega(u, v) in {0, 1}."""
    return ((u[0] & v[1]) ^ (u[1] & v[0])).bit_count() & 1


def pauli_mul_commuting_sign(u: Pauli, v: Pauli) -> tuple[int, Pauli]:
    """
    For commuting Hermitian Paulis P(u), P(v) (omega(u, v) == 0), return (sigma, u ^ v)
    where sigma in {-1, +1} satisfies P(u) P(v) = sigma * P(u ^ v) identically.
    Raises ValueError if u and v anticommute.
    """
    if symp_inner(u, v) != 0:
        raise ValueError(f"Paulis {u} and {v} anticommute; product is not Hermitian +/- P")
    ux, uz = u
    vx, vz = v
    wx, wz = ux ^ vx, uz ^ vz
    # H(x, z) = i^{|x & z|} X^x Z^z
    # H(u) H(v) = i^{|ux&uz| + |vx&vz| - |wx&wz| + 2*|uz&vx|} H(u^v)
    exp_i = (
        (ux & uz).bit_count()
        + (vx & vz).bit_count()
        - (wx & wz).bit_count()
        + 2 * (uz & vx).bit_count()
    ) % 4
    if exp_i == 0:
        return 1, (wx, wz)
    if exp_i == 2:
        return -1, (wx, wz)
    raise RuntimeError(f"Unexpected phase i^{exp_i} for commuting Paulis {u}, {v}")


def pauli_to_str(u: Pauli, n_qubits: int) -> str:
    """Convert bitmask Pauli (x, z) to human-readable string of length n_qubits (qubit 0 on left)."""
    ux, uz = u
    chars = []
    for k in range(n_qubits):
        bx = (ux >> k) & 1
        bz = (uz >> k) & 1
        chars.append(["I", "Z", "X", "Y"][2 * bx + bz])
    return "".join(chars)


def str_to_pauli(s: str) -> Pauli:
    """Parse Pauli string like 'XYZI' (qubit 0 on left) into exact (x, z) bitmask pair."""
    x = 0
    z = 0
    for k, ch in enumerate(s):
        if ch == "I":
            pass
        elif ch == "X":
            x |= 1 << k
        elif ch == "Z":
            z |= 1 << k
        elif ch == "Y":
            x |= 1 << k
            z |= 1 << k
        else:
            raise ValueError(f"Invalid Pauli character {ch!r} in {s!r}")
    return (x, z)


def pauli_weight(u: Pauli) -> int:
    """Number of qubits on which P(u) acts nontrivially."""
    return (u[0] | u[1]).bit_count()


def jordan_wigner_clique(n_qubits: int) -> list[Pauli]:
    """
    Construct the canonical maximal anticommuting Jordan-Wigner Majorana clique K_{2n+1}
    on n_qubits >= 1 qubits:
      Gamma_{2k}   = Z_0 ... Z_{k-1} X_k         (k = 0..n-1)
      Gamma_{2k+1} = Z_0 ... Z_{k-1} Y_k         (k = 0..n-1)
      Gamma_{2n}   = Z_0 ... Z_{n-1}
    """
    if n_qubits < 1:
        raise ValueError("n_qubits must be >= 1")
    clq: list[Pauli] = []
    for k in range(n_qubits):
        z_pref = (1 << k) - 1
        clq.append((1 << k, z_pref))
        clq.append((1 << k, z_pref | (1 << k)))
    clq.append((0, (1 << n_qubits) - 1))
    return clq


def enumerate_stabilizer_states_exact(n_qubits: int) -> list[dict[Pauli, int]]:
    """
    Exact enumeration of ALL pure n-qubit stabilizer states |STAB_n| for n in {1, 2, 3}
    (6, 60, and 1080 states respectively). Each state is returned as a dict mapping every
    Pauli u in its abelian stabilizer subgroup S_phi (|S_phi| = 2^n) to its exact sign in {-1, +1}.
    """
    if n_qubits not in (1, 2, 3):
        raise ValueError("Exact enumeration intended for n_qubits in {1, 2, 3}")
    N = 1 << (2 * n_qubits)
    mask = (1 << n_qubits) - 1
    all_v = [(k >> n_qubits, k & mask) for k in range(1, N)]
    subspaces: set[frozenset[Pauli]] = set()

    def rec(basis: list[Pauli], span: set[Pauli]) -> None:
        if len(basis) == n_qubits:
            subspaces.add(frozenset(span))
            return
        for v in all_v:
            if v in span:
                continue
            if basis and v < basis[-1]:
                continue
            if any(symp_inner(v, g) for g in basis):
                continue
            ns = set(span)
            for s in span:
                ns.add((s[0] ^ v[0], s[1] ^ v[1]))
            rec(basis + [v], ns)

    rec([], {(0, 0)})
    states: list[dict[Pauli, int]] = []
    for S in sorted(subspaces, key=lambda fs: sorted(fs)):
        basis: list[Pauli] = []
        span = {(0, 0)}
        for v in sorted(S):
            if v not in span:
                basis.append(v)
                span = span | {(s[0] ^ v[0], s[1] ^ v[1]) for s in span}
        for signs in itertools.product((1, -1), repeat=n_qubits):
            elems: dict[Pauli, int] = {(0, 0): 1}
            for g, sg in zip(basis, signs):
                ne = dict(elems)
                for u, su in elems.items():
                    sig, uv = pauli_mul_commuting_sign(u, g)
                    ne[uv] = su * sg * sig
                elems = ne
            states.append(elems)
    return states
