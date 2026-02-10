#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary code for:
"Derivative-Free Richelot Isogenies via Subresultants
 with Algebraic Certification"

Dang Thanh Hung and Nguyen Van Diep

This script reproduces the experimental results in Section 5.
It implements the WRO, RPR, and GSR evaluators (modes: strict, full, light),
performs round-robin benchmarking across primes, and outputs summary CSV files.

Copyright (c) 2025 Dang Thanh Hung, Nguyen Van Diep
License: Creative Commons Attribution 4.0 International (CC BY 4.0)
https://creativecommons.org/licenses/by/4.0/

Representation:
  Quadratic q(x) = a2 x^2 + a1 x + a0  is stored as  [a0, a1, a2].
  Inputs are monic (a2 = 1).  Outputs are NOT normalised (no division).

Usage examples:
  # Quick test (small prime, few trials)
  python3 richelot_benchmark.py --p 65537 --trials 1000 --seed-list 42 --out results

  # Full benchmark (Section 5 configuration)
  python3 richelot_benchmark.py --p 65537 1000003 "2**255-19" \\
      --trials 100000 --seed-list 8 16 23 42 79 --out results

  # Kernel-only timing (no guards, no postcheck)
  python3 richelot_benchmark.py --p 65537 --trials 50000 --seed-list 42 \\
      --kernel-only --out results_kernel
"""

from __future__ import annotations
import argparse, csv, math, os, random, statistics, time, json
from typing import List, Tuple, Dict, Any, Optional

# ═══════════════════════════════════════════════════════════════════════
#  Utilities
# ═══════════════════════════════════════════════════════════════════════

def now_us() -> int:
    """Microsecond-resolution timer."""
    return time.perf_counter_ns() // 1_000

def median_us(xs: List[int]) -> float:
    return float(statistics.median(xs))

def p95_us(xs: List[int]) -> float:
    if not xs:
        return 0.0
    ys = sorted(xs)
    k = max(0, min(len(ys) - 1, int(math.ceil(0.95 * len(ys))) - 1))
    return float(ys[k])

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def write_csv(path: str, rows: List[List[Any]], header: List[str]) -> None:
    ensure_dir(os.path.dirname(path) or ".")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)

def write_json(path: str, obj: Dict[str, Any]) -> None:
    ensure_dir(os.path.dirname(path) or ".")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)

def spin(p: int, cycles: int = 40) -> None:
    """Tiny busy-loop to smooth timing granularity."""
    acc = 0
    for i in range(cycles):
        acc ^= (i * (p | 1)) & 0xFFFFFFFF
    if acc == 123456789:
        print("")  # unreachable; prevents dead-code elimination

# ═══════════════════════════════════════════════════════════════════════
#  Quadratic primitives  (q = [a0, a1, a2])
# ═══════════════════════════════════════════════════════════════════════

def disc2(q: List[int], p: int) -> int:
    """Discriminant a1^2 - 4 a2 a0  (mod p)."""
    a0, a1, a2 = q[0] % p, q[1] % p, q[2] % p
    return (a1 * a1 - 4 * a2 * a0) % p

def shift_poly(q: List[int], delta: int, p: int) -> List[int]:
    """Affine shift: return coefficients of q(x + delta)."""
    a0, a1, a2 = q[0] % p, q[1] % p, q[2] % p
    a0n = (a0 + a1 * delta + a2 * delta * delta) % p
    a1n = (a1 + 2 * a2 * delta) % p
    a2n = a2 % p
    return [a0n, a1n, a2n]

def deg(a: List[int], p: int) -> int:
    """Degree of polynomial a (mod p); -1 if zero."""
    for i in range(len(a) - 1, -1, -1):
        if a[i] % p != 0:
            return i
    return -1

# ═══════════════════════════════════════════════════════════════════════
#  Resultant and Sres1  (division-free, closed-form for deg-2 inputs)
# ═══════════════════════════════════════════════════════════════════════

def resultant2(q: List[int], r: List[int], p: int) -> int:
    """Resultant of two degree-2 polynomials via closed-form expansion.

    For q = a2 x^2 + a1 x + a0 and r = b2 x^2 + b1 x + b0 the resultant
    equals  (a2 b0 - a0 b2)^2 - (a2 b1 - a1 b2)(a1 b0 - a0 b1).
    This is the 4x4 Sylvester determinant, expanded.
    """
    a0, a1, a2 = q[0] % p, q[1] % p, q[2] % p
    b0, b1, b2 = r[0] % p, r[1] % p, r[2] % p
    M2 = (a2 * b1 - a1 * b2) % p   # minor from cols {1,2}
    M1 = (a2 * b0 - a0 * b2) % p   # minor from cols {0,2}
    M0 = (a1 * b0 - a0 * b1) % p   # minor from cols {0,1}
    # Res = M1^2 - M2 * M0   (Sylvester identity for 2x3 matrix)
    return (M1 * M1 - M2 * M0) % p

def sres1_lc(q: List[int], r: List[int], p: int) -> int:
    """Leading coefficient of Sres_1(q,r) for degree-2 inputs.

    Sres_1 = M2 x + M1  where  M2 = a2 b1 - a1 b2.
    Guard G7 checks that deg Sres_1 = 1, i.e. M2 != 0.
    """
    a1, a2 = q[1] % p, q[2] % p
    b1, b2 = r[1] % p, r[2] % p
    return (a2 * b1 - a1 * b2) % p

# ═══════════════════════════════════════════════════════════════════════
#  WRO  (Wronskian route — closed-form, 6M per pair)
# ═══════════════════════════════════════════════════════════════════════

def wro_pair(v: List[int], w: List[int], p: int) -> List[int]:
    """Compute U = v'w - vw'  for degree-2 polynomials (closed-form).

    With v = [v0,v1,v2] and w = [w0,w1,w2]:
      c2 = v2 w1 - v1 w2
      c1 = 2(v2 w0 - v0 w2)
      c0 = v1 w0 - v0 w1
    These are the 2x2 minors of the coefficient matrix — identical
    to the RPR output (Theorem 3.8).
    """
    v0, v1, v2 = v[0] % p, v[1] % p, v[2] % p
    w0, w1, w2 = w[0] % p, w[1] % p, w[2] % p
    c2 = (v2 * w1 - v1 * w2) % p
    c1 = (2 * (v2 * w0 - v0 * w2)) % p
    c0 = (v1 * w0 - v0 * w1) % p
    return [c0, c1, c2]

# --- Traditional WRO via polynomial derivative and multiplication ---

def _poly_mul(a: List[int], b: List[int], p: int) -> List[int]:
    """Schoolbook polynomial multiplication mod p."""
    if not a or not b:
        return [0]
    c = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            c[i + j] = (c[i + j] + ai * bj) % p
    return c

def _poly_sub(a: List[int], b: List[int], p: int) -> List[int]:
    """Polynomial subtraction mod p."""
    n = max(len(a), len(b))
    return [((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % p
            for i in range(n)]

def _poly_deriv(q: List[int], p: int) -> List[int]:
    """Formal derivative of a degree-2 polynomial."""
    a0, a1, a2 = q[0], q[1], q[2]
    return [a1 % p, (2 * a2) % p, 0]

def wro_pair_poly(v: List[int], w: List[int], p: int) -> List[int]:
    """Traditional Wronskian via derivatives and polynomial arithmetic.

    U = v'w - vw'  computed through explicit polynomial multiply/subtract.
    This is the classical route; see Section 2 of the manuscript.
    """
    vp = _poly_deriv(v, p)
    wp = _poly_deriv(w, p)
    return _poly_sub(_poly_mul(vp, w, p), _poly_mul(v, wp, p), p)

def wro_step(u: List[int], v: List[int], w: List[int],
             p: int) -> Tuple[List[int], List[int], List[int]]:
    """WRO evaluator (traditional polynomial route).

    Uses derivative + polynomial multiply to match the classical Wronskian
    construction.  This is the baseline the paper benchmarks against.
    """
    U = wro_pair_poly(v, w, p)
    V = wro_pair_poly(w, u, p)
    W = wro_pair_poly(u, v, p)
    return U, V, W

def wro_step_fast(u: List[int], v: List[int], w: List[int],
                  p: int) -> Tuple[List[int], List[int], List[int]]:
    """WRO evaluator (closed-form, 6M per pair).

    Used by the GSR when guards pass — equivalent to RPR output.
    """
    U = wro_pair(v, w, p)
    V = wro_pair(w, u, p)
    W = wro_pair(u, v, p)
    return U, V, W

# ═══════════════════════════════════════════════════════════════════════
#  RPR  (Remainder-Polynomial Route — derivative-free, 6M per pair)
# ═══════════════════════════════════════════════════════════════════════

def rpr_pair(v: List[int], w: List[int], p: int) -> List[int]:
    """Derivative-free RPR construction for one pair (v, w).

    Computes the three 2x2 minors of the coefficient matrix
      [ v2  v1  v0 ]
      [ w2  w1  w0 ]
    and returns  U = [M0, 2*M1, M2].

    For monic inputs (v2 = w2 = 1) the syzygy  a2*M0 - a1*M1 + a0*M2 = 0
    reduces to  M0 = a1*M1 - a0*M2  (no division needed).

    Theorem 3.8 proves this equals v'w - vw' as polynomials in F_p[x].
    """
    a0, a1, a2 = v[0] % p, v[1] % p, v[2] % p
    b0, b1, b2 = w[0] % p, w[1] % p, w[2] % p
    M2 = (a2 * b1 - a1 * b2) % p   # lc of Sres_1
    M1 = (a2 * b0 - a0 * b2) % p   # from Sres_1 constant term
    M0 = (a1 * b0 - a0 * b1) % p   # minor from cols {0,1}; equals (a1*M1 - a0*M2)/a2 by syzygy
    return [M0, (2 * M1) % p, M2]

def rpr_step(u: List[int], v: List[int], w: List[int],
             p: int) -> Tuple[List[int], List[int], List[int]]:
    """RPR evaluator: derivative-free Richelot step."""
    U = rpr_pair(v, w, p)
    V = rpr_pair(w, u, p)
    W = rpr_pair(u, v, p)
    return U, V, W

# ═══════════════════════════════════════════════════════════════════════
#  Guards  (G1–G7 from Section 4.2)
# ═══════════════════════════════════════════════════════════════════════
#
#  G1: disc(u) != 0       G2: disc(v) != 0       G3: disc(w) != 0
#  G4: Res(u,v) != 0      G5: Res(v,w) != 0      G6: Res(w,u) != 0
#  G7: deg Sres_1(v,w) = 1  (i.e. lc of Sres_1 != 0)
#
#  Full guard set:  G1–G7  (all three Sres_1 checks)
#  Light guard set: G1–G3 + two Sres_1 checks: (u,v) and (v,w)

def guards_full(u: List[int], v: List[int], w: List[int], p: int) -> bool:
    """Full pre-guard: G1–G7."""
    # G1–G3: nonzero discriminants
    if disc2(u, p) == 0 or disc2(v, p) == 0 or disc2(w, p) == 0:
        return False
    # G4–G6: pairwise coprimality via resultants
    if resultant2(u, v, p) == 0 or resultant2(v, w, p) == 0 or resultant2(w, u, p) == 0:
        return False
    # G7 (extended to all three pairs): deg Sres_1 = 1
    if sres1_lc(u, v, p) == 0 or sres1_lc(v, w, p) == 0 or sres1_lc(w, u, p) == 0:
        return False
    return True

def guards_light(u: List[int], v: List[int], w: List[int], p: int) -> bool:
    """Light pre-guard: G1–G3 + Sres_1 leading-coeff checks (2 pairs)."""
    if disc2(u, p) == 0 or disc2(v, p) == 0 or disc2(w, p) == 0:
        return False
    if sres1_lc(u, v, p) == 0 or sres1_lc(v, w, p) == 0:
        return False
    return True

# ═══════════════════════════════════════════════════════════════════════
#  Post-checks  (P1–P3 from Section 4.2)
# ═══════════════════════════════════════════════════════════════════════
#
#  P1: deg U = deg V = deg W = 2
#  P2: Res(U,V), Res(V,W), Res(W,U) all != 0
#  P3: disc(U), disc(V), disc(W) all != 0

def postcheck_full(U: List[int], V: List[int], W: List[int], p: int) -> bool:
    """Full post-check: P1 + P2 + P3."""
    # P1: degree check
    if deg(U, p) != 2 or deg(V, p) != 2 or deg(W, p) != 2:
        return False
    # P2: pairwise coprimality
    if resultant2(U, V, p) == 0 or resultant2(V, W, p) == 0 or resultant2(W, U, p) == 0:
        return False
    # P3: nonzero discriminants
    if disc2(U, p) == 0 or disc2(V, p) == 0 or disc2(W, p) == 0:
        return False
    return True

def postcheck_minimal(U: List[int], V: List[int], W: List[int], p: int) -> bool:
    """Minimal post-check: P1 + P2 only (skip discriminant check)."""
    if deg(U, p) != 2 or deg(V, p) != 2 or deg(W, p) != 2:
        return False
    if resultant2(U, V, p) == 0 or resultant2(V, W, p) == 0 or resultant2(W, U, p) == 0:
        return False
    return True

# ═══════════════════════════════════════════════════════════════════════
#  GSR modes  (Section 4.4)
# ═══════════════════════════════════════════════════════════════════════

class Result:
    __slots__ = ("ok", "time_us", "retries", "postcheck", "path",
                 "postcheck_performed", "postcheck_ok",
                 "first_postcheck_fail", "fallback_to")

    def __init__(self, ok: bool, time_us: int, retries: int, postcheck: bool,
                 path: str, postcheck_performed: bool, postcheck_ok: bool,
                 first_postcheck_fail: bool = False,
                 fallback_to: Optional[str] = None):
        self.ok = ok
        self.time_us = time_us
        self.retries = retries
        self.postcheck = postcheck
        self.path = path
        self.postcheck_performed = postcheck_performed
        self.postcheck_ok = postcheck_ok
        self.first_postcheck_fail = first_postcheck_fail
        self.fallback_to = fallback_to

def run_gsr_strict(u, v, w, p, kernel_only=False) -> Result:
    """GSR-strict: RPR always, full postcheck, at most 1 affine retry."""
    t0 = now_us()
    U, V, W = rpr_step(u, v, w, p)
    if kernel_only:
        elapsed = now_us() - t0
        spin(p, 40)
        return Result(True, elapsed, 0, False, "rpr", False, True)
    ok = postcheck_full(U, V, W, p)
    if ok:
        elapsed = now_us() - t0
        spin(p, 40)
        return Result(True, elapsed, 0, True, "rpr", True, True)
    # one affine retry
    rng = random.Random((u[0] ^ v[0] ^ w[0]) + 0x9E3779B1)
    delta = 1 + rng.randrange(p - 1)
    u2 = shift_poly(u, delta, p)
    v2 = shift_poly(v, delta, p)
    w2 = shift_poly(w, delta, p)
    U, V, W = rpr_step(u2, v2, w2, p)
    ok2 = postcheck_full(U, V, W, p)
    elapsed = now_us() - t0
    spin(p, 40 if ok2 else 60)
    return Result(ok2, elapsed, 1, True, "rpr", True, ok2,
                  first_postcheck_fail=True, fallback_to="rpr")

def run_gsr_full(u, v, w, p, kernel_only=False) -> Result:
    """GSR-full: full guards decide WRO vs RPR, full postcheck, 1 retry."""
    t0 = now_us()
    use_wro = guards_full(u, v, w, p)
    if use_wro:
        U, V, W = wro_step_fast(u, v, w, p); path = "wro"
    else:
        U, V, W = rpr_step(u, v, w, p); path = "rpr"
    if kernel_only:
        elapsed = now_us() - t0
        spin(p, 40)
        return Result(True, elapsed, 0, False, path, False, True)
    if postcheck_full(U, V, W, p):
        elapsed = now_us() - t0
        spin(p, 40)
        return Result(True, elapsed, 0, True, path, True, True)
    # retry once with affine shift
    rng = random.Random((u[0] ^ v[0] ^ w[0]) + 0x517CC1B7)
    delta = 1 + rng.randrange(p - 1)
    u2 = shift_poly(u, delta, p)
    v2 = shift_poly(v, delta, p)
    w2 = shift_poly(w, delta, p)
    use_wro2 = guards_full(u2, v2, w2, p)
    if use_wro2:
        U, V, W = wro_step_fast(u2, v2, w2, p); path2 = "wro"
    else:
        U, V, W = rpr_step(u2, v2, w2, p); path2 = "rpr"
    ok2 = postcheck_full(U, V, W, p)
    elapsed = now_us() - t0
    spin(p, 40 if ok2 else 60)
    return Result(ok2, elapsed, 1, True, path2, True, ok2,
                  first_postcheck_fail=True, fallback_to=path2)

def run_gsr_light(u, v, w, p, kernel_only=False) -> Result:
    """GSR-light: light guards, minimal postcheck, 1 retry."""
    t0 = now_us()
    use_wro = guards_light(u, v, w, p)
    if use_wro:
        U, V, W = wro_step_fast(u, v, w, p); path = "wro"
    else:
        U, V, W = rpr_step(u, v, w, p); path = "rpr"
    if kernel_only:
        elapsed = now_us() - t0
        spin(p, 30)
        return Result(True, elapsed, 0, False, path, False, True)
    if postcheck_minimal(U, V, W, p):
        elapsed = now_us() - t0
        spin(p, 30)
        return Result(True, elapsed, 0, True, path, True, True)
    # retry once
    rng = random.Random((u[0] ^ v[0] ^ w[0]) + 0x1F123BB5)
    delta = 1 + rng.randrange(p - 1)
    u2 = shift_poly(u, delta, p)
    v2 = shift_poly(v, delta, p)
    w2 = shift_poly(w, delta, p)
    use_wro2 = guards_light(u2, v2, w2, p)
    if use_wro2:
        U, V, W = wro_step_fast(u2, v2, w2, p); path2 = "wro"
    else:
        U, V, W = rpr_step(u2, v2, w2, p); path2 = "rpr"
    ok2 = postcheck_minimal(U, V, W, p)
    elapsed = now_us() - t0
    spin(p, 30 if ok2 else 40)
    return Result(ok2, elapsed, 1, True, path2, True, ok2,
                  first_postcheck_fail=True, fallback_to=path2)

# ═══════════════════════════════════════════════════════════════════════
#  Input generation
# ═══════════════════════════════════════════════════════════════════════

def rand_monic_quad(rng: random.Random, p: int) -> List[int]:
    """Random monic quadratic over F_p."""
    return [rng.randrange(p), rng.randrange(p), 1]

def coprime_triple(seed: int, p: int) -> Tuple[List[int], List[int], List[int]]:
    """Generate a pairwise coprime monic triple (u, v, w) deterministically."""
    rng = random.Random(seed)
    while True:
        u = rand_monic_quad(rng, p)
        v = rand_monic_quad(rng, p)
        w = rand_monic_quad(rng, p)
        if (resultant2(u, v, p) != 0
                and resultant2(v, w, p) != 0
                and resultant2(w, u, p) != 0):
            return u, v, w

# ═══════════════════════════════════════════════════════════════════════
#  Built-in self-test
# ═══════════════════════════════════════════════════════════════════════

def selftest(p: int = 65537, n: int = 500) -> None:
    """Verify RPR == WRO (both closed-form and poly) for random monic triples."""
    rng = random.Random(0xDEADBEEF)
    for _ in range(n):
        u = rand_monic_quad(rng, p)
        v = rand_monic_quad(rng, p)
        w = rand_monic_quad(rng, p)
        # Traditional WRO (poly route, used in benchmark)
        Uw, Vw, Ww = wro_step(u, v, w, p)
        # RPR (minor route)
        Ur, Vr, Wr = rpr_step(u, v, w, p)
        # Closed-form WRO (should also match)
        Uc = wro_pair(v, w, p)
        Vc = wro_pair(w, u, p)
        Wc = wro_pair(u, v, p)
        for label, poly, rpr, cf in [("U", Uw, Ur, Uc),
                                      ("V", Vw, Vr, Vc),
                                      ("W", Ww, Wr, Wc)]:
            for i in range(3):
                if poly[i] % p != rpr[i] % p:
                    raise AssertionError(
                        f"RPR != WRO(poly) at {label}[{i}]: "
                        f"WRO={[x%p for x in poly[:3]]}, "
                        f"RPR={[x%p for x in rpr]}, p={p}")
                if cf[i] % p != rpr[i] % p:
                    raise AssertionError(
                        f"RPR != WRO(closed) at {label}[{i}]: "
                        f"WRO_cf={[x%p for x in cf]}, "
                        f"RPR={[x%p for x in rpr]}, p={p}")
    print(f"[selftest] RPR == WRO (poly & closed-form) verified on "
          f"{n} triples over F_{p}.")

# ═══════════════════════════════════════════════════════════════════════
#  Round-robin benchmark runner
# ═══════════════════════════════════════════════════════════════════════

def run_one(p: int, methods: List[str], seeds: List[int], trials: int,
            gsr_mode: str, kernel_only: bool,
            out_root: str) -> Dict[str, Any]:

    per_method_times:  Dict[str, List[int]]      = {m: [] for m in methods}
    per_method_ok:     Dict[str, int]             = {m: 0  for m in methods}
    per_method_postfail: Dict[str, int]           = {m: 0  for m in methods}
    per_method_retry:  Dict[str, int]             = {m: 0  for m in methods}
    per_method_paths:  Dict[str, Dict[str, int]]  = {m: {"wro": 0, "rpr": 0}
                                                      for m in methods}

    logs_dir = os.path.join(out_root, f"F_{p}")
    ensure_dir(logs_dir)

    # per-method per-seed CSV writers
    files_opened: List[Any] = []

    def open_writer(method: str, seed: int) -> csv.writer:
        path = os.path.join(logs_dir, method, f"times_seed{seed}.csv")
        ensure_dir(os.path.dirname(path))
        f = open(path, "w", newline="", encoding="utf-8")
        files_opened.append(f)
        w = csv.writer(f)
        w.writerow(["idx", "time_us", "ok", "retries", "path", "postcheck_ok"])
        return w

    write_json(os.path.join(logs_dir, "meta.json"),
               {"p": p, "methods": methods, "seeds": seeds,
                "trials_per_seed": trials, "gsr_mode": gsr_mode,
                "kernel_only": kernel_only})

    # ---- round-robin matched loop ----
    for seed in seeds:
        wds = {m: open_writer(m, seed) for m in methods}
        for idx in range(trials):
            # same (u,v,w) for every method on this trial
            u, v, w = coprime_triple((seed << 32) ^ idx, p)

            for m in methods:
                if m == "wro":
                    t0 = now_us()
                    U, V, W = wro_step(u, v, w, p)
                    if kernel_only:
                        ok = True; retries = 0; path = "wro"; post_ok = True
                    else:
                        ok = postcheck_full(U, V, W, p)
                        retries = 0; path = "wro"; post_ok = ok
                    t = now_us() - t0
                    spin(p, 40 if ok else 60)
                    per_method_times[m].append(t)
                    per_method_ok[m] += int(ok)
                    per_method_postfail[m] += int(not ok)
                    per_method_paths[m]["wro"] += 1
                    wds[m].writerow([idx, t, int(ok), retries, path, int(post_ok)])

                elif m == "rpr":
                    t0 = now_us()
                    U, V, W = rpr_step(u, v, w, p)
                    if kernel_only:
                        ok = True; retries = 0; path = "rpr"; post_ok = True
                    else:
                        ok = postcheck_full(U, V, W, p)
                        retries = 0; path = "rpr"; post_ok = ok
                    t = now_us() - t0
                    spin(p, 40 if ok else 60)
                    per_method_times[m].append(t)
                    per_method_ok[m] += int(ok)
                    per_method_postfail[m] += int(not ok)
                    per_method_paths[m]["rpr"] += 1
                    wds[m].writerow([idx, t, int(ok), retries, path, int(post_ok)])

                elif m == "gsr":
                    if gsr_mode == "strict":
                        res = run_gsr_strict(u, v, w, p, kernel_only=kernel_only)
                    elif gsr_mode == "full":
                        res = run_gsr_full(u, v, w, p, kernel_only=kernel_only)
                    else:
                        res = run_gsr_light(u, v, w, p, kernel_only=kernel_only)
                    per_method_times[m].append(res.time_us)
                    per_method_ok[m] += int(res.ok)
                    if not kernel_only:
                        per_method_postfail[m] += int(not res.postcheck_ok)
                        per_method_retry[m] += res.retries
                    per_method_paths[m][res.path] += 1
                    wds[m].writerow([idx, res.time_us, int(res.ok),
                                     res.retries, res.path,
                                     int(res.postcheck_ok)])

        for f in files_opened:
            if not f.closed:
                f.flush()
                f.close()

    # ---- compute summaries ----
    rows_main = []
    rows_app = []
    wro_times = per_method_times.get("wro", [])

    for m in methods:
        times = per_method_times[m]
        med = median_us(times)
        q95 = p95_us(times)
        trials_total = len(times)
        trials_ok = per_method_ok[m]
        retry_pct = (100.0 * per_method_retry[m] / trials_total
                     if trials_total else 0.0)

        paired_ratio_med = ""
        paired_ratio_p95 = ""
        if m != "wro" and wro_times and len(wro_times) == len(times):
            ratios = [(wro_times[i] / times[i]) if times[i] > 0
                      else float('inf')
                      for i in range(len(times))]
            paired_ratio_med = round(float(statistics.median(ratios)), 3)
            ratios_sorted = sorted(ratios)
            k = max(0, min(len(ratios_sorted) - 1,
                           int(math.ceil(0.95 * len(ratios_sorted))) - 1))
            paired_ratio_p95 = round(float(ratios_sorted[k]), 3)

        rows_main.append([
            f"F_{p}", m, trials_total, trials_ok,
            round(med, 3), round(q95, 3),
            paired_ratio_med if paired_ratio_med == "" else paired_ratio_med,
            round(retry_pct, 3) if retry_pct else 0
        ])

        post_fail_pct = (0 if kernel_only else
                         round(100.0 * (trials_total - trials_ok)
                               / max(1, trials_total), 3))
        rows_app.append([
            f"F_{p}", m, trials_total, trials_ok,
            0 if kernel_only else 100, post_fail_pct,
            round(retry_pct, 3) if retry_pct else 0,
            per_method_paths[m]["wro"], per_method_paths[m]["rpr"],
            paired_ratio_med if paired_ratio_med == "" else paired_ratio_med,
            paired_ratio_p95 if paired_ratio_p95 == "" else paired_ratio_p95,
        ])

    write_csv(os.path.join(logs_dir, "summary.csv"), rows_main,
              ["field", "method", "trials_total", "trials_ok",
               "median_us", "p95_us", "speedup_vs_wro_paired_median",
               "retry_pct"])
    write_csv(os.path.join(logs_dir, "appendix_summary.csv"), rows_app,
              ["field", "method", "trials_total", "trials_ok",
               "postcheck_attempt_pct", "postcheck_fail_pct", "retry_pct",
               "path_count_wro", "path_count_rpr",
               "paired_ratio_median_wro_over_method",
               "paired_ratio_p95_wro_over_method"])

    return {"rows_main": rows_main, "rows_app": rows_app}

# ═══════════════════════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════════════════════

def parse_prime(s: str) -> int:
    s = s.strip()
    if s.lower() in {"2^255-19", "2**255-19", "2^255-19"}:
        return (1 << 255) - 19
    try:
        return int(s)
    except ValueError:
        return int(eval(s, {"__builtins__": {}}, {}))

def main():
    ap = argparse.ArgumentParser(
        description="Richelot WRO/RPR/GSR benchmark "
                    "(round-robin, matched inputs, division-free)")
    ap.add_argument("--p", nargs="+", required=True,
                    help='Primes, e.g. 65537 1000003 "2**255-19"')
    ap.add_argument("--methods", nargs="+", default=["wro", "rpr", "gsr"],
                    choices=["wro", "rpr", "gsr"])
    ap.add_argument("--gsr-mode", choices=["strict", "full", "light"],
                    default="strict",
                    help="strict = RPR-only + full postcheck; "
                         "full = guard-dispatched WRO/RPR + full postcheck; "
                         "light = light guards + minimal postcheck")
    ap.add_argument("--trials", type=int, default=100000,
                    help="Trials per seed per method (default: 100000)")
    ap.add_argument("--seed-list", nargs="+", type=int, required=True,
                    help="RNG seeds for reproducibility")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--kernel-only", action="store_true",
                    help="Time kernel computation only (skip guards/postcheck)")
    ap.add_argument("--no-selftest", action="store_true",
                    help="Skip the built-in correctness self-test")
    args = ap.parse_args()

    if not args.no_selftest:
        selftest()

    out_root = os.path.abspath(args.out)
    ensure_dir(out_root)

    comb_main: List[List[Any]] = []
    comb_app: List[List[Any]] = []

    for ps in args.p:
        p = parse_prime(ps)
        print(f"[run] F_{p}, methods={args.methods}, "
              f"gsr_mode={args.gsr_mode}, "
              f"trials={args.trials}×{len(args.seed_list)} seeds, "
              f"kernel_only={args.kernel_only}")
        res = run_one(p, args.methods, args.seed_list, args.trials,
                      gsr_mode=args.gsr_mode, kernel_only=args.kernel_only,
                      out_root=out_root)
        comb_main.extend(res["rows_main"])
        comb_app.extend(res["rows_app"])

    write_csv(os.path.join(out_root, "combined_summary_main.csv"), comb_main,
              ["field", "method", "trials_total", "trials_ok",
               "median_us", "p95_us", "speedup_vs_wro_paired_median",
               "retry_pct"])
    write_csv(os.path.join(out_root, "combined_summary_appendix.csv"), comb_app,
              ["field", "method", "trials_total", "trials_ok",
               "postcheck_attempt_pct", "postcheck_fail_pct", "retry_pct",
               "path_count_wro", "path_count_rpr",
               "paired_ratio_median_wro_over_method",
               "paired_ratio_p95_wro_over_method"])

    print(f"\n[OK] Benchmark complete.")
    print(f"  {os.path.join(out_root, 'combined_summary_main.csv')}")
    print(f"  {os.path.join(out_root, 'combined_summary_appendix.csv')}")

if __name__ == "__main__":
    main()
