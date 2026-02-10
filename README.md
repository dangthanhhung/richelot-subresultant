# richelot-subresultant

Python and SageMath code accompanying the paper:

> **Derivative-Free Richelot Isogenies via Subresultants with Algebraic Certification**
> Hung T. Dang, Diep V. Nguyen
> Submitted to *Finite Fields and Their Applications*, 2026.

## Contents

- `benchmark.py` — Benchmark runner for WRO, RPR, and GSR routes ($10^6$ matched trials per prime).
- `verify.py` — Correctness verification: coefficient-wise RPR vs WRO comparison and double-Richelot involution test.
- `cross_validate.sage` — Independent SageMath cross-validation script.

## Requirements

- Python ≥ 3.10
- SageMath ≥ 9.5 (for `cross_validate.sage` only)

No external Python packages are required; all arithmetic uses built-in integers and the `random` module.

## Usage

```bash
# Run benchmarks (default: p = 65537, 10^6+3, 2^255-19)
python benchmark.py

# Run correctness tests (five primes, 2.5 × 10^5 triples)
python verify.py

# SageMath cross-validation
sage cross_validate.sage
```

## Overview

The code implements three routes for a single Richelot (2,2)-isogeny step on monic quadratic triples over F_p (p > 2):

| Route | Description |
|-------|-------------|
| **WRO** | Classical Wronskian: $U = v'w - vw'$ via explicit derivatives |
| **RPR** | Remainder-Polynomial Route: derivative-free, uses pseudo-remainder + minor syzygy |
| **GSR** | Guarded Subresultant Route: RPR core with algebraic guards and post-check |

See the paper for full details on the algebraic equivalence (Theorem 3.8) and the certified evaluator (Algorithm 1).

## License

MIT. See [LICENSE](LICENSE).

## Citation

If you use this code, please cite:

```bibtex
@article{DangNguyen2026,
  title   = {Derivative-Free Richelot Isogenies via Subresultants
             with Algebraic Certification},
  author  = {Dang, Hung T. and Nguyen, Diep V.},
  journal = {Finite Fields and Their Applications},
  year    = {2026},
  note    = {Submitted}
}
```
