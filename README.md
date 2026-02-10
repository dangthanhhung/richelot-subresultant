# richelot-subresultant

Supplementary code for the paper:

> **Derivative-Free Richelot Isogenies via Subresultants with Algebraic Certification**
> Hung T. Dang, Diep V. Nguyen
> Submitted to *Finite Fields and Their Applications*, 2026.

## Contents

| File | Description |
|------|-------------|
| `benchmark.py` | Benchmark runner: kernel-only, GSR-strict, GSR-full, GSR-light, and multi-prime scaling configurations (§5.2–5.6) |
| `verify.py` | Correctness verification: coefficient-wise RPR vs WRO comparison and double-Richelot involution test (§5.1) |
| `cross_validate.sage` | Independent SageMath cross-validation of the Richelot step |

## Requirements

- Python ≥ 3.13
- SageMath ≥ 9.5 (for `cross_validate.sage` only)

No external Python packages are required.

## Usage

```bash
# Benchmarks: WRO, RPR, GSR over three primes
# (p = 65537, 10^6+3, 2^255-19; 10^6 matched trials each)
python benchmark.py

# Correctness tests: 2.5 × 10^5 random triples across five primes
# (p = 101, 65537, 10^6+3, 2^31-1, 2^61-1)
python verify.py

# SageMath cross-validation
sage cross_validate.sage
```

## Routes

The code implements three routes for a single Richelot (2,2)-isogeny step on monic quadratic triples over F_p, p > 2:

| Route | Description |
|-------|-------------|
| **WRO** | Classical Wronskian: U = v'w − vw' via explicit derivatives |
| **RPR** | Remainder-Polynomial Route: derivative-free via pseudo-remainder and minor syzygy |
| **GSR** | Guarded Subresultant Route: RPR core with algebraic guards (G1–G7) and post-check (P1–P3) |

## License

MIT. See [LICENSE](LICENSE).

## Citation

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
