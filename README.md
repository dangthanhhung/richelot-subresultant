# Supplementary Code

**Paper:** Derivative-Free Richelot Isogenies via Subresultants
with Algebraic Certification

**Authors:** Hung T. Dang, Diep V. Nguyen

## Contents

| File | Description |
|---|---|
| `richelot_benchmark.py` | Benchmark runner for WRO, RPR, GSR (Section 5) |
| `verify_correctness.py` | Standalone correctness verification (no SageMath) |
| `verify_sage.sage` | SageMath cross-validation script |
| `run_all_benchmarks.sh` | Run all benchmarks (Linux/macOS) |
| `run_all_benchmarks.bat` | Run all benchmarks (Windows) |

## Requirements

- Python 3.10+ (standard library only; no external packages)
- SageMath 9.0+ (only for `verify_sage.sage`)

## Quick Start

### Correctness verification

```bash
python3 verify_correctness.py
```

Verifies three properties over random monic coprime triples across
five primes (p = 101, 65537, 10^6+3, 2^31-1, 2^61-1):

1. **RPR == WRO** (Theorem 3.8): exact polynomial identity
2. **Double-Richelot involution**: Richelot(Richelot(u,v,w)) is proportional
   to (u,v,w), proving the construction computes a valid (2,2)-isogeny
3. **Output admissibility**: degree, coprimality, discriminant checks

### SageMath cross-validation

```bash
sage verify_sage.sage
```

Cross-validates our closed-form WRO and RPR against SageMath's
symbolic polynomial arithmetic.

### Full benchmark suite

```bash
# Linux/macOS
chmod +x run_all_benchmarks.sh
./run_all_benchmarks.sh

# Windows
run_all_benchmarks.bat
```

Runs all five benchmark configurations from Section 5:
1. Kernel-only (Table 3)
2. GSR-strict (Table 4, GSR-strict rows)
3. GSR-full (Table 4, GSR-full rows)
4. GSR-light (Table 4, GSR-light rows)
5. Multi-prime scaling (Table 6)

### Single benchmark example

```bash
# Quick test
python3 richelot_benchmark.py --p 65537 --trials 1000 --seed-list 42 --out results

# Kernel-only timing
python3 richelot_benchmark.py --p 65537 1000003 "2**255-19" \
    --trials 100000 --seed-list 8 16 23 42 79 \
    --kernel-only --out results_kernel
```

## Output Structure

```
bench_kernel/
  combined_summary_main.csv       # Aggregate results (Tables 3, 4, 6)
  combined_summary_appendix.csv   # Retry rates, path dispatch counts
  F_65537/
    meta.json                     # Run configuration
    summary.csv                   # Per-prime summary
    wro/times_seed42.csv          # Per-trial raw timings
    rpr/times_seed42.csv
    gsr/times_seed42.csv
```

## License

Creative Commons Attribution 4.0 International (CC BY 4.0)
