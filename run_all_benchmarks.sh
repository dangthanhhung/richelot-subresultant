#!/usr/bin/env bash
# ================================================================
#  run_all_benchmarks.sh
#  Runs all 5 benchmark configurations for Section 5 of the paper.
#  Total: ~3M trials x 5 configs. Expect 30-90 minutes depending on CPU.
#  Run from the directory containing richelot_benchmark.py.
#
#  Usage:  chmod +x run_all_benchmarks.sh && ./run_all_benchmarks.sh
# ================================================================

set -e

echo "============================================================"
echo "  Richelot Benchmark Suite"
echo "  Paper: Derivative-Free Richelot Isogenies via Subresultants"
echo "         with Algebraic Certification"
echo "============================================================"
echo

# --- Selftest (run once) ---
echo "[0/5] Running selftest..."
python3 richelot_benchmark.py --p 65537 --trials 1 --seed-list 1 --out _selftest_tmp
rm -rf _selftest_tmp
echo

# --- Benchmark 1: Kernel-only (core WRO vs RPR speedup) ---
echo "[1/5] Kernel-only: WRO vs RPR vs GSR (no guards/postcheck)"
echo "      Primes: 65537, 1000003, 2^255-19"
echo "      Trials: 200000 x 5 seeds = 1M per prime"
python3 richelot_benchmark.py --p 65537 1000003 "2**255-19" \
    --trials 200000 --seed-list 8 16 23 42 79 \
    --kernel-only --no-selftest --out bench_kernel
echo

# --- Benchmark 2: GSR strict (RPR-only, full postcheck) ---
echo "[2/5] GSR strict: RPR path + full postcheck + retry"
python3 richelot_benchmark.py --p 65537 1000003 "2**255-19" \
    --trials 200000 --seed-list 8 16 23 42 79 \
    --gsr-mode strict --no-selftest --out bench_strict
echo

# --- Benchmark 3: GSR full (guard-dispatched) ---
echo "[3/5] GSR full: full guards dispatch WRO/RPR + full postcheck"
python3 richelot_benchmark.py --p 65537 1000003 "2**255-19" \
    --trials 200000 --seed-list 8 16 23 42 79 \
    --gsr-mode full --no-selftest --out bench_full
echo

# --- Benchmark 4: GSR light (minimal guards) ---
echo "[4/5] GSR light: light guards + minimal postcheck"
python3 richelot_benchmark.py --p 65537 1000003 "2**255-19" \
    --trials 200000 --seed-list 8 16 23 42 79 \
    --gsr-mode light --no-selftest --out bench_light
echo

# --- Benchmark 5: Multi-prime scaling (kernel-only) ---
echo "[5/5] Multi-prime scaling: 5 primes, kernel-only"
python3 richelot_benchmark.py --p 101 65537 1000003 "2**127-1" "2**255-19" \
    --trials 100000 --seed-list 8 16 23 42 79 \
    --kernel-only --no-selftest --out bench_primes
echo

echo "============================================================"
echo "  All benchmarks complete!  Output directories:"
echo
echo "    bench_kernel/combined_summary_main.csv     -- Table 3 (kernel speedup)"
echo "    bench_strict/combined_summary_main.csv     -- Table 4a (GSR strict)"
echo "    bench_full/combined_summary_main.csv       -- Table 4b (GSR full)"
echo "    bench_light/combined_summary_main.csv      -- Table 4c (GSR light)"
echo "    bench_primes/combined_summary_main.csv     -- Table 5 (prime scaling)"
echo
echo "  Appendix data (retry rates, path counts, postcheck failures):"
echo "    bench_*/combined_summary_appendix.csv"
echo
echo "  Per-trial raw data:"
echo "    bench_*/F_*/{wro,rpr,gsr}/times_seed*.csv"
echo "============================================================"
