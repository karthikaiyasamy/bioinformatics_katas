#!/usr/bin/env bash
# ==============================================================================
# Bioinformatics High-Throughput Benchmarking Suite
# Reproducible comparison of Python vs Rust on FASTQ quality & GC analytics
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

FASTQ_FILE="data/benchmark_500k.fastq"

echo "======================================================================"
echo " BIOINFORMATICS HIGH-THROUGHPUT BENCHMARK REPRODUCIBILITY SUITE"
echo " Comparing Python 3 vs Rust (Rayon) on 75 Million Genomic Bases"
echo "======================================================================"
echo ""

# 1. Verify / Generate Benchmark Dataset
if [[ ! -f "$FASTQ_FILE" ]]; then
    echo ">> [1/5] Generating 500,000-read benchmark dataset..."
    /usr/bin/python3 generate_benchmark_fastq.py
else
    echo ">> [1/5] Benchmark dataset already present: $FASTQ_FILE"
fi
echo ""

# 2. Run Test Suite
echo ">> [2/5] Running algorithm verification test suite..."
/usr/bin/python3 run_tests.py
echo ""

# 3. Build Rust Analyzer
echo ">> [3/5] Building Rust release binary with LTO..."
cd fastq_summary_rs
cargo build --release
cd ..
echo ""

# 4. Run Benchmarks
echo ">> [4/5] Executing Python Sequential Benchmark..."
/usr/bin/time -l /usr/bin/python3 run_python_bench.py "$FASTQ_FILE"
echo ""

echo ">> [4/5] Executing Python Concurrent Benchmark (4 Workers)..."
/usr/bin/time -l /usr/bin/python3 run_python_bench.py "$FASTQ_FILE" --concurrent
echo ""

echo ">> [5/5] Executing Rust Sequential Benchmark (Zero-Allocation)..."
/usr/bin/time -l ./fastq_summary_rs/target/release/fastq_summary_rs "$FASTQ_FILE"
echo ""

echo ">> [5/5] Executing Rust Parallel Benchmark (Rayon)..."
/usr/bin/time -l ./fastq_summary_rs/target/release/fastq_summary_rs "$FASTQ_FILE" --parallel
echo ""

echo "======================================================================"
echo " Benchmarks completed successfully. See COMPARATIVE_STUDY_PYTHON_VS_RUST.md"
echo "======================================================================"
