#!/usr/bin/env python3
"""
Python FASTQ Benchmarking Runner.
Runs sequential and concurrent quality & GC summarization on FASTQ files.
"""

import sys
import time
from pathlib import Path

# Add project root to Python path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.fastq_summary.fastq_parser import stream_fastq
from src.fastq_summary.concurrent_summarizer import summarize_fastq_concurrent


def run_sequential(fastq_path: Path):
    print("=== Python FASTQ Quality & GC Analyzer ===")
    print(f"File: {fastq_path}")
    print("Mode: Sequential (Pure Python Stream)")

    start = time.perf_counter()

    total_reads = 0
    total_bases = 0
    gc_bases = 0
    q20_bases = 0
    q30_bases = 0
    qual_sum = 0

    for read in stream_fastq(fastq_path):
        seq = read.sequence
        qual = read.quality
        r_len = len(seq)

        total_reads += 1
        total_bases += r_len

        for i in range(r_len):
            base = seq[i]
            if base in ("G", "C"):
                gc_bases += 1

            q_score = ord(qual[i]) - 33
            qual_sum += q_score

            if q_score >= 20:
                q20_bases += 1
            if q_score >= 30:
                q30_bases += 1

    duration = time.perf_counter() - start

    mean_qual = (qual_sum / total_bases) if total_bases > 0 else 0.0
    gc_pct = (gc_bases / total_bases * 100.0) if total_bases > 0 else 0.0
    q20_pct = (q20_bases / total_bases * 100.0) if total_bases > 0 else 0.0
    q30_pct = (q30_bases / total_bases * 100.0) if total_bases > 0 else 0.0

    print("\n--- Results ---")
    print(f"Total Reads:        {total_reads:>12,}")
    print(f"Total Bases:        {total_bases:>12,}")
    print(f"Mean Quality (Q):   {mean_qual:>12.2f}")
    print(f"GC Percentage:      {gc_pct:>11.2f}%")
    print(f"Q20 Percentage:     {q20_pct:>11.2f}%")
    print(f"Q30 Percentage:     {q30_pct:>11.2f}%")
    print(f"Elapsed Time:       {duration:>12.3f}s ({duration * 1000:.1f}ms)")

    throughput = (total_bases / 1_000_000.0) / duration if duration > 0 else 0.0
    print(f"Throughput:         {throughput:>11.2f} Megabases/sec")


def run_concurrent(fastq_path: Path):
    print("=== Python FASTQ Quality & GC Analyzer ===")
    print(f"File: {fastq_path}")
    print("Mode: Concurrent (ProcessPoolExecutor - 4 Workers)")

    start = time.perf_counter()
    report = summarize_fastq_concurrent(fastq_path, chunk_size=10_000, max_workers=4)
    duration = time.perf_counter() - start

    print("\n--- Results ---")
    print(f"Total Reads:        {report.total_reads:>12,}")
    print(f"Total Bases:        {report.total_bases:>12,}")
    print(f"Mean Quality (Q):   {report.mean_quality:>12.2f}")
    print(f"GC Percentage:      {report.gc_percentage:>11.2f}%")
    print(f"Q20 Percentage:     {report.q20_percentage:>11.2f}%")
    print(f"Q30 Percentage:     {report.q30_percentage:>11.2f}%")
    print(f"Elapsed Time:       {duration:>12.3f}s ({duration * 1000:.1f}ms)")

    throughput = (report.total_bases / 1_000_000.0) / duration if duration > 0 else 0.0
    print(f"Throughput:         {throughput:>11.2f} Megabases/sec")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python run_python_bench.py <path_to_fastq> [--concurrent]")
        sys.exit(1)

    path = Path(sys.argv[1])
    if "--concurrent" in sys.argv:
        run_concurrent(path)
    else:
        run_sequential(path)
