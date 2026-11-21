#!/usr/bin/env python3
"""
Synthetic FASTQ Generator for Genomics Benchmarking.
Generates realistic Illumina NovaSeq-style 150bp paired-end reads with
decaying Phred quality profiles toward the 3-prime end.
"""

import sys
import random
from pathlib import Path

BASES = ("A", "C", "G", "T")
# Weighted toward realistic human genome GC content (~42-45%)
BASE_WEIGHTS = (0.28, 0.22, 0.22, 0.28)


def generate_benchmark_fastq(output_path: Path, num_reads: int = 1_000_000, read_length: int = 150):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Generating {num_reads:,} synthetic FASTQ reads ({read_length} bp each) into {output_path}...")

    # Pre-compute quality distributions simulating Illumina cycle degradation
    # Cycles 1-100: High Q37-Q40 (ASCII 'F' to 'I')
    # Cycles 101-150: Slight degradation Q28-Q34 (ASCII '=' to 'C')
    high_qual_chars = [chr(q + 33) for q in range(35, 41)]
    low_qual_chars = [chr(q + 33) for q in range(25, 35)]

    chunk_size = 50_000
    written = 0

    with open(output_path, "w", encoding="utf-8", buffering=1024 * 1024) as f:
        while written < num_reads:
            batch_count = min(chunk_size, num_reads - written)
            lines = []

            for i in range(batch_count):
                read_num = written + i + 1
                header = f"@NOVASEQ:7:1101:{read_num}:1000 1:N:0:TCGA"

                # Random sequence with realistic GC
                seq = "".join(random.choices(BASES, weights=BASE_WEIGHTS, k=read_length))

                # Quality string with realistic decay
                q_start = "".join(random.choices(high_qual_chars, k=100))
                q_end = "".join(random.choices(low_qual_chars, k=read_length - 100))
                qual = q_start + q_end

                lines.append(f"{header}\n{seq}\n+\n{qual}\n")

            f.writelines(lines)
            written += batch_count
            percent = (written / num_reads) * 100
            print(f"Progress: {written:,} / {num_reads:,} reads ({percent:.1f}%)")

    file_size_mb = output_path.stat().st_size / (1024 * 1024)
    print(f"Done. File created: {output_path} ({file_size_mb:.2f} MB)")


if __name__ == "__main__":
    count = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
    dest = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("data/benchmark_1m.fastq")
    generate_benchmark_fastq(dest, num_reads=count)
