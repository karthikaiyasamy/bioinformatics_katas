"""Concurrent quality score and GC-content summarizer for high-throughput FASTQ data."""

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import field
from pathlib import Path
from typing import Dict, List, Optional
import os

from ..compat import slotted_dataclass
from .fastq_parser import FastqRead, stream_fastq_chunks


@slotted_dataclass(slots=True)
class FastqChunkMetrics:
    """Intermediate metrics accumulated from a single chunk of reads."""
    total_reads: int = 0
    total_bases: int = 0
    gc_bases: int = 0
    q20_bases: int = 0
    q30_bases: int = 0
    total_quality_sum: int = 0
    cycle_quality_sums: List[int] = field(default_factory=list)
    cycle_counts: List[int] = field(default_factory=list)


@slotted_dataclass(slots=True, frozen=True)
class FastqSummaryReport:
    """Final aggregated quality report for a FASTQ dataset."""
    total_reads: int
    total_bases: int
    mean_quality: float
    gc_percentage: float
    q20_percentage: float
    q30_percentage: float
    per_cycle_mean_qualities: List[float]


def analyze_chunk(reads: List[FastqRead]) -> FastqChunkMetrics:
    """Worker function: Computes Phred metrics for a list of reads without GIL blocking."""
    metrics = FastqChunkMetrics()

    for read in reads:
        seq = read.sequence
        qual = read.quality
        r_len = len(seq)

        metrics.total_reads += 1
        metrics.total_bases += r_len

        # Grow cycle arrays if this read is longer than previous ones
        if len(metrics.cycle_quality_sums) < r_len:
            extend_len = r_len - len(metrics.cycle_quality_sums)
            metrics.cycle_quality_sums.extend([0] * extend_len)
            metrics.cycle_counts.extend([0] * extend_len)

        for idx, (base, q_char) in enumerate(zip(seq, qual)):
            if base in ("G", "C"):
                metrics.gc_bases += 1

            # Phred+33 quality decoding
            q_score = ord(q_char) - 33
            metrics.total_quality_sum += q_score

            if q_score >= 20:
                metrics.q20_bases += 1
            if q_score >= 30:
                metrics.q30_bases += 1

            metrics.cycle_quality_sums[idx] += q_score
            metrics.cycle_counts[idx] += 1

    return metrics


def summarize_fastq_concurrent(
    fastq_path: Path,
    chunk_size: int = 2500,
    max_workers: Optional[int] = None,
) -> FastqSummaryReport:
    """Summarizes FASTQ statistics across multiple CPU processes for maximum throughput."""
    fastq_path = Path(fastq_path)
    if not fastq_path.exists():
        raise FileNotFoundError(f"File not found: {fastq_path}")

    # Default to available CPUs (capped at 8 for standard jobs)
    workers = max_workers or min(os.cpu_count() or 2, 8)

    total_reads = 0
    total_bases = 0
    total_gc = 0
    total_q20 = 0
    total_q30 = 0
    total_qual_sum = 0
    master_cycle_sums: List[int] = []
    master_cycle_counts: List[int] = []

    chunks = list(stream_fastq_chunks(fastq_path, chunk_size=chunk_size))
    if not chunks:
        return FastqSummaryReport(
            total_reads=0,
            total_bases=0,
            mean_quality=0.0,
            gc_percentage=0.0,
            q20_percentage=0.0,
            q30_percentage=0.0,
            per_cycle_mean_qualities=[],
        )

    # Process chunks in parallel using ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_to_chunk = {executor.submit(analyze_chunk, chunk): chunk for chunk in chunks}

        for future in as_completed(future_to_chunk):
            res = future.result()
            total_reads += res.total_reads
            total_bases += res.total_bases
            total_gc += res.gc_bases
            total_q20 += res.q20_bases
            total_q30 += res.q30_bases
            total_qual_sum += res.total_quality_sum

            # Merge per-cycle arrays
            if len(master_cycle_sums) < len(res.cycle_quality_sums):
                ext = len(res.cycle_quality_sums) - len(master_cycle_sums)
                master_cycle_sums.extend([0] * ext)
                master_cycle_counts.extend([0] * ext)

            for i in range(len(res.cycle_quality_sums)):
                master_cycle_sums[i] += res.cycle_quality_sums[i]
                master_cycle_counts[i] += res.cycle_counts[i]

    # Calculate final aggregate metrics
    mean_qual = (total_qual_sum / total_bases) if total_bases > 0 else 0.0
    gc_pct = (total_gc / total_bases * 100.0) if total_bases > 0 else 0.0
    q20_pct = (total_q20 / total_bases * 100.0) if total_bases > 0 else 0.0
    q30_pct = (total_q30 / total_bases * 100.0) if total_bases > 0 else 0.0

    cycle_means = [
        round(master_cycle_sums[i] / master_cycle_counts[i], 2)
        if master_cycle_counts[i] > 0
        else 0.0
        for i in range(len(master_cycle_sums))
    ]

    return FastqSummaryReport(
        total_reads=total_reads,
        total_bases=total_bases,
        mean_quality=round(mean_qual, 2),
        gc_percentage=round(gc_pct, 2),
        q20_percentage=round(q20_pct, 2),
        q30_percentage=round(q30_pct, 2),
        per_cycle_mean_qualities=cycle_means,
    )
