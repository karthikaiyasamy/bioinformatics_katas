"""Tests for Kata 3: Multi-Threaded FASTQ Quality Score Summarizer."""

from pathlib import Path
try:
    import pytest
except ImportError:
    pytest = None

from src.fastq_summary.fastq_parser import FastqRead, stream_fastq, stream_fastq_chunks
from src.fastq_summary.concurrent_summarizer import (
    analyze_chunk,
    summarize_fastq_concurrent,
)


def test_fastq_read_and_phred_decoding():
    # 'I' is ASCII 73 -> Phred quality = 73 - 33 = 40
    # '#' is ASCII 35 -> Phred quality = 35 - 33 = 2
    read = FastqRead(
        header="@SEQ_ID_01",
        sequence="ACGT",
        quality="I#I#",
    )
    assert read.length == 4

    metrics = analyze_chunk([read])
    assert metrics.total_reads == 1
    assert metrics.total_bases == 4
    assert metrics.gc_bases == 2  # C and G
    assert metrics.q20_bases == 2  # Two 'I' characters (Q40 >= 20)
    assert metrics.q30_bases == 2  # Two 'I' characters (Q40 >= 30)
    assert metrics.total_quality_sum == (40 + 2 + 40 + 2)  # 84
    assert metrics.cycle_quality_sums == [40, 2, 40, 2]


def test_stream_fastq(sample_fastq_path: Path, sample_fastq_gz_path: Path):
    plain_reads = list(stream_fastq(sample_fastq_path))
    gz_reads = list(stream_fastq(sample_fastq_gz_path))

    assert len(plain_reads) == 1000
    assert len(gz_reads) == 1000

    r0 = plain_reads[0]
    assert r0.header.startswith("@A00123")
    assert r0.length == 75
    assert len(r0.quality) == 75


def test_stream_fastq_chunks(sample_fastq_path: Path):
    chunks = list(stream_fastq_chunks(sample_fastq_path, chunk_size=300))
    # 1000 reads in chunks of 300 -> 300, 300, 300, 100 (4 chunks)
    assert len(chunks) == 4
    assert len(chunks[0]) == 300
    assert len(chunks[3]) == 100


def test_summarize_fastq_concurrent(sample_fastq_gz_path: Path):
    report = summarize_fastq_concurrent(
        fastq_path=sample_fastq_gz_path,
        chunk_size=250,
        max_workers=2,
    )

    assert report.total_reads == 1000
    assert report.total_bases == 75000  # 1000 reads * 75 bp
    assert 40.0 <= report.gc_percentage <= 60.0  # Random 4 bases is ~50%
    assert report.mean_quality > 25.0
    assert report.q20_percentage > 80.0
    assert len(report.per_cycle_mean_qualities) == 75

    # Confirm 5' quality is higher than 3' quality as simulated
    cycle_5_prime = report.per_cycle_mean_qualities[0]
    cycle_3_prime = report.per_cycle_mean_qualities[-1]
    assert cycle_5_prime > cycle_3_prime
