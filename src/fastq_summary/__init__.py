"""FASTQ Concurrent Summarizer Package."""

from .fastq_parser import FastqRead, stream_fastq, stream_fastq_chunks
from .concurrent_summarizer import (
    FastqChunkMetrics,
    FastqSummaryReport,
    analyze_chunk,
    summarize_fastq_concurrent,
)

__all__ = [
    "FastqRead",
    "FastqChunkMetrics",
    "FastqSummaryReport",
    "analyze_chunk",
    "stream_fastq",
    "stream_fastq_chunks",
    "summarize_fastq_concurrent",
]

