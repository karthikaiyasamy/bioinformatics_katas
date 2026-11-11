#!/usr/bin/env python3
"""Zero-dependency test runner for bioinformatics katas.

Executes all test functions across the suite without requiring third-party libraries.
If pytest is installed, 'pytest' can also be run directly.
"""

from pathlib import Path
import os
import sys
import tempfile
import time

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

DATA_DIR = PROJECT_ROOT / "data"
SAMPLE_VCF = DATA_DIR / "sample.vcf"
SAMPLE_VCF_GZ = DATA_DIR / "sample.vcf.gz"
REGIONS_A = DATA_DIR / "regions_a.bed"
REGIONS_B = DATA_DIR / "regions_b.bed"
SAMPLE_FASTQ = DATA_DIR / "sample_reads.fastq"
SAMPLE_FASTQ_GZ = DATA_DIR / "sample_reads.fastq.gz"


def run_all():
    print("=" * 70)
    print("BIOINFORMATICS HIGH-THROUGHPUT KATAS: TEST RUNNER")
    print("=" * 70)

    start_time = time.time()
    total_passed = 0

    # Kata 1: VCF Filter & snpEff Extractor
    print("\n[Kata 1] VCF Streaming Filter & snpEff Extractor:")
    import tests.test_vcf_filter as tvcf

    tvcf.test_parse_info_string()
    print("  ✓ test_parse_info_string PASSED")
    total_passed += 1

    tvcf.test_parse_snpeff_ann_tag()
    print("  ✓ test_parse_snpeff_ann_tag PASSED")
    total_passed += 1

    tvcf.test_stream_vcf_plain_and_gzipped(SAMPLE_VCF, SAMPLE_VCF_GZ)
    print("  ✓ test_stream_vcf_plain_and_gzipped PASSED")
    total_passed += 1

    tvcf.test_filter_vcf_stream(SAMPLE_VCF)
    print("  ✓ test_filter_vcf_stream PASSED")
    total_passed += 1

    with tempfile.TemporaryDirectory() as tmp:
        tvcf.test_export_filtered_report(SAMPLE_VCF_GZ, Path(tmp))
    print("  ✓ test_export_filtered_report PASSED")
    total_passed += 1

    # Kata 2: Interval Overlap Engine
    print("\n[Kata 2] Zero-Dependency Genomic Interval Overlap Engine:")
    import tests.test_bed_intersection as tbed

    tbed.test_bed_interval_mechanics()
    print("  ✓ test_bed_interval_mechanics PASSED")
    total_passed += 1

    tbed.test_parse_bed_line()
    print("  ✓ test_parse_bed_line PASSED")
    total_passed += 1

    tbed.test_sweep_line_intersect(REGIONS_A, REGIONS_B)
    print("  ✓ test_sweep_line_intersect PASSED")
    total_passed += 1

    with tempfile.TemporaryDirectory() as tmp:
        tbed.test_export_bed_intersection(REGIONS_A, REGIONS_B, Path(tmp))
    print("  ✓ test_export_bed_intersection PASSED")
    total_passed += 1

    # Kata 3: Multi-Threaded FASTQ Summarizer
    print("\n[Kata 3] Multi-Threaded FASTQ Quality Score Summarizer:")
    import tests.test_fastq_summary as tfq

    tfq.test_fastq_read_and_phred_decoding()
    print("  ✓ test_fastq_read_and_phred_decoding PASSED")
    total_passed += 1

    tfq.test_stream_fastq(SAMPLE_FASTQ, SAMPLE_FASTQ_GZ)
    print("  ✓ test_stream_fastq PASSED")
    total_passed += 1

    tfq.test_stream_fastq_chunks(SAMPLE_FASTQ)
    print("  ✓ test_stream_fastq_chunks PASSED")
    total_passed += 1

    tfq.test_summarize_fastq_concurrent(SAMPLE_FASTQ_GZ)
    print("  ✓ test_summarize_fastq_concurrent PASSED")
    total_passed += 1

    elapsed = time.time() - start_time
    print("\n" + "=" * 70)
    print(f"SUCCESS: ALL {total_passed} TESTS PASSED in {elapsed:.3f} seconds!")
    print("=" * 70)


if __name__ == "__main__":
    run_all()

