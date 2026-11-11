"""Tests for Kata 2: Zero-Dependency Interval Overlap Engine."""

from pathlib import Path
try:
    import pytest
except ImportError:
    pytest = None

from src.bed_intersection.interval import BedInterval
from src.bed_intersection.intersect_engine import (
    export_bed_intersection,
    load_bed_intervals,
    parse_bed_line,
    sweep_line_intersect,
)


def test_bed_interval_mechanics():
    # 0-based half-open coordinates: [100, 200) length is 100
    iv1 = BedInterval("chr1", 100, 200, name="A")
    assert iv1.length == 100

    # Overlapping interval: [150, 250)
    iv2 = BedInterval("chr1", 150, 250, name="B")
    assert iv1.overlaps_with(iv2)

    inter = iv1.intersect(iv2)
    assert inter is not None
    assert inter.start == 150
    assert inter.end == 200
    assert inter.length == 50

    # Abutting intervals: [100, 200) and [200, 300) do NOT overlap in half-open logic
    iv3 = BedInterval("chr1", 200, 300, name="C")
    assert not iv1.overlaps_with(iv3)
    assert iv1.intersect(iv3) is None

    # Different chromosomes do NOT overlap
    iv_diff = BedInterval("chr2", 150, 250)
    assert not iv1.overlaps_with(iv_diff)


def test_parse_bed_line():
    line = "chr7\t55240000\t55250000\tEGFR_Exon19\t1000\t+\n"
    iv = parse_bed_line(line)
    assert iv is not None
    assert iv.chrom == "chr7"
    assert iv.start == 55240000
    assert iv.end == 55250000
    assert iv.name == "EGFR_Exon19"
    assert iv.score == 1000.0
    assert iv.strand == "+"

    # Comments and track lines return None
    assert parse_bed_line("# Chromosome start end") is None
    assert parse_bed_line("track name=probes") is None


def test_sweep_line_intersect(regions_a_bed_path: Path, regions_b_bed_path: Path):
    intervals_a = load_bed_intervals(regions_a_bed_path)
    intervals_b = load_bed_intervals(regions_b_bed_path)

    assert len(intervals_a) == 8
    assert len(intervals_b) == 8

    overlaps = sweep_line_intersect(intervals_a, intervals_b, min_overlap=1)

    # In our sample data:
    # 1. chr7 EGFR_Exon19 [55240000, 55250000) & Probe_EGFR_Hotspot [55242000, 55245000) -> 3000 bp
    # 2. chr7 BRAF_Exon15 [140450000, 140460000) & Probe_BRAF_V600E [140452000, 140455000) -> 3000 bp
    # 3. chr12 KRAS_Exon2 [25390000, 25400000) & Probe_KRAS_Spanning [25395000, 25405000) -> 5000 bp (25395000 to 25400000)
    # 4. chr17 TP53_Core [7570000, 7580000) & Probe_TP53_DNA_Binding [7575000, 7578000) -> 3000 bp
    # 5. chr17 BRCA1_RegionA [41240000, 41250000) & Probe_BRCA1_Splice [41243000, 41246000) -> 3000 bp
    # 6. chr3 PIK3CA_Exon9 [178930000, 178940000) & Probe_PIK3CA_Helical [178935000, 178938000) -> 3000 bp
    # Non-overlapping:
    # chr13 BRCA2_Exon11 [32900000, 32910000) vs Probe_BRCA2 [32950000, 32960000) -> disjoint!
    # chr2 ALK vs chr1 NRAS -> disjoint chroms!
    assert len(overlaps) == 6

    # Verify KRAS spanning overlap coordinates
    kras_ov = [ov for ov in overlaps if ov.interval_a.chrom == "chr12"][0]
    assert kras_ov.overlap_start == 25395000
    assert kras_ov.overlap_end == 25400000
    assert kras_ov.overlap_bases == 5000


def test_export_bed_intersection(regions_a_bed_path: Path, regions_b_bed_path: Path, tmp_path: Path):
    out_bed = tmp_path / "intersected_targets.bed"
    count = export_bed_intersection(regions_a_bed_path, regions_b_bed_path, out_bed)
    assert count == 6
    assert out_bed.exists()

    content = out_bed.read_text(encoding="utf-8")
    lines = content.strip().split("\n")
    assert lines[0].startswith("#CHROM\tSTART\tEND")
    assert len(lines) == 7  # Header + 6 overlaps
