"""Tests for Kata 1: Constant-Memory VCF Stream Parser & snpEff Extractor."""

from pathlib import Path
try:
    import pytest
except ImportError:
    pytest = None

from src.vcf_filter.models import VcfRecord
from src.vcf_filter.snpeff_parser import parse_snpeff_ann_tag
from src.vcf_filter.vcf_stream import (
    export_filtered_report,
    filter_vcf_stream,
    parse_info_string,
    stream_vcf,
)


def test_parse_info_string():
    raw_info = "DP=150;AF=0.45;SOMATIC;ANN=A|missense|MODERATE|TP53"
    info_map = parse_info_string(raw_info)
    assert info_map["DP"] == "150"
    assert info_map["AF"] == "0.45"
    assert info_map["SOMATIC"] == "True"
    assert info_map["ANN"].startswith("A|missense")


def test_parse_snpeff_ann_tag():
    raw_ann = "A|stop_gained|HIGH|TP53|ENSG00000141510|transcript|ENST00000269305|protein_coding|8/11|c.659G>A|p.Arg220Gln|850/2576|659/1182|220/393||"
    annotations = parse_snpeff_ann_tag(raw_ann)
    assert len(annotations) == 1
    ann = annotations[0]
    assert ann.allele == "A"
    assert ann.annotation == "stop_gained"
    assert ann.impact == "HIGH"
    assert ann.gene_name == "TP53"
    assert ann.hgvs_c == "c.659G>A"
    assert ann.hgvs_p == "p.Arg220Gln"


def test_stream_vcf_plain_and_gzipped(sample_vcf_path: Path, sample_vcf_gz_path: Path):
    plain_records = list(stream_vcf(sample_vcf_path))
    gz_records = list(stream_vcf(sample_vcf_gz_path))

    assert len(plain_records) == 8
    assert len(gz_records) == 8

    # Verify field extraction on first record (EGFR)
    r0 = plain_records[0]
    assert r0.chrom == "chr7"
    assert r0.pos == 55242465
    assert r0.ref == "G"
    assert r0.alt == "A"
    assert r0.qual == 99.0
    assert r0.filter_status == "PASS"
    assert r0.depth == 142
    assert r0.allele_freq == 0.42
    assert len(r0.annotations) == 1
    assert r0.annotations[0].gene_name == "EGFR"


def test_filter_vcf_stream(sample_vcf_path: Path):
    records = stream_vcf(sample_vcf_path)

    # Filter for PASS variants with QUAL >= 30, DP >= 50, and HIGH or MODERATE impact
    filtered = list(
        filter_vcf_stream(
            records,
            min_qual=30.0,
            min_depth=50,
            pass_only=True,
            target_impacts={"HIGH", "MODERATE"},
        )
    )

    # In our sample data:
    # chr7 EGFR (PASS, DP 142, QUAL 99, MODERATE) -> Keep
    # chr17 TP53 (PASS, DP 210, QUAL 150, HIGH) -> Keep
    # chr12 KRAS (PASS, DP 98, QUAL 120, MODERATE) -> Keep
    # chr7 BRAF (PASS, DP 320, QUAL 200, MODERATE) -> Keep
    # chr3 PIK3CA (PASS, DP 75, QUAL 85, MODERATE) -> Keep
    # chr17 BRCA1 (LowQual, DP 12, QUAL 18.5) -> Filtered out
    # chr13 BRCA2 (LowQual, DP 8, QUAL 12) -> Filtered out
    # chr17 TP53 (PASS, DP 115, QUAL 95, MODERATE) -> Keep
    assert len(filtered) == 6

    # Test filtering specifically for TP53
    records_tp53 = stream_vcf(sample_vcf_path)
    tp53_only = list(
        filter_vcf_stream(
            records_tp53,
            target_genes={"TP53"},
        )
    )
    assert len(tp53_only) == 2
    for r in tp53_only:
        assert r.annotations[0].gene_name == "TP53"


def test_export_filtered_report(sample_vcf_gz_path: Path, tmp_path: Path):
    out_tsv = tmp_path / "filtered_clinical_report.tsv"
    count = export_filtered_report(
        input_vcf=sample_vcf_gz_path,
        output_path=out_tsv,
        min_qual=50.0,
        min_depth=50,
        pass_only=True,
        target_impacts={"HIGH", "MODERATE"},
    )
    assert count == 6
    assert out_tsv.exists()

    content = out_tsv.read_text(encoding="utf-8")
    lines = content.strip().split("\n")
    assert lines[0].startswith("CHROM\tPOS\tID\tREF\tALT")
    assert len(lines) == 7  # Header + 6 records
