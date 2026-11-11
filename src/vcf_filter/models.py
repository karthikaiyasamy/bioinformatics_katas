"""Data models for high-throughput VCF records and snpEff annotations."""

from typing import Dict, List, Optional
import sys
from ..compat import slotted_dataclass


@slotted_dataclass(slots=True, frozen=True)
class SnpEffAnnotation:
    """Structured snpEff functional effect annotation from the ANN tag."""
    allele: str
    annotation: str
    impact: str
    gene_name: str
    gene_id: str
    feature_type: str
    feature_id: str
    transcript_biotype: str
    rank: str
    hgvs_c: str
    hgvs_p: str
    cdna_pos: str
    cds_pos: str
    aa_pos: str
    distance: str
    warnings: str


@slotted_dataclass(slots=True, frozen=True)
class VcfRecord:
    """Memory-optimized representation of a parsed VCF variant record."""
    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: float
    filter_status: str
    info_raw: str
    info_map: Dict[str, str]
    annotations: List[SnpEffAnnotation]
    depth: int = 0
    allele_freq: float = 0.0
