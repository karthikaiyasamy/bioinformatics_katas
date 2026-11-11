"""Parser for snpEff functional annotation (ANN) INFO tags."""

from typing import List
import sys
from .models import SnpEffAnnotation


def parse_snpeff_ann_tag(ann_raw: str) -> List[SnpEffAnnotation]:
    """Parse a raw snpEff ANN tag into a list of structured SnpEffAnnotation models.

    Format specification:
    Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type |
    Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length |
    CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO
    """
    if not ann_raw:
        return []

    results: List[SnpEffAnnotation] = []
    transcripts = ann_raw.split(",")

    for tr in transcripts:
        parts = [p.strip() for p in tr.split("|")]
        # Pad with empty strings if fewer than 16 fields present
        while len(parts) < 16:
            parts.append("")

        results.append(
            SnpEffAnnotation(
                allele=sys.intern(parts[0]) if parts[0] else "",
                annotation=sys.intern(parts[1]) if parts[1] else "",
                impact=sys.intern(parts[2]) if parts[2] else "",
                gene_name=sys.intern(parts[3]) if parts[3] else "",
                gene_id=parts[4],
                feature_type=parts[5],
                feature_id=parts[6],
                transcript_biotype=parts[7],
                rank=parts[8],
                hgvs_c=parts[9],
                hgvs_p=parts[10],
                cdna_pos=parts[11],
                cds_pos=parts[12],
                aa_pos=parts[13],
                distance=parts[14],
                warnings=parts[15],
            )
        )

    return results

