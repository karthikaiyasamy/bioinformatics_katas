"""Constant-memory streaming VCF parser and filter."""

from contextlib import contextmanager
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set
import gzip
import logging
import os
import sys
import tempfile

from .models import SnpEffAnnotation, VcfRecord
from .snpeff_parser import parse_snpeff_ann_tag

logger = logging.getLogger("VcfStream")


@contextmanager
def atomic_file_writer(target_path: Path, mode: str = "w", encoding: str = "utf-8"):
    """Atomically write file via a temporary file on the same mount."""
    target_path = Path(target_path)
    target_path.parent.mkdir(parents=True, exist_ok=True)

    tmp = tempfile.NamedTemporaryFile(
        mode=mode,
        dir=target_path.parent,
        prefix=f".{target_path.name}.tmp_",
        delete=False,
        encoding=encoding if "b" not in mode else None,
    )
    tmp_path = Path(tmp.name)

    try:
        yield tmp
        tmp.flush()
        os.fsync(tmp.fileno())
        tmp.close()
        os.replace(tmp_path, target_path)
    except Exception:
        tmp.close()
        if tmp_path.exists():
            tmp_path.unlink()
        raise


def parse_info_string(info_str: str) -> Dict[str, str]:
    """Parse key=value pairs from the VCF INFO column."""
    info_dict: Dict[str, str] = {}
    if not info_str or info_str == ".":
        return info_dict

    items = info_str.split(";")
    for item in items:
        if "=" in item:
            key, val = item.split("=", 1)
            info_dict[key.strip()] = val.strip()
        else:
            info_dict[item.strip()] = "True"
    return info_dict


def stream_vcf(vcf_path: Path) -> Iterator[VcfRecord]:
    """Stream VCF records one by one in O(1) auxiliary space."""
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    is_gz = vcf_path.suffix == ".gz" or vcf_path.name.endswith(".vcf.gz")
    open_fn = gzip.open if is_gz else open
    mode = "rt" if is_gz else "r"

    with open_fn(vcf_path, mode=mode, encoding="utf-8", errors="replace") as stream:
        for line_num, line in enumerate(stream, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 8:
                logger.warning("Line %d malformed: fewer than 8 required VCF columns", line_num)
                continue

            try:
                pos = int(fields[1])
                qual = float(fields[5]) if fields[5] != "." else 0.0
            except ValueError:
                logger.warning("Line %d malformed: non-numeric POS or QUAL", line_num)
                continue

            info_raw = fields[7]
            info_map = parse_info_string(info_raw)

            # Parse DP and AF if present
            depth = 0
            if "DP" in info_map:
                try:
                    depth = int(info_map["DP"])
                except ValueError:
                    pass

            allele_freq = 0.0
            if "AF" in info_map:
                try:
                    allele_freq = float(info_map["AF"])
                except ValueError:
                    pass

            # Parse snpEff annotations from ANN tag if present
            annotations: List[SnpEffAnnotation] = []
            if "ANN" in info_map:
                annotations = parse_snpeff_ann_tag(info_map["ANN"])

            yield VcfRecord(
                chrom=sys.intern(fields[0]),
                pos=pos,
                id=fields[2],
                ref=fields[3],
                alt=fields[4],
                qual=qual,
                filter_status=sys.intern(fields[6]),
                info_raw=info_raw,
                info_map=info_map,
                annotations=annotations,
                depth=depth,
                allele_freq=allele_freq,
            )


def filter_vcf_stream(
    records: Iterator[VcfRecord],
    min_qual: float = 30.0,
    min_depth: int = 10,
    pass_only: bool = True,
    target_impacts: Optional[Set[str]] = None,
    target_genes: Optional[Set[str]] = None,
) -> Iterator[VcfRecord]:
    """Yield records that satisfy all filtering criteria."""
    for record in records:
        if pass_only and record.filter_status != "PASS":
            continue
        if record.qual < min_qual:
            continue
        if record.depth < min_depth:
            continue

        if target_impacts or target_genes:
            # Check if any annotation matches impact and gene filters
            has_match = False
            for ann in record.annotations:
                impact_match = not target_impacts or (ann.impact in target_impacts)
                gene_match = not target_genes or (ann.gene_name in target_genes)
                if impact_match and gene_match:
                    has_match = True
                    break
            if not has_match:
                continue

        yield record


def export_filtered_report(
    input_vcf: Path,
    output_path: Path,
    min_qual: float = 30.0,
    min_depth: int = 10,
    pass_only: bool = True,
    target_impacts: Optional[Set[str]] = None,
    target_genes: Optional[Set[str]] = None,
) -> int:
    """Filter records and write a TSV report with atomic write guarantees."""
    exported_count = 0
    records = stream_vcf(input_vcf)
    filtered = filter_vcf_stream(
        records,
        min_qual=min_qual,
        min_depth=min_depth,
        pass_only=pass_only,
        target_impacts=target_impacts,
        target_genes=target_genes,
    )

    with atomic_file_writer(output_path) as out:
        out.write("CHROM\tPOS\tID\tREF\tALT\tQUAL\tDEPTH\tAF\tGENE\tEFFECT\tIMPACT\tHGVS_C\tHGVS_P\n")
        for rec in filtered:
            # Primary annotation
            ann = rec.annotations[0] if rec.annotations else None
            gene = ann.gene_name if ann else "."
            effect = ann.annotation if ann else "."
            impact = ann.impact if ann else "."
            hgvs_c = ann.hgvs_c if ann else "."
            hgvs_p = ann.hgvs_p if ann else "."

            out.write(
                f"{rec.chrom}\t{rec.pos}\t{rec.id}\t{rec.ref}\t{rec.alt}\t"
                f"{rec.qual:.1f}\t{rec.depth}\t{rec.allele_freq:.3f}\t"
                f"{gene}\t{effect}\t{impact}\t{hgvs_c}\t{hgvs_p}\n"
            )
            exported_count += 1

    return exported_count

