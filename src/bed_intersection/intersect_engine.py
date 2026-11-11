"""High-performance pure Python genomic interval intersection engine."""

from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Optional
import sys
import tempfile
import os

from .interval import BedInterval, OverlapResult


def parse_bed_line(line: str) -> Optional[BedInterval]:
    """Parse a single BED line into a BedInterval, ignoring headers and track definitions."""
    line = line.strip()
    if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
        return None

    fields = line.split("\t")
    if len(fields) < 3:
        # Fall back to whitespace splitting if not strictly tab-delimited
        fields = line.split()
        if len(fields) < 3:
            return None

    chrom = sys.intern(fields[0])
    start = int(fields[1])
    end = int(fields[2])
    name = fields[3] if len(fields) > 3 else None
    score = float(fields[4]) if len(fields) > 4 and fields[4] != "." else None
    strand = fields[5] if len(fields) > 5 and fields[5] in ("+", "-", ".") else None

    return BedInterval(chrom=chrom, start=start, end=end, name=name, score=score, strand=strand)


def load_bed_intervals(bed_path: Path) -> List[BedInterval]:
    """Load all valid intervals from a BED file."""
    bed_path = Path(bed_path)
    if not bed_path.exists():
        raise FileNotFoundError(f"BED file not found: {bed_path}")

    intervals: List[BedInterval] = []
    with open(bed_path, "r", encoding="utf-8") as f:
        for line in f:
            interval = parse_bed_line(line)
            if interval:
                intervals.append(interval)
    return intervals


def sweep_line_intersect(
    intervals_a: List[BedInterval],
    intervals_b: List[BedInterval],
    min_overlap: int = 1,
) -> List[OverlapResult]:
    """Find all intersecting intervals using a sweep-line algorithm with chromosome grouping.

    Complexity: O(N log N + M log M + K) where K is number of overlaps.
    Does not require any third-party C-extensions or libraries.
    """
    # Group intervals by chromosome
    chrom_map_a: Dict[str, List[BedInterval]] = defaultdict(list)
    chrom_map_b: Dict[str, List[BedInterval]] = defaultdict(list)

    for a in intervals_a:
        chrom_map_a[a.chrom].append(a)
    for b in intervals_b:
        chrom_map_b[b.chrom].append(b)

    results: List[OverlapResult] = []

    # Process intersecting chromosomes
    common_chroms = set(chrom_map_a.keys()) & set(chrom_map_b.keys())

    for chrom in sorted(common_chroms):
        # Sort by start coordinate ascending
        sorted_a = sorted(chrom_map_a[chrom], key=lambda x: (x.start, x.end))
        sorted_b = sorted(chrom_map_b[chrom], key=lambda x: (x.start, x.end))

        b_idx_start = 0
        len_b = len(sorted_b)

        for a in sorted_a:
            # Advance b_idx_start while b ends before or at a.start
            while b_idx_start < len_b and sorted_b[b_idx_start].end <= a.start:
                b_idx_start += 1

            # Check all intervals in B that could overlap with A
            curr_b = b_idx_start
            while curr_b < len_b:
                b = sorted_b[curr_b]
                # If b starts at or after a ends, no subsequent B can overlap this A
                if b.start >= a.end:
                    break

                # Overlap check
                o_start = max(a.start, b.start)
                o_end = min(a.end, b.end)
                overlap_len = o_end - o_start

                if overlap_len >= min_overlap:
                    results.append(
                        OverlapResult(
                            interval_a=a,
                            interval_b=b,
                            overlap_start=o_start,
                            overlap_end=o_end,
                            overlap_bases=overlap_len,
                        )
                    )
                curr_b += 1

    return results


def export_bed_intersection(
    bed_a_path: Path,
    bed_b_path: Path,
    output_path: Path,
    min_overlap: int = 1,
) -> int:
    """Intersects two BED files and writes the intersecting intervals with atomic rename."""
    intervals_a = load_bed_intervals(bed_a_path)
    intervals_b = load_bed_intervals(bed_b_path)

    overlaps = sweep_line_intersect(intervals_a, intervals_b, min_overlap=min_overlap)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    tmp = tempfile.NamedTemporaryFile(
        mode="w",
        dir=output_path.parent,
        prefix=f".{output_path.name}.tmp_",
        delete=False,
        encoding="utf-8",
    )
    tmp_path = Path(tmp.name)

    try:
        tmp.write("#CHROM\tSTART\tEND\tOVERLAP_NAME\tOVERLAP_BP\tA_NAME\tB_NAME\n")
        for ov in overlaps:
            name_a = ov.interval_a.name or "."
            name_b = ov.interval_b.name or "."
            tmp.write(
                f"{ov.interval_a.chrom}\t{ov.overlap_start}\t{ov.overlap_end}\t"
                f"{name_a}&{name_b}\t{ov.overlap_bases}\t{name_a}\t{name_b}\n"
            )
        tmp.flush()
        os.fsync(tmp.fileno())
        tmp.close()
        os.replace(tmp_path, output_path)
    except Exception:
        tmp.close()
        if tmp_path.exists():
            tmp_path.unlink()
        raise

    return len(overlaps)

