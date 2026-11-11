"""Genomic interval data model with 0-based half-open coordinate mechanics."""

from typing import Optional
import sys
from ..compat import slotted_dataclass


@slotted_dataclass(slots=True, frozen=True)
class BedInterval:
    """Represents a standard BED genomic interval [start, end) with 0-based half-open coordinates."""
    chrom: str
    start: int
    end: int
    name: Optional[str] = None
    score: Optional[float] = None
    strand: Optional[str] = None

    def __post_init__(self):
        if self.start < 0:
            raise ValueError(f"Start coordinate cannot be negative: {self.start}")
        if self.end < self.start:
            raise ValueError(f"End coordinate ({self.end}) cannot precede start ({self.start})")

    @property
    def length(self) -> int:
        """Length of interval in base pairs."""
        return self.end - self.start

    def overlaps_with(self, other: "BedInterval") -> bool:
        """Returns True if two intervals overlap on the same chromosome."""
        if self.chrom != other.chrom:
            return False
        return max(self.start, other.start) < min(self.end, other.end)

    def intersect(self, other: "BedInterval") -> Optional["BedInterval"]:
        """Compute the intersecting interval, or None if no overlap exists."""
        if not self.overlaps_with(other):
            return None

        o_start = max(self.start, other.start)
        o_end = min(self.end, other.end)
        name_combined = f"{self.name or '.'}&{other.name or '.'}"
        return BedInterval(chrom=self.chrom, start=o_start, end=o_end, name=name_combined)


@slotted_dataclass(slots=True, frozen=True)
class OverlapResult:
    """Detailed overlap result between a query interval (A) and database interval (B)."""
    interval_a: BedInterval
    interval_b: BedInterval
    overlap_start: int
    overlap_end: int
    overlap_bases: int
