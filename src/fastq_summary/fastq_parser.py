"""Memory-efficient streaming FASTQ reader."""

from pathlib import Path
from typing import Iterator, List
import gzip
from ..compat import slotted_dataclass


@slotted_dataclass(slots=True, frozen=True)
class FastqRead:
    """Represents a single sequencing read from a FASTQ file."""
    header: str
    sequence: str
    quality: str

    @property
    def length(self) -> int:
        return len(self.sequence)


def stream_fastq(fastq_path: Path) -> Iterator[FastqRead]:
    """Streams 4-line FASTQ records one by one in O(1) memory space."""
    fastq_path = Path(fastq_path)
    if not fastq_path.exists():
        raise FileNotFoundError(f"FASTQ file not found: {fastq_path}")

    is_gz = fastq_path.suffix == ".gz" or fastq_path.name.endswith(".fastq.gz") or fastq_path.name.endswith(".fq.gz")
    open_fn = gzip.open if is_gz else open
    mode = "rt" if is_gz else "r"

    with open_fn(fastq_path, mode=mode, encoding="utf-8", errors="replace") as stream:
        while True:
            header = stream.readline()
            if not header:
                break  # Clean EOF

            seq = stream.readline()
            plus = stream.readline()
            qual = stream.readline()

            if not (seq and plus and qual):
                # Truncated record at end of file
                break

            header = header.strip()
            seq = seq.strip().upper()
            qual = qual.strip()

            if header.startswith("@") and plus.startswith("+") and len(seq) == len(qual):
                yield FastqRead(header=header, sequence=seq, quality=qual)


def stream_fastq_chunks(fastq_path: Path, chunk_size: int = 5000) -> Iterator[List[FastqRead]]:
    """Yields batches of FastqRead objects for parallel worker consumption."""
    batch: List[FastqRead] = []
    for read in stream_fastq(fastq_path):
        batch.append(read)
        if len(batch) >= chunk_size:
            yield batch
            batch = []
    if batch:
        yield batch
