"""Pytest shared fixtures and data paths for bioinformatics katas."""

from pathlib import Path
import pytest


@pytest.fixture(scope="session")
def kata_root() -> Path:
    """Returns absolute path to the bioinformatics_katas root directory."""
    return Path(__file__).resolve().parent.parent


@pytest.fixture(scope="session")
def data_dir(kata_root: Path) -> Path:
    """Returns absolute path to sample test datasets."""
    return kata_root / "data"


@pytest.fixture(scope="session")
def sample_vcf_path(data_dir: Path) -> Path:
    return data_dir / "sample.vcf"


@pytest.fixture(scope="session")
def sample_vcf_gz_path(data_dir: Path) -> Path:
    return data_dir / "sample.vcf.gz"


@pytest.fixture(scope="session")
def regions_a_bed_path(data_dir: Path) -> Path:
    return data_dir / "regions_a.bed"


@pytest.fixture(scope="session")
def regions_b_bed_path(data_dir: Path) -> Path:
    return data_dir / "regions_b.bed"


@pytest.fixture(scope="session")
def sample_fastq_path(data_dir: Path) -> Path:
    return data_dir / "sample_reads.fastq"


@pytest.fixture(scope="session")
def sample_fastq_gz_path(data_dir: Path) -> Path:
    return data_dir / "sample_reads.fastq.gz"

