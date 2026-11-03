from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class PrimerProfile:
    id: str
    name: str
    marker: str
    forward_primer: str
    reverse_primer: str
    min_insert_bp: int
    max_insert_bp: int
    expected_insert_bp: int
    max_primer_error_rate: float
    taxonomic_scope: tuple[str, ...]
    doi: str


@dataclass
class FastqStats:
    reads_seen: int = 0
    bases_seen: int = 0
    mean_quality_sum: float = 0.0
    min_length: int | None = None
    max_length: int | None = None
    files: list[str] = field(default_factory=list)

    def add(self, length: int, mean_quality: float) -> None:
        self.reads_seen += 1
        self.bases_seen += length
        self.mean_quality_sum += mean_quality
        self.min_length = length if self.min_length is None else min(self.min_length, length)
        self.max_length = length if self.max_length is None else max(self.max_length, length)

    @property
    def mean_length(self) -> float:
        return self.bases_seen / self.reads_seen if self.reads_seen else 0.0

    @property
    def mean_quality(self) -> float:
        return self.mean_quality_sum / self.reads_seen if self.reads_seen else 0.0


@dataclass(frozen=True)
class PrimerDetection:
    profile: PrimerProfile
    sampled_reads: int
    paired_hits: int
    support_fraction: float
    confidence: str
    median_insert_bp: int | None


@dataclass(frozen=True)
class Cluster:
    id: str
    representative: str
    reads: int
    fraction: float
    unique_variants: int


@dataclass(frozen=True)
class ReferenceRecord:
    accession: str
    taxon: str
    rank: str
    marker: str
    sequence: str
    source: str = "unknown"
    audit_level: str = "unreviewed"


@dataclass(frozen=True)
class TaxonomicHit:
    cluster_id: str
    taxon: str
    rank: str
    accession: str
    identity: float
    query_coverage: float
    margin_to_second: float | None
    assignment: str
    confidence: str
    reference_source: str = "unknown"
    audit_level: str = "unreviewed"


@dataclass
class AnalysisConfig:
    sample_name: str
    inputs: list[Path]
    reference_path: Path
    output_dir: Path
    max_reads: int | None = None
    detection_sample_reads: int = 10_000
    min_mean_q: float = 10.0
    min_cluster_reads: int = 3
    min_cluster_fraction: float = 0.001
    cluster_identity: float = 0.94

