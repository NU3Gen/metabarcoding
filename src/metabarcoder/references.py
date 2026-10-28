from __future__ import annotations

from pathlib import Path

from .models import ReferenceRecord, TaxonomicHit
from .sequence import levenshtein


def load_reference_fasta(path: Path) -> list[ReferenceRecord]:
    records: list[ReferenceRecord] = []
    header: str | None = None
    sequence: list[str] = []

    def emit() -> None:
        if header is None:
            return
        fields = dict(part.split("=", 1) for part in header.split("|") if "=" in part)
        records.append(ReferenceRecord(
            accession=fields.get("accession", header.split("|")[0]),
            taxon=fields.get("taxon", "unclassified"), rank=fields.get("rank", "species"),
            marker=fields.get("marker", "unknown"), sequence="".join(sequence).upper(),
        ))

    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                emit()
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
        emit()
    if not records:
        raise ValueError(f"Referenčná databáza je prázdna: {path}")
    return records


def classify(sequence: str, cluster_id: str, references: list[ReferenceRecord],
             min_identity: float = 0.90, species_identity: float = 0.97,
             species_margin: float = 0.01) -> TaxonomicHit:
    scored: list[tuple[float, ReferenceRecord]] = []
    for record in references:
        distance = levenshtein(sequence, record.sequence)
        identity = 1.0 - distance / max(len(sequence), len(record.sequence))
        scored.append((identity, record))
    scored.sort(key=lambda item: item[0], reverse=True)
    best_identity, best = scored[0]
    second_identity = scored[1][0] if len(scored) > 1 else None
    margin = best_identity - second_identity if second_identity is not None else None
    coverage = min(len(sequence), len(best.sequence)) / max(1, len(sequence))
    if best_identity < min_identity:
        assignment, confidence = "unclassified", "none"
    elif best_identity >= species_identity and margin is not None and margin >= species_margin:
        assignment, confidence = "species", "high" if best_identity >= 0.99 else "medium"
    elif best_identity >= species_identity and margin is None:
        assignment, confidence = "provisional", "medium"
    else:
        assignment, confidence = "ambiguous", "low"
    return TaxonomicHit(
        cluster_id=cluster_id, taxon=best.taxon, rank=best.rank, accession=best.accession,
        identity=best_identity, query_coverage=coverage, margin_to_second=margin,
        assignment=assignment, confidence=confidence,
    )

