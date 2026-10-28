from __future__ import annotations

import json
import statistics
from importlib.resources import files
from pathlib import Path

from .models import PrimerDetection, PrimerProfile
from .sequence import find_amplicon


def load_primer_profiles(path: Path | None = None) -> list[PrimerProfile]:
    source = path if path else files("metabarcoder").joinpath("data/primer_profiles.json")
    data = json.loads(source.read_text(encoding="utf-8"))
    return [PrimerProfile(
        id=row["id"], name=row["name"], marker=row["marker"],
        forward_primer=row["forward_primer"], reverse_primer=row["reverse_primer"],
        min_insert_bp=row["min_insert_bp"], max_insert_bp=row["max_insert_bp"],
        expected_insert_bp=row["expected_insert_bp"],
        max_primer_error_rate=row.get("max_primer_error_rate", 0.15),
        taxonomic_scope=tuple(row.get("taxonomic_scope", [])), doi=row.get("doi", ""),
    ) for row in data]


def detect_primers(sequences: list[str], profiles: list[PrimerProfile]) -> list[PrimerDetection]:
    detections: list[PrimerDetection] = []
    for profile in profiles:
        max_errors = max(1, round(max(len(profile.forward_primer), len(profile.reverse_primer)) *
                                  profile.max_primer_error_rate))
        lengths: list[int] = []
        for sequence in sequences:
            found = find_amplicon(sequence, profile.forward_primer, profile.reverse_primer,
                                  max_errors, profile.min_insert_bp, profile.max_insert_bp)
            if found:
                lengths.append(len(found[0]))
        support = len(lengths) / len(sequences) if sequences else 0.0
        confidence = "high" if support >= 0.30 else "medium" if support >= 0.05 else "low" if support > 0 else "none"
        detections.append(PrimerDetection(
            profile=profile, sampled_reads=len(sequences), paired_hits=len(lengths),
            support_fraction=support, confidence=confidence,
            median_insert_bp=round(statistics.median(lengths)) if lengths else None,
        ))
    return sorted(detections, key=lambda item: (item.support_fraction, item.paired_hits), reverse=True)

