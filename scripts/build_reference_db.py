#!/usr/bin/env python3
"""Build a traceable animal mitochondrial reference panel.

The default build retrieves complete RefSeq mitochondrial records and extracts
annotated MT-CYB, COI, 12S and 16S features. GenBank fallback is opt-in. A
supplemental FASTA is accepted only when every record carries explicit source,
accession, taxon, marker and audit metadata in the standard header format.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Seq import UndefinedSequenceError

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from metabarcoder.sequence import find_amplicon  # noqa: E402


def marker_for(feature) -> str | None:
    genes = {value.upper() for value in feature.qualifiers.get("gene", [])}
    products = {value.upper() for value in feature.qualifiers.get("product", [])}
    notes = " ".join(feature.qualifiers.get("note", [])).upper()
    if genes & {"CYTB", "MT-CYB", "COB"} or "CYTOCHROME B" in products:
        return "MT-CYB"
    if genes & {"COI", "CO1", "COX1", "MT-CO1"} or products & {
        "CYTOCHROME C OXIDASE SUBUNIT I", "CYTOCHROME C OXIDASE SUBUNIT 1"
    }:
        return "COI"
    if genes & {"12S", "MT-RNR1", "RNR1"} or products & {"12S RIBOSOMAL RNA", "S-RRNA"} or "12S RIBOSOMAL RNA" in notes:
        return "12S"
    if genes & {"16S", "MT-RNR2", "RNR2"} or products & {"16S RIBOSOMAL RNA", "L-RRNA"} or "16S RIBOSOMAL RNA" in notes:
        return "16S"
    return None


def search_accessions(species: str, allow_genbank: bool) -> tuple[list[str], str]:
    curated = f'"{species}"[Organism] AND mitochondrion[filter] AND refseq[filter]'
    with Entrez.esearch(db="nuccore", term=curated, retmax=5) as response:
        ids = Entrez.read(response)["IdList"]
    if ids:
        return ids, "curated"
    if not allow_genbank:
        return [], "missing"
    fallback = f'"{species}"[Organism] AND mitochondrion[filter] AND "complete genome"[Title]'
    with Entrez.esearch(db="nuccore", term=fallback, retmax=5) as response:
        return Entrez.read(response)["IdList"], "traceable"


def fetch_records(ids: list[str]):
    if not ids:
        return []
    with Entrez.efetch(db="nuccore", id=",".join(ids), rettype="gb", retmode="text") as response:
        return list(SeqIO.parse(response, "genbank"))


def load_profiles(path: Path) -> list[dict]:
    return json.loads(path.read_text(encoding="utf-8"))


def amplicons_for(sequence: str, marker: str, profiles: list[dict]):
    for profile in profiles:
        if profile["marker"] != marker:
            continue
        longest_primer = max(len(profile["forward_primer"]), len(profile["reverse_primer"]))
        max_errors = max(1, round(longest_primer * profile["max_primer_error_rate"]))
        result = find_amplicon(sequence, profile["forward_primer"], profile["reverse_primer"],
                               max_errors, profile["min_insert_bp"], profile["max_insert_bp"])
        if result:
            yield profile["id"], result[0]


def extract_record(record, species: str, audit: str, profiles: list[dict]):
    try:
        if len(record.seq) < 10_000:
            return []
        str(record.seq)
    except UndefinedSequenceError:
        return []
    extracted = []
    accession = record.id
    source = "NCBI RefSeq" if audit == "curated" else "NCBI GenBank"
    for feature in record.features:
        marker = marker_for(feature)
        if marker is None:
            continue
        expected_type = "CDS" if marker in {"MT-CYB", "COI"} else "rRNA"
        if feature.type != expected_type:
            continue
        try:
            feature_sequence = str(feature.extract(record.seq)).upper()
        except UndefinedSequenceError:
            continue
        for profile_id, sequence in amplicons_for(feature_sequence, marker, profiles):
            fields = {"accession": accession, "taxon": species, "marker": marker,
                      "primer_profile": profile_id, "source": source, "audit": audit,
                      "length": str(len(sequence))}
            header = standard_header(accession, species, marker, source, audit)
            extracted.append((f"{header}|primer_profile={profile_id}", sequence, fields))
    return extracted


def standard_header(accession: str, taxon: str, marker: str, source: str, audit: str) -> str:
    return f"{accession}|accession={accession}|taxon={taxon}|rank=species|marker={marker}|source={source}|audit={audit}"


def parse_standard_header(header: str) -> dict[str, str]:
    return dict(part.split("=", 1) for part in header.split("|") if "=" in part)


def import_supplement(path: Path) -> list[tuple[str, str, dict]]:
    imported = []
    for record in SeqIO.parse(path, "fasta"):
        fields = parse_standard_header(record.description)
        required = {"accession", "taxon", "marker", "source", "audit"}
        missing = required - fields.keys()
        if missing:
            raise ValueError(f"Supplement {record.id} lacks metadata: {sorted(missing)}")
        if fields["source"] not in {"MIDORI2", "BOLD"}:
            raise ValueError(f"Unsupported supplement source: {fields['source']}")
        imported.append((record.description, str(record.seq).upper(), fields))
    return imported


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", type=Path, default=Path("references/manifests/target_species.tsv"))
    parser.add_argument("--output", type=Path, default=Path("references/releases/nu3gen_animals_v1.fasta"))
    parser.add_argument("--manifest", type=Path, default=Path("references/releases/nu3gen_animals_v1.records.tsv"))
    parser.add_argument("--profiles", type=Path, default=Path("src/metabarcoder/data/primer_profiles.json"))
    parser.add_argument("--fallback-accessions", type=Path,
                        default=Path("references/manifests/fallback_accessions.tsv"))
    parser.add_argument("--email", default=os.getenv("NCBI_EMAIL"), help="Required by NCBI E-utilities")
    parser.add_argument("--api-key", default=os.getenv("NCBI_API_KEY"))
    parser.add_argument("--allow-genbank-fallback", action="store_true")
    parser.add_argument("--supplement-fasta", action="append", type=Path, default=[])
    args = parser.parse_args()
    if not args.email:
        parser.error("Provide --email or NCBI_EMAIL")
    Entrez.email = args.email
    Entrez.api_key = args.api_key

    targets = list(csv.DictReader(args.species.open(encoding="utf-8"), delimiter="\t"))
    profiles = load_profiles(args.profiles)
    sequences: list[tuple[str, str, dict]] = []
    for target in targets:
        species = target["scientific_name"]
        ids, audit = search_accessions(species, args.allow_genbank_fallback)
        for record in fetch_records(ids):
            sequences.extend(extract_record(record, species, audit, profiles))
        time.sleep(0.11 if args.api_key else 0.35)

    present = {row[2]["taxon"] for row in sequences}
    if args.fallback_accessions.exists():
        for fallback in csv.DictReader(args.fallback_accessions.open(encoding="utf-8"), delimiter="\t"):
            species = fallback["scientific_name"]
            if species in present:
                continue
            for record in fetch_records([fallback["accession"]]):
                sequences.extend(extract_record(record, species, "traceable", profiles))
            time.sleep(0.11 if args.api_key else 0.35)

    for supplement in args.supplement_fasta:
        sequences.extend(import_supplement(supplement))

    unique = {}
    for header, sequence, fields in sequences:
        unique[(fields["taxon"], fields["marker"], sequence)] = (header, sequence, fields)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as fasta, args.manifest.open("w", newline="", encoding="utf-8") as manifest:
        columns = ["accession", "taxon", "marker", "primer_profile", "source", "audit", "length"]
        writer = csv.DictWriter(manifest, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for header, sequence, fields in sorted(unique.values(), key=lambda row: (row[2]["taxon"], row[2]["marker"], row[2]["accession"])):
            fasta.write(f">{header}\n{sequence}\n")
            fields["length"] = str(len(sequence))
            writer.writerow(fields)
    summary = {"records": len(unique), "species": len({row[2]["taxon"] for row in unique.values()}),
               "sources": sorted({row[2]["source"] for row in unique.values()})}
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
