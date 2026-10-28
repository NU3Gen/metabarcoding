from __future__ import annotations

import csv
import html
import json
from collections import Counter
from dataclasses import asdict
from pathlib import Path

from .models import AnalysisConfig, Cluster, FastqStats
from .primers import detect_primers, load_primer_profiles
from .references import classify, load_reference_fasta
from .sequence import discover_fastq, find_amplicon, iter_fastq, levenshtein, mean_phred


def _cluster_variants(counts: Counter[str], total: int, identity_threshold: float,
                      min_reads: int, min_fraction: float) -> list[Cluster]:
    buckets: list[dict] = []
    for sequence, count in counts.most_common():
        target = None
        for bucket in buckets:
            max_len = max(len(sequence), len(bucket["representative"]))
            distance = levenshtein(sequence, bucket["representative"],
                                   max_distance=max(1, round(max_len * (1 - identity_threshold))))
            if 1 - distance / max_len >= identity_threshold:
                target = bucket
                break
        if target is None:
            buckets.append({"representative": sequence, "reads": count, "variants": 1})
        else:
            target["reads"] += count
            target["variants"] += 1
    clusters: list[Cluster] = []
    for index, bucket in enumerate(sorted(buckets, key=lambda x: x["reads"], reverse=True), 1):
        fraction = bucket["reads"] / total if total else 0.0
        if bucket["reads"] >= min_reads and fraction >= min_fraction:
            clusters.append(Cluster(
                id=f"cluster_{index}", representative=bucket["representative"],
                reads=bucket["reads"], fraction=fraction, unique_variants=bucket["variants"],
            ))
    return clusters


def run_analysis(config: AnalysisConfig) -> dict:
    paths = discover_fastq(config.inputs)
    profiles = load_primer_profiles()
    sampled_sequences: list[str] = []
    stats = FastqStats(files=[str(path) for path in paths])
    for _, sequence, quality in iter_fastq(paths):
        q = mean_phred(quality)
        stats.add(len(sequence), q)
        if len(sampled_sequences) < config.detection_sample_reads and q >= config.min_mean_q:
            sampled_sequences.append(sequence)
        if config.max_reads and stats.reads_seen >= config.max_reads:
            break

    detections = detect_primers(sampled_sequences, profiles)
    if not detections or detections[0].confidence == "none":
        raise ValueError("V čítaniach sa nenašiel žiadny podporovaný primerový pár")
    selected = detections[0]
    profile = selected.profile
    max_errors = max(1, round(max(len(profile.forward_primer), len(profile.reverse_primer)) *
                              profile.max_primer_error_rate))

    inserts: Counter[str] = Counter()
    reads_examined = 0
    reads_below_quality = 0
    reads_without_pair = 0
    for _, sequence, quality in iter_fastq(paths):
        reads_examined += 1
        if mean_phred(quality) < config.min_mean_q:
            reads_below_quality += 1
            continue
        found = find_amplicon(sequence, profile.forward_primer, profile.reverse_primer,
                              max_errors, profile.min_insert_bp, profile.max_insert_bp)
        if found is None:
            reads_without_pair += 1
        else:
            inserts[found[0]] += 1
        if config.max_reads and reads_examined >= config.max_reads:
            break

    filtered_total = sum(inserts.values())
    if filtered_total == 0:
        raise ValueError("Po primerovom a kvalitatívnom filtri nezostali žiadne čítania")
    clusters = _cluster_variants(inserts, filtered_total, config.cluster_identity,
                                 config.min_cluster_reads, config.min_cluster_fraction)
    references = load_reference_fasta(config.reference_path)
    hits = [classify(cluster.representative, cluster.id, references) for cluster in clusters]

    result = {
        "schema_version": "1.0",
        "pipeline_version": "0.1.0",
        "sample": config.sample_name,
        "status": "completed",
        "input": {
            "files": stats.files, "reads": stats.reads_seen, "bases": stats.bases_seen,
            "mean_read_length": round(stats.mean_length, 2),
            "mean_read_quality": round(stats.mean_quality, 2),
            "min_read_length": stats.min_length, "max_read_length": stats.max_length,
        },
        "primer_detection": {
            "profile_id": profile.id, "name": profile.name, "marker": profile.marker,
            "sampled_reads": selected.sampled_reads, "paired_hits": selected.paired_hits,
            "support_fraction": round(selected.support_fraction, 6),
            "confidence": selected.confidence, "median_insert_bp": selected.median_insert_bp,
        },
        "filtering": {
            "reads_examined": reads_examined, "reads_below_quality": reads_below_quality,
            "reads_without_primer_pair": reads_without_pair,
            "reads_with_valid_insert": filtered_total,
        },
        "clusters": [asdict(cluster) for cluster in clusters],
        "taxonomy": [asdict(hit) for hit in hits],
        "warnings": [
            "Taxonomické priradenie je platné iba v rozsahu použitej referenčnej databázy.",
            "MVP nepoužíva negatívnu kontrolu ani produkčný chiméra/NUMT filter.",
        ],
    }
    write_outputs(result, config.output_dir)
    return result


def write_outputs(result: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "report.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    with (output_dir / "species.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cluster", "reads", "fraction", "taxon", "rank", "identity", "coverage", "assignment", "confidence"])
        clusters = {row["id"]: row for row in result["clusters"]}
        for hit in result["taxonomy"]:
            cluster = clusters[hit["cluster_id"]]
            writer.writerow([hit["cluster_id"], cluster["reads"], f'{cluster["fraction"]:.6f}',
                             hit["taxon"], hit["rank"], f'{hit["identity"]:.6f}',
                             f'{hit["query_coverage"]:.6f}', hit["assignment"], hit["confidence"]])
    rows = []
    clusters = {row["id"]: row for row in result["clusters"]}
    for hit in result["taxonomy"]:
        cluster = clusters[hit["cluster_id"]]
        rows.append("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in (
            hit["cluster_id"], cluster["reads"], f'{cluster["fraction"]:.2%}', hit["taxon"],
            f'{hit["identity"]:.2%}', f'{hit["query_coverage"]:.2%}', hit["assignment"], hit["confidence"],
        )) + "</tr>")
    report_html = f"""<!doctype html><html lang='sk'><meta charset='utf-8'>
<meta name='viewport' content='width=device-width,initial-scale=1'><title>Metabarcoding – {html.escape(result['sample'])}</title>
<style>body{{font:16px system-ui;max-width:1100px;margin:40px auto;padding:0 20px;color:#17202a}}table{{border-collapse:collapse;width:100%}}th,td{{padding:10px;border-bottom:1px solid #ddd;text-align:left}}.card{{background:#f5f8fa;padding:20px;border-radius:14px;margin:18px 0}}code{{background:#eef2f4;padding:2px 5px}}</style>
<h1>Výsledok: {html.escape(result['sample'])}</h1>
<div class='card'><b>Marker:</b> {html.escape(result['primer_detection']['marker'])}<br>
<b>Primer:</b> {html.escape(result['primer_detection']['name'])}<br>
<b>Podpora detekcie:</b> {result['primer_detection']['support_fraction']:.2%}<br>
<b>Platné inserty:</b> {result['filtering']['reads_with_valid_insert']}</div>
<table><thead><tr><th>Zhluk</th><th>Čítania</th><th>Podiel</th><th>Taxón</th><th>Identita</th><th>Pokrytie</th><th>Priradenie</th><th>Istota</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
<p><small>{html.escape(' '.join(result['warnings']))}</small></p></html>"""
    (output_dir / "report.html").write_text(report_html, encoding="utf-8")

