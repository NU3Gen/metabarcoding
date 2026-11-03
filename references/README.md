# NU3Gen Animal Reference Database

The production reference is generated, not hand-edited. `target_species.tsv`
defines the initial food-relevant and Slovak wildlife panel. The builder uses
NCBI RefSeq complete mitochondrial records as the primary source. From their
annotated marker features, the builder exports only amplicons delimited by the
primer profiles currently supported by the pipeline. This prevents a short
consensus from being compared with an entire mitochondrial gene.

MIDORI2 and BOLD records may be added through `--supplement-fasta`, but only in
the standardized header format containing `accession`, `taxon`, `marker`,
`source` and `audit`. This keeps every call traceable to its source record.

```bash
NCBI_EMAIL=name@example.org python scripts/build_reference_db.py
```

Pinned GenBank fallbacks in `manifests/fallback_accessions.tsv` are used only
when a target has no usable RefSeq amplicon. Such records are labelled
`traceable`, never `curated`. NCBI `nt` is available
from each report for independent BLAST verification, but its mixed-review
content is not automatically promoted into the local curated database.

The generated release FASTA and record manifest should be versioned together
with a checksum and the exact MIDORI2 release DOI/BOLD query when supplements
are used.

Version 1 contains 22 amplicons for 16 target taxa. Nine target taxa do not
produce an in-range amplicon with the currently implemented L15601/H15748
profile and remain listed in the target manifest for a future primer profile.
Absence from this panel must not be reported as proof of absence from a sample.
