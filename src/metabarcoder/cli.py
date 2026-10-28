from __future__ import annotations

import argparse
import json
from pathlib import Path

from .models import AnalysisConfig
from .pipeline import run_analysis


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="metabarcoder")
    subparsers = parser.add_subparsers(dest="command", required=True)
    analyze = subparsers.add_parser("analyze", help="Analyzovať FASTQ/FASTQ.gz")
    analyze.add_argument("--input", action="append", type=Path, required=True,
                         help="FASTQ súbor alebo priečinok; možno opakovať")
    analyze.add_argument("--reference", type=Path, required=True)
    analyze.add_argument("--output", type=Path, required=True)
    analyze.add_argument("--sample", default="sample")
    analyze.add_argument("--max-reads", type=int)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "analyze":
        result = run_analysis(AnalysisConfig(
            sample_name=args.sample, inputs=args.input, reference_path=args.reference,
            output_dir=args.output, max_reads=args.max_reads,
        ))
        print(json.dumps({
            "sample": result["sample"], "status": result["status"],
            "primer": result["primer_detection"]["profile_id"],
            "results": result["taxonomy"], "report": str(args.output / "report.html"),
        }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

