from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

from metabarcoder.models import AnalysisConfig
from metabarcoder.pipeline import run_analysis
from metabarcoder.primers import detect_primers, load_primer_profiles
from metabarcoder.sequence import find_amplicon, reverse_complement


INSERT = "TAACAAACTAGGAGGAGTACTAGCCCTAATCTCATCAATCCTAATCTTGATCCTTATACCCCTCCTCCATACATCTAAACAACGCAGTATAATGTTCCGGCCATTTAGTCAATGCCTATTCTGAATCCTAGTAGCTGACCTATTAACACT"
FORWARD = "TACGCAATCCTACGATCAATTCC"
REVERSE = "GGTTGTCCTCCAATTCATGTTAG"
READ = "ACGTACGT" + FORWARD + INSERT + reverse_complement(REVERSE) + "TTAA"


class CoreTests(unittest.TestCase):
    def test_finds_amplicon_in_both_orientations(self):
        for sequence in (READ, reverse_complement(READ)):
            found = find_amplicon(sequence, FORWARD, REVERSE, 3, 125, 175)
            self.assertIsNotNone(found)
            self.assertEqual(found[0], INSERT)

    def test_detects_known_profile(self):
        result = detect_primers([READ] * 20, load_primer_profiles())[0]
        self.assertEqual(result.profile.id, "cytb_l15601_h15748")
        self.assertEqual(result.confidence, "high")
        self.assertEqual(result.median_insert_bp, 150)

    def test_end_to_end_gz_fastq(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fastq = root / "sample.fastq.gz"
            with gzip.open(fastq, "wt") as handle:
                for index in range(12):
                    handle.write(f"@read{index}\n{READ}\n+\n{'I' * len(READ)}\n")
            reference = root / "reference.fasta"
            reference.write_text(">test|taxon=Capreolus capreolus|rank=species|marker=MT-CYB\n" + INSERT + "\n")
            result = run_analysis(AnalysisConfig(
                sample_name="AMP25-test", inputs=[fastq], reference_path=reference,
                output_dir=root / "results", min_cluster_fraction=0,
            ))
            self.assertEqual(result["primer_detection"]["profile_id"], "cytb_l15601_h15748")
            self.assertEqual(result["taxonomy"][0]["taxon"], "Capreolus capreolus")
            self.assertEqual(result["taxonomy"][0]["assignment"], "provisional")
            self.assertTrue((root / "results/report.html").exists())


if __name__ == "__main__":
    unittest.main()

