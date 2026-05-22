import tempfile
import unittest
from pathlib import Path

from scripts.annotate_codon_positions import annotate_alignment, find_cds_in_seq
from scripts.validate_translations import (
    compare_proteins,
    extract_annotated_cds_from_alignment,
    translate_cds,
    translated_ref_start,
)


class AnnotateCoordinateTests(unittest.TestCase):
    def test_find_cds_in_seq_trimmed_returns_offset(self):
        mode, start, end = find_cds_in_seq("AAACCC", "ATGAAACCCTAA")

        self.assertEqual(mode, "trimmed")
        self.assertEqual(start, 3)
        self.assertEqual(end, 9)

    def test_annotate_alignment_threads_codon_positions_across_gaps(self):
        annotations = annotate_alignment("AA--ACCC", cds_start=0, cds_end=6, codon_offset=3)

        self.assertEqual(annotations, ["1", "2", "-", "-", "3", "1", "2", "3"])


class ValidationCoordinateTests(unittest.TestCase):
    def test_translated_reference_start_advances_after_partial_codon_skip(self):
        self.assertEqual(translated_ref_start(0), 0)
        self.assertEqual(translated_ref_start(1), 1)
        self.assertEqual(translated_ref_start(2), 1)
        self.assertEqual(translated_ref_start(3), 1)

    def test_compare_proteins_uses_start_after_leading_skip(self):
        translated = translate_cds("CCATG", codon_offset=1)
        result = compare_proteins(translated, "XM", codon_offset=1)

        self.assertEqual(translated, "M")
        self.assertEqual(result["match"], "exact")

    def test_extract_annotated_cds_reads_tsv_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "gene.dna.aln.fasta"
            tsv = tmp_path / "gene_symbol_codon_pos.tsv"
            fasta.write_text(
                ">1 Gene_Name:SYMBOL Species:Arabidopsis_thaliana\n"
                "ATG--AAAT\n"
                ">2 Gene_Name:SYMBOL Species:Other_species\n"
                "ATGCCAAAT\n"
            )
            tsv.write_text(
                "# Gene: 1 (SYMBOL)\n"
                "col_index\tcodon_pos\n"
                "0\t1\n"
                "1\t2\n"
                "2\t3\n"
                "3\t-\n"
                "4\t-\n"
                "5\t1\n"
                "6\t2\n"
                "7\t3\n"
                "8\tU\n"
            )

            seq, diagnostics = extract_annotated_cds_from_alignment(fasta, tsv)

        self.assertEqual(seq, "ATGAAA")
        self.assertEqual(diagnostics["coding_columns"], 6)
        self.assertEqual(diagnostics["coding_gaps"], 0)
        self.assertEqual(diagnostics["phase_errors"], [])


if __name__ == "__main__":
    unittest.main()
