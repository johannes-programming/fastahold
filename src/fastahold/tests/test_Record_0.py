import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from fastahold.core.Record import Record


class TestRecord(unittest.TestCase):

    def test_init_fields(self):
        record = Record(description="Test Description", seq="ATGC")
        self.assertEqual(record.description, "Test Description")
        self.assertEqual(str(record.seq), "ATGC")

    def test_init_text(self):
        fasta_text = ">Test Description\nATGC"
        record = Record(text=fasta_text)
        self.assertEqual(record.description, "Test Description")
        self.assertEqual(str(record.seq), "ATGC")

    def test_init_bio(self):
        bio_seq_record = SeqRecord(
            Seq("ATGC"), id="TestID", description="Test Description"
        )
        record = Record(bio=bio_seq_record)
        self.assertEqual(record.description, "Test Description")
        self.assertEqual(str(record.seq), "ATGC")

    def test_copy(self):
        record = Record(description="Test Description", seq="ATGC")
        copied_record = record.copy()
        self.assertEqual(copied_record.description, "Test Description")
        self.assertEqual(str(copied_record.seq), "ATGC")
        self.assertIsNot(record, copied_record)  # Ensure they are different objects

    def test_id_property(self):
        record = Record(description="TestID Description", seq="ATGC")
        self.assertEqual(record.id, "TestID")
        record.id = "NewID"
        self.assertEqual(record.id, "NewID")

    def test_supplementary_property(self):
        record = Record(description="TestID Supplementary", seq="ATGC")
        self.assertEqual(record.supplementary, "Supplementary")
        record.supplementary = "New Info"
        self.assertEqual(record.supplementary, "New Info")

    def test_bio_property(self):
        record = Record(description="Test Description", seq="ATGC")
        bio_record = record.bio
        self.assertEqual(bio_record.description, "Test Description")
        self.assertEqual(str(bio_record.seq), "ATGC")

    def test_text_property(self):
        record = Record(description="Test Description", seq="ATGC")
        expected_text = ">Test Description\nATGC\n"
        self.assertEqual(record.text, expected_text)

    def test_text_setter(self):
        record = Record()
        record.text = ">Test Description\nATGC"
        self.assertEqual(record.description, "Test Description")
        self.assertEqual(str(record.seq), "ATGC")

    def test_seq_setter(self):
        record = Record()
        record.seq = "ATGC"
        self.assertEqual(str(record.seq), "ATGC")

    def test_invalid_init(self):
        with self.assertRaises(TypeError):
            Record(123)

    def test_invalid_copy(self):
        with self.assertRaises(TypeError):
            Record().copy("invalid argument")


if __name__ == "__main__":
    unittest.main()
