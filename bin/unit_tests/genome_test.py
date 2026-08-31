import os
import pathlib
import sys
import unittest
from typing import Generator

from Bio import SeqIO

sys.path.append(str(pathlib.Path(__file__).parent.parent.parent))

from bin.Genome import Genome

# import the module instead of the class so that unittest does not collect and
# run ResultsTest as part of this module
import results_test


class GenomeTest(unittest.TestCase):
    """class for testing the Genome class"""

    # constants
    TEST_DIR = results_test.ResultsTest.TEST_DIR
    GENBANK = "genbank"
    FASTA = "fasta"

    # keys for looking up each of the three input files
    GBK = "gbk"
    FNA = "fna"
    EMPTY = "empty"

    # basenames for each of the three input files
    BASENAMES = {
        GBK: results_test.ResultsTest.INGROUP_FILES[0],
        FNA: "i1.fna",
        EMPTY: "empty.txt",
    }

    # the ftp path of the genbank file to download (Mycoplasma mycoides)
    GBK_FTP = results_test.ResultsTest.DEFAULT_INGROUP[
        results_test.ResultsTest.INGROUP_FILES[0]
    ]

    @classmethod
    def setUpClass(cls) -> None:
        """downloads a genbank file, derives a fasta file and an empty file from
        it, then builds a persistent Genome object for each file both with and
        without preloading
        """
        # make the test directory if it doesn't exist
        if not os.path.isdir(GenomeTest.TEST_DIR):
            os.mkdir(GenomeTest.TEST_DIR)

        # get the path of each input file
        cls.files: dict[str, pathlib.Path] = {
            key: pathlib.Path(os.path.join(GenomeTest.TEST_DIR, basename))
            for key, basename in GenomeTest.BASENAMES.items()
        }

        # download the genbank file
        results_test.ResultsTest._downloadOneGenome(
            GenomeTest.GBK_FTP, str(cls.files[GenomeTest.GBK])
        )

        # convert the genbank file to a fasta file
        SeqIO.convert(
            str(cls.files[GenomeTest.GBK]),
            GenomeTest.GENBANK,
            str(cls.files[GenomeTest.FNA]),
            GenomeTest.FASTA,
        )

        # make an empty file
        with open(cls.files[GenomeTest.EMPTY], "w"):
            pass

        # the expected format of each file; an empty file has no format
        cls.formats: dict[str, str | None] = {
            GenomeTest.GBK: GenomeTest.GENBANK,
            GenomeTest.FNA: GenomeTest.FASTA,
            GenomeTest.EMPTY: None,
        }

        # get the expected number of records and total length directly from
        # biopython so that the expectations are independent of the Genome class;
        # the fasta file was derived from the genbank file so both files share
        # the same records
        with open(cls.files[GenomeTest.GBK], "r") as fh:
            records = list(SeqIO.parse(fh, GenomeTest.GENBANK))

        cls.numSeqs: dict[str, int] = {
            GenomeTest.GBK: len(records),
            GenomeTest.FNA: len(records),
            GenomeTest.EMPTY: 0,
        }

        cls.lengths: dict[str, int] = {
            GenomeTest.GBK: sum(map(len, records)),
            GenomeTest.FNA: sum(map(len, records)),
            GenomeTest.EMPTY: 0,
        }

        # build a persistent Genome object for each file with and without preloading
        cls.genomes: dict[tuple[str, bool], Genome] = {
            (key, preload): Genome(fn, preload)
            for key, fn in cls.files.items()
            for preload in (False, True)
        }

    @classmethod
    def tearDownClass(cls) -> None:
        """removes the files that were created for these tests"""
        for fn in GenomeTest.files.values():
            if fn.exists():
                os.remove(fn)

    # helper functions
    def _iterGenomes(self) -> Generator[tuple[str, bool, Genome], None, None]:
        """iterates over all six Genome objects

        Yields:
            Generator[tuple[str,bool,Genome],None,None]: the file key, the
            preload status, and the Genome object
        """
        for key in (GenomeTest.GBK, GenomeTest.FNA, GenomeTest.EMPTY):
            for preload in (False, True):
                yield key, preload, self.genomes[(key, preload)]

    @staticmethod
    def _getMsg(prefix: str, key: str, preload: bool) -> str:
        """builds a failure message that identifies one Genome object

        Args:
            prefix (str): the start of the message
            key (str): the file key
            preload (bool): the preload status

        Returns:
            str: the failure message
        """
        return f"{prefix} for {GenomeTest.BASENAMES[key]} (preload={preload})"

    # test cases
    def testA_check_format(self) -> None:
        """does the sniffed file format match the expected format"""
        # constant
        FAIL_MSG = "wrong format"

        # the format is sniffed from the first line of the file
        for key, preload, genome in self._iterGenomes():
            self.assertEqual(
                genome._format,
                self.formats[key],
                GenomeTest._getMsg(FAIL_MSG, key, preload),
            )

    def testB_check_name(self) -> None:
        """does the name attribute match the filename without its extension"""
        # constant
        FAIL_MSG = "wrong name"

        # the name is the stem of the file
        for key, preload, genome in self._iterGenomes():
            self.assertEqual(
                genome.name,
                self.files[key].name,
                GenomeTest._getMsg(FAIL_MSG, key, preload),
            )

    def testC_check_filename(self) -> None:
        """does the fn attribute match the file the Genome was built from"""
        # constant
        FAIL_MSG = "wrong filename"

        # the fn attribute is the file passed to the constructor
        for key, preload, genome in self._iterGenomes():
            self.assertEqual(
                genome.fn,
                self.files[key],
                GenomeTest._getMsg(FAIL_MSG, key, preload),
            )

    def testD_check_len(self) -> None:
        """does __len__ return the number of records in the file"""
        # constant
        FAIL_MSG = "wrong number of records"

        # an empty file has no records
        for key, preload, genome in self._iterGenomes():
            self.assertEqual(
                len(genome),
                self.numSeqs[key],
                GenomeTest._getMsg(FAIL_MSG, key, preload),
            )

    def testE_check_length(self) -> None:
        """does the length attribute equal the total length of all records"""
        # constant
        FAIL_MSG = "wrong length"

        # an empty file has no sequence
        for key, preload, genome in self._iterGenomes():
            self.assertEqual(
                genome.length,
                self.lengths[key],
                GenomeTest._getMsg(FAIL_MSG, key, preload),
            )

    def testF_check_equality(self) -> None:
        """are two Genomes equal only when they were built from the same file"""
        # constants
        FAIL_EQ = "genomes built from the same file are not equal"
        FAIL_NE = "genomes built from different files are equal"
        FAIL_HASH = "genomes built from the same file have different hashes"

        for key, preload, genome in self._iterGenomes():
            # equality depends on the file, so the preload status is irrelevant
            same = self.genomes[(key, not preload)]
            self.assertEqual(genome, same, GenomeTest._getMsg(FAIL_EQ, key, preload))
            self.assertFalse(
                genome != same, GenomeTest._getMsg(FAIL_EQ, key, preload)
            )

            # equal genomes must hash identically so they can share a dictionary key
            self.assertEqual(
                hash(genome), hash(same), GenomeTest._getMsg(FAIL_HASH, key, preload)
            )

            # a genome built from any other file is not equal
            for otherKey, otherPreload, other in self._iterGenomes():
                if otherKey == key:
                    continue

                self.assertNotEqual(
                    genome, other, GenomeTest._getMsg(FAIL_NE, key, preload)
                )
                self.assertTrue(
                    genome != other, GenomeTest._getMsg(FAIL_NE, key, preload)
                )

    def testG_check_gt(self) -> None:
        """does __gt__ compare the total length of the genomes"""
        # constants
        FAIL_SAME = "a genome is greater than a genome of equal length"
        FAIL_GT = "a genome is not greater than an empty genome"
        FAIL_LT = "an empty genome is greater than a genome with sequence"

        for preload in (False, True):
            gbk = self.genomes[(GenomeTest.GBK, preload)]
            fna = self.genomes[(GenomeTest.FNA, preload)]
            empty = self.genomes[(GenomeTest.EMPTY, preload)]

            # the fasta file was derived from the genbank file, so the two
            # genomes have the same length and neither one is greater
            self.assertFalse(
                gbk > fna, GenomeTest._getMsg(FAIL_SAME, GenomeTest.GBK, preload)
            )
            self.assertFalse(
                fna > gbk, GenomeTest._getMsg(FAIL_SAME, GenomeTest.FNA, preload)
            )

            # an empty file has no sequence, so it is shorter than both genomes
            self.assertTrue(
                gbk > empty, GenomeTest._getMsg(FAIL_GT, GenomeTest.GBK, preload)
            )
            self.assertTrue(
                fna > empty, GenomeTest._getMsg(FAIL_GT, GenomeTest.FNA, preload)
            )
            self.assertFalse(
                empty > gbk, GenomeTest._getMsg(FAIL_LT, GenomeTest.EMPTY, preload)
            )
            self.assertFalse(
                empty > fna, GenomeTest._getMsg(FAIL_LT, GenomeTest.EMPTY, preload)
            )

    def testH_check_preload(self) -> None:
        """is _seqs populated only when the sequences were preloaded"""
        # constants
        FAIL_LOADED = "sequences were not preloaded"
        FAIL_STREAM = "sequences were loaded into memory instead of streamed"
        FAIL_RECORD = "wrong record preloaded"

        for key, preload, genome in self._iterGenomes():
            # a file without a format has no sequences to preload
            if preload and self.formats[key] is not None:
                self.assertEqual(
                    len(genome._seqs),
                    self.numSeqs[key],
                    GenomeTest._getMsg(FAIL_LOADED, key, preload),
                )

                # the preloaded records must match those in the file
                with open(self.files[key], "r") as fh:
                    expected = list(SeqIO.parse(fh, self.formats[key]))

                for obs, exp in zip(genome._seqs, expected):
                    self.assertEqual(
                        obs.id, exp.id, GenomeTest._getMsg(FAIL_RECORD, key, preload)
                    )
                    self.assertEqual(
                        obs.seq, exp.seq, GenomeTest._getMsg(FAIL_RECORD, key, preload)
                    )

            # without preloading the sequences are always streamed from file
            else:
                self.assertEqual(
                    genome._seqs, [], GenomeTest._getMsg(FAIL_STREAM, key, preload)
                )

    def testI_check_clear(self) -> None:
        """does _clear remove the sequences from memory without breaking iteration"""
        # constants
        FAIL_CLEAR = "sequences were not removed from memory"
        FAIL_STREAM = "wrong number of records streamed after clearing"
        FAIL_RECORD = "wrong record streamed after clearing"

        # build new objects so that the persistent genomes are not modified
        for key, fn in self.files.items():
            for preload in (False, True):
                genome = Genome(fn, preload)

                # clearing always empties the sequences in memory
                genome._clear()
                self.assertEqual(
                    genome._seqs, [], GenomeTest._getMsg(FAIL_CLEAR, key, preload)
                )

                # the records must still be streamed from file after clearing
                observed = list(genome)
                self.assertEqual(
                    len(observed),
                    self.numSeqs[key],
                    GenomeTest._getMsg(FAIL_STREAM, key, preload),
                )

                # an empty file has no records to compare
                if self.formats[key] is None:
                    continue

                with open(fn, "r") as fh:
                    expected = list(SeqIO.parse(fh, self.formats[key]))

                for obs, exp in zip(observed, expected):
                    self.assertEqual(
                        obs.id, exp.id, GenomeTest._getMsg(FAIL_RECORD, key, preload)
                    )
                    self.assertEqual(
                        obs.seq, exp.seq, GenomeTest._getMsg(FAIL_RECORD, key, preload)
                    )


if __name__ == "__main__":
    unittest.main()
