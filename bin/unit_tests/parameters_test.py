import os
import pathlib
import sys
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.append(str(pathlib.Path(__file__).parent.parent.parent))

from bin.Genome import InvalidFileFormat
from bin.Parameters import DefaultArgs, Parameters, Log
from bin.cli import _parseArgs


class ParametersTest(unittest.TestCase):
    # default arguments
    DEFAULT_ARGS = DefaultArgs()

    # values for inputs
    IG_FNS_GB = ("itmp1.gb", "itmp2.gb", "itmp3.gb")
    IG_FNS_FA = ("itmp1.fa", "itmp2.fa", "itmp3.fa")
    OG_FNS_GB = ("otmp1.gb", "otmp2.gb", "otmp3.gb")
    OG_FNS_FA = ("otmp1.fa", "otmp2.fa", "otmp3.fa")

    # a mixture of genbank and fasta files; the format of each is auto-detected
    IG_FNS_MIX = (IG_FNS_GB[0], IG_FNS_FA[1], IG_FNS_GB[2])
    OG_FNS_MIX = (OG_FNS_FA[0], OG_FNS_GB[1], OG_FNS_FA[2])

    # genbank files whose sequences differ in length; listed longest first so
    # that the input order does not match the expected sorted order
    SIZE_FNS = ("stmp1.gb", "stmp2.gb", "stmp3.gb")
    SIZE_LENS = (40, 20, 60)

    # a file that is neither genbank nor fasta
    BAD_FMT_FN = "btmp1.txt"

    OUT_FN = "result.testfile"
    BAD_SIZE = "64,160"
    FORMAT_GB = "genbank"
    FORMAT_FA = "fasta"
    PRIMER_LEN = "15,25"
    GC_RANGE = "42,64"
    TM_RANGE = "68,72"
    PCR_SIZE = "80,128"
    TM_DIFF = "4.2"
    THREADS = "8"
    DUMP_FN = "tmp.p"
    VERSION = "parameters_test"
    AUTHOR = "null"
    TEST_SEQ = SeqRecord(
        Seq("atcg"),
        id="id",
        name="name",
        description="description",
        annotations={"molecule_type": "DNA"},
    )

    def setUp(self) -> None:
        """set up for each test"""
        # silence prints
        self.stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")

        # make the dummy files needed for testing
        ParametersTest._makeDummyFiles()

        # save the directory
        self.dir = pathlib.Path().cwd()

        # sys.argv for default values using short args
        self.basic1 = [
            "primerForge.py",
            "-i",
            *ParametersTest.IG_FNS_GB,
            "-o",
            ParametersTest.OUT_FN,
        ]

        # sys.argv for default values using long args
        self.basic2 = [
            "primerForge.py",
            "--ingroup",
            *ParametersTest.IG_FNS_GB,
            "--out",
            ParametersTest.OUT_FN,
        ]

        # sys.argv with custom values using short args for genbank file format
        self.short1 = [
            "primerForge.py",
            "-i",
            *ParametersTest.IG_FNS_GB,
            "-o",
            ParametersTest.OUT_FN,
            "-u",
            *ParametersTest.OG_FNS_GB,
            "-b",
            ParametersTest.BAD_SIZE,
            "-p",
            ParametersTest.PRIMER_LEN,
            "-g",
            ParametersTest.GC_RANGE,
            "-t",
            ParametersTest.TM_RANGE,
            "-r",
            ParametersTest.PCR_SIZE,
            "-d",
            ParametersTest.TM_DIFF,
            "-n",
            ParametersTest.THREADS,
        ]

        # sys.argv with custom values using short args for fasta file format
        self.short2 = [
            "primerForge.py",
            "-i",
            *ParametersTest.IG_FNS_FA,
            "-o",
            ParametersTest.OUT_FN,
            "-u",
            *ParametersTest.OG_FNS_FA,
            "-b",
            ParametersTest.BAD_SIZE,
            "-p",
            ParametersTest.PRIMER_LEN,
            "-g",
            ParametersTest.GC_RANGE,
            "-t",
            ParametersTest.TM_RANGE,
            "-r",
            ParametersTest.PCR_SIZE,
            "-d",
            ParametersTest.TM_DIFF,
            "-n",
            ParametersTest.THREADS,
        ]

        # sys.argv with custom values using long args for genbank file format
        self.long1 = [
            "primerForge.py",
            "--ingroup",
            *ParametersTest.IG_FNS_GB,
            "--out",
            ParametersTest.OUT_FN,
            "--outgroup",
            *ParametersTest.OG_FNS_GB,
            "--bad_sizes",
            ParametersTest.BAD_SIZE,
            "--primer_len",
            ParametersTest.PRIMER_LEN,
            "--gc_range",
            ParametersTest.GC_RANGE,
            "--tm_range",
            ParametersTest.TM_RANGE,
            "--pcr_prod",
            ParametersTest.PCR_SIZE,
            "--tm_diff",
            ParametersTest.TM_DIFF,
            "--num_threads",
            ParametersTest.THREADS,
        ]

        # sys.argv with custom values using long args for fasta file format
        self.long2 = [
            "primerForge.py",
            "--ingroup",
            *ParametersTest.IG_FNS_FA,
            "--out",
            ParametersTest.OUT_FN,
            "--outgroup",
            *ParametersTest.OG_FNS_FA,
            "--bad_sizes",
            ParametersTest.BAD_SIZE,
            "--primer_len",
            ParametersTest.PRIMER_LEN,
            "--gc_range",
            ParametersTest.GC_RANGE,
            "--tm_range",
            ParametersTest.TM_RANGE,
            "--pcr_prod",
            ParametersTest.PCR_SIZE,
            "--tm_diff",
            ParametersTest.TM_DIFF,
            "--num_threads",
            ParametersTest.THREADS,
        ]

        # sys.argv with a mixture of genbank and fasta files
        self.mixed = [
            "primerForge.py",
            "-i",
            *ParametersTest.IG_FNS_MIX,
            "-o",
            ParametersTest.OUT_FN,
            "-u",
            *ParametersTest.OG_FNS_MIX,
        ]

        # sys.argv with ingroup files whose sequences differ in length
        self.sizes = [
            "primerForge.py",
            "-i",
            *ParametersTest.SIZE_FNS,
            "-o",
            ParametersTest.OUT_FN,
        ]

        # sys.argv with a file that is neither genbank nor fasta
        self.badFmt = [
            "primerForge.py",
            "-i",
            *ParametersTest.IG_FNS_GB,
            ParametersTest.BAD_FMT_FN,
            "-o",
            ParametersTest.OUT_FN,
        ]

        # additional sys.argv to test
        self.help1 = ["primerForge.py", "-h"]
        self.help2 = ["primerForge.py", "--help"]
        self.vers1 = ["primerForge.py", "-v"]
        self.vers2 = ["primerForge.py", "--version"]
        self.debug1 = self.basic1 + ["--debug"]
        self.debug2 = self.basic2 + ["--debug"]
        self.debug3 = self.short1 + ["--debug"]
        self.debug4 = self.long1 + ["--debug"]

    def tearDown(self) -> None:
        """tear downs after each test"""
        # reset sys.stdout
        sys.stdout.close()
        sys.stdout = self.stdout

        # remove dummy files
        ParametersTest._removeDummyFiles()

    @staticmethod
    def _makeSeqRecord(length: int) -> SeqRecord:
        """builds a record whose sequence is the requested length

        Args:
            length (int): the length of the sequence; must be a multiple of four

        Returns:
            SeqRecord: a record that can be written as genbank or fasta
        """
        return SeqRecord(
            Seq("atcg" * (length // 4)),
            id="id",
            name="name",
            description="description",
            annotations={"molecule_type": "DNA"},
        )

    @staticmethod
    def _expectedFormat(fn: str) -> str:
        """determines the format a file should be detected as

        Args:
            fn (str): the filename

        Returns:
            str: the expected format
        """
        if fn.endswith(".gb"):
            return ParametersTest.FORMAT_GB

        return ParametersTest.FORMAT_FA

    @staticmethod
    def _makeDummyFiles() -> None:
        """creates the dummy files"""
        for fn in ParametersTest.IG_FNS_GB + ParametersTest.OG_FNS_GB:
            SeqIO.write(ParametersTest.TEST_SEQ, fn, ParametersTest.FORMAT_GB)

        for fn in ParametersTest.IG_FNS_FA + ParametersTest.OG_FNS_FA:
            SeqIO.write(ParametersTest.TEST_SEQ, fn, ParametersTest.FORMAT_FA)

        # make genbank files whose sequences differ in length
        for fn, length in zip(ParametersTest.SIZE_FNS, ParametersTest.SIZE_LENS):
            SeqIO.write(
                ParametersTest._makeSeqRecord(length), fn, ParametersTest.FORMAT_GB
            )

        # make a file that cannot be detected as genbank or fasta
        with open(ParametersTest.BAD_FMT_FN, "w") as fh:
            fh.write("this is not a sequence file\n")

    @staticmethod
    def _removeDummyFiles() -> None:
        """removes the dummy files"""
        for fn in (
            ParametersTest.IG_FNS_GB
            + ParametersTest.OG_FNS_GB
            + ParametersTest.IG_FNS_FA
            + ParametersTest.OG_FNS_FA
            + ParametersTest.SIZE_FNS
            + (ParametersTest.BAD_FMT_FN,)
        ):
            os.remove(fn)
        if os.path.exists(ParametersTest.DUMP_FN):
            os.remove(ParametersTest.DUMP_FN)

    def _checkDefaultValues(self, params: Parameters) -> None:
        """evaluates that params has the appropriate values when default

        Args:
            params (Parameters): Parameters object
        """
        # make sure the ingroup files were correctly parsed
        for fn in self.IG_FNS_GB:
            self.assertIn(fn, [x.fn.name for x in params.ingroup])

        # make sure the outfile is correct
        self.assertEqual(
            params.resultsFn.absolute(), self.dir.joinpath(ParametersTest.OUT_FN).absolute()
        )

        # check optional arguments match default values
        self.assertEqual(params.outgroup, ParametersTest.DEFAULT_ARGS.OUTGROUP)
        self.assertEqual(params.minLen, ParametersTest.DEFAULT_ARGS.MIN_LEN)
        self.assertEqual(params.maxLen, ParametersTest.DEFAULT_ARGS.MAX_LEN)
        self.assertEqual(params.minGc, ParametersTest.DEFAULT_ARGS.MIN_GC)
        self.assertEqual(params.maxGc, ParametersTest.DEFAULT_ARGS.MAX_GC)
        self.assertEqual(params.minTm, ParametersTest.DEFAULT_ARGS.MIN_TM)
        self.assertEqual(params.maxTm, ParametersTest.DEFAULT_ARGS.MAX_TM)
        self.assertEqual(params.minPcr, ParametersTest.DEFAULT_ARGS.MIN_PCR)
        self.assertEqual(params.maxPcr, ParametersTest.DEFAULT_ARGS.MAX_PCR)
        self.assertEqual(params.maxTmDiff, ParametersTest.DEFAULT_ARGS.MAX_TM_DIFF)
        self.assertEqual(params.numThreads, ParametersTest.DEFAULT_ARGS.NUM_THREADS)
        self.assertEqual(params.debug, ParametersTest.DEFAULT_ARGS.DEBUG)
        self.assertEqual(
            params.disallowedLens,
            range(DefaultArgs.MIN_PCR, DefaultArgs.MAX_PCR + 1),
        )

    def _checkGenomeFilesPresent(self, params: Parameters, frmt: str) -> None:
        """checks if the genome files are present

        Args:
            params (Parameters): Parameters object
            frmt (str): the format of the genome files
        """
        if frmt == ParametersTest.FORMAT_GB:
            expectedIngroupFns = ParametersTest.IG_FNS_GB
            expectedOutgroupFns = ParametersTest.OG_FNS_GB
        else:
            expectedIngroupFns = ParametersTest.IG_FNS_FA
            expectedOutgroupFns = ParametersTest.OG_FNS_FA

        for genome in params.ingroup:
            self.assertIn(genome.fn.name, expectedIngroupFns)

        # check for all outgroup files
        for genome in params.outgroup:
            self.assertIn(genome.fn.name, expectedOutgroupFns)

    def _checkCustomValues(self, params: Parameters, frmt: str) -> None:
        """evaluates that params has the appropriate values when custom values are specified

        Args:
            params (Parameters): a Parameters object
        """
        # check for all genome files
        self._checkGenomeFilesPresent(params, frmt)

        # parse values from the inputs
        minLen, maxLen = map(int, ParametersTest.PRIMER_LEN.split(","))
        minGc, maxGc = map(float, ParametersTest.GC_RANGE.split(","))
        minTm, maxTm = map(float, ParametersTest.TM_RANGE.split(","))
        minPcr, maxPcr = map(int, ParametersTest.PCR_SIZE.split(","))
        tmDiff = float(ParametersTest.TM_DIFF)
        threads = int(ParametersTest.THREADS)

        # get the disallowed outgroup sizes
        m, n = map(int, ParametersTest.BAD_SIZE.split(","))
        badSizes = range(m, n + 1)

        # make sure the parameters are correct
        self.assertEqual(
            params.resultsFn.absolute(), self.dir.joinpath(ParametersTest.OUT_FN).absolute()
        )
        self.assertEqual(params.minLen, minLen)
        self.assertEqual(params.maxLen, maxLen)
        self.assertEqual(params.minGc, minGc)
        self.assertEqual(params.maxGc, maxGc)
        self.assertEqual(params.minTm, minTm)
        self.assertEqual(params.maxTm, maxTm)
        self.assertEqual(params.minPcr, minPcr)
        self.assertEqual(params.maxPcr, maxPcr)
        self.assertEqual(params.maxTmDiff, tmDiff)
        self.assertEqual(params.numThreads, threads)
        self.assertEqual(params.disallowedLens, badSizes)
        self.assertEqual(params.debug, ParametersTest.DEFAULT_ARGS.DEBUG)

    def _dumpLoadTest(self, params: Parameters, obj) -> None:
        """evaluates if the dumpObj method is working

        Args:
            params (Parameters): Parameters object
            obj (Any): an object to dump
        """
        # dump the object to the dummy file
        params.dumpObj(obj, ParametersTest.DUMP_FN, "test")

        # load the dumped object and remove the file
        imported = params.loadObj(
            os.path.join(params.log.debugDir, ParametersTest.DUMP_FN)
        )
        os.remove(os.path.join(params.log.debugDir, ParametersTest.DUMP_FN))

        # make sure the original object matches the loaded one
        self.assertEqual(obj, imported)

    def testA_parseBasic1(self) -> None:
        """are args parsed with short flags and default values"""
        sys.argv = self.basic1
        params = _parseArgs()
        self._checkDefaultValues(params)

    def testB_parseBasic2(self) -> None:
        """are args parsed with long flags and default values"""
        sys.argv = self.basic2
        params = _parseArgs()
        self._checkDefaultValues(params)

    def testC_parseShort1(self) -> None:
        """are args parsed with short flags and custom values for genbank files"""
        sys.argv = self.short1
        params = _parseArgs()
        self._checkCustomValues(params, ParametersTest.FORMAT_GB)

    def testD_parseShort2(self) -> None:
        """are args parsed with short flags and custom values for fasta files"""
        sys.argv = self.short2
        params = _parseArgs()
        self._checkCustomValues(params, ParametersTest.FORMAT_FA)

    def testE_parseLong1(self) -> None:
        """are args parsed with long flags and custom args for genbank files"""
        sys.argv = self.long1
        params = _parseArgs()
        self._checkCustomValues(params, ParametersTest.FORMAT_GB)

    def testF_parseLong2(self) -> None:
        """are args parsed with long flags and custom args for genbank files"""
        sys.argv = self.long2
        params = _parseArgs()
        self._checkCustomValues(params, ParametersTest.FORMAT_FA)

    def testG_parseHelp(self) -> None:
        """Do we reach a sys.exit(0) when -h/--help is requested"""
        # check short flag
        with self.assertRaises(SystemExit) as e:
            sys.argv = self.help1
            params = _parseArgs()
            self.assertEqual(e.exception.code, 0)

        # check long flag
        with self.assertRaises(SystemExit) as e:
            sys.argv = self.help2
            params = _parseArgs()
            self.assertEqual(e.exception.code, 0)

    def testH_version(self) -> None:
        """Do we reach a sys.exit(0) when -v/--version is requested"""
        # check short flag
        with self.assertRaises(SystemExit) as e:
            sys.argv = self.vers1
            params = _parseArgs()
            self.assertEqual(e.exception.code, 0)

        # check long flag
        with self.assertRaises(SystemExit) as e:
            sys.argv = self.vers2
            params = _parseArgs()
            self.assertEqual(e.exception.code, 0)

    def testI_debug1(self) -> None:
        """checks if params.debug is True and has default values when in debug mode"""
        # check short flags with default args
        sys.argv = self.debug1
        params = _parseArgs()
        self.assertTrue(params.debug)
        params.debug = False
        self._checkDefaultValues(params)

        # check long flags with default args
        sys.argv = self.debug2
        params = _parseArgs()
        self.assertTrue(params.debug)
        params.debug = False
        self._checkDefaultValues(params)

        # check short flags with custom args
        sys.argv = self.debug3
        params = _parseArgs()
        self.assertTrue(params.debug)
        params.debug = False
        self._checkCustomValues(params, ParametersTest.FORMAT_GB)

        # check long flags with custom args
        sys.argv = self.debug4
        params = _parseArgs()
        self.assertTrue(params.debug)
        params.debug = False
        self._checkCustomValues(params, ParametersTest.FORMAT_GB)

    def testJ_debug2(self) -> None:
        """is the logger working"""
        # create the params object
        sys.argv = self.debug1
        params = _parseArgs()

        # replace the current Log object with one that references this directory
        params.log = Log(os.getcwd(), debug=True)

        # rename the log
        params.log.rename(ParametersTest.testJ_debug2.__name__)

        # verify that the log file exists
        self.assertTrue(os.path.exists(params.log.logFn))

        # make sure each writer works
        params.log.critical("")
        params.log.debug("")
        params.log.error("")
        params.log.info("")

        # remove the log file
        os.remove(params.log.logFn)

    def testK_dumpObjects(self) -> None:
        """evaluate Parameters.dumpObj"""
        # create a parameters object
        sys.argv = self.basic1
        params = _parseArgs()

        # initialize the log object
        params.log = Log(os.getcwd())
        params.log.rename(ParametersTest.testK_dumpObjects.__name__)

        # check sets
        obj = {1, 2, 3, 4, 5}
        self._dumpLoadTest(params, obj)

        # check lists
        obj = ["asdf", "jkl;"]
        self._dumpLoadTest(params, obj)

        # check dictionaries
        obj = {1: "one", 2: "two"}
        self._dumpLoadTest(params, obj)

    def testL_equality(self) -> None:
        """make sure equality overload works"""
        # create parameters objects for comparing
        sys.argv = self.basic1
        basic1 = _parseArgs()

        sys.argv = self.basic2
        basic2 = _parseArgs()

        sys.argv = self.short1
        short1 = _parseArgs()

        sys.argv = self.short2
        short2 = _parseArgs()

        sys.argv = self.long1
        long1 = _parseArgs()

        sys.argv = self.long2
        long2 = _parseArgs()

        sys.argv = self.debug1
        debugBasic1 = _parseArgs()

        sys.argv = self.debug2
        debugBasic2 = _parseArgs()

        sys.argv = self.debug3
        debugShort1 = _parseArgs()

        sys.argv = self.debug4
        debugLong1 = _parseArgs()

        sys.argv = self.long1[:-2] + ["--num_threads", "1"]
        oneThreadLong1 = _parseArgs()

        # the same objects should be equal
        self.assertEqual(basic1, basic1)
        self.assertEqual(basic2, basic2)
        self.assertEqual(short1, short1)
        self.assertEqual(short2, short2)
        self.assertEqual(long1, long1)
        self.assertEqual(long2, long2)
        self.assertEqual(debugBasic1, debugBasic1)
        self.assertEqual(debugBasic2, debugBasic2)
        self.assertEqual(debugShort1, debugShort1)
        self.assertEqual(debugLong1, debugLong1)
        self.assertEqual(oneThreadLong1, oneThreadLong1)

        # all basics and debug basics should be equal
        self.assertEqual(basic1, basic2)
        self.assertEqual(basic1, debugBasic1)
        self.assertEqual(basic1, debugBasic2)

        # long1, short1, debug variants, and thread variants should be equal
        self.assertEqual(long1, short1)
        self.assertEqual(long1, debugLong1)
        self.assertEqual(long1, debugShort1)
        self.assertEqual(long1, oneThreadLong1)

        # different files should not be equal
        self.assertNotEqual(long1, long2)
        self.assertNotEqual(short1, short2)
        self.assertNotEqual(long1, short2)
        self.assertNotEqual(short1, long2)

        # different input values with same files should not be equal
        self.assertNotEqual(basic1, short1)
        self.assertNotEqual(basic2, short1)
        self.assertNotEqual(basic1, long1)
        self.assertNotEqual(basic2, long1)

    def testM_sniffedFormats(self) -> None:
        """is the format of each genome detected from the file contents"""
        # constant
        FAIL_MSG = "wrong detected format for "

        # the format is no longer supplied on the command line, so each genome
        # must determine it from the contents of its own file
        for argv, frmt in (
            (self.short1, ParametersTest.FORMAT_GB),
            (self.short2, ParametersTest.FORMAT_FA),
        ):
            sys.argv = argv
            params = _parseArgs()

            for genome in params.ingroup + params.outgroup:
                self.assertEqual(
                    genome._format, frmt, f"{FAIL_MSG}{genome.fn.name}"
                )

    def testN_mixedFormats(self) -> None:
        """can genbank and fasta files be used together as inputs"""
        # constants
        FAIL_FMT = "wrong detected format for "
        FAIL_LEN = "wrong number of genomes parsed"

        # a mixture of formats is allowed because each file is sniffed separately
        sys.argv = self.mixed
        params = _parseArgs()

        # every input file must be present
        self.assertEqual(
            len(params.ingroup), len(ParametersTest.IG_FNS_MIX), FAIL_LEN
        )
        self.assertEqual(
            len(params.outgroup), len(ParametersTest.OG_FNS_MIX), FAIL_LEN
        )

        # each genome must be detected as the format matching its extension
        for genome in params.ingroup + params.outgroup:
            self.assertEqual(
                genome._format,
                ParametersTest._expectedFormat(genome.fn.name),
                f"{FAIL_FMT}{genome.fn.name}",
            )

        # the ingroup and outgroup files must be the ones that were requested
        self.assertEqual(
            sorted(x.fn.name for x in params.ingroup),
            sorted(ParametersTest.IG_FNS_MIX),
        )
        self.assertEqual(
            sorted(x.fn.name for x in params.outgroup),
            sorted(ParametersTest.OG_FNS_MIX),
        )

    def testO_sortIngroup(self) -> None:
        """are the ingroup genomes sorted from shortest to longest"""
        # constant
        FAIL_MSG = "ingroup genomes are not sorted by length"

        sys.argv = self.sizes
        params = _parseArgs()

        # the ingroup is sorted so that the smallest genome is processed first
        observed = [x.length for x in params.ingroup]
        self.assertEqual(observed, sorted(ParametersTest.SIZE_LENS), FAIL_MSG)

    def testP_preLoad(self) -> None:
        """are sequences streamed from file unless preloading was requested"""
        # constants
        FAIL_FLAG = "wrong pre_load value"
        FAIL_SEQS = "sequences were loaded into memory instead of streamed"

        # preloading is off by default
        sys.argv = self.basic1
        params = _parseArgs()
        self.assertEqual(
            params.pre_load, ParametersTest.DEFAULT_ARGS.PRE_LOAD, FAIL_FLAG
        )

        # without preloading the sequences are always streamed from file
        for genome in params.ingroup + params.outgroup:
            self.assertEqual(genome._seqs, [], f"{FAIL_SEQS}: {genome.fn.name}")

    def testQ_invalidFormat(self) -> None:
        """is an error raised when a file is neither genbank nor fasta"""
        # a file that cannot be sniffed is rejected after the genomes are built
        sys.argv = self.badFmt
        with self.assertRaises(InvalidFileFormat):
            _parseArgs()


if __name__ == "__main__":
    unittest.main()
