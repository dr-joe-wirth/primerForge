import pathlib, sys
sys.path.append(str(pathlib.Path(__file__).parent.parent.parent))

import os, unittest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters, Log


class ParametersTest(unittest.TestCase):
    # values for inputs
    IG_FNS_GB = ('itmp1.gb', 'itmp2.gb', 'itmp3.gb')
    IG_FNS_FA = ('itmp1.fa', 'itmp2.fa', 'itmp3.fa')
    OG_FNS_GB = ('otmp1.gb', 'otmp2.gb', 'otmp3.gb')
    OG_FNS_FA = ('otmp1.fa', 'otmp2.fa', 'otmp3.fa')
    IG_PATTERN_GB = "itmp*gb"
    IG_PATTERN_FA = "itmp*fa"
    OG_PATTERN_GB = "otmp*gb"
    OG_PATTERN_FA = "otmp*fa"
    OUT_FN = "result.testfile"
    BAD_SIZE = "64,160"
    FORMAT_GB = "genbank"
    FORMAT_FA = "fasta"
    PRIMER_LEN = "14,25"
    GC_RANGE = "42,64"
    TM_RANGE = "68,72"
    PCR_SIZE = "80,128"
    TM_DIFF = "4.2"
    THREADS = "8"
    DUMP_FN = "tmp.p"
    VERSION = 'parameters_test'
    AUTHOR = 'null'
    TEST_SEQ = SeqRecord(Seq('atcg'),
                         id='id',
                         name='name',
                         description='description',
                         annotations={'molecule_type':'DNA'})
    
    
    def setUp(self) -> None:
        """set up for each test
        """
        # silence prints
        self.stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        
        # make the dummy files needed for testing
        ParametersTest._makeDummyFiles()
        
        # save the directory
        self.dir = os.getcwd()
        
        # sys.argv for default values using short args
        self.basic1 = ['primerForge.py',
                       '-i', ParametersTest.IG_PATTERN_GB,
                       '-o', ParametersTest.OUT_FN]
        
        # sys.argv for default values using long args
        self.basic2 = ['primerForge.py',
                       '--ingroup', ParametersTest.IG_PATTERN_GB,
                       '--out', ParametersTest.OUT_FN]
        
        # sys.argv with custom values using short args for genbank file format
        self.short1  = ['primerForge.py',
                        '-i', ParametersTest.IG_PATTERN_GB,
                        '-o', ParametersTest.OUT_FN,
                        '-u', ParametersTest.OG_PATTERN_GB,
                        '-b', ParametersTest.BAD_SIZE,
                        '-f', ParametersTest.FORMAT_GB,
                        '-p', ParametersTest.PRIMER_LEN,
                        '-g', ParametersTest.GC_RANGE,
                        '-t', ParametersTest.TM_RANGE,
                        '-r', ParametersTest.PCR_SIZE,
                        '-d', ParametersTest.TM_DIFF,
                        '-n', ParametersTest.THREADS]
        
        # sys.argv with custom values using short args for fasta file format
        self.short2  = ['primerForge.py',
                        '-i', ParametersTest.IG_PATTERN_FA,
                        '-o', ParametersTest.OUT_FN,
                        '-u', ParametersTest.OG_PATTERN_FA,
                        '-b', ParametersTest.BAD_SIZE,
                        '-f', ParametersTest.FORMAT_FA,
                        '-p', ParametersTest.PRIMER_LEN,
                        '-g', ParametersTest.GC_RANGE,
                        '-t', ParametersTest.TM_RANGE,
                        '-r', ParametersTest.PCR_SIZE,
                        '-d', ParametersTest.TM_DIFF,
                        '-n', ParametersTest.THREADS]
        
        # sys.argv with custom values using long args for genbank file format
        self.long1   =  ['primerForge.py',
                        '--ingroup', ParametersTest.IG_PATTERN_GB,
                        '--out', ParametersTest.OUT_FN,
                        '--outgroup', ParametersTest.OG_PATTERN_GB,
                        '--bad_sizes', ParametersTest.BAD_SIZE,
                        '--format', ParametersTest.FORMAT_GB,
                        '--primer_len', ParametersTest.PRIMER_LEN,
                        '--gc_range', ParametersTest.GC_RANGE,
                        '--tm_range', ParametersTest.TM_RANGE,
                        '--pcr_prod', ParametersTest.PCR_SIZE,
                        '--tm_diff', ParametersTest.TM_DIFF,
                        '--num_threads', ParametersTest.THREADS]
        
        # sys.argv with custom values using long args for fasta file format
        self.long2   =  ['primerForge.py',
                        '--ingroup', ParametersTest.IG_PATTERN_FA,
                        '--out', ParametersTest.OUT_FN,
                        '--outgroup', ParametersTest.OG_PATTERN_FA,
                        '--bad_sizes', ParametersTest.BAD_SIZE,
                        '--format', ParametersTest.FORMAT_FA,
                        '--primer_len', ParametersTest.PRIMER_LEN,
                        '--gc_range', ParametersTest.GC_RANGE,
                        '--tm_range', ParametersTest.TM_RANGE,
                        '--pcr_prod', ParametersTest.PCR_SIZE,
                        '--tm_diff', ParametersTest.TM_DIFF,
                        '--num_threads', ParametersTest.THREADS]
        
        # additional sys.argv to test
        self.help1  = ['primerForge.py', '-h']
        self.help2  = ['primerForge.py', '--help']
        self.vers1  = ['primerForge.py', '-v']
        self.vers2  = ['primerForge.py', '--version']
        self.debug1 = self.basic1 + ['--debug']
        self.debug2 = self.basic2 + ['--debug']
        self.debug3 = self.short1 + ['--debug']
        self.debug4 = self.long1  + ['--debug']
    
    def tearDown(self) -> None:
        """tear downs after each test
        """
        # reset sys.stdout
        sys.stdout.close()
        sys.stdout = self.stdout
        
        # remove dummy files
        ParametersTest._removeDummyFiles()

    def _makeDummyFiles() -> None:
        """creates the dummy files
        """
        for fn in ParametersTest.IG_FNS_GB + ParametersTest.OG_FNS_GB:
            SeqIO.write(ParametersTest.TEST_SEQ, fn, ParametersTest.FORMAT_GB)
        
        for fn in ParametersTest.IG_FNS_FA + ParametersTest.OG_FNS_FA:
            SeqIO.write(ParametersTest.TEST_SEQ, fn, ParametersTest.FORMAT_FA)

    def _removeDummyFiles() -> None:
        """removes the dummy files
        """
        for fn in ParametersTest.IG_FNS_GB + ParametersTest.OG_FNS_GB + ParametersTest.IG_FNS_FA + ParametersTest.OG_FNS_FA:
            os.remove(fn)
        if os.path.exists(ParametersTest.DUMP_FN):
            os.remove(ParametersTest.DUMP_FN)

    def _checkDefaultValues(self, params:Parameters) -> None:
        """evaluates that params has the appropriate values when default

        Args:
            params (Parameters): Parameters object
        """
        # make sure the ingroup files were correctly parsed
        for fn in self.IG_FNS_GB:
            self.assertIn(fn, params.ingroupFns)
        
        # make sure the outfile is correct
        self.assertEqual(params.resultsFn, os.path.join(self.dir, ParametersTest.OUT_FN))
        
        # check optional arguments match default values
        self.assertEqual(params.outgroupFns, Parameters._DEF_OUTGROUP)
        self.assertEqual(params.format, Parameters._DEF_FRMT)
        self.assertEqual(params.minLen, Parameters._DEF_MIN_LEN)
        self.assertEqual(params.maxLen, Parameters._DEF_MAX_LEN)
        self.assertEqual(params.minGc, Parameters._DEF_MIN_GC)
        self.assertEqual(params.maxGc, Parameters._DEF_MAX_GC)
        self.assertEqual(params.minTm, Parameters._DEF_MIN_TM)
        self.assertEqual(params.maxTm, Parameters._DEF_MAX_TM)
        self.assertEqual(params.minPcr, Parameters._DEF_MIN_PCR)
        self.assertEqual(params.maxPcr, Parameters._DEF_MAX_PCR)
        self.assertEqual(params.maxTmDiff, Parameters._DEF_MAX_TM_DIFF)
        self.assertEqual(params.numThreads, Parameters._DEF_NUM_THREADS)
        self.assertEqual(params.helpRequested, Parameters._DEF_HELP)
        self.assertEqual(params.debug, Parameters._DEF_DEBUG)
        self.assertEqual(params.disallowedLens, range(Parameters._DEF_MIN_PCR, Parameters._DEF_MAX_PCR+1))

    def _checkGenomeFilesPresent(self, params:Parameters, frmt:str) -> None:
        """checks if the genome files are present

        Args:
            params (Parameters): Parameters object
            frmt (str): the format of the genome files
        """
        if frmt == ParametersTest.FORMAT_GB:
            expectedIngroupFns  = ParametersTest.IG_FNS_GB
            expectedOutgroupFns = ParametersTest.OG_FNS_GB
        else:
            expectedIngroupFns  = ParametersTest.IG_FNS_FA
            expectedOutgroupFns = ParametersTest.OG_FNS_FA
            
        for fn in params.ingroupFns:
            self.assertIn(fn, expectedIngroupFns)
        
        # check for all outgroup files
        for fn in params.outgroupFns:
            self.assertIn(fn, expectedOutgroupFns)
    
    def _checkCustomValues(self, params:Parameters, frmt:str) -> None:
        """evaluates that params has the appropriate values when custom values are specified

        Args:
            params (Parameters): a Parameters object
        """
        # check for all genome files
        self._checkGenomeFilesPresent(params, frmt)
        
        # parse values from the inputs
        minLen,maxLen = map(int, ParametersTest.PRIMER_LEN.split(','))
        minGc,maxGc = map(float, ParametersTest.GC_RANGE.split(","))
        minTm,maxTm = map(float, ParametersTest.TM_RANGE.split(','))
        minPcr,maxPcr = map(int, ParametersTest.PCR_SIZE.split(','))
        tmDiff = float(ParametersTest.TM_DIFF)
        threads = int(ParametersTest.THREADS)
        
        # get the disallowed outgroup sizes
        m,n = map(int, ParametersTest.BAD_SIZE.split(','))
        badSizes = range(m,n+1)
        
        # make sure the parameters are correct
        self.assertEqual(params.resultsFn, os.path.join(self.dir, ParametersTest.OUT_FN))
        self.assertEqual(params.format, frmt)
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
        self.assertEqual(params.helpRequested, Parameters._DEF_HELP)
        self.assertEqual(params.debug, Parameters._DEF_DEBUG)
    
    def _dumpLoadTest(self, params:Parameters, obj) -> None:
        """evaluates if the dumpObj method is working

        Args:
            params (Parameters): Parameters object
            obj (Any): an object to dump
        """
        # dump the object to the dummy file
        params.dumpObj(obj, ParametersTest.DUMP_FN, 'test')
        
        # load the dumped object and remove the file
        imported = params.loadObj(os.path.join(params.log.debugDir, ParametersTest.DUMP_FN))
        os.remove(os.path.join(params.log.debugDir, ParametersTest.DUMP_FN))
        
        # make sure the original object matches the loaded one
        self.assertEqual(obj, imported)
    
    def testA_parseBasic1(self) -> None:
        """are args parsed with short flags and default values
        """
        sys.argv = self.basic1
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        self._checkDefaultValues(params)
    
    def testB_parseBasic2(self) -> None:
        """are args parsed with long flags and default values
        """
        sys.argv = self.basic2
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        self._checkDefaultValues(params)

    def testC_parseShort1(self) -> None:
        """are args parsed with short flags and custom values for genbank files
        """
        sys.argv = self.short1
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        self._checkCustomValues(params, ParametersTest.FORMAT_GB)
    
    def testD_parseShort2(self) -> None:
        """are args parsed with short flags and custom values for fasta files
        """
        sys.argv = self.short2
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        self._checkCustomValues(params, ParametersTest.FORMAT_FA)
    
    def testE_parseLong1(self) -> None:
        """are args parsed with long flags and custom args for genbank files
        """
        sys.argv = self.long1
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        self._checkCustomValues(params, ParametersTest.FORMAT_GB)
    
    def testF_parseLong2(self) -> None:
        """are args parsed with long flags and custom args for genbank files
        """
        sys.argv = self.long2
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        self._checkCustomValues(params, ParametersTest.FORMAT_FA)
    
    def testE_parseHelp(self) -> None:
        """is params.helpRequested True when help is requested
        """
        # check short flag
        sys.argv = self.help1
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        self.assertTrue(params.helpRequested)
        
        # check long flag
        sys.argv = self.help2
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        self.assertTrue(params.helpRequested)
    
    def testF_version(self) -> None:
        """is params.helpRequested True when version is requested
        """
        # check short flag
        sys.argv = self.vers1
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        self.assertTrue(params.helpRequested)
        
        # check long flag
        sys.argv = self.vers2
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        self.assertTrue(params.helpRequested)
    
    def testG_debug1(self) -> None:
        """checks if params.debug is True and has default values when in debug mode
        """
        # check short flags with default args
        sys.argv = self.debug1
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        self.assertTrue(params.debug)
        params.debug = False
        self._checkDefaultValues(params)
        
        # check long flags with default args
        sys.argv = self.debug2
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        self.assertTrue(params.debug)
        params.debug = False
        self._checkDefaultValues(params)
        
        # check short flags with custom args
        sys.argv = self.debug3
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        self.assertTrue(params.debug)
        params.debug = False
        self._checkCustomValues(params, Parameters._DEF_FRMT)
        
        # check long flags with custom args
        sys.argv = self.debug4
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        self.assertTrue(params.debug)
        params.debug = False
        self._checkCustomValues(params, Parameters._DEF_FRMT)
    
    def testH_debug2(self) -> None:
        """is the logger working
        """
        # create the params object
        sys.argv = self.debug1
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False, initializeLog=False)
        
        # replace the current Log object with one that references this directory
        params.log = Log(os.getcwd(), debug=True)
        
        # rename the log
        params.log.rename(ParametersTest.testH_debug2.__name__)
        
        # verify that the log file exists
        self.assertTrue(os.path.exists(params.log.logFn))
        
        # make sure each writer works
        params.log.critical('')
        params.log.debug('')
        params.log.error('')
        params.log.info('')
        
        # remove the log file
        os.remove(params.log.logFn)
    
    def testI_dumpObjects(self) -> None:
        """evaluate Parameters.dumpObj
        """
        # create a parameters object
        sys.argv = self.basic1
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False, initializeLog=False)
        
        # initialize the log object
        params.log = Log(os.getcwd())
        params.log.rename(ParametersTest.testI_dumpObjects.__name__)
        
        # check sets
        obj  = {1,2,3,4,5}
        self._dumpLoadTest(params, obj)
        
        # check lists
        obj  = ['asdf', 'jkl;']
        self._dumpLoadTest(params, obj)
        
        # check dictionaries
        obj  = {1:'one', 2:'two'}
        self._dumpLoadTest(params, obj)

    def testJ_equality(self) -> None:
        """make sure equality overload works
        """
        # create parameters objects for comparing
        sys.argv = self.basic1
        basic1 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        sys.argv = self.basic2
        basic2 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        sys.argv = self.short1
        short1 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        sys.argv = self.short2
        short2 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        sys.argv = self.long1
        long1 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        sys.argv = self.long2
        long2 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        sys.argv = self.help1
        help1 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        sys.argv = self.help2
        help2 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        sys.argv = self.vers1
        vers1 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        sys.argv = self.vers2
        vers2 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        sys.argv = self.debug1
        debugBasic1 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        sys.argv = self.debug2
        debugBasic2 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        sys.argv = self.debug3
        debugShort1 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        sys.argv = self.debug4
        debugLong1 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        sys.argv = self.long1[:-2] + ['--num_threads', '1']
        oneThreadLong1 = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION, makeWd=False)
        
        # the same objects should be equal
        self.assertEqual(basic1, basic1)
        self.assertEqual(basic2, basic2)
        self.assertEqual(short1, short1)
        self.assertEqual(short2, short2)
        self.assertEqual(long1, long1)
        self.assertEqual(long2, long2)
        self.assertEqual(help1, help1)
        self.assertEqual(help2, help2)
        self.assertEqual(vers1, vers1)
        self.assertEqual(vers2, vers2)
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
        
        # help and vers should be equal
        self.assertEqual(help1, help2)
        self.assertEqual(help1, vers1)
        self.assertEqual(help2, vers2)
        
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


if __name__ == "__main__":
    unittest.main()
