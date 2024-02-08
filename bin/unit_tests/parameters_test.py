import os, pickle, sys, unittest
from bin.Parameters import Parameters, Log


class ParametersTest(unittest.TestCase):
    # values for inputs
    IG_FNS = ('itmp1.testfile', 'itmp2.testfile', 'itmp3.testfile')
    OG_FNS = ('otmp4.testfile', 'otmp5.testfile', 'otmp6.testfile')
    IG_PATTERN = "itmp*testfile"
    OG_PATTERN = "otmp*testfile"
    OUT_FN = "result.testfile"
    BAD_SIZE = "64,160"
    FORMAT = "genbank"
    PRIMER_LEN = "14,25"
    GC_RANGE = "42,64"
    TM_RANGE = "68,72"
    PCR_SIZE = "80,128"
    TM_DIFF = "4.2"
    THREADS = "8"
    DUMP_FN = "tmp.p"
    VERSION = 'null'
    AUTHOR = 'null'
    
    
    def setUp(self) -> None:
        """set up for each test
        """
        # silence prints
        self.stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        
        # make the dummy files needed for testing
        ParametersTest._makeDummyFiles()
        
        # sys.argv for default values using short args
        self.basic1 = ['primerForge.py',
                       '-i', ParametersTest.IG_PATTERN,
                       '-o', ParametersTest.OUT_FN]
        
        # sys.argv for default values using long args
        self.basic2 = ['primerForge.py',
                       '--ingroup', ParametersTest.IG_PATTERN,
                       '--out', ParametersTest.OUT_FN]
        
        # sys.argv with custom values using short args
        self.short  = ['primerForge.py',
                       '-i', ParametersTest.IG_PATTERN,
                       '-o', ParametersTest.OUT_FN,
                       '-u', ParametersTest.OG_PATTERN,
                       '-b', ParametersTest.BAD_SIZE,
                       '-f', ParametersTest.FORMAT,
                       '-p', ParametersTest.PRIMER_LEN,
                       '-g', ParametersTest.GC_RANGE,
                       '-t', ParametersTest.TM_RANGE,
                       '-r', ParametersTest.PCR_SIZE,
                       '-d', ParametersTest.TM_DIFF,
                       '-n', ParametersTest.THREADS]

        # sys.argv with custom values using long args
        self.long   =  ['primerForge.py',
                       '--ingroup', ParametersTest.IG_PATTERN,
                       '--out', ParametersTest.OUT_FN,
                       '--outgroup', ParametersTest.OG_PATTERN,
                       '--bad_sizes', ParametersTest.BAD_SIZE,
                       '--format', ParametersTest.FORMAT,
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
        self.debug3 = self.short + ['--debug']
        self.debug4 = self.long + ['--debug']
    
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
        for fn in ParametersTest.IG_FNS + ParametersTest.OG_FNS:
            fh = open(fn, 'w')
            fh.close()
    
    def _removeDummyFiles() -> None:
        """removes the dummy files
        """
        for fn in ParametersTest.IG_FNS + ParametersTest.OG_FNS:
            os.remove(fn)
        if os.path.exists(ParametersTest.DUMP_FN):
            os.remove(ParametersTest.DUMP_FN)

    def _checkDefaultValues(self, params:Parameters) -> None:
        """evaluates that params has the appropriate values when default

        Args:
            params (Parameters): Parameters object
        """
        # make sure the ingroup files were correctly parsed
        for fn in self.IG_FNS:
            self.assertIn(fn, params.ingroupFns)
        
        # make sure the outfile is correct
        self.assertEqual(params.outFn, ParametersTest.OUT_FN)
        
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
    
    def _checkCustomValues(self, params:Parameters) -> None:
        """evaluates that params has the appropriate values when custom values are specified

        Args:
            params (Parameters): a Parameters object
        """
        # check for all ingroup files
        for fn in params.ingroupFns:
            self.assertIn(fn, ParametersTest.IG_FNS)
        
        # check for all outgroup files
        for fn in params.outgroupFns:
            self.assertIn(fn, ParametersTest.OG_FNS)
        
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
        self.assertEqual(params.outFn, ParametersTest.OUT_FN)
        self.assertEqual(params.format, ParametersTest.FORMAT)
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
        with open(os.path.join(params.log.debugDir, ParametersTest.DUMP_FN), 'rb') as fh:
            imported = pickle.load(fh)
        os.remove(os.path.join(params.log.debugDir, ParametersTest.DUMP_FN))
        
        # make sure the original object matches the loaded one
        self.assertEqual(obj, imported)
    
    def testA_parseBasic1(self) -> None:
        """are args parsed with short flags and default values
        """
        sys.argv = self.basic1
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION)
        self._checkDefaultValues(params)
    
    def testB_parseBasic2(self) -> None:
        """are args parsed with long flags and default values
        """
        sys.argv = self.basic2
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION)
        self._checkDefaultValues(params)

    def testC_parseShort(self) -> None:
        """are args parsed with short flags and custom values
        """
        sys.argv = self.short
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION)
        self._checkCustomValues(params)
    
    def testD_parseLong(self) -> None:
        """are args parsed with long flags and custom args
        """
        sys.argv = self.long
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION)
        self._checkCustomValues(params)
    
    def testE_parseHelp(self) -> None:
        """is params.helpRequested True when help is requested
        """
        # check short flag
        sys.argv = self.help1
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION)
        self.assertTrue(params.helpRequested)
        
        # check long flag
        sys.argv = self.help2
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION)
        self.assertTrue(params.helpRequested)
    
    def testF_version(self) -> None:
        """is params.helpRequested True when version is requested
        """
        # check short flag
        sys.argv = self.vers1
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION)
        self.assertTrue(params.helpRequested)
        
        # check long flag
        sys.argv = self.vers2
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION)
        self.assertTrue(params.helpRequested)
    
    def testG_debug1(self) -> None:
        """checks if params.debug is True and has default values when in debug mode
        """
        # check short flags with default args
        sys.argv = self.debug1
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION)
        self.assertTrue(params.debug)
        params.debug = False
        self._checkDefaultValues(params)
        
        # check long flags with default args
        sys.argv = self.debug2
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION)
        self.assertTrue(params.debug)
        params.debug = False
        self._checkDefaultValues(params)
        
        # check short flags with custom args
        sys.argv = self.debug3
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION)
        self.assertTrue(params.debug)
        params.debug = False
        self._checkCustomValues(params)
        
        # check long flags with custom args
        sys.argv = self.debug4
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION)
        self.assertTrue(params.debug)
        params.debug = False
        self._checkCustomValues(params)
    
    def testH_debug2(self) -> None:
        """is the logger working
        """
        # create the params object
        sys.argv = self.debug1
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION)
        
        # replace the current Log object with one that references this directory
        params.log = Log(os.getcwd())
        
        # start up the log
        params.log.initialize(ParametersTest.testH_debug2.__name__)
        
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
        params = Parameters(ParametersTest.AUTHOR, ParametersTest.VERSION)
        
        # initialize the log object
        params.log = Log(os.getcwd())
        params.log.initialize('tmp')
        
        # check sets
        obj  = {1,2,3,4,5}
        self._dumpLoadTest(params, obj)
        
        # check lists
        obj  = ['asdf', 'jkl;']
        self._dumpLoadTest(params, obj)
        
        # check dictionaries
        obj  = {1:'one', 2:'two'}
        self._dumpLoadTest(params, obj)
