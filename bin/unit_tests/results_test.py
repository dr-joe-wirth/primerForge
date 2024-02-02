from __future__ import annotations
from bin.Clock import Clock, _printDone, _printStart
from bin.getCandidateKmers import _kmpSearch
import gzip, os, subprocess, sys, unittest
from bin.Parameters import Parameters
from Bio.SeqUtils import MeltingTemp
from bin.main import _main
from Bio.Seq import Seq
from bin.Log import Log

class Result():
    """class to save results for easy lookup
    """
    def __init__(self, fwdTm:float, fwdGc:float, revTm:float, revGc:float, additional:dict[str]) -> Result:
        self.fwdTm:float = fwdTm
        self.fwdGc:float = fwdGc
        self.revTm:float = revTm
        self.revGc:float = revGc
        self.additional:dict[str] = additional
    
    def __repr__(self) -> str:
        return str(vars(self))

class ResultsTest(unittest.TestCase):
    """class for testing the results of primerForge
    """
    # constants
    NUM_THREADS = 24
    TEST_DIR = os.path.join(os.getcwd(), "test_dir")
    RESULT_FN = os.path.join(TEST_DIR, "results.tsv")
    GENOMES = {'i1.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/650/295/GCF_001650295.1_ASM165029v1/GCF_001650295.1_ASM165029v1_genomic.gbff.gz',
               'i2.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/727/175/GCF_008727175.1_ASM872717v1/GCF_008727175.1_ASM872717v1_genomic.gbff.gz',
               'i3.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/208/865/GCF_002208865.2_ASM220886v2/GCF_002208865.2_ASM220886v2_genomic.gbff.gz',
               'o1.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz',
               'o2.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gbff.gz'}
    
    @classmethod
    def setUpClass(cls) -> None:
        """sets up the class for all tests
        """
        # initialize clocks
        total = Clock()
        clock = Clock()
        
        # get the parameters
        cls.params:Parameters = ResultsTest._getParameters()
        cls.params.log.initialize(ResultsTest.setUpClass.__name__)
        
        # load existing results if present
        if os.path.exists(ResultsTest.RESULT_FN):
            print(f'running tests on existing results file: {ResultsTest.RESULT_FN}')
        
        # otherwise run primerForge
        else:
            _printStart(total, 'getting results for testing', end=' ...\n')
            
            _printStart(clock, 'downloading genomes')
            ResultsTest._downloadTestData()
            _printDone(clock)
            
            _printStart(clock, 'running primerForge', end=' ...\n')
            _main(cls.params)
            _printDone(clock)
            print()
            _printDone(total)
        
        # load the results into memory
        _printStart(clock, 'reading results into memory')
        cls.results:dict[tuple[Seq,Seq],Result] = ResultsTest._parseResultsFile(ResultsTest.RESULT_FN)
        _printDone(clock)
    
    def _downloadOneGenome(ftp:str, fn:str) -> None:
        """downloads a genome using wget

        Args:
            ftp (str): the ftp path of the file to download
            fn (str): the file to save the genome

        Raises:
            RuntimeError: failed to download a genome
        """
        # constant
        CMD = ('wget', '-q', '-O')
        
        # get the gzip filename
        gz = fn + ".gz"
        
        # build the command to run
        cmd = list(CMD)
        cmd.append(gz)
        cmd.append(ftp)
        
        # download the genome
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError:
            raise RuntimeError(f"download failed for {os.path.basename(fn)}")
        
        # extract the gzipped file
        with gzip.open(gz, 'rt') as zh:
            with open(fn, 'w') as fh:
                for line in zh:
                    fh.write(line)
        
        # remove the gzipped file
        os.remove(gz)
    
    def _downloadTestData() -> None:
        """downloads all the genomes for testing
        """
        # make the test directory if it doesn't exist
        if not os.path.isdir(ResultsTest.TEST_DIR):
            os.mkdir(ResultsTest.TEST_DIR)
        
        # download each genome from ncbi
        for fn in ResultsTest.GENOMES.keys():
            # get the ftp path and the filename
            ftp = ResultsTest.GENOMES[fn]
            fn = os.path.join(ResultsTest.TEST_DIR, fn)
            
            # download the genome
            ResultsTest._downloadOneGenome(ftp, fn)
    
    def _getParameters() -> Parameters:
        """creates a Parameters object for testing

        Returns:
            Parameters: a Parameters object
        """
        # make the parameters object
        sys.argv = ['primerForge.py',
                    '-i', os.path.join(ResultsTest.TEST_DIR, "i[123].gbff"),
                    '-u', os.path.join(ResultsTest.TEST_DIR, "o[12].gbff"),
                    '-r', '70,120',
                    '-b', '45,150',
                    '-n', ResultsTest.NUM_THREADS,
                    '-o', ResultsTest.RESULT_FN,
                    '--debug']
        params = Parameters('', '')
        
        # use a different logger pointing at the test directory
        params.log = Log(ResultsTest.TEST_DIR)
        
        return params

    def _parseResultsFile(fn:str) -> dict[tuple[Seq,Seq],Result]:
        """parses the results file produced by running primerForge

        Args:
            fn (str): the filename of the results file

        Returns:
            dict[tuple[Seq,Seq],Result]: key=primer pair; value=Result object
        """
        # constants
        SEP_1 = "\t"
        SEP_2 = ","
        FWD_SEQ = 0
        FWD_TM  = 1
        FWD_GC  = 2
        REV_SEQ = 3
        REV_TM  = 4
        REV_GC  = 5
        FIRST   = 6
        LEN_SUFFIX = "_length"
        CON_SUFFIX = "_contig"
        
        # initialize variables
        out = dict()
        key = dict()
        header = True
        
        # go through each line in the file
        with open(fn, 'r') as fh:
            for line in fh:
                # conver the line to a row of columns
                line = line.rstrip()
                row = line.split(SEP_1)
                
                # find out the header names for the genome-specific names
                if header:
                    for idx in range(FIRST, len(row)):
                        key[idx] = row[idx]
                    header = False
                
                # extract data if not a  header
                else:
                    # get known expected columns
                    fwd = Seq(row[FWD_SEQ])
                    ftm = float(row[FWD_TM])
                    fgc = float(row[FWD_GC])
                    rev = Seq(row[REV_SEQ])
                    rtm = float(row[REV_TM])
                    rgc = float(row[REV_GC])
                    
                    # import genome-specific data
                    additional = dict()
                    for idx in key.keys():
                        # if the `_length` string is in the column name
                        if LEN_SUFFIX == key[idx][-len(LEN_SUFFIX):]:
                            # get the name of the sequence
                            name = key[idx][:-len(LEN_SUFFIX)]
                            
                            # initialize a dictionary for the name if one doesn't exist
                            additional[name] = additional.get(name, dict())
                            
                            # then save a list of product sizes
                            additional[name]['length'] = list(map(int, row[idx].split(SEP_2)))
                        
                        # otherwise it is a contig name
                        else:
                            # get the name of the sequence
                            name = key[idx][:-len(CON_SUFFIX)]
                            
                            # initialize a dictionary for the name if one deosn't exist
                            additional[name] = additional.get(name, dict())
                            
                            # then save the a list of contig names
                            additional[name]['contig'] = row[idx].split(SEP_2)
                    
                    # store the result for this primer pair
                    out[(fwd,rev)] = Result(ftm,fgc,rtm,rgc,additional)
        
        return out

    def _getGc(seq:Seq) -> float:
        """gets the percent G+C for the provided sequence

        Args:
            seq (Seq): the sequence to calculate

        Returns:
            float: the percent G+C
        """
        seq = str(seq.upper())
        numGc = len(seq.replace('A', '').replace('T', ''))
        return numGc / len(seq) * 100
    
    def _noLongRepeats(seq:Seq) -> bool:
        """determines if a sequence does not have long homopolymers

        Args:
            seq (Seq): a sequence to evaluate

        Returns:
            bool: indicates if the sequences does not have long repeats
        """
        # constants
        MAX_LEN = 4
        REPEATS = map(Seq, ("A"*MAX_LEN, "T"*MAX_LEN, "C"*MAX_LEN, "G"*MAX_LEN))

        # check for each repeat in the primer
        for repeat in REPEATS:
            if _kmpSearch(seq, repeat)[0]:
                return False
        return True
    
    def _getBindingSites(self, seq:Seq) -> None:
        pass
    
    def _isRepeatedInIngroup(self, fbind, rbind) -> None:
        pass
    
    def _isProductLengthCorrect(self, fbind, rbind) -> None:
        pass
    
    def _onSeparateStrands(self, fbind, rbind) -> None:
        pass
    
    def _noDisallowedOutgroupProducts(self, fbind, rbind) -> None:
        pass
    
    def testA_checkTm(self) -> None:
        """does the Tm match the expected value
        """
        # for each pair, make sure the saved Tm matches the expected value
        for fwd,rev in self.results.keys():
            self.assertEqual(self.results[(fwd,rev)].fwdTm, round(MeltingTemp.Tm_Wallace(fwd), 1))
            self.assertEqual(self.results[(fwd,rev)].revTm, round(MeltingTemp.Tm_Wallace(rev), 1))
    
    def testB_tmWithinRange(self) -> None:
        """is the Tm within the specified range
        """
        # for each pair, make sure the Tms are within the specified ranges
        for result in self.results.values():
            self.assertGreaterEqual(result.fwdTm, self.params.minTm)
            self.assertGreaterEqual(result.revTm, self.params.minTm)
            self.assertLessEqual(result.fwdTm, self.params.maxTm)
            self.assertLessEqual(result.revTm, self.params.maxTm)
    
    def testC_tmDiffWithinRange(self) -> None:
        """is the Tm difference below the specified threshold
        """
        # for each pair, make sure the Tm difference is below the threshold
        for result in self.results.values():
            tmDiff = abs(result.fwdTm - result.revTm)
            self.assertLessEqual(tmDiff, self.params.maxTmDiff)
    
    def testD_checkGcPercent(self) -> None:
        """is the G+C percentage the expected value
        """
        # for each pair, make sure the GC percent is accurate
        for fwd,rev in self.results.keys():
            self.assertEqual(self.results[(fwd,rev)].fwdGc, round(ResultsTest._getGc(fwd), 1))
            self.assertEqual(self.results[(fwd,rev)].revGc, round(ResultsTest._getGc(rev), 1))
    
    def testE_gcPercentInRange(self) -> None:
        """is the G+C percentage within the specified range
        """
        # for each pair, make sure the GC percent is within the specified range
        for result in self.results.values():
            self.assertGreaterEqual(result.fwdGc, self.params.minGc)
            self.assertGreaterEqual(result.revGc, self.params.minGc)
            self.assertLessEqual(result.fwdGc, self.params.maxGc)
            self.assertLessEqual(result.revGc, self.params.maxGc)
    
    def testF_threePrimeIsGc(self) -> None:
        """is the 3' end of the primer a G or C
        """
        # constant
        GC = ("G", "C")
        
        # for each pair, make sure the 3' end is a G or a C
        for fwd,rev in self.results.keys():
            self.assertIn(fwd[-1], GC)
            self.assertIn(rev[-1], GC)

    def testG_noHomoPolymers(self) -> None:
        """does the primer contain homopolymers
        """
        # for each pair, make sure there are no homopolymers
        for fwd,rev in self.results.keys():
            self.assertTrue(ResultsTest._noLongRepeats(fwd))
            self.assertTrue(ResultsTest._noLongRepeats(rev))

    def testH_ingroupHasOneProduct(self) -> None:
        """do all pairs produce one product size in the ingroup
        """
        for pair in self.results.keys():
            pass

    def testH_sequenceTests(self) -> None:
        """evaluate several additional tests
        """
        for fwd,rev in self.results.keys():
            fwdIngroupBind,fwdOutgroupBind = self._getBindingSites(fwd)
            revIngroupBind,revOutgroupBind = self._getBindingSites(rev)
            self._isRepeatedInIngroup(fwdIngroupBind,revIngroupBind)
            self._isProductLengthCorrect(fwdIngroupBind, revIngroupBind)
            self._onSeparateStrands(fwdIngroupBind, revIngroupBind)
            self._noDisallowedOutgroupProducts(fwdOutgroupBind, revOutgroupBind)
