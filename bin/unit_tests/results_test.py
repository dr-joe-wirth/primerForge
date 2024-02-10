from __future__ import annotations
from bin.Clock import Clock, _printDone, _printStart
from bin.getPrimerPairs import _formsDimers
import gzip, os, subprocess, sys, unittest
from bin.Parameters import Parameters
from Bio.SeqUtils import MeltingTemp
from Bio.SeqRecord import SeqRecord
from bin.Primer import Primer
from bin.main import _main
from Bio.Seq import Seq
from bin.Log import Log
from Bio import SeqIO

class Result():
    """class to save results for easy lookup
    """
    def __init__(self, fwdTm:float, fwdGc:float, revTm:float, revGc:float, additional:dict[str,dict[str,list]]) -> Result:
        self.fwdTm:float = fwdTm
        self.fwdGc:float = fwdGc
        self.revTm:float = revTm
        self.revGc:float = revGc
        self.additional:dict[str,dict[str,list]] = additional
    
    def __repr__(self) -> str:
        return str(vars(self))

class ResultsTest(unittest.TestCase):
    """class for testing the results of primerForge
    """
    # constants
    NUM_THREADS = 24
    TEST_DIR = os.path.join(os.getcwd(), "test_dir")
    RESULT_FN = os.path.join(TEST_DIR, "results.tsv")
    INGROUP  = {'i1.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/650/295/GCF_001650295.1_ASM165029v1/GCF_001650295.1_ASM165029v1_genomic.gbff.gz',
                'i2.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/727/175/GCF_008727175.1_ASM872717v1/GCF_008727175.1_ASM872717v1_genomic.gbff.gz',
                'i3.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/208/865/GCF_002208865.2_ASM220886v2/GCF_002208865.2_ASM220886v2_genomic.gbff.gz'}
    OUTGROUP = {'o1.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz',
                'o2.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gbff.gz'}
    PCRLEN = "length"
    CONTIG = "contig"
    PLS = "+"
    MNS = "-"
    
    @classmethod
    def setUpClass(cls) -> None:
        """sets up the class for all tests
        """
        # initialize clocks
        total = Clock()
        clock = Clock()
        
        # make the test directory if it doesn't exist
        if not os.path.isdir(ResultsTest.TEST_DIR):
            os.mkdir(ResultsTest.TEST_DIR)
        
        # get a list of all the genome files
        allFiles = list(ResultsTest.INGROUP.keys()) + list(ResultsTest.OUTGROUP.keys())
        allFiles = map(os.path.join, [ResultsTest.TEST_DIR]*len(allFiles), allFiles)
        
        # download genome files if necessary
        if not all(map(os.path.exists, allFiles)):
            _printStart(clock, 'downloading genomes')
            ResultsTest._downloadGenomesForGroup(ResultsTest.INGROUP)
            ResultsTest._downloadGenomesForGroup(ResultsTest.OUTGROUP)
            _printDone(clock)
        
        # load existing results if present
        if os.path.exists(ResultsTest.RESULT_FN):
            # get the parameters
            cls.params:Parameters = ResultsTest._getParameters()
            cls.params.log.initialize(ResultsTest.setUpClass.__name__)
            
            print(f'running tests on existing results file: {ResultsTest.RESULT_FN}')
        
        # otherwise run primerForge
        else:
            _printStart(total, 'getting results for testing', end=' ...\n')
            
            # get the parameters
            cls.params:Parameters = ResultsTest._getParameters()
            cls.params.log.initialize(ResultsTest.setUpClass.__name__)
            
            
            # run primerForge
            _printStart(clock, 'running primerForge', end=' ...\n')
            _main(cls.params)
            _printDone(clock)
            print()
            _printDone(total)
        
        # load the results file into memory
        _printStart(clock, 'reading results into memory')
        cls.results:dict[tuple[Seq,Seq],Result] = ResultsTest._parseResultsFile(ResultsTest.RESULT_FN)
        _printDone(clock)
    
        # load the genomic sequences into memory
        _printStart(clock, 'reading genomic sequences into memory')
        cls.sequences:dict[str,dict[str,dict[str,Seq]]] = ResultsTest._loadGenomeSequences()
        _printDone(clock)

        # get all the primer binding sites
        _printStart(clock, 'calculating primer binding sites in each genome')
        cls.bindingSites = ResultsTest._getAllBindingSites(cls.results, cls.sequences, cls.params.minLen, cls.params.maxLen)
        _printDone(clock)
    
    # functions for setting up the class
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
        
        # don't overwrite an existing results file
        if os.path.exists(ResultsTest.RESULT_FN):
            sys.argv[-2] = 'fakefile'
            params = Parameters('', '')
            params.outFn = ResultsTest.RESULT_FN
        
        # proceed like normal if the file doesn not yet exist
        else:
            params = Parameters('', '')
        
        # use a different logger pointing at the test directory
        params.log = Log(ResultsTest.TEST_DIR)
        
        return params

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
    
    def _downloadGenomesForGroup(group:dict[str,str]) -> None:
        """downloads a group of genomes

        Args:
            group (dict[str,str]): key=file basename; val=ftp path
        """
        # download each genome from ncbi
        for fn in group.keys():
            # get the ftp path and the filename
            ftp = group[fn]
            fn = os.path.join(ResultsTest.TEST_DIR, fn)
            
            # download the genome
            ResultsTest._downloadOneGenome(ftp, fn)

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
        LEN_SUFFIX = "_" + ResultsTest.PCRLEN
        CON_SUFFIX = "_" + ResultsTest.CONTIG
        
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
                            additional[name][ResultsTest.PCRLEN] = list(map(int, row[idx].split(SEP_2)))
                        
                        # otherwise it is a contig name
                        else:
                            # get the name of the sequence
                            name = key[idx][:-len(CON_SUFFIX)]
                            
                            # initialize a dictionary for the name if one deosn't exist
                            additional[name] = additional.get(name, dict())
                            
                            # then save the a list of contig names
                            additional[name][ResultsTest.CONTIG] = row[idx].split(SEP_2)
                    
                    # store the result for this primer pair
                    out[(fwd,rev)] = Result(ftm,fgc,rtm,rgc,additional)
        
        return out

    def _loadOneGenomeSequence(fn:str) -> dict[str,dict[str,Seq]]:
        """loads a single genome sequence into memory

        Args:
            fn (str): the filename (basename)

        Returns:
            dict[str,SeqRecord]: key=contig name; val=dict: key=strand; val=sequence
        """
        # initialize variables
        out = dict()
        contig:SeqRecord
        
        # store each contig underneath its name
        for contig in SeqIO.parse(os.path.join(ResultsTest.TEST_DIR, fn), 'genbank'):
            out[contig.id] = {ResultsTest.PLS: contig.seq.upper(),
                              ResultsTest.MNS: contig.seq.reverse_complement().upper()}
        
        return out

    def _loadGenomeSequences() -> dict[str,dict[str,dict[str,Seq]]]:
        """loads genome sequences into memory

        Returns:
            dict[str,list[SeqRecord]]: key=genome name; val=dict: key=contig name; val=dict: key=strand; val=sequence
        """
        # initialize output
        out = dict()
        
        # load ingroup sequences
        for fn in ResultsTest.INGROUP.keys():
            out[fn] = ResultsTest._loadOneGenomeSequence(fn)
        
        # load outgroup sequences
        for fn in ResultsTest.OUTGROUP.keys():
            out[fn] = ResultsTest._loadOneGenomeSequence(fn)
        
        return out
    
    def _getKmersForOneContig(seq:Seq, minLen:int, maxLen:int) -> dict[Seq,list[int]]:
        """gets all the kmers for a single contig strand

        Args:
            seq (Seq): the contig sequence from which to retrieve kmers
            minLen (int): the minimum kmer length
            maxLen (int): the maximum kmer length

        Returns:
            dict[Seq,list[int]]: key=kmers; val=list of start positions
        """
        # initialize variables
        krange = range(minLen, maxLen + 1)
        smallest = min(krange)
        seqLen = len(seq)
        done = False
        out = dict()
        
        # get every possible kmer start position
        for start in range(seqLen):
            # for each allowed kmer length
            for klen in krange:
                # stop looping through the contig once we're past the smallest kmer
                if start+smallest > seqLen:
                    done = True
                    break
                
                # stop looping through the kmers once the length is too long
                elif start+klen > seqLen:
                    break
            
                # proceed if the extracted kmer length is good
                else:
                    # extract the kmer sequences
                    kmer = seq[start:start+klen]
                    
                    # create dictionaries stored under this kmer
                    out[kmer] = out.get(kmer, list())
                    out[kmer].append(start)
            
            # stop iterating through the contig when we're done with it
            if done:
                break
        
        return out

    def _getAllBindingSites(results:dict[tuple[Seq,Seq],Result], sequences:dict[str,dict[str,dict[str,Seq]]], minLen:int, maxLen:int) -> dict[Seq,dict[str,dict[str,dict[str,list[int]]]]]:
        """get all of the primer binding sites in the provided genomes

        Args:
            results (dict[tuple[Seq,Seq],Result]): the dictionary produced by _parseResultsFile
            sequences (dict[str,dict[str,dict[str,Seq]]]): the dictionary produced by _loadGenomeSequences
            minLen (int): the minimum primer size
            maxLen (int): the maximum primer size

        Returns:
            dict[Seq,dict[str,dict[str,dict[str,list[int]]]]]: key=primer; val=dict:key=genome name; val=dict: key=contig; val=dict: key=strand; val=list of start positions
        """
        # initialize variables
        bindingSites = dict()
        allPrimers = {p for pair in results.keys() for p in pair}
        
        # for each genome
        for name in sequences.keys():
            # for each contig
            for contig in sequences[name].keys():
                # for each strand
                for strand in sequences[name][contig].keys():
                    # get the all the kmers
                    kmers = ResultsTest._getKmersForOneContig(sequences[name][contig][strand], minLen, maxLen)

                    # for each primer
                    for primer in allPrimers:
                        # build the data structure: seq, name, contig, strand, list of start positions
                        bindingSites[primer] = bindingSites.get(primer, dict())
                        bindingSites[primer][name] = bindingSites[primer].get(name, dict())
                        bindingSites[primer][name][contig] = bindingSites[primer][name].get(contig, dict())
                        bindingSites[primer][name][contig][strand] = bindingSites[primer][name][contig].get(strand, list())
                        
                        # only save data for the primers present in the results
                        try:
                            bindingSites[primer][name][contig][strand].extend(kmers[primer])
                        except KeyError:
                            pass
        
        return bindingSites

    # functions for testing
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
            if repeat in seq:
                return False
        return True
    
    # test cases
    def testA_isSequenceUpper(self) -> None:
        """is primer sequence uppercase
        """
        # constants
        FAIL_MSG = " is not uppercase"
        
        # make sure each primer saved is uppercase
        for fwd,rev in self.results.keys():
            self.assertEqual(fwd, fwd.upper(), f"{fwd}{FAIL_MSG}")
            self.assertEqual(rev, rev.upper(), f"{rev}{FAIL_MSG}")
    
    def testB_checkTm(self) -> None:
        """does the Tm match the expected value
        """
        # constant
        FAIL_MSG = "wrong tm for "
        
        # for each pair, make sure the saved Tm matches the expected value
        for fwd,rev in self.results.keys():
            self.assertEqual(self.results[(fwd,rev)].fwdTm, round(MeltingTemp.Tm_Wallace(fwd), 1), f"{FAIL_MSG}{fwd}")
            self.assertEqual(self.results[(fwd,rev)].revTm, round(MeltingTemp.Tm_Wallace(rev), 1), f"{FAIL_MSG}{rev}")
    
    def testC_tmWithinRange(self) -> None:
        """is the Tm within the specified range
        """
        # constant
        FAIL_MSG = "tm out of range for "
        
        # for each pair, make sure the Tms are within the specified ranges
        for fwd,rev in self.results.keys():
            self.assertGreaterEqual(self.results[(fwd,rev)].fwdTm, self.params.minTm, f"{FAIL_MSG}{fwd}")
            self.assertGreaterEqual(self.results[(fwd,rev)].revTm, self.params.minTm, f"{FAIL_MSG}{rev}")
            self.assertLessEqual(self.results[(fwd,rev)].fwdTm, self.params.maxTm, f"{FAIL_MSG}{fwd}")
            self.assertLessEqual(self.results[(fwd,rev)].revTm, self.params.maxTm, f"{FAIL_MSG}{rev}")
    
    def testD_tmDiffWithinRange(self) -> None:
        """is the Tm difference below the specified threshold
        """
        # for each pair, make sure the Tm difference is below the threshold
        for fwd,rev in self.results.keys():
            tmDiff = abs(self.results[(fwd,rev)].fwdTm - self.results[(fwd,rev)].revTm)
            self.assertLessEqual(tmDiff, self.params.maxTmDiff, f"tm diff to large for {fwd}, {rev}")
    
    def testE_checkGcPercent(self) -> None:
        """is the G+C percentage the expected value
        """
        # constant
        FAIL_MSG = "GC percent incorrect for "
        
        # for each pair, make sure the GC percent is accurate
        for fwd,rev in self.results.keys():
            self.assertEqual(self.results[(fwd,rev)].fwdGc, round(ResultsTest._getGc(fwd), 1), f"{FAIL_MSG}{fwd}")
            self.assertEqual(self.results[(fwd,rev)].revGc, round(ResultsTest._getGc(rev), 1), f"{FAIL_MSG}{fwd}")
    
    def testF_gcPercentInRange(self) -> None:
        """is the G+C percentage within the specified range
        """
        # constant
        FAIL_MSG = "GC percent out of range for "
        
        # for each pair, make sure the GC percent is within the specified range
        for fwd,rev in self.results.keys():
            self.assertGreaterEqual(self.results[(fwd,rev)].fwdGc, self.params.minGc, f"{FAIL_MSG}{fwd}")
            self.assertGreaterEqual(self.results[(fwd,rev)].revGc, self.params.minGc, f"{FAIL_MSG}{rev}")
            self.assertLessEqual(self.results[(fwd,rev)].fwdGc, self.params.maxGc, f"{FAIL_MSG}{fwd}")
            self.assertLessEqual(self.results[(fwd,rev)].revGc, self.params.maxGc, f"{FAIL_MSG}{rev}")
    
    def testG_threePrimeIsGc(self) -> None:
        """is the 3' end of the primer a G or C
        """
        # constants
        GC = ("G", "C")
        FAIL_MSG = "3' end not G|C in "
        
        # for each pair, make sure the 3' end is a G or a C
        for fwd,rev in self.results.keys():
            self.assertIn(fwd[-1], GC, f"{FAIL_MSG}{fwd}")
            self.assertIn(rev[-1], GC, f"{FAIL_MSG}{rev}")

    def testH_noHomoPolymers(self) -> None:
        """does the primer contain homopolymers
        """
        # constant
        FAIL_MSG = "homopolymers present in "
        
        # for each pair, make sure there are no homopolymers
        for fwd,rev in self.results.keys():
            self.assertTrue(ResultsTest._noLongRepeats(fwd), f"{FAIL_MSG}{fwd}")
            self.assertTrue(ResultsTest._noLongRepeats(rev), f"{FAIL_MSG}{rev}")

    def testI_noPrimerDimers(self) -> None:
        """do all pairs not produce primer dimers
        """
        for fwd,rev in self.results.keys():
            # create mock primers
            fwd = Primer(fwd, '', 100, len(fwd), Primer.PLUS)
            rev = Primer(rev, '', 100, len(rev), Primer.MINUS)

            # check for primer dimers; both primers need to be on the same strand
            self.assertFalse(_formsDimers(fwd,rev.reverseComplement()), f'primer dimers present in {fwd}, {rev}')

    def testJ_ingroupHasOneProduct(self) -> None:
        """do all pairs produce one product size in the ingroup
        """
        # constants
        FAIL_MSG_A = "multiple product sizes for "
        FAIL_MSG_B = " in "
        
        # for each pair in the ingroup
        for pair in self.results.keys():
            for name in ResultsTest.INGROUP.keys():
                # there should be exactly one pcr product size for each ingroup genome
                self.assertEqual(len(self.results[pair].additional[name][ResultsTest.PCRLEN]), 1, f"{FAIL_MSG_A}{pair}{FAIL_MSG_B}{name}")
                self.assertEqual(len(self.results[pair].additional[name][ResultsTest.CONTIG]), 1, f"{FAIL_MSG_A}{pair}{FAIL_MSG_B}{name}")
    
    def testK_ingroupProductsWithinRange(self) -> None:
        """is ingroup pcr product size in the expected range
        """
        # constants
        FAIL_MSG_A = "product size out of range for "
        FAIL_MSG_B = " in "
        
        # for each pari in the ingroup
        for pair in self.results.keys():
            for name in ResultsTest.INGROUP.keys():
                # make sure the pcr products are within the allowed range
                self.assertGreaterEqual(self.results[pair].additional[name][ResultsTest.PCRLEN][0], self.params.minPcr, f"{FAIL_MSG_A}{pair}{FAIL_MSG_B}{name}")
                self.assertLessEqual(self.results[pair].additional[name][ResultsTest.PCRLEN][0], self.params.maxPcr, f"{FAIL_MSG_A}{pair}{FAIL_MSG_B}{name}")
    
    def testL_outgroupProductsSavedCorrectly(self) -> None:
        """do outgroup pcr products look valid
        """
        # constants
        FAIL_MSG_A = "outgroup products saved incorrectly for "
        FAIL_MSG_B = " in "
        
        # for each pair in the outgroup
        for pair in self.results.keys():
            for name in ResultsTest.OUTGROUP.keys():
                # the length of the contig and pcr size lists should be equal
                numContigs = len(self.results[pair].additional[name][ResultsTest.CONTIG])
                numPcrLens = len(self.results[pair].additional[name][ResultsTest.PCRLEN])
                self.assertEqual(numContigs, numPcrLens, f"{FAIL_MSG_A}{pair}{FAIL_MSG_B}{name}")
    
    def testM_outgroupProductsNaAreZero(self) -> None:
        """do outgroup products of 0 have NA as the contig and vice-versa
        """
        # constants
        FAIL_MSG_1A = "NA is not 0 for "
        FAIL_MSG_2A = "0 is not NA for "
        FAIL_MSG_B  = " in "
        
        # for each pair in the outgroup
        for pair in self.results.keys():
            for name in ResultsTest.OUTGROUP.keys():
                # for each pcr product saved
                for idx in range(len(self.results[pair].additional[name][ResultsTest.CONTIG])):
                    # if the contig is NA, then the pcr size should be zero
                    if self.results[pair].additional[name][ResultsTest.CONTIG][idx] == 'NA':
                        self.assertEqual(self.results[pair].additional[name][ResultsTest.PCRLEN][idx], 0, f"{FAIL_MSG_1A}{pair}{FAIL_MSG_B}{name}")
                    
                    # if the pcr size is zero, then the contig should be NA
                    if self.results[pair].additional[name][ResultsTest.PCRLEN][idx] == 0:
                        self.assertEqual(self.results[pair].additional[name][ResultsTest.CONTIG][idx], 'NA', f"{FAIL_MSG_2A}{pair}{FAIL_MSG_B}{name}")
                
    def testN_outgroupProductsNotWithinDisallowedRange(self) -> None:
        """are outgroup pcr products outside the disallowed range
        """
        # constants
        FAIL_MSG_A = "outgroup product(s) within the disallowed range for "
        FAIL_MSG_B = " in "
        
        # for each pair in the outgroup
        for pair in self.results.keys():
            for name in ResultsTest.OUTGROUP.keys():
                # each pcr product size should not be in the disallowed range
                for pcrLen in self.results[pair].additional[name][ResultsTest.PCRLEN]:
                    self.assertNotIn(pcrLen, self.params.disallowedLens, f"{FAIL_MSG_A}{pair}{FAIL_MSG_B}{name}")

    def testO_isRepeatedInIngroup(self) -> None:
        """does the primer pair have exactly one binding site
        """
        # constants
        FAIL_MSG_A = " has "
        FAIL_MSG_B = " binding sites in "
        
        for fwd,rev in self.results.keys():
            fbind = self.bindingSites[fwd]
            rbind = self.bindingSites[rev]
            
            # for each ingroup genome
            for name in ResultsTest.INGROUP.keys():
                # count the number of binding sites in the genome for each primer
                fCount = 0
                rCount = 0
                for contig in fbind[name].keys():
                    fCount += len(fbind[name][contig][ResultsTest.PLS]) + len(fbind[name][contig][ResultsTest.MNS])
                    rCount += len(rbind[name][contig][ResultsTest.PLS]) + len(rbind[name][contig][ResultsTest.MNS])
            
                # each ingroup genome should have exactly one binding site for each primer
                self.assertEqual(rCount, 1, f"{rev}{FAIL_MSG_A}{rCount}{FAIL_MSG_B}{name}")
                self.assertEqual(fCount, 1, f"{fwd}{FAIL_MSG_A}{fCount}{FAIL_MSG_B}{name}")

    def testP_onSeparateStrandsInIngroup(self) -> None:
        """evaluates if the forward and reverse primers bind on separate strands in the ingroup
        """
        # for each primer pair
        for fwd,rev in self.results.keys():
            # get the binding sites
            fbind = self.bindingSites[fwd]
            rbind = self.bindingSites[rev]
            
            # only evaluate for the ingroup
            for name in ResultsTest.INGROUP.keys():
                for contig in fbind[name].keys():
                    # if the forward primer is on the (+) strand
                    if fbind[name][contig][ResultsTest.PLS] != []:
                        # then the reverse primer should be on the (-) strand
                        self.assertNotEqual(rbind[name][contig][ResultsTest.MNS], [], f"{rev} is not on opposite strand of {fwd} in {name}|{contig}")
                        
                        # the forward primer should not be on the (-) strand
                        self.assertEqual(fbind[name][contig][ResultsTest.MNS], [], f"{fwd} is on both strands in {name}|{contig}")
                        
                        # the reverse primer should not be on the (+) strand
                        self.assertEqual(rbind[name][contig][ResultsTest.PLS], [], f"{rev} is on the same strand as {fwd} in {name}|{contig}")
                    
                    # if the forward primer is on the (-) strand
                    elif fbind[name][contig][ResultsTest.MNS] != []:
                        # then the reverse primer should be on the (+) strand
                        self.assertNotEqual(rbind[name][contig][ResultsTest.PLS], [], f"{rev} is ont on opposite strand of {fwd} in {name}|{contig}")
                        
                        # the forward primer should not be on the (+) strand
                        self.assertEqual(fbind[name][contig][ResultsTest.PLS], [], f"{fwd} is on both strands in {name}|{contig}")
                        
                        # the reverse primer should not be on the (-) strand
                        self.assertEqual(rbind[name][contig][ResultsTest.MNS], [], f"{rev} is on the same strand as {fwd} in {name}|{contig}")
    
    def testQ_isProductLengthCorrect(self) -> None:
        """determines if the product length is correctly saved
        """
        # for each primer pair
        for fwd,rev in self.results.keys():
            # extract the binding sites
            fbind = self.bindingSites[fwd]
            rbind = self.bindingSites[rev]
            
            # for each ingroup genome
            for name in ResultsTest.INGROUP.keys():
                # extract the stored pcr length and the contig
                pcrLen = self.results[(fwd,rev)].additional[name][ResultsTest.PCRLEN][0]
                contig = self.results[(fwd,rev)].additional[name][ResultsTest.CONTIG][0]
                
                # find where the forward and reverse primers bind
                if fbind[name][contig][ResultsTest.PLS] != []:
                    fstart = fbind[name][contig][ResultsTest.PLS][0]
                    rstart = rbind[name][contig][ResultsTest.MNS][0]
                else:
                    fstart = fbind[name][contig][ResultsTest.MNS][0]
                    rstart = rbind[name][contig][ResultsTest.PLS][0]
                
                # the true length is just the space between the two ends
                truLen = len(self.sequences[name][contig][ResultsTest.PLS]) - fstart - rstart
                
                # make sure the saved length and the true length match
                self.assertEqual(pcrLen, truLen, f"bad pcr product sizes in {name} for {fwd}, {rev}")
            
            # for each outgroup genome
            for name in ResultsTest.OUTGROUP.keys():
                # get the list of of pcr products and their contigs
                pcrLens = self.results[(fwd,rev)].additional[name][ResultsTest.PCRLEN]
                contigs = self.results[(fwd,rev)].additional[name][ResultsTest.CONTIG]
                
                # key the product sizes under the contig names
                res = dict()
                for idx in range(len(pcrLens)):
                    res[contigs[idx]] = res.get(contigs[idx], set())
                    res[contigs[idx]].add(pcrLens[idx])

                # initialize a dictionary to store the true values
                tru = dict()
                
                # for each contig
                for contig in fbind[name].keys():
                    # only proceed if the reverse primer also binds to this contig
                    if contig in rbind[name].keys():
                        # extrac the forward and reverse binding sites on both strands
                        fwdPlsStarts = fbind[name][contig][ResultsTest.PLS]
                        fwdMnsStarts = fbind[name][contig][ResultsTest.MNS]
                        revPlsStarts = rbind[name][contig][ResultsTest.PLS]
                        revMnsStarts = rbind[name][contig][ResultsTest.MNS]
                        
                        # extract the length of the contig
                        contigLen = len(self.sequences[name][contig][ResultsTest.PLS])
                        
                        # for each fwd/rev pair on the (+)/(-) strands
                        for fstart in fwdPlsStarts:
                            for rstart in revMnsStarts:
                                # calculate the pcr lengths
                                pcrLen = contigLen - fstart - rstart
                                
                                # save the pcr length if it is positive
                                if pcrLen > 0:
                                    tru[contig] = tru.get(contig, set())
                                    tru[contig].add(pcrLen)
                        
                        # for each fwd/rev pair on the (-)/(+) strands
                        for fstart in fwdMnsStarts:
                            for rstart in revPlsStarts:
                                # calculate the pcr lengths
                                pcrLen = contigLen - fstart - rstart
                                
                                # only save the pcr length if it is positive
                                if pcrLen > 0:
                                    tru[contig] = tru.get(contig, set())
                                    tru[contig].add(pcrLen)
                
                # if the value was "NA", then there should be no saved data
                if contigs == ["NA"]:
                    self.assertEqual(tru, dict(), f"'NA' in {name} but products produced by {fwd},{rev}")
                
                else:
                    # for each contig
                    for contig in res.keys():
                        # the pcr lengths should be equal
                        self.assertEqual(res[contig], tru[contig], f"unexpected sizes in {name} using {fwd}, {rev}")
                        
                        # and the pcr lengths should not be disallowed
                        for plen in tru[contig]:
                            self.assertNotIn(plen, self.params.disallowedLens, f"disallowed sizes in {name} with {fwd},{rev}")
