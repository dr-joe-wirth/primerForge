from __future__ import annotations
from Bio import SeqIO
from Bio.Seq import Seq
from bin.Log import Log
from bin.main import _runner
from bin.Clock import Clock
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters
from bin.AnalysisData import AnalysisData
from bin.getPrimerPairs import _formsDimers
import gzip, os, pickle, primer3, re, subprocess, sys, unittest

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
    TEST_DIR = os.path.join(os.getcwd(), "test_dir")
    RESULT_FN = os.path.join(TEST_DIR, "results.tsv")
    ANALYSIS_BASENAME = os.path.join(TEST_DIR, "distribution")
    FAKE_FN = 'fakefile'
    BINDING_SITES_FN = os.path.join(TEST_DIR, "bindingSites.p")
    INGROUP  = {'i1.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/650/295/GCF_001650295.1_ASM165029v1/GCF_001650295.1_ASM165029v1_genomic.gbff.gz',
                'i2.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/727/175/GCF_008727175.1_ASM872717v1/GCF_008727175.1_ASM872717v1_genomic.gbff.gz',
                'i3.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/208/865/GCF_002208865.2_ASM220886v2/GCF_002208865.2_ASM220886v2_genomic.gbff.gz'}
    OUTGROUP = {'o1.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz',
                'o2.gbff': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gbff.gz'}
    PCRLEN = "length"
    CONTIG = "contig"
    
    @classmethod
    def setUpClass(cls) -> None:
        """sets up the class for all tests
        """
        # initialize clocks
        clock = Clock()
        
        # get the number of threads to use
        numThreads = ResultsTest._getNumThreads()
        
        # make the test directory if it doesn't exist
        if not os.path.isdir(ResultsTest.TEST_DIR):
            os.mkdir(ResultsTest.TEST_DIR)
        
        # get a list of all the genome files
        allFiles = list(ResultsTest.INGROUP.keys()) + list(ResultsTest.OUTGROUP.keys())
        allFiles = map(os.path.join, [ResultsTest.TEST_DIR]*len(allFiles), allFiles)
        
        # download genome files if necessary
        if not all(map(os.path.exists, allFiles)):
            clock.printStart('downloading genomes')
            ResultsTest._downloadGenomesForGroup(ResultsTest.INGROUP)
            ResultsTest._downloadGenomesForGroup(ResultsTest.OUTGROUP)
            clock.printDone()
        
        # get the parameters
        cls.params:Parameters = ResultsTest._getParameters(numThreads)
        cls.params.log.rename(ResultsTest.setUpClass.__name__)
        
        # run primerForge if results file does not exist
        if not os.path.exists(ResultsTest.params.resultsFn) or not os.path.exists(ResultsTest.params.plotDataFn):  
            # run primerForge
            clock.printStart('running primerForge', end=' ...\n', spin=False)
            _runner(cls.params)
            clock.printDone()
        
        # otherwise using existing file
        else:
            print('running tests on existing files')
        
        # load the results file into memory
        clock.printStart('reading results into memory')
        cls.results:dict[tuple[Seq,Seq],Result] = ResultsTest._parseResultsFile(ResultsTest.params.resultsFn)
        clock.printDone()
        
        # load the analysis file into memory
        clock.printStart('reading analysis data into memory')
        cls.analysis:dict[tuple[Seq,str],AnalysisData] = ResultsTest._parseAnalysisData(ResultsTest.params.plotDataFn)
        clock.printDone()
    
        # load the genomic sequences into memory
        clock.printStart('reading genomic sequences into memory')
        cls.sequences:dict[str,dict[str,dict[str,Seq]]] = ResultsTest._loadGenomeSequences()
        clock.printDone()

        # load the binding sites if the file already exists
        if os.path.exists(ResultsTest.BINDING_SITES_FN):
            clock.printStart('loading binding sites from file')
            with open(ResultsTest.BINDING_SITES_FN, 'rb') as fh:
                cls.bindingSites:dict[Seq,dict[str,dict[str,dict[str,list[int]]]]] = pickle.load(fh)
            clock.printDone()
        
        # otherwise calculate and dump the binding sites
        else:
            # calculate the binding sites
            clock.printStart('getting binding sites for all kmers and primers', end='...\n', spin=False)
            cls.bindingSites:dict[Seq,dict[str,dict[str,dict[str,list[int]]]]] = ResultsTest._getKmerBindingSites(cls.results, cls.analysis, cls.sequences, cls.params.minLen, cls.params.maxLen)
            clock.printDone()
            
            # dump the binding sites
            clock.printStart('dumping binding sites to file')
            with open(ResultsTest.BINDING_SITES_FN, 'wb') as fh:
                pickle.dump(cls.bindingSites, fh)
            clock.printDone()
    
    # functions for setting up the class
    def _getNumThreads() -> int:
        # message
        MSG = "\nplease specify a number of threads to use: "
        
        # get the number of threads
        numThreads = input(MSG)

        # make sure the threads is an int
        while not type(numThreads) is int:
            try:
                numThreads = int(numThreads)
            except:
                print('invalid number of threads')
                numThreads = input(MSG)
        
        return numThreads
    
    def _getParameters(numThreads:int) -> Parameters:
        """creates a Parameters object for testing
        
        Args:
            numThreads (int): the number of threads to use

        Returns:
            Parameters: a Parameters object
        """ 
        # make the parameters object
        sys.argv = ['primerForge.py',
                    '-i', os.path.join(ResultsTest.TEST_DIR, "i[123].gbff"),
                    '-u', os.path.join(ResultsTest.TEST_DIR, "o[12].gbff"),
                    '-r', '70,120',
                    '-b', '45,150',
                    '-n', numThreads,
                    '-o', ResultsTest.RESULT_FN,
                    '-a', ResultsTest.ANALYSIS_BASENAME,
                    '--debug']
        
        # don't overwrite existing results or analysis files
        if os.path.exists(ResultsTest.RESULT_FN):
            sys.argv[12] = ResultsTest.FAKE_FN
            sys.argv[14] = ResultsTest.FAKE_FN
            params = Parameters('', '')
            params.plotDataFn = ResultsTest.ANALYSIS_BASENAME + Parameters._DATA_EXT
            params.plotsFn = ResultsTest.ANALYSIS_BASENAME + Parameters._PLOT_EXT
            params.resultsFn = ResultsTest.RESULT_FN
        
        # proceed like normal if the files do not yet exist
        else:
            params = Parameters('', '', initializeLog=False)
        
        # move the log file to the test directory
        params.log = Log(debugDir=ResultsTest.TEST_DIR, debug=True)
        
        # move the pickle directory into the test directory
        oldDir = os.path.dirname(next(iter(params.pickles.values())))
        newDir = os.path.join(ResultsTest.TEST_DIR, os.path.basename(oldDir))
        
        # move the directory, but allow it to fail for loading parameters outside the tests
        try:
            os.rename(oldDir, newDir)
        except:
            pass
        
        # update the location of the pickle files
        for key in params.pickles.keys():
            params.pickles[key] = os.path.join(newDir, os.path.basename(params.pickles[key]))
        
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
                # convert the line to a row of columns
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

    def _parseAnalysisData(fn:str) -> dict[tuple[Seq,str], AnalysisData]:
        """parses an analysis data file

        Args:
            fn (str): the file to parse

        Returns:
            dict[tuple[Seq,str], AnalysisData]: key=(kmer seq, contig name); val=corresponding AnalysisData object
        """
        # constants
        SEP = "\t"
        G_COORD = 'genome coord'
        GREP_FIND = r'^(.+) \((\d+), (\d+)\)$'
        GREP_REPL = r'\1\t\2\t\3'
        
        # initialize variables
        header = True
        indices = dict()
        recNum = 0
        out = dict()
        
        # for each line in the file
        with open(fn, 'r') as fh:
            for line in fh:
                # convert to a row of columns
                line = line.rstrip()
                row = line.split(SEP)
                
                # save the index for each genome's kmer coordinates
                if header:
                    # key = index; val = genome name
                    indices = {idx:row[idx][:-len(G_COORD)].rstrip() for idx in range(2, len(row)) if G_COORD in row[idx]}
                    header  = False
                
                # process data rows
                else:
                    # extract the sequence and the level
                    seq = Seq(row[0])
                    level = row[1]

                    # for each genome
                    for idx in indices.keys():
                        # get the name
                        name = indices[idx]

                        # extract the positional data
                        contig, start, end = re.sub(GREP_FIND, GREP_REPL, row[idx]).split(SEP)
                        
                        # cast coordinates as integers
                        start = int(start)
                        end = int(end)
                        
                        # determine the strand; `start` needs to be lowest value
                        if start < end:
                            strand = Primer.PLUS
                        else:
                            start = end
                            strand = Primer.MINUS
                        
                        # create a primer object
                        primer = Primer(seq, contig, start, len(seq), strand)
                        
                        # create and save an analysis data object; use the row as the index
                        data = AnalysisData(primer, recNum, name)
                        data.setLevel(level)
                        out[(seq,contig)] = data
                        
                        # update the record number for the next entry
                        recNum += 1
        
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
            out[contig.id] = {Primer.PLUS:  contig.seq.upper(),
                              Primer.MINUS: contig.seq.reverse_complement().upper()}
        
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
    
    def _getKmersForOneContig(seq:dict[str,Seq], minLen:int, maxLen:int) -> dict[str,dict[Seq,list[int]]]:
        """gets all the kmers for a single contig

        Args:
            seq (dict[str,Seq]): key=strand; val=contig sequence
            minLen (int): the minimum kmer length
            maxLen (int): the maximum kmer length

        Returns:
            dict[str,dict[Seq,list[int]]]: key=strand; val=dict: key=kmers; val=list of start positions
        """
        # initialize variables
        krange = range(minLen, maxLen + 1)
        smallest = min(krange)
        seqLen = len(seq[Primer.PLUS])
        done = False
        out = {Primer.PLUS:  dict(),
               Primer.MINUS: dict()}
        
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
                    # calculate the end
                    end = start + klen
                    
                    # extract the forward sequence
                    fKmer = seq[Primer.PLUS][start:end]
                    
                    # extract the reverse sequence
                    if start == 0:
                        rKmer = seq[Primer.MINUS][-end:]
                    else:
                        rKmer = seq[Primer.MINUS][-end:-start]
                    
                    # initialize lists stored under the kmers
                    out[Primer.PLUS][fKmer] = out[Primer.PLUS].get(fKmer, list())
                    out[Primer.MINUS][rKmer] = out[Primer.MINUS].get(rKmer, list())
                    
                    # save the start positions
                    out[Primer.PLUS][fKmer].append(start)
                    out[Primer.MINUS][rKmer].append(end-1)
            
            # stop iterating through the contig when we're done with it
            if done:
                break
        
        return out

    def _getKmerBindingSites(results:dict[tuple[Seq,Seq],Result], analysis:dict[tuple[Seq,str],AnalysisData], sequences:dict[str,dict[str,dict[str,Seq]]], minLen:int, maxLen:int) -> dict[Seq,dict[str,dict[str,dict[str,list[int]]]]]:
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
        savedKmers = {k for k,c in analysis.keys()}
        savedKmers.update(allPrimers)
        clock = Clock()
        
        # for each genome
        for name in sequences.keys():
            # print status
            clock.printStart(f'    getting kmer binding sites for {name}')
            
            # for each contig
            for contig in sequences[name].keys():
                # get the all the kmers
                kmers = ResultsTest._getKmersForOneContig(sequences[name][contig], minLen, maxLen)

                # determine which kmers actually need to be evaluated
                relevantKmers = set(kmers[Primer.PLUS])
                relevantKmers.update(set(kmers[Primer.MINUS]))                
                relevantKmers.intersection_update(savedKmers)
                
                # for each kmer
                for kmer in savedKmers:
                    # build the data structure: seq, name, contig, strand, list of start positions
                    bindingSites[kmer] = bindingSites.get(kmer, dict())
                    bindingSites[kmer][name] = bindingSites[kmer].get(name, dict())
                    bindingSites[kmer][name][contig] = bindingSites[kmer][name].get(contig, dict())
                    bindingSites[kmer][name][contig][Primer.PLUS] = bindingSites[kmer][name][contig].get(Primer.PLUS, list())
                    bindingSites[kmer][name][contig][Primer.MINUS] = bindingSites[kmer][name][contig].get(Primer.MINUS, list())
                    
                    # only add start positions for kmers in this contig
                    if kmer in relevantKmers:
                        # make sure these lists are empty
                        plsStarts = list()
                        mnsStarts = list()
                        
                        # try to get the plus start positions
                        try:
                            plsStarts = kmers[Primer.PLUS][kmer]
                        except KeyError:
                            pass
                        
                        # try to get the minus start positions
                        try:
                            mnsStarts = kmers[Primer.MINUS][kmer]
                        except KeyError:
                            pass
                        
                        # save the start positions for this kmer
                        try:
                            bindingSites[kmer][name][contig][Primer.PLUS].extend(plsStarts)
                            bindingSites[kmer][name][contig][Primer.MINUS].extend(mnsStarts)
                        
                        # the kmers will not necessarily be present in every contig
                        except KeyError:
                            pass
                
            # print status
            clock.printDone()
        
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
    
    def _noHomoDimers(seq:Seq, minTm:float) -> bool:
        """evaluates sequences for the presence of homodimers

        Args:
            seq (Seq): a DNA sequence
            minTm (float): the minimum primer Tm

        Returns:
            bool: indicates if the sequence does not form homodimers
        """
        # constant
        FIVE_DEGREES = 5
        
        # homodimers don't exist as long as the tm for formation is 5° below the min Tm
        return primer3.calc_homodimer_tm(str(seq)) < (minTm - FIVE_DEGREES)
    
    def _noHairpins(seq:Seq, minTm:float) -> bool:
        """evaluates sequences for the presence of hairpins

        Args:
            seq (Seq): a DNA sequence
            minTm (float): the minimum primer Tm

        Returns:
            bool: indicates if the sequence does not form hairpins
        """
        # constant
        FIVE_DEGREES = 5
        
        # hairpins don't exist as long as the tm for formation is 5° below the min Tm
        return primer3.calc_hairpin_tm(str(seq)) < (minTm - FIVE_DEGREES)
    
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
            self.assertEqual(self.results[(fwd,rev)].fwdTm, round(primer3.calc_tm(str(fwd)), 1), f"{FAIL_MSG}{fwd}")
            self.assertEqual(self.results[(fwd,rev)].revTm, round(primer3.calc_tm(str(rev)), 1), f"{FAIL_MSG}{rev}")
    
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

    def testI_noHomoDimers(self) -> None:
        # constant
        FAIL_MSG = "homodimer present in "
        
        # for each pair, make sure there are no homodimers
        for fwd,rev in self.results.keys():
            self.assertTrue(ResultsTest._noHomoDimers(fwd, self.params.minTm), f"{FAIL_MSG}{fwd}")
            self.assertTrue(ResultsTest._noHomoDimers(rev, self.params.minTm), f"{FAIL_MSG}{rev}")
    
    def testJ_noHairpins(self) -> None:
        # constant
        FAIL_MSG = "hairpin present in "
        
        # for each pair, make sure there are no homodimers
        for fwd,rev in self.results.keys():
            self.assertTrue(ResultsTest._noHairpins(fwd, self.params.minTm), f"{FAIL_MSG}{fwd}")
            self.assertTrue(ResultsTest._noHairpins(rev, self.params.minTm), f"{FAIL_MSG}{rev}")

    def testK_noPrimerDimers(self) -> None:
        """do all pairs not produce primer dimers
        """
        for fwd,rev in self.results.keys():
            # create mock primers
            fwd = Primer(fwd, '', 100, len(fwd), Primer.PLUS)
            rev = Primer(rev, '', 100, len(rev), Primer.MINUS)

            # check for primer dimers; both primers need to be on the same strand
            self.assertFalse(_formsDimers(fwd,rev.reverseComplement()), f'primer dimers present in {fwd}, {rev}')

    def testL_ingroupHasOneProduct(self) -> None:
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
    
    def testM_ingroupProductsWithinRange(self) -> None:
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
    
    def testN_outgroupProductsSavedCorrectly(self) -> None:
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
    
    def testO_outgroupProductsNaAreZero(self) -> None:
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

    def testP_outgroupProductsNotWithinDisallowedRange(self) -> None:
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

    def testQ_isRepeatedInIngroup(self) -> None:
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
                    fCount += len(fbind[name][contig][Primer.PLUS]) + len(fbind[name][contig][Primer.MINUS])
                    rCount += len(rbind[name][contig][Primer.PLUS]) + len(rbind[name][contig][Primer.MINUS])
            
                # each ingroup genome should have exactly one binding site for each primer
                self.assertEqual(rCount, 1, f"{rev}{FAIL_MSG_A}{rCount}{FAIL_MSG_B}{name}")
                self.assertEqual(fCount, 1, f"{fwd}{FAIL_MSG_A}{fCount}{FAIL_MSG_B}{name}")

    def testR_onSeparateStrandsInIngroup(self) -> None:
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
                    if fbind[name][contig][Primer.PLUS] != []:
                        # then the reverse primer should be on the (-) strand
                        self.assertNotEqual(rbind[name][contig][Primer.MINUS], [], f"{rev} is not on opposite strand of {fwd} in {name}|{contig}")
                        
                        # the forward primer should not be on the (-) strand
                        self.assertEqual(fbind[name][contig][Primer.MINUS], [], f"{fwd} is on both strands in {name}|{contig}")
                        
                        # the reverse primer should not be on the (+) strand
                        self.assertEqual(rbind[name][contig][Primer.PLUS], [], f"{rev} is on the same strand as {fwd} in {name}|{contig}")
                    
                    # if the forward primer is on the (-) strand
                    elif fbind[name][contig][Primer.MINUS] != []:
                        # then the reverse primer should be on the (+) strand
                        self.assertNotEqual(rbind[name][contig][Primer.PLUS], [], f"{rev} is ont on opposite strand of {fwd} in {name}|{contig}")
                        
                        # the forward primer should not be on the (+) strand
                        self.assertEqual(fbind[name][contig][Primer.PLUS], [], f"{fwd} is on both strands in {name}|{contig}")
                        
                        # the reverse primer should not be on the (-) strand
                        self.assertEqual(rbind[name][contig][Primer.MINUS], [], f"{rev} is on the same strand as {fwd} in {name}|{contig}")
    
    def testS_isProductLengthCorrect(self) -> None:
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
                if fbind[name][contig][Primer.PLUS] != []:
                    fstart = fbind[name][contig][Primer.PLUS][0]
                    rstart = rbind[name][contig][Primer.MINUS][0]
                    plus = True
                else:
                    fstart = fbind[name][contig][Primer.MINUS][0]
                    rstart = rbind[name][contig][Primer.PLUS][0]
                    plus = False
                
                # the true length is just the space between the two ends
                if plus:
                    truLen = rstart - fstart + 1
                else:
                    truLen = fstart - rstart + 1
                
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
                    # extract the forward and reverse binding sites on both strands
                    fwdPlsStarts = fbind[name][contig][Primer.PLUS]
                    fwdMnsStarts = fbind[name][contig][Primer.MINUS]
                    revPlsStarts = rbind[name][contig][Primer.PLUS]
                    revMnsStarts = rbind[name][contig][Primer.MINUS]
                    
                    # for each fwd/rev pair on the (+)/(-) strands
                    for fstart in fwdPlsStarts:
                        for rstart in revMnsStarts:
                            # calculate the pcr lengths
                            pcrLen = rstart - fstart + 1
                            
                            # save the pcr length if it is positive
                            if pcrLen > 0:
                                tru[contig] = tru.get(contig, set())
                                tru[contig].add(pcrLen)
                    
                    # for each fwd/rev pair on the (-)/(+) strands
                    for fstart in fwdMnsStarts:
                        for rstart in revPlsStarts:
                            # calculate the pcr lengths
                            pcrLen = fstart - rstart + 1
                            
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
    
    def testT_analysisDataMatchesKmerPositions(self) -> None:
        """check that the analysis data are congruent with the binding sites
        """
        # for each kmer in the analysis file
        for (kmer,contig),analysis in self.analysis.items():
            # make sure the kmer is present in the genomes
            self.assertIn(kmer, self.bindingSites.keys())
            
            # make sure the kmer binds exactly one strand
            bindsPlus  = self.bindingSites[kmer][analysis.name][contig][Primer.PLUS]  != []
            bindsMinus = self.bindingSites[kmer][analysis.name][contig][Primer.MINUS] != []
            
            self.assertFalse(bindsPlus and bindsMinus, f"{kmer} binds multiple strands")
            self.assertTrue(bindsPlus  or  bindsMinus, f'{kmer} does not bind either strand')
        
            # extract the start position of the kmer
            if bindsPlus:
                starts = self.bindingSites[kmer][analysis.name][contig][Primer.PLUS]
            elif bindsMinus:
                starts = self.bindingSites[kmer][analysis.name][contig][Primer.MINUS]
            
            # should have at least one start position
            self.assertNotEqual(starts, [], f"{(kmer,contig)} has no start positions")
            
            # only one start position should exist for all kmers
            self.assertEqual(len(starts), 1, f"{(kmer,contig)} has multiple binding sites")
            
            # the start positions should match
            self.assertEqual(analysis.primer.start, starts[0], f"{(kmer,contig)}: start position does not match")
    
    def testU_resultsMatchesAnalysis(self) -> None:
        """checks that the results data and the analysis data are congruent
        """
        # track seen kmers to prevent redundant work
        seen = set()
        
        # for each pair
        for pair in self.results.keys():
            # for each genome
            for name in self.results[pair].additional.keys():
                # for each contig
                for contig in self.results[pair].additional[name][ResultsTest.CONTIG]:
                    # for each kmer in the pair
                    for kmer in pair:
                        # only process each kmer once
                        if kmer not in seen:
                            # mark the kmer as seen
                            seen.add(kmer)
                            
                            # the kmer (or its rev comp) should be in analysis data
                            try:
                                analysis = self.analysis[(kmer, contig)]
                                
                            except KeyError:
                                try:
                                    analysis = self.analysis[(kmer.reverse_complement(), contig)]
                                
                                except KeyError:
                                    self.fail(f'{(kmer, contig)} from {pair} not in analysis data')

                            # everything in the results should have the maximum level
                            self.assertEqual(analysis.getLevel(), max(AnalysisData.LEVELS))
