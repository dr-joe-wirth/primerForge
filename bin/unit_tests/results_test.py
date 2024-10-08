from __future__ import annotations

import pathlib, sys
sys.path.append(str(pathlib.Path(__file__).parent.parent.parent))

from Bio import SeqIO
from Bio.Seq import Seq
from bin.Log import Log
from bin.Clock import Clock
from bin.main import _runner
from khmer import Countgraph
from bin.Primer import Primer
from bin.Product import Product
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters
from bin.sortPrimerPairs import _sortPairs
from bin.getPrimerPairs import _formsDimers
import glob, gzip, os, pickle, primer3, signal, subprocess, time, unittest

class Result():
    """class to save results for easy lookup
    """
    def __init__(self, rowNum:int, fwdTm:float, fwdGc:float, revTm:float, revGc:float, additional:dict[str,dict[str,list]]) -> Result:
        self.rowNum:int  = rowNum
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
    TEST_DIR = os.path.join(str(pathlib.Path(__file__).parent.parent.parent), "test_dir")
    RESULT_FN = os.path.join(TEST_DIR, Parameters._DEF_RESULTS_FN)
    BED_FN = os.path.join(TEST_DIR, Parameters._DEF_BED_FN)
    FAKE_FN = 'fakefile'
    PCRLEN = "length"
    CONTIG = "contig"
    INGROUP_FILES = ('i1.gbff', 'i2.gbff', 'i3.gbff')
    OUTGROUP_FILES = ('o1.gbff', 'o2.gbff')
    
    # default ingroup and outgroup genomes (Mycoplasma mycoides)
    DEFAULT_INGROUP  = {INGROUP_FILES[0]:  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/034/305/GCF_003034305.1_ASM303430v1/GCF_003034305.1_ASM303430v1_genomic.gbff.gz',
                        INGROUP_FILES[1]:  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/034/275/GCF_003034275.1_ASM303427v1/GCF_003034275.1_ASM303427v1_genomic.gbff.gz',
                        INGROUP_FILES[2]:  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/034/345/GCF_003034345.1_ASM303434v1/GCF_003034345.1_ASM303434v1_genomic.gbff.gz'}
    DEFAULT_OUTGROUP = {OUTGROUP_FILES[0]: 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/489/555/GCF_900489555.1_MMC68/GCF_900489555.1_MMC68_genomic.gbff.gz',
                        OUTGROUP_FILES[1]: 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/389/745/GCF_018389745.1_ASM1838974v1/GCF_018389745.1_ASM1838974v1_genomic.gbff.gz'}
    
    # alternate ingroup and outgroup genomes (Escherichia coli)
    ALTRNTE_INGROUP  = {INGROUP_FILES[0]:  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/650/295/GCF_001650295.1_ASM165029v1/GCF_001650295.1_ASM165029v1_genomic.gbff.gz',
                        INGROUP_FILES[1]:  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/727/175/GCF_008727175.1_ASM872717v1/GCF_008727175.1_ASM872717v1_genomic.gbff.gz',
                        INGROUP_FILES[2]:  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/208/865/GCF_002208865.2_ASM220886v2/GCF_002208865.2_ASM220886v2_genomic.gbff.gz'}
    ALTRNTE_OUTGROUP = {OUTGROUP_FILES[0]: 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz',
                        OUTGROUP_FILES[1]: 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gbff.gz'}
    
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
        
        # get a list of the genome filenames
        cls.seqFiles:list[str] = list(map(lambda x: os.path.join(ResultsTest.TEST_DIR, x),
                                          ResultsTest.INGROUP_FILES + ResultsTest.OUTGROUP_FILES))
        
        # download genome files if necessary
        if not all(map(os.path.exists, ResultsTest.seqFiles)):
            # determine which dataset to use
            ingroupFtp, outgroupFtp = ResultsTest._selectDataset()
            
            # download the files
            clock.printStart('downloading genomes')
            ResultsTest._downloadGenomesForGroup(ingroupFtp)
            ResultsTest._downloadGenomesForGroup(outgroupFtp)
            clock.printDone()
        
        # get the parameters
        cls.params:Parameters = ResultsTest._getParameters(numThreads)
        cls.params.log.rename(ResultsTest.setUpClass.__name__)
        
        # run primerForge if results file or bed file does not exist
        if not os.path.exists(ResultsTest.params.resultsFn) or not os.path.exists(ResultsTest.params.bedFn):
            # run primerForge
            clock.printStart('running primerForge', end=' ...\n\n', spin=False)
            _runner(ResultsTest.params)
            clock.printDone()
        
        # otherwise using existing file
        else:
            print('running tests on existing files\n')
        
        # load the results file into memory
        clock.printStart('reading results into memory')
        cls.results:dict[tuple[Seq,Seq],Result] = ResultsTest._parseResultsFile(ResultsTest.params.resultsFn)
        clock.printDone()
        
        # load the bed file into memory
        clock.printStart('reading BED file into memory')
        cls.bed:dict[tuple[Seq,int],dict[str,tuple[int,int,str]]] = ResultsTest._parseBedFile(ResultsTest.params.bedFn)
        clock.printDone()

        # load the genomic sequences into memory
        clock.printStart('reading genomic sequences into memory')
        cls.sequences:dict[str,dict[str,dict[str,Seq]]] = ResultsTest._loadGenomeSequences()
        clock.printDone()
        
        # load the pairs from the pickle
        clock.printStart('unpickling pairs')
        with open(ResultsTest.params.pickles[Parameters._PAIR_3], 'rb') as fh:
            cls.pairs:dict[tuple[Primer,Primer],dict[str,Product]] = pickle.load(fh)
        clock.printDone()
        
        # calculate the binding sites
        clock.printStart('getting binding sites for all kmers and primers', end='...\n', spin=False)
        cls.bindingSites:dict[Seq,dict[str,dict[str,dict[str,list[int]]]]] = ResultsTest._getKmerBindingSites(cls.results, cls.sequences, cls.params.minLen, cls.params.maxLen)
        clock.printDone()
    
    # functions for setting up the class
    def _getNumThreads() -> int:
        """gets the number of threads

        Returns:
            int: the number of cpus; default value used after a timeout
        """
        # constants
        TIMEOUT = 10
        DEFAULT_NUM = 4
        
        # messages
        MSG = f"\nplease specify a number of cpus to use (defaults to {DEFAULT_NUM} after {TIMEOUT} seconds): "
        INVALID = '\ninvalid number of threads'
        
        # helper function for handling signals
        def handler(signum, frame):
            raise TimeoutError
        
        # initialize a signal handler and numThreads
        signal.signal(signal.SIGALRM, handler)
        numThreads = None
        
        # return default value if not running interactively
        if not sys.stdin.isatty():
            return DEFAULT_NUM
        
        # make sure the threads is an int
        while not type(numThreads) is int:
            # start the timeout clock
            signal.alarm(TIMEOUT)
            
            # get the user input
            try:
                numThreads = int(input(MSG))
            
            # print a message then try again if invalid input
            except ValueError:
                print(INVALID)
                numThreads = input(MSG)
            
            # use the default value after the timeout
            except TimeoutError:
                print(f"timeout; using {DEFAULT_NUM} threads")
                numThreads = int(DEFAULT_NUM)
        
        # turn off the alarm
        signal.alarm(0)
        
        return numThreads
    
    def _selectDataset() -> tuple[dict[str,str],dict[str,str]]:
        """allows the user to interactively select the test dataset

        Returns:
            tuple[dict[str,str],dict[str,str]]: an ingroup and an outgroup (see class attributes for details)
        """
        # constants
        TIMEOUT = 10
        ALLOWED = ('mycoplasma', 'ecoli')
        
        # messages
        MSG = f"\nplease specify a dataset {str(ALLOWED)} (default: {ALLOWED[0]}): "
        INVALID = '\ninvalid dataset'
        
        # helper function for handling signals
        def handler(signum, frame):
            raise TimeoutError
        
        # initialize a signal handler and numThreads
        signal.signal(signal.SIGALRM, handler)
        dataset = None
        
        # run interactively
        if sys.stdin.isatty():
            # make sure the threads is an int
            while True:
                # start the timeout clock
                signal.alarm(TIMEOUT)
                
                try:
                    # get the user input
                    dataset = input(MSG)
                    
                # use the default value after the timeout
                except TimeoutError:
                    print(f"timeout; using {ALLOWED[0]} dataset")
                    dataset = ALLOWED[0]

                # default value
                if dataset == ALLOWED[0]:
                    out = ResultsTest.DEFAULT_INGROUP, ResultsTest.DEFAULT_OUTGROUP
                    break
                
                # alternate value
                elif dataset == ALLOWED[1]:
                    out = ResultsTest.ALTRNTE_INGROUP, ResultsTest.ALTRNTE_OUTGROUP
                    break
                
                # otherwise it is invalid
                else:
                    print(INVALID)

            # turn off the alarm
            signal.alarm(0)
        
        # default value if not running interactively
        else:
            out = ResultsTest.DEFAULT_INGROUP, ResultsTest.DEFAULT_OUTGROUP
        
        return out
    
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
                    '-B', ResultsTest.BED_FN,
                    '--debug']
        
        # get the parameters
        params = Parameters('', '', initializeLog=False)
        
        # move the log file to the test directory
        params.log = Log(debugDir=ResultsTest.TEST_DIR, debug=True)
        
        # check if a run has already been done
        existing = glob.glob(os.path.join(ResultsTest.TEST_DIR, Parameters._WORKDIR_PREFIX + '*'))
        if existing != []:
            oldDir = os.path.basename(existing[0])
        
        # otherwise get the old directory from the newly created parameters object
        else:
            oldDir = os.path.dirname(next(iter(params.pickles.values())))
        
        # move the working directory into the test directory
        newDir = os.path.join(ResultsTest.TEST_DIR, os.path.basename(oldDir))
        
        # move the directory, but allow it to fail for loading parameters outside the tests
        try:
            os.rename(oldDir, newDir)
        except:
            pass
        
        # update the location of the pickle files
        for key in params.pickles.keys():
            params.pickles[key] = os.path.join(newDir, os.path.basename(params.pickles[key]))
        
        # update the location of the all_contigs.fna file
        params.allContigsFna = os.path.join(newDir, os.path.basename(params.allContigsFna))
        
        return params

    def _downloadOneGenome(ftp:str, fn:str) -> None:
        """downloads a genome using wget

        Args:
            ftp (str): the ftp path of the file to download
            fn (str): the file to save the genome

        Raises:
            RuntimeError: failed to download a genome
        """
        # constants
        CMD = ('wget', '-q', '-O')
        MAX_TRIES = 5
        SLEEP_SEC = 1
        
        # get the gzip filename
        gz = fn + ".gz"
        
        # build the command to run
        cmd = list(CMD)
        cmd.append(gz)
        cmd.append(ftp)
        
        # try to download the genome
        numTries = 0
        while numTries < MAX_TRIES:
            try:
                subprocess.run(cmd, check=True)
                break
            
            # keep trying; sleep between attempts
            except subprocess.CalledProcessError:
                numTries += 1
                time.sleep(SLEEP_SEC)
        
        # if the genome failed to be downloaded, then raise an error
        if numTries == MAX_TRIES:
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
        rowNum = -1
        
        # go through each line in the file
        with open(fn, 'r') as fh:
            for line in fh:
                # increment the row number
                rowNum += 1
                
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
                    out[(fwd,rev)] = Result(rowNum,ftm,fgc,rtm,rgc,additional)
        
        return out

    def _parseBedFile(fn:str) -> dict[tuple[Seq,int],dict[str,tuple[int,int,str]]]:
        """parses a file in BED format

        Args:
            fn (str): a file in BED format

        Returns:
            dict[tuple[Seq,int],dict[str,tuple[int,int,str]]]: key=(primer, pair num); val=dict:key=contig; val=(start position, end position, strand)
        """
        # constant
        SEP = "\t"
        
        # initialize output
        out = dict()
        
        # for each line in the file
        with open(fn, 'r') as fh:
            for line in fh:
                # parse the line
                contig,start,end,seq,num,strand = line.rstrip().split(SEP)
                
                # cast non-strings appropriately
                start = int(start)
                end = int(end)
                seq = Seq(seq)
                num = int(num)
                
                # save output
                out[(seq,num)] = out.get((seq,num), dict())
                out[(seq,num)][contig] = (start, end, strand)

        return out

    def _loadOneGenomeSequence(fn:str) -> dict[str,dict[str,Seq]]:
        """loads a single genome sequence into memory

        Args:
            fn (str): the filename (basename)

        Returns:
            dict[str,dict[str,Seq]]: key=contig name; val=dict: key=strand; val=sequence
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
            dict[str,dict[str,dict[str,Seq]]]: key=genome name; val=dict: key=contig name; val=dict: key=strand; val=sequence
        """
        return {os.path.basename(fn):ResultsTest._loadOneGenomeSequence(fn) for fn in ResultsTest.seqFiles}
    
    def _kmerListToDict(klist:list[tuple[str,int]]) -> dict[str,list[int]]:
        """converts a list of kmers to a dictionary

        Args:
            klist (list[tuple[str,int]]): a list of kmers paired with their start positions

        Returns:
            dict[str,list[int]]: key=kmer; val=list of start positions
        """
        # initialize output
        out = dict()
        
        # for each kmer
        for kmer,start in klist:
            # create the list (if necessary) and append the start position to it
            out[kmer] = out.get(kmer, list())
            out[kmer].append(start)
        
        return out
    
    def _getKmersForOneContig(seq:dict[str,Seq], minLen:int, maxLen:int, primers:set[str]) -> dict[str,dict[Seq,list[int]]]:
        """gets all the kmers for a single contig

        Args:
            seq (dict[str,Seq]): key=strand; val=contig sequence
            minLen (int): the minimum kmer length
            maxLen (int): the maximum kmer length
            primers (set[str]): a set of all the primers being evaluated

        Returns:
            dict[str,dict[Seq,list[int]]]: key=strand; val=dict: key=kmers; val=list of start positions
        """
        # initialize output
        out = {Primer.PLUS:  dict(),
               Primer.MINUS: dict()}
        
        # extract the sequences from the dictionary as well as its length
        fwd = str(seq[Primer.PLUS])
        rev = str(seq[Primer.MINUS])
        length = len(fwd)
        
        # for each kmer length
        for k in range(minLen, maxLen+1):
            # create countgraphs for each strand
            fkh = Countgraph(k, 1e7, 1)
            rkh = Countgraph(k, 1e7, 1)
            
            # get all the kmers and their start positions for each strand
            fmers = [(kmer, start) for start,kmer in enumerate(fkh.get_kmers(fwd)) if kmer in primers]
            rmers = [(kmer, length - start - 1) for start,kmer in enumerate(rkh.get_kmers(rev)) if kmer in primers]
            
            # convert the lists to dictionaries and save them
            out[Primer.PLUS].update(ResultsTest._kmerListToDict(fmers))
            out[Primer.MINUS].update(ResultsTest._kmerListToDict(rmers))
        
        return out

    def _getKmerBindingSites(results:dict[tuple[Seq,Seq],Result], sequences:dict[str,dict[str,dict[str,Seq]]], minLen:int, maxLen:int) -> dict[Seq,dict[str,dict[str,dict[str,list[int]]]]]:
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
        allPrimers = {str(p) for pair in results.keys() for p in pair}
        clock = Clock()
        
        # for each genome
        for name in sequences.keys():
            # print status
            clock.printStart(f'    getting kmer binding sites for {name}')
            
            # for each contig
            for contig in sequences[name].keys():
                # get the all the kmers
                kmers = ResultsTest._getKmersForOneContig(sequences[name][contig], minLen, maxLen, allPrimers)
                
                # for each kmer
                for kmer in allPrimers:
                    # build the data structure: seq, name, contig, strand, list of start positions
                    bindingSites[kmer] = bindingSites.get(kmer, dict())
                    bindingSites[kmer][name] = bindingSites[kmer].get(name, dict())
                    bindingSites[kmer][name][contig] = bindingSites[kmer][name].get(contig, dict())
                    bindingSites[kmer][name][contig][Primer.PLUS] = bindingSites[kmer][name][contig].get(Primer.PLUS, list())
                    bindingSites[kmer][name][contig][Primer.MINUS] = bindingSites[kmer][name][contig].get(Primer.MINUS, list())
                    
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
    
    def _noHomoDimers(seq:Seq, minTm:float, tempTolerance:float) -> bool:
        """evaluates sequences for the presence of homodimers

        Args:
            seq (Seq): a DNA sequence
            minTm (float): the minimum primer Tm
            tempTolerance (float): temp tolerance for secondary structures

        Returns:
            bool: indicates if the sequence does not form homodimers
        """
        # homodimers don't exist as long as the tm for formation is 5° below the min Tm
        return primer3.calc_homodimer_tm(str(seq)) < (minTm - tempTolerance)
    
    def _noHairpins(seq:Seq, minTm:float, tempTolerance:float) -> bool:
        """evaluates sequences for the presence of hairpins

        Args:
            seq (Seq): a DNA sequence
            minTm (float): the minimum primer Tm
            tempTolerance (float): temp tolerance for secondary structures

        Returns:
            bool: indicates if the sequence does not form hairpins
        """
        # hairpins don't exist as long as the tm for formation is 5° below the min Tm
        return primer3.calc_hairpin_tm(str(seq)) < (minTm - tempTolerance)
    
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
            self.assertTrue(ResultsTest._noHomoDimers(fwd, self.params.minTm, self.params.tempTolerance), f"{FAIL_MSG}{fwd}")
            self.assertTrue(ResultsTest._noHomoDimers(rev, self.params.minTm, self.params.tempTolerance), f"{FAIL_MSG}{rev}")
    
    def testJ_noHairpins(self) -> None:
        # constant
        FAIL_MSG = "hairpin present in "
        
        # for each pair, make sure there are no homodimers
        for fwd,rev in self.results.keys():
            self.assertTrue(ResultsTest._noHairpins(fwd, self.params.minTm, self.params.tempTolerance), f"{FAIL_MSG}{fwd}")
            self.assertTrue(ResultsTest._noHairpins(rev, self.params.minTm, self.params.tempTolerance), f"{FAIL_MSG}{rev}")

    def testK_noPrimerDimers(self) -> None:
        """do all pairs not produce primer dimers
        """
        for fwd,rev in self.results.keys():
            # create mock primers
            fwd = Primer(fwd, '', 100, len(fwd), Primer.PLUS)
            rev = Primer(rev, '', 100, len(rev), Primer.MINUS)

            # check dimer Tm
            forms,tm = _formsDimers(str(fwd),
                                    str(rev.reverseComplement()),
                                    fwd.Tm,
                                    rev.Tm,
                                    self.params.p3_mvConc,
                                    self.params.p3_dvConc,
                                    self.params.p3_dntpConc,
                                    self.params.p3_dnaConc,
                                    self.params.p3_tempC,
                                    self.params.p3_maxLoop,
                                    self.params.tempTolerance)
            
            # check for primer dimers; both primers need to be on the same strand
            self.assertFalse(forms, f'primer dimers present in {fwd}, {rev}')

    def testL_ingroupHasOneProduct(self) -> None:
        """do all pairs produce one product size in the ingroup
        """
        # constants
        FAIL_MSG_A = "multiple product sizes for "
        FAIL_MSG_B = " in "
        
        # for each pair in the ingroup
        for pair in self.results.keys():
            for name in ResultsTest.INGROUP_FILES:
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
            for name in ResultsTest.INGROUP_FILES:
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
            for name in ResultsTest.OUTGROUP_FILES:
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
            for name in ResultsTest.OUTGROUP_FILES:
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
            for name in ResultsTest.OUTGROUP_FILES:
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
            for name in ResultsTest.INGROUP_FILES:
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
            for name in ResultsTest.INGROUP_FILES:
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
            for name in ResultsTest.INGROUP_FILES:
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
            for name in ResultsTest.OUTGROUP_FILES:
                # get the list of pcr products and their contigs
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
    
    def testT_allPairsWritten(self) -> None:
        """verifies that all primer pairs were written to file
        """
        # unpickle the pairs
        with open(self.params.pickles[Parameters._PAIR_3], 'rb') as fh:
            pairs:dict[tuple[Primer,Primer]] = pickle.load(fh)
        
        # convert the pairs to a set of pairs as Seq
        expected = {(f.seq,r.seq) for f,r in pairs.keys()}
        
        # get the pairs saved in the results as a set
        observed = set(self.results.keys())
        
        # make sure that the pairs match
        self.assertEqual(len(expected), len(observed), f'expected {len(expected)} primer pairs but observed {len(observed)}')
        self.assertSetEqual(expected, observed, 'expected pairs do not match those written to file')
    
    def testU_arePairsSortedCorrectly(self) -> None:
        """verifies that the pairs were written in the correct order
        """
        # get the expected sort order as a tuples of Seq
        expected = [(f.seq,r.seq) for f,r in _sortPairs(self.params, self.pairs)]
        
        # get the observed order that was written
        observed  = sorted(self.results.keys(), key=lambda x: self.results[x].rowNum)
        
        # check that the order is conserved
        self.assertListEqual(expected, observed, 'expected sort order does not match the order written to file')
    
    def testV_primerTmsWereCorrect(self) -> None:
        """verifies that the homodimer Tms and hairpin Tms were stored correctly
        """
        # get a set of the non-redundant primer sequences
        primers = {x for p in self.pairs.keys() for x in p}
        
        # for each primer
        for primer in primers:
            # check hairpin Tms
            self.assertAlmostEqual(primer.hairpinTm, primer3.calc_hairpin_tm(str(primer)), places=3)
            self.assertAlmostEqual(primer.rcHairpin, primer3.calc_hairpin_tm(str(primer.reverseComplement())), places=3)
            
            # check homodimer Tms
            self.assertAlmostEqual(primer.homodimerTm, primer3.calc_homodimer_tm(str(primer)), places=3)
            self.assertAlmostEqual(primer.rcHomodimer, primer3.calc_homodimer_tm(str(primer.reverseComplement())), places=3)
    
    def testW_productTmsWereCorrect(self) -> None:
        """verifies that heterodimer Tm was correct for all PCR products
        """
        # for each pair
        for fwd,rev in self.pairs.keys():
            # calculate the heterodimer Tm for this pair
            dimerTm = primer3.calc_heterodimer_tm(str(fwd), str(rev))
            
            # for each product
            for name in self.pairs[(fwd,rev)].keys():
                product = self.pairs[(fwd,rev)][name]
                
                # check that the dimer Tm was saved correctly
                self.assertAlmostEqual(dimerTm, product.dimerTm, places=3)
    
    def testX_bedFileHasAllPrimers(self) -> None:
        """checks that the BED file has the expected sequences
        """
        # get the expected set of the non-redundant primer sequences
        expected = {x.seq for p in self.pairs.keys() for x in p}
        
        # get the primer sequences observed in the bed file
        observed = {x[0] for x in self.bed.keys()}
        
        self.assertSetEqual(observed, expected, "BED file primers do not match expected values")
    
    def testY_bedFileCorrect(self) -> None:
        """checks that the data saved in each row of the BED file is accurate"""
        # get the expected sort order of pairs as a tuples of strings
        sortOrder = [(str(f.seq),str(r.seq)) for f,r in _sortPairs(self.params, self.pairs)]
        
        # for each primer
        for primer,num in self.bed.keys():
            # check the pair number matches the sort order index
            self.assertIn(str(primer), sortOrder[num], f"unexpected sort order for {primer}: {num}")
            
            # for each contig
            for contig in self.bed[(primer,num)].keys():
                # extract data
                start, end, strand = self.bed[(primer,num)][contig]
                
                # migrate to the correct genome name
                for name in self.bindingSites[primer]:
                    if contig in self.bindingSites[primer][name].keys():
                        break
                
                # make sure the strand matches; an empty list indicates the wrong strand
                self.assertNotEqual(self.bindingSites[primer][name][contig][strand], [], f"strand mismatch for {primer} on {name}: {contig}")
                
                # end positions should be start + primer length if on the plus strand
                self.assertEqual(end, start + len(primer), f"end position is not correct for {primer} on {name}: {contig}")
                
                # the end will be the start position in the binding sites; it will be 1bp ahead of the actual site (exclusive ends)
                if strand == Primer.MINUS:
                    start = end - 1
                
                # start positions need to match the binding sites
                self.assertEqual(start, self.bindingSites[primer][name][contig][strand][0], f"start positions for {name}: {contig} do not match for {primer}")

# entrypoint
if __name__ == "__main__":
    unittest.main()
