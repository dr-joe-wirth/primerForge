import getopt, os, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# define a few global constants
PLUS_STRAND = "+"
MINUS_STRAND = "-"

def parseArgs() -> tuple[str,str,str,int,int,float,float,float,float,int,int,float,int,bool]:
    """parses command line arguments

    Raises:
        ValueError: invalid input file
        ValueError: invalid file format
        ValueError: invalid primer length
        ValueError: primer lengths not integers
        ValueError: specify range of GC
        ValueError: GC values non-numeric
        ValueError: specify range of Tm
        ValueError: Tm values non-numeric
        ValueError: invalid PCR product length
        ValueError: PCR product lengths not integers
        ValueError: Tm difference is non-numeric
        ValueError: num threads not an integer
        ValueError: must specify an input file

    Returns:
        tuple[str,str,str,int,int,float,float,float,float,int,int,float,int,bool]:
            inFN,outFN,format,minPrimerLen,maxPrimerLen,minGc,maxGc,minTm,maxTm,minPcrLen,maxPcrLen,maxTmDiff,numThreads,helpRequested
    """
    # constants
    ALLOWED_FORMATS = ('genbank', 'fasta')
    SEP = ","
    
    # flags
    IN_FLAGS = ('-i', '--in')
    OUT_FLAGS = ('-o', '--out')
    FMT_FLAGS = ('-f', '--format')
    PRIMER_LEN_FLAGS = ('-p', '--primer_len')
    GC_FLAGS = ('-g', '--gc_range')
    TM_FLAGS = ('-t', '--tm_range')
    THREADS_FLAGS = ('-n', '--num_threads')
    PCR_LEN_FLAGS = ('-r', '--pcr_prod_len')
    TM_DIFF_FLAGS = ('-d', '--tm_diff')
    HELP_FLAGS = ('-h', '--help')
    SHORT_OPTS = IN_FLAGS[0][-1] + ":" + \
                 OUT_FLAGS[0][-1] + ":" + \
                 FMT_FLAGS[0][-1] + ":" + \
                 PRIMER_LEN_FLAGS[0][-1] + ":" + \
                 GC_FLAGS[0][-1] + ":" + \
                 TM_FLAGS[0][-1] + ":" + \
                 PCR_LEN_FLAGS[0][-1] + ":" + \
                 TM_DIFF_FLAGS[0][-1] + ":" + \
                 THREADS_FLAGS[0][-1] + ":" + \
                 HELP_FLAGS[0][-1]
    LONG_OPTS = (IN_FLAGS[1][2:] + "=",
                 OUT_FLAGS[1][2:] + "=",
                 FMT_FLAGS[1][2:] + "=",
                 PRIMER_LEN_FLAGS[1][2:] + "=",
                 GC_FLAGS[1][2:] + "=",
                 TM_FLAGS[1][2:] + "=",
                 PCR_LEN_FLAGS[1][2:] + "=",
                 TM_DIFF_FLAGS[1][2:] + "=",
                 THREADS_FLAGS[1][2:] + "=",
                 HELP_FLAGS[1][2:])

    # default values
    DEF_FRMT = ALLOWED_FORMATS[0]
    DEF_MIN_LEN = 16
    DEF_MAX_LEN = 20
    DEF_MIN_GC = 40.0
    DEF_MAX_GC = 60.0
    DEF_MIN_TM = 55.0
    DEF_MAX_TM = 68.0
    DEF_MIN_PCR = 120
    DEF_MAX_PCR = 2400
    DEF_MAX_TM_DIFF = 5.0
    DEF_NUM_THREADS = 1

    # messages
    IGNORE_MSG = 'ignoring unused argument: '
    ERR_MSG_1  = 'invalid or missing input file'
    ERR_MSG_2  = 'invalid format'
    ERR_MSG_3  = 'can only specify one primer length or a range (min,max)'
    ERR_MSG_4  = 'primer lengths are not integers'
    ERR_MSG_5  = 'must specify a range of GC values (min,max)'
    ERR_MSG_6  = 'gc values are not numeric'
    ERR_MSG_7  = 'must specify a range of Tm values (min, max)'
    ERR_MSG_8  = 'Tm values are not numeric'
    ERR_MSG_9  = 'can only specify one PCR product length or a range (min,max)'
    ERR_MSG_10 = 'PCR product lengths are not integers'
    ERR_MSG_11 = 'max Tm difference is not numeric'
    ERR_MSG_12 = 'num threads is not an integer'
    ERR_MSG_13 = 'must specify an input file'
    ERR_MSG_14 = 'must specify an output file'

    def printHelp():
        GAP = " "*4
        EOL = "\n"
        SEP_1 = ", "
        SEP_2 = "|"
        SEP_3 = ","
        DEF_OPEN = ' (default: '
        CLOSE = ')'
        HELP_MSG = EOL + "Finds pairs of primers suitable for an input genome." + EOL + \
                   GAP + "Joseph S. Wirth, 2023" + EOL*2 + \
                   "usage:" + EOL + \
                   GAP + "primerDesign.py [-iofpgtnrdh]" + EOL*2 + \
                   "required arguments:" + EOL + \
                   GAP + f"{IN_FLAGS[0] + SEP_1 + IN_FLAGS[1]:<22}{'[file] input filename'}" + EOL + \
                   GAP + f"{OUT_FLAGS[0] + SEP_1 + OUT_FLAGS[1]:<22}{'[file] output filename'}" + EOL*2 + \
                   "optional arguments:" + EOL + \
                   GAP + f"{FMT_FLAGS[0] + SEP_1 + FMT_FLAGS[1]:<22}{'[str] input file format '}{{{ALLOWED_FORMATS[0] + SEP_2 + ALLOWED_FORMATS[1]}}}{DEF_OPEN + DEF_FRMT + CLOSE}" + EOL + \
                   GAP + f"{PRIMER_LEN_FLAGS[0] + SEP_1 + PRIMER_LEN_FLAGS[1]:<22}{'[int(s)] a single primer length or a range specified as '}" + "'min,max'" + f"{DEF_OPEN + str(DEF_MIN_LEN) + SEP_3 + str(DEF_MAX_LEN) + CLOSE}" + EOL + \
                   GAP + f"{GC_FLAGS[0] + SEP_1 + GC_FLAGS[1]:<22}{'[float,float] a min and max percent GC specified as a comma separated list' + DEF_OPEN + str(DEF_MIN_GC) + SEP_3 + str(DEF_MAX_GC) + CLOSE}" + EOL + \
                   GAP + f"{TM_FLAGS[0] + SEP_1 + TM_FLAGS[1]:<22}{'[float,float] a min and max melting temp (Tm) specified as a comma separated list' + DEF_OPEN + str(DEF_MIN_TM) + SEP_3 + str(DEF_MAX_TM) + CLOSE}" + EOL + \
                   GAP + f"{PCR_LEN_FLAGS[0] + SEP_1 + PCR_LEN_FLAGS[1]:<22}{'[int(s)] a single PCR product length or a range specified as '}" + "'min,max'" + f"{DEF_OPEN + str(DEF_MIN_PCR) + SEP_3 + str(DEF_MAX_PCR) + CLOSE}" + EOL + \
                   GAP + f"{TM_DIFF_FLAGS[0] + SEP_1 + TM_DIFF_FLAGS[1]:<22}{'[float] the maximum allowable Tm difference between a pair of primers' + DEF_OPEN + str(DEF_MAX_TM_DIFF) + CLOSE}" + EOL + \
                   GAP + f"{THREADS_FLAGS[0] + SEP_1 + THREADS_FLAGS[1]:<22}{'[int] the number of threads for parallel processing' + DEF_OPEN + str(DEF_NUM_THREADS) + CLOSE}" + EOL + \
                   GAP + f"{HELP_FLAGS[0] + SEP_1 + HELP_FLAGS[1]:<22}{'print this message'}" + EOL*2
        
        print(HELP_MSG)
        
    # set default values
    inFN = None
    outFN = None
    frmt = DEF_FRMT
    minLen = DEF_MIN_LEN
    maxLen = DEF_MAX_LEN
    minGc = DEF_MIN_GC
    maxGc = DEF_MAX_GC
    minTm = DEF_MIN_TM
    maxTm = DEF_MAX_TM
    minPcr = DEF_MIN_PCR
    maxPcr = DEF_MAX_PCR
    maxTmDiff = DEF_MAX_TM_DIFF
    numThreads = DEF_NUM_THREADS
    helpRequested = False
    
    # give help if requested
    if HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv or len(sys.argv) == 1:
        helpRequested = True
        printHelp()
    
    # parse command line arguments
    else:
        opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
        for opt,arg in opts:
            # get the input file
            if opt in IN_FLAGS:
                if not os.path.isfile(arg):
                    raise ValueError(ERR_MSG_1)
                inFN = arg
            
            # get output filehandle
            elif opt in OUT_FLAGS:
                outFN = arg
            
            # get the input file format
            elif opt in FMT_FLAGS:
                if arg not in ALLOWED_FORMATS:
                    raise ValueError(ERR_MSG_2)
                frmt = arg
            
            # get the primer lengths
            elif opt in PRIMER_LEN_FLAGS:
                # split comma-separated list
                primerRange = arg.split(SEP)
                
                # make sure at one or two primers specified
                if len(primerRange) not in {1,2}:
                    raise ValueError(ERR_MSG_3)
                
                # coerce to lengths to ints
                try:
                    primerRange = [int(x) for x in primerRange]
                except:
                    raise ValueError(ERR_MSG_4)
                
                # save values
                minLen = min(primerRange)
                maxLen = max(primerRange)
            
            # get the allowed GC range
            elif opt in GC_FLAGS:
                # expecting two values separated by a comma
                gcRange = arg.split(SEP)
                if len(gcRange) != 2:
                    raise ValueError(ERR_MSG_5)
                
                # make sure the values are numeric
                try:
                    gcRange = [float(x) for x in gcRange]
                except:
                    raise ValueError(ERR_MSG_6)
            
                # save values
                minGc = min(gcRange)
                maxGc = max(gcRange)
            
            # get the allowed Tm range
            elif opt in TM_FLAGS:
                # expecting two values separated by a comma
                tmRange = arg.split(SEP)
                if len(tmRange) != 2:
                    raise ValueError(ERR_MSG_7)
            
                # make sure the values are numeric
                try:
                    tmRange = [float(x) for x in tmRange]
                except:
                    raise ValueError(ERR_MSG_8)
            
                # save values
                minTm = min(tmRange)
                maxTm = max(tmRange)
        
            # get the allowed PCR lengths
            elif opt in PCR_LEN_FLAGS:
                # expecting one or two values
                pcrRange = arg.split(SEP)
                if len(pcrRange) not in {1,2}:
                    raise ValueError(ERR_MSG_9)
                
                # coerce to integers
                try:
                    pcrRange = [int(x) for x in pcrRange]
                except:
                    raise ValueError(ERR_MSG_10)
            
                # save values
                minPcr = min(pcrRange)
                maxPcr = max(pcrRange)
            
            # get the allowed Tm difference between primer pairs
            elif opt in TM_DIFF_FLAGS:
                # make sure input is numeric
                try:
                    maxTmDiff = float(arg)
                except:
                    raise ValueError(ERR_MSG_11)
            
            # get the number of threads to use
            elif opt in THREADS_FLAGS:
                # make sure input is an integer
                try:
                    numThreads = int(arg)
                except:
                    raise ValueError(ERR_MSG_12)
            
            else:
                print(IGNORE_MSG + opt + " " + arg)
        
        # make sure an input file was specified
        if inFN is None:
            raise ValueError(ERR_MSG_13)
        
        # make sure an output file was specified
        if outFN is None:
            raise ValueError(ERR_MSG_14)
    
    return inFN,outFN,frmt,minLen,maxLen,minGc,maxGc,minTm,maxTm,minPcr,maxPcr,maxTmDiff,numThreads,helpRequested


def readSequenceData(seqFiles:list[str], frmt:str) -> dict[str, list[SeqRecord]]:
    """reads sequence data into file

    Args:
        seqFiles (list[str]): a list of sequence files to read
        frmt (str): the format of the sequence files

    Returns:
        dict[str, list[SeqRecord]]: key=genome name; val=list of contigs as SeqRecords
    """
    # initialize output
    out = dict()
    
    # for each file in the list
    for fn in seqFiles:
        # get the genome name and use it as a key to store the list of parsed contigs
        name = os.path.splitext(os.path.basename(fn))[0]
        out[name] = list(SeqIO.parse(fn, frmt))
    
    return out


def kmpSearch(text:str, pattern:str) -> bool:
    """implementation of the KMP string search algorithm: O(n+m)

    Args:
        text (str): text to search
        pattern (str): pattern to search for

    Returns:
        bool: indicates if the pattern is found in the text
    """

    # helper function
    def computeLpsArray(patternLen:int) -> list[int]:
        """precomputes the longest prefix that is also a suffix array

        Args:
            patternLen (int): the length of the pattern

        Returns:
            list[int]: the precomputed lps array
        """
        # initialize the lps array
        lps = [0] * patternLen
        
        # initialize the indices; lps[0] always == 0 so start idx at 1
        matchLen = 0 
        idx = 1
        
        # go through each position in the pattern 
        while idx < patternLen:
            # if a match occurs update the match length, store in the array, and move to the next position
            if pattern[idx] == pattern[matchLen]:
                matchLen += 1
                lps[idx] = matchLen
                idx += 1
            
            # if a mismatch
            else:
                # if at previous string also did not match, then no lps; move to next position
                if matchLen == 0:
                    lps[idx] = 0
                    idx += 1 
                
                # if the previous string matched, we need to start at the beginning of the last match
                else:
                    matchLen = lps[matchLen - 1]
        
        return lps

    # get the length of the text and the pattern
    textLen = len(text)
    patternLen = len(pattern)
    
    lps = computeLpsArray(patternLen)
    
    # initialize indices
    idx = 0
    jdx = 0
    
    # keep evaluating the text until at the end
    while idx < textLen:
        # keep comparing strings while they match
        if text[idx] == pattern[jdx]:
            idx += 1
            jdx += 1
        
        # if a mismatch occurs
        else:
            # if at the start of the pattern, cannot reset; move to next text character
            if jdx == 0:
                idx += 1
            
            # otherwise, move to the previous character as defined by the lps array
            else:
                jdx = lps[jdx -1]
    
        # match found if at the end of the pattern
        if jdx == patternLen:
            return True
    
    return False


def getUniqueKmers(name:str, seqs:list[SeqRecord], minLen:int, maxLen:int, sharedDict) -> None:
    """designed to run in parallel. gets a collection of kmers that are:
        * not repeated anywhere in the genome
        * at least one end is a G or C

    Args:
        name (str): the name of the sequence
        seqs (list[SeqRecord]): a list of contigs as SeqRecord objects
        minLen (int): the minimum kmer length
        maxLen (int): the maximum kmer length
        sharedDict (DictProxy): a shared dictionary for parallel processing
    
    Returns:
        does not return; saves result in the shared dictionary
    """
    # helper functions to clarify boolean expressions
    
    def isOneEndGc(seq:Seq) -> bool:
        """ evaluates if one end of a sequence is a GC
        """
        # constants
        GC = {"G", "C", 'g', 'c'}
        
        return seq[-1] in GC or seq[0] in GC

    # initialize variables
    kmers = dict()
    kmers[PLUS_STRAND] = dict()
    kmers[MINUS_STRAND] = dict()
    bad = set()
    
    # go through each contig
    for contig in seqs:
        fwdSeq:Seq = contig.seq
        revSeq:Seq = fwdSeq.reverse_complement()
        
        # go through each kmer length
        for primerLen in range(minLen, maxLen+1):
            # get every possible kmer start position
            for start in range(len(contig)- (primerLen - 1)):
                # extract the kmer sequence
                fwdKmer = fwdSeq[start:start+primerLen]
                revKmer = revSeq[-(start+primerLen):-start]
                
                # mark duplicate primers for removal
                if fwdKmer in kmers[PLUS_STRAND].keys() or fwdKmer in kmers[MINUS_STRAND].keys():
                    bad.add(fwdKmer)
                if revKmer in kmers[PLUS_STRAND].keys() or revKmer in kmers[MINUS_STRAND].keys():
                    bad.add(revKmer)
                
                # only save primers that have GC at one end, no long repeats, and no complements
                if isOneEndGc(fwdSeq):
                    kmers[PLUS_STRAND][fwdKmer] = (contig.id, start, primerLen)
                    kmers[MINUS_STRAND][revKmer] = (contig.id, start, primerLen)

    # discard any sequences that were duplicated
    for seq in bad:
        try: kmers[PLUS_STRAND].pop(seq)
        except KeyError: pass
        
        try: kmers[MINUS_STRAND].pop(seq)
        except KeyError: pass
    
    sharedDict[name] = kmers


def getAllKmers(contig:SeqRecord, minLen:int, maxLen:int) -> set[Seq]:
    """gets all of the kmers in a given contig

    Args:
        contig (SeqRecord): the contig as a SeqRecord object
        minLen (int): the minimum kmer length
        maxLen (int): the maximum kmer length

    Returns:
        set[Seq]: a set of kmers
    """
    # initialize output
    out = set()
    
    # for each allowed primer length
    for primerLen in range(minLen, maxLen+1):
        # get every possible primer start position
        for start in range(len(contig)- (primerLen - 1)):
            # extract the primer sequence and save it
            seq:Seq = contig.seq[start:start+primerLen]
            out.add(seq)
            out.add(seq.reverse_complement())

    return out
