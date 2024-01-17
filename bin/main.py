from Bio import SeqIO
import getopt, os, shutil, sys
from Bio.SeqRecord import SeqRecord
from bin.getPrimerPairs import _getPrimerPairs
from bin.getCandidateKmers import _getAllCandidateKmers


def __parseArgs() -> tuple[list[str],list[str],str,str,int,int,float,float,float,float,int,int,float,int,bool]:
    """parses command line arguments

    Raises:
        ValueError: invalid ingroup file(s)
        ValueError: invalid outgroup file(s)
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
        ValueError: must specify an ingroup file
        ValueError: must specify an output file

    Returns:
        tuple[list[str],list[str],str,str,int,int,float,float,float,float,int,int,float,int,bool,bool]:
            ingroupFns,outgroupFns,outFN,format,minPrimerLen,maxPrimerLen,minGc,maxGc,minTm,maxTm,minPcrLen,maxPcrLen,maxTmDiff,numThreads,keep,helpRequested
    """
    # constants
    ALLOWED_FORMATS = ('genbank', 'fasta')
    SEP = ","
    
    # flags
    INGROUP_FLAGS = ('-i', '--ingroup')
    OUT_FLAGS = ('-o', '--out')
    OUTGROUP_FLAGS = ('-u', '--outgroup')
    FMT_FLAGS = ('-f', '--format')
    PRIMER_LEN_FLAGS = ('-p', '--primer_len')
    GC_FLAGS = ('-g', '--gc_range')
    TM_FLAGS = ('-t', '--tm_range')
    THREADS_FLAGS = ('-n', '--num_threads')
    PCR_LEN_FLAGS = ('-r', '--pcr_prod_len')
    TM_DIFF_FLAGS = ('-d', '--tm_diff')
    KEEP_FLAGS = ('-k', '--keep')
    HELP_FLAGS = ('-h', '--help')
    SHORT_OPTS = INGROUP_FLAGS[0][-1] + ":" + \
                 OUT_FLAGS[0][-1] + ":" + \
                 OUTGROUP_FLAGS[0][-1] + ":" + \
                 FMT_FLAGS[0][-1] + ":" + \
                 PRIMER_LEN_FLAGS[0][-1] + ":" + \
                 GC_FLAGS[0][-1] + ":" + \
                 TM_FLAGS[0][-1] + ":" + \
                 PCR_LEN_FLAGS[0][-1] + ":" + \
                 TM_DIFF_FLAGS[0][-1] + ":" + \
                 THREADS_FLAGS[0][-1] + ":" + \
                 KEEP_FLAGS[0][-1] +\
                 HELP_FLAGS[0][-1]
    LONG_OPTS = (INGROUP_FLAGS[1][2:] + "=",
                 OUT_FLAGS[1][2:] + "=",
                 OUTGROUP_FLAGS[1][2:] + "=",
                 FMT_FLAGS[1][2:] + "=",
                 PRIMER_LEN_FLAGS[1][2:] + "=",
                 GC_FLAGS[1][2:] + "=",
                 TM_FLAGS[1][2:] + "=",
                 PCR_LEN_FLAGS[1][2:] + "=",
                 TM_DIFF_FLAGS[1][2:] + "=",
                 THREADS_FLAGS[1][2:] + "=",
                 KEEP_FLAGS[1][2:],
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
    DEF_KEEP = False

    # messages
    IGNORE_MSG = 'ignoring unused argument: '
    ERR_MSG_1  = 'invalid or missing ingroup file(s)'
    ERR_MSG_2  = 'invalid or missing outgroup file(s)'
    ERR_MSG_3  = 'invalid format'
    ERR_MSG_4  = 'can only specify one primer length or a range (min,max)'
    ERR_MSG_5  = 'primer lengths are not integers'
    ERR_MSG_6  = 'must specify a range of GC values (min,max)'
    ERR_MSG_7  = 'gc values are not numeric'
    ERR_MSG_8  = 'must specify a range of Tm values (min, max)'
    ERR_MSG_9  = 'Tm values are not numeric'
    ERR_MSG_10 = 'can only specify one PCR product length or a range (min,max)'
    ERR_MSG_11 = 'PCR product lengths are not integers'
    ERR_MSG_12 = 'max Tm difference is not numeric'
    ERR_MSG_13 = 'num threads is not an integer'
    ERR_MSG_14 = 'must specify one or more ingroup files'
    ERR_MSG_15 = 'must specify an output file'

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
                   GAP + f"{INGROUP_FLAGS[0] + SEP_1 + INGROUP_FLAGS[1]:<22}{'[file(s)] ingroup filename(s); comma-separated list'}" + EOL + \
                   GAP + f"{OUT_FLAGS[0] + SEP_1 + OUT_FLAGS[1]:<22}{'[file] output filename'}" + EOL*2 + \
                   "optional arguments:" + EOL + \
                   GAP + f"{OUTGROUP_FLAGS[0] + SEP_1 + OUTGROUP_FLAGS[1]:<22}{'[file(s)] outgroup filename(s); comma-separated list'}" + EOL + \
                   GAP + f"{FMT_FLAGS[0] + SEP_1 + FMT_FLAGS[1]:<22}{'[str] file format of the ingroup and outgroup'}{{{ALLOWED_FORMATS[0] + SEP_2 + ALLOWED_FORMATS[1]}}}{DEF_OPEN + DEF_FRMT + CLOSE}" + EOL + \
                   GAP + f"{PRIMER_LEN_FLAGS[0] + SEP_1 + PRIMER_LEN_FLAGS[1]:<22}{'[int(s)] a single primer length or a range specified as '}" + "'min,max'" + f"{DEF_OPEN + str(DEF_MIN_LEN) + SEP_3 + str(DEF_MAX_LEN) + CLOSE}" + EOL + \
                   GAP + f"{GC_FLAGS[0] + SEP_1 + GC_FLAGS[1]:<22}{'[float,float] a min and max percent GC specified as a comma separated list' + DEF_OPEN + str(DEF_MIN_GC) + SEP_3 + str(DEF_MAX_GC) + CLOSE}" + EOL + \
                   GAP + f"{TM_FLAGS[0] + SEP_1 + TM_FLAGS[1]:<22}{'[float,float] a min and max melting temp (Tm) specified as a comma separated list' + DEF_OPEN + str(DEF_MIN_TM) + SEP_3 + str(DEF_MAX_TM) + CLOSE}" + EOL + \
                   GAP + f"{PCR_LEN_FLAGS[0] + SEP_1 + PCR_LEN_FLAGS[1]:<22}{'[int(s)] a single PCR product length or a range specified as '}" + "'min,max'" + f"{DEF_OPEN + str(DEF_MIN_PCR) + SEP_3 + str(DEF_MAX_PCR) + CLOSE}" + EOL + \
                   GAP + f"{TM_DIFF_FLAGS[0] + SEP_1 + TM_DIFF_FLAGS[1]:<22}{'[float] the maximum allowable Tm difference between a pair of primers' + DEF_OPEN + str(DEF_MAX_TM_DIFF) + CLOSE}" + EOL + \
                   GAP + f"{THREADS_FLAGS[0] + SEP_1 + THREADS_FLAGS[1]:<22}{'[int] the number of threads for parallel processing' + DEF_OPEN + str(DEF_NUM_THREADS) + CLOSE}" + EOL + \
                   GAP + f"{KEEP_FLAGS[0] + SEP_1 + KEEP_FLAGS[1]:<22}{'keep intermediate files' + DEF_OPEN + str(DEF_KEEP) + CLOSE}" + EOL + \
                   GAP + f"{HELP_FLAGS[0] + SEP_1 + HELP_FLAGS[1]:<22}{'print this message'}" + EOL*2
        
        print(HELP_MSG)
        
    # set default values
    ingroupFns = None
    outFN = None
    outgroupFns = list()
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
    keep = DEF_KEEP
    helpRequested = False
    
    # give help if requested
    if HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv or len(sys.argv) == 1:
        helpRequested = True
        printHelp()
    
    # parse command line arguments
    else:
        opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
        for opt,arg in opts:
            # get the ingroup filenames
            if opt in INGROUP_FLAGS:
                arg = arg.split(SEP)
                for fn in arg:
                    if not os.path.isfile(fn):
                        raise ValueError(ERR_MSG_1)
                ingroupFns = arg
            
            # get output filehandle
            elif opt in OUT_FLAGS:
                outFN = arg
            
            # get the outgroup filenames
            elif opt in OUTGROUP_FLAGS:
                arg = arg.split(SEP)
                for fn in arg:
                    if not os.path.isfile(fn):
                        raise ValueError(ERR_MSG_2)
                outgroupFns = arg
            
            # get the file format
            elif opt in FMT_FLAGS:
                if arg not in ALLOWED_FORMATS:
                    raise ValueError(ERR_MSG_3)
                frmt = arg
            
            # get the primer lengths
            elif opt in PRIMER_LEN_FLAGS:
                # split comma-separated list
                primerRange = arg.split(SEP)
                
                # make sure at one or two primers specified
                if len(primerRange) not in {1,2}:
                    raise ValueError(ERR_MSG_4)
                
                # coerce to lengths to ints
                try:
                    primerRange = [int(x) for x in primerRange]
                except:
                    raise ValueError(ERR_MSG_5)
                
                # save values
                minLen = min(primerRange)
                maxLen = max(primerRange)
            
            # get the allowed GC range
            elif opt in GC_FLAGS:
                # expecting two values separated by a comma
                gcRange = arg.split(SEP)
                if len(gcRange) != 2:
                    raise ValueError(ERR_MSG_6)
                
                # make sure the values are numeric
                try:
                    gcRange = [float(x) for x in gcRange]
                except:
                    raise ValueError(ERR_MSG_7)
            
                # save values
                minGc = min(gcRange)
                maxGc = max(gcRange)
            
            # get the allowed Tm range
            elif opt in TM_FLAGS:
                # expecting two values separated by a comma
                tmRange = arg.split(SEP)
                if len(tmRange) != 2:
                    raise ValueError(ERR_MSG_8)
            
                # make sure the values are numeric
                try:
                    tmRange = [float(x) for x in tmRange]
                except:
                    raise ValueError(ERR_MSG_9)
            
                # save values
                minTm = min(tmRange)
                maxTm = max(tmRange)
        
            # get the allowed PCR lengths
            elif opt in PCR_LEN_FLAGS:
                # expecting one or two values
                pcrRange = arg.split(SEP)
                if len(pcrRange) not in {1,2}:
                    raise ValueError(ERR_MSG_10)
                
                # coerce to integers
                try:
                    pcrRange = [int(x) for x in pcrRange]
                except:
                    raise ValueError(ERR_MSG_11)
            
                # save values
                minPcr = min(pcrRange)
                maxPcr = max(pcrRange)
            
            # get the allowed Tm difference between primer pairs
            elif opt in TM_DIFF_FLAGS:
                # make sure input is numeric
                try:
                    maxTmDiff = float(arg)
                except:
                    raise ValueError(ERR_MSG_12)
            
            # get the number of threads to use
            elif opt in THREADS_FLAGS:
                # make sure input is an integer
                try:
                    numThreads = int(arg)
                except:
                    raise ValueError(ERR_MSG_13)
            
            # determine if intermediate files should be kept
            elif opt in KEEP_FLAGS:
                keep = True
            
            else:
                print(IGNORE_MSG + opt + " " + arg)
        
        # make sure an input file was specified
        if ingroupFns is None:
            raise ValueError(ERR_MSG_14)
        
        # make sure an output file was specified
        if outFN is None:
            raise ValueError(ERR_MSG_15)
    
    return ingroupFns,outgroupFns,outFN,frmt,minLen,maxLen,minGc,maxGc,minTm,maxTm,minPcr,maxPcr,maxTmDiff,numThreads,keep,helpRequested


def __readSequenceData(seqFiles:list[str], frmt:str) -> dict[str, list[SeqRecord]]:
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


def _main() -> None:
    """main runner function:
        * reads ingroup and outgroup sequences into memory
        * gets candidate kmers to use to search for primer pairs
        * gets primer pairs that are present in all ingroup and excluded from all outgroup
    """
    # constant
    TEMP_DIR = os.path.join(".", "kmerwd")

    # parse command line arguments
    ingroupFiles,outgroupFiles,outFn,frmt,minPrimerLen,maxPrimerLen,minGc,maxGc,minTm,maxTm,minPcrLen,maxPcrLen,maxTmDiff,numThreads,keep,helpRequested = __parseArgs()
    
    # only do work if help was not requested
    if not helpRequested:
        # read the ingroup and outgroup sequences into memory
        ingroupSeqs = __readSequenceData(ingroupFiles, frmt)
        outgroupSeqs = __readSequenceData(outgroupFiles, frmt)
        
        # get the candidate kmers for the ingroup
        candidateKmers = _getAllCandidateKmers(ingroupSeqs, outgroupSeqs, minPrimerLen, maxPrimerLen, minGc, maxGc, minTm, maxTm, numThreads)
        
        # create the temporary directory
        if not os.path.exists(TEMP_DIR):
            os.mkdir(TEMP_DIR)
        
        # get the suitable primer pairs for the ingroup
        _getPrimerPairs(candidateKmers, minPrimerLen, minPcrLen, maxPcrLen, maxTmDiff, outFn, numThreads, TEMP_DIR)

        # remove the temporary directory unless user asked to keep it
        if not keep:
            shutil.rmtree(TEMP_DIR)
