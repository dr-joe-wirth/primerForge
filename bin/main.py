import os
from Bio import SeqIO
from Bio.Seq import Seq
from bin.Clock import Clock
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters
from bin.AnalysisData import AnalysisData
from bin.getCandidateKmers import _getAllCandidateKmers
from bin.removeOutgroupPrimers import _removeOutgroupPrimers
from bin.getPrimerPairs import _getPrimerPairs, _keepOnePairPerBinPair
from bin.analysis import _initializeAnalysisData, _updateAnalysisData, _plotAnalysisData

# global constants
__version__ = "0.7.4"
__author__ = "Joseph S. Wirth"
__SHARED_NUM = 0
__CAND_NUM   = 1
__PAIR_1_NUM = 2
__PAIR_2_NUM = 3
__PAIR_3_NUM = 4
__ANAL_NUM   = 5


def __getCheckpoint(params:Parameters) -> tuple[bool,bool,bool,bool,bool,bool]:
    """gets booleans for checkpointing

    Args:
        params (Parameters): a Parameters object

    Returns:
        tuple[bool,bool,bool,bool,bool,bool]: each boolean indicates if the pickle exists: shared kmers; candidate kmers; unfiltered pairs; filtered pairs; final pairs; analysis data
    """
    
    # determine if each pickle exists
    sharedExists = os.path.exists(params.pickles[__SHARED_NUM])
    candidExists = os.path.exists(params.pickles[__CAND_NUM])
    unfiltExists = os.path.exists(params.pickles[__PAIR_1_NUM])
    filterExists = os.path.exists(params.pickles[__PAIR_2_NUM])
    finalsExists = os.path.exists(params.pickles[__PAIR_3_NUM])
    analysExists = os.path.exists(params.pickles[__ANAL_NUM])
    
    return sharedExists, candidExists, unfiltExists, filterExists, finalsExists, analysExists


def __readSequenceData(seqFiles:list[str], frmt:str) -> dict[str, list[SeqRecord]]:
    """reads sequence data into file

    Args:
        seqFiles (list[str]): a list of sequence files to read
        frmt (str): the format of the sequence files

    Returns:
        dict[str, list[SeqRecord]]: key=genome name; val=a sequence iterator; functionally a list of SeqRecords
    """
    # initialize output
    out = dict()
    
    # for each file in the list
    for fn in seqFiles:
        # get the genome name and use it as a key to store the list of parsed contigs
        name = os.path.basename(fn)
        out[name] = SeqIO.parse(fn, frmt)
    
    return out


def __getCandidates(params:Parameters, sharedExists:bool, clock:Clock) -> dict[str,dict[str,list[Primer]]]:
    """gets the candidate kmers

    Args:
        params (Parameters): a Parameters object
        sharedExists (bool): indicates if the shared kmers pickle exists
        clock (Clock): a Clock object

    Returns:
        dict[str,dict[str,list[Primer]]]: key=genome name; val=dict: key=contig name; val=list of Primers
    """
    # messages
    MSG_1A = "identifying kmers suitable for use as primers in all "
    MSG_1B = " ingroup genome sequences"
    
    # log data
    params.log.rename(__getCandidates.__name__)
    
    # read the ingroup sequences into memory
    ingroupSeqs = __readSequenceData(params.ingroupFns, params.format)
    
    # print status
    params.log.info(f"{MSG_1A}{len(ingroupSeqs)}{MSG_1B}")
    clock.printStart(f"{MSG_1A}{len(ingroupSeqs)}{MSG_1B}", end='\n', spin=False)
    
    # get the candidate kmers
    candidateKmers = _getAllCandidateKmers(ingroupSeqs, params, sharedExists)
    
    # print status
    clock.printDone()
    params.log.rename(__getCandidates.__name__)
    params.log.info(f"done {clock.getTimeString()}")
    
    # dump the candidates to file
    params.dumpObj(candidateKmers, params.pickles[__CAND_NUM], "candidate kmers")
    
    return candidateKmers


def __initializeAnalysis(params:Parameters, candidateKmers:dict[str,dict[str,list[Primer]]], clock:Clock) -> dict[tuple[Seq,str],AnalysisData]:
    """initialize the analysis data

    Args:
        params (Parameters): a Parameters object
        candidateKmers (dict[str,dict[str,list[Primer]]]): the dictionary produced by __getCandidates
        clock (Clock): a Clock object

    Returns:
        dict[tuple[Seq,str],AnalysisData]: key=(kmer sequence, contig name); val=AnalysisData object
    """
    
    # message
    MSG  = "initializing analysis data"
    
    # log data
    params.log.rename(__initializeAnalysis.__name__)
    
    # print status
    params.log.info(MSG)
    clock.printStart(MSG)
    
    # initialize analysis data
    analysisData = _initializeAnalysisData(candidateKmers)
    
    # print status
    clock.printDone()
    params.log.info(f"done {clock.getTimeString()}")
    
    # dump the analysis to file
    params.dumpObj(analysisData, params.pickles[__ANAL_NUM], "initial analysis data")
    
    return analysisData


def __getUnfilteredPairs(params:Parameters, candidateKmers:dict[str,dict[str,list[Primer]]], analysisData:dict[tuple[Seq,str],AnalysisData], clock:Clock) -> dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]:
    """gets the unfiltered pairs

    Args:
        params (Parameters): a Parameters object
        candidateKmers (dict[str,dict[str,list[Primer]]]): the dictionary produced by __getCandidates
        analysisData (dict[tuple[Seq,str],AnalysisData]): the dictionary produced by __initializeAnalysis
        clock (Clock): a Clock object

    Returns:
        dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]: key=Primer pairs; val=dict: key=genome name; val=tuple: contig, pcr product size, bin pair (contig, bin1, bin2)
    """
    # messages
    MSG_1  = "identifying pairs of primers found in all ingroup sequences"
    MSG_2  = "updating analysis data"
    MSG_3A = "    identified "
    MSG_3B = " primer pairs shared in all ingroup sequences"
    
    # log data
    params.log.rename(__getUnfilteredPairs.__name__)
    
    # get the suitable primer pairs for the ingroup
    params.log.info(MSG_1)
    clock.printStart(MSG_1)
    pairs = _getPrimerPairs(candidateKmers, params)
    clock.printDone()
    params.log.info(f"done {clock.getTimeString()}")
    
    # remove the candidate kmers to free up memory
    del candidateKmers
    
    # dump the pairs to file
    params.dumpObj(pairs, params.pickles[__PAIR_1_NUM], "unfiltered pairs")
    
    # update the analysis data
    params.log.info(MSG_2)
    clock.printStart(MSG_2)
    _updateAnalysisData(analysisData, pairs.keys())
    clock.printDone()
    params.log.info(f"done {clock.getTimeString()}")
    
    # dump the analysis to file
    params.dumpObj(analysisData, params.pickles[__ANAL_NUM], 'updated analysis data')
    
    # print the number of candidate primer pairs
    print(f"{MSG_3A}{len(pairs)}{MSG_3B}")
    params.log.info(f"{MSG_3A}{len(pairs)}{MSG_3B}")
    
    return pairs


def __removeOutgroup(params:Parameters, pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]], analysisData:dict[tuple[Seq,str],AnalysisData], clock:Clock) -> None:
    """removes outgroup primers from the pairs dictionary. does not return

    Args:
        params (Parameters): a Parameters object
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]): the dictionary produced by __getUnfilteredPairs
        analysisData (dict[tuple[Seq,str],AnalysisData]): the dictionary produced by __initializeAnalysis
        clock (Clock): a Clock object
    """
    # message
    MSG = "updating analysis data"
    
    # log data
    params.log.rename(__removeOutgroup.__name__)
    
    # only have work to do if there are outgroup genomes
    if params.outgroupFns != []:
        # load genome sequences
        outgroupSeqs = __readSequenceData(params.outgroupFns, params.format)
        
        # remove the outgroup primers
        _removeOutgroupPrimers(outgroupSeqs, pairs, params)
        
        # dump the results
        params.log.rename(__removeOutgroup.__name__)
        params.dumpObj(pairs, params.pickles[__PAIR_2_NUM], "filtered pairs")
    
        # print status
        params.log.info(MSG)
        clock.printStart(MSG)
        
        # update analysis data
        _updateAnalysisData(analysisData, pairs.keys())
        
        # print status
        clock.printDone()
        params.log.info(f'done {clock.getTimeString()}')
    
        # dump the updated analysis
        params.dumpObj(analysisData, params.pickles[__ANAL_NUM], 'updated analysis data')


def __getFinalPairs(params:Parameters, pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]], analysisData:dict[tuple[Seq,str],AnalysisData], clock:Clock) -> None:
    """gets the final pairs by keeping one per bin pair. does not return

    Args:
        params (Parameters): a Parameters object
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]): the dictionary produced by __getUnfilteredPairs
        analysisData (dict[tuple[Seq,str],AnalysisData]): the dictionary produced by __initializeAnalysis
        clock (Clock): a Clock object
    """
    # messages
    MSG_1  = "keeping one primer pair per locus"
    MSG_2 = "updating analysis data"
    
    # print status
    params.log.rename(__getFinalPairs.__name__)
    params.log.info(MSG_1)
    clock.printStart(MSG_1)
    
    # only keep one pair per bin pair (only process ingroup bins)
    for name in map(os.path.basename, params.ingroupFns):
        _keepOnePairPerBinPair(pairs, name)
    
    # print status
    clock.printDone()
    params.log.info(f"done {clock.getTimeString()}")
    
    # print status
    params.log.info(MSG_2)
    clock.printStart(MSG_2)
    
    # update analysis data (twice if no outgroup)
    _updateAnalysisData(analysisData, pairs.keys())
    if params.outgroupFns == []:
        _updateAnalysisData(analysisData, pairs.keys())
    
    # print status
    clock.printDone()
    params.log.info(f'done {clock.getTimeString()}')
    
    # dump analysis data to file
    params.dumpObj(analysisData, params.pickles[__ANAL_NUM], 'updated analysis data')
    
    # dump pairs to file
    params.dumpObj(pairs, params.pickles[__PAIR_3_NUM], 'final pairs')


def __writePrimerPairs(fn:str, pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]) -> None:
    """writes pairs of primers to file

    Args:
        fn (str): the filename to write
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]])): key=primer pair; val=dict: key=genome name; val=tuple: contig; pcr len;  bin pair
    """
    # contants
    EOL = "\n"
    SEP = "\t"
    NUM_DEC = 1
    
    # helper function to create the headers
    def getHeaders(names) -> list[str]:
        # constants
        HEADERS = ('fwd_seq',
                   'fwd_Tm',
                   'fwd_GC',
                   'rev_seq',
                   'rev_Tm',
                   'rev_GC')
        CONTIG = "_contig"
        LENGTH = "_length"
        
        # each name will have a contig and a length
        headers = list(HEADERS)
        for name in names:
            headers.append(name + CONTIG)
            headers.append(name + LENGTH)
        
        return headers
    
    # get the names of genomes, and create headers with the same order
    names = list(next(iter(pairs.values())).keys())
    headers = getHeaders(names)
    
    # open the file
    with open(fn, 'w') as fh:
        # write the headers
        fh.write(SEP.join(headers) + EOL)
        fh.flush()
        
        # for each primer pair
        for fwd,rev in pairs.keys():
            # save the primer pair data
            row = [fwd.seq,
                   round(fwd.Tm, NUM_DEC),
                   round(fwd.gcPer, NUM_DEC),
                   rev.seq,
                   round(rev.Tm, NUM_DEC),
                   round(rev.gcPer, NUM_DEC)]
            
            # then save the contig name and PCR product length for each genome
            for name in names:
                contig = pairs[(fwd,rev)][name][0]
                pcrLen = pairs[(fwd,rev)][name][1]
                row.append(contig)
                row.append(pcrLen)
            
            # write the data to the file
            fh.write(f"{SEP.join(map(str, row))}{EOL}")
            fh.flush()


def __plotAndWrite(params:Parameters, pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]], analysisData:dict[tuple[Seq,str],AnalysisData], clock:Clock) -> None:
    """writes and plots final data

    Args:
        params (Parameters): a Parameters object
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]): the dictionary produced by __getUnfilteredPairs
        analysisData (dict[tuple[Seq,str],AnalysisData]): the dictionary produced by __initializeAnalysis
        clock (Clock): a Clock object
    """
    # messages
    MSG_1A = "writing "
    MSG_1B = " primer pairs to "
    MSG_2  = "plotting kmer distributions"
    
    # print status
    params.log.rename(__plotAndWrite.__name__)
    params.log.info(f"{MSG_1A}{len(pairs)}{MSG_1B}{params.resultsFn}")
    clock.printStart(f"{MSG_1A}{len(pairs)}{MSG_1B}{params.resultsFn}")
    
    # write primer pairs to file
    __writePrimerPairs(params.resultsFn, pairs)
    
    # print status
    clock.printDone()
    params.log.info(f"done {clock.getTimeString()}")
    params.log.info(MSG_2)
    clock.printStart(MSG_2, end=' ...\n', spin=False)
    
    # plot analysis data
    _plotAnalysisData(analysisData, params)
    
    # print status
    clock.printDone()
    params.log.rename(__plotAndWrite.__name__)
    params.log.info(f"done {clock.getTimeString()}")


def __removePickles(params:Parameters) -> None:
    """removes pickle files

    Args:
        params (Parameters): a Parameters object
    """
    # remove all the pickle files
    for fn in params.pickles.values():
        os.remove(fn)
    
    # remove the pickle directory
    os.rmdir(os.path.dirname(fn))


def _runner(params:Parameters) -> None:
    """main runner function:
        * reads ingroup and outgroup sequences into memory
        * gets candidate kmers to use to search for primer pairs
        * gets primer pairs that are present in all ingroup
        * removes primer pairs that are present in the outgroup
        * writes data to file
    """
    # messages
    MSG  = 'total runtime: '
    
    # start the timers
    totalClock = Clock()
    clock = Clock()
    
    # start the logger
    params.log.rename(_runner.__name__)
    
    # get the checkpoint status
    sharedExists, candExists, unfiltExists, filtExists, finalExists, analExists = __getCheckpoint(params)
    
    # jump straight to the end if possible
    if finalExists and analExists:
        # load data
        pairs = params.loadObj(params.pickles[__PAIR_3_NUM])
        analysisData = params.loadObj(params.pickles[__ANAL_NUM])
        
        # plot and write data
        __plotAndWrite(params, pairs, analysisData, clock)
    
    # otherwise jump to final pair retrieval if possible
    elif filtExists and analExists:
        # load data
        pairs = params.loadObj(params.pickles[__PAIR_2_NUM])
        analysisData = params.loadObj(params.pickles[__ANAL_NUM])
        
        # get final pairs
        __getFinalPairs(params, pairs, analysisData, clock)
        
        # plot and write data
        __plotAndWrite(params, pairs, analysisData, clock)
    
    # otherwise jump to outgroup removal if possible
    elif unfiltExists and analExists:
        # load data
        pairs = params.loadObj(params.pickles[__PAIR_1_NUM])
        analysisData = params.loadObj(params.pickles[__ANAL_NUM])
        
        # remove outgroup
        __removeOutgroup(params, pairs, analysisData, clock)
        
        # get final pairs
        __getFinalPairs(params, pairs, analysisData, clock)
        
        # plot and write data
        __plotAndWrite(params, pairs, analysisData, clock)
    
    # otherwise jump to pair retrieval
    elif candExists:
        # load data
        candidateKmers = params.loadObj(params.pickles[__CAND_NUM])
        
        # get analysis data
        if analExists:
            analysisData = params.loadObj(params.pickles[__ANAL_NUM])
        else:
            analysisData = __initializeAnalysis(params, candidateKmers, clock)
        
        # get unfiltered pairs
        pairs = __getUnfilteredPairs(params, candidateKmers, analysisData, clock)
        
        # remove outgroup
        __removeOutgroup(params, pairs, analysisData, clock)
        
        # get final pairs
        __getFinalPairs(params, pairs, analysisData, clock)
        
        # plot and write data
        __plotAndWrite(params, pairs, analysisData, clock)
    
    # otherwise run from the beginning (_getAllCandidates handles sharedExists)
    else:
        # get the candidate kmers for the ingroup
        candidateKmers = __getCandidates(params, sharedExists, clock)
        
        # initialize the analysis data
        analysisData = __initializeAnalysis(params, candidateKmers, clock)
        
        # get the pairs shared in the ingroup
        pairs = __getUnfilteredPairs(params, candidateKmers, analysisData, clock)
        
        # remove primer pairs that make products in the outgroup
        __removeOutgroup(params, pairs, analysisData, clock)
        
        # keep one pair per bin pair
        __getFinalPairs(params, pairs, analysisData, clock)
        
        # make plots and write data
        __plotAndWrite(params, pairs, analysisData, clock)
    
    # move the logger back
    params.log.rename(_runner.__name__)
    
    # remove the pickles unless keeping them
    if not params.keepPickles:
        __removePickles(params)
    
    # print the total runtime
    print(f"\n{MSG}{totalClock.getTimeString()}\n")
    params.log.info(f"{MSG}{totalClock.getTimeString()}\n\n")


def main() -> None:
    """main runner function for primerForge

    Raises:
        Exception: catches all downstream exceptions
    """
    # parse command line arguments
    params = Parameters(__author__, __version__)
    
    if not params.helpRequested:
        try:
            # start up the logger
            params.log.rename(__name__)
            params.logRunDetails()
            
            # run primerForge
            _runner(params)
        
        # catch all error messages
        except Exception as e:
            # save the error message if in debug mode
            params.log.rename(__name__)
            params.log.critical(e)

            # terminate        
            raise Exception(e)
