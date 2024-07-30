import os
from Bio import SeqIO
from bin.Clock import Clock
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters
from bin.getPrimerPairs import _getPrimerPairs
from bin.validatePrimers import _validatePrimerPairs
from bin.getCandidateKmers import _getAllCandidateKmers
from bin.removeOutgroupPrimers import _removeOutgroupPrimers

# global constants
__version__ = "1.3.2"
__author__ = "Joseph S. Wirth"
__SHARED_NUM = 0
__CAND_NUM   = 1
__PAIR_1_NUM = 2
__PAIR_2_NUM = 3
__PAIR_3_NUM = 4


def __getCheckpoint(params:Parameters) -> tuple[bool,bool,bool,bool,bool]:
    """gets booleans for checkpointing

    Args:
        params (Parameters): a Parameters object

    Returns:
        tuple[bool,bool,bool,bool,bool,bool]: each boolean indicates if the pickle exists: shared kmers; candidate kmers; unfiltered pairs; filtered pairs; validated pairs
    """
    # determine if each pickle exists
    sharedExists = os.path.exists(params.pickles[__SHARED_NUM])
    candidExists = os.path.exists(params.pickles[__CAND_NUM])
    unfiltExists = os.path.exists(params.pickles[__PAIR_1_NUM])
    filterExists = os.path.exists(params.pickles[__PAIR_2_NUM])
    validsExists = os.path.exists(params.pickles[__PAIR_3_NUM])
    
    return sharedExists, candidExists, unfiltExists, filterExists, validsExists


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
    
    # print status
    params.log.info(f"{MSG_1A}{len(params.ingroupFns)}{MSG_1B}")
    clock.printStart(f"{MSG_1A}{len(params.ingroupFns)}{MSG_1B}", end='\n', spin=False)
    
    # get the candidate kmers
    candidateKmers = _getAllCandidateKmers(params, sharedExists)
    
    # print status
    clock.printDone()
    params.log.rename(__getCandidates.__name__)
    params.log.info(f"done {clock.getTimeString()}")
    
    return candidateKmers


def __getUnfilteredPairs(params:Parameters, candidateKmers:dict[str,dict[str,list[Primer]]], clock:Clock) -> dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]:
    """gets the unfiltered pairs

    Args:
        params (Parameters): a Parameters object
        candidateKmers (dict[str,dict[str,list[Primer]]]): the dictionary produced by __getCandidates
        clock (Clock): a Clock object

    Returns:
        dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]: key=Primer pairs; val=dict: key=genome name; val=tuple: contig, pcr product size, bin pair (contig, bin1, bin2)
    """
    # messages
    MSG_1  = "identifying pairs of primers found in all ingroup sequences"
    MSG_2A = "    identified "
    MSG_2B = " primer pairs shared in all ingroup sequences"
    
    # log data
    params.log.rename(__getUnfilteredPairs.__name__)
    
    # get the suitable primer pairs for the ingroup
    params.log.info(MSG_1)
    clock.printStart(MSG_1)
    pairs = _getPrimerPairs(candidateKmers, params)
    clock.printDone()
    params.log.info(f"done {clock.getTimeString()}")
    
    # print the number of candidate primer pairs
    print(f"{MSG_2A}{len(pairs)}{MSG_2B}")
    params.log.info(f"{MSG_2A}{len(pairs)}{MSG_2B}")
    
    # dump the pairs to file
    params.dumpObj(pairs, params.pickles[__PAIR_1_NUM], "unfiltered pairs", prefix=" "*4)
    
    return pairs


def __removeOutgroup(params:Parameters, pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]) -> None:
    """removes outgroup primers from the pairs dictionary. does not return

    Args:
        params (Parameters): a Parameters object
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]): the dictionary produced by __getUnfilteredPairs
    """
    # messages
    GAP = " "*4
    MSG_1A = "removed "
    MSG_1B = " pairs present in the outgroup ("
    MSG_1C = " pairs remaining)"
    
    # log data
    params.log.rename(__removeOutgroup.__name__)
    
    # only have work to do if there are outgroup genomes
    if params.outgroupFns != []:
        # load genome sequences
        outgroupSeqs = __readSequenceData(params.outgroupFns, params.format)
        
        # determine the starting number of pairs
        startingNum = len(pairs)
        
        # remove the outgroup primers
        _removeOutgroupPrimers(outgroupSeqs, pairs, params)
        
        # move the log back
        params.log.rename(__removeOutgroup.__name__)
        
        # print the number of pairs removed
        finalNum = len(pairs)
        print(f"{GAP}{MSG_1A}{startingNum - finalNum}{MSG_1B}{finalNum}{MSG_1C}")
        params.log.info(f"{GAP}{MSG_1A}{startingNum - finalNum}{MSG_1B}{finalNum}{MSG_1C}")
        
        # dump the results
        params.dumpObj(pairs, params.pickles[__PAIR_2_NUM], "filtered pairs", prefix=GAP)


def __validatePairs(params:Parameters, pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]) -> None:
    """validates primer pairs with in silico PCR and removes invalid pairs

    Args:
        params (Parameters): a Parameters object
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]): the dictionary produced by __getUnfilteredPairs
    """
    # messages
    GAP = ' '*4
    MSG_1A = 'removed '
    MSG_1B = ' pairs not validated by isPcr ('
    MSG_1C = ' pairs remaining)'
    
    # get the starting number of pairs
    startingNum = len(pairs)
    
    # validate pairs
    _validatePrimerPairs(params, pairs)
    
    # move the log back
    params.log.rename(__validatePairs.__name__)
    
    # print and log the number of pairs removed and remaining
    finalNum = len(pairs)
    print(f'{GAP}{MSG_1A}{startingNum - finalNum}{MSG_1B}{finalNum}{MSG_1C}')
    params.log.info(f'{GAP}{MSG_1A}{startingNum - finalNum}{MSG_1B}{finalNum}{MSG_1C}')
    
    # dump the results
    params.dumpObj(pairs, params.pickles[__PAIR_3_NUM], "validated pairs", prefix=GAP)


def __writePrimerPairs(params:Parameters, pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]], clock:Clock) -> None:
    """writes pairs of primers to file

    Args:
        params (Parameters): the Parameters object for the run
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]])): key=primer pair; val=dict: key=genome name; val=tuple: contig; pcr len;  bin pair
        clock (Clock): a Clock object
    """
    # messages
    MSG_A = "writing "
    MSG_B = " primer pairs to "
    
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
    
    # print status
    params.log.rename(__writePrimerPairs.__name__)
    params.log.info(f"{MSG_A}{len(pairs)}{MSG_B}'{os.path.relpath(params.resultsFn)}'")
    clock.printStart(f"{MSG_A}{len(pairs)}{MSG_B}'{os.path.relpath(params.resultsFn)}'")
    
    # get the names of genomes, and create headers with the same order
    names = list(next(iter(pairs.values())).keys())
    headers = getHeaders(names)
    
    # open the file
    with open(params.resultsFn, 'w') as fh:
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
    
    # print status
    clock.printDone()
    params.log.info(f"done {clock.getTimeString()}")    


def __removeIntermediateFiles(params:Parameters) -> None:
    """removes intermediate files

    Args:
        params (Parameters): a Parameters object
    """
    # remove all the pickle files
    for fn in params.pickles.values():
        os.remove(fn)
    
    # remove the pickle directory
    os.rmdir(os.path.dirname(fn))
    
    # remove the isPcr files
    os.remove(params.allContigsFna)
    os.remove(params.queryFn)


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
    sharedExists, candExists, unfiltExists, filtExists, validExists = __getCheckpoint(params)
    
    # jump straight to the end if possible
    if validExists:
        # load data
        pairs = params.loadObj(params.pickles[__PAIR_3_NUM])
        
        # write data to file
        __writePrimerPairs(params, pairs, clock)
    
    # otherwise jump to validatoin if possible
    elif filtExists:
        # load data
        pairs = params.loadObj(params.pickles[__PAIR_2_NUM])
        
        # validate pairs
        __validatePairs(params, pairs)
        
        #  write data to file
        __writePrimerPairs(params, pairs, clock)
    
    # otherwise jump to outgroup removal if possible
    elif unfiltExists:
        # load data
        pairs = params.loadObj(params.pickles[__PAIR_1_NUM])
        
        # remove outgroup
        __removeOutgroup(params, pairs)
        
        # validate pairs
        __validatePairs(params, pairs)
        
        # write data to file
        __writePrimerPairs(params, pairs, clock)
    
    # otherwise jump to pair retrieval
    elif candExists:
        # load data
        candidateKmers = params.loadObj(params.pickles[__CAND_NUM])
        
        # get unfiltered pairs
        pairs = __getUnfilteredPairs(params, candidateKmers, clock)
        
        # remove the candidate kmers to free up memory
        del candidateKmers
        
        # remove outgroup
        __removeOutgroup(params, pairs)
        
        # validate pairs
        __validatePairs(params, pairs)
        
        # plot and write data
        __writePrimerPairs(params, pairs, clock)
    
    # otherwise run from the beginning (_getAllCandidates handles sharedExists)
    else:
        # get the candidate kmers for the ingroup
        candidateKmers = __getCandidates(params, sharedExists, clock)
        
        # get the pairs shared in the ingroup
        pairs = __getUnfilteredPairs(params, candidateKmers, clock)
        
        # remove the candidate kmers to free up memory
        del candidateKmers
        
        # remove primer pairs that make products in the outgroup
        __removeOutgroup(params, pairs)
        
        # validate pairs
        __validatePairs(params, pairs)
        
        # make plots and write data
        __writePrimerPairs(params, pairs, clock)
    
    # move the logger back
    params.log.rename(_runner.__name__)
    
    # remove the intermediate files unless keeping them
    if not params.keepIntermediateFiles:
        __removeIntermediateFiles(params)
    
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
            params.log.critical(e)
            
            Clock._killWheel()
            
            # terminate        
            raise Exception(e)
