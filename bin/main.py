import os
from Bio import SeqIO
from bin.Clock import Clock
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters
from bin.getCandidateKmers import _getAllCandidateKmers
from bin.removeOutgroupPrimers import _removeOutgroupPrimers
from bin.getPrimerPairs import _getPrimerPairs, _keepOnePairPerBinPair
from bin.analysis import _initializeAnalysisData, _updateAnalysisData, _plotAnalysisData


def __getCheckpoint(params:Parameters) -> tuple[bool,bool,bool,bool,bool]:
    sharedExists = os.path.exists(params.pickles[0])
    candidExists  = os.path.exists(params.pickles[1])
    unfiltExists = os.path.exists(params.pickles[2])
    filterExists = os.path.exists(params.pickles[3])
    analysExists = os.path.exists(params.pickles[4])
    
    return sharedExists, candidExists, unfiltExists, filterExists, analysExists
        


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


def _main(params:Parameters) -> None:
    """main runner function:
        * reads ingroup and outgroup sequences into memory
        * gets candidate kmers to use to search for primer pairs
        * gets primer pairs that are present in all ingroup
        * removes primer pairs that are present in the outgroup
        * writes data to file
    """
    # pickle constants
    CAND_NUM = 1
    PAIR_1_NUM = 2
    PAIR_2_NUM = 3
    ANAL_NUM = 4
    
    # debugging constants
    DONE = "done"
    
    # messages
    MSG_1A = "identifying kmers suitable for use as primers in all "
    MSG_1B = " ingroup genome sequences"
    MSG_X  = "initializing analysis data"
    MSG_2  = "identifying pairs of primers found in all ingroup sequences"
    MSG_J  = "updating analysis data"
    MSG_3A = "    identified "
    MSG_3B = " primer pairs shared in all ingroup sequences"
    MSG_4  = "keeping one primer pair per locus"
    MSG_5A = "writing "
    MSG_5B = " primer pairs to "
    MSG_6  = "plotting kmer distributions"
    MSG_7  = '\ntotal runtime: '
    
    # start the timers
    totalClock = Clock()
    clock = Clock()
    
    # read the ingroup sequences into memory
    ingroupSeqs = __readSequenceData(params.ingroupFns, params.format)
    
    # log data
    params.log.rename(_main.__name__)
    
    
    params.log.info(f"{MSG_1A}{len(ingroupSeqs)}{MSG_1B}")
    
    # get the candidate kmers for the ingroup
    clock.printStart(f"{MSG_1A}{len(ingroupSeqs)}{MSG_1B}", end='\n', spin=False)
    candidateKmers = _getAllCandidateKmers(ingroupSeqs, params, taskStart==0)
    clock.printDone()
    params.log.rename(_main.__name__)
    params.log.info(f"{DONE} {clock.getTimeString()}")
    
    # dump the candidates to file
    params.dumpObj(candidateKmers, params.pickles[CAND_NUM], "candidate kmers")
    
    # initialize the analysis data
    params.log.info(MSG_X)
    clock.printStart(MSG_X)
    analysisData = _initializeAnalysisData(candidateKmers)
    clock.printDone()
    params.log.info(f"{DONE} {clock.getTimeString()}")
    
    # dump the analysis to file
    params.dumpObj(analysisData, params.pickles[ANAL_NUM], "initial analysis data")
    
    # get the suitable primer pairs for the ingroup
    params.log.info(MSG_2)
    clock.printStart(MSG_2)
    pairs = _getPrimerPairs(candidateKmers, params)
    clock.printDone()
    params.log.info(f"{DONE} {clock.getTimeString()}")
    
    # remove the candidate kmers to free up memory
    del candidateKmers
    
    # dump the pairs to file
    params.dumpObj(pairs, params.pickles[PAIR_1_NUM], "unfiltered pairs")
    
    # update the analysis data
    params.log.info(MSG_J)
    clock.printStart(MSG_J)
    _updateAnalysisData(analysisData, pairs.keys())
    clock.printDone()
    params.log.info(f"{DONE} {clock.getTimeString()}")

    # dump the analysis to file
    params.dumpObj(analysisData, params.pickles[ANAL_NUM], 'updated analysis data')
    
    # print the number of candidate primer pairs
    print(f"{MSG_3A}{len(pairs)}{MSG_3B}")
    params.log.info(f"{MSG_3A}{len(pairs)}{MSG_3B}")
    
    # remove primer pairs that make products in the outgroup
    if params.outgroupFns != []:
        outgroupSeqs = __readSequenceData(params.outgroupFns, params.format)
        _removeOutgroupPrimers(outgroupSeqs, pairs, params)
        params.log.rename(_main.__name__)
        params.dumpObj(pairs, params.pickles[PAIR_2_NUM], "filtered pairs")
    
    # update analysis data
    params.log.info(MSG_J)
    clock.printStart(MSG_J)
    _updateAnalysisData(analysisData, pairs.keys())
    clock.printDone()
    params.log.info(f'{DONE} {clock.getTimeString()}')
    
    # dump the updated analysis
    params.dumpObj(analysisData, params.pickles[ANAL_NUM], 'updated analysis data')
    
    # only keep one pair per bin pair (only process ingroup bins)
    clock.printStart(MSG_4)
    params.log.info(MSG_4)
    for name in map(os.path.basename, params.ingroupFns):
        _keepOnePairPerBinPair(pairs, name)
    params.log.info(f"{DONE} {clock.getTimeString()}")
    
    # update analysis data
    params.log.info(MSG_J)
    clock.printStart(MSG_J)
    _updateAnalysisData(analysisData, pairs.keys())
    clock.printDone()
    params.log.info(f'{DONE} {clock.getTimeString()}')
    
    # write pair results to file
    params.log.info(f"{MSG_5A}{len(pairs)}{MSG_5B}{params.resultsFn}")
    clock.printStart(f"{MSG_5A}{len(pairs)}{MSG_5B}{params.resultsFn}")
    __writePrimerPairs(params.resultsFn, pairs)
    clock.printDone()
    params.log.info(f"{DONE} {clock.getTimeString()}")
    
    # save analysis data to file if in debug mode
    params.dumpObj(analysisData, params.pickles[ANAL_NUM], 'updated analysis data')
    
    # make analysis plots
    params.log.info(MSG_6)
    clock.printStart(MSG_6, end=' ...\n', spin=False)
    _plotAnalysisData(analysisData, params)
    clock.printDone()
    params.log.rename(_main.__name__)
    params.log.info(f"{DONE} {clock.getTimeString()}")
    
    # print the total runtime
    print(MSG_7, end='', flush=True)
    totalClock.printTime()
    params.log.info(f"{MSG_7} {totalClock.getTimeString()}")
