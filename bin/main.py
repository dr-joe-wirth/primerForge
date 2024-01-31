__version__ = "0.2.0"
__author__ = "Joseph S. Wirth"

import os
from Bio import SeqIO
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters
from bin.getPrimerPairs import _getPrimerPairs
from bin.Clock import Clock, _printStart, _printDone
from bin.getCandidateKmers import _getAllCandidateKmers
from bin.removeOutgroupPrimers import _removeOutgroupPrimers


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
    # debugging constants
    CAND_FN = "candidates.p"
    PAIR_1_FN = "pairs.p"
    PAIR_2_FN = "pairs_noOutgroup.p"
    DONE = "done"
    
    # messages
    MSG_1A = "identifying kmers suitable for use as primers in all "
    MSG_1B = " ingroup genome sequences"
    MSG_2 = "identifying pairs of primers found in all ingroup sequences"
    MSG_3A = "    identified "
    MSG_3B = " primer pairs shared in all ingroup sequences"
    MSG_4  = "removing primer pairs present in the outgroup sequences"
    MSG_5A = "writing "
    MSG_5B = " primer pairs to "
    MSG_6  = '\ntotal runtime: '
    
    # start the timers
    totalClock = Clock()
    clock = Clock()
    
    # only do work if help was not requested
    if not params.helpRequested:
        if params.debug:
            params.log.setLogger(_main.__name__)

        # read the ingroup sequences into memory
        ingroupSeqs = __readSequenceData(params.ingroupFns, params.format)

        # get the candidate kmers for the ingroup
        _printStart(clock, f"{MSG_1A}{len(ingroupSeqs)}{MSG_1B}", '\n')
        if params.debug: params.log.writeInfoMsg(f"{MSG_1A}{len(ingroupSeqs)}{MSG_1B}")
        candidateKmers = _getAllCandidateKmers(ingroupSeqs, params)
        _printDone(clock)
        
        # save the candidates to file if in debug mode
        if params.debug:
            params.log.setLogger(_main.__name__)
            params.log.writeInfoMsg(f"{DONE} {clock.getTime()}")
            params.dumpObj(candidateKmers, CAND_FN, "candidate kmers")
        
        # get the suitable primer pairs for the ingroup
        _printStart(clock, MSG_2)
        if params.debug: params.log.writeInfoMsg(MSG_2)
        pairs = _getPrimerPairs(candidateKmers, params)
        _printDone(clock)
        
        # remove the candidate kmers to free up memory
        del candidateKmers
        
        # save the pairs to file if in debug mode
        if params.debug:
            params.log.writeInfoMsg(f"{DONE} {clock.getTime()}")
            params.dumpObj(pairs, PAIR_1_FN, "unfiltered pairs")
        
        # print the number of candidate primer pairs
        print(f"{MSG_3A}{len(pairs)}{MSG_3B}")
        if params.debug: params.log.writeInfoMsg(f"{MSG_3A}{len(pairs)}{MSG_3B}")
        
        # remove primer pairs that make products in the outgroup
        if params.outgroupFns != []:
            _printStart(clock, MSG_4)
            if params.debug: params.log.writeInfoMsg(MSG_4)
            outgroupSeqs = __readSequenceData(params.outgroupFns, params.format)
            _removeOutgroupPrimers(outgroupSeqs, pairs, params)   ######## what is that final parameter??
            _printDone(clock)
            
            # save pairs to file if in debug mode
            if params.debug:
                params.log.writeInfoMsg(f"{DONE} {clock.getTime()}")
                params.dumpObj(pairs, PAIR_2_FN, "filterd pairs")
        
        # write results to file
        _printStart(clock, f"{MSG_5A}{len(pairs)}{MSG_5B}{params.outFn}")
        if params.debug: params.log.writeInfoMsg(f"{MSG_5A}{len(pairs)}{MSG_5B}{params.outFn}")
        __writePrimerPairs(params.outFn, pairs)
        if params.debug: params.log.writeInfoMsg(f"{DONE} {clock.getTime()}")
        _printDone(clock)
        
        # print the total runtime
        print(MSG_6, end='', flush=True)
        totalClock.printTime()
        if params.debug: params.log.writeInfoMsg(f"{MSG_6} {totalClock.getTime()}")
