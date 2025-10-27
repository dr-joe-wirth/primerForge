from Bio import SeqIO
from Bio.Seq import Seq
from typing import Iterator
from bin.Clock import Clock
from bin.Primer import Primer
from ahocorasick import Automaton
import multiprocessing, os, primer3
from bin.Parameters import Parameters


# functions
def __getAllKmers(seq:str, k:int) -> Iterator[str]:
    """gets all possible kmers from a sequence

    Args:
        seq (str): the sequence to process
        k (int): the size of the kmers

    Yields:
        Iterator[str]: individual kmers
    """
    for i in range(len(seq) -k + 1):
        yield seq[i:i+k]


def __getAllowedPlusStrandKmers(fn:str, frmt:str, k:int) -> tuple[dict[str,dict[str,tuple[str,int,str]]], Automaton]:
    """gets all the allowed kmers from only the plus strand

    Args:
        fn (str): the filename of a genome
        frmt (str): the filename format
        k (int): the kmer size

    Returns:
        tuple[dict[str,dict[str,tuple[str,int,str]]], Automaton]: {kmer: {name: (contig, start, strand)}}, Aho-Corasick automaton
    """
    # helper function to identify suitable kmers
    def isSuitableKmer(seq:str) -> bool:
        """determines if a kmer is suitable

        Args:
            seq (str): the kmer sequence to evaluate

        Returns:
            bool: indicates suitability
        """
        GC = frozenset(("G", "C", 'g', 'c'))

        # check two conditions
        isOneEndGc  = seq[-1] in GC or seq[0] in GC
        palindromic = seq == Seq(seq).reverse_complement()

        return isOneEndGc and not palindromic

    # initialize outputs
    kmers = dict()
    auto = Automaton()

    # save the name
    name = os.path.basename(fn)
    
    # initialize some sets to store kmers
    bad = set()
    good = set()

    # categorize each kmer
    for contig in SeqIO.parse(fn, frmt):
        for kmer in __getAllKmers(str(contig.seq), k):
            if isSuitableKmer(kmer):
                # bad kmers appear more than once
                if kmer in good:
                    bad.add(kmer)
                    good.remove(kmer)
                    
                # good kmers are not bad
                elif kmer not in bad:
                    good.add(kmer)
            
    # don't need bad kmers anymore
    del bad

    # build the automaton from the good kmers
    for kmer in good:
        auto.add_word(kmer, kmer)
    auto.make_automaton()

    # extract the coordinate data for this genome
    for contig in SeqIO.parse(fn, frmt):
        for end,kmer in auto.iter(str(contig.seq)):
            kmers[kmer] = {name: (contig.id, end-k+1, Primer.PLUS)}

    return kmers, auto


def __updateAllowedKmers(fn:str, frmt:str, auto:Automaton, k:int, kmers:dict[str,dict[str,tuple[str,int,str]]]) -> Automaton:
    """updates the kmer dictionary by removing kmers that appears more than once

    Args:
        fn (str): the filename of the genome to evaluate
        frmt (str): the file format of the genome
        auto (Automaton): an Aho-Corasick automaton containing kmers of interest
        k (int): kmer size
        kmers (dict[str,dict[str,tuple[str,str,int]]]): {kmer: {genome: (contig, start, strand)}}

    Returns:
        Automaton: an updated Aho-Corasick automaton with only the kmers that passed filtering
    """
    def evaluateOneStrand(s:str, strand:str) -> None:
        """evaluates the kmers on a single strand

        Args:
            s (str): the sequence
            strand (str): the strand
        """
        # for each kmer appearance
        for end,kmer in auto.iter(s):
            # flag kmers that have been seen more than once
            if fn in kmers[kmer].keys():
                bad.add(kmer)
                del kmers[kmer]
            
            # save the coordinate for a kmers when seen the first time
            elif kmer not in bad:
                # determine the start position
                if strand == Primer.PLUS:
                    start = end - k + 1
                elif strand == Primer.MINUS:
                    start = len(s) - end - 1
                
                # store the coordinates
                kmers[kmer][name] = (contig.id, start, strand)
            
    # initialize a collection of bad kmers and get the name
    bad = set()
    name = os.path.basename(fn)

    # evaluate all kmers on both strands
    for contig in SeqIO.parse(fn, frmt):
        evaluateOneStrand(str(contig.seq), Primer.PLUS)
        evaluateOneStrand(str(contig.seq.reverse_complement()), Primer.MINUS)

    del bad

    # replace the existing automaton with a new one
    auto = Automaton()

    # filter kmers
    for kmer in list(kmers.keys()):        
        # kmers that did not appear in this genome must be removed
        if name not in kmers[kmer].keys():
            del kmers[kmer]
        
        # all other kmers go into the automaton
        else:
            auto.add_word(kmer, kmer)
    
    # finalize automaton
    auto.make_automaton()

    return auto


def __getSharedKmersOneK(arguments:tuple[list[str],str,int]) -> dict[str,dict[str,tuple[str,int,str]]]:
    """gets the shared kmers from genomes given a single kmer size
    designed for parallel processing with imap

    Args:
        arguments (tuple[list[str],str,int]): ingroup files, file format, kmer size

    Returns:
        dict[str,dict[str,tuple[str,str,int]]]: {kmer: {name: (contig, start, strand)}}
    """
    # parse the arguments into its individual components
    ingroupFns,frmt,k = arguments

    # get the kmers and an automaton for the smallest (first) genome
    kmers,auto = __getAllowedPlusStrandKmers(ingroupFns[0], frmt, k)

    # for each remaining genome, update the allowed kmer
    for fn in ingroupFns[1:]:
        auto = __updateAllowedKmers(fn, frmt, auto, k, kmers)
    
    # don't need the automaton anymore
    del auto

    return kmers


def __getSharedKmers(params:Parameters) -> dict[str,dict[str,tuple[str,int,str]]]:
    """retrieves all the kmers that are shared between the input genomes

    Args:
        params (Parameters): a Parameters object

    Returns:
        dict[str,dict[str,tuple[str,int,str]]]: key=kmer; val=dict: key=genome name: val=tuple: contig, start, strand
    """
    # helper function to generate arguments for __getSharedKmersOneK
    def genArgs() -> Iterator[tuple[list[str],str,int]]:
        for k in range(params.minLen, params.maxLen + 1):
            yield (params.ingroupFns, params.format, k)
        
    # initialize output
    out = dict()

    # # open the pool
    # with multiprocessing.Pool(params.numThreads) as pool:
    #     # save results as they become available
    #     for result in pool.imap_unordered(__getSharedKmersOneK, genArgs()):
    #         out.update(result)

    for k in range(params.minLen, params.maxLen + 1):
        out.update(__getSharedKmersOneK((params.ingroupFns, params.format, k)))
    
    return out


def __reorganizeDataByPosition(name:str, kmers:dict[str,dict[str,tuple[str,int,str]]]) -> dict[str,dict[int,list[tuple[str,str]]]]:
    """reorganizes data from __getSharedKmers by its genomic position

    Args:
        name (str): the name of the genome to be processed
        kmers (dict[str,dict[str,tuple[str,int,str]]]): the dictionary produced by __getSharedKmers

    Returns:
        dict[str,dict[int,list[tuple[str,str]]]]: key=contig; val=dict: key=start position; val=list of tuples: kmer, strand
    """
    # initialize output
    out = dict()
    
    # for each kmer
    for kmer in kmers.keys():
        # extract data from the dictionary
        contig, start, strand = kmers[kmer][name]
        
        # contig = top level key; start position = second level key; val = list
        out[contig] = out.get(contig, dict())
        out[contig][start] = out[contig].get(start, list())
        
        # add the sequence and its length to the list
        out[contig][start].append((kmer, strand))
    
    return out


def __removeRedundantKmerGroups(positions:dict[str,dict[int,list[tuple[str,str]]]], seen:set[tuple[str]]) -> None:
    """removes redundant groups of kmers from the positions dictionary

    Args:
        positions (dict[str,dict[int,list[tuple[str,str]]]]): the dictionary produced by __reorganizeDataByPosition
        seen (set[tuple[str]]): a set of kmer groups that have already been processed
    """
    # for each contig
    for contig in positions.keys():
        # get a list of kmer positions to allow on-the-fly deleting
        starts = list(positions[contig].keys())
        
        # for each start position
        for start in starts:
            # extract the group of kmers for this position
            group = tuple(sorted(x[0] for x in positions[contig][start]))
            
            # remove this from the dictionary if its kmer group has been seen
            if group in seen:
                del positions[contig][start]
            
            # mark previously unseen groups as seen
            else:
                seen.add(group)


def __evaluateKmersAtOnePosition(contig:str, start:int, positions:list[tuple[str,str]], minGc:float, maxGc:float, minTm:float, maxTm:float, \
                                 mvConc:float, dvConc:float, dntpConc:float, dnaConc:float, tempC:float, maxLoop:int, tempTolerance:float, repeats:Automaton) -> Primer:
    """evaluates all the primers at a single position in the genome; designed for parallel calls

    Args:
        contig (str): the name of the contig
        start (int): the start position in the sequence
        positions (list[tuple[str,int]]): a list of positions (kmer, strand)
        minGc (float): the minimum percent GC allowed
        maxGc (float): the maximum percent GC allowed
        minTm (float): the minimum melting temperature allowed
        maxTm (float): the maximum melting temperature allowed
        mvConc (float): primer3 mv_conc
        dvConc (float): primer3 dv_conc
        dntpConc (float): primer3 dntp_conc
        dnaConc (float): primer3 dna_conc
        tempC (float): primer3 temp_c
        maxLoop (int): primer3 max_loop
        tempTolerance (float): the minimum degrees below primer Tm allowed for secondary structures
        repeats (Automaton): an Aho-Corasick Automaton consisting of only the disallowed homopolymers
    
    Returns:
        Primer: a suitable primer at the given position
    """
    # define helper functions to make booleans below more readable
    def isGcWithinRange(primer:Primer) -> bool:
        """is the percent GC within the acceptable range?"""
        return primer.gcPer >= minGc and primer.gcPer <= maxGc

    def isTmWithinRange(primer:Primer) -> bool:
        """is the Tm within the acceptable range?"""
        return primer.Tm >= minTm and primer.Tm <= maxTm
    
    def noLongRepeats(primer:Primer) -> bool:
        """verifies that a primer does not have long repeats
        """
        # check for any repeats in the primer
        try:
            # this will work if a repeat is present
            next(iter(repeats.iter(str(primer))))
            return False
        
        # absence of repeats will not be iterable
        except StopIteration:
            return True

    def noHairpins(primer:Primer) -> bool:
        """verifies that the primer does not form hairpins
        """
        # save the hairpin Tms
        primer.hairpinTm = primer3.calc_hairpin_tm(str(primer),
                                                   mv_conc=mvConc,
                                                   dv_conc=dvConc,
                                                   dntp_conc=dntpConc,
                                                   dna_conc=dnaConc,
                                                   temp_c=tempC,
                                                   max_loop=maxLoop)
        primer.rcHairpin = primer3.calc_hairpin_tm(str(primer.reverseComplement()),
                                                   mv_conc=mvConc,
                                                   dv_conc=dvConc,
                                                   dntp_conc=dntpConc,
                                                   dna_conc=dnaConc,
                                                   temp_c=tempC,
                                                   max_loop=maxLoop)
        
        # hairpin tm should be less than (minTm - 5°); need to check both strands
        fwdOk = primer.hairpinTm < (minTm - tempTolerance)
        revOk = primer.rcHairpin < (minTm - tempTolerance)
        
        return fwdOk and revOk

    def noHomodimers(primer:Primer) -> bool:
        """verifies that the primer does not form homodimers
        """
        # save the homodimer Tms
        primer.homodimerTm = primer3.calc_homodimer_tm(str(primer),
                                                       mv_conc=mvConc,
                                                       dv_conc=dvConc,
                                                       dntp_conc=dntpConc,
                                                       dna_conc=dnaConc,
                                                       temp_c=tempC,
                                                       max_loop=maxLoop)
        primer.rcHomodimer = primer3.calc_homodimer_tm(str(primer.reverseComplement()),
                                                       mv_conc=mvConc,
                                                       dv_conc=dvConc,
                                                       dntp_conc=dntpConc,
                                                       dna_conc=dnaConc,
                                                       temp_c=tempC,
                                                       max_loop=maxLoop)
        
        # homodimer tm should be less than (minTm - 5°); need to check both strands
        fwdOk = primer.homodimerTm < (minTm - tempTolerance)
        revOk = primer.rcHomodimer < (minTm - tempTolerance)
        
        return fwdOk and revOk
    
    # initialize values for the while loop
    idx = 0
    
    # continue to iterate through each primer in the list until a primer is found 
    for idx in range(len(positions)):
        # extract data from the list
        seq,strand = positions[idx]
        
        # create a Primer object
        primer = Primer(seq, contig, start, len(seq), strand)
        
        # evaluate the primer's percent GC, Tm, hairpin potential, and homodimer potential; save if passes
        if isGcWithinRange(primer) and isTmWithinRange(primer): # O(1)
            if noLongRepeats(primer):
                if noHairpins(primer):
                    if noHomodimers(primer):
                        return primer


def __evaluateAllKmers(kmers:dict[str,dict[int,list[tuple[str,str]]]], params:Parameters) -> list[Primer]:
    """evaluates kmers at each position for their suitability as primers

    Args:
        kmers (dict[str,dict[int,list[tuple[str,str]]]]): the dictionary produced by __reorganizeDataByPosition
        params (Parameters): a Parameters object

    Returns:
        list[Primer]: a list of suitable primers as Primer objects
    """
    # generator function for getting arguments
    def generateArgs() -> Iterator[tuple[str,int,list[tuple[str,str]],float,float,float,float]]:
        """ generates arguments for __evaluateKmersAtOnePosition
        """
        # create an Aho-Corasick Automaton for detecting homopolymers in primer sequences
        repeats = Automaton()
        for ch in 'ATCG':
            repeats.add_word(ch * params.maxRepeatLen, ch)
        repeats.make_automaton()
        
        # each contig needs to be evalutated
        for contig in kmers.keys():
            # each start position within the contig needs to be evaluated
            for start in kmers[contig].keys():
                # save arguments to pass in parallel
                yield (contig,
                       start,
                       kmers[contig][start],
                       params.minGc,
                       params.maxGc,
                       params.minTm,
                       params.maxTm,
                       params.p3_mvConc,
                       params.p3_dvConc,
                       params.p3_dntpConc,
                       params.p3_dnaConc,
                       params.p3_tempC,
                       params.p3_maxLoop,
                       params.tempTolerance,
                       repeats)

    # parallelize primer evaluations
    pool = multiprocessing.Pool(processes=params.numThreads)
    results = pool.starmap(__evaluateKmersAtOnePosition, generateArgs())
    pool.close()
    pool.join()

    # remove failed searches before returning
    return [x for x in results if x is not None]


def __buildOutput(kmers:dict[str,dict[str,tuple[str,int,str]]], candidates:list[Primer]) -> dict[str,dict[str,list[Primer]]]:
    """builds the candidate primer output

    Args:
        kmers (dict[str,dict[str,tuple[str,int,str]]]): the dictionary produced by __getSharedKmers
        candidates (list[Primer]): the list produced by __evaluateAllKmers

    Returns:
        dict[str,dict[str,list[Primer]]]: key=genome name; val=dict: key=contig; val=list of Primers
    """
    # initialize output
    out = dict()
    
    # for each candidate primer
    for cand in candidates:
        # identify which sequence is present in the genome
        try:
            entry = kmers[cand.seq]
            hairpinTm = cand.hairpinTm
            rcHairpin = cand.rcHairpin
            homodimerTm = cand.homodimerTm
            rcHomodimer = cand.rcHomodimer
        
        except:
            entry = kmers[cand.seq.reverse_complement()]
            hairpinTm = cand.rcHairpin
            rcHairpin = cand.hairpinTm
            homodimerTm = cand.rcHomodimer
            rcHomodimer = cand.homodimerTm
        
        # for each genome
        for name in entry.keys():
            # extract the data from the entry
            contig, start, strand = entry[name]

            # create the new primer for this genome
            primer = Primer(cand.seq, contig, start, len(cand.seq), strand)
            primer.hairpinTm = hairpinTm
            primer.rcHairpin = rcHairpin
            primer.homodimerTm = homodimerTm
            primer.rcHomodimer = rcHomodimer
            
            # save the Primer in this contig's list
            out[name] = out.get(name, dict())
            out[name][contig] = out[name].get(contig, list())
            out[name][contig].append(primer)
    
    return out


def __getCandidatesForOneGenome(name:str, kmers:dict[str,dict[str,tuple[str,int,str]]], seen:set[tuple[str]], params:Parameters) -> dict[str,dict[str,list[Primer]]]:
    """gets candidate primers for a single genome

    Args:
        name (str): the name of the genome to get primers for
        kmers (dict[str,dict[str,tuple[str,int,str]]]): the dictionary produced by __getSharedKmers
        seen (set[tuple[str]]): a set of kmer groups (tuple of strings) that have already been processed
        params (Parameters): a Parameters object

    Returns:
        dict[str.dict[str,list[Primer]]]: key=genome name; val=dict: key=contig; val=list of Primers
    """
    # reorganize data by each unique start positions for one genome
    positions = __reorganizeDataByPosition(name, kmers)
    
    # remove groups of kmers that have already been seen
    __removeRedundantKmerGroups(positions, seen)
    
    # get a list of the kmers that pass the evaulation
    candidates = __evaluateAllKmers(positions, params)
    
    # create a dictionary whose keys are contigs and values are lists of candidate primers
    return __buildOutput(kmers, candidates)


def _getAllCandidateKmers(params:Parameters, sharedExists:bool) -> dict[str,dict[str,list[Primer]]]:
    """gets all the candidate kmer sequences for a given ingroup

    Args:
        params (Parameters): a Parameters object
        sharedExists (bool): indicates if the shared kmers are already pickled

    Raises:
        RuntimeError: no shared ingroup kmers
        RuntimeError: shared ingroup kmers are not suitable for use as primers

    Returns:
        dict[str,dict[str,list[Primer]]]: key=genome name; val=dict: key=contig; val=list of Primers
    """
    # messages
    GAP = " "*4
    MSG_1 = f"{GAP}getting shared ingroup kmers that appear once in each genome"
    MSG_2A = f"{GAP}evaluating "
    MSG_2B = " kmers"
    MSG_3A = f"{GAP}identified "
    MSG_3B = " candidate kmers"
    ERR_MSG_1 = "failed to identify a set of kmers shared between the ingroup genomes"
    ERR_MSG_2 = "none of the ingroup kmers are suitable for use as a primer"
    
    # initialize variables
    clock = Clock()
    numCand = 0
    seen = set()
    out = dict()
    
    # setup debugger
    params.log.rename(_getAllCandidateKmers.__name__)
    
    if sharedExists:
        # load existing shared kmers from file
        kmers = params.loadObj(params.pickles[Parameters._SHARED])

    else:
        # get all non-duplicated kmers that are shared in the ingroup
        params.log.info(MSG_1)
        clock.printStart(MSG_1)
        kmers = __getSharedKmers(params)
        clock.printDone()
        
        # move log back to this function
        params.log.rename(_getAllCandidateKmers.__name__)
        params.log.info(f'{GAP}done {clock.getTimeString()}')
    
        # make sure that ingroup kmers were identified
        if kmers == dict():
            params.log.error(ERR_MSG_1)
            raise RuntimeError(ERR_MSG_1)
        
        # dump the shared kmers to file
        params.dumpObj(kmers, params.pickles[Parameters._SHARED], "shared kmers", prefix=GAP)
    
    # print status
    clock.printStart(f'{MSG_2A}{len(kmers)}{MSG_2B}')
    params.log.info(f'{MSG_2A}{len(kmers)}{MSG_2B}')
    
    # go through each genome name
    names = list(next(iter(kmers.values())).keys())
    for name in names:
        # get the candidate kmers for the genome
        candidates = __getCandidatesForOneGenome(name, kmers, seen, params)
        
        # for each genome in the candidates
        for genome in candidates.keys():
            # create a sub dictionary if one does not already exist
            out[genome] = out.get(genome, dict())
            
            # for each contig in the genome
            for contig in candidates[genome].keys():
                # store a set of all the candidates for this contig
                out[genome][contig] = out[genome].get(contig, set())
                out[genome][contig].update(candidates[genome][contig])
    
    # done with seen; remove it
    del seen
    
    # for each genome
    for name in out.keys():
        # for each contig
        for contig in out[name].keys():
            # sort kmers by their start position on the (+) strand
            out[name][contig] = sorted(out[name][contig], key=lambda x: min(x.start, x.end))
            
            # count the number of candidate kmers for the first genome only
            if name == names[0]:
                numCand += len(out[name][contig])

    # make sure candidates were found
    if numCand == 0:
        params.log.error(ERR_MSG_2)
        raise RuntimeError(ERR_MSG_2)
    
    # print status
    clock.printDone()
    print(f"{MSG_3A}{numCand}{MSG_3B}")
    
    # log status
    params.log.info(f'{GAP}done {clock.getTimeString()}')
    params.log.info(f'{MSG_3A}{numCand}{MSG_3B}')
    
    # dump the candidate kmers to file
    params.dumpObj(out, params.pickles[Parameters._CAND], "candidate kmers", prefix=GAP)

    return out
