import os
from Bio import SeqIO
from Bio.Seq import Seq
from bin.Clock import Clock
from khmer import Countgraph
from bin.Primer import Primer
import multiprocessing, primer3
from ahocorasick import Automaton
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters

# global constant
__JUNCTION_CHAR = "~"

# functions
def __getAllowedKmers(fwd:str, rev:str, k:int, first:bool, name:str) -> tuple[str,set[str]]:
    """gets all the kmers that are allowed

    Args:
        fwd (str): the entire (concatenated) forward sequence
        rev (str): the entire (concatenated) reverse sequence
        k (int): the length of the kmers
        first (bool): indicates if this is the first genome to be processed
        name (str): the name of the input genome

    Returns:
        set[str]: a set of the allowed kmers
    """
    # helper functions
    def isOneEndGc(seq:str) -> bool:
        """ evaluates if one end of a sequence is a GC
        """
        GC = {"G", "C", 'g', 'c'}
        return seq[-1] in GC or seq[0] in GC
    
    def appearsOnce(seq:str, kh:Countgraph) -> bool:
        return kh.get(seq) == 1
    
    def atJunction(seq:str) -> bool:
        return __JUNCTION_CHAR in seq
    
    # create count graphs for both strands
    fkh = Countgraph(k, 1e7, 1)
    rkh = Countgraph(k, 1e7, 1)
    
    # consume the sequences so kmers can be counted
    fkh.consume(fwd)
    rkh.consume(rev)
    
    # extract the allowed kmers for both strands
    out    = {x for x in fkh.get_kmers(fwd) if isOneEndGc(x) and appearsOnce(x, fkh) and not atJunction(x)}
    rKmers = {x for x in rkh.get_kmers(rev) if isOneEndGc(x) and appearsOnce(x, rkh) and not atJunction(x)}
    
    # remove the reverse kmers if this is the first genome
    if first:
        out.difference_update(rKmers)
    
    # otherwise remove kmers that are shared between strands
    else:
        out.symmetric_difference_update(rKmers)
    
    return name,out


def __getAllAllowedKmers(params:Parameters) -> Automaton:
    """gets all the allowed kmers for the entire ingroup

    Args:
        params (Parameters): a Parameters object

    Returns:
        Automaton: an Aho-Corasick automaton for substring searching
    """
    # message
    MSG = " shared kmers after processing "
    
    # initialize variables
    tmp = dict()
    rec:SeqRecord
    args = list()
    out = Automaton()
    firstGenome = True
    
    # for each file
    for fn in params.ingroupFns:
        # create lists of the contig sequences as strings
        fwd = list()
        rev = list()
        for rec in SeqIO.parse(fn, params.format):
            fwd.append(str(rec.seq))
            rev.append(str(rec.seq.reverse_complement()))
        
        # concatenate contig sequences separated by junction characters
        # this will prevent the creation of "chimeric junction kmers"
        fwd = (__JUNCTION_CHAR * params.maxLen).join(fwd)
        rev = (__JUNCTION_CHAR * params.maxLen).join(rev)
        
        # for each kmer length
        for k in range(params.minLen, params.maxLen+1):
            # create a list of arguments for __getAllowedKmers
            args.append((fwd, rev, k, firstGenome, os.path.basename(fn)))
        
        # reset first genome boolean
        firstGenome = False
    
    # identify the allowed kmers for each k and each genome in parallel
    pool = multiprocessing.Pool(params.numThreads)
    allowed = pool.starmap(__getAllowedKmers, args)
    pool.close()
    pool.join()
    
    # sort the list by genome name
    allowed.sort(key=lambda x: x[0], reverse=True)
    
    # save the first name that will be examined
    prevName = allowed[-1][0]
    
    # save the shared kmers for each k
    while allowed != []:
        # remove the last item from the list
        name,a = allowed.pop()
        
        # log the number of shared kmers after processing a genome
        if name != prevName:
            params.log.debug(f'{sum(map(len, tmp.values()))}{MSG}{prevName}')
            prevName = name
        
        # save the intersection for this k if it already exists
        try:
            tmp[len(next(iter(a)))].intersection_update(a)
        
        # otherwise save the entire set for this k
        except KeyError:
            tmp[len(next(iter(a)))] = a
    
    # log the number of shared kmers after processing the last genome
    params.log.debug(f'{sum(map(len, tmp.values()))}{MSG}{prevName}')
    
    # collapse the dictionary into an Automaton for substring searching
    for k,kmers in tmp.items():
        # use pop to prevent data duplication
        while kmers != set():
            kmer = kmers.pop()
            out.add_word(kmer, kmer)
    
    # finish automaton creation
    out.make_automaton()
    
    return out


def __extractKmerData(name:str, contig:str, strand:str, seq:str, allowed:Automaton, shared:dict[str,dict[str,tuple[str,int,str]]]) -> None:
    """extracts positional data from a contig for the allowed kmers. modifies the input; does not return.

    Args:
        name (str): the name of the sequence file
        contig (str): the name of the contig
        strand (str): the strand of the sequence
        seq (str): the dna sequence
        allowed (Automaton): an Aho-Corasick automaton for substring searching
        shared (dict[str,dict[str,tuple[str,int,str]]]): the dictionary where results will be stored
    """
    for end,kmer in allowed.iter(seq):        
        # extract the start position if on the (+) strand
        if strand == Primer.PLUS:
            start = end - len(kmer) + 1
        
        # extract the start position if on the (-) strand
        else:
            start = len(seq) - end - 1
        
        # extract the positional data
        entry = {name: (contig, start, strand)}
        
        # save the data in the shared object
        try:
            shared[kmer].update(entry)
        except KeyError:
            shared[kmer] = entry


def __getSharedKmers(params:Parameters) -> dict[str,dict[str,tuple[str,int,str]]]:
    """retrieves all the kmers that are shared between the input genomes

    Args:
        params (Parameters): a Parameters object

    Returns:
        dict[str,dict[str,tuple[str,int,str]]]: key=kmer; val=dict: key=genome name: val=tuple: contig, start, strand
    """
    # intialize variables
    contig:SeqRecord
    out = dict()
    
    # get all of the kmers that are allowed for the ingroup genomes
    allAllowed = __getAllAllowedKmers(params)
    
    # for each file
    for fn in params.ingroupFns:
        # extract the name of the file
        name = os.path.basename(fn)
        
        # for each contig
        for contig in SeqIO.parse(fn, params.format):
            # extract data for the (+) strand
            __extractKmerData(name, contig.id, Primer.PLUS, str(contig.seq), allAllowed, out)
            
            # extract data for the (-) strand
            __extractKmerData(name, contig.id, Primer.MINUS, str(contig.seq.reverse_complement()), allAllowed, out)
    
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


def __evaluateKmersAtOnePosition(contig:str, start:int, posL:list[tuple[str,str]], minGc:float, maxGc:float, minTm:float, maxTm:float) -> Primer:
    """evaluates all the primers at a single position in the genome; designed for parallel calls

    Args:
        contig (str): the name of the contig
        start (int): the start position in the sequence
        posL (list[tuple[str,int]]): a list of tuples: kmer, strand
        minGc (float): the minimum percent GC allowed
        maxGc (float): the maximum percent GC allowed
        minTm (float): the minimum melting temperature allowed
        maxTm (float): the maximum melting temperature allowed
    
    Returns:
        Primer: a suitable primer at the given position
    """
    # constant
    FIVE_DEGREES = 5
    
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
        # constants
        MAX_LEN = 4
        REPEATS = ("A"*MAX_LEN, "T"*MAX_LEN, "C"*MAX_LEN, "G"*MAX_LEN)

        # check for each repeat in the primer
        for repeat in REPEATS:
            if repeat in primer.seq:
                return False
        return True

    def noHairpins(primer:Primer) -> bool:
        """verifies that the primer does not form hairpins
        """
        # hairpin tm should be less than (minTm - 5°); need to check both strands
        fwdOk = primer3.calc_hairpin_tm(str(primer)) < (minTm - FIVE_DEGREES)
        revOk = primer3.calc_hairpin_tm(str(primer.reverseComplement())) < (minTm - FIVE_DEGREES)
        
        return fwdOk and revOk

    def noHomodimers(primer:Primer) -> bool:
        """verifies that the primer does not form homodimers
        """
        # homodimer tm should be less than (minTm - 5°); need to check both strands
        fwdOk = primer3.calc_homodimer_tm(str(primer)) < (minTm - FIVE_DEGREES)
        revOk = primer3.calc_homodimer_tm(str(primer.reverseComplement())) < (minTm - FIVE_DEGREES)
        
        return fwdOk and revOk
    
    # initialize values for the while loop
    idx = 0
    
    # continue to iterate through each primer in the list until a primer is found 
    for idx in range(len(posL)):
        # extract data from the list
        seq,strand = posL[idx]
        
        # create a Primer object
        primer = Primer(seq, contig, start, len(seq), strand)
        
        # evaluate the primer's percent GC, Tm, hairpin potential, and homodimer potential; save if passes
        if isGcWithinRange(primer) and isTmWithinRange(primer): # O(1)
            if noLongRepeats(primer): # O(1)
                if noHairpins(primer):
                    if noHomodimers(primer):
                        return primer


def __evaluateAllKmers(kmers:dict[str,dict[int,list[tuple[str,str]]]], minGc:float, maxGc:float, minTm:float, maxTm:float, numThreads:int) -> list[Primer]:
    """evaluates kmers at each position for their suitability as primers

    Args:
        kmers (dict[str,dict[int,list[tuple[str,str]]]]): the dictionary produced by __reorganizeDataByPosition
        minGc (float): the minimum percent G+C
        maxGc (float): the maximum percent G+C
        minTm (float): the minimum melting temperature
        maxTm (float): the maximum melting temperature
        numThreads (int): number of threads for parallel processing

    Returns:
        list[Primer]: a list of suitable primers as Primer objects
    """
    # initialize a list of arguments
    args = list()

    # each contig needs to be evalutated
    for contig in kmers.keys():
        # each start position within the contig needs to be evaluated
        for start in kmers[contig].keys():
            # save arguments to pass in parallel
            args.append((contig, start, kmers[contig][start], minGc, maxGc, minTm, maxTm))

    # parallelize primer evaluations
    pool = multiprocessing.Pool(processes=numThreads)
    results = pool.starmap(__evaluateKmersAtOnePosition, args)
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
        except:
            entry = kmers[cand.seq.reverse_complement()]
        
        # for each genome
        for name in entry.keys():
            # extract the data from the entry
            contig, start, strand = entry[name]
        
            # save a Primer object in this contig's list
            out[name] = out.get(name, dict())
            out[name][contig] = out[name].get(contig, list())
            out[name][contig].append(Primer(cand.seq, contig, start, len(cand.seq), strand))
    
    return out


def __getCandidatesForOneGenome(name:str, kmers:dict[str,dict[str,tuple[str,int,str]]], params:Parameters) -> dict[str,dict[str,list[Primer]]]:
    """gets candidate primers for a single genome

    Args:
        name (str): the name of the genome to get primers for
        kmers (dict[Seq,dict[str,tuple[str,int,int]]]): the dictionary produced by __getSharedKmers
        params (Parameters): a Parameters object

    Returns:
        dict[str.dict[str,list[Primer]]]: key=genome name; val=dict: key=contig; val=list of Primers
    """
    # reorganize data by each unique start positions for one genome
    positions = __reorganizeDataByPosition(name, kmers)
    
    # get a list of the kmers that pass the evaulation
    candidates = __evaluateAllKmers(positions, params.minGc, params.maxGc, params.minTm, params.maxTm, params.numThreads)
    
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
    # constants
    SHARED_KMER_NUM = 0
    CANDID_KMER_NUM = 1
    
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
    out = dict()
    
    # setup debugger
    params.log.rename(_getAllCandidateKmers.__name__)
    
    if sharedExists:
        # load existing shared kmers from file
        kmers = params.loadObj(params.pickles[SHARED_KMER_NUM])

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
        params.dumpObj(kmers, params.pickles[SHARED_KMER_NUM], "shared kmers", prefix=GAP)
    
    # print status
    clock.printStart(f'{MSG_2A}{len(kmers)}{MSG_2B}')
    params.log.info(f'{MSG_2A}{len(kmers)}{MSG_2B}')
    
    # go through each genome name
    names = list(next(iter(kmers.values())).keys())
    for name in names:
        # get the candidate kmers for the genome
        candidates = __getCandidatesForOneGenome(name, kmers, params)
        
        # for each genome in the candidates
        for genome in candidates.keys():
            # create a sub dictionary if one does not already exist
            out[genome] = out.get(genome, dict())
            
            # for each contig in the genome
            for contig in candidates[genome].keys():
                # store a set of all the candidates for this contig
                out[genome][contig] = out[genome].get(contig, set())
                out[genome][contig].update(candidates[genome][contig])
    
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
    params.dumpObj(out, params.pickles[CANDID_KMER_NUM], "candidate kmers", prefix=GAP)

    return out
