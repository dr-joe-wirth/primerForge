import multiprocessing
from Bio.Seq import Seq
from bin.Clock import Clock, _printStart, _printDone
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from multiprocessing.managers import ListProxy

# global constants
__PLUS_STRAND = "+"
__MINUS_STRAND = "-"

# functions
def _kmpSearch(text:str, pattern:str) -> tuple[bool, int]:
    """implementation of the KMP string search algorithm: O(n+m)

    Args:
        text (str): text to search
        pattern (str): pattern to search for

    Returns:
        tuple[bool, int]: indicates if the pattern is found in the text; the start position of the match
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
            return True, (idx - patternLen)
    
    return False, -1


def __getUniqueKmers(seqs:list[SeqRecord], minLen:int, maxLen:int) -> dict[str,dict[Seq,tuple[str,int,int]]]:
    """gets a collection of kmers from a genome that are:
        * not repeated anywhere in the genome
        * at least one end is a G or C

    Args:
        seqs (list[SeqRecord]): a list of contigs as SeqRecord objects
        minLen (int): the minimum kmer length
        maxLen (int): the maximum kmer length
    
    Returns:
        dict[str,dict[Seq,tuple[str,int,int]]]: key=strand; val=dict: key=kmer seq; val=contig id, start position, kmer length
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
    kmers[__PLUS_STRAND] = dict()
    kmers[__MINUS_STRAND] = dict()
    bad = set()
    
    # go through each contig
    for contig in seqs:
        fwdSeq:Seq = contig.seq
        revSeq:Seq = fwdSeq.reverse_complement()
        
        # go through each kmer length
        for primerLen in range(minLen, maxLen+1):
            # get every possible kmer start position
            for start in range(len(contig)- (primerLen - 1)):
                # extract the kmer sequences
                fwdKmer = fwdSeq[start:start+primerLen]
                revKmer = revSeq[-(start+primerLen):-start]
                
                # mark duplicate kmers for removal
                if fwdKmer in kmers[__PLUS_STRAND].keys() or fwdKmer in kmers[__MINUS_STRAND].keys():
                    bad.add(fwdKmer)
                
                # only save kmers that have GC at one end
                if isOneEndGc(fwdSeq):
                    kmers[__PLUS_STRAND][fwdKmer] = (contig.id, start, primerLen)
                    kmers[__MINUS_STRAND][revKmer] = (contig.id, start, primerLen)

    # discard any sequences that were duplicated
    for seq in bad:
        try: kmers[__PLUS_STRAND].pop(seq)
        except KeyError: pass
        
        try: kmers[__MINUS_STRAND].pop(seq)
        except KeyError: pass
    
    return kmers
            

def __getSharedKmers(seqs:dict[str,list[SeqRecord]], minLen:int, maxLen:int) -> dict[str,dict[Seq,tuple[str,int,int]]]:
    """retrieves all the kmers that are shared between the input genomes

    Args:
        seqs (dict[str, list[SeqRecord]]): key=genome name; val=list of contigs as SeqRecord objects
        minLen (int): the minimum kmer length
        maxLen (int): the maximum kmer length

    Returns:
        dict[str,dict[Seq,tuple[str,int,int]]]: key=genome name; val=dict: key=kmer sequence; val=(contig, start, length)
    """
    # intialize variables
    sharedKmers = set()
    kmers = dict()
    
    # get the unique kmers for each genome
    for name in seqs.keys():
        kmers[name] = __getUniqueKmers(seqs[name], minLen, maxLen)
    
    # identify the shared kmers
    # for each genome
    for name in kmers.keys():
        # keep only the (+) strand kmers if this is the first genome
        if sharedKmers == set():
            sharedKmers.update(set(kmers[name][__PLUS_STRAND].keys()))
        
        # otherwise keep the shared kmers (both + and -)
        else:
            # get all the kmers for this genome (both + and - strands)
            thisGenome = set(kmers[name][__PLUS_STRAND].keys()).union(kmers[name][__MINUS_STRAND].keys())
            
            # update the shared set with the kmers found in this genome
            sharedKmers.intersection_update(thisGenome)
    
    # remove the unshared kmers
    # for each genome
    for name in kmers.keys():
        # handle each strand separately
        for strand in kmers[name].keys():
            # get the set of kmers in this genome/strand that are not shared
            bad = set(kmers[name][strand].keys()).difference(sharedKmers)
        
            # remove all the unshared kmers from this dictionary
            for kmer in bad:
                kmers[name][strand].pop(kmer)
    
    # restructure the data
    # for each genome
    for name in kmers.keys():
        # save (-) strand in temp variable
        tmp = kmers[name][__MINUS_STRAND]
        
        # save all the (+) strand kmers under the genome name
        kmers[name] = kmers[name][__PLUS_STRAND]
        
        # save the rev-comp of the (-) strand kmers
        for seq in set(tmp.keys()):
            kmers[name][seq.reverse_complement()] = tmp.pop(seq)
    
    return kmers


def __reorganizeDataByPosition(kmers:dict[str,dict[Seq,tuple[str,int,int]]]) -> dict[str,dict[int,list[tuple[Seq,int]]]]:
    """reorganizes data from __getSharedKmers by its genomic position

    Args:
        kmers (dict[str,dict[Seq,tuple[str,int,int]]]): key=contig; val=dict: key=kmer sequence; val=contig,start,length

    Returns:
        dict[str,dict[str,dict[int,list[tuple[Seq,int]]]]]: key=contig; val=dict: key=start position; val=list of tuples: (kmer sequence, kmer length)
    """
    # initialize output
    out = dict()
    
    # for each kmer
    for kmer in kmers.keys():
        # extract data from the dictionary
        contig, start, length = kmers[kmer]
        
        # contig = top level key; start position = second level key; val = list
        out[contig] = out.get(contig, dict())
        out[contig][start] = out[contig].get(start, list())
        
        # add the sequence and its length to the list
        out[contig][start].append((kmer, length))
    
    return out


def __evaluateKmersAtOnePosition(contig:str, start:int, posL:list[tuple[Seq,int]], minGc:float, maxGc:float, minTm:float, maxTm:float) -> Primer:
    """evaluates all the primers at a single position in the genome; designed for parallel calls

    Args:
        contig (str): the name of the contig
        start (int): the start position in the sequence
        posL (list[tuple[Seq,int]]): a list of tuples; primer sequence and primer length
        minGc (float): the minimum percent GC allowed
        maxGc (float): the maximum percent GC allowed
        minTm (float): the minimum melting temperature allowed
        maxTm (float): the maximum melting temperature allowed
    
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
        # constants
        MAX_LEN = 4
        REPEATS = ("A"*MAX_LEN, "T"*MAX_LEN, "C"*MAX_LEN, "G"*MAX_LEN)

        # check for each repeat in the primer
        for repeat in REPEATS:
            if _kmpSearch(primer.seq, repeat)[0]:
                return False
        return True

    def noIntraPrimerComplements(primer:Primer) -> bool:
        """verifies that the primer does not have hairpin potential
        """
        # constants
        MAX_LEN = 3
        LEN_TO_CHECK = MAX_LEN + 1
        
        # get the reverse complement of the primer
        revComp = primer.seq.reverse_complement()
        
        for idx in range(len(primer)-MAX_LEN):
            # see if the reverse complement exists downstream
            if _kmpSearch(primer.seq, revComp[idx:idx+LEN_TO_CHECK])[0]:
                return False
        
        return True
    
    # initialize values for the while loop
    idx = 0
    
    # continue to iterate through each primer in the list until a primer is found 
    for idx in range(len(posL)):
        # extract data from the list
        seq,length = posL[idx]
        
        # create a Primer object
        primer = Primer(seq, contig, start, length)
        
        # evaluate the primer's percent GC, Tm, and homology; save if found
        if isGcWithinRange(primer) and isTmWithinRange(primer): # O(1)
            if noLongRepeats(primer): # O(1)
                if noIntraPrimerComplements(primer): # this runtime is the worst O(len(primer)); evaluate last
                    return primer


def __evaluateAllKmers(kmers:dict[str,dict[int,list[tuple[Seq,int]]]], minGc:float, maxGc:float, minTm:float, maxTm:float, numThreads:int) -> list[Primer]:
    """evaluates kmers at each position for their suitability as primers

    Args:
        kmers (dict[str,dict[int,list[tuple[Seq,int]]]]): the dictionary produced by __reorganizeDataByPosition
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
 

def __buildOutput(kmers:dict[str,dict[Seq,tuple[str,int,int]]], candidates:list[Primer]) -> dict[str,dict[str,list[Primer]]]:
    """builds the candidate primer output

    Args:
        kmers (dict[str,dict[Seq,tuple[str,int,int]]]): the dictionary produced by __getSharedKmers
        candidates (list[Primer]): the list produced by __evaluateAllKmers

    Returns:
        dict[str,dict[str,list[Primer]]]: key=genome name; val=dict: key=contig; val=list of Primers
    """
    # initialize output
    out = dict()
    
    # for each genome
    for name in kmers:
        out[name] = dict()
        
        # for each candidate primer
        for cand in candidates:
            # try to retrieve the data for the (+) strand
            try:
                contig, start, length = kmers[name][cand.seq]
            
            # otherwise retrieve the data for the (-) strand
            except:
                contig, start, length = kmers[name][cand.seq.reverse_complement()]
            
            # save the Primer object in this contig's list
            out[name][contig] = out[name].get(contig, list())
            out[name][contig].append(Primer(cand.seq, contig, start, length))
    
    # sort each list of primers by start their position
    for name in out.keys():
        for contig in out[name].keys():
            out[name][contig] = sorted(out[name][contig], key=lambda x: x.start)
    
    return out


def _getAllCandidateKmers(ingroup:dict[str,list[SeqRecord]], minLen:int, maxLen:int, minGc:float, maxGc:float, minTm:float, maxTm:float, numThreads:int) -> dict[str,dict[str,list[Primer]]]:
    """gets all the candidate kmer sequences for a given ingroup with respect to a given outgroup

    Args:
        ingroup (dict[str,list[SeqRecord]]): the ingroup sequences: key=genome name; val=contigs as SeqRecords
        minLen (int): minimum primer length
        maxLen (int): maximum primer length
        minGc (float): minimum primer G+C percent
        maxGc (float): maximum primer G+C percent
        minTm (float): minimum primer melting temp
        maxTm (float): maximum primer melting temp
        numThreads (int): number threads available for parallel processing

    Returns:
        dict[str,dict[str,list[Primer]]]: key=genome name; val=dict: key=contig; val=list of Primers
    """
    GAP = " "*4
    MSG_1 = GAP + "getting shared ingroup kmers that appear once in each genome"
    MSG_2 = GAP + "evaluating candidate ingroup kmers"
    MSG_3A = GAP*2 + "identified "
    MSG_3B = " candidate kmers suitable for use as primers"
    ERR_MSG_1 = "failed to identify a set of kmers shared between the ingroup genomes"
    ERR_MSG_2 = "could not find ingroup kmers that are absent in the outgroup"
    ERR_MSG_3 = "none of the ingroup kmers are suitable for use as a primer"
    
    # start the timer
    clock = Clock()
    
    # get all non-duplicated kmers that are shared in the ingroup
    _printStart(clock, MSG_1)
    ingroupKmers = __getSharedKmers(ingroup, minLen, maxLen)
    _printDone(clock)
    
    # make sure that ingroup kmers were identified
    if len(ingroupKmers.values()) == 0:
        raise RuntimeError(ERR_MSG_1)
    
    # make sure that there are still kmers for the ingroup
    if len(ingroupKmers.values()) == 0:
        raise RuntimeError(ERR_MSG_2)
    
    # reorganize data by each unique start positions for one genome
    positions = __reorganizeDataByPosition(next(iter(ingroupKmers.values())))
    
    # get a list of the kmers that pass the evaulation
    _printStart(clock, MSG_2)
    candidates = __evaluateAllKmers(positions, minGc, maxGc, minTm, maxTm, numThreads)
    _printDone(clock)
    
    # make sure there are candidate kmers
    if candidates == []:
        raise RuntimeError(ERR_MSG_3)
    
    # print number of candidate kmers found
    print(f"{MSG_3A}{len(candidates)}{MSG_3B}")
    
    # create a dictionary whose keys are contigs and values are the primers
    return __buildOutput(ingroupKmers, candidates)
