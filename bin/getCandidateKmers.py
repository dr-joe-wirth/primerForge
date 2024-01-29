import multiprocessing
from Bio.Seq import Seq
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters
from bin.Clock import Clock, _printStart, _printDone

# global constants
__PLUS  = "+"
__MINUS = "-"

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


def __getUniqueKmers(seqs:list[SeqRecord], minLen:int, maxLen:int, name:str) -> dict[str,dict[Seq,dict[str,tuple[str,int,int]]]]:
    """gets a collection of kmers from a genome that are:
        * not repeated anywhere in the genome
        * at least one end is a G or C

    Args:
        seqs (list[SeqRecord]): a list of contigs as SeqRecord objects
        minLen (int): the minimum kmer length
        maxLen (int): the maximum kmer length
    
    Returns:
        dict[str,dict[Seq,dict[str,tuple[str,int,int]]]]: key=strand; val=dict: key=kmer seq; val=dict: key=name; val=tuple: contig, start, klen
    """
    # constant
    GC = {"G", "C", 'g', 'c'}
    
    # helper functions to clarify boolean expressions
    def isOneEndGc(seq:Seq) -> bool:
        """ evaluates if one end of a sequence is a GC
        """
        return seq[-1] in GC or seq[0] in GC

    def isDuplicated(seq:Seq) -> bool:
        return seq in kmers[__PLUS].keys() or seq in kmers[__MINUS].keys()

    # initialize variables
    kmers = dict()
    kmers[__PLUS] = dict()
    kmers[__MINUS] = dict()
    bad = set()
    krange = range(minLen, maxLen + 1)
    smallest = min(krange)
    
    # go through each contig
    for contig in seqs:
        # extract the contig sequences
        fwdSeq:Seq = contig.seq
        revSeq:Seq = fwdSeq.reverse_complement()
        
        # get the length of the contig
        contigLen = len(contig)
        done = False
        
        # get every possible kmer start position
        for start in range(contigLen):
            # for each allowed kmer length
            for klen in krange:
                # stop looping through the contig once we're past the smallest kmer
                if start+smallest > contigLen:
                    done = True
                    break
                
                # stop looping through the kmers once the length is too long
                elif start+klen > contigLen:
                    break
            
                # proceed if the extracted kmer length is good
                else:
                    # extract the kmer sequences
                    fwdKmer = fwdSeq[start:start+klen]
                    revKmer = revSeq[-(start+klen):-start]
                    
                    # mark duplicate kmers for removal
                    if isDuplicated(fwdKmer):
                        bad.add(fwdKmer)
                        bad.add(revKmer)
                    
                    # only save kmers that have GC at one end
                    elif isOneEndGc(fwdSeq):
                        kmers[__PLUS][fwdKmer]  = {name: (contig.id, start, klen)}
                        kmers[__MINUS][revKmer] = {name: (contig.id, start, klen)}
            
            # stop iterating through the contig when we're done with it
            if done:
                break

    # discard any sequences that were duplicated
    for seq in bad:
        try: kmers[__PLUS].pop(seq)
        except KeyError: pass
        
        try: kmers[__MINUS].pop(seq)
        except KeyError: pass
    
    return kmers
            

def __getSharedKmers(seqs:dict[str,list[SeqRecord]], params:Parameters) -> dict[Seq,dict[str,tuple[str,int,int]]]:
    """retrieves all the kmers that are shared between the input genomes

    Args:
        seqs (dict[str, list[SeqRecord]]): key=genome name; val=list of contigs as SeqRecord objects
        params (Parameters): a Parameters object

    Returns:
        dict[Seq,dict[str,tuple[str,int,int]]]: key=kmer sequence; val=dict: key=genome name; val=(contig, start, length)
    """
    # messages
    ERR_MSG_1 = 'failed to extract kmers from '
    
    # intialize variables
    sharedKmers = dict()
    firstGenome = True

    # for each genome
    for name in seqs.keys():
        # get the unique kmers for this genome/kmer length pair
        kmers = __getUniqueKmers(seqs[name], params.minLen, params.maxLen, name)
        
        # make sure there are kmers for this genome
        if kmers == dict():
            if params.debug:
                params.log.setLogger(__getSharedKmers.__name__)
                params.log.writeErrorMsg(f"{ERR_MSG_1}{name}")
            raise RuntimeError(f"{ERR_MSG_1}{name}")
        
        # keep only the (+) strand kmers if this is the first genome
        if firstGenome:
            sharedKmers.update(kmers[__PLUS])
            firstGenome = False
        
        # otherwise keep the shared kmers (both + and -)
        else:
            # get all the kmers for this genome (both + and - strands)
            thisGenome = set(kmers[__PLUS].keys()).union(kmers[__MINUS].keys())
            
            # go through each kmer that is currently shared by all processed genomes
            shared = set(sharedKmers.keys())
            for kmer in shared:
                # save data for kmers shared with this genome
                if kmer in thisGenome:
                    # strand information is now lost; try both
                    try:
                        sharedKmers[kmer].update(kmers[__PLUS][kmer])
                    except KeyError:
                        sharedKmers[kmer].update(kmers[__MINUS][kmer])
                
                # remove any kmers that are not shared with this genome
                else:
                    sharedKmers.pop(kmer)
        
        # log the number of shared kmers if debugging
        if params.debug:
            params.log.writeDebugMsg(f"{' '*8}{len(sharedKmers)} shared kmers after processing {name}")
    
    return sharedKmers


def __reorganizeDataByPosition(name:str, kmers:dict[Seq,dict[str,tuple[str,int,int]]]) -> dict[str,dict[int,list[tuple[Seq,int]]]]:
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
        contig, start, length = kmers[kmer][name]
        
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
            contig, start, length = entry[name]
        
            # save a Primer object in this contig's list
            out[name] = out.get(name, dict())
            out[name][contig] = out[name].get(contig, list())
            out[name][contig].append(Primer(cand.seq, contig, start, length))
    
    return out


def __getCandidatesForOneGenome(name:str, kmers:dict[Seq,dict[str,tuple[str,int,int]]], params:Parameters) -> dict[str,dict[str,list[Primer]]]:
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


def _getAllCandidateKmers(ingroup:dict[str,list[SeqRecord]], params:Parameters) -> dict[str,dict[str,list[Primer]]]:
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
        debug (bool, optional): indicates if it should run in debug mode. Defaults to False.

    Raises:
        RuntimeError: no shared ingroup kmers
        RuntimeError: shared ingroup kmers are not suitable for use as primers

    Returns:
        dict[str,dict[str,list[Primer]]]: key=genome name; val=dict: key=contig; val=list of Primers
    """
    # constant for debugging
    KMER_FN = "sharedKmers.p"

    # messages
    GAP = " "*4
    MSG_1 = f"{GAP}getting shared ingroup kmers that appear once in each genome"
    MSG_2 = f"{GAP}evaluating kmers"
    MSG_3A = f"{GAP}identified "
    MSG_3B = " candidate kmers"
    ERR_MSG_1 = "failed to identify a set of kmers shared between the ingroup genomes"
    ERR_MSG_2 = "none of the ingroup kmers are suitable for use as a primer"
    
    # initialize variables
    clock = Clock()
    numCand = 0
    shared = set()
    allseen = set()
    out = dict()
    
    # setup debugger
    if params.debug:
        params.log.setLogger(_getAllCandidateKmers.__name__)
    
    # get all non-duplicated kmers that are shared in the ingroup
    _printStart(clock, MSG_1)
    if params.debug: params.log.writeInfoMsg(f'{MSG_1}')
    kmers = __getSharedKmers(ingroup, params)
    
    # make sure that ingroup kmers were identified
    if kmers == dict():
        if params.debug:
            params.log.writeErrorMsg(ERR_MSG_1)
        raise RuntimeError(ERR_MSG_1)
    
    # print status
    _printDone(clock)
    
    # save data if debugging
    if params.debug:
        params.log.writeInfoMsg(f'{GAP}done {clock.getTime()}')
        params.dumpObj(kmers, KMER_FN, "shared kmers")
    
    # print status
    _printStart(clock, MSG_2)
    if params.debug: params.log.writeInfoMsg(MSG_2)
    
    # go through each genome name
    names = list(next(iter(kmers.values())).keys())
    for name in names:
        # initialize a set to track the kmers seen in this genome
        seen = set()
        
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
                
                # keep track of the kmers seen in this genome
                seen.update(candidates[genome][contig])
        
        # keep track of all the kmers seen
        allseen.update(seen)
        
        # calculate the kmers shared in all genomes
        if shared == set():
            shared.update(seen)
        else:
            shared.intersection_update(seen)
    
    # get all unshared kmers
    bad = allseen.difference(shared)
    
    # for each genome
    for name in out.keys():
        # for each genome
        for contig in out[name].keys():
            # remove unshared contigs and sort kmers by start position
            out[name][contig].difference_update(bad)
            out[name][contig] = sorted(out[name][contig], key=lambda x: x.start)
            
            # count the number of candidate kmers for the first genome only
            if name == names[0]:
                numCand += len(out[name][contig])

    # make sure candidates were found
    if numCand == 0:
        if params.debug: params.log.writeErrorMsg(ERR_MSG_2)
        raise RuntimeError(ERR_MSG_2)
    
    # print status
    _printDone(clock)
    print(f"{MSG_3A}{numCand}{MSG_3B}")
    
    # save data if debugging
    if params.debug:
        params.log.writeInfoMsg(f'{GAP}done {clock.getTime()}')
        params.log.writeInfoMsg(f'{MSG_3A}{numCand}{MSG_3B}')

    return out
