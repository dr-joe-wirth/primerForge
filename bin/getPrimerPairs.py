import multiprocessing
from Bio.Seq import Seq
from bin.Primer import Primer
from bin.Parameters import Parameters


def __binOverlappingPrimers(candidates:dict[str,list[Primer]]) -> dict[str,dict[int,list[Primer]]]:
    """bins primers that overlap to form a contiguous sequence based on positional data

    Args:
        candidates (dict[str,list[Primer]]): key=contig; val=list of Primers (sorted by start position)

    Returns:
        dict[str,dict[int,list[Primer]]]: key=contig; val=dict: key=bin number; val=list of overlapping Primer objects
    """
    # initialize output
    out = dict()
    
    # for each contig in the dictionary
    for contig in candidates:
        # initialize variables
        currentBin = 0
        prevEnd = None
        out[contig] = dict()
        
        # for each candidate primer in the contig
        for cand in candidates[contig]:
            # the first candidate needs to go into its own bin
            if prevEnd is None:
                out[contig][currentBin] = [cand]
                prevEnd = cand.end
            
            # if the candidate overlaps the previous one, then add it to the bin
            elif cand.start < prevEnd:
                prevEnd = cand.end
                out[contig][currentBin].append(cand)
            
            # otherwise the candidate belongs in a new bin
            else:
                prevEnd = cand.end
                currentBin += 1
                out[contig][currentBin] = [cand]
    
    return out


def __minimizeOverlaps(bins:dict[str,dict[int,list[Primer]]], minimizerLen:int) -> None:
    """splits excessively long overlap bins by the minimizer sequences of the binned primers

    Args:
        bins (dict[str, dict[int,list[Primer]]]): the dictionary produced by __binOverlappingPrimers
        minimizerLen (int): the length of the minimizer sequence
    
    Returns:
        does not return. modifies input dictionary
    """
    # helper function to evaluate overlap length
    def overlapIsTooLong(primers:list[Primer]) -> bool:
        """determines if the overlap covered in the bin is too long"""
        MAX_LEN = 64
        length = primers[-1].end - primers[0].start
        return length > MAX_LEN
    
    # for each contig
    for contig in bins:
        # get the next empty bin number
        curBin = max(bins[contig].keys()) + 1
        
        # for each existing bin number
        for binNum in set(bins[contig].keys()):
            # if the length of the overlap is too long
            if overlapIsTooLong(bins[contig][binNum]):
                # get the primers from that bin and remove it from the dictionary
                primers = bins[contig].pop(binNum)
                
                # cluster each Primer by their minimizer sequence
                clusters = dict()
                for primer in primers:
                    minimizer = primer.getMinimizer(minimizerLen)
                    
                    clusters[minimizer] = clusters.get(minimizer, list())
                    clusters[minimizer].append(primer)
                
                # add each cluster as a new bin in the input dictionary
                for clust in clusters.values():
                    bins[contig][curBin] = clust
                    curBin += 1


def __binCandidateKmers(candidates:dict[str,list[Primer]], minPrimerLen:int) -> dict[str,dict[int,list[Primer]]]:
    """bins overlapping candidate primers; splits long overlaps on minimizer sequences

    Args:
        candidates (dict[str,list[Primer]]): key=contig; val=list of candidate Primers
        minPrimerSize (int): the minimum length of a primer

    Returns:
        dict[str,dict[int,list[Primer]]]: key=contig; val=dict: key=bin number; val=list of candidate Primers
    """
    # first bin overlapping primers
    bins = __binOverlappingPrimers(candidates)
    
    # next divide excessively long overlapping primers by minimizers (window size 1/2 min primer len)
    __minimizeOverlaps(bins, int(minPrimerLen / 2))
      
    return bins


def __getBinPairs(binned:dict[str,dict[int,list[Primer]]], minPrimerLen:int, minProdLen:int, maxProdLen:int) -> list[tuple[str,int,int]]:
    # initialize a list of arguments for __evaluateOneBinPair
    out = list()

    # for each contig in the genome
    for contig in binned.keys():
        # sort bins by their start position
        sortedBins = sorted(binned[contig].keys(), key=lambda x: binned[contig][x][0].start)
        
        # go through pairs of bins in order of their start positions
        for idx in range(len(sortedBins)-1):
            bin1 = binned[contig][sortedBins[idx]]
            for jdx in range(idx+1, len(sortedBins)):
                bin2 = binned[contig][sortedBins[jdx]]
                
                # calculate the smallest possible pcr product length
                smallest = (bin2[0].start + minPrimerLen) - (bin1[-1].end - minPrimerLen)
                
                # calculate the largest posisble pcr product length
                largest = bin2[-1].end - bin1[0].start
    
                # move to the next bin1 if the smallest product is too large
                if smallest > maxProdLen:
                    break
                
                # move to the next bin2 if largest possible product is too small
                elif largest < minProdLen:
                    continue
                
                # otherwise, these two bins could produce viable PCR product lengths; compare them
                else:
                    out.append((contig, sortedBins[idx], sortedBins[jdx]))
    return out


def __isPairSuitable(fwd:Primer, rev:Primer, minPcr:int, maxPcr:int, maxTmDiff:float) -> tuple[bool,int]:
    """determines if a pair of primers is suitable for PCR

    Args:
        fwd (Primer): the forward primer
        rev (Primer): the reverse primer
        minPcr (int): minimum PCR product size
        maxPcr (int): maximum PCR product size
        maxTmDiff (float): maximum Tm difference bewteen Primers

    Returns:
        tuple[bool,int]: bool indicating if the pair is suitable, int indicating the PCR product length
    """
    # calculate the pcr product length
    pcrLen = (rev.end - fwd.start + 1)

    # unsuitable if the product length is not within the allowed range
    if not minPcr <= pcrLen <= maxPcr: return False, pcrLen
    
    # unsuitable if the Tm difference is not within the allowed range
    if not abs(fwd.Tm - rev.Tm) <= maxTmDiff: return False, pcrLen
    
    return True, pcrLen


def __evaluateOnePair(fwd:Primer, rev:Primer, minPcr:int, maxPcr:int, maxTmDiff:float, binPair:tuple[int,int]) -> tuple[Primer,Primer,int,tuple[int,int]]:
    """evaluates a primer pair; designed to work in parallel; only returns if the pair passes evaluation

    Args:
        fwd (Primer): the forward Primer
        rev (Primer): the reverse Primer
        minPcr (int): the minimum PCR product size
        maxPcr (int): the maximum PCR product size
        maxTmDiff (float): the maximum difference between primer Tm
        binPair (tuple): 

    Returns:
        tuple[Primer,Primer,int,int,int]: forward primer, reverse primer, pcr product size, bin pair
    """
    # determine the product size and if the pair is suitable for PCR
    acceptablePair, pcrLen = __isPairSuitable(fwd, rev, minPcr, maxPcr, maxTmDiff)
    
    # return acceptable pairs and their pcr product sizes
    if acceptablePair: return fwd,rev,pcrLen,binPair


def __getCandidatePrimerPairs(binPairs:list[tuple[str,int,int]], bins:dict[str,dict[int,list[Primer]]], params:Parameters) -> list[tuple[Primer,Primer,int,tuple[int,int]]]:
    """gets candidate primer pairs from pairs of bins

    Args:
        binPairs (list[tuple[str,int,int]]): the list produced by __getBinPairs
        bins (dict[str,dict[int,list[Primer]]]): the list produced by __binCandidateKmers
        params (Parameters): a Parameters object

    Returns:
        list[tuple[Primer,Primer,int,int,int]]: a list of primer pairs and the corresponding pcr product size and the bin pair
    """
    # constants
    FWD = 'forward'
    REV = 'reverse'
    GC = {"G", "C"}
    
    # helper functions for evaluating primers
    def isThreePrimeGc(primer:Primer, direction:str=FWD) -> bool:
        """checks if a primer has a G or C at its 3' end
        """
        # different check depending if forward or reverse primer
        if direction == FWD:
            return primer.seq[-1] in GC
        if direction == REV:
            return primer.seq[0] in GC

    # initialize a list of arguments for parallel processing
    args = list()
    
    # for each pair of bins
    for contig,num1,num2 in binPairs:
        # extract the lists of primers for these bins
        bin1 = bins[contig][num1]
        bin2 = bins[contig][num2]

        # for each forward primer
        for fwd in bin1:
            # only evaluate if the 3' end is GC
            if isThreePrimeGc(fwd):
                # for each reverse primer
                for rev in bin2:
                    # only evaluate if the 3' end is GC
                    if isThreePrimeGc(rev, REV):
                        # save the arguments to the list to evaluate in parallel
                        args.append((fwd, rev, params.minPcr, params.maxPcr, params.maxTmDiff, (num1, num2)))
    
    # evaluate pairs of primers in parallel
    with multiprocessing.Pool(params.numThreads) as pool:
        out = pool.starmap(__evaluateOnePair, args)
    
    # remove failed results from the list before returning
    return [x for x in out if x is not None]


def __restructureCandidateKmerData(candidates:dict[str,list[Primer]]) -> dict[Seq,Primer]:
    """restructures candidate kmer data to allow O(1) Primer look-up by its sequence

    Args:
        candidates (dict[str,list[Primer]]): key=contig; val=list of Primer objects

    Returns:
        dict[Seq,Primer]: key=primer sequence; val=Primer object
    """
    # initialize the output
    out = dict()
    
    # for each primer in each contig
    for primers in candidates.values():
        for primer in primers:
            # store the primer under its sequence
            out[primer.seq] = primer
    
    return out


def __getAllSharedPrimerPairs(firstName:str, candidateKmers:dict[str,dict[str,list[Primer]]], candidatePairs:list[tuple[Primer,Primer,int]], minPcr:int, maxPcr:int) -> dict[tuple[Primer,Primer],dict[str,tuple[str,int]]]:
    """gets all the primer pairs that are shared in all the genomes

    Args:
        firstName (str): the name of the genome that has already been evaluated
        candidateKmers (dict[str,dict[str,list[Primer]]]): key=genome name; val=dict: key=contig name; val=list of candidate Primers
        candidatePairs (list[tuple[Primer,Primer,int]]): the list produced by __evaluateBinPairs
        minPcr (int): the minimum PCR product length
        maxPcr (int): the maximum PCR product length

    Returns:
        dict[tuple[Primer,Primer],dict[str,tuple[str,int]]]: key=pair of Primers; val=dict: key=genome name; val=tuple: contig name, PCR product length
    """
    # initialize variables
    k1:Primer
    k2:Primer
    out = dict()
    allowedLengths = range(minPcr, maxPcr+1)
    
    # get a set of genomes that need to be evaluated (the first name has already been evaluated)
    remaining = set(candidateKmers.keys())
    remaining.remove(firstName)

    # for each remaining genome
    kmers = dict()
    for name in remaining: 
        # restructure the data for O(1) lookup of primers by their sequences
        kmers[name] = __restructureCandidateKmerData(candidateKmers[name])
        
    # for each pair of primers in the candidate pairs
    for p1,p2,length in candidatePairs:
        # store the data for the genome (firstName) that has already been evaluated
        out[(p1,p2)] = dict()
        out[(p1,p2)][firstName] = (p1.contig,length)
        
        # for each unprocessed genome
        for name in remaining:
            # get the first kmer (may be rev comp)
            try:
                k1 = kmers[name][p1.seq]
                k1_rev = False
                
            except KeyError:
                k1 = kmers[name][p1.seq.reverse_complement()]
                k1_rev = True
            
            # get the second kmer (may be rev comp)
            try:
                k2 = kmers[name][p2.seq]
                k2_rev = False
    
            except KeyError:
                k2 = kmers[name][p2.seq.reverse_complement()]
                k2_rev = True
            
            # both primers need to be on the same contig
            if k1.contig == k2.contig:
                # only proceed if on opposite strands (p2 is already rev comped)
                if k1_rev == k2_rev:
                    # save values when k1 is the foward primer and k2 is the reverse primer
                    if k1.start < k2.start:
                        fwd = k1
                        rev = k2
                    
                    # save values when k1 is the reverse primer and k2 is the forward primer
                    else:
                        fwd = k2
                        rev = k1
                    
                    # only save the pair if the length falls within the acceptable range
                    length = rev.end - fwd.start + 1
                    if length in allowedLengths:
                        out[(p1,p2)][name] = (fwd.contig, length)

    # remove any pairs that are not universally suitable in all genomes
    for pair in set(out.keys()):
        if len(out[pair]) < len(candidateKmers.keys()):
            out.pop(pair)
    
    return out


def __keepOnePairPerBinPair(pairs:list[tuple[Primer,Primer,int,tuple[int,int]]]) -> None:
    """keeps only one primer pair per bin pair

    Args:
        pairs (list[tuple[Primer,Primer,int,tuple[int,int]]]): the list created by __getAllSharedPrimerPairs
    
    Returns:
        does not return; modifies input
    """
    # initialize a set to track the pairs that have already been seen
    seen = set()
    
    # go through the list indices in reverse order for on-the-fly popping
    for idx in range(len(pairs)-1,-1,-1):
        # extract the bin pair for this primer pair
        binPair = pairs[idx][-1]
        
        # remove any pairs that are already represented
        if binPair in seen:
            pairs.pop(idx)
        
        # otherwise this is the first time this bin pair has been seen
        else:
            seen.add(binPair)
            

def _getPrimerPairs(candidateKmers:dict[str,dict[str,list[Primer]]], params:Parameters) -> dict[tuple[Primer,Primer],dict[str,int]]:
    """gets primer pairs found in all the ingroup genomes

    Args:
        candidateKmers (dict[str,dict[str,list[Primer]]]): key=genome name; val=dict: key=contig; val=list of Primers
        params (Parameters): a Parameters object

    Raises:
        RuntimeError: unable to identify suitable primer pairs from the candidate kmers
        RuntimeError: unable to identify primer pairs in every ingroup genome

    Returns:
        dict[tuple[Primer,Primer],dict[str,int]]: key=Primer pairs; val=dict: key=genome name; val=tuple: contig, pcr product size
    """
    # messages
    ERR_MSG_1 = "could not identify suitable primer pairs from the candidate kmers"
    ERR_MSG_2 = "could not identify primer pairs present in every ingroup genome"
    
    # initialize variables
    allCand = list()
    out = dict()
    
    # process each genome as different results may be obtained
    for name in candidateKmers.keys():
        # bin kmers to reduce time complexity
        binnedCandidateKmers = __binCandidateKmers(candidateKmers[name], params.minLen)
        
        # get the bins that could work
        binPairs = __getBinPairs(binnedCandidateKmers, params.minLen, params.minPcr, params.maxPcr)

        # get candidate primer pairs
        candidatePairs = __getCandidatePrimerPairs(binPairs, binnedCandidateKmers, params)
        
        # get the primer pairs that are shared in all genomes
        pairs = __getAllSharedPrimerPairs(name, candidateKmers, candidatePairs, params.minPcr, params.maxPcr)
        
        # wittle down the pairs; keep only one pair per bin
        __keepOnePairPerBinPair(pairs)
        
        # keep track of all the candidate pairs and final pairs
        allCand.extend(candidatePairs)
        out.update(pairs)
    
    # make sure that some candidate pairs were identified
    if allCand == []:
        if params.debug:
            params.log.setLogger(_getPrimerPairs.__name__)
            params.log.writeErrorMsg(ERR_MSG_1)
        raise RuntimeError(ERR_MSG_1)
    
    # make sure that some candidate pairs were universal
    elif out == dict():
        if params.debug:
            params.log.setLogger(_getPrimerPairs.__name__)
            params.log.writeErrorMsg(ERR_MSG_2)
        raise RuntimeError(ERR_MSG_2)
    
    return out
