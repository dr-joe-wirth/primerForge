from Bio.Seq import Seq
import multiprocessing, os
from bin.Primer import Primer
from bin.Parameters import Parameters
from bin.Clock import Clock, _printDone, _printStart


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


def __binCandidatePrimers(candidates:dict[str,list[Primer]], minPrimerLen:int) -> dict[str,dict[int,list[Primer]]]:
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


def __evaluateOneBinPair(bin1:list[Primer], bin2:list[Primer], maxTmDiff:float, minProdLen:int, maxProdLen:int) -> tuple[Primer, Primer, int]:
    """evaluates a pair of bins of primers to find a single suitable pair; designed to run in parallel

    Args:
        bin1 (list[Primer]): the first (upstream) bin of primers
        bin2 (list[Primer]): the second (downstream) bin of primers
        maxTmDiff (float): the maximum difference in melting temps
        minProdLen (int): the minimum PCR product length
        maxProdLen (int): the maximum PCR product length
    
    Returns:
        tuple[Primer,Primer,int]: a pair of primers and their corresponding pcr product length
    """
    # constants
    FWD = 'forward'
    REV = 'reverse'
    
    # helper functions for evaluating primers
    def isThreePrimeGc(primer:Primer, direction:str=FWD) -> bool:
        """checks if a primer has a G or C at its 3' end
        """
        # constant
        GC = {"G", "C"}
        
        # different check depending if forward or reverse primer
        if direction == FWD:
            return primer.seq[-1] in GC
        if direction == REV:
            return primer.seq[0] in GC
    
    def isTmDiffWithinRange(p1:Primer, p2:Primer) -> bool:
        """ checks if the difference between melting temps is within the specified range
        """
        return abs(p1.Tm - p2.Tm) <= maxTmDiff

    def noPrimerDimer(p1:Primer, p2:Primer) -> bool:
        """verifies that the primer pair will not form primer dimers
        """
        # constant
        MAX_PID = 0.9
        
        # create the aligner; cost-free end gaps; no internal gaps
        from Bio.Align import PairwiseAligner, Alignment
        aligner = PairwiseAligner(mode='global',
                                  end_gap_score=0,
                                  end_open_gap_score=0,
                                  internal_open_gap_score=-float('inf'),
                                  match_score=2,
                                  mismatch_score=-1)
        
        # go through the alignments and check for high percent identities
        aln:Alignment
        for aln in aligner.align(p1.seq, p2.seq):
            matches = aln.counts().identities
            pid1 = matches / len(p1)
            pid2 = matches / len(p2)
            
            if max(pid1,pid2) > MAX_PID:
                return False
        
        return True
    
    # for each primer in the first bin
    for p1 in bin1:
        # only proceed if three prime end is GC
        if isThreePrimeGc(p1):
            # for each primer in the second bin
            for p2 in bin2:
                # get the pcr product length for this pair
                productLen = p2.end - p1.start + 1
                
                # only proceed if the product length is within allowable limits
                if minProdLen <= productLen <= maxProdLen:
                    # only proceed if the reverse primer's end is GC
                    if isThreePrimeGc(p2, REV):
                        # only proceed if the Tm difference is allowable
                        if isTmDiffWithinRange(p1, p2):
                            # only proceed if primer dimers won't form
                            if noPrimerDimer(p1, p2):
                                # found a suitable pair; stop comparing other primers in these bins
                                return p1, p2.reverseComplement(), productLen


def __evaluateBinPairs(binned:dict[str,dict[int,list[Primer]]], minPrimerLen:int, minProdLen:int, maxProdLen:int, maxTmDiff:float, numThreads:int) -> list[tuple[Primer,Primer,int]]:
    """evaluate pairs of bins to identify candidate primer pairs for a single genome

    Args:
        binned (dict[str,dict[int,list[Primer]]]): the dictionary produced by __binCandidatePrimers
        minPrimerLen (int): the minimum primer length
        minProdLen (int): the minimum PCR product length
        maxProdLen (int): the maximum PCR product length
        maxTmDiff (float): the maximum Tm difference bewteen two primers
        numThreads (int): the number of threads available for parallel processing

    Returns:
        list[tuple[Primer,Primer,int]]: a list of tuples produced by __evaluateOneBinPair
    """
    # initialize a list of arguments for __evaluateOneBinPair
    args = list()

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
                    args.append((bin1, bin2, maxTmDiff, minProdLen, maxProdLen))
    
    # process pairs of bins in parallel
    with multiprocessing.Pool(numThreads) as pool:
        out = pool.starmap(__evaluateOneBinPair, args)
    
    # remove failed bin pairs (None) before returning
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


def __getAllSharedPrimerPairs(firstName:str, candidateKmers:dict[str,dict[str,list[Primer]]], candidatePairs:list[tuple[Primer,Primer,int]], minProdLen:int, maxProdLen:int) -> dict[tuple[Primer,Primer],dict[str,tuple[str,int]]]:
    """gets all the primer pairs that are shared in all the genomes

    Args:
        firstName (str): the name of the genome that has already been evaluated
        candidateKmers (dict[str,dict[str,list[Primer]]]): key=genome name; val=dict: key=contig name; val=list of candidate Primers
        candidatePairs (list[tuple[Primer,Primer,int]]): the list produced by __evaluateBinPairs
        minProdLen (int): the minimum PCR product length
        maxProdLen (int): the maximum PCR product length

    Returns:
        dict[tuple[Primer,Primer],dict[str,tuple[str,int]]]: key=pair of Primers; val=dict: key=genome name; val=tuple: contig name, PCR product length
    """
    # initialize variables    
    out = dict()
    allowedLengths = range(minProdLen, maxProdLen+1)
    
    # get a set of genomes that need to be evaluated (the first name has already been evaluated)
    remaining = set(candidateKmers.keys())
    remaining.remove(firstName)

    # for each remaining genome
    kmers = dict()
    for name in remaining: 
        # restructure the data for O(1) lookup of primers by their sequences
        kmers[name] = __restructureCandidateKmerData(candidateKmers[name])
        
    # for each pair of primers in the candidate pairs
    k1:Primer
    k2:Primer
    for p1,p2,length in candidatePairs:
        # store the data for the genome (firstName) that has already been evaluated
        out[(p1,p2)] = dict()
        out[(p1,p2)][firstName] = (p1.contig,length)
        
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


def _getPrimerPairs(candidateKmers:dict[str,dict[str,list[Primer]]], params:Parameters) -> dict[tuple[Primer,Primer],dict[str,int]]:
    """gets pairs of primers suitable for use in all ingroup genomes

    Args:
        candidateKmers (dict[str,dict[str,list[Primer]]]): key=genome name; val=dict: key=contig name; val=list of candidate primers
        minPrimerLen (int): minimum primer length
        minProdLen (int): minimum PCR product length
        maxProdLen (int): maximum PCR product length
        maxTmDiff (float): maximum Tm difference between a pair of primers
        numThreads (int): the number of threads available for multiprocessing

    Raises:
        RuntimeError: unable to find any suitable primer pairs 
        RuntimeError: unable to find suitable primer pairs shared in all ingroup genomes

    Returns:
        dict[tuple[Primer,Primer],dict[str,int]]: the dictionary produced by __getAllSharedPrimerPairs
    """
    # messages
    ERR_MSG_1 = "could not identify suitable primer pairs from the candidate kmers"
    ERR_MSG_2 = "could not identify primer pairs present in every ingroup genome"
    
    # bin kmers for `firstName` only
    firstName = next(iter(candidateKmers.keys()))
    binnedCandidateKmers = __binCandidatePrimers(candidateKmers[firstName], params.minLen)
    
    # evaluate the pairs of bins in parallel
    candidatePairs = __evaluateBinPairs(binnedCandidateKmers, params.minLen, params.minPcr, params.maxPcr, params.maxTmDiff, params.numThreads)
    
    # make sure that some candidate pairs were identified
    if candidatePairs == []:
        if params.debug:
            params.debugger.setLogger(_getPrimerPairs.__name__)
            params.debugger.writeErrorMsg(ERR_MSG_1)
        raise RuntimeError(ERR_MSG_1)
    
    # find all the pairs that are shared between every ingroup genome
    pairs = __getAllSharedPrimerPairs(firstName, candidateKmers, candidatePairs, params.minPcr, params.maxPcr)
    
    # make sure that some candidate pairs still exist
    if len(pairs) == 0:
        if params.debug:
            params.debugger.setLogger(_getPrimerPairs.__name__)
            params.debugger.writeErrorMsg(ERR_MSG_1)
        raise RuntimeError(ERR_MSG_2)
    
    return pairs
