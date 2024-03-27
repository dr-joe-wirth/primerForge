from Bio.Seq import Seq
from bin.Primer import Primer
import multiprocessing, primer3
from bin.Parameters import Parameters


def __binCandidateKmers(candidates:dict[str,list[Primer]]) -> dict[str,dict[int,list[Primer]]]:
    """bins primers that overlap to form a contiguous sequence based on positional data

    Args:
        candidates (dict[str,list[Primer]]): key=contig; val=list of Primers (sorted by start position)

    Returns:
        dict[str,dict[int,list[Primer]]]: key=contig; val=dict: key=bin number; val=list of overlapping Primer objects
    """
    # constant
    MAX_BIN_LEN = 64
    
    # initialize output
    out = dict()
    
    # for each contig in the dictionary
    for contig in candidates.keys():
        # initialize variables
        currentBin = 0
        prevEnd = None
        binStart = 0
        out[contig] = dict()
        
        # for each candidate primer in the contig
        for cand in candidates[contig]:
            # get the current start/end on the (+) strand
            curStart = min(cand.start, cand.end)
            curEnd = max(cand.start, cand.end)
            
            # get the length of the current bin
            binLen = curEnd - binStart
            
            # the first candidate needs to go into its own bin
            if prevEnd is None:
                out[contig][currentBin] = [cand]
                prevEnd = curEnd
                binStart = curStart
            
            # if the candidate overlaps the previous one and the bin isn't too long, then add it to the bin
            elif curStart < prevEnd and binLen <= MAX_BIN_LEN:
                prevEnd = curEnd
                out[contig][currentBin].append(cand)
            
            # otherwise the candidate belongs in a new bin
            else:
                prevEnd = curEnd
                binStart = curStart
                currentBin += 1
                
                out[contig][currentBin] = [cand]
    
    return out


def __getBinPairs(binned:dict[str,dict[int,list[Primer]]], minPrimerLen:int, minProdLen:int, maxProdLen:int) -> list[tuple[str,int,int]]:
    """gets a list of bin pairs that can be used to find primer pairs

    Args:
        binned (dict[str,dict[int,list[Primer]]]): the dictionary produced by __binCandidateKmers
        minPrimerLen (int): the minimum primer length
        minProdLen (int): the minimum pcr product size
        maxProdLen (int): the maximum pcr product size

    Returns:
        list[tuple[str,int,int]]: list of tuples: contig, bin1 number, bin2 number
    """
    # initialize a list of arguments for __evaluateOneBinPair
    out = list()

    # for each contig in the genome
    for contig in binned.keys():
        # sort bins by their start position
        sortedBins = sorted(binned[contig].keys(), key=lambda x: min(binned[contig][x][0].start, binned[contig][x][0].end))
        
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


def _formsDimers(fwd:Primer, rev:Primer) -> bool:
    """evaluates if two primers may form primer dimers

    Args:
        fwd (Primer): the forward primer
        rev (Primer): the reverse primer

    Returns:
        bool: indicates if primer dimer formation is possible
    """
    # constant
    FIVE_DEGREES = 5

    # dimer formation is a concern if the Tm is sufficiently high
    return primer3.calc_heterodimer_tm(str(fwd), str(rev)) >= (min(fwd.Tm, rev.Tm) - FIVE_DEGREES)


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

    # unsuitable if primer dimers can form
    if _formsDimers(fwd, rev): return False, pcrLen
    
    return True, pcrLen


def __evaluateOnePair(fwd:Primer, rev:Primer, minPcr:int, maxPcr:int, maxTmDiff:float, binPair:tuple[str,int,int]) -> tuple[Primer,Primer,int,tuple[str,int,int]]:
    """evaluates a primer pair; designed to work in parallel; only returns if the pair passes evaluation

    Args:
        fwd (Primer): the forward Primer
        rev (Primer): the reverse Primer
        minPcr (int): the minimum PCR product size
        maxPcr (int): the maximum PCR product size
        maxTmDiff (float): the maximum difference between primer Tm
        binPair (tuple[str,int,int]): contig, bin1 number, bin2 number

    Returns:
        tuple[Primer,Primer,int,int,int]: forward primer, reverse primer, pcr product size, bin pair
    """
    # determine the product size and if the pair is suitable for PCR
    acceptablePair, pcrLen = __isPairSuitable(fwd, rev, minPcr, maxPcr, maxTmDiff)
    
    # return acceptable pairs and their pcr product sizes
    if acceptablePair: return fwd,rev.reverseComplement(),pcrLen,binPair


def __getCandidatePrimerPairs(binPairs:list[tuple[str,int,int]], bins:dict[str,dict[int,list[Primer]]], params:Parameters) -> list[tuple[Primer,Primer,int,tuple[str,int,int]]]:
    """gets candidate primer pairs from pairs of bins

    Args:
        binPairs (list[tuple[str,int,int]]): the list produced by __getBinPairs
        bins (dict[str,dict[int,list[Primer]]]): the dictionary produced by __binCandidateKmers
        params (Parameters): a Parameters object

    Returns:
        list[tuple[Primer,Primer,int,tuple[str,int,int]]]: a list of primer pairs, the corresponding pcr product size, and the bin pair (contig, bin1, bin2)
    """
    # constants
    FWD = 'forward'
    REV = 'reverse'
    GC = {"G", "C"}
    
    # helper function for evaluating primers
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
            # make sure the primer is on the plus strand
            if fwd.strand == Primer.MINUS:
                fwd = fwd.reverseComplement()
                
            # only evaluate if the 3' end is GC
            if isThreePrimeGc(fwd):
                # for each reverse primer
                for rev in bin2:
                    # make sure the primer is on the plus strand
                    if rev.strand == Primer.MINUS:
                        rev = rev.reverseComplement()

                    # only evaluate if the 3' end is GC
                    if isThreePrimeGc(rev, REV):
                        # save the pair for evaluation
                        args.append((fwd, rev, params.minPcr, params.maxPcr, params.maxTmDiff, (contig, num1, num2)))
    
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


def __getAllSharedPrimerPairs(firstName:str, candidateKmers:dict[str,dict[str,list[Primer]]], candidatePairs:list[tuple[Primer,Primer,int,tuple[str,int,int]]], params:Parameters) -> dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]:
    """gets all the primer pairs that are shared in all the genomes

    Args:
        firstName (str): the name of the genome that has already been evaluated
        candidateKmers (dict[str,dict[str,list[Primer]]]): key=genome name; val=dict: key=contig name; val=list of candidate kmers
        candidatePairs (list[tuple[Primer,Primer,int,tuple[str,int,int]]]): the list produced by __evaluateBinPairs
        params (Parameters): a Parameters object

    Returns:
        dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]: key=pair of Primers; val=dict: key=genome name; val=tuple: contig name, PCR product length, bin pair (contig, bin1, bin2)
    """
    # initialize variables
    k1:Primer
    k2:Primer
    out = dict()
    allowedLengths = range(params.minPcr, params.maxPcr+1)
    
    # get a set of genomes that need to be evaluated (the first name has already been evaluated)
    remaining = set(candidateKmers.keys())
    remaining.remove(firstName)

    # for each remaining genome
    kmers = dict()
    for name in remaining: 
        # restructure the data for O(1) lookup of primers by their sequences
        kmers[name] = __restructureCandidateKmerData(candidateKmers[name])
        
    # for each pair of primers in the candidate pairs
    for p1,p2,length,binPair in candidatePairs:
        # store the data for the genome (firstName) that has already been evaluated
        out[(p1,p2)] = dict()
        out[(p1,p2)][firstName] = (p1.contig,length,binPair)
        
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
                if k1_rev != k2_rev:
                    # save values when k1 is the foward primer and k2 is the reverse primer
                    if k1.start < k2.start:
                        fwd = k1
                        rev = k2
                    
                    # save values when k1 is the reverse primer and k2 is the forward primer
                    else:
                        fwd = k2
                        rev = k1
                    
                    # strand mismatches only occur if primer pairs in firstName are not suitable for the current genome; skip
                    if fwd.strand == rev.strand:
                        # make sure both primers are on the (+) strand before calculating the product length
                        if fwd.strand == Primer.MINUS:
                            fwd = fwd.reverseComplement()
                            rev = rev.reverseComplement()
                        
                        # only save the pair if the length falls within the acceptable range
                        length = rev.end - fwd.start + 1
                        if length in allowedLengths:
                            out[(p1,p2)][name] = (fwd.contig, length, binPair)

    # remove any pairs that are not universally suitable in all genomes
    for pair in set(out.keys()):
        if len(out[pair]) < len(candidateKmers.keys()):
            out.pop(pair)
    
    return out


def _getPrimerPairs(candidateKmers:dict[str,dict[str,list[Primer]]], params:Parameters) -> dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]:
    """gets primer pairs found in all the ingroup genomes

    Args:
        candidateKmers (dict[str,dict[str,list[Primer]]]): key=genome name; val=dict: key=contig; val=list of Primers
        params (Parameters): a Parameters object

    Raises:
        RuntimeError: unable to identify suitable primer pairs from the candidate kmers
        RuntimeError: unable to identify primer pairs in every ingroup genome

    Returns:
        dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]: key=Primer pairs; val=dict: key=genome name; val=tuple: contig, pcr product size, bin pair (contig, bin1, bin2)
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
        binnedCandidateKmers = __binCandidateKmers(candidateKmers[name])
        
        # get the bins that could work
        binPairs = __getBinPairs(binnedCandidateKmers, params.minLen, params.minPcr, params.maxPcr)

        # get candidate primer pairs
        candidatePairs = __getCandidatePrimerPairs(binPairs, binnedCandidateKmers, params)
        
        # get the primer pairs that are shared in all genomes
        pairs = __getAllSharedPrimerPairs(name, candidateKmers, candidatePairs, params)
        
        # keep track of all the candidate pairs and final pairs
        allCand.extend(candidatePairs)
        out.update(pairs)
    
    # make sure that some candidate pairs were identified
    if allCand == []:
        params.log.rename(_getPrimerPairs.__name__)
        params.log.error(ERR_MSG_1)
        raise RuntimeError(ERR_MSG_1)
    
    # make sure that some candidate pairs were universal
    elif out == dict():
        params.log.rename(_getPrimerPairs.__name__)
        params.log.error(ERR_MSG_2)
        raise RuntimeError(ERR_MSG_2)
    
    return out


def _keepOnePairPerBinPair(pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]], name:str) -> None:
    """keeps only one primer pair per bin pair

    Args:
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]): the dictionary created by __getAllSharedPrimerPairs
        name (str): the name of the genome used to bin kmers
    
    Returns:
        does not return; modifies input
    """
    # initialize a set to track the pairs that have already been seen
    seen = set()
    
    # go through each pair in the dictionary
    for pair in set(pairs.keys()):
        # extract the bin pair for the reference genome
        binPair = pairs[pair][name][-1]
        
        # if the bin pair has already been seen, then remove this primer pair
        if binPair in seen and binPair != (): # outgroup binpairs are empty tuples
            pairs.pop(pair)
        
        # otherwise mark the bin pair as seen
        else:
            seen.add(binPair)
