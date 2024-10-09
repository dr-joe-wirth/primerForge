from Bio.Seq import Seq
from bin.Primer import Primer
from bin.Product import Product
from itertools import combinations
import multiprocessing, os, primer3
from typing import Generator, Union
from bin.Parameters import Parameters


def __reduceBinSize(bins:dict[str,dict[int,list[Primer]]]) -> None:
    """reduces the number of primers in each bin; modifies input

    Args:
        bins (dict[str,dict[int,list[Primer]]]): the dictionary being created within __binCandidateKmers
    """
    # for each contig
    for contig in bins.keys():
        # for each bin
        for binNum in bins[contig].keys():
            # sort the primers by their lengths
            lengths = dict()
            for primer in bins[contig][binNum]:
                lengths[len(primer)] = lengths.get(len(primer), list())
                lengths[len(primer)].append(primer)

            # remove smaller primers that are sub strings
            for small,big in combinations(sorted(lengths.keys()), 2):
                # for each smaller primer
                for p1 in lengths[small]:
                    # for each larger primer
                    for p2 in lengths[big]:
                        # remove any primer that is a substring of an existing primer
                        if p1.seq in p2.seq:
                            bins[contig][binNum].remove(p1)
                            lengths[small].remove(p1)
                            break


def __binCandidateKmers(candidates:dict[str,list[Primer]], maxBinSize:int, isFirstGenome=False) -> dict[str,dict[int,list[Primer]]]:
    """bins primers that overlap to form a contiguous sequence based on positional data

    Args:
        candidates (dict[str,list[Primer]]): key=contig; val=list of Primers (sorted by start position)
        maxBinSize (int): the maximum length (bp) of bins
        isFirstGenome (bool, optional): indicates if this is the first genome. Defaults to False.

    Returns:
        dict[str,dict[int,list[Primer]]]: key=contig; val=dict: key=bin number; val=list of overlapping Primer objects
    """
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
            # force the candidate to the plus strand
            if cand.strand == Primer.MINUS:
                cand = cand.reverseComplement()
            
            # get the current start/end on the (+) strand
            curStart = cand.start
            curEnd = cand.end
            
            # get the length of the current bin
            binLen = curEnd - binStart
            
            # the first candidate needs to go into its own bin
            if prevEnd is None:
                out[contig][currentBin] = [cand]
                prevEnd = curEnd
                binStart = curStart
            
            # if the candidate overlaps the previous one and the bin isn't too long, then add it to the bin
            elif curStart < prevEnd and binLen <= maxBinSize:
                prevEnd = curEnd
                out[contig][currentBin].append(cand)
            
            # otherwise the candidate belongs in a new bin
            else:
                prevEnd = curEnd
                binStart = curStart
                currentBin += 1
                
                out[contig][currentBin] = [cand]
    
    # reduce the number of kmers in a bin for the first genome only
    if isFirstGenome:
        __reduceBinSize(out)
    
    return out


def __getBinPairs(binned:dict[str,dict[int,list[Primer]]], minPrimerLen:int, minProdLen:int, maxProdLen:int) -> list[tuple[str,int,int]]:
    """gets a list of bin pairs that can be used to find primer pairs

    Args:
        binned (dict[str,dict[int,list[Primer]]]): the dictionary produced by __binCandidateKmers
        minPrimerLen (int): the minimum primer length
        minProdLen (int): the minimum pcr product size
        maxProdLen (int): the maximum pcr product size

    Returns:
        list[tuple[str,int,int]]: list of bin pairs (contig, bin1 number, bin2 number)
    """
    # initialize output
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
                
                # otherwise, these two bins could produce viable PCR product lengths
                else:
                    out.append((contig, sortedBins[idx], sortedBins[jdx]))
    return out


def _formsDimers(fwd:str, rev:str, fwdTm:float, revTm:float, mvConc:float, dvConc:float, \
                 dntpConc:float, dnaConc:float, tempC:float, maxLoop:int, tempTolerance:float) -> tuple[bool, float]:
    """evaluates if two primers may form primer dimers

    Args:
        fwd (str): the forward primer sequence
        rev (str): the reverse primer sequence
        fwdTm (float): the forward primer melting temperature
        revTm (float): the reverse primer melting temperature
        mvConc (float): primer3 mv_conc
        dvConc (float): primer3 dv_conc
        dntpConc (float): primer3 dntp_conc
        dnaConc (float): primer3 dna_conc
        tempC (float): primer3 temp_c
        maxLoop (int): primer3 max_loop
        tempTolerance (float): the minimum degrees below primer Tm allowed for secondary structures

    Returns:
        tuple[bool, float]: indicates if primer dimer formation is possible; the dimer Tm
    """
    dimerTm = primer3.calc_heterodimer_tm(fwd,
                                          rev,
                                          mv_conc=mvConc,
                                          dv_conc=dvConc,
                                          dntp_conc=dntpConc,
                                          dna_conc=dnaConc,
                                          temp_c=tempC,
                                          max_loop=maxLoop)
    
    # dimer formation is a concern if the Tm is sufficiently high
    highEnough = dimerTm >= (min(fwdTm, revTm) - tempTolerance)

    return highEnough, dimerTm


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


def __evaluateOnePair(fwd:str, rev:str, binPair:str, fTm:float, rTm:float, mvConc:float, dvConc:float, \
                      dntpConc:float, dnaConc:float, tempC:float, maxLoop:int, tempTolerance:float) -> Union[tuple[str, float],None]:
    """evaluates a primer pair; designed to work in parallel; only returns if the pair passes evaluation

    Args:
        fwd (str): the forward Primer
        rev (str): the reverse Primer
        binPair (str): the bin pair for this pair
        fTm (float): the forward melting temperature
        rTm (float): the reverse melting temperature
        mvConc (float): primer3 mv_conc
        dvConc (float): primer3 dv_conc
        dntpConc (float): primer3 dntp_conc
        dnaConc (float): primer3 dna_conc
        tempC (float): primer3 temp_c
        maxLoop (int): primer3 max_loop
        tempTolerance (float): the minimum degrees below primer Tm allowed for secondary structures

    Returns:
        tuple[str, float] | None: the bin pair for this pair, the heterodimer Tm; None if the pair does not pass evaluation
    """
    # get the dimer Tm and if it forms dimers
    forms, dimerTm = _formsDimers(fwd, rev, fTm, rTm, mvConc, dvConc, dntpConc, dnaConc, tempC, maxLoop, tempTolerance)
    
    # return acceptable pairs and their pcr product sizes
    if not forms: return binPair, dimerTm


def __getCandidatePrimerPairs(binPairs:list[tuple[str,int,int]], bins:dict[str,dict[int,list[Primer]]], params:Parameters) -> list[tuple[Primer,Primer,Product]]:
    """gets candidate primer pairs from pairs of bins

    Args:
        binPairs (list[tuple[str,int,int]]): the list produced by __getBinPairs
        bins (dict[str,dict[int,list[Primer]]]): the dictionary produced by __binCandidateKmers
        params (Parameters): a Parameters object

    Returns:
        list[tuple[Primer,Primer,Product]]: a list of (forward primer, reverse primer, pcr product)
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
    
    # generator function for getting arguments
    def generateArgs() -> Generator[tuple[str, str, float, float, str],None,None]:
        """ generates arguments for __evaluateOnePair
        """
        # for each bin pair
        for contig,num1,num2 in binPairs:
            # indicate if a suitable primer pair has been found
            found = False
            
            # extract the lists of primers for these bins
            bin1 = bins[contig][num1]
            bin2 = bins[contig][num2]
            
            # go through each possible forward primer
            for idx,fwd in enumerate(bin1):
                # only consider fwd primers with 3' GC
                if isThreePrimeGc(fwd):
                    # go through each possible reverse primer
                    for jdx,rev in enumerate(bin2):
                        # only consider rev primers with 3' GC
                        if isThreePrimeGc(rev, REV):
                            # determine if the pair is suitable and its PCR product length
                            suitable,pcrLen = __isPairSuitable(fwd, rev, params.minPcr, params.maxPcr, params.maxTmDiff)
                            
                            if suitable:
                                # cast the bin pair as a string to reduce memory footprint
                                binPair = str((contig, num1, num2))
                                
                                # store the indices of these primers within their respective bins
                                lookup[binPair] = (idx, jdx, pcrLen)
                                
                                # yield data to the child process
                                yield (str(fwd),
                                       str(rev.reverseComplement()),
                                       binPair,
                                       fwd.Tm,
                                       rev.Tm,
                                       params.p3_mvConc,
                                       params.p3_dvConc,
                                       params.p3_dntpConc,
                                       params.p3_dnaConc,
                                       params.p3_tempC,
                                       params.p3_maxLoop,
                                       params.tempTolerance)
                                
                                # indicate that a suitable primer pair has been found and break out of the loop
                                found = True
                                break
                    
                    # move to the next bin pair if a suitable primer pair has been found
                    if found:
                        break
    
    # initialize variables
    out = list()
    lookup = dict()
    
    # evaluate pairs of primers for heterodimer potential in parallel
    pool = multiprocessing.Pool(params.numThreads)
    results = pool.starmap(__evaluateOnePair, generateArgs())
    pool.close()
    pool.join()
    
    # for each valid pair
    for pair,dimerTm in [r for r in results if r is not None]:
        # recast the bin pair from a string to its component parts
        contig, fbin, rbin = eval(pair)
        
        # extract the indices and pcr length for this pair
        idx, jdx, pcrLen = lookup[pair]
        
        # retrieve the actual primer objects for this pair
        fwd:Primer = bins[contig][fbin][idx]
        rev:Primer = bins[contig][rbin][jdx]
        
        # reverse complement the primer to save appropriate data in the product
        rev = rev.reverseComplement()
        
        # save this pair in the output
        out.append((fwd, rev, Product(contig, pcrLen, fbin, rbin, dimerTm, fwd.start, fwd.end, fwd.strand, rev.start, rev.end, rev.strand)))
    
    return out


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


def __getAllSharedPrimerPairs(firstName:str, candidateKmers:dict[str,dict[str,list[Primer]]], candidatePairs:list[tuple[Primer,Primer,Product]], params:Parameters) -> dict[tuple[Primer,Primer],dict[str,Product]]:
    """gets all the primer pairs that are shared in all the genomes

    Args:
        firstName (str): the name of the genome that has already been evaluated
        candidateKmers (dict[str,dict[str,list[Primer]]]): key=genome name; val=dict: key=contig name; val=list of candidate kmers
        candidatePairs (list[tuple[Primer,Primer,Product]): the list produced by __getCandidatePrimerPairs
        params (Parameters): a Parameters object

    Returns:
        dict[tuple[Primer,Primer],dict[str,Product]]: key=pair of Primers; val=dict: key=genome name; val=Product object
    """
    # constant
    TEMP_BIN = -1
    
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
    for p1,p2,product in candidatePairs:
        # store the data for the genome (firstName) that has already been evaluated
        out[(p1,p2)] = dict()
        out[(p1,p2)][firstName] = product
        
        # for each unprocessed genome
        for name in remaining:
            # get the first kmer (may be rev comp) and its Tms
            try:
                k1:Primer = kmers[name][p1.seq]
                k1_rev = False
                k1.hairpinTm = p1.hairpinTm
                k1.rcHairpin = p1.rcHairpin
                k1.homodimerTm = p1.homodimerTm
                k1.rcHomodimer = p1.rcHomodimer
                
            except KeyError:
                k1 = kmers[name][p1.seq.reverse_complement()]
                k1_rev = True
                k1.hairpinTm = p1.rcHairpin
                k1.rcHairpin = p1.hairpinTm
                k1.homodimerTm = p1.rcHomodimer
                k1.rcHomodimer = p1.homodimerTm
            
            # get the second kmer (may be rev comp) and its Tms
            try:
                k2 = kmers[name][p2.seq]
                k2_rev = False
                k2.hairpinTm = p2.hairpinTm
                k2.rcHairpin = p2.rcHairpin
                k2.homodimerTm = p2.homodimerTm
                k2.rcHomodimer = p2.rcHomodimer
    
            except KeyError:
                k2 = kmers[name][p2.seq.reverse_complement()]
                k2_rev = True
                k2.hairpinTm = p2.rcHairpin
                k2.rcHairpin = p2.hairpinTm
                k2.homodimerTm = p2.rcHomodimer
                k2.rcHomodimer = p2.homodimerTm
            
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
                            # put k2 on the opposite strand of k1 before saving product information
                            k2 = k2.reverseComplement()
                            
                            # save product information
                            ## use k1 and k2 so the sequences matche p1 and p2
                            ## the bin numbers are incorrect and will need to be fixed later
                            out[(p1,p2)][name] = Product(k1.contig, length, TEMP_BIN, TEMP_BIN, product.dimerTm, k1.start, k1.end, k1.strand, k2.start, k2.end, k2.strand)

    # remove any pairs that are not universally suitable in all genomes
    for pair in set(out.keys()):
        if len(out[pair]) < len(candidateKmers.keys()):
            del out[pair]
    
    return out


def __updateBinsForUnprocessedGenomes(name:str, maxBinSize:int, kmers:dict[str,list[Primer]], pairs:dict[tuple[Primer,Primer],dict[str,Product]]) -> None:
    """updates the bin pairs for the unprocessed genomes to ensure that all redundant primers are removed

    Args:
        name (str): the name of an unprocessed genome
        maxBinSize (int): the maximum length (bp) of bins
        kmers (dict[str,list[Primer]]): candidate kmers for this genome (key=contig; val=list of candidate kmers)
        pairs (dict[tuple[Primer,Primer],dict[str,Product]]): the dictionary produced by __getAllSharedPrimerPairs
    """
    # bin the candidate kmers for this genome
    binned = __binCandidateKmers(kmers, maxBinSize)
    
    # restructure the data for O(1) lookup of bins by primer sequence
    primerToBin = dict()
    for contig in binned.keys():
        for bin in binned[contig].keys():
            for primer in binned[contig][bin]:
                primerToBin[primer] = (contig, bin)
    
    # for each primer pair (primers are not on the same strand)
    for fwd,rev in pairs.keys():
        # check for presence on forward strand
        try:
            contig, fbin = primerToBin[fwd]
            contig, rbin = primerToBin[rev.reverseComplement()]
        
        # check for presence on the reverse strand
        except KeyError:
            contig, fbin = primerToBin[fwd.reverseComplement()]
            contig, rbin = primerToBin[rev]
        
        # update the bin numbers for this pair
        pairs[(fwd,rev)][name].fbin = fbin
        pairs[(fwd,rev)][name].rbin = rbin


def _getPrimerPairs(candidateKmers:dict[str,dict[str,list[Primer]]], params:Parameters) -> dict[tuple[Primer,Primer],dict[str,Product]]:
    """gets primer pairs found in all the ingroup genomes

    Args:
        candidateKmers (dict[str,dict[str,list[Primer]]]): key=genome name; val=dict: key=contig; val=list of Primers
        params (Parameters): a Parameters object

    Raises:
        RuntimeError: unable to identify suitable primer pairs from the candidate kmers
        RuntimeError: unable to identify primer pairs in every ingroup genome

    Returns:
        dict[tuple[Primer,Primer],dict[str,Product]]: key=Primer pairs; val=dict: key=genome name; val=Product
    """
    # messages
    ERR_MSG_1 = "could not identify suitable primer pairs from the candidate kmers"
    ERR_MSG_2 = "could not identify primer pairs present in every ingroup genome"
    
    # only need to get pairs for one genome
    firstName = os.path.basename(params.ingroupFns[0])
    
    # bin kmers to reduce time complexity
    binnedCandidateKmers = __binCandidateKmers(candidateKmers[firstName], params.maxBinSize, isFirstGenome=True)
    
    # get the bins that could work
    binPairs = __getBinPairs(binnedCandidateKmers, params.minLen, params.minPcr, params.maxPcr)

    # get candidate primer pairs
    candidatePairs = __getCandidatePrimerPairs(binPairs, binnedCandidateKmers, params)
    
    # make sure that some candidate pairs were identified
    if candidatePairs == []:
        params.log.rename(_getPrimerPairs.__name__)
        params.log.error(ERR_MSG_1)
        raise RuntimeError(ERR_MSG_1)
    
    # get the primer pairs that are shared in all genomes
    pairs = __getAllSharedPrimerPairs(firstName, candidateKmers, candidatePairs, params)
    
    # make sure that some candidate pairs were universal
    if pairs == dict():
        params.log.rename(_getPrimerPairs.__name__)
        params.log.error(ERR_MSG_2)
        raise RuntimeError(ERR_MSG_2)
    
    # bin the primers from the unprocessed genomes
    unprocessed = set(candidateKmers.keys()).difference({firstName})
    for name in unprocessed:
        __updateBinsForUnprocessedGenomes(name, params.maxBinSize, candidateKmers[name], pairs)
    
    return pairs
