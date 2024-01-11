import multiprocessing
from Bio.Seq import Seq
from bin.Primer import Primer
from bin.utilities import parseArgs
from bin.getCandidateKmers import getAllCandidatePrimers



def __binOverlappingPrimers(candidates:dict[str,list[Primer]]) -> dict[str, dict[int, list[Primer]]]:
    """bins primers that overlap to form a contiguous sequence based on positional data

    Args:
        candidates (dict[str,list[Primer]]): the dictionary produced by __getAllCandidatePrimers

    Returns:
        dict[str, dict[int, list[Primer]]]: key = contig; val=dict[key = bin number; val = list of overlapping Primer objects]
    """
    # initialize output
    bins = dict()
    
    # for each contig in the dictionary
    for contig in candidates:
        currentBin = 0
        prevEnd = None
        bins[contig] = dict()
        
        # for each candidate primer in the contig
        for cand in candidates[contig]:
            # the first candidate needs to go into its own bin
            if prevEnd is None:
                bins[contig][currentBin] = [cand]
                prevEnd = cand.end
            
            # if the candidate overlaps the previous one, then add it to the bin
            elif cand.start < prevEnd:
                prevEnd = cand.end
                bins[contig][currentBin].append(cand)
            
            # otherwise the candidate belongs in a new bin
            else:
                prevEnd = cand.end
                currentBin += 1
                bins[contig][currentBin] = [cand]
    
    return bins


def __minimizeOverlaps(bins:dict[str, dict[int,list[Primer]]], minimizerSize:int) -> None:
    """splits excessively long overlap bins by the minimizer sequences of the binned primers
    mutates input dictionary, does not return

    Args:
        bins (dict[str, dict[int,list[Primer]]]): the dictionary produced by __binOverlappingPrimers
        minimizerSize (int): the length of the minimizer sequence
    """
    def overlapIsTooLong(primers:list[Primer]) -> bool:
        """determines if the overlap covered in the bin is too long"""
        MAX_LEN = 64
        length = primers[-1].end - primers[0].start
        return length > MAX_LEN
    
    # for each bin number
    for contig in bins:
        curBin = max(bins[contig].keys()) + 1
        for binNum in list(bins[contig].keys()):
            if overlapIsTooLong(bins[contig][binNum]):
                
                primers = bins[contig].pop(binNum)
                
                # for each primer; cluster them by minimizer sequence
                clusters = dict()
                for primer in primers:
                    minimizer = primer.getMinimizer(minimizerSize)
                    
                    clusters[minimizer] = clusters.get(minimizer, list())
                    clusters[minimizer].append(primer)
                
                # add each cluster as a new bin in the input dictionary
                for clust in clusters.values():
                    bins[contig][curBin] = clust
                    curBin += 1


def _binCandidatePrimers(candidates:dict[str,list[Primer]], minPrimerSize:int) -> dict[str, dict[int, list[Primer]]]:
    """bins candidate primers based on overlapping positions and minimizers if the overlap is too long

    Args:
        candidates (dict[str,list[Primer]]): the dictionary produced by __getAllCandidatePrimers
        minPrimerSize (int): the minimum primer length

    Returns:
        dict[str, dict[int, list[Primer]]]: key = contig; val = dict[key = bin number; val = list of Primer objects]
    """
    # first bin overlapping primers
    binned = __binOverlappingPrimers(candidates)
    
    # next divide the overlapping primers by minimizers (window size 1/2 min primer len)
    __minimizeOverlaps(binned, int(minPrimerSize / 2))
    
    return binned


def __evaluateOneBinPair(bin1:list[Primer], bin2:list[Primer], maxTmDiff:float, minProdLen:int, maxProdLen:int, queue) -> None:
    """evaluates a pair of bins of primers to find a single suitable pair
    designed to run in parallel

    Args:
        bin1 (list[Primer]): the first (upstream) bin of primers
        bin2 (list[Primer]): the second (downstream) bin of primers
        maxTmDiff (float): the maximum difference in melting temps
        minProdLen (int): the minimum PCR product length
        maxProdLen (int): the maximum PCR product length
        queue (_type_): a multiprocessing.Manager().Queue() object to store results for writing
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
                # only proceed if the product length is within allowable limits
                productLen = p2.end - p1.start
                if minProdLen <= productLen <= maxProdLen:
                    # only proceed if the reverse primer's end is GC
                    if isThreePrimeGc(p2, REV):
                        # only proceed if the Tm difference is allowable
                        if isTmDiffWithinRange(p1, p2):
                            # only proceed if primer dimers won't form
                            if noPrimerDimer(p1, p2):
                                # write the suitable pair; stop comparing other primers in these bins
                                queue.put((p1,p2.reverseComplement(),productLen))
                                return


def __writePairsToFile(fn:str, terminator:str, queue) -> None:
    """writes primer pairs to file in parallel

    Args:
        fn (str): the filename for writing the results
        terminator (str): the string that will tell this process to terminiate
        queue: a multiprocessing.Manager().Queue() object to store results for writing
    """
    # constants
    SEP = '\t'
    EOL = '\n'
    NUM_DEC = 1
    HEADER = ['contig',
              'fwd_seq',
              'fwd_start',
              'fwd_end',
              'fwd_Tm',
              'fwd_GC',
              'rev_seq',
              'rev_start',
              'rev_end',
              'rev_Tm',
              'rev_GC',
              'product_len']
    
    with open(fn, 'w') as fh:
        # write the header
        fh.write(SEP.join(HEADER) + EOL)
        
        # keep checking the queue
        while True:
            # extract the value from the queue
            val = queue.get()
            
            # stop loooping when the terminate signal is received
            if val == terminator:
                break
            
            # extract the primers and PCR product length from the queue
            p1:Primer
            p2:Primer
            length:int
            p1,p2,length = queue.get()
            
            # write the row to file
            row = [p1.contig,
                   str(p1.seq),
                   str(p1.start),
                   str(p1.end),
                   str(round(p1.Tm, NUM_DEC)),
                   str(round(p1.gcPer, NUM_DEC)),
                   str(p2.seq),
                   str(p2.start),
                   str(p2.end),
                   str(round(p2.Tm, NUM_DEC)),
                   str(round(p2.gcPer, NUM_DEC)),
                   str(length)]
            fh.write(SEP.join(row) + EOL)
            fh.flush()


def _getPrimerPairs(primersD:dict[str, dict[int, list[Primer]]], minPrimerLen:int, minProdLen:int, maxProdLen:int, maxTmDiff:float, outFN:str, numThreads:int) -> None:
    """picks primer pairs from a collection of binned candidate primers and write them to file

    Args:
        primersD (dict): the dictionary produced by __binCandidatePrimers
        minPrimerLen (int): minimum primer length
        minProdLen (int): min PCR product length
        maxProdLen (int): max PCR product length
        maxTmDiff (float): maximum difference in melting temps between primers
        outFN (str): filename to write the results
        numThreads (int): number of threads for parallel processing
    """
    # constants
    TERMINATOR = 'END'
    
    # messages
    MSG_1 = "choosing pairs of primer clusters to compare"
    MSG_2 = "evaluating pairs of primer clusters in parallel"
    END = ' ... '
    DONE = 'done.'
    
    # make a queue so results can be written to the file in parallel
    queue = multiprocessing.Manager().Queue()

    # intitialize the pool and start the writer process
    pool = multiprocessing.Pool(processes=numThreads)
    pool.apply_async(__writePairsToFile, (outFN, TERMINATOR, queue))

    # initialize a list of arguments for __evaluateOneBinPair
    argsL = list()
    
    # print status
    print(MSG_1, end=END, flush=True)
    
    # build a list of arguments to process in parallel
    # process each contig separately
    for contig in primersD.keys():
        # a list of the bin numbers for this contig
        bins = list(primersD[contig].keys())
        
        # get the first primer bin
        for idx in range(len(bins)-1):
            bin1 = primersD[contig][bins[idx]]
            
            # get the second primer bin
            for jdx in range(idx+1,len(bins)):
                bin2 = primersD[contig][bins[jdx]]
                
                # calculate the smallest possible product length
                # smallest primer at beginning of b2 minus smallest primer at end of b1
                smallest = (bin2[0].start + minPrimerLen) - (bin1[-1].end - minPrimerLen)
                
                # calculate the largest posisble product length
                # difference between the two extremes
                largest = bin2[-1].end - bin1[0].start
                
                # done with all subsequent bins if smallest possible product is too large
                if smallest > maxProdLen:
                    break
                
                # move to next bin2 if largest possible product is too small
                elif largest < minProdLen:
                    continue
                
                # otherwise, these two bins could produce viable PCR product lengths; compare them
                else:
                    argsL.append((bin1, bin2, maxTmDiff, minProdLen, maxProdLen, queue))

    # print status
    print(DONE)
    print(MSG_2, end=END, flush=True)

    # evaluate pairs of binned primers and write them in parallel
    pool.starmap(__evaluateOneBinPair, argsL)
    
    # stop the writer and close the pool
    queue.put(TERMINATOR)
    pool.close()
    pool.join()
    
    # print status
    print(DONE)


def __main():
    """main runner function
        * parses command line arguments
        * retrieves all the candidate primers
        * bins the candidate primers to reduce unnecessary comparisons
        * finds suitable primer pairs
        * writes suitable primer pairs to file
    """
    # parse command line arguments
    inFN,outFN,frmt,minPrimerLen,maxPrimerLen,minGc,maxGc,minTm,maxTm,minPcrLen,maxPcrLen,maxTmDiff,numThreads,helpRequested = parseArgs()
    
    # no work to do if help was requested
    if not helpRequested:
        # get the candidate primer sequences
        candidatesD = _getAllCandidatePrimers(inFN,
                                            frmt,
                                            minPrimerLen,
                                            maxPrimerLen,
                                            minGc,
                                            maxGc,
                                            minTm,
                                            maxTm,
                                            numThreads)
        
        # bin the primers by overlap and (sometimes) maximizer sequences
        # maximizer window size is 1/2 the minimum primer length
        binned = _binCandidatePrimers(candidatesD, int(minPrimerLen/2))

        # evaluate binned primer pairs and write suitable primer pairs to file in parallel
        _getPrimerPairs(binned,
                         minPrimerLen,
                         minPcrLen,
                         maxPcrLen,
                         maxTmDiff,
                         outFN,
                         numThreads)
