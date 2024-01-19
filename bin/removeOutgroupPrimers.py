import multiprocessing
from bin.Clock import Clock
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from bin.getCandidateKmers import _kmpSearch


def __getPcrLen(p1:Primer, p2:Primer, contig:SeqRecord) -> tuple[Primer,Primer,int]:
    """gets the pcr length of a primer pair for a given contig using substring searches

    Args:
        p1 (Primer): the forward primer
        p2 (Primer): the reverse primer
        contig (SeqRecord): the contig to amplify using the primer pairs

    Returns:
        tuple[Primer,Primer,int]: a pair of primers and the corresponding pcr product length
    """
    # search the contig for the first primer sequence
    found1,start1 = _kmpSearch(contig.seq, p1.seq)
    
    # nothing to do unless first primer found
    if found1:
        # search the contig for the second primer (it is rc; need to undo that)
        found2,start2 = _kmpSearch(contig.seq, p2.seq.reverse_complement())
        
        # if a second primer is found
        if found2:
            # return the length of the pcr product
            if start1 > start2:
                pcrLen = start1 - start2 + len(p1)
            else:
                pcrLen = start2 - start1 + len(p2)
            
            return p1,p2,pcrLen


def _removeOutgroupPrimers(outgroup:dict[str,list[SeqRecord]], pairs:dict[tuple[Primer,Primer],dict[str,int]], disallowedLens:range, numThreads:int) -> None:
    """removes any primer pairs that produce PCR products of the disallowed lengths

    Args:
        outgroup (dict[str,list[SeqRecord]]): key=genome name; val=list of contigs as SeqRecord objects
        pairs (dict[tuple[Primer,Primer],dict[str,int]]): key=Primer pairs; val=dict: key=genome name; val=tuple: contig, PCR product length
        disallowedLens (range): the range of pcr product lengths that are not allowed for the outgroup
        numThreads (int): the number of threads available for parallel processing
    """
    # for each contig
    for name in outgroup.keys():
        for contig in outgroup[name]:
            # make sure the contig sequence is capitalized
            contig.seq = contig.seq.upper()
            
            # for each primer pair (reevaluate each iteration as it may change)
            pairsToCheck = set(pairs.keys())
            
            # add each pair to check to a list of arguments
            args = list()
            for p1,p2 in pairsToCheck:
                args.append((p1, p2, contig))
            
            # process the pairs in parallel for this contig
            with multiprocessing.Pool(numThreads) as pool:
                results = pool.starmap(__getPcrLen, args)
            
            # evaluate each result
            for p1,p2,pcrLen in [r for r in results if r is not None]:
                # remove the pair from the dictionary if the length is not allowed
                if pcrLen in disallowedLens:
                    pairs.pop((p1,p2))
                
                # otherwise save length in the dictionary
                elif pcrLen > 0:
                    pairs[(p1,p2)][name] = (contig,pcrLen)
