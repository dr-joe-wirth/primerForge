from bin.Clock import Clock
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from bin.getCandidateKmers import _kmpSearch


def __getPcrLen(p1:Primer, p2:Primer, contig:SeqRecord) -> int:
    # search the contig for the first primer sequence
    found,start1 = _kmpSearch(contig.seq, p1.seq)
    
    # nothing to do unless first primer found
    if found:
        # search the contig for the second primer (it is rc; need to undo that)
        found,start2 = _kmpSearch(contig.seq, p2.seq.reverse_complement())
        
        # if a second primer is found
        if found:
            # return the length of the pcr product
            if start1 > start2:
                return abs(start1 - start2) + len(p2)
            else:
                return abs(start1 - start2) + len(p1)
    
    # return a 
    return -1


def _removeOutgroupPrimers(outgroup:dict[str,list[SeqRecord]], pairs:dict[tuple[Primer,Primer],dict[str,int]], disallowedLens:range) -> None:
    # for each contig
    for name in outgroup.keys():
        for contig in outgroup[name]:
            # for each pair in the dict (reevaluate each iteration)
            pairsToCheck = set(pairs.keys())
            for p1,p2 in pairsToCheck:
                # get the pcr product length
                pcrLen = __getPcrLen(p1, p2, contig)
                
                # remove the pair from the dictionary if the length is bad
                if pcrLen in disallowedLens:  # maybe parallelize this function?  # if so, have to parallelize each time I check contigs
                    pairs.pop((p1,p2))
                
                # otherwise save length in the dictionary
                else:
                    pairs[(p1,p2)][name] = (contig,pcrLen)
