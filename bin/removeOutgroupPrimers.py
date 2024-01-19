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


def _removeOutgroupPrimers(outgroup:dict[str,list[SeqRecord]], pairs:dict[tuple[Primer,Primer],dict[str,int]], bad) -> None:
    # for each contig
    for contigs in outgroup.values():
        for contig in contigs:
            # for each pair in the dict (reevaluate each iteration)
            pairsToCheck = set(pairs.keys())
            for p1,p2 in pairsToCheck:
                # remove the pair from the dictionar
                if __getPcrLen(p1, p2, contig) is bad:  # maybe parallelize this function?  # what is bad??
                    pairs.pop((p1,p2))
