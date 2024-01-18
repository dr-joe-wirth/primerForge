from bin.getCandidateKmers import _kmpSearch
from bin.Primer import Primer
from bin.Clock import Clock
from Bio.SeqRecord import SeqRecord


def __getPcrLen(p1:Primer, p2:Primer, contig:SeqRecord) -> int:
    found,start1 = _kmpSearch(contig.seq, p1.seq)
    if found:
        found,start2 = _kmpSearch(contig.seq, p2.seq.reverse_complement())
        if found:
            return abs(start1 - start2) + len(p1)


def _removeOutgroupPrimers(outgroup:dict[str,list[SeqRecord]], pairs:dict[tuple[Primer,Primer],dict[str,int]], disallowedLens:range) -> None:
    for contigs in outgroup.values():
        for contig in contigs:
            pairsToCheck = set(pairs.keys())
            for p1,p2 in pairsToCheck:
                if __getPcrLen(p1, p2, contig) in disallowedLens:
                    pairs.pop((p1,p2))


