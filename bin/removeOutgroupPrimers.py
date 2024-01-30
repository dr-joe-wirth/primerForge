from Bio.Seq import Seq
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters


def __getAllKmers(contig:SeqRecord, minLen:int, maxLen:int) -> dict[Seq,list[int]]:
    """gets all the kmers and their start positions from a contig

    Args:
        contig (SeqRecord): a contig as a SeqRecord object
        minLen (int): the minimum kmer length
        maxLen (int): the maximum kmer length

    Returns:
        dict[Seq,list[int]]: key=kmer sequence; val=list of start positions
    """
    # initialize variables
    kmers = dict()
    krange = range(minLen, maxLen + 1)
    smallest = min(krange)
    
    # extract the contig sequences
    fwdSeq:Seq = contig.seq
    
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
                kmer = fwdSeq[start:start+klen]
                kmers[kmer] = kmers.get(kmer, list())
                kmers[kmer].append(start)
        
        # stop iterating through the contig when we're done with it
        if done:
            break

    return kmers


def __getOutgroupProductSizes(kmers:dict[Seq,list[int]], fwd:Primer, rev:Primer) -> set[int]:
    """gets a set of pcr product sizes for a primer pair

    Args:
        kmers (dict[Seq,list[int]]): the dictionary produced by __getAllKmers
        fwd (Primer): the forward primer
        rev (Primer): the reverse primer

    Returns:
        set[int]: a set of pcr product sizes
    """
    # initialize the output
    out = set()
    
    # try to get the primer binding sites from the kmers
    try:
        # (+) strand if forward and reverse sequence are both present
        fStarts = kmers[fwd.seq]
        rStarts = kmers[rev.seq.reverse_complement()]
        reversed = False
    
    except KeyError:
        try:
            # (-) strand if the reverse complements are present
            fStarts = kmers[fwd.seq.reverse_complement()]
            rStarts = kmers[rev.seq]
            reversed = True
        
        # no binding sites ==> empty set
        except KeyError:
            return out
    
    # we found the binding sites; calculate the pcr product sizes
    for f in fStarts:
        for r in rStarts:
            # equation if the binding sites are on the (+) strand
            if not reversed:
                pcrLen = r + len(rev) - f
            
            # equation if the binding sites are on the (-) strand
            else:
                pcrLen = f + len(fwd) - r
            
            # negative products mean the primers are oriented opposite from each other
            if pcrLen > 0:
                out.add(pcrLen)
    
    return out


def _removeOutgroupPrimers(outgroup:dict[str,list[SeqRecord]], pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]], params:Parameters) -> None:
    """removes primers found in the outgroup that produce disallowed product sizes

    Args:
        outgroup (dict[str,list[SeqRecord]]): key=genome name; val=list of contigs
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]): key=Primer pair; dict:key=genome name; val=tuple:contig,pcr product size,bin pair(contig, num1, num2)
        params (Parameters): a Parameters object

    Raises:
        RuntimeError: all primer pairs were present in the outgroup
    """
    # message
    ERR_MSG = "failed to find primer pairs that are absent in the outgroup"
    
    # for each outgroup genome
    for name in outgroup.keys():
        # for each contig in the genome
        for contig in outgroup[name]:
            # get all the kmers in this genome
            kmers = __getAllKmers(contig, params.minLen, params.maxLen)
            
            # recalculate which pairs need to be evaluated
            pairsToCheck = set(pairs.keys())
            for fwd,rev in pairsToCheck:
                # get the outgroup products for this primer pair
                outgroupProducts = __getOutgroupProductSizes(kmers, fwd, rev)
                
                # if there are no products, then the size is 0
                if outgroupProducts == set():
                    pairs[(fwd,rev)][name] = ("NA", 0, ())
                
                # if there are products
                else:
                    # determine if multiple product lengths need to be reported
                    count = 0
                    
                    # for each pcr product length
                    for pcrLen in outgroupProducts:
                        # remove any pairs that produce disallowed product sizes
                        if pcrLen in params.disallowedLens:
                            pairs.pop((fwd,rev))
                        
                        # otherwise count the number of sizes seen
                        else:
                            count += 1
                    
                    # if more than one product size, then join them as a comma-separated list
                    if count > 1:
                        pcrLen = ",".join(map(str, outgroupProducts))
                    
                    # save the pcr product size for this genome (binpair is empty tuple)
                    pairs[(fwd,rev)][name] = (contig.id, pcrLen, ())

    # if the pairs is now empty, then raise an error
    if pairs == dict():
        # save details if debugging
        if params.debug:
            params.log.setLogger(_removeOutgroupPrimers.__name__)
            params.log.writeErrorMsg(ERR_MSG)
        raise RuntimeError(ERR_MSG)
