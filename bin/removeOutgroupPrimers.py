from Bio.Seq import Seq
from bin.Clock import Clock
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters

# global constant
__NULL_PRODUCT = ("NA", 0, ())


def __getAllKmers(contig:SeqRecord, minLen:int, maxLen:int) -> dict[str,dict[Seq,list[int]]]:
    """gets all the kmers and their start positions from a contig

    Args:
        contig (SeqRecord): a contig as a SeqRecord object
        minLen (int): the minimum kmer length
        maxLen (int): the maximum kmer length

    Returns:
        dict[str,dict[Seq,list[int]]]: key=strand; val=dict: key=kmer sequence; val=list of start positions
    """
    # initialize variables
    kmers = {Primer.PLUS:  dict(),
             Primer.MINUS: dict()}
    krange = range(minLen, maxLen + 1)
    smallest = min(krange)
    
    # extract the contig sequences
    fwdSeq:Seq = contig.seq
    revSeq:Seq = contig.seq.reverse_complement()
    
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
                # extract and save the forward kmer sequence
                fKmer = fwdSeq[start:start+klen]
                kmers[Primer.PLUS][fKmer] = kmers[Primer.PLUS].get(fKmer, list())
                kmers[Primer.PLUS][fKmer].append(start)
                
                # calculate the end position (rev start position)
                end = start + klen
                
                # extract the reverse kmer sequence
                if start == 0:
                    rKmer = revSeq[-end:]
                else:
                    rKmer = revSeq[-end:-start]
                
                # save the reverse kmer sequence
                kmers[Primer.MINUS][rKmer] = kmers[Primer.MINUS].get(rKmer, list())
                kmers[Primer.MINUS][rKmer].append(end-1)
        
        # stop iterating through the contig when we're done with it
        if done:
            break

    return kmers


def __productSizesFromStartPositions(plusStarts:list[int], minusStarts:list[int]) -> set[int]:
    """calculates the product sizes from lists of primer start positions

    Args:
        plusStarts (list[int]): the start positions on the plus strand
        minusStarts (list[int]): the start positions on the minus strand

    Returns:
        set[int]: a set of pcr product sizes
    """
    # initialize output
    out = set()
    
    # for each pair of start positions on opposite strands
    for fStart in plusStarts:
        for rStart in minusStarts:
            # calculate the PCR product length
            pcrLen = rStart - fStart + 1
            
            # negative values indicate primers that are facing away from each other
            if pcrLen > 0:
                out.add(pcrLen)
    
    return out


def __getOutgroupProductSizes(kmers:dict[str,dict[Seq,list[int]]], fwd:Primer, rev:Primer) -> set[int]:
    """gets a set of pcr product sizes for a primer pair

    Args:
        kmers (dict[str,dict[Seq,list[int]]]): the dictionary produced by __getAllKmers
        fwd (Primer): the forward primer
        rev (Primer): the reverse primer

    Returns:
        set[int]: a set of pcr product sizes
    """
    # initialize the output
    out = set()
    
    # get the sizes when fwd (+) and rev (-)
    try:
        productSizes = __productSizesFromStartPositions(kmers[Primer.PLUS][fwd.seq], kmers[Primer.MINUS][rev.seq])
        out.update(productSizes)

    # primers may not bind those strands
    except KeyError:
        pass
    
    # get the sizes when fwd (-) and rev (+)
    try:
        productSizes = __productSizesFromStartPositions(kmers[Primer.PLUS][rev.seq], kmers[Primer.MINUS][fwd.seq])
        out.update(productSizes)
    
    # primers may not bind those strands
    except KeyError:
        pass
    
    return out


def __processOutgroupResults(outgroupProducts:dict[str,dict[tuple[Primer,Primer],set[tuple[str,int,tuple]]]], pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]) -> None:
    """adds the outgroup results to the pairs dictionary

    Args:
        outgroupProducts (dict[str,dict[tuple[Primer,Primer],set[tuple[str,int,tuple]]]]): key=genome name; val=dict: key=primer pair; val=set of tuples: contig, pcrLen, binpair
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]): key=Primer pair; val=dict: key=genome name; val=tuple: contig, pcrLen, binpair
    
    Returns:
        does not return. modifies the pairs dictionary
    """
    # for each pair remaining to process
    for pair in pairs.keys():
        # for each outgroup genome
        for name in outgroupProducts.keys():
            # extract the outgroup product sizes for this pair
            result = outgroupProducts[name][pair]

            # if there is only one primer seize, then save it
            if len(result) == 1:
                pairs[pair][name] = result.pop()
            
            # otherwise
            else:
                # remove any null products from the set
                try: result.remove(__NULL_PRODUCT)
                except KeyError: pass
                
                # if there is only one primer size, then save it
                if len(result) == 1:
                    pairs[pair][name] = result.pop()
                
                # otherwise
                else:
                    # combine all contigs and pcrLens into separate lists
                    contigs = list()
                    pcrLens = list()
                    for contig,pcrLen,binPair in result:
                        contigs.append(contig)
                        pcrLens.append(pcrLen)
                    
                    # convert the contigs and lengths to comma-separated strings; add empty bin pair
                    pairs[pair][name] = (",".join(contigs), ",".join(map(str,pcrLens)), ())


def _removeOutgroupPrimers(outgroup:dict[str,list[SeqRecord]], pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]], params:Parameters) -> None:
    """removes primers found in the outgroup that produce disallowed product sizes

    Args:
        outgroup (dict[str,list[SeqRecord]]): key=genome name; val=list of contigs
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]): key=Primer pair; dict:key=genome name; val=tuple: contig, pcr product size, bin pair (contig, num1, num2)
        params (Parameters): a Parameters object

    Raises:
        RuntimeError: all candidate primer pairs were present in the outgroup
    """
    # messages
    MSG_1   = "removing primer pairs present in the outgroup sequences"
    MSG_2   = "processing outgroup results"
    MSG_3A  = "removed "
    MSG_3B  = " pairs after processing "
    MSG_3C  = " ("
    MSG_3D  = " pairs remaining)"
    ERR_MSG = "failed to find primer pairs that are absent in the outgroup"
    
    # initialize variables
    clock = Clock()
    outgroupProducts = dict()
    
    # print status and log
    params.log.rename(_removeOutgroupPrimers.__name__)
    params.log.info(MSG_1)
    clock.printStart(MSG_1)
    
    # for each outgroup genome
    for name in outgroup.keys():
        # save the current number of primer pairs
        startNumPairs = len(pairs)
        
        # stop looping if the number of pairs is 0
        if startNumPairs == 0:
            break
        
        # initialize a dictionary for the current outgroup genome
        outgroupProducts[name] = dict()
        
        # for each contig in the genome
        for contig in outgroup[name]:
            # get all the kmers
            kmers = __getAllKmers(contig, params.minLen, params.maxLen)
            
            # evaluate only the pairs that are still present in the dictionary
            for fwd,rev in set(pairs.keys()):
                # initialize an empty set if one does not already exist
                outgroupProducts[name][(fwd,rev)] = outgroupProducts[name].get((fwd,rev), set())

                # get the outgroup products for this primer pair
                products = __getOutgroupProductSizes(kmers, fwd, rev)
                
                # if there are no products, then the size is 0
                if products == set():
                    outgroupProducts[name][(fwd,rev)].update({__NULL_PRODUCT})
                
                # if there are products
                else:
                    # initialize variable to determine if this product needs to be processed further
                    done = False
                    
                    # for each pcr product length
                    for pcrLen in products:
                        # remove any pairs that produce disallowed product sizes
                        if pcrLen in params.disallowedLens:
                            pairs.pop((fwd,rev))
                            done = True
                            break
                    
                    # if the pcr product lengths are not disallowed, then save them in the dictionary
                    if not done:
                        outgroupProducts[name][(fwd,rev)].update({(contig.id, x, ()) for x in products})
        
        # log the number of pairs removed and remaining if debugging
        params.log.debug(f"{MSG_3A}{startNumPairs - len(pairs)}{MSG_3B}{name}{MSG_3C}{len(pairs)}{MSG_3D}")
    
    # print status; log if debugging
    clock.printDone()
    params.log.info(f"done {clock.getTimeString()}")
    
    # if the pairs dictionary is now empty, then raise an error
    if pairs == dict():
        params.log.error(ERR_MSG)
        raise RuntimeError(ERR_MSG)
    
    # print status and log
    params.log.info(MSG_2)
    clock.printStart(MSG_2)
    
    # process the outgroup results and add them to the pairs dictionary
    __processOutgroupResults(outgroupProducts, pairs)
    
    # print status; log if debugging
    clock.printDone()
    params.log.info(f"done {clock.getTimeString()}")
