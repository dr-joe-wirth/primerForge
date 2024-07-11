from bin.Clock import Clock
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from ahocorasick import Automaton
from bin.Parameters import Parameters

# global constant
__NULL_PRODUCT = ("NA", 0, ())


# functions
def __getAutomaton(pairs:list[tuple[Primer,Primer]]) -> Automaton:
    """creates an Aho-Corasick Automaton for substring lookup

    Args:
        pairs (list[tuple[Primer,Primer]]): a list of primer pairs

    Returns:
        Automaton: an Aho-Corasick Automaton
    """
    # initialize output
    out = Automaton()
    
    # add each kmer in the pairs to the automaton once and make the automaton
    for kmer in {str(x) for y in pairs for x in y}:
        out.add_word(kmer, kmer)
    out.make_automaton()

    return out


def __extractKmerData(seq:str, strand:str, kmers:Automaton) -> dict[str, list[int]]:
    """extracts data from a contig sequence for a collection of kmers

    Args:
        seq (str): the contig sequence
        strand (str): the strand of the sequence
        kmers (Automaton): the collection of kmers as an Aho-Corasick Automaton

    Returns:
        dict[str, list[int]]: key=kmer; val=list of start positions
    """
    # initialize output
    out = dict()
    
    # for each kmer in the sequence
    for end,kmer in kmers.iter(seq):
        # extract the start position if on the 
        start = end - len(kmer) + 1
        
        # modify the start position if on the (-) strand
        if strand == Primer.MINUS:
            start = len(seq) - start - 1
        
        # build the output
        out[kmer] = out.get(kmer, list())
        out[kmer].append(start)
    
    return out


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


def __getOutgroupProductSizes(kmers:dict[str,dict[str,list[int]]], fwd:str, rev:str) -> set[int]:
    """gets a set of pcr product sizes for a primer pair

    Args:
        kmers (dict[str,dict[str,list[int]]]): key=strand; val=dict: key=kmer; val=list of start positions
        fwd (str): the forward primer
        rev (str): the reverse primer

    Returns:
        set[int]: a set of pcr product sizes
    """
    # initialize the output
    out = set()
    
    # get the sizes when fwd (+) and rev (-)
    try:
        productSizes = __productSizesFromStartPositions(kmers[Primer.PLUS][fwd], kmers[Primer.MINUS][rev])
        out.update(productSizes)

    # primers may not bind those strands
    except KeyError:
        pass
    
    # get the sizes when fwd (-) and rev (+)
    try:
        productSizes = __productSizesFromStartPositions(kmers[Primer.PLUS][rev], kmers[Primer.MINUS][fwd])
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

            # if there is only one primer size, then save it
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
        outgroup (dict[str,list[SeqRecord]]): key=genome name; val=list of contigs (or a generator)
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]): key=Primer pair; dict:key=genome name; val=tuple: contig, pcr product size, bin pair (contig, num1, num2)
        params (Parameters): a Parameters object

    Raises:
        RuntimeError: all candidate primer pairs were present in the outgroup
    """
    # messages
    GAP = " "*4
    DONE = "done "
    MSG_1   = "removing primer pairs present in the outgroup sequences"
    MSG_2   = f"{GAP}getting outgroup PCR products"
    MSG_3   = f"{GAP}filtering primer pairs"
    MSG_4   = f"{GAP}processing outgroup results"
    MSG_5A  = "removed "
    MSG_5B  = " pairs after processing "
    MSG_5C  = " ("
    MSG_5D  = " pairs remaining)"
    ERR_MSG = "failed to find primer pairs that are absent in the outgroup"
    
    # initialize variables
    clock = Clock()
    prevName = None
    outgroupKmers = dict()
    outgroupProducts = dict()
    startNumPairs = len(pairs)
    
    # print status and log
    params.log.rename(_removeOutgroupPrimers.__name__)
    params.log.info(MSG_1)
    print(MSG_1)
    clock.printStart(MSG_2)
    params.log.info(MSG_2)
    
    # create the automaton from the pairs
    kmers = __getAutomaton(pairs.keys())
    
    # for each outgroup genome
    for name in outgroup.keys():
        # initialize a dictionary for the current outgroup genome
        outgroupKmers[name] = dict()
        
        # add each contig in the genome to the argument list
        for contig in outgroup[name]:
            outgroupKmers[name][contig.id] = {Primer.PLUS:  dict(),
                                              Primer.MINUS: dict()}
            
            # extract the forward and reverse sequences
            fwd = str(contig.seq)
            rev = str(contig.seq.reverse_complement())
            
            # extract kmer data
            outgroupKmers[name][contig.id][Primer.PLUS].update((__extractKmerData(fwd, Primer.PLUS, kmers)))
            outgroupKmers[name][contig.id][Primer.MINUS].update(__extractKmerData(rev, Primer.MINUS, kmers))
    
    # print status and log
    clock.printDone()
    params.log.info(f"{DONE}{clock.getTimeString()}")
    clock.printStart(MSG_3)
    params.log.info(MSG_3)
    
    # for each outgroup genome
    prevName = next(iter(outgroupKmers.keys()))
    for name in outgroupKmers.keys():
        # keep track of how the pairs are changing after each outgroup genome
        if name != prevName:
            # log the number of pairs removed and remaining if debugging
            params.log.debug(f"{MSG_5A}{startNumPairs - len(pairs)}{MSG_5B}{prevName}{MSG_5C}{len(pairs)}{MSG_5D}")
            
            # reset the starting number and update the prev name
            startNumPairs = len(pairs)
            prevName = name
        
        # create the sub dictionary
        outgroupProducts[name] = dict()
        
        # for each contig
        for contig in outgroupKmers[name].keys():
            # evaluate the pairs present in the dictionary (list allows on-the-fly removal)
            for fwd,rev in list(pairs.keys()):
                # initialize an empty set if one does not already exist
                outgroupProducts[name][(fwd,rev)] = outgroupProducts[name].get((fwd,rev), set())

                # get the outgroup products for this primer pair
                products = __getOutgroupProductSizes(outgroupKmers[name][contig], str(fwd), str(rev))
                
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
                            del pairs[(fwd,rev)]
                            done = True
                            break
                    
                    # if the pcr product lengths are not disallowed, then save them in the dictionary
                    if not done:
                        outgroupProducts[name][(fwd,rev)].update({(contig, x, ()) for x in products})
        
    # log the number of pairs removed and remaining from the last genome if debugging
    params.log.debug(f"{MSG_5A}{startNumPairs - len(pairs)}{MSG_5B}{name}{MSG_5C}{len(pairs)}{MSG_5D}")
    
    # print status; log if debugging
    clock.printDone()
    params.log.info(f"{DONE}{clock.getTimeString()}")
    
    # if the pairs dictionary is now empty, then raise an error
    if pairs == dict():
        params.log.error(ERR_MSG)
        raise RuntimeError(ERR_MSG)
    
    # print status and log
    params.log.info(MSG_4)
    clock.printStart(MSG_4)
    
    # process the outgroup results and add them to the pairs dictionary
    __processOutgroupResults(outgroupProducts, pairs)
    
    # print status; log if debugging
    clock.printDone()
    params.log.info(f"{DONE}{clock.getTimeString()}")
