import os
from Bio import SeqIO
from bin.Clock import Clock
from typing import Iterator
from bin.Primer import Primer
from bin.Product import Product
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters
from bin.kmer_counting.kmer_counter import _encodeKmer, _getAllStartPositionsAndDecodeAllowedEncodings


# global constant
__NULL_PRODUCT = ("NA", 0)


# functions
def __getKmerEncodings(pairs:Iterator[tuple[Primer,Primer]]) -> dict[int,set[int]]:
    """gets the kmer encodings from a collection of primer pairs

    Args:
        pairs (Iterator[tuple[Primer,Primer]]): a collection of primer pairs (eg dictionary keys)

    Returns:
        dict[int,set[int]]: {kmer length: {encodings}}
    """
    # initialize the output
    out = defaultdict(set)

    # create a set of kmer encodings stored under their lengths
    for kmer in {p for pair in pairs for p in pair}:
        out[len(kmer)].add(_encodeKmer(str(kmer)))
    
    return dict(out)


def __extractKmerPositions(seq:str, encodings:dict[int,set[int]]) -> dict[str,dict[str,list[int]]]:
    """extracts kmer positions from a collection of encodings

    Args:
        seq (str): the sequence to evaluate
        encodings (dict[int,set[int]]): {kmer length: {encodings}}

    Returns:
        dict[str,dict[str,list[int]]]: {strand: {kmer: [start positions]}}
    """
    # initialize output
    out = {Primer.PLUS: defaultdict(list),
           Primer.MINUS: defaultdict(list)}

    # for each kmer length
    for k in encodings.keys():
        # decode encodings and get their start positions
        positions = _getAllStartPositionsAndDecodeAllowedEncodings(encodings[k], k, seq)

        # for each kmer
        for kmer,starts in positions.items():
            # for each start position, save the start positions under the correct strand
            for start in starts:
                # positive start indicates plus strand
                if start > 0:
                    out[Primer.PLUS][kmer].append(start)
                
                # negatave start indicates minus strand
                elif start < 0:
                    out[Primer.MINUS][kmer].append(-start)
                
                # zero start that matches beginning of seq indicates plus strand
                elif seq[:len(kmer)] == kmer:
                    out[Primer.PLUS][kmer].append(start)
                
                # zero start that doesn't match beginning of seq indicates minus strand
                else:
                    out[Primer.MINUS][kmer].append(start)
    
    # recast defaultdict to dict
    out[Primer.PLUS] = dict(out[Primer.PLUS])
    out[Primer.MINUS] = dict(out[Primer.MINUS])

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


def __processOutgroupResults(outgroupProducts:dict[str,dict[tuple[Primer,Primer],set[tuple[str,int]]]], pairs:dict[tuple[Primer,Primer],dict[str,Product]]) -> None:
    """adds the outgroup results to the pairs dictionary

    Args:
        outgroupProducts (dict[str,dict[tuple[Primer,Primer],set[tuple[str,int]]]]): key=genome name; val=dict: key=primer pair; val=set of tuples: contig, pcrLen
        pairs (dict[tuple[Primer,Primer],dict[str,Product]]): key=Primer pair; val=dict: key=genome name; val=Product
    
    Returns:
        does not return. modifies the pairs dictionary
    """
    # constant
    FAKE_BIN = -1
    FAKE_COORD = -1
    FAKE_STRAND = ''
    
    # get the name of an ingroup genome
    ingroupName = next(iter([n for v in pairs.values() for n in v.keys() if n not in outgroupProducts.keys()]))
    
    # for each pair remaining to process
    for pair in pairs.keys():
        # for each outgroup genome
        for name in outgroupProducts.keys():
            # extract the outgroup product sizes for this pair
            result = outgroupProducts[name][pair]

            # if there is only one primer size, then save it
            if len(result) == 1:
                contig, size = result.pop()
            
            # otherwise
            else:
                # remove any null products from the set
                try: result.remove(__NULL_PRODUCT)
                except KeyError: pass
                
                # if there is only one primer size, then save it
                if len(result) == 1:
                    contig, size = result.pop()
                
                # otherwise
                else:
                    # combine all contigs and pcrLens into separate lists
                    contigs = list()
                    pcrLens = list()
                    for contig,pcrLen in result:
                        contigs.append(contig)
                        pcrLens.append(pcrLen)
                    
                    # convert the contigs and lengths to comma-separated strings
                    contig = ",".join(contigs)
                    size = ",".join(map(str, pcrLens))
                
            # create and save the Product for this pair (add fake bins and fake coordinates)
            pairs[pair][name] = Product(contig, size, FAKE_BIN, FAKE_BIN, pairs[pair][ingroupName].dimerTm, FAKE_COORD, FAKE_COORD, FAKE_STRAND, FAKE_COORD, FAKE_COORD, FAKE_STRAND)


def _removeOutgroupPrimers(pairs:dict[tuple[Primer,Primer],dict[str,Product]], params:Parameters) -> None:
    """removes primers found in the outgroup that produce disallowed product sizes

    Args:
        pairs (dict[tuple[Primer,Primer],dict[str,Product]]): key=Primer pair; dict:key=genome name; val=Product
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
    
    # get kmer encodings
    encodings = __getKmerEncodings(pairs.keys())
    
    # for each outgroup genome
    for fn in params.outgroupFns:
        # get the name 
        name = os.path.basename(fn)
        
        # initialize a dictionary for the current outgroup genome
        outgroupKmers[name] = dict()

        with open(fn, 'r') as fh:
            # add each contig in the genome to the argument list
            for contig in SeqIO.parse(fh, params.format):
                # extract kmer data
                outgroupKmers[name][contig.id] = __extractKmerPositions(str(contig.seq), encodings)
        
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
                        outgroupProducts[name][(fwd,rev)].update({(contig, x) for x in products})
        
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
