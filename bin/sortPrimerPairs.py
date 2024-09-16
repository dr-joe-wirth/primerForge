import os, numpy
from typing import Union
from bin.Primer import Primer
from bin.Product import Product
from bin.Parameters import Parameters


def __parseOutgroupProductSizes(product:Union[str,int]) -> list[int]:
    """parses the PCR product sizes for an outgroup genome

    Args:
        product (Union[str,int]): a comma-separated list of ints (str) or a single int (int)

    Returns:
        list[int]: a list of PCR product sizes
    """
    # constant
    SEP = ","
    
    # convert comma-separated list (str) to a list of ints
    if type(product) is str:
        out = [int(x) for x in product.split(SEP)]
    
    # if it just an int, then make it a list
    elif type(product) is int:
        out = [product]
    
    return out


def __minDiffFromRange(intgr:int, rng:range) -> int:
    """calculates the minimum distance of an integer from a range

    Args:
        intgr (int): an integer
        rng (range): a range

    Returns:
        int: the minimum distance of the integer from the range
    """
    # calculate the differnces from the upper and lower limits
    lower = abs(intgr - min(rng))
    upper = abs(intgr - max(rng))
    
    # return the minimum distance from the disallowed range
    return min((lower, upper))


def __sortOutgroupProducts(metadata:dict[str,Product], outgroup:list[str], ingroupPcrLens:range) -> tuple[int,int]:
    """evaluates the outgroup products for a primer pair in order to sort the pair

    Args:
        metadata (dict[str,Product]): the metadata for a primer pair
        outgroup (list[str]): a list of outgroup filenames
        ingroupPcrLens (range): the range that ingroup PCR products are within

    Returns:
        tuple[int,int]: minimum distance from the range (negated); the number of outgroup products
    """
    # initialize variables
    products = list()
    
    # determine the genome names from the ingroup file names
    outgroup = frozenset(map(os.path.basename, outgroup))
    
    # get the outgroup products that are >0bp
    products = [v.size for k,v in metadata.items() if k in outgroup and v.size != 0]
    
    # handle situations where there are no products (best case)
    if products == []:
        minDiff = float('inf') # largest number possible
    
    # otherwise, need to calculate the minimum distance from the disallowed range
    else:
        # parse the outgroup products so that we have a list ints
        products = list(map(__parseOutgroupProductSizes, products))
        products = [y for x in products for y in x]
        
        # get the smallest value distance from the disallowed range
        minDiff = min(map(lambda x:__minDiffFromRange(x, ingroupPcrLens), products))
    
    # return the minimum distance (negated so larger differences are sorted first) and the number of products
    return -minDiff, len(products)


def __getHomopolymerTemps(pair:tuple[Primer,Primer]) -> float:
    """returns the sum of the homodimer melting temps for a primer pair

    Args:
        pair (tuple[Primer,Primer]): a pair of primers

    Returns:
        float: the sum of both primer's homodimer melting temperatures
    """
    # parse the pair as strings
    fwd,rev = pair
    
    # return the sum of the dimer temps
    return fwd.homodimerTm + rev.homodimerTm


def __getHairpinTemps(pair:tuple[Primer,Primer]) -> float:
    """gets the sume of hairpin melting tempertaures for a primer pair

    Args:
        pair (tuple[Primer,Primer]): a pair of primer

    Returns:
        float: the sum of the hairpin melting tempertures
    """
    # parse the pair as strings
    fwd,rev = pair
    
    # return the sum of the hairpin temps
    return fwd.hairpinTm + rev.hairpinTm


def __getHeteroDimerTemp(products:dict[str,Product]) -> float:
    """gets the heterodimer melting temperature for a pair of primers

    Args:
        pair (tuple[Primer,Primer]): a pair of primers

    Returns:
        float: the heterodimer melting temperature
    """
    # all products have the same dimerTm; extract it from the first one
    return next(iter(products.values())).dimerTm


def __getGcDiff(pair:tuple[Primer,Primer]) -> float:
    """finds the G+C content (mol %) difference between a pair of primers

    Args:
        pair (tuple[Primer,Primer]): a pair of primers

    Returns:
        float: the difference in G+C content (mol %)
    """
    # parse the pair
    fwd,rev = pair
    
    return abs(fwd.gcPer - rev.gcPer)


def __getTmDiff(pair:tuple[Primer,Primer]) -> float:
    """finds the melting temperature difference between a pair of primers

    Args:
        pair (tuple[Primer,Primer]): a pair of primers

    Returns:
        float: the difference between melting temperatures
    """
    # parse the pair
    fwd,rev = pair

    return abs(fwd.Tm - rev.Tm)


def __sortIngroupProducts(metadata:dict[str,Product], ingroup:list[str]) -> tuple[float,float]:
    """calculates the median product size for the ingroup genomes

    Args:
        metadata (dict[str,Product]): the metadata for a primer pair

    Returns:
        tuple[float,float]: variance of product sizes; median product size (negated)
    """
    # determine the genome names from the ingroup file names
    ingroup = frozenset(map(os.path.basename, ingroup))
    
    # extract the sizes from the metadata for the pair
    products = [prod.size for name,prod in metadata.items() if name in ingroup]
    
    # calculate the variance and the median of the 
    variance = numpy.var(products)
    median = numpy.var(products)
    
    # return the variance and the median (negated to prioritize larger product sizes)
    return variance, -median


def _sortPairs(params:Parameters, pairs:dict[tuple[Primer,Primer],dict[str,Product]]) -> list[tuple[Primer,Primer]]:
    """sorts pairs of primers based on several characteristics:
          * first by difference in GC
          * next by difference in Tm
          * next by primer-dimer Tm
          * next by homopolymer Tm
          * next by hairpin Tm
          * finally by product size (largest to smallest)
    
    Args:
        params (Parameters): a Parameters object
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]): key=primer pair; val=dict; key=genome name; val=(contig, pcr size, binpair)

    Returns:
        list[tuple[Primer,Primer]]: a list of sorted primer pairs
    """
    # sort the results using the following criteria (in order):
    #  * largest difference from the ingroup product range for outgroup PCR products
    #  * fewest number of outgroup PCR products
    #  * lowest G+C difference between the pair
    #  * lowest Tm difference between the pair
    #  * lowest Tm for heterodimers
    #  * lowest Tm (sum) for homodimers
    #  * lowest Tm (sum) for hairpins
    #  * lowest variance in ingroup PCR product sizes
    #  * largest median ingroup PCR product size
    return sorted(pairs.keys(), key=lambda p: (__sortOutgroupProducts(pairs[p],
                                                                      params.outgroupFns,
                                                                      range(params.minPcr, params.maxPcr+1)),
                                               __getGcDiff(p),
                                               __getTmDiff(p),
                                               __getHeteroDimerTemp(pairs[p]),
                                               __getHomopolymerTemps(p),
                                               __getHairpinTemps(p),
                                               __sortIngroupProducts(pairs[p], params.ingroupFns)))
