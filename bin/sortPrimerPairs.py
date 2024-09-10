import os, numpy, primer3
from bin.Primer import Primer
from bin.Parameters import Parameters


def __getHomopolymerTemps(pair:tuple[Primer,Primer]) -> float:
    """returns the sum of the homodimer melting temps for a primer pair

    Args:
        pair (tuple[Primer,Primer]): a pair of primers

    Returns:
        float: the sum of both primer's homodimer melting temperatures
    """
    # parse the pair as strings
    fwd,rev = map(str, pair)
    
    # calculate dimer tms
    fDimer = primer3.calc_homodimer_tm(fwd)
    rDimer = primer3.calc_homodimer_tm(rev)
    
    # return the sum of these temps
    return fDimer + rDimer


def __getHairpinTemps(pair:tuple[Primer,Primer]) -> float:
    """gets the sume of hairpin melting tempertaures for a primer pair

    Args:
        pair (tuple[Primer,Primer]): a pair of primer

    Returns:
        float: the sum of the hairpin melting tempertures
    """
    # parse the pair as strings
    fwd,rev = map(str, pair)

    # calculate hairpin tms
    fHairpin = primer3.calc_hairpin_tm(fwd)
    rHairpin = primer3.calc_hairpin_tm(rev)
    
    # return the sum of these temps
    return fHairpin + rHairpin


def __getHeteroDimerTemp(pair:tuple[Primer,Primer]) -> float:
    """gets the heterodimer melting temperature for a pair of primers

    Args:
        pair (tuple[Primer,Primer]): a pair of primers

    Returns:
        float: the heterodimer melting temperature
    """
    
    # parse the pair as strings
    fwd,rev = map(str, pair)

    return primer3.calc_heterodimer_tm(fwd, rev)


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


def __getMedianIngroupProductSize(metadata:dict[str,tuple[str,int,tuple[str,int,int]]], ingroup:list[str]) -> float:
    """calculates the median product size for the ingroup genomes

    Args:
        metadata (dict[str,tuple[str,int,tuple[str,int,int]]]): the metadata for a primer pair

    Returns:
        float: the median PCR product size for the ingroup genomes
    """
    # determine the genome names from the ingroup file names
    ingroup = set(map(os.path.basename, ingroup))
    
    # extract the sizes from the metadata for the pair
    sizes = [size for name,(contig,size,binpair) in metadata.items() if name in ingroup]
    
    # negate the median so that pairs will be sorted largest to smallest
    return -numpy.median(sizes)


def _sortPairs(params:Parameters, pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]) -> list[tuple[Primer,Primer]]:
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
    # sort the results:
    #  * first by difference in GC
    #  * next by difference in Tm
    #  * next by primer-dimer Tm
    #  * next by homopolymer Tm
    #  * next by hairpin Tm
    #  * finally by product size (largest to smallest)
    return sorted(pairs.keys(), key=lambda p: (__getGcDiff(p),
                                               __getTmDiff(p),
                                               __getHeteroDimerTemp(p),
                                               __getHomopolymerTemps(p),
                                               __getHairpinTemps(p),
                                               __getMedianIngroupProductSize(pairs[p], params.ingroupFns)))
