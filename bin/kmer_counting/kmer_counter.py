from Bio import SeqIO
from bin.kmer_counting._kmer_counter import (count_kmers_rolling_encoding,
                                             count_allowlist_kmers_rolling_encoding,
                                             decode_kmer_encoding,
                                             gc_percentage_kmer_encoding,
                                             is_palindrome_kmer_encoding,
                                             has_gc_clamp_kmer_encoding)


def __getSingletonEncodings(counts:dict[int,int]) -> set[int]:
    """gets kmer encodings that appeared exactly once

    Args:
        counts (dict[int,int]): {kmer encoding: num appearances}

    Returns:
        set[int]: kmer encodings that appear once
    """
    return {encoding for encoding,count in counts.items() if count == 1}


def _getFilteredKmerEncodings(seq:str, k:int, minGc:float, maxGc:float) -> set[int]:
    """gets kmer encodings that pass a collection of filters

    Args:
        seq (str): the sequence for which to get kmer encodings
        k (int): the size of the kmer
        minGc (float): the minimum allowed GC (%)
        maxGc (float): the maximum allowed GC (%)

    Returns:
        set[int]: a set of kmer encodings that:
                    * appear once
                    * are not palindromic
                    * fall within (inclusive) the GC range
                    * has a GC clamp on at least one end
    """
    # helper functions
    def isNotPalindromic(enc:int) -> bool:
        return not is_palindrome_kmer_encoding(enc, k)
    
    def isGcWithinRange(enc:int) -> bool:
        gc = gc_percentage_kmer_encoding(enc, k)
        return minGc <= gc <= maxGc
    
    def hasGcClamp(enc:int) -> bool:
        return has_gc_clamp_kmer_encoding(enc, k)

    # count all the kmers in the sequence
    kmerCounts:dict[int,int] = count_kmers_rolling_encoding(seq, k)

    # only keep singleton kmers
    out = __getSingletonEncodings(kmerCounts)

    # keep non-palindromes, those within GC range, and those with a GC clamp
    return {x for x in out if isNotPalindromic(x) and isGcWithinRange(x) and hasGcClamp(x)}


def _getAllowedKmerEncodings(seq:str, k:int, allowed:set[int]) -> set[int]:
    """gets kmer encodings that appear exactly once and are contained in the allowed set

    Args:
        seq (str): the sequence for which to get kmer encodings
        k (int): the kmer size
        allowed (set[int]): a set of kmer encodings that are allowed to be retained

    Returns:
        set[int]: the kmer encodings that appear exactly once and are present in the allowed set
    """
    # count only the kmer encodings that are allowed
    kmerCounts = count_allowlist_kmers_rolling_encoding(seq, k, allowed)

    # return only those kmer ecodings that appear exactly once
    return __getSingletonEncodings(kmerCounts)


def _decodeKmerEncoding(encoding:int, k:int) -> str:
    """decodes a kmer encoding to its corresponding string

    Args:
        encoding (int): the kmer encoding
        k (int): the kmer size

    Returns:
        str: the decoded kmer
    """
    return decode_kmer_encoding(encoding, k)
