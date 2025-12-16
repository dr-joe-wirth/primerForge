# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False

import cython
from libc.stdint cimport uint64_t, uint32_t
from libc.stdlib cimport malloc, free


# boundscheck(False) and wraparound(False) are safe here because we perform no array indexing
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef inline int encode_base(char base) nogil:
    """Convert nucleotide base to 2-bit integer.
    
    :param base: a nucleotide base to get a 2-bit encoding for
    
    :return: an int encoding of the base
    """
    if base == 'A' or base == 'a':
        return 0
    elif base == 'C' or base == 'c':
        return 1
    elif base == 'G' or base == 'g':
        return 2
    elif base == 'T' or base == 't':
        return 3
    else:
        return 4  # Invalid base (N, etc.)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef inline uint64_t reverse_complement_kmer_encoding(uint64_t kmer_encoding, int k) nogil:
    """
    Reverse-complement a 2-bit encoded kmer.

    Encoding: A=0, C=1, G=2, T=3
    NOTE: k must be <= 32
    """
    cdef uint64_t revcomp = 0
    cdef int i
    cdef uint64_t base

    for i in range(k):
        # extract lowest 2 bits
        base = kmer_encoding & 3

        # complement (flip both bits)
        base = (~base) & 3

        # shift into result
        revcomp = (revcomp << 2) | base

        # shift input down
        kmer_encoding >>= 2

    return revcomp


# boundscheck(False) is safe here because looping over i in range(seq_len), so indices are guaranteed to be in-bounds
# wraparound(False) is safe here because no negative indexing is utilized
@cython.boundscheck(False)
@cython.wraparound(False)
def count_kmers_rolling_encoding(str sequence, int k) -> dict[int,int]:
    """Fast kmer counting using rolling encoding with 2-bit encoding.
    NOTE: kmers larger than k=32 will exhibit undefined behavior.

    :param sequence: a nucleotide sequence to count kmers for
    :param k: the size of kmers to count

    :return: a dictionary: {kmer_hash: count}
    """
    # enfore maximum kmer sizes of 32
    if k > 32:
        raise OverflowError(f'kmer size {k} is too large to encode (maximum 32)')
    
    # For masking out high bits irrelevant to kmers smaller than 32 bits
    # `mask` starts as all 1s, then we shift right until only kmer-relevant bits are 1s
    cdef uint64_t mask = (~0ULL) >> (64 - (2 * k))
    cdef int seq_len = len(sequence)
    cdef uint64_t kmer_encoding = 0
    cdef int valid_bases = 0
    cdef uint32_t base_val
    cdef int i
    cdef bytes seq_bytes = sequence.encode('utf-8')
    cdef char* seq_ptr = seq_bytes
    
    # Dictionary to store {hash: count} mappings
    kmer_counts = {}
    
    # for each position in the sequence
    for i in range(seq_len):
        # cast the base as an integer
        base_val = encode_base(seq_ptr[i])
        
        # Invalid base resets the hash and base count
        if base_val == 4:
            kmer_encoding = 0
            valid_bases = 0
            continue
        
        # SHIFT the existing kmer_hash leftwards 2 bits
        # OR with base_value to load the current base encoding into the bottom bits
        # AND with mask to blank out the high bits not relevant to this kmer size
        kmer_encoding = ((kmer_encoding << 2) | base_val) & mask

        # increment the number of consecutive valid bases
        valid_bases += 1

        # If we've encountered enough sequential valid bases that this kmer_hash corresponds to a valid kmer, count it
        if valid_bases >= k:
            if kmer_encoding in kmer_counts:
                kmer_counts[kmer_encoding] += 1
            else:
                kmer_counts[kmer_encoding] = 1
    
    return kmer_counts


# boundscheck(False) is safe here because looping over i in range(seq_len), so indices are guaranteed to be in-bounds
# wraparound(False) is safe here because no negative indexing is utilized
@cython.boundscheck(False)
@cython.wraparound(False)
def count_allowlist_kmers_rolling_encoding(str sequence, int k, set allowed_encoded_kmers) -> dict[int,int]:
    """Fast k-mer counting with allowlist filtering using rolling hash.

    :param sequence: a nucleotide sequence to count kmers for
    :param k: the size of kmers to count (must be <=32)
    :param allowed_encoded_kmers: set of k-mer encodings to count (ignore others)

    :return: a dictionary: {kmer_encoding: count}
    """
    # enfore maximum kmer sizes of 32
    if k > 32:
        raise OverflowError(f'kmer size {k} is too large to encode (maximum 32)')
    
    # For masking out high bits irrelevant to kmers smaller than 32 bits
    # `mask` starts as all 1s, then we shift right until only kmer-relevant bits are 1s
    cdef uint64_t mask = (~0ULL) >> (64 - (2 * k))
    cdef int seq_len = len(sequence)
    cdef uint64_t kmer_encoding = 0
    cdef int valid_bases = 0
    cdef uint32_t base_val
    cdef int i
    cdef bytes seq_bytes = sequence.encode('utf-8')
    cdef char* seq_ptr = seq_bytes

    # Dictionary to store {hash: count} mappings
    kmer_counts = {}

    # for each position in the sequence
    for i in range(seq_len):
        # cast the base as an integer
        base_val = encode_base(seq_ptr[i])

        # Invalid base resets the hash and base count
        if base_val == 4:
            kmer_encoding = 0
            valid_bases = 0
            continue
        
        # SHIFT the existing kmer_hash leftwards 2 bits
        # OR with base_value to load the current base encoding into the bottom bits
        # AND with mask to blank out the high bits not relevant to this kmer size
        kmer_encoding = ((kmer_encoding << 2) | base_val) & mask

        # increment the number of consecutive valid bases
        valid_bases += 1

        # If we've encountered enough sequential valid bases that this kmer_hash corresponds to a valid kmer,
        # AND the kmer is one we've been told to care about, count it
        if valid_bases >= k and kmer_encoding in allowed_encoded_kmers:
            if kmer_encoding in kmer_counts:
                kmer_counts[kmer_encoding] += 1
            else:
                kmer_counts[kmer_encoding] = 1
    
    return kmer_counts


# boundscheck(False) is safe here because looping over i in range(k), so indices are guaranteed to be in-bounds
# wraparound(False) is safe here because no negative indexing is utilized
@cython.boundscheck(False)
@cython.wraparound(False)
def decode_kmer_encoding(uint64_t kmer_encoding, int k) -> str:
    """Decode a kmer hash back to nucleotide string.

    :param kmer_encoding: a 2-bit encoded kmer 
    :param k: the size of the encoded kmer

    :return: the string representation of the kmer hash
    """
    cdef char* bases = "ACGT"
    cdef char* result = <char*>malloc((k + 1) * sizeof(char))       # C strings terminate with null byte --> extra char
    cdef int i

    # if k > 32:
    #     raise OverflowError(f'{k} exceeds maximum kmer size (32)')

    # for each position in the kmer
    for i in range(k):
        # SHIFT kmer_encoding
        ## make a new int
        ### slide current base hash to end
        ### look up its value in the string array (`bases`)
        result[i] = bases[(kmer_encoding >> (2 * (k - 1 - i))) & 3] # Isolate the 2 bits for i-th base and decode it
    
    result[k] = 0                                                   # Null terminator
    py_result = result[:k].decode('ascii')                          # Cast to python string

    # no memory leaks allowed
    free(result)

    return py_result


# boundscheck(False) is safe here because looping over i in range(k), so indices are guaranteed to be in-bounds
# wraparound(False) is safe here because no negative indexing is utilized
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef inline int count_gc_kmer_encoding(uint64_t kmer_encoding, int k) nogil:
    """Calculate a kmer's GC percentage from its hash value.

    :param kemr_encoding: a 2-bit encoded kmer to calculate GC content for

    :return: the number of G or C in the encoded kmer
    """
    cdef int gc_count = 0
    cdef int i
    cdef uint32_t base_val

    # for each position in the kmer
    for i in range(k):
        # Isolate i-th kmer bases from right to left (order doesn't matter for counting G's and C's
        base_val = (kmer_encoding >> (2 * i)) & 3       # next line equivalent but processes bases left-to-right instead

        # found a C or G -- count it
        if base_val == 1 or base_val == 2:
            gc_count += 1
    
    return gc_count


# boundscheck(False) is safe here because looping over i in range(k), so indices are guaranteed to be in-bounds
# wraparound(False) is safe here because no negative indexing is utilized
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef inline double gc_percentage_kmer_encoding(uint64_t kmer_encoding, int k) nogil:
    """Calculate a kmer's GC percentage from its hash value.

    :param kmer_encoding: a 2-bit encoded kmer
    :param k: the size of the encoded kmer

    :return: gc_content as a percentage
    """
    cdef int gc_count
    
    gc_count = count_gc_kmer_encoding(kmer_encoding, k)

    return <double>(100.0 * gc_count) / k


# boundscheck(False) is safe here because looping over i in range(k // 2), so indices are guaranteed to be in-bounds
# wraparound(False) is safe here because no negative indexing is utilized
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef inline bint is_palindrome_kmer_encoding(uint64_t kmer_encoding, int k) nogil:
    """Determine whether a k-mer is a palindrome from its hash value.
    :param kmer_encoding: a 2-bit encoded kmer
    :param k: the size of the encoded kmer

    :return: a boolean indicating if the kmer is palindromic
    """
    cdef int i
    cdef uint32_t left_base, right_base

    # k=1 cannot be a palindrome
    if k % 2 == 1:
        return False

    # go through first half of the kmer
    for i in range(k // 2):
        # Isolate the left base and right base
        left_base = (kmer_encoding >> (2 * (k - 1 - i))) & 3
        right_base = (kmer_encoding >> (2 * i)) & 3

        # check if left_base == complement(right_base)
        if left_base != (3 - right_base):
            return False
    
    return True


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef inline bint has_gc_clamp_kmer_encoding(uint64_t kmer_encoding, int k) nogil:
    """Determine whether a kmer has a GC clamp at either end from its hash value.
    
    Checks if either the leftmost 5 bases or rightmost 5 bases have 1-3 GC
    (good for PCR primer design).
    
    :param kmer_encoding: a 2-bit encoded kmer hash
    :param k: the size of the encoded kmer (must be >= 5)

    :return: True if either end has appropriate GC content
    """
    # constants
    cdef int MIN_GC = 1
    cdef int MAX_GC = 3

    cdef double left_gc, right_gc
    cdef uint64_t left_end, right_end

    # Sanity check -- if k < 5 we can't possibly have a GC clamp
    if k < 5:
        return False
    
    # Extract leftmost 5 bases (10 bits) of the kmer
    left_end = kmer_encoding >> (2 * (k - 5))
    left_gc = count_gc_kmer_encoding(left_end, k)
    
    # Extract rightmost 5 bases (10 bits) of the kmer
    right_end = kmer_encoding & 0x3FF  # 0x3FF = 1023 = 10 bits of 1s
    right_gc = count_gc_kmer_encoding(right_end, k)
    
    # Check if either end has the allowed number of G and C
    return (MIN_GC <= left_gc <= MAX_GC) or (MIN_GC <= right_gc <= MAX_GC)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef has_long_homopolymer_in_kmer_encoding(uint64_t kmer_encoding, int k, int min_len):
    """Determine whether a kmer has a long repeating homopolymer
    
    :param kmer_encoding: a 2-bit encoded kmer
    :param k: the size of the encoded kmer (must be >= 5)
    :param max_len: the minimum length of a homopolymer to be detected

    :return: True if either end has appropriate GC content
    """
    cdef int count = 0
    cdef int last_seen = 4  # have to start with an known invalid base

    # for each position in the kmer
    for i in range(k):
        # Isolate i-th kmer bases from right to left
        base_val = (kmer_encoding >> (2 * i)) & 3       # next line equivalent but processes bases left-to-right instead

        if base_val == last_seen:
            count += 1
        else:
            last_seen = base_val
            count = 1
        
        if count > min_len:
            return True

    return False


# boundscheck(False) is safe here because looping over i in range(seq_len), so indices are guaranteed to be in-bounds
# wraparound(False) is safe here because no negative indexing is utilized
@cython.boundscheck(False)
@cython.wraparound(False)
def first_start_position_from_kmer_encodings(str sequence, int k, dict kmer_start_positions):
    """Fast k-mer counting with allowlist filtering using rolling hash.
    Only finds the first start position. This is ok if the kmer appears exactly once

    :param sequence: a nucleotide sequence to count kmers for
    :param k: the size of kmers to count (must be <=32)
    :param kmer_start_positions: a dictionary of {encoding: None} for allowed kmers

    :return: {kmer_str: index} # index is the first seen start position (strand is indicated by positive or negative start positions)
    """
    # For masking out high bits irrelevant to kmers smaller than 32 bits
    # `mask` starts as all 1s, then we shift right until only kmer-relevant bits are 1s
    cdef uint64_t mask = (~0ULL) >> (64 - (2 * k))
    cdef int seq_len = len(sequence)
    cdef uint64_t kmer_encoding = 0
    cdef int valid_bases = 0
    cdef uint32_t base_val
    cdef int end
    cdef int start
    cdef bytes seq_bytes = sequence.encode('utf-8')
    cdef char* seq_ptr = seq_bytes

    # go through the sequence
    for end in range(seq_len):
        # encode the current base
        base_val = encode_base(seq_ptr[end])

        # Invalid base resets the encoding and base count
        if base_val == 4:
            kmer_encoding = 0
            valid_bases = 0
            continue
        
        # SHIFT the existing encoding leftwards 2 bits
        # OR with base_value to load the current base encoding into the bottom bits
        # AND with mask to blank out the high bits not relevant to this kmer size
        kmer_encoding = ((kmer_encoding << 2) | base_val) & mask

        valid_bases += 1

        # If we've encountered enough sequential valid bases that this encoding corresponds to a valid kmer
        # AND the kmer is one we've been told to care about, count it
        if valid_bases >= k:
            # determine the reverse complement of the kmer
            rev_comp_encoding = reverse_complement_kmer_encoding(kmer_encoding, k)
            start = end - k + 1

            # save the current encoding and its start position if found on the plus strand
            if kmer_encoding in kmer_start_positions.keys():
                desired_encoding = kmer_encoding

            
            # calculate the start position for minus strand
            elif rev_comp_encoding in kmer_start_positions.keys():
                start = -start
                desired_encoding = rev_comp_encoding
            
            else:
                continue
            
            # remove the encoding from the dictionary (it will be replaced with its decoded string)
            del kmer_start_positions[desired_encoding]

            # save the decoded kmer and its start position
            kmer_start_positions[decode_kmer_encoding(desired_encoding, k)] = start

    return kmer_start_positions


# boundscheck(False) is safe here because looping over i in range(seq_len), so indices are guaranteed to be in-bounds
# wraparound(False) is safe here because no negative indexing is utilized
@cython.boundscheck(False)
@cython.wraparound(False)
def all_start_positions_from_kmer_encodings(str sequence, int k, set allowed_kmer_encodings):
    """Fast k-mer counting with allowlist filtering using rolling hash.
    Only finds the last start position. This is ok if the kmer appears exactly once

    :param sequence: a nucleotide sequence to count kmers for
    :param k: the size of kmers to count (must be <=32)
    :param allowed_kmer_encodings: set of k-mer hashes to count (ignore others)

    :return: {kmer_str: [indices]} # indices are the start positions (strand is indicated by positive or negative start positions)
    """
    # For masking out high bits irrelevant to kmers smaller than 32 bits
    # `mask` starts as all 1s, then we shift right until only kmer-relevant bits are 1s
    cdef uint64_t mask = (~0ULL) >> (64 - (2 * k))
    cdef int seq_len = len(sequence)
    cdef uint64_t kmer_encoding = 0
    cdef int valid_bases = 0
    cdef uint32_t base_val
    cdef int end
    cdef int start
    cdef bytes seq_bytes = sequence.encode('utf-8')
    cdef char* seq_ptr = seq_bytes
    cdef dict kmer_start_positions = dict()

    # go through the sequence
    for end in range(seq_len):
        # encode the current base
        base_val = encode_base(seq_ptr[end])

        # Invalid base resets the encoding and base count
        if base_val == 4:
            kmer_encoding = 0
            valid_bases = 0
            continue
        
        # SHIFT the existing encoding leftwards 2 bits
        # OR with base_value to load the current base encoding into the bottom bits
        # AND with mask to blank out the high bits not relevant to this kmer size
        kmer_encoding = ((kmer_encoding << 2) | base_val) & mask

        valid_bases += 1

        # If we've encountered enough sequential valid bases that this encoding corresponds to a valid kmer
        # AND the kmer is one we've been told to care about, count it
        if valid_bases >= k:
            # determine the reverse complement of the kmer
            rev_comp_encoding = reverse_complement_kmer_encoding(kmer_encoding, k)
            start = end - k + 1

            # save the current encoding and its start position if found on the plus strand
            if kmer_encoding in allowed_kmer_encodings:
                desired_encoding = kmer_encoding

            
            # calculate the start position for minus strand
            elif rev_comp_encoding in allowed_kmer_encodings:
                start = -start
                desired_encoding = rev_comp_encoding
            
            else:
                continue

            # decode the kmer
            kmer = decode_kmer_encoding(desired_encoding, k)
            
            # initialize a list if this is the first time this kmer has been seen
            if kmer not in kmer_start_positions.keys():
                kmer_start_positions[kmer] = list()

            # save the kmer and its start position
            kmer_start_positions[kmer].append(start)

    return kmer_start_positions
