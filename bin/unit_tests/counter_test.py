from Bio.Seq import Seq
from collections import defaultdict
import pathlib, unittest, string, sys

sys.path.append(str(pathlib.Path(__file__).parent.parent.parent))

from bin.getCandidateKmers import __JUNCTION_CHAR as JC
from bin.kmer_counting.kmer_counter import (_decodeKmerEncoding,
                                            _getAllowedKmerEncodings,
                                            _getFilteredKmerEncodings)
from bin.kmer_counting._kmer_counter import (count_allowlist_kmers_rolling_encoding,
                                             count_kmers_rolling_encoding,
                                             count_gc_kmer_encoding,
                                             decode_kmer_encoding,
                                             encode_base,
                                             gc_percentage_kmer_encoding,
                                             has_gc_clamp_kmer_encoding,
                                             has_long_homopolymer_in_kmer_encoding,
                                             is_palindrome_kmer_encoding)

class CounterTest(unittest.TestCase):
    """test class for evaluating kmer counting functionality
    """
    # characters that are not expected to appear in sequences
    BAD_CHARS = [x for x in string.ascii_letters + string.digits + string.punctuation if x not in f'AaCcGgTt{JC}']

    # expected key-value pairs for single-base encodings
    SINGLE_BASES = {'a': 0,
                    'A': 0,
                    'c': 1,
                    'C': 1,
                    'g': 2,
                    'G': 2,
                    't': 3,
                    'T': 3}

    # a collection of kmers to test with
    KMERS = ['ATTACGCGGAAGCCG',
             'TTAACGCGACGTCA',
             'GATCGGACGGCGCGCGAAATATATCGA',
             'GTCAACCATGGTTGAC',   # palindrome
             'cgatgctgtgtcagtagtgcaattttaat',
             'TAGGCCTAGCTAGCTATACGTACG',
             'CGTACGTATAGCTAGCTAGGCCTA', # revcomp of prev seq
             'ATCG' * 10,
             'A',
             'C',
             'G',
             'T'
             'A' * 32,
             'a' * 32,
             'C' * 32,
             'c' * 32,
             'G' * 32,
             'g' * 32,
             'T' * 32,
             't' * 32]
    
    # a sequence to test with
    SEQ = 'TAATCAATTACGCGGAAGCCGTCGAACTTAACGCGACGTCAAGATTGGATCGGACGGCGCGCGAAATATATCGAAGTCAAGTCAACCATGGTTGACTTCGGACGATGCTGTGTCAGTAGTGCAATTTTAATTGGTACTAGGCCTAGCTAGCTATACGTACGTTTCAA'

    # python implementations to test functionality
    def encodeKmer(kmer:str) -> int:
        """encodes a kmer to its integer

        Args:
            kmer (str): the kmer to encode

        Returns:
            int: the encoding
        """
        # initialize a string
        binaryString = ''

        # for each base
        for base in kmer:
            # get the binary encoding of the base
            binary = bin(encode_base(ord(base)))

            # ensure the string is 2 characters long and drop the preceding 0b
            binaryString += binary.removeprefix('0b').rjust(2, '0')
        
        # cast the string as an integer (base 2 counting)
        return int(binaryString, 2)
    
    def countKmerEncodings(seq:str, k:int) -> dict[int,int]:
        """counts kmer encodings in a given sequence

        Args:
            seq (str): the sequence to count kmers in
            k (int): the kmer size

        Returns:
            dict[int,int]: {kmer encoding: count}
        """
        # initialize the output
        out = defaultdict(int)

        # count each kmer (key is the encoding)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            out[CounterTest.encodeKmer(kmer)] += 1
        
        return out
    
    def countKmerEncodingsAllowList(seq:str, k:int, allowed:set[int]) -> dict[int,int]:
        """counts kmer encodings for those encodings in an allowed set

        Args:
            seq (str): the sequence for which to get kmer encodings
            k (int): the kmer size
            allowed (set[int]): a set of kmer encodings that are allowed

        Returns:
            dict[int,int]: {kmer encoding: count}
        """
        # count the kmers in the sequence
        kmerCounts = CounterTest.countKmerEncodings(seq, k)

        # return the counts for those encodings in the allowed set
        return {k:v for k,v in kmerCounts.items() if k in allowed}
    
    def countKmerGc(seq:str) -> int:
        """counts the number of G/C in a sequence

        Args:
            seq (str): the sequence to count

        Returns:
            int: the number G/C
        """
        return seq.upper().count('C') + seq.upper().count('G')
    
    def kmerGcPercent(seq:str) -> float:
        """calculates the GC % for a given sequence

        Args:
            seq (str): the sequence to evaluate

        Returns:
            float: the G+C %
        """
        return CounterTest.countKmerGc(seq) / len(seq) * 100.0

    def kmerHasGcClamp(seq:str) -> bool:
        """determines if a sequence has a GC clamp on at least one side

        Args:
            seq (str): the sequence to evaluate

        Returns:
            bool: indicates if the sequence has a clamp
        """
        # constants
        MIN_GC = 1
        MAX_GC = 3

        # minimum 5 bases required
        if len(seq) < 5:
            return False

        # count the number of G+C on the lefhand side and righthand size
        lhs = CounterTest.countKmerGc(seq[:5])
        rhs = CounterTest.countKmerGc(seq[-5:])

        # report if at least one side has a clamp
        return MIN_GC <= lhs <= MAX_GC or MIN_GC <= rhs <= MAX_GC

    def kmerIsPalindrome(seq:str) -> bool:
        """reports if the sequence is palindromic

        Args:
            seq (str): the sequence to evaluate

        Returns:
            bool: indicates if the sequence is a palindrome
        """
        return Seq(seq.upper()) == Seq(seq.upper()).reverse_complement()

    def getFilteredKmerEncodings(seq:str, k:int, minGc:float, maxGc:float, maxRepeat:int) -> set[int]:
        """gets a set of kmer encodings that pass a filter

        Args:
            seq (str): the sequence to evaluate
            k (int): the size of the kmer
            minGc (float): the minimum G+C (%)
            maxGc (float): the maximum G+C (%)

        Returns:
            set[int]: a set of kmer encodings that:
                        * appear once
                        * are not palindromic
                        * fall within (inclusive) the GC range
                        * has a GC clamp on at least one end
        """
        # helper function to filter kmers
        def passesFilter(s:str, c:int) -> bool:
            appearsOnce = c == 1
            isNotPalindromic = not CounterTest.kmerIsPalindrome(s)
            gcInRange = minGc <= CounterTest.kmerGcPercent(s) <= maxGc
            hasClamp = CounterTest.kmerHasGcClamp(s)
            hasNoLongRepeat = not CounterTest.hasLongHomopolymer(s, maxRepeat)

            return appearsOnce and isNotPalindromic and gcInRange and hasClamp and hasNoLongRepeat

        # count number of occurences for all the kmers
        kmers = defaultdict(int)
        for i in range(len(seq) - k + 1):
            kmers[seq[i:i+k]] += 1
        
        # keep only kmers that pass the filter; return the encodings
        return {CounterTest.encodeKmer(kmer) for kmer,count in kmers.items() if passesFilter(kmer, count)}

    def getAllowedKmerEncodings(seq:str, k:int, allowed:set[str]) -> set[int]:
        """gets all the kmer encodings that appear once and are in the allowed list

        Args:
            seq (str): the sequence for which to get kmers
            k (int): the kmer size
            allowed (set[str]): a set of allowed kmer sequences

        Returns:
            set[int]: the encodings of the kmers that appear once and are allowed
        """
        # ensure all the allowed kmer are uppercase
        allowed = {x.upper() for x in allowed}

        # count the number of times each kmer appears
        kmers = defaultdict(int)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k].upper()
            if kmer in allowed:
                kmers[kmer] += 1
        
        # return the encodings for the kmers that appear exactly once
        return {CounterTest.encodeKmer(kmer) for kmer,count in kmers.items() if count == 1}

    def hasLongHomopolymer(seq:str, maxLen:int) -> bool:
        """detects long homopolymers in a sequence

        Args:
            seq (str): the sequence to evaluate
            maxLen (int): the maximum allowed homopolymer length

        Returns:
            bool: indicates if the sequence contain a long homopolymer
        """
        # ensure upper case
        seq = seq.upper()

        # search the sequence for each homopolymer that exceeds the maximum
        for homopoly in [x * (maxLen + 1) for x in 'ATCG']:
            if homopoly in seq:
                return True
        
        return False

    # tests
    def testA_encoding(self) -> None:
        """test if the encoding functionality is working
        """
        # encoding of all invalid characters
        INVALID_VAL = 4

        # check that valid bases get encoded properly
        for base,val in CounterTest.SINGLE_BASES.items():
            self.assertEqual(encode_base(ord(base)), val)
            self.assertEqual(encode_base(ord(base)), CounterTest.encodeKmer(base))
        
        # check that the junction character is 4
        self.assertEqual(encode_base(ord(JC)), INVALID_VAL)

        # check all other invalid characters
        for char in CounterTest.BAD_CHARS:
            if char not in CounterTest.SINGLE_BASES.keys() and char != JC:
                self.assertEqual(encode_base(ord(char)), INVALID_VAL)
        
        # check that the rolling encoder is properly encoding kmers
        for kmer in CounterTest.KMERS:
            # kmers less than 32bp are allowed
            if len(kmer) <= 32:
                counts = count_kmers_rolling_encoding(kmer, len(kmer))
                encoding = next(iter(counts.keys()))
                self.assertEqual(encoding, CounterTest.encodeKmer(kmer))
    
    def testB_decoding(self) -> None:
        """test if the decoding functionality is working
        """
        # make sure all single bases work
        for base,val in CounterTest.SINGLE_BASES.items():
            self.assertEqual(decode_kmer_encoding(val, 1), base.upper())
            self.assertEqual(_decodeKmerEncoding(val, 1), base.upper())
        
        # decode the test kmers
        for kmer in CounterTest.KMERS:
            # encode the kmer
            encoding = CounterTest.encodeKmer(kmer)

            # encodings of 32bp or less should work
            if len(kmer) <= 32:
                self.assertEqual(decode_kmer_encoding(encoding, len(kmer)), kmer.upper())
                self.assertEqual(_decodeKmerEncoding(encoding, len(kmer)), kmer.upper())
            
            # encodings greater than 32bp should fail bc the integers are too large
            # NOTE: this class's encoder can generate integers that are too big but
            #       the rolling counter encoder will generate undefined 32bit ints that
            #       can be decoded...but the decoding will not match the original kmer
            if len(kmer) > 32:
                with self.assertRaises(OverflowError):
                    decode_kmer_encoding(encoding, len(kmer))
                with self.assertRaises(OverflowError):
                    _decodeKmerEncoding(encoding, len(kmer))

    def testC_countKmers(self) -> None:
        """test if the rolling kmer counter works
        """
        for k in range(1,81):
            # all counts should match for kmers of size 1-32
            if k <= 32:
                observed = count_kmers_rolling_encoding(CounterTest.SEQ, k)
                expected = CounterTest.countKmerEncodings(CounterTest.SEQ, k)
                self.assertDictEqual(observed, expected)
            
            # k greater than 32 should fail
            else:
                with self.assertRaises(OverflowError):
                    count_kmers_rolling_encoding(CounterTest.SEQ, k)

    def testD_countKmersAllowList(self) -> None:
        """tests if rolling allowlist kmer counter
        """
        for k in range(1,81):
            # only keep the first six kmers
            allowed = {CounterTest.encodeKmer(x) for x in CounterTest.KMERS[:6] if len(x) == k}

            # all counts should match for kmers of size 1-32
            if k <= 32:
                observed = count_allowlist_kmers_rolling_encoding(CounterTest.SEQ, k, allowed)
                expected = CounterTest.countKmerEncodingsAllowList(CounterTest.SEQ, k, allowed)
                self.assertDictEqual(observed, expected)
            
            # k greater than 32 should fail
            else:
                with self.assertRaises(OverflowError):
                    count_allowlist_kmers_rolling_encoding(CounterTest.SEQ, k, allowed)

    def testE_countKmerGc(self) -> None:
        """tests the G+C counting functionality
        """
        for kmer in CounterTest.KMERS:
            # encode the kmer and count the G+C
            encoding = CounterTest.encodeKmer(kmer)
            expected = CounterTest.countKmerGc(kmer)

            # encodings of 32bp or less should work
            if len(kmer) <= 32:
                observed = count_gc_kmer_encoding(encoding, len(kmer))
                self.assertEqual(observed, expected)
            
            # encodings greater than 32bp should fail bc the integers are too large
            # NOTE: this class's encoder can generate integers that are too big but
            #       the rolling counter encoder will generate undefined 32bit ints that
            #       can be decoded...but the decoding will not match the original kmer
            else:
                with self.assertRaises(OverflowError):
                    count_gc_kmer_encoding(encoding, len(kmer))

    def testF_kmerGcPercentage(self) -> None:
        """tests the G+C percentage calculator
        """
        for kmer in CounterTest.KMERS:
            # encode the kmer and get the G+C %
            encoding = CounterTest.encodeKmer(kmer)
            expected = CounterTest.kmerGcPercent(kmer)

            # encodings of 32bp or less should work
            # equal to 10 decimal points is good enough
            if len(kmer) <= 32:
                observed = gc_percentage_kmer_encoding(encoding, len(kmer))
                self.assertAlmostEqual(observed, expected, 10)
            
            # encodings greater than 32bp should fail bc the integers are too large
            # NOTE: this class's encoder can generate integers that are too big but
            #       the rolling counter encoder will generate undefined 32bit ints that
            #       can be decoded...but the decoding will not match the original kmer
            else:
                with self.assertRaises(OverflowError):
                    gc_percentage_kmer_encoding(encoding, len(kmer))

    def testG_kmerGcClamp(self) -> None:
        """tests GC clamp detection
        """
        for kmer in CounterTest.KMERS:
            # encode the kmer and detect a clamp
            encoding = CounterTest.encodeKmer(kmer)
            expected = CounterTest.kmerHasGcClamp(kmer)

            # encodings of 32bp or less should work
            # equal to 10 decimal points is good enough
            if len(kmer) <= 32:
                observed = has_gc_clamp_kmer_encoding(encoding, len(kmer))
                self.assertAlmostEqual(observed, expected, 10)
            
            # encodings greater than 32bp should fail bc the integers are too large
            # NOTE: this class's encoder can generate integers that are too big but
            #       the rolling counter encoder will generate undefined 32bit ints that
            #       can be decoded...but the decoding will not match the original kmer
            else:
                with self.assertRaises(OverflowError):
                    has_gc_clamp_kmer_encoding(encoding, len(kmer))

    def testH_kmerPalindrome(self) -> None:
        """tests palindrome detection
        """
        for kmer in CounterTest.KMERS:
            # encode the kmer and detect a palindrome
            encoding = CounterTest.encodeKmer(kmer)
            expected = CounterTest.kmerIsPalindrome(kmer)

            # encodings of 32bp or less should work
            if len(kmer) <= 32:
                observed = is_palindrome_kmer_encoding(encoding, len(kmer))
                self.assertEqual(observed, expected)
            
            # encodings greater than 32bp should fail bc the integers are too large
            # NOTE: this class's encoder can generate integers that are too big but
            #       the rolling counter encoder will generate undefined 32bit ints that
            #       can be decoded...but the decoding will not match the original kmer
            else:
                with self.assertRaises(OverflowError):
                    is_palindrome_kmer_encoding(encoding, len(kmer))

    def testI_filteredKmers(self) -> None:
        """tests filtered kmer retrieval
        """
        # set the G+C percentage range
        MIN_GC = 40.0
        MAX_GC = 60.0
        MAX_REPEAT = 3

        for k in range(1,81):
            # kmers of 32bp or less should work
            if k <= 32:
                if k == 5:
                    observed = _getFilteredKmerEncodings(CounterTest.SEQ, k, MIN_GC, MAX_GC, MAX_REPEAT)
                    expected = CounterTest.getFilteredKmerEncodings(CounterTest.SEQ, k, MIN_GC, MAX_GC, MAX_REPEAT)

                    self.assertSetEqual(observed, expected)
            
            # kmers more than 32bp should fail
            else:
                with self.assertRaises(OverflowError):
                    _getFilteredKmerEncodings(CounterTest.SEQ, k, MIN_GC, MAX_GC, MAX_REPEAT)
        
    def testJ_allowedKmers(self) -> None:
        """tests kmer allowlist retrieval
        """
        for k in range(1,81):
            # only keep the first six kmers
            allowedSeqs = {x for x in CounterTest.KMERS[:6] if len(x) == k}
            allowedEncodings = {CounterTest.encodeKmer(x) for x in allowedSeqs}

            # any k less than or equal to 32 should work
            if k <= 32:
                observed = _getAllowedKmerEncodings(CounterTest.SEQ, k, allowedEncodings)
                expected = CounterTest.getAllowedKmerEncodings(CounterTest.SEQ, k, allowedSeqs)

                self.assertSetEqual(observed, expected)
            
            # k larger than 32 should throw an overflow error
            else:
                with self.assertRaises(OverflowError):
                    _getAllowedKmerEncodings(CounterTest.SEQ, k, allowedEncodings)
    
    def testK_noLongHomopolymers(self) -> None:
        """tests long homopolymer detection
        """
        # for lengths 2 through 31
        for length in range(2,31):
            for kmer in CounterTest.KMERS:
                # encode each kmer
                encoding = CounterTest.encodeKmer(kmer)

                # should work for 32bp or
                if len(kmer) <= 32:
                    observed = has_long_homopolymer_in_kmer_encoding(encoding, len(kmer), length)
                    expected = CounterTest.hasLongHomopolymer(kmer, length)

                    self.assertEqual(observed, expected)
                
                # encodings greater than 32bp should fail bc the integers are too large
                # NOTE: this class's encoder can generate integers that are too big but
                #       the rolling counter encoder will generate undefined 32bit ints that
                #       can be decoded without raising an error...but the decoding will not
                #       match the original kmer
                else:
                    with self.assertRaises(OverflowError):
                        has_long_homopolymer_in_kmer_encoding(encoding, len(kmer), length)


if __name__ == "__main__":
    unittest.main()
