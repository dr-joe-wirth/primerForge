import multiprocessing
import os
import primer3
from collections import defaultdict
from typing import Iterator, Union

from Bio import SeqIO

from bin.Clock import Clock
from bin.Genome import Genome
from bin.Parameters import Parameters
from bin.Primer import Primer
from bin.kmer_counting.kmer_counter import (
    _getAllKmerEncodings,
    _getAllowedKmerEncodings,
    _getFirstStartPositionsAndDecodeAllowedEncodings,
    _getFilteredKmerEncodings,
)

# constant
__JUNCTION_CHAR = "~"


# functions
def __getConcatenatedSequence(genome: Genome, strand: str) -> str:
    """gets the entire sequence from a file, concatenated with __JUNCTION_CHAR

    Args:
        genome (Genome): the genome sequence
        strand (str): the strand to retrieve

    Raises:
        ValueError: invalid strand specified

    Returns:
        str: the concatenated sequence
    """
    # import the forward sequences if requested
    if strand == Primer.PLUS:
        seqs = [str(r.seq) for r in genome]

    # import the reverse sequences if requested
    elif strand == Primer.MINUS:
        seqs = [str(r.seq.reverse_complement()) for r in genome]

    # fail on invalid strand
    else:
        raise ValueError(f"invalid strand specified: {strand}")

    # concatenate sequences
    return __JUNCTION_CHAR.join(seqs)


def __updateAllowedKmerEncodings(
    genome:Genome, k: int, sharedEncodings: set[int]
) -> None:
    """updates the shared encodings set with those that are found in the provided sequence

    Args:
        genome (Genome): the genome sequence
        k (int): the kmer length
        sharedEncodings (set[int]): a collection of allowed kmer encodings
    """
    # get the forward and reverse sequences
    fwd = __getConcatenatedSequence(genome, Primer.PLUS)
    rev = __getConcatenatedSequence(genome, Primer.MINUS)

    # combine sequences into a single string
    seq = __JUNCTION_CHAR.join((fwd, rev))

    # get the kmer encodings that appear in this genome
    newEncodings = _getAllowedKmerEncodings(seq, k, sharedEncodings)

    # update the shared encodings to only keep those that are truly shared
    sharedEncodings.intersection_update(newEncodings)


def __convertKmerEncodingsToPrimers(
    genomes: list[Genome], encodings: set[int], k: int
) -> dict[str, list[Primer]]:
    """converts kmer encodings to a list of Primer objects

    Args:
        genomes (list[Genome]): a list of genome sequences
        encodings (set[int]): a collection of kmer encodings
        k (int): the kmer length

    Returns:
        dict[str,list[Primer]]: a list of Primer objects
    """
    # initialize the output
    out = defaultdict(list)

    # for each file
    for genome in genomes:
        # for each contig
        for rec in genome:
            # get the start positions and decode the encoding
            positions = _getFirstStartPositionsAndDecodeAllowedEncodings(
                encodings, k, str(rec.seq)
            )

            # for each kmer and its start position
            for kmer, start in positions.items():
                # positive coordinate indicates plus strand
                if start > 0:
                    strand = Primer.PLUS

                # negative coordinate indicates minus strand
                elif start < 0:
                    start = abs(start)
                    strand = Primer.MINUS

                # if the start is 0, then does it match the beginning of the sequence?
                elif kmer == rec.seq[: len(kmer)]:
                    strand = Primer.PLUS

                # if not, then it must be minus strand
                else:
                    strand = Primer.MINUS

                # create a Primer object and save it in the list
                out[genome.name].append(Primer(kmer, rec.id, start, k, strand))

    return dict(out)


def __getSharedPrimersOneK(
    ingroup: list[Genome],
    k: int,
    minGc: float,
    maxGc: float,
    maxRepeatLen: int,
) -> dict[str, list[Primer]]:
    """gets shared primers for a single kmer length

    Args:
        ingroup (list[Genome]): the genome sequences
        k (int): the kmer length
        minGc (float): the minimum allowed GC percent
        maxGc (float): the maximum allowed GC percent
        maxRepeatLen (int): the maximum allowed repeat length

    Returns:
        list[Primer]: a list of shared Primer objects
    """
    # get the first genome's forward and reverse sequences
    fwd = __getConcatenatedSequence(ingroup[0], Primer.PLUS)
    rev = __getConcatenatedSequence(ingroup[0], Primer.MINUS)

    # get the kmer encodings for the first (smallest) genome's plus strand
    sharedEncodings = _getFilteredKmerEncodings(fwd, k, minGc, maxGc, maxRepeatLen)

    # get ALL the kmer encodings for first genome's minus strand
    minusEncodings = _getAllKmerEncodings(rev, k)

    # remove any minus strand encodings that are present
    sharedEncodings.difference_update(minusEncodings)

    # discard unused items
    del minusEncodings
    del fwd
    del rev

    # for each remaining genome, update the shared kmer encodings
    for genome in ingroup[1:]:
        __updateAllowedKmerEncodings(genome, k, sharedEncodings)

    # convert the shared kmers to Primer objects
    return __convertKmerEncodingsToPrimers(ingroup, sharedEncodings, k)


def __getSharedPrimers(params: Parameters) -> dict[str, list[Primer]]:
    """retrieves all the primers that are shared between the input genomes

    Args:
        params (Parameters): a Parameters object

    Returns:
        dict[str,list[Primer]]: {genome name: [shared Primers]}
    """

    # helper function to generate arguments for __getSharedKmersOneK
    def genArgs() -> Iterator[tuple[list[str], str, int, float, float]]:
        for k in range(params.minLen, params.maxLen + 1):
            yield params.ingroup, k, params.minGc, params.maxGc, params.maxRepeatLen

    # initialize output
    out = defaultdict(list)

    # get the shared primers for each kmer length in parallel
    with multiprocessing.Pool(params.numThreads) as pool:
        for result in pool.starmap(__getSharedPrimersOneK, genArgs()):
            for name, primers in result.items():
                out[name].extend(primers)

    # sort the lists of primers alphabetically
    for name in out.keys():
        out[name].sort()

    return dict(out)


def __evaluateOnePrimerSequence(
    args: tuple[
        Primer, int, float, float, float, float, float, float, float, int, float
    ],
) -> tuple[Union[Primer, None], int]:
    """evaluates one primer; designed for parallel calls

    Args:
        args (tuple):
            primer (Primer): the Primer to evaluate
            idx (int): the index of the primer in the list of the calling function
            minTm (float): the minimum melting temperature allowed
            maxTm (float): the maximum melting temperature allowed
            mvConc (float): primer3 mv_conc
            dvConc (float): primer3 dv_conc
            dntpConc (float): primer3 dntp_conc
            dnaConc (float): primer3 dna_conc
            tempC (float): primer3 temp_c
            maxLoop (int): primer3 max_loop
            tempTolerance (float): the minimum degrees below primer Tm allowed for secondary structures

    Returns:
        tuple[Union[Primer,None],int]: a Primer object (or None if the eval failed) and its index
    """

    # define helper functions to make booleans below more readable
    def isTmWithinRange(p: Primer) -> bool:
        """is the Tm within the acceptable range?"""
        return p.Tm >= minTm and p.Tm <= maxTm

    def noHairpins(p: Primer) -> bool:
        """verifies that the primer does not form hairpins"""
        # calculate the hairpin Tms
        p.hairpinTm = primer3.calc_hairpin_tm(
            str(p),
            mv_conc=mvConc,
            dv_conc=dvConc,
            dntp_conc=dntpConc,
            dna_conc=dnaConc,
            temp_c=tempC,
            max_loop=maxLoop,
        )
        p.rcHairpin = primer3.calc_hairpin_tm(
            str(p.reverseComplement()),
            mv_conc=mvConc,
            dv_conc=dvConc,
            dntp_conc=dntpConc,
            dna_conc=dnaConc,
            temp_c=tempC,
            max_loop=maxLoop,
        )

        # hairpin tm should be less than (minTm - tolerance°); need to check both strands
        fwdOk = p.hairpinTm < (minTm - tempTolerance)
        revOk = p.rcHairpin < (minTm - tempTolerance)

        return fwdOk and revOk

    def noHomodimers(p: Primer) -> bool:
        """verifies that the primer does not form homodimers"""
        # calculate the homodimer Tms
        p.homodimerTm = primer3.calc_homodimer_tm(
            str(p),
            mv_conc=mvConc,
            dv_conc=dvConc,
            dntp_conc=dntpConc,
            dna_conc=dnaConc,
            temp_c=tempC,
            max_loop=maxLoop,
        )
        p.rcHomodimer = primer3.calc_homodimer_tm(
            str(p.reverseComplement()),
            mv_conc=mvConc,
            dv_conc=dvConc,
            dntp_conc=dntpConc,
            dna_conc=dnaConc,
            temp_c=tempC,
            max_loop=maxLoop,
        )

        # homodimer tm should be less than (minTm - tolerance°); need to check both strands
        fwdOk = p.homodimerTm < (minTm - tempTolerance)
        revOk = p.rcHomodimer < (minTm - tempTolerance)

        return fwdOk and revOk

    # parse the arguments
    (
        primer,
        idx,
        minTm,
        maxTm,
        mvConc,
        dvConc,
        dntpConc,
        dnaConc,
        tempC,
        maxLoop,
        tempTolerance,
    ) = args

    # evaluate the primer's percent GC, Tm, hairpin potential, and homodimer potential
    if isTmWithinRange(primer):
        if noHairpins(primer):
            if noHomodimers(primer):
                # don't return this index if the primer is good
                return primer, idx

    return None, idx


def __removeBadPrimers(primers: dict[str, list[Primer]], params: Parameters) -> None:
    """evaluates primers and removes those that are not suitable; updates hairpin and homodimer Tm

    Args:
        primers (list[Primer]): the list produced by __getSharedPrimers; sorted alphabetically
        params (Parameters): a Parameters object
    """

    # generator function for getting arguments
    def generateArgs(
        n: str,
    ) -> Iterator[
        tuple[Primer, int, float, float, float, float, float, float, float, int, float]
    ]:
        """generates arguments for __evaluateOnePrimerSequence"""
        # emit arguments for each primer for the first genome only (primer lists are equivalent)
        for idx in range(len(primers[n])):
            yield (
                primers[n][idx],
                idx,
                params.minTm,
                params.maxTm,
                params.p3_mvConc,
                params.p3_dvConc,
                params.p3_dntpConc,
                params.p3_dnaConc,
                params.p3_tempC,
                params.p3_maxLoop,
                params.tempTolerance,
            )

    # get the names
    names = list(primers.keys())

    # initialize a list of bad indices
    badIndices = list()

    # parallelize primer evaluations
    with multiprocessing.Pool(processes=params.numThreads) as pool:
        # only need to evaluate the first genome; primer lists are equivalent
        # process results as they become available
        for primer, index in pool.imap_unordered(
            __evaluateOnePrimerSequence, generateArgs(names[0])
        ):
            # track which indices failed the evaluation
            if primer is None:
                badIndices.append(index)

            # replace the existing primer with one that has hairpin and homodimer Tm
            else:
                primers[names[0]][index] = primer

    # remove bad primers from the list
    badIndices.sort(reverse=True)
    for name in names:
        for idx in badIndices:
            primers[name].pop(idx)

    # update the hairpin and homodimer melting temps for the other genomes
    for name in names[1:]:
        for idx in range(len(primers[names[0]])):
            primers[name][idx].hairpinTm = primers[names[0]][idx].hairpinTm
            primers[name][idx].rcHairpin = primers[names[0]][idx].rcHairpin
            primers[name][idx].homodimerTm = primers[names[0]][idx].homodimerTm
            primers[name][idx].rcHomodimer = primers[names[0]][idx].rcHomodimer


def __buildOutput(
    primers: dict[str, list[Primer]],
) -> dict[str, dict[str, list[Primer]]]:
    """builds the datastructure needed for downstream applications

    Args:
        primers (dict[str,list[Primer]]): the dictionary produced by __getSharedPrimers

    Returns:
        dict[str,dict[str,list[Primer]]]: {name: {contig: [Primer objects]}}
    """
    # initialize the output
    out = dict()

    # for each name
    for name in primers.keys():
        # create the subordinate dictionary
        out[name] = defaultdict(list)

        # store each primer under its respective contig
        for primer in primers[name]:
            out[name][primer.contig].append(primer)

        # recast defaultdict to dict
        out[name] = dict(out[name])

    return out


def _getAllCandidateKmers(
    params: Parameters, sharedExists: bool
) -> dict[str, dict[str, list[Primer]]]:
    """gets all the candidate kmer sequences for a given ingroup

    Args:
        params (Parameters): a Parameters object
        sharedExists (bool): indicates if the shared kmers are already pickled

    Raises:
        RuntimeError: no shared ingroup kmers
        RuntimeError: shared ingroup kmers are not suitable for use as primers

    Returns:
        dict[str,dict[str,list[Primer]]]: key=genome name; val=dict: key=contig; val=list of Primers
    """
    # messages
    GAP = " " * 4
    MSG_1 = f"{GAP}getting shared ingroup kmers that appear once in each genome"
    MSG_2A = f"{GAP}evaluating "
    MSG_2B = " kmers"
    MSG_3A = f"{GAP}identified "
    MSG_3B = " candidate kmers"
    ERR_MSG_1 = "failed to identify a set of kmers shared between the ingroup genomes"
    ERR_MSG_2 = "none of the ingroup kmers are suitable for use as a primer"

    # initialize clock
    clock = Clock()

    # setup debugger
    params.log.rename(_getAllCandidateKmers.__name__)

    if sharedExists:
        # load existing shared kmers from file
        primers = params.loadObj(params.pickles[Parameters._SHARED])

        # determine the number of kmers that were identified
        numCand = len(next(iter(primers.values())))

    else:
        # get all non-duplicated kmers that are shared in the ingroup
        params.log.info(MSG_1)
        clock.printStart(MSG_1)
        primers = __getSharedPrimers(params)
        clock.printDone()

        # move log back to this function
        params.log.rename(_getAllCandidateKmers.__name__)
        params.log.info(f"{GAP}done {clock.getTimeString()}")

        # determine the number of kmers that were identified
        numCand = len(next(iter(primers.values())))

        # make sure that ingroup kmers were identified
        if numCand == 0:
            params.log.error(ERR_MSG_1)
            raise RuntimeError(ERR_MSG_1)

        # dump the shared kmers to file
        params.dumpObj(
            primers, params.pickles[Parameters._SHARED], "shared kmers", prefix=GAP
        )

    # print status
    clock.printStart(f"{MSG_2A}{numCand}{MSG_2B}")
    params.log.info(f"{MSG_2A}{numCand}{MSG_2B}")

    # evaluate the candidates to remove bad primers
    __removeBadPrimers(primers, params)

    # make sure candidates were found
    numCand = len(next(iter(primers.values())))
    if numCand == 0:
        params.log.error(ERR_MSG_2)
        raise RuntimeError(ERR_MSG_2)

    # build the output
    out = __buildOutput(primers)

    # print status
    clock.printDone()
    print(f"{MSG_3A}{numCand}{MSG_3B}")

    # log status
    params.log.info(f"{GAP}done {clock.getTimeString()}")
    params.log.info(f"{MSG_3A}{numCand}{MSG_3B}")

    # dump the candidate kmers to file
    params.dumpObj(
        primers, params.pickles[Parameters._CAND], "candidate kmers", prefix=GAP
    )

    return out
