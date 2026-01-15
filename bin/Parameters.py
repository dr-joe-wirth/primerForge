from __future__ import annotations

import argparse
import os
import pathlib
import pickle
import sys
import uuid
from dataclasses import dataclass
from typing import Any

from Bio import SeqIO
from primer3.bindings import DEFAULT_P3_ARGS

from bin import __author__, __version__
from bin.Clock import Clock
from bin.Log import Log


@dataclass
class DefaultArgs:
    """class to store default arguments"""
    out: str = "results.tsv"
    bed_file: str = "primers.bed"
    format: str = "genbank"
    min_len: int = 16
    max_len: int = 20
    min_gc: float = 40.0
    max_gc: float = 60.0
    min_tm: float = 55.0
    max_tm: float = 68.0
    min_pcr: int = 120
    max_pcr: int = 2400
    tm_diff: float = 5.0
    num_threads: int = 1
    primer3_mv_conc: float = DEFAULT_P3_ARGS.mv_conc
    primer3_dv_conc: float = DEFAULT_P3_ARGS.dv_conc
    primer3_dntp_conc: float = DEFAULT_P3_ARGS.dntp_conc
    primer3_dna_conc: float = DEFAULT_P3_ARGS.dna_conc
    primer3_temp_c: float = DEFAULT_P3_ARGS.temp_c
    primer3_max_loop: int = DEFAULT_P3_ARGS.max_loop
    ispcr_min_good: int = 6
    ispcr_min_perfect: int = 8
    ispcr_tile_size: int = 10
    temp_tolerance: float = 5.0
    max_repeats: int = 3
    bin_size: int = 64
    keep: bool = False
    debug: bool = False


class Parameters:
    """class to store arguments and debug utilities"""
    # constants
    _MIN_LEN = 15
    _MAX_LEN = 32
    __ALLOWED_FORMATS = ("genbank", "fasta")
    __ALL_CONTIGS_FNA = "all_contigs.fna"
    _WORKDIR_PREFIX = "primerforge_"
    _PARAMS = 0
    _SHARED = 1
    _CAND = 2
    _PAIR_1 = 3
    _PAIR_2 = 4
    _PAIR_3 = 5
    __PICKLE_FNS = {
        _PARAMS: "parameters.p",
        _SHARED: "sharedKmers.p",
        _CAND: "candidates.p",
        _PAIR_1: "pairs.p",
        _PAIR_2: "pairs_noOutgroup.p",
        _PAIR_3: "pairs_noOutgroup_validated.p",
    }

    # default values
    _DEF_RESULTS_FN = "results.tsv"
    _DEF_BED_FN = "primers.bed"
    _DEF_OUTGROUP = list()
    _DEF_FRMT = __ALLOWED_FORMATS[0]
    _DEF_MIN_LEN = 16
    _DEF_MAX_LEN = 20
    _DEF_MIN_GC = 40.0
    _DEF_MAX_GC = 60.0
    _DEF_MIN_TM = 55.0
    _DEF_MAX_TM = 68.0
    _DEF_MIN_PCR = 120
    _DEF_MAX_PCR = 2400
    _DEF_MAX_TM_DIFF = 5.0
    _DEF_NUM_THREADS = 1
    _DEF_ISPCR_MIN_GOOD = 6
    _DEF_ISPCR_MIN_PERFECT = 8
    _DEF_ISPCR_TILE_SIZE = 10
    _DEF_DEGREES = 5.0
    _DEF_REPEATS = 3
    _DEF_BINSIZE = 64
    _DEF_KEEP = False
    _DEF_DEBUG = False
    _DEF_HELP = False

    # overloads
    def __init__(self, argparse_ns: argparse.Namespace, initializeLog: bool = True):
        """Initialize a new Parameters object.

        Args:
            argparse_ns: parsed command-line arguments
            initializeLog: if True, initialize the log. Defaults to True.
        """
        self.ingroupFns: list[pathlib.Path] = argparse_ns.ingroup
        self.outgroupFns: list[pathlib.Path] = argparse_ns.outgroup
        self.resultsFn: pathlib.Path = argparse_ns.out
        self.bedFn: pathlib.Path = argparse_ns.bed_file
        self.format: str = argparse_ns.format
        self.minLen: int = argparse_ns.min_len
        self.maxLen: int = argparse_ns.max_len
        self.minGc: float = argparse_ns.min_gc
        self.maxGc: float = argparse_ns.max_gc
        self.minTm: float = argparse_ns.min_tm
        self.maxTm: float = argparse_ns.max_tm
        self.minPcr: int = argparse_ns.min_pcr
        self.maxPcr: int = argparse_ns.max_pcr
        self.maxTmDiff: float = argparse_ns.tm_diff
        self.disallowedLens: range = range(argparse_ns.min_bad, argparse_ns.max_bad + 1)
        self.numThreads: int = argparse_ns.num_threads
        self.debug: bool = argparse_ns.debug
        self.log: Log = Log(debug=self.debug, initialize=initializeLog)
        self.__workdir: pathlib.Path
        self.pickles: dict[int, pathlib.Path] = {}
        self.allContigsFna: pathlib.Path

        self.keepIntermediateFiles: bool = argparse_ns.keep

        # advanced attributes for primer3
        self.p3_mvConc: float = argparse_ns.primer3_mv_conc
        self.p3_dvConc: float = argparse_ns.primer3_dv_conc
        self.p3_dntpConc: float = argparse_ns.primer3_dntp_conc
        self.p3_dnaConc: float = argparse_ns.primer3_dna_conc
        self.p3_tempC: float = argparse_ns.primer3_temp_c
        self.p3_maxLoop: int = argparse_ns.primer3_max_loop

        # advanced attributes for isPcr
        self.isPcr_minGood: int = argparse_ns.isPcr_minGood
        self.isPcr_minPerfect: int = argparse_ns.isPcr_minPerfect
        self.isPcr_tileSize: int = argparse_ns.isPcr_tileSize

        # additional advanced attributes
        self.tempTolerance: float = argparse_ns.temp_tolerance
        self.maxRepeatLen: int = argparse_ns.max_repeats
        self.maxBinSize: int = argparse_ns.bin_size

        # save author and version as private attributes
        self.__author: list[str] = __author__
        self.__version: str = __version__

        # Sort the ingroup files by size (smallest first)
        self.ingroupFns.sort(key=lambda x: x.stat().st_size)

        # Now populate working directory + its linked attributes
        self.workdir = self.__getIntermediateDirname()

        # Finally, perform sanity checks
        # Check 1 -- all ingroup files must exist
        nonexistent = [str(x) for x in self.ingroupFns if not x.exists()]
        if len(nonexistent) > 0:
            raise FileNotFoundError(f"the following ingroup files do not exist: {', '.join(nonexistent)}")
        # Check 2 -- all outgroup files must exist
        nonexistent = [str(x) for x in self.outgroupFns if not x.is_file()]
        if len(nonexistent) > 0:
            raise FileNotFoundError(f"the following ingroup files do not exist: {', '.join(nonexistent)}")
        # Check 3 -- all ingroup and outgroup files must be in the specified format
        self.__checkGenomeFormat()
        # Check 4 -- the output files must be writable
        for output_file in (self.bedFn, self.resultsFn):
            self.__checkOutputFile(output_file)
        # Check 5 -- tileSize must be sanely set (must be less than min primer length)
        if self.isPcr_tileSize > Parameters._MIN_LEN:
            raise ValueError(f"maximum tileSize is {Parameters._MIN_LEN} bp")

    def __eq__(self, other: Parameters) -> bool:
        """equality overload

        Args:
            other (Parameters): another Parameters object

        Raises:
            TypeError: can only compare Parameters to Parameters

        Returns:
            bool: are objects equal?
        """
        # only compare to Parameters
        if not type(other) is Parameters:
            raise TypeError(f"cannot compare Parameters object to type '{type(other)}'")

        # determine if ingroup files match
        if self.ingroupFns is None:
            sameIngroup = self.ingroupFns == other.ingroupFns
        else:
            sameIngroup = set(map(os.path.abspath, self.ingroupFns)) == set(
                map(os.path.abspath, other.ingroupFns)
            )

        # determine if outgroup files match
        if self.outgroupFns is None:
            sameOutgroup = self.outgroupFns == other.outgroupFns
        else:
            sameOutgroup = set(map(os.path.abspath, self.outgroupFns)) == set(
                map(os.path.abspath, other.outgroupFns)
            )

        # determine if other important attributes match
        samePrimerLens = self.minLen == other.minLen and self.maxLen == other.maxLen
        samePrimerGc = self.minGc == other.minGc and self.maxGc == other.maxGc
        samePrimerTm = self.minTm == other.minTm and self.maxTm == other.maxTm
        samePcrLen = self.minPcr == other.minPcr and self.maxPcr == other.maxPcr
        sameTmDiff = self.maxTmDiff == other.maxTmDiff
        sameBadLens = self.disallowedLens == other.disallowedLens

        # determine if primer3 attributes match
        sameP3 = (
            self.p3_dnaConc == other.p3_dnaConc
            and self.p3_dntpConc == other.p3_dntpConc
            and self.p3_dvConc == other.p3_dvConc
            and self.p3_mvConc == other.p3_mvConc
            and self.p3_tempC == other.p3_tempC
            and self.p3_maxLoop == other.p3_maxLoop
        )

        # determine if isPcr attributes match
        sameIsPcr = (
            self.isPcr_minGood == other.isPcr_minGood
            and self.isPcr_minPerfect == other.isPcr_minPerfect
            and self.isPcr_tileSize == other.isPcr_tileSize
        )

        # determine if additional advanced attributes match
        sameAddtnl = (
            self.tempTolerance == other.tempTolerance
            and self.maxRepeatLen == other.maxRepeatLen
            and self.maxBinSize == other.maxBinSize
        )

        # evaluate important attributes; all must be equivalent
        return all(
            (
                sameIngroup,
                sameOutgroup,
                samePrimerLens,
                samePrimerGc,
                samePrimerTm,
                samePcrLen,
                sameTmDiff,
                sameBadLens,
                sameP3,
                sameIsPcr,
                sameAddtnl,
            )
        )

    def __ne__(self, other: Parameters) -> bool:
        """inequality overload

        Args:
            other (Parameters): another Parameters object

        Returns:
            bool: are objects not equal?
        """
        return not self == other

    # private methods
    @staticmethod
    def __checkOutputFile(fn: pathlib.Path) -> None:
        """checks if an output file is valid

        Args:
            fn (pathlib.Path): the filename to check

        Raises:
            FileExistsError: file already exists
            ValueError: file cannot be written to
        """
        # constants
        YN = ["y", "n"]
        WARN_MSG_A = "\nfile '"
        WARN_MSG_B = "' already exists."
        PROCEED_MSG = f"overwrite existing file? [{'/'.join(YN)}] "
        INVALID_SELECTION = "invalid selection"
        ERR_MSG = "cannot write to "

        # make sure the file doesn't already exist
        if fn.exists():
            # warn the user and abort if the file already exists
            print(f"{WARN_MSG_A}{fn}{WARN_MSG_B}")

            # ask the user if they wish to proceed
            proceed = input(PROCEED_MSG)
            while proceed not in YN:
                print(INVALID_SELECTION)
                proceed = input(PROCEED_MSG)

            # if the user declined, then abort
            if proceed == YN[1]:
                raise FileExistsError(f"{WARN_MSG_A}{fn}{WARN_MSG_B}")

        # make sure we can write the file to the output directory
        try:
            os.access(fn.parent, os.W_OK)
        except:
            raise ValueError(f"{ERR_MSG}{fn}")

    def __checkGenomeFormat(self) -> None:
        """checks the file format of the input genome files

        Raises:
            ValueError: empty or improperly formatted file encountered
        """
        # error message
        ERR_MSG = f" is empty or an improperly formatted {self.format} file"

        # initialize boolean to track status
        fail = False

        # for each genome file
        for fn in self.ingroupFns + self.outgroupFns:
            # open the file
            with open(fn, "r") as fh:
                # attempt to extract the first record from the generator
                try:
                    next(iter(SeqIO.parse(fh, self.format)))

                # failure indicates empty file or improperly formatted
                except StopIteration:
                    fail = True

            # raise an error only after the file is closed
            if fail:
                raise ValueError(f"{fn}{ERR_MSG}")

    def __getIntermediateDirname(self) -> pathlib.Path:
        """determines the intermediate directory name and handles checkpointing

        Returns:
            pathlib.Path: the absolute path to the intermediate directory
        """
        # constants
        CP_MSG_A = "\ncheckpoint detected. using cached data found in '"
        CP_MSG_B = "'\nresume previous run (yes; no; abort)? [y/n/a] "
        KILL_MSG = "\naborting\n"
        IGNORE_MSG = "\nignoring checkpoint; starting new run\n"
        INVALID_MSG = "invalid response"

        # initialize a variable to determine if an existing directory was found
        found = False

        uid = str(uuid.uuid4()).replace("-", "")  # don't keep hypens

        # find all existing Parameters pickles
        allParamFns = [x.joinpath(Parameters.__PICKLE_FNS[Parameters._PARAMS])
                       for x in pathlib.Path().cwd().glob(f"{Parameters._WORKDIR_PREFIX}*")]

        # if there are existing parameters, then process them
        if allParamFns != []:
            # sort them by time; most recently created files first
            allParamFns.sort(key=lambda x: os.path.getctime(x), reverse=True)

            # for each existing Parameter file
            for fn in allParamFns:
                # load the existing Parameters into memory
                with open(fn, "rb") as fh:
                    other: Parameters = pickle.load(fh)

                # check if it is compatible with the current run; use existing folder if possible
                if other == self:
                    workdir = os.path.dirname(fn)
                    found = True
                    break

        # ask the user if they want to checkpoint
        if found:
            while True:
                # print message that a checkpoint was found
                proceed = input(f"{CP_MSG_A}{workdir}{CP_MSG_B}")

                # let the user decide
                if proceed.startswith('n'):
                    print(IGNORE_MSG)
                    found = False
                    break

                elif proceed.startswith('y'):
                    print()
                    break

                elif proceed.startswith('a'):
                    print(KILL_MSG)
                    self.helpRequested = True
                    break

                else:
                    print(INVALID_MSG)


        # create new directory if not found or requested to skip existing
        if not found:
            workdir = Parameters._WORKDIR_PREFIX + uid

        return pathlib.Path(workdir).absolute()

    # public methods
    def logRunDetails(self) -> None:
        """saves the details for the current instance of the program"""
        # constant
        WIDTH = 34

        # save the command
        self.log.info(f'{"command:":<{WIDTH}}{" ".join(sys.argv)}')

        # write the parameters to the log file
        self.log.info(f'{"version:":<{WIDTH}}{self.__version}')
        self.log.info(f'{"ingroup:":<{WIDTH}}{",".join([str(x) for x in self.ingroupFns])}')
        self.log.info(f'{"outgroup:":<{WIDTH}}{",".join([str(x) for x in self.outgroupFns])}')
        self.log.info(f'{"results filename:":{WIDTH}}{self.resultsFn}')
        self.log.info(f'{"file format:":{WIDTH}}{self.format}')
        self.log.info(f'{"min kmer len:":{WIDTH}}{self.minLen}')
        self.log.info(f'{"max kmer len:":<{WIDTH}}{self.maxLen}')
        self.log.info(f'{"min % G+C":<{WIDTH}}{self.minGc}')
        self.log.info(f'{"max % G+C":<{WIDTH}}{self.maxGc}')
        self.log.info(f'{"min Tm:":<{WIDTH}}{self.minTm}')
        self.log.info(f'{"max Tm:":<{WIDTH}}{self.maxTm}')
        self.log.info(f'{"max Tm difference:":<{WIDTH}}{self.maxTmDiff}')
        self.log.info(f'{"min PCR size:":<{WIDTH}}{self.minPcr}')
        self.log.info(f'{"max PCR size:":<{WIDTH}}{self.maxPcr}')
        self.log.info(
            f'{"disallowed outgroup PCR sizes:":<{WIDTH}}{"-".join(map(str,[min(self.disallowedLens),max(self.disallowedLens)]))}'
        )
        self.log.info(f'{"num threads:":<{WIDTH}}{self.numThreads}')

    def dumpObj(self, obj: Any, fn: str, objName: str, prefix: str = "") -> None:
        """dumps an object in memory to file as a pickle

        Args:
            obj (any): the object to dump
            fn (str): the filename where object will be dumped
            objName (str): the name of the dumped object
            prefix (str, optional): a prefix for the printed message. Defaults to ''
        """
        # messages
        MSG_A = "dumping "
        MSG_B = " to '"
        MSG_C = "'"

        # start the timer
        clock = Clock()

        # get the filename
        fn = os.path.join(self.log.debugDir, fn)

        # determine which filename to print
        printedFn = fn[len(os.getcwd()) + 1 :]

        # print status
        self.log.info(MSG_A + objName + MSG_B + printedFn + MSG_C)
        clock.printStart(MSG_A + objName + MSG_B + printedFn + MSG_C, prefix=prefix)

        # dump the object to file
        with open(fn, "wb") as fh:
            pickle.dump(obj, fh)

        # print status
        clock.printDone()
        self.log.info(f"done {clock.getTimeString()}")

    def loadObj(self, fn: str):
        """loads an object from a pickle file

        Args:
            fn (str): a pickle filename

        Returns:
            any: the unpickled object
        """
        # constants
        MSG_A = "loading pickle from '"
        MSG_B = "'"

        # determine which filename to print
        printedFn = os.path.basename(fn)

        # start clock
        clock = Clock()

        # print status
        self.log.info(MSG_A + printedFn + MSG_B)
        clock.printStart(MSG_A + printedFn + MSG_B)

        # load the pickle
        with open(fn, "rb") as fh:
            out = pickle.load(fh)

        # print status
        clock.printDone()
        self.log.info("done " + clock.getTimeString())

        return out

    @property
    def workdir(self) -> pathlib.Path:
        """getter for the work directory

        Returns:
            pathlib.Path: the work directory
        """
        return self.__workdir

    @workdir.setter
    def workdir(self, value: pathlib.Path) -> None:
        """Update the working directory and all attributes that rely on it"""
        value = value.absolute()
        self.pickles = {x: value.joinpath(y) for x, y in self.__PICKLE_FNS.items()}
        self.allContigsFna = value.joinpath(self.__ALL_CONTIGS_FNA)
        self.__workdir = value
