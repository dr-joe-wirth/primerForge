from __future__ import annotations

import argparse
import os
import pathlib
import pickle
import sys
import uuid
from dataclasses import dataclass, field
from typing import Any

from Bio import SeqIO
from primer3.bindings import DEFAULT_P3_ARGS

from bin import __author__, __version__
from bin.Clock import Clock
from bin.Genome import Genome, InvalidFileFormat
from bin.Log import Log


@dataclass
class DefaultArgs:
    """class to store default arguments"""
    OUT_FN:str = "results.tsv"
    BED_FN:str = "primers.bed"
    FORMAT:str = "genbank"
    OUTGROUP:list[str] = field(default_factory=list) 
    MIN_LEN:int = 16
    MAX_LEN:int = 20
    MIN_GC:float = 40.0
    MAX_GC:float = 60.0
    MIN_TM:float = 55.0
    MAX_TM:float = 68.0
    MIN_PCR:int = 120
    MAX_PCR:int = 2400
    MAX_TM_DIFF:float = 5.0
    NUM_THREADS:int = 1
    PRIMER3_MV_CONC:float = DEFAULT_P3_ARGS.mv_conc
    PRIMER3_DV_CONC:float = DEFAULT_P3_ARGS.dv_conc
    PRIMER3_DNTP_CONC:float = DEFAULT_P3_ARGS.dntp_conc
    PRIMER3_DNA_CONC:float = DEFAULT_P3_ARGS.dna_conc
    PRIMER3_TEMP_C:float = DEFAULT_P3_ARGS.temp_c
    PRIMER3_MAX_LOOP:int = DEFAULT_P3_ARGS.max_loop
    ISPCR_MIN_GOOD:int = 6
    ISPCR_MIN_PERFECT:int = 8
    ISPCR_TILE_SIZE:int = 10
    TEMP_TOLERANCE:float = 5.0
    MAX_REPEATS:int = 3
    BIN_SIZE:int = 64
    PRE_LOAD:bool = False
    KEEP:bool = False
    DEBUG:bool = False


class Parameters:
    """class to store arguments and debug utilities"""
    # constants
    _MIN_ALLOWED_LEN = 15
    _MAX_ALLOWED_LEN = 32
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

    # overloads
    def __init__(self, parsedArgs:argparse.Namespace, initializeLog: bool = True):
        """Initialize a new Parameters object.

        Args:
            parsedArgs: parsed command-line arguments
            initializeLog: if True, initialize the log. Defaults to True.
        """
        self.ingroup:list[Genome] = [Genome(x, parsedArgs.pre_load) for x in parsedArgs.ingroup]
        self.outgroup:list[Genome] = [Genome(x, parsedArgs.pre_load) for x in parsedArgs.outgroup]
        self.resultsFn:pathlib.Path = parsedArgs.out
        self.bedFn:pathlib.Path = parsedArgs.bed_file
        self.minLen:int = parsedArgs.min_len
        self.maxLen:int = parsedArgs.max_len
        self.minGc:float = parsedArgs.min_gc
        self.maxGc:float = parsedArgs.max_gc
        self.minTm:float = parsedArgs.min_tm
        self.maxTm:float = parsedArgs.max_tm
        self.minPcr:int = parsedArgs.min_pcr
        self.maxPcr:int = parsedArgs.max_pcr
        self.maxTmDiff:float = parsedArgs.tm_diff
        self.disallowedLens:range = range(parsedArgs.min_bad, parsedArgs.max_bad + 1)
        self.pre_load:bool = parsedArgs.pre_load
        self.numThreads:int = parsedArgs.num_threads
        self.debug:bool = parsedArgs.debug
        self.log:Log = Log(debug=self.debug, initialize=initializeLog)
        self.__workdir:pathlib.Path
        self.pickles:dict[int, pathlib.Path] = dict()
        self.allContigsFna:pathlib.Path

        self.keepIntermediateFiles: bool = parsedArgs.keep

        # advanced attributes for primer3
        self.p3_mvConc:float = parsedArgs.primer3_mv_conc
        self.p3_dvConc:float = parsedArgs.primer3_dv_conc
        self.p3_dntpConc:float = parsedArgs.primer3_dntp_conc
        self.p3_dnaConc:float = parsedArgs.primer3_dna_conc
        self.p3_tempC:float = parsedArgs.primer3_temp_c
        self.p3_maxLoop:int = parsedArgs.primer3_max_loop

        # advanced attributes for isPcr
        self.isPcr_minGood:int = parsedArgs.isPcr_minGood
        self.isPcr_minPerfect:int = parsedArgs.isPcr_minPerfect
        self.isPcr_tileSize:int = parsedArgs.isPcr_tileSize

        # additional advanced attributes
        self.tempTolerance:float = parsedArgs.temp_tolerance
        self.maxRepeatLen:int = parsedArgs.max_repeats
        self.maxBinSize:int = parsedArgs.bin_size

        # Sort the ingroup files by size (smallest first)
        self.ingroup.sort()

        # Now populate working directory + its linked attributes
        self.workdir = self.__getIntermediateDirname()

        # Finally, perform sanity checks
        # Check 1 -- all ingroup files must exist
        nonexistent = [str(x.fn) for x in self.ingroup if not x.fn.exists()]
        if len(nonexistent) > 0:
            raise FileNotFoundError(f"the following ingroup files do not exist: {', '.join(nonexistent)}")
        
        # Check 2 -- all outgroup files must exist
        nonexistent = [str(x.fn) for x in self.outgroup if not x.fn.is_file()]
        if len(nonexistent) > 0:
            raise FileNotFoundError(f"the following ingroup files do not exist: {', '.join(nonexistent)}")
        
        # Check 3 -- all ingroup and outgroup files must be in an allowed format
        bad_files = [str(x.fn) for x in self.ingroup + self.outgroup if not x._format]
        if bad_files:
            raise InvalidFileFormat(f"the following files are invalidly formatted: {', '.join(bad_files)}")

        # Check 4 -- the output files must be writable
        for output_file in (self.bedFn, self.resultsFn):
            self.__checkOutputFile(output_file)
        
        # Check 5 -- primer lengths must be within range
        if self.minLen < Parameters._MIN_ALLOWED_LEN or self.maxLen > Parameters._MAX_ALLOWED_LEN:
            raise ValueError(f"invalid primer sizes: '{self.minLen}-{self.maxLen}'. primer lengths must be {Parameters._MIN_ALLOWED_LEN}-{Parameters._MAX_ALLOWED_LEN}bp")

        # Check 6 -- tileSize must be sanely set (must be less than min primer length)
        if self.isPcr_tileSize > Parameters._MIN_ALLOWED_LEN:
            raise ValueError(f"maximum tileSize is {Parameters._MIN_ALLOWED_LEN} bp")

    def __eq__(self, other:Parameters) -> bool:
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
        sameIngroup = set(self.ingroup) == set(other.ingroup)

        # determine if outgroup files match
        sameOutgroup = set(self.outgroup) == set(other.outgroup)

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

    def __ne__(self, other:Parameters) -> bool:
        """inequality overload

        Args:
            other (Parameters): another Parameters object

        Returns:
            bool: are objects not equal?
        """
        return not self == other

    # private methods
    @staticmethod
    def __checkOutputFile(fn:pathlib.Path) -> None:
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
        self.log.info(f'{"ingroup:":<{WIDTH}}{",".join([str(x) for x in self.ingroup])}')
        self.log.info(f'{"outgroup:":<{WIDTH}}{",".join([str(x) for x in self.outgroup])}')
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

    def dumpObj(self, obj:Any, fn:str, objName:str, prefix:str="") -> None:
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

    def loadObj(self, fn:str):
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
    def workdir(self, value:pathlib.Path) -> None:
        """Update the working directory and all attributes that rely on it"""
        value = value.absolute()
        self.pickles = {x: value.joinpath(y) for x, y in self.__PICKLE_FNS.items()}
        self.allContigsFna = value.joinpath(self.__ALL_CONTIGS_FNA)
        self.__workdir = value
