from __future__ import annotations
from Bio import SeqIO
from bin.Log import Log
from bin.Clock import Clock
from primer3.bindings import DEFAULT_P3_ARGS
import getopt, glob, os, pickle, shutil, sys, uuid

class Parameters():
    """class to store arguments and debug utilities
    """
    # constants
    _MIN_LEN = 10
    __ALLOWED_FORMATS = ('genbank', 'fasta')
    __ALL_CONTIGS_FNA = "all_contigs.fna"
    _WORKDIR_PREFIX = "primerforge_"
    _PARAMS = 0
    _SHARED = 1
    _CAND   = 2
    _PAIR_1 = 3
    _PAIR_2 = 4
    _PAIR_3 = 5
    __PICKLE_FNS = {_PARAMS: "parameters.p",
                    _SHARED: "sharedKmers.p",
                    _CAND:   "candidates.p",
                    _PAIR_1: "pairs.p",
                    _PAIR_2: "pairs_noOutgroup.p",
                    _PAIR_3: "pairs_noOutgroup_validated.p"}
    
    # default values
    _DEF_RESULTS_FN = 'results.tsv'
    _DEF_BED_FN = 'primers.bed'
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
    _DEF_REPEATS = 4
    _DEF_BINSIZE = 64
    _DEF_KEEP = False
    _DEF_DEBUG = False
    _DEF_HELP = False
    
    # overloads
    def __init__(self, author:str, version:str, initializeLog:bool=True, makeWd:bool=True) -> Parameters:
        # type hint attributes
        self.ingroupFns:list[str]
        self.outgroupFns:list[str]
        self.resultsFn:str
        self.bedFn:str
        self.format:str
        self.minLen:int
        self.maxLen:int
        self.minGc:float
        self.maxGc:float
        self.minTm:float
        self.maxTm:float
        self.minPcr:int
        self.maxPcr:int
        self.maxTmDiff:float
        self.disallowedLens:range
        self.numThreads:int
        self.debug:bool
        self.helpRequested:bool
        self.log:Log
        self.pickles:dict[int,str]
        self.keepIntermediateFiles:bool
        self.allContigsFna:str
        
        # advanced attributes for primer3
        self.p3_mvConc:float
        self.p3_dvConc:float
        self.p3_dntpConc:float
        self.p3_dnaConc:float
        self.p3_tempC:float
        self.p3_maxLoop:int
        
        # advanced attributes for isPcr
        self.isPcr_minGood:int
        self.isPcr_minPerfect:int
        self.isPcr_tileSize:int
        
        # additional advanced attributes
        self.tempTolerance:float
        self.maxRepeatLen:int
        self.maxBinSize:int
        
        # save author and version as private attributes
        self.__author:str = author
        self.__version:str = version
        
        # parse the command line arguments (populates attributes)
        self.__parseArgs()
        
        # handle a few upkeep tasks if running
        if not self.helpRequested:
            # initialize a logger
            self.log = Log(debug=self.debug, initialize=initializeLog)
            
            # create the intermediate directory, or find an existing one for checkpointing
            workdir = self.__getIntermediateDirname(makeWd)
            
            # get the pickle filenames
            self.pickles = {x:os.path.join(workdir, y) for x,y in self.__PICKLE_FNS.items()}
            
            # get the all_contigs fasta filename
            self.allContigsFna = os.path.join(workdir, Parameters.__ALL_CONTIGS_FNA)

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
        if self.ingroupFns is None:
            sameIngroup = self.ingroupFns == other.ingroupFns
        else:
            sameIngroup = set(map(os.path.abspath, self.ingroupFns)) == set(map(os.path.abspath, other.ingroupFns))
        
        # determine if outgroup files match
        if self.outgroupFns is None:
            sameOutgroup = self.outgroupFns == other.outgroupFns
        else:  
            sameOutgroup = set(map(os.path.abspath, self.outgroupFns)) == set(map(os.path.abspath, other.outgroupFns))
        
        # determine if other important attributes match
        samePrimerLens = self.minLen == other.minLen and self.maxLen == other.maxLen
        samePrimerGc = self.minGc == other.minGc and self.maxGc == other.maxGc
        samePrimerTm = self.minTm == other.minTm and self.maxTm == other.maxTm
        samePcrLen = self.minPcr == other.minPcr and self.maxPcr == other.maxPcr
        sameTmDiff = self.maxTmDiff == other.maxTmDiff
        sameBadLens = self.disallowedLens == other.disallowedLens
        
        # determine if primer3 attributes match
        sameP3 = self.p3_dnaConc == other.p3_dnaConc and \
                 self.p3_dntpConc == other.p3_dntpConc and \
                 self.p3_dvConc == other.p3_dvConc and \
                 self.p3_mvConc == other.p3_mvConc and \
                 self.p3_tempC == other.p3_tempC and \
                 self.p3_maxLoop == other.p3_maxLoop
        
        # determine if isPcr attributes match
        sameIsPcr = self.isPcr_minGood == other.isPcr_minGood and \
                    self.isPcr_minPerfect == other.isPcr_minPerfect and \
                    self.isPcr_tileSize == other.isPcr_tileSize
        
        # determine if additional advanced attributes match
        sameAddtnl = self.tempTolerance == other.tempTolerance and\
                     self.maxRepeatLen == other.maxRepeatLen and \
                     self.maxBinSize == other.maxBinSize

        # evaluate important attributes; all must be equivalent
        return all((sameIngroup, sameOutgroup, samePrimerLens, samePrimerGc, samePrimerTm,
                    samePcrLen, sameTmDiff, sameBadLens, sameP3, sameIsPcr, sameAddtnl))
    
    def __ne__(self, other:Parameters) -> bool:
        """inequality overload

        Args:
            other (Parameters): another Parameters object

        Returns:
            bool: are objects not equal?
        """
        return not self == other
    
    # private methods
    def __checkOutputFile(fn:str) -> None:
        """checks if an output file is valid

        Args:
            fn (str): the filename to check

        Raises:
            FileExistsError: file already exists
            ValueError: file cannot be written to
        """
        # constants
        YN = ['y', 'n']
        WARN_MSG_A = "\nfile '"
        WARN_MSG_B = "' already exists."
        PROCEED_MSG = f"overwrite existing file? [{'/'.join(YN)}] "
        INVALID_SELECTION = 'invalid selection'
        ERR_MSG  = 'cannot write to '
        
        # make sure the file doesn't already exist
        if os.path.exists(fn):
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
            os.access(os.path.dirname(fn), os.W_OK)
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
            with open(fn, 'r') as fh:
                # attempt to extract the first record from the generator
                try:
                    next(iter(SeqIO.parse(fh, self.format)))
                
                # failure indicates empty file or improperly formatted
                except StopIteration:
                    fail = True
            
            # raise an error only after the file is closed
            if fail:
                raise ValueError(f"{fn}{ERR_MSG}")
    
    def __isPcrInstalled() -> None:
        """determines if isPcr is installed and in the PATH

        Raises:
            BaseException: isPcr is not installed or not in the PATH
        """
        # constant
        IS_PCR = "isPcr"
        ERR_MSG = ' is not installed or not in the PATH'
        
        # make sure that isPcr is in the path
        if not shutil.which(IS_PCR):
            raise BaseException(f"'{IS_PCR}'{ERR_MSG}")
    
    def __checkInstallation(self) -> None:
        """checks the installation of all primerforge dependencies

        Raises:
            BaseException: incompatible python version
            BaseException: biopython not installed
            BaseException: biopython version is bad
            BaseException: numpy not installed
            BaseException: primer3-py not installed
            BaseException: primer3-py version is bad
            BaseException: scipy not installed
            BaseException: scipy version is bad
        """
        # constants
        PY_MAJOR = 3
        PY_MINOR = (9, 11)
        BIO_VER = (1, 81)
        KMR_VER = (2, 1)
        P3_VER = 2
        SCI_VER = (1, 10)
        
        # messages
        BAD_VER = ' version is incompatible (requires '
        NOT_INS = ' package is not installed'
        SUCCESS = f'primerForge v{self.__version} is properly installed'
        
        # check python version
        if sys.version_info.major != PY_MAJOR or not (PY_MINOR[0] <= sys.version_info.minor <= PY_MINOR[1]):
            raise BaseException(f'incompatible python version (requires {PY_MAJOR}.{PY_MINOR[0]} to {PY_MAJOR}.{PY_MINOR[1]})')
    
        # check that all dependencies exist
        # check Bio installation
        try:
            import Bio
        except:
            raise BaseException(f"'Bio'{NOT_INS}")
        
         # check bio version
        vers = tuple(map(int, Bio.__version__.split('.')))
        if vers[0] < BIO_VER[0] or (vers[0] == BIO_VER[0] and vers[1] < BIO_VER[1]):
            raise BaseException(f"'Bio'{BAD_VER}{'.'.join(map(str, BIO_VER))})")
        
        # check khmer installation
        try:
            import khmer
        except:
            raise BaseException(f"'khmer'{NOT_INS}")
        
        # check khmer version
        vers = tuple(map(int, khmer.__version__.split('.')))
        if vers[0] < KMR_VER[0] or (vers[0] == KMR_VER[0] and vers[1] < KMR_VER[1]):
            raise BaseException(f"'khmer'{KMR_VER}{'.'.join(map(str, KMR_VER))}")
        
        # check numpy installation
        try:
            import numpy
        except:
            raise BaseException(f"'numpy'{NOT_INS}")
        
        # check primer3-py installation
        try:
            import primer3
        except:
            raise BaseException(f"'primer3-py'{NOT_INS}")
        
        # check primer3-py version
        if int(primer3.__version__.split('.')[0]) < P3_VER:
            raise BaseException(f"'primer3-py'{BAD_VER}{P3_VER} or above)")
        
        # check pyahocorasick installation
        try:
            import ahocorasick
        except:
            raise BaseException(f"'pyahocorasick'{NOT_INS}")
        
        # check scipy installation
        try:
            import scipy
        except:
            raise BaseException(f"'scipy'{NOT_INS}")
        
        # check scipy version
        vers = tuple(map(int, scipy.__version__.split('.')))
        if vers[0] < SCI_VER[0] or (vers[0] >= SCI_VER[0] and vers[1] < SCI_VER[1]):
            raise BaseException(f"'scipy'{BAD_VER}{'.'.join(map(str,SCI_VER))} or above)")
        
        # make sure isPcr is installed
        Parameters.__isPcrInstalled()
        
        # print success message
        print(SUCCESS)
    
    def __parseArgs(self) -> None:
        """parses command line arguments

        Raises:
            BaseException: cannot use wildcards without enclosing them in quotes
            ValueError: invalid/missing ingroup files
            ValueError: invalid/missing ingroup files
            ValueError: invalid/missing outgroup files
            ValueError: invalid/missing outgroup files
            ValueError: invalid disallowed sizes: not two values
            ValueError: invalid disallowed sizes: not int
            ValueError: invalid file format
            ValueError: invalid primer length: wrong num arguments
            ValueError: invalid primer length: not int
            ValueError: invalid primer length: below minimum threshold
            ValueError: invalid GC: not two values
            ValueError: invalid GC: not numeric
            ValueError: invalid Tm: not two values
            ValueError: invalid Tm: not numeric
            ValueError: invalid PCR size: not two values
            ValueError: invalid PCR size: not int
            ValueError: Tm difference is not numeric
            ValueError: num threads is not an integer
            ValueError: primer3: mv_conc is not numeric
            ValueError: primer3: dv_conc is not numeric
            ValueError: primer3: dntp_conc is not numeric
            ValueError: primer3: dna_conc is not numeric
            ValueError: primer3: temp_c is not numeric
            ValueError: primer3: max_loop is not int
            ValueError: isPcr: minGood is not int
            ValueError: isPcr: minPerfect is not int
            ValueError: isPcr: tileSize is not int
            ValueError: isPcr: tileSize cannot exceed _MIN_LEN
            ValueError: temp tolerance is not int
            ValueError: max repeats is not int
            ValueError: bin size is not int
            ValueError: tileSize cannot exceed minimum of input primer size
            ValueError: ingroup file argument is missing
        """
        # constants
        SEP = ","
        
        # flags
        INGROUP_FLAGS = ('-i', '--ingroup')
        OUT_FLAGS = ('-o', '--out')
        BED_FLAGS = ('-B', '--bed_file')
        OUTGROUP_FLAGS = ('-u', '--outgroup')
        DISALLOW_FLAGS = ("-b", "--bad_sizes")
        FMT_FLAGS = ('-f', '--format')
        PRIMER_LEN_FLAGS = ('-p', '--primer_len')
        GC_FLAGS = ('-g', '--gc_range')
        TM_FLAGS = ('-t', '--tm_range')
        THREADS_FLAGS = ('-n', '--num_threads')
        PCR_LEN_FLAGS = ('-r', '--pcr_prod')
        TM_DIFF_FLAGS = ('-d', '--tm_diff')
        KEEP_FLAGS = ('-k', '--keep')
        VERSION_FLAGS = ('-v', '--version')
        HELP_FLAGS = ('-h', '--help')
        CHECK_FLAGS = ('--check_install',)
        DEBUG_FLAGS = ("--debug",)
        ADVANCED_FLAGS = ("--advanced",)
        P3_MV_CONC_FLAGS = ("--primer3_mv_conc",)
        P3_DV_CONC_FLAGS = ("--primer3_dv_conc",)
        P3_DNTP_CONC_FLAGS = ("--primer3_dntp_conc",)
        P3_DNA_CONC_FLAGS = ("--primer3_dna_conc",)
        P3_TEMP_C_FLAGS = ("--primer3_temp_c",)
        P3_MAX_LOOP_FLAGS = ("--primer3_max_loop",)
        IP_MIN_GOOD_FLAGS = ("--isPcr_minGood",)
        IP_MIN_PERF_FLAGS = ("--isPcr_minPerfect",)
        IP_TILESIZE_FLAGS = ("--isPcr_tileSize",)
        ADDTNL_DEGREES_FLAGS = ("--temp_tolerance",)
        ADDTNL_REPEATS_FLAGS = ("--max_repeats",)
        ADDTNL_BINSIZE_FLAGS = ("--bin_size",)
        SHORT_OPTS = INGROUP_FLAGS[0][-1] + ":" + \
                     OUT_FLAGS[0][-1] + ":" + \
                     BED_FLAGS[0][-1] + ":" + \
                     OUTGROUP_FLAGS[0][-1] + ":" + \
                     DISALLOW_FLAGS[0][-1] + ":" + \
                     FMT_FLAGS[0][-1] + ":" + \
                     PRIMER_LEN_FLAGS[0][-1] + ":" + \
                     GC_FLAGS[0][-1] + ":" + \
                     TM_FLAGS[0][-1] + ":" + \
                     PCR_LEN_FLAGS[0][-1] + ":" + \
                     TM_DIFF_FLAGS[0][-1] + ":" + \
                     THREADS_FLAGS[0][-1] + ":" + \
                     KEEP_FLAGS[0][-1] + \
                     VERSION_FLAGS[0][-1] + \
                     HELP_FLAGS[0][-1]
        LONG_OPTS = (INGROUP_FLAGS[1][2:] + "=",
                     OUT_FLAGS[1][2:] + "=",
                     BED_FLAGS[1][2:] + "=",
                     OUTGROUP_FLAGS[1][2:] + "=",
                     DISALLOW_FLAGS[1][2:] + "=",
                     FMT_FLAGS[1][2:] + "=",
                     PRIMER_LEN_FLAGS[1][2:] + "=",
                     GC_FLAGS[1][2:] + "=",
                     TM_FLAGS[1][2:] + "=",
                     PCR_LEN_FLAGS[1][2:] + "=",
                     TM_DIFF_FLAGS[1][2:] + "=",
                     KEEP_FLAGS[1][2:],
                     THREADS_FLAGS[1][2:] + "=",
                     VERSION_FLAGS[1][2:],
                     HELP_FLAGS[1][2:],
                     CHECK_FLAGS[0][2:],
                     DEBUG_FLAGS[0][2:],
                     ADVANCED_FLAGS[0][2:],
                     P3_MV_CONC_FLAGS[0][2:] + "=",
                     P3_DV_CONC_FLAGS[0][2:] + "=",
                     P3_DNTP_CONC_FLAGS[0][2:] + "=",
                     P3_DNA_CONC_FLAGS[0][2:] + "=",
                     P3_TEMP_C_FLAGS[0][2:] + "=",
                     P3_MAX_LOOP_FLAGS[0][2:] + "=",
                     IP_MIN_GOOD_FLAGS[0][2:] + "=",
                     IP_MIN_PERF_FLAGS[0][2:] + "=",
                     IP_TILESIZE_FLAGS[0][2:] + "=",
                     ADDTNL_DEGREES_FLAGS[0][2:] + "=",
                     ADDTNL_REPEATS_FLAGS[0][2:] + "=",
                     ADDTNL_BINSIZE_FLAGS[0][2:] + "=")

        # messages
        ERR_MSG_0  = 'detected wildcards that are not enclosed in quotes; aborting'
        ERR_MSG_1  = 'invalid or missing ingroup file(s)'
        ERR_MSG_2  = 'invalid or missing outgroup file(s)'
        ERR_MSG_3  = 'must specify exactly two integers for bad outgroup PCR product length'
        ERR_MSG_4  = 'bad outgroup PCR product lengths are not integers'
        ERR_MSG_5  = 'invalid format'
        ERR_MSG_6  = 'can only specify one primer length or a range (min,max)'
        ERR_MSG_7  = 'primer lengths are not integers'
        ERR_MSG_8  = 'minimum primer length is'
        ERR_MSG_9  = 'must specify a range of GC values (min,max)'
        ERR_MSG_10 = 'gc values are not numeric'
        ERR_MSG_11 = 'must specify a range of Tm values (min, max)'
        ERR_MSG_12 = 'Tm values are not numeric'
        ERR_MSG_13 = 'can only specify one PCR product length or a range (min,max)'
        ERR_MSG_14 = 'PCR product lengths are not integers'
        ERR_MSG_15 = 'max Tm difference is not numeric'
        ERR_MSG_16 = 'num threads is not an integer'
        ERR_MSG_17 = "expected 'float' for '"
        ERR_MSG_18 = "exepcted 'int' for '"
        ERR_MSG_19 = 'minimum tile size is'
        ERR_MSG_20 = 'tile size cannot exceed smallest primer length'
        ERR_MSG_21 = 'must specify one or more ingroup files'

        # helper function
        def printHelp(advanced:bool):
            # constants
            GAP = " "*4
            EOL = "\n"
            SEP_1 = ", "
            SEP_2 = "|"
            SEP_3 = ","
            DEF_OPEN = ' (default: '
            CLOSE = ')'
            WIDTH = max(map(lambda x: len('--' + x.replace("=", '')), LONG_OPTS)) + len(GAP) # argument list determines column width
            
            # primary help message
            HELP_MSG = f"{EOL}Finds pairs of primers suitable for a group of input genomes{EOL}" + \
                       f"{GAP}{self.__author}, 2024{EOL*2}" + \
                       f"If you use this software, please cite our article:{EOL}" + \
                       f"{GAP}https://doi.org/10.21105/joss.06850{EOL*2}" + \
                       f"usage:{EOL}" + \
                       f"{GAP}primerForge [-{SHORT_OPTS.replace(':','')}]{EOL*2}" + \
                       f"required arguments:{EOL}" + \
                       f'{GAP}{INGROUP_FLAGS[0] + SEP_1 + INGROUP_FLAGS[1]:<{WIDTH}}[file] ingroup filename or a file pattern inside double-quotes (eg."*.gbff"){EOL*2}' + \
                       f"optional arguments:{EOL}" + \
                       f"{GAP}{OUT_FLAGS[0] + SEP_1 + OUT_FLAGS[1]:<{WIDTH}}[file] output filename for primer pair data{DEF_OPEN}{Parameters._DEF_RESULTS_FN}{CLOSE}{EOL}" + \
                       f"{GAP}{BED_FLAGS[0] + SEP_1 + BED_FLAGS[1]:<{WIDTH}}[file] output filename for primer data in BED file format{DEF_OPEN}{Parameters._DEF_BED_FN}{CLOSE}{EOL}" + \
                       f'{GAP}{OUTGROUP_FLAGS[0] + SEP_1 + OUTGROUP_FLAGS[1]:<{WIDTH}}[file] outgroup filename or a file pattern inside double-quotes (eg."*.gbff"){EOL}' + \
                       f"{GAP}{DISALLOW_FLAGS[0] + SEP_1 + DISALLOW_FLAGS[1]:<{WIDTH}}[int,int] a range of PCR product lengths that the outgroup cannot produce{DEF_OPEN}same as '{PCR_LEN_FLAGS[1]}'{CLOSE}{EOL}" + \
                       f"{GAP}{FMT_FLAGS[0] + SEP_1 + FMT_FLAGS[1]:<{WIDTH}}[str] file format of the ingroup and outgroup [{Parameters.__ALLOWED_FORMATS[0]}{SEP_2}{Parameters.__ALLOWED_FORMATS[1]}]{DEF_OPEN}{Parameters._DEF_FRMT}{CLOSE}{EOL}" + \
                       f"{GAP}{PRIMER_LEN_FLAGS[0] + SEP_1 + PRIMER_LEN_FLAGS[1]:<{WIDTH}}[int(s)] a single primer length or a range specified as 'min,max'; (minimum {Parameters._MIN_LEN}){DEF_OPEN}{Parameters._DEF_MIN_LEN}{SEP_3}{Parameters._DEF_MAX_LEN}{CLOSE}{EOL}" + \
                       f"{GAP}{GC_FLAGS[0] + SEP_1 + GC_FLAGS[1]:<{WIDTH}}[float,float] a min and max percent GC specified as a comma separated list{DEF_OPEN}{Parameters._DEF_MIN_GC}{SEP_3}{Parameters._DEF_MAX_GC}{CLOSE}{EOL}" + \
                       f"{GAP}{TM_FLAGS[0] + SEP_1 + TM_FLAGS[1]:<{WIDTH}}[float,float] a min and max melting temp (°C) specified as a comma separated list{DEF_OPEN}{Parameters._DEF_MIN_TM}{SEP_3}{Parameters._DEF_MAX_TM}{CLOSE}{EOL}" + \
                       f"{GAP}{PCR_LEN_FLAGS[0] + SEP_1 + PCR_LEN_FLAGS[1]:<{WIDTH}}[int(s)] a single PCR product length or a range specified as 'min,max'{DEF_OPEN}{Parameters._DEF_MIN_PCR}{SEP_3}{Parameters._DEF_MAX_PCR}{CLOSE}{EOL}" + \
                       f"{GAP}{TM_DIFF_FLAGS[0] + SEP_1 + TM_DIFF_FLAGS[1]:<{WIDTH}}[float] the maximum allowable Tm difference °C between a pair of primers{DEF_OPEN}{Parameters._DEF_MAX_TM_DIFF}{CLOSE}{EOL}" + \
                       f"{GAP}{THREADS_FLAGS[0] + SEP_1 + THREADS_FLAGS[1]:<{WIDTH}}[int] the number of threads for parallel processing{DEF_OPEN}{Parameters._DEF_NUM_THREADS}{CLOSE}{EOL}" + \
                       f"{GAP}{KEEP_FLAGS[0] + SEP_1 + KEEP_FLAGS[1]:<{WIDTH}}keep intermediate files{DEF_OPEN}{Parameters._DEF_KEEP}{CLOSE}{EOL}" + \
                       f"{GAP}{VERSION_FLAGS[0] + SEP_1 + VERSION_FLAGS[1]:<{WIDTH}}print the version{EOL}" + \
                       f"{GAP}{HELP_FLAGS[0] + SEP_1 + HELP_FLAGS[1]:<{WIDTH}}print this message{EOL*2}" + \
                       f"{GAP}{CHECK_FLAGS[0]:<{WIDTH}}check installation{EOL}" + \
                       f"{GAP}{DEBUG_FLAGS[0]:<{WIDTH}}run in debug mode{EOL}" + \
                       f"{GAP}{ADVANCED_FLAGS[0]:<{WIDTH}}print advanced options{EOL}"
            
            # advanced help message
            ADV_MSG  = f"primer3 parameters:{EOL}" + \
                       f"{GAP}{P3_MV_CONC_FLAGS[0]:<{WIDTH}}[float] monovalent cation concentration (mM){DEF_OPEN}{DEFAULT_P3_ARGS.mv_conc}{CLOSE}{EOL}" + \
                       f"{GAP}{P3_DV_CONC_FLAGS[0]:<{WIDTH}}[float] divalent cation concentration (mM){DEF_OPEN}{DEFAULT_P3_ARGS.dv_conc}{CLOSE}{EOL}" + \
                       f"{GAP}{P3_DNTP_CONC_FLAGS[0]:<{WIDTH}}[float] dNTP concentration (mM){DEF_OPEN}{DEFAULT_P3_ARGS.dntp_conc}{CLOSE}{EOL}" + \
                       f"{GAP}{P3_DNA_CONC_FLAGS[0]:<{WIDTH}}[float] template DNA concentration (nM){DEF_OPEN}{DEFAULT_P3_ARGS.dna_conc}{CLOSE}{EOL}" + \
                       f"{GAP}{P3_TEMP_C_FLAGS[0]:<{WIDTH}}[float] simulation temp (°C) for ΔG{DEF_OPEN}{DEFAULT_P3_ARGS.temp_c}{CLOSE}{EOL}" + \
                       f"{GAP}{P3_MAX_LOOP_FLAGS[0]:<{WIDTH}}[int] maximum size (bp) of loops in primer secondary structures{DEF_OPEN}{DEFAULT_P3_ARGS.max_loop}{CLOSE}{EOL*2}" + \
                       f"isPcr parameters:{EOL}" + \
                       f"{GAP}{IP_MIN_GOOD_FLAGS[0]:<{WIDTH}}[int] minimum size (bp) where there must be 2 matches for each mismatch{DEF_OPEN}{Parameters._DEF_ISPCR_MIN_GOOD}{CLOSE}{EOL}" + \
                       f"{GAP}{IP_MIN_PERF_FLAGS[0]:<{WIDTH}}[int] minimum size (bp) of perfect match at 3' end of primer{DEF_OPEN}{Parameters._DEF_ISPCR_MIN_PERFECT}{CLOSE}{EOL}" + \
                       f"{GAP}{IP_TILESIZE_FLAGS[0]:<{WIDTH}}[int] the size of match that triggers an alignment{DEF_OPEN}{Parameters._DEF_ISPCR_TILE_SIZE}{CLOSE}{EOL*2}" + \
                       f"additional parameters:{EOL}" + \
                       f"{GAP}{ADDTNL_DEGREES_FLAGS[0]:<{WIDTH}}[float] minimum number of degrees (°C) below primer Tm allowed for secondary structure Tm{DEF_OPEN}{Parameters._DEF_DEGREES}{CLOSE}{EOL}" + \
                       f"{GAP}{ADDTNL_REPEATS_FLAGS[0]:<{WIDTH}}[int] maximum allowed length (bp) of homopolymers (repeats) in primer sequences{DEF_OPEN}{Parameters._DEF_REPEATS}{CLOSE}{EOL}" + \
                       f"{GAP}{ADDTNL_BINSIZE_FLAGS[0]:<{WIDTH}}[int] maximum allowed length (bp) of contiguous regions of overlapping primers (bins){DEF_OPEN}{Parameters._DEF_BINSIZE}{CLOSE}{EOL}"
            
            # print the help message
            print(HELP_MSG)
            
            # print the advanced help message only if requested
            if advanced:
                print(ADV_MSG)
            
        # set default values
        self.ingroupFns = None
        self.resultsFn = os.path.join(os.getcwd(), Parameters._DEF_RESULTS_FN)
        self.bedFn = os.path.join(os.getcwd(), Parameters._DEF_BED_FN)
        self.outgroupFns = Parameters._DEF_OUTGROUP
        self.disallowedLens = None # start as None; update at the end
        self.format = Parameters._DEF_FRMT
        self.minLen = Parameters._DEF_MIN_LEN
        self.maxLen = Parameters._DEF_MAX_LEN
        self.minGc = Parameters._DEF_MIN_GC
        self.maxGc = Parameters._DEF_MAX_GC
        self.minTm = Parameters._DEF_MIN_TM
        self.maxTm = Parameters._DEF_MAX_TM
        self.minPcr = Parameters._DEF_MIN_PCR
        self.maxPcr = Parameters._DEF_MAX_PCR
        self.maxTmDiff = Parameters._DEF_MAX_TM_DIFF
        self.numThreads = Parameters._DEF_NUM_THREADS
        self.keepIntermediateFiles = Parameters._DEF_KEEP
        self.debug = Parameters._DEF_DEBUG
        self.helpRequested = Parameters._DEF_HELP
        self.isPcr_minGood = Parameters._DEF_ISPCR_MIN_GOOD
        self.isPcr_minPerfect = Parameters._DEF_ISPCR_MIN_PERFECT
        self.isPcr_tileSize = Parameters._DEF_ISPCR_TILE_SIZE
        self.tempTolerance = Parameters._DEF_DEGREES
        self.maxRepeatLen = Parameters._DEF_REPEATS
        self.maxBinSize = Parameters._DEF_BINSIZE
        
        # import default primer3 values
        self.p3_mvConc = DEFAULT_P3_ARGS.mv_conc
        self.p3_dvConc = DEFAULT_P3_ARGS.dv_conc
        self.p3_dntpConc = DEFAULT_P3_ARGS.dntp_conc
        self.p3_dnaConc = DEFAULT_P3_ARGS.dna_conc
        self.p3_tempC = DEFAULT_P3_ARGS.temp_c
        self.p3_maxLoop = DEFAULT_P3_ARGS.max_loop
        
        # check for the advanced option flag
        if ADVANCED_FLAGS[0] in sys.argv:
            self.helpRequested = True
            printHelp(True)
        
        # give help if requested
        elif HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv or len(sys.argv) == 1:
            self.helpRequested = True
            printHelp(False)
        
        # print the version if requested
        elif VERSION_FLAGS[0] in sys.argv or VERSION_FLAGS[1] in sys.argv:
            self.helpRequested = True
            print(self.__version)
        
        # check the installation if requested
        elif CHECK_FLAGS[0] in sys.argv:
            self.helpRequested = True
            self.__checkInstallation()
        
        else:
            # make sure isPcr is installed
            Parameters.__isPcrInstalled()
            
            # parse command line arguments
            opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
            
            # abort if args is not empty (happens if wildcards are not enclosed in quotes)
            if args != []:
                raise BaseException(ERR_MSG_0)
            
            for opt,arg in opts:
                # get the ingroup filenames
                if opt in INGROUP_FLAGS:
                    self.ingroupFns = glob.glob(arg)
                    for fn in self.ingroupFns:
                        if not os.path.isfile(fn):
                            raise ValueError(ERR_MSG_1)
                    if self.ingroupFns == []:
                        raise ValueError(ERR_MSG_1)
                
                # get output filename
                elif opt in OUT_FLAGS:                    
                    self.resultsFn = os.path.abspath(arg)
                
                # get the BED filename
                elif opt in BED_FLAGS:
                    self.bedFn = os.path.abspath(arg)
                
                # get the outgroup filenames
                elif opt in OUTGROUP_FLAGS:
                    self.outgroupFns = glob.glob(arg)
                    for fn in self.outgroupFns:
                        if not os.path.isfile(fn):
                            raise ValueError(ERR_MSG_2)
                    if self.outgroupFns == []:
                        raise ValueError(ERR_MSG_2)
                
                # get the disallowed outgroup pcr product lengths
                elif opt in DISALLOW_FLAGS:
                    # split comma-separated list
                    disallowed = arg.split(SEP)
                    
                    # make sure exactly two values provided
                    if len(disallowed) != 2:
                        raise ValueError(ERR_MSG_3)
                    
                    # coerce lengths to ints
                    try:
                        disallowed = [int(x) for x in disallowed]
                    except:
                        raise ValueError(ERR_MSG_4)
                    
                    # save the values
                    self.disallowedLens = range(min(disallowed), max(disallowed)+1)
                
                # get the file format
                elif opt in FMT_FLAGS:
                    if arg not in Parameters.__ALLOWED_FORMATS:
                        raise ValueError(ERR_MSG_5)
                    self.format = arg
                
                # get the primer lengths
                elif opt in PRIMER_LEN_FLAGS:
                    # split comma-separated list
                    primerRange = arg.split(SEP)
                    
                    # make sure at one or two primers specified
                    if len(primerRange) not in {1,2}:
                        raise ValueError(ERR_MSG_6)
                    
                    # coerce to lengths to ints
                    try:
                        primerRange = [int(x) for x in primerRange]
                    except:
                        raise ValueError(ERR_MSG_7)
                    
                    # make sure that the minimum is within the allowed range
                    if min(primerRange) < Parameters._MIN_LEN:
                        raise ValueError(f"{ERR_MSG_8} {Parameters._MIN_LEN} bp")
                    
                    # save values
                    self.minLen = min(primerRange)
                    self.maxLen = max(primerRange)
                
                # get the allowed GC range
                elif opt in GC_FLAGS:
                    # expecting two values separated by a comma
                    gcRange = arg.split(SEP)
                    if len(gcRange) != 2:
                        raise ValueError(ERR_MSG_9)
                    
                    # make sure the values are numeric
                    try:
                        gcRange = [float(x) for x in gcRange]
                    except:
                        raise ValueError(ERR_MSG_10)
                
                    # save values
                    self.minGc = min(gcRange)
                    self.maxGc = max(gcRange)
                
                # get the allowed Tm range
                elif opt in TM_FLAGS:
                    # expecting two values separated by a comma
                    tmRange = arg.split(SEP)
                    if len(tmRange) != 2:
                        raise ValueError(ERR_MSG_11)
                
                    # make sure the values are numeric
                    try:
                        tmRange = [float(x) for x in tmRange]
                    except:
                        raise ValueError(ERR_MSG_12)
                
                    # save values
                    self.minTm = min(tmRange)
                    self.maxTm = max(tmRange)
            
                # get the allowed PCR lengths
                elif opt in PCR_LEN_FLAGS:
                    # expecting one or two values
                    pcrRange = arg.split(SEP)
                    if len(pcrRange) not in {1,2}:
                        raise ValueError(ERR_MSG_13)
                    
                    # coerce to integers
                    try:
                        pcrRange = [int(x) for x in pcrRange]
                    except:
                        raise ValueError(ERR_MSG_14)
                
                    # save values
                    self.minPcr = min(pcrRange)
                    self.maxPcr = max(pcrRange)
                    
                    # see if the disallowed needs to be changed
                    if self.disallowedLens is None:
                        self.disallowedLens = range(self.minPcr,self.maxPcr+1)
                
                # get the allowed Tm difference between primer pairs
                elif opt in TM_DIFF_FLAGS:
                    # make sure input is numeric
                    try:
                        self.maxTmDiff = float(arg)
                    except:
                        raise ValueError(ERR_MSG_15)
                
                # get the number of threads to use
                elif opt in THREADS_FLAGS:
                    # make sure input is an integer
                    try:
                        self.numThreads = int(arg)
                    except:
                        raise ValueError(ERR_MSG_16)
                
                # update keep to True if requested
                if opt in KEEP_FLAGS:
                    self.keepIntermediateFiles = True
                
                # update debug to True if requested
                elif opt in DEBUG_FLAGS:
                    self.debug = True
                    self.keepIntermediateFiles = True
                
                # get the primer3 parameters
                elif opt in P3_MV_CONC_FLAGS:
                    try:
                        self.p3_mvConc = float(arg)
                    except ValueError:
                        raise ValueError(f"{ERR_MSG_17}{opt}'")
                
                elif opt in P3_DV_CONC_FLAGS:
                    try:
                        self.p3_dvConc = float(arg)
                    except ValueError:
                        raise ValueError(f"{ERR_MSG_17}{opt}'")
                
                elif opt in P3_DNTP_CONC_FLAGS:
                    try:
                        self.p3_dntpConc = float(arg)
                    except ValueError:
                        raise ValueError(f"{ERR_MSG_17}{opt}'")
                
                elif opt in P3_DNA_CONC_FLAGS:
                    try:
                        self.p3_dnaConc = float(arg)
                    except ValueError:
                        raise ValueError(f"{ERR_MSG_17}{opt}'")
                
                elif opt in P3_TEMP_C_FLAGS:
                    try:
                        self.p3_tempC = float(arg)
                    except ValueError:
                        raise ValueError(f"{ERR_MSG_17}{opt}'")
                
                elif opt in P3_MAX_LOOP_FLAGS:
                    try:
                        self.p3_maxLoop = int(arg)
                    except ValueError:
                        raise ValueError(f"{ERR_MSG_18}{opt}'")
                
                # get isPcr parameters
                elif opt in IP_MIN_GOOD_FLAGS:
                    try:
                        self.isPcr_minGood = int(arg)
                    except ValueError:
                        raise ValueError(f"{ERR_MSG_18}{opt}'")
                
                elif opt in IP_MIN_PERF_FLAGS:
                    try:
                        self.isPcr_minPerfect = int(arg)
                    except ValueError:
                        raise ValueError(f"{ERR_MSG_18}{opt}'")
                
                elif opt in IP_TILESIZE_FLAGS:
                    try:
                        self.isPcr_tileSize = int(arg)
                    except ValueError:
                        raise ValueError(f"{ERR_MSG_18}{opt}'")
                    
                    # tilesize cannot be lower than the minimum primer length
                    if self.isPcr_tileSize < Parameters._MIN_LEN:
                        raise ValueError(f"{ERR_MSG_19} {Parameters._MIN_LEN} bp")
                
                # get additional advanced parameters
                elif opt in ADDTNL_DEGREES_FLAGS:
                    try:
                        self.tempTolerance = float(arg)
                    except ValueError:
                        raise ValueError(f"{ERR_MSG_17}{opt}")
                
                elif opt in ADDTNL_REPEATS_FLAGS:
                    try:
                        self.maxRepeatLen = int(arg)
                    except ValueError:
                        raise ValueError(f"{ERR_MSG_18}{opt}")
                
                elif opt in ADDTNL_BINSIZE_FLAGS:
                    try:
                        self.maxBinSize = int(arg)
                    except ValueError:
                        raise ValueError(f"{ERR_MSG_18}{opt}")
            
            # make sure that the tileSize is not larger than the smallest primer length
            if self.isPcr_tileSize > self.minLen:
                raise ValueError(ERR_MSG_20)
            
            # update disallowed to match pcr parameters unless it was already specified
            if self.disallowedLens is None:
                self.disallowedLens = range(self.minPcr, self.maxPcr + 1)
            
            # make sure an input file was specified
            if self.ingroupFns is None:
                raise ValueError(ERR_MSG_21)
            
            # make sure that all genomes are formatted correctly
            self.__checkGenomeFormat()
            
            # check for existing files
            Parameters.__checkOutputFile(self.resultsFn)
            Parameters.__checkOutputFile(self.bedFn)
    
    def __getIntermediateDirname(self, makeWd:bool) -> str:
        """determines the intermediate directory name and handles checkpointing

        Args:
            makeWd (bool): should the intermediate directory me made? (exists for testing)

        Returns:
            str: the name of the intermediate directory
        """
        # constants
        CP_MSG_A = "\ncheckpoint detected. using cached data found in '"
        CP_MSG_B = "'\nresume previous run (yes; no; abort)? [y/n/a] "
        KILL_MSG = "\naborting\n"
        IGNORE_MSG = "\nignoring checkpoint; starting new run\n"
        INVALID_MSG = "invalid response"
        
        # initialize a variable to determine if an existing directory was found
        found = False
        
        # find all existing Parameters pickles
        allParamFns = glob.glob(os.path.join(os.curdir, f'{Parameters._WORKDIR_PREFIX}*', Parameters.__PICKLE_FNS[Parameters._PARAMS]))
        
        # if there are existing parameters, then process them
        if allParamFns != []:
            # sort them by time; most recently created files first
            allParamFns.sort(key=lambda x: os.path.getctime(x), reverse=True)
            
            # for each existing Parameter file
            for fn in allParamFns:
                # load the existing Parameters into memory
                with open(fn, 'rb') as fh:
                    other:Parameters = pickle.load(fh)
                
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
                if proceed == 'n' or proceed == 'no':
                    print(IGNORE_MSG)
                    found = False
                    break
                
                elif proceed == 'y' or proceed == 'yes':
                    print()
                    break
                
                elif proceed == 'a' or proceed == 'abort':
                    print(KILL_MSG)
                    self.helpRequested = True
                    break
                
                else:
                    print(INVALID_MSG)
        
        # create new directory if not found or requested to skip existing
        if not found:
            uid = str(uuid.uuid4()).replace('-','') # don't keep hypens
            workdir = Parameters._WORKDIR_PREFIX + uid
            
            if makeWd:
                os.mkdir(workdir)
        
        return workdir
    
    # public methods
    def logRunDetails(self) -> None:
        """saves the details for the current instance of the program
        """
        # constant
        WIDTH = 34
        
        # save the command
        self.log.info(f'{"command:":<{WIDTH}}{" ".join(sys.argv)}')
        
        # write the parameters to the log file
        self.log.info(f'{"version:":<{WIDTH}}{self.__version}')
        self.log.info(f'{"ingroup:":<{WIDTH}}{",".join(self.ingroupFns)}')
        self.log.info(f'{"outgroup:":<{WIDTH}}{",".join(self.outgroupFns)}')
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
        self.log.info(f'{"disallowed outgroup PCR sizes:":<{WIDTH}}{"-".join(map(str,[min(self.disallowedLens),max(self.disallowedLens)]))}')
        self.log.info(f'{"num threads:":<{WIDTH}}{self.numThreads}')
    
    def dumpObj(self, obj:any, fn:str, objName:str, prefix:str='') -> None:
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
        printedFn = fn[len(os.getcwd())+1:]
        
        # print status
        self.log.info(MSG_A + objName + MSG_B + printedFn + MSG_C)
        clock.printStart(MSG_A + objName + MSG_B + printedFn + MSG_C, prefix=prefix)
        
        # dump the object to file
        with open(fn, 'wb') as fh:
            pickle.dump(obj, fh)
        
        # print status
        clock.printDone()
        self.log.info(f'done {clock.getTimeString()}')

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
        with open(fn, 'rb') as fh:
            out = pickle.load(fh)
            
        # print status
        clock.printDone()
        self.log.info('done ' + clock.getTimeString())
        
        return out
