#!/usr/bin/env python3

import argparse
import multiprocessing as mp
import os
import pathlib
import shutil
import sys

from bin import __author__, __version__
from bin.Parameters import Parameters, DefaultArgs

DEFAULT_ARGS = DefaultArgs()

try:
    # Try to handle the fact that we may be in an HPC environment by using only scheduled CPU cores
    MAX_CPUS = len(os.sched_getaffinity(0))
except AttributeError:
    # Fall back to the standard method to allow all available cores
    MAX_CPUS = mp.cpu_count()


def _check_installation() -> None:
    """Check installation of all primerforge dependencies.

    Checks:
        - python version
        - biopython, numpy, primer3, and scipy versions
        - isPcr present
        - new primerForge kmer counter is compiled

    Raises:
        BaseException: if any dependency is missing or incompatible
    """
    # Constants
    PY_VER = (3, 9)
    BIO_VER = (1, 81)
    P3_VER = 2
    SCI_VER = (1, 10)

    # Check python version
    if sys.version_info.major < PY_VER[0] or (
            sys.version_info.major == PY_VER[0] and sys.version_info.minor < PY_VER[1]
    ):
        raise BaseException(
            f"incompatible python version (requires {PY_VER[0]}.{PY_VER[1]} or above)"
        )

    # Check dependencies
    try:
        import Bio
        vers = tuple(map(int, Bio.__version__.split(".")[:2]))
        if vers[0] < BIO_VER[0] or (vers[0] == BIO_VER[0] and vers[1] < BIO_VER[1]):
            raise BaseException(f"'Bio' version is incompatible (requires {'.'.join(map(str, BIO_VER))})")
    except ImportError:
        raise BaseException("'Bio' package is not installed")

    try:
        import numpy
    except ImportError:
        raise BaseException("'numpy' package is not installed")

    try:
        import primer3
        if int(primer3.__version__.split(".")[0]) < P3_VER:
            raise BaseException(f"'primer3-py' version is incompatible (requires {P3_VER} or above)")
    except ImportError:
        raise BaseException("'primer3-py' package is not installed")

    try:
        import scipy
        vers = tuple(map(int, scipy.__version__.split(".")[:2]))
        if vers[0] < SCI_VER[0] or (vers[0] >= SCI_VER[0] and vers[1] < SCI_VER[1]):
            raise BaseException(
                f"'scipy' version is incompatible (requires {'.'.join(map(str, SCI_VER))} or above)"
            )
    except ImportError:
        raise BaseException("'scipy' package is not installed")

    # Check isPcr
    if shutil.which("isPcr") is None:
        raise BaseException("'isPcr' is not installed or not in the PATH")

    # Check kmer counting package
    try:
        from bin.kmer_counting.kmer_counter import _decodeKmerEncoding
    except ImportError:
        raise BaseException("kmer counting package was not properly built")

    print(f"\nprimerForge is properly installed\n")


def _validate_min_max_ranges(arg_name: str, arg_value: str, must_have_two: bool = False) -> tuple[float, float]:
    """Take a comma-separated string of two numeric values corresponding to a range
    and return them cast as floats, sorted as (min_value, max_value).

    Args:
        arg_name: name of the argument (for error messages)
        arg_value: comma-separated string of two numbers (e.g., '15,21')
        must_have_two: if True, raise an error if only one value is given

    Returns:
        tuple[float, float]: (min_value, max_value)
    """
    # Split on commas
    values = [float(x) for x in arg_value.split(",")]

    # If only received one value and we were not told that two are required, use the one value in both positions
    if len(values) == 1 and not must_have_two:
        min_value, max_value = values[0], values[0]
    # If we didn't get exactly two values and we were told to expect two, raise a ValueError with helpful message
    elif len(values) != 2 and must_have_two:
        raise ValueError(f"must supply exactly two values for argument --{arg_name} (got: {values})")
    # If we got two, this is the normal case
    elif len(values) == 2:
        min_value, max_value = min(values), max(values)
    # Anything else is unexpected behavior
    else:
        raise ValueError(f"unexpected number of values for argument --{arg_name} (got: {values}")

    return min_value, max_value


class UnifiedHelpFormatter(argparse.HelpFormatter):
    def __init__(self, prog, indent_increment=2, max_help_position=24, width=None):
        super().__init__(prog, indent_increment, max_help_position, width)
        self.show_advanced = '--advanced' in sys.argv

    def _format_action_invocation(self, action):
        """Make sure metavar is GONE. It's ugly."""
        if not action.option_strings:
            return super()._format_action_invocation(action)
        return ', '.join(action.option_strings)

    def format_help(self):
        # Get basic help
        help_text = super().format_help()
        
        if self.show_advanced:
            # Add advanced sections
            help_text += self._format_advanced_help()
        
        return help_text

    @staticmethod
    def _format_advanced_help():
        GAP = " " * 2
        EOL = "\n"
        WIDTH = 21
        DEF_OPEN = " (default: "
        CLOSE = ")"
        
        return (
            f"{EOL}primer3 parameters:{EOL}"
            f"{GAP}{'--primer3_mv_conc':<{WIDTH}}[float] monovalent cation concentration (mM){DEF_OPEN}{DEFAULT_ARGS.primer3_mv_conc}{CLOSE}{EOL}"
            f"{GAP}{'--primer3_dv_conc':<{WIDTH}}[float] divalent cation concentration (mM){DEF_OPEN}{DEFAULT_ARGS.primer3_dv_conc}{CLOSE}{EOL}"
            f"{GAP}{'--primer3_dntp_conc':<{WIDTH}}[float] dNTP concentration (mM){DEF_OPEN}{DEFAULT_ARGS.primer3_dntp_conc}{CLOSE}{EOL}"
            f"{GAP}{'--primer3_dna_conc':<{WIDTH}}[float] template DNA concentration (nM){DEF_OPEN}{DEFAULT_ARGS.primer3_dna_conc}{CLOSE}{EOL}"
            f"{GAP}{'--primer3_temp_c':<{WIDTH}}[float] simulation temp (°C) for ΔG calculation{DEF_OPEN}{DEFAULT_ARGS.primer3_temp_c}{CLOSE}{EOL}"
            f"{GAP}{'--primer3_max_loop':<{WIDTH}}[int] maximum size (bp) of loops in primer secondary structures{DEF_OPEN}{DEFAULT_ARGS.primer3_max_loop}{CLOSE}{EOL}"
            f"{EOL}isPcr parameters:{EOL}"
            f"{GAP}{'--isPcr_minGood':<{WIDTH}}[int] minimum size (bp) where there must be 2 matches for each mismatch{DEF_OPEN}{DEFAULT_ARGS.ispcr_min_good}{CLOSE}{EOL}"
            f"{GAP}{'--isPcr_minPerfect':<{WIDTH}}[int] minimum size (bp) of perfect match at 3' end of primer{DEF_OPEN}{DEFAULT_ARGS.ispcr_min_perfect}{CLOSE}{EOL}"
            f"{GAP}{'--isPcr_tileSize':<{WIDTH}}[int] the size of match that triggers an alignment{DEF_OPEN}{DEFAULT_ARGS.ispcr_tile_size}{CLOSE}{EOL}"
            f"{EOL}additional parameters:{EOL}"
            f"{GAP}{'--temp_tolerance':<{WIDTH}}[float] minimum number of degrees (°C) below primer Tm allowed for secondary structure Tm{DEF_OPEN}{DEFAULT_ARGS.temp_tolerance}{CLOSE}{EOL}"
            f"{GAP}{'--max_repeats':<{WIDTH}}[int] maximum allowed length (bp) of homopolymers (repeats) in primer sequences{DEF_OPEN}{DEFAULT_ARGS.max_repeats}{CLOSE}{EOL}"
            f"{GAP}{'--bin_size':<{WIDTH}}[int] maximum allowed length (bp) of contiguous regions of overlapping primers (bins){DEF_OPEN}{DEFAULT_ARGS.bin_size}{CLOSE}{EOL}"
        )


def create_parser() -> argparse.ArgumentParser:
    """Create and return the ArgumentParser to be used for commandline argument parsing.

    Returns:
        an ArgumentParser configured for primerForge
    """
    parser = argparse.ArgumentParser(
        prog='primerForge',
        description='Finds pairs of primers suitable for a group of input genomes',
        usage='primerForge -i [-ouBubfpgtrdnkvh] [--check_install] [--debug] [--advanced]',
        formatter_class=UnifiedHelpFormatter,
        add_help=False  # We'll handle help manually
    )

    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--ingroup', required=True, nargs="+", type=pathlib.Path,
                          help='[file] ingroup filename or a file pattern (e.g. ingroup/*.gbff)')
    
    # Optional arguments  
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-o', '--out',  default=DEFAULT_ARGS.out, type=pathlib.Path,
                          help='[file] output filename for primer pair data (default: %(default)s)')
    optional.add_argument('-B', '--bed_file', default=DEFAULT_ARGS.bed_file, type=pathlib.Path,
                          help='[file] output filename for primer data in BED file format (default: %(default)s)')
    optional.add_argument('-u', '--outgroup', default=list(), nargs="+", type=pathlib.Path,
                          help='[file] outgroup filename or a file pattern (e.g. outgroup/*.gbff)')
    optional.add_argument('-b', '--bad_sizes',
                          help="[int,int] a range of PCR product lengths that the outgroup cannot produce (default: same as '--pcr_prod')")
    optional.add_argument('-f', '--format', default=DEFAULT_ARGS.format, choices=['genbank', 'fasta'],
                          help='[str] file format of the ingroup and outgroup [genbank|fasta] (default: %(default)s)')
    optional.add_argument('-p', '--primer_len', default=f"{DEFAULT_ARGS.min_len},{DEFAULT_ARGS.max_len}",
                          help='[int(s)] a single primer length or a range specified as \'min,max\'; (minimum 10; maximum 32) (default: %(default)s)')
    optional.add_argument('-g', '--gc_range', default=f"{DEFAULT_ARGS.min_gc},{DEFAULT_ARGS.max_gc}",
                          help='[float,float] a min and max percent GC specified as a comma separated list (default: %(default)s)')
    optional.add_argument('-t', '--tm_range', default=f"{DEFAULT_ARGS.min_tm},{DEFAULT_ARGS.max_tm}",
                          help='[float,float] a min and max melting temp (°C) specified as a comma separated list (default: %(default)s)')
    optional.add_argument('-r', '--pcr_prod', default='120,2400',
                          help='[int(s)] a single PCR product length or a range specified as \'min,max\' (default: 120,2400)')
    optional.add_argument('-d', '--tm_diff', type=float, default=5.0,
                          help='[float] the maximum allowable Tm difference °C between a pair of primers (default: 5.0)')
    optional.add_argument('-n', '--num_threads', type=int, default=1,
                          help='[int] the number of threads for parallel processing (default: 1)')
    optional.add_argument('-k', '--keep', action='store_true',
                          help='keep intermediate files (default: False)')
    optional.add_argument('-v', '--version', action='store_true',
                          help='print the version')
    optional.add_argument('-h', '--help', action='store_true',
                          help='print this message')
    
    # Special flags
    special = parser.add_argument_group('special arguments')
    special.add_argument('--check_install', action='store_true',
                         help='check that all dependencies are available')
    special.add_argument('--debug', action='store_true',
                         help='run in debug mode')
    special.add_argument('--advanced', action='store_true',
                         help='print advanced help options')
    
    # Advanced arguments (only shown when --advanced is used)
    # help=argparse.SUPPRESS --> argparse will accept and handle these args but never display them (we do that manually)
    # primer3 parameters
    parser.add_argument('--primer3_mv_conc', type=float, default=DEFAULT_ARGS.primer3_mv_conc, help=argparse.SUPPRESS)
    parser.add_argument('--primer3_dv_conc', type=float, default=DEFAULT_ARGS.primer3_dv_conc, help=argparse.SUPPRESS)
    parser.add_argument('--primer3_dntp_conc', type=float, default=DEFAULT_ARGS.primer3_dntp_conc, help=argparse.SUPPRESS)
    parser.add_argument('--primer3_dna_conc', type=float, default=DEFAULT_ARGS.primer3_dna_conc, help=argparse.SUPPRESS)
    parser.add_argument('--primer3_temp_c', type=float, default=DEFAULT_ARGS.primer3_temp_c, help=argparse.SUPPRESS)
    parser.add_argument('--primer3_max_loop', type=int, default=DEFAULT_ARGS.primer3_max_loop, help=argparse.SUPPRESS)

    # isPcr parameters
    parser.add_argument('--isPcr_minGood', type=int, default=DEFAULT_ARGS.ispcr_min_good, help=argparse.SUPPRESS)
    parser.add_argument('--isPcr_minPerfect', type=int, default=DEFAULT_ARGS.ispcr_min_perfect, help=argparse.SUPPRESS)
    parser.add_argument('--isPcr_tileSize', type=int, default=DEFAULT_ARGS.ispcr_tile_size, help=argparse.SUPPRESS)

    # additional parameters
    parser.add_argument('--temp_tolerance', type=float, default=DEFAULT_ARGS.temp_tolerance, help=argparse.SUPPRESS)
    parser.add_argument('--max_repeats', type=int, default=DEFAULT_ARGS.max_repeats, help=argparse.SUPPRESS)
    parser.add_argument('--bin_size', type=int, default=DEFAULT_ARGS.bin_size, help=argparse.SUPPRESS)
    
    return parser


def parse_args() -> Parameters:
    """Parse commandline arguments and construct a Parameters object from parsed args.

    Returns
        a loaded Parameters object
    """
    parser = create_parser()

    # If help or advanced, print the help menu with/without the advanced options (as appropriate) and exit
    if '--help' in sys.argv or '-h' in sys.argv or '--advanced' in sys.argv:
        print()
        parser.print_help()
        print()
        print(f"If you use this software, please cite our article:")
        print(f"    https://doi.org/10.21105/joss.06850")
        print()
        print(f"{', '.join(__author__)}; 2026")
        print()
        sys.exit(0)
    # If version, print the primerForge version and exit
    elif '--version' in sys.argv or '-v' in sys.argv:
        print(f"primerForge v{__version__}")
        sys.exit(0)
    # If check install, check that all dependencies are available and exit
    elif '--check_install' in sys.argv:
        _check_installation()
        sys.exit(0)
    # Otherwise, parse the arguments normally
    else:
        args = parser.parse_args()

    # Enforce connection between linked args
    if not args.bad_sizes:
        args.bad_sizes = args.pcr_prod
    if args.debug:
        args.keep = True    # keep intermediate files in debug mode

    # Use MAX_CPUS if threads < 1
    if args.num_threads < 1:
        args.num_threads = MAX_CPUS

    # Validate range arguments
    # primer length, pcr product, and bad sizes all tolerate single-values
    args.min_len, args.max_len = [
        int(x) for x in _validate_min_max_ranges("primer_len", args.primer_len, must_have_two=False)]
    args.min_pcr, args.max_pcr = [
        int(x) for x in _validate_min_max_ranges("pcr_prod", args.pcr_prod, must_have_two=False)]
    args.min_bad, args.max_bad = [
        int(x) for x in _validate_min_max_ranges("bad_sizes", args.bad_sizes, must_have_two=False)]

    # GC content and TM range always require two values
    args.min_gc, args.max_gc = _validate_min_max_ranges("gc_range", args.gc_range, must_have_two=True)
    args.min_tm, args.max_tm = _validate_min_max_ranges("tm_range", args.tm_range, must_have_two=True)

    return Parameters(args, initializeLog=True)
