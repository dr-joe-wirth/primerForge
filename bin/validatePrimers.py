import os, subprocess
from Bio import SeqIO
from bin.Clock import Clock
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters


# global constants
__KEY_SEP = "|"
__ROW_SEP = "\t"


# functions
def __makeFastaFile(params:Parameters) -> None:
    """creates a fasta file containing all the contigs in the data set

    Args:
        params (Parameters): a Parameters object
    """
    # constant
    OUT_FORMAT = "fasta"
    
    # initialize variable
    contig:SeqRecord
    
    # remove the fasta file if it already exists
    if os.path.exists(params.allContigsFna):
        os.remove(params.allContigsFna)
    
    # open the file in append mode
    with open(params.allContigsFna, 'a') as fh:
        # for each sequence file
        for fn in params.ingroupFns + params.outgroupFns:
            # write each contig to file
            for contig in SeqIO.parse(fn, params.format):
                # make sure the sequence name matches what was saved in other data structures
                SeqIO.write(SeqRecord(contig.seq, contig.id, '', ''), fh, OUT_FORMAT)


def __makeQueryFile(params:Parameters, pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]) -> None:
    """makes the query file for 

    Args:
        params (Parameters): a Parameters object
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]): key=Primer pair; dict:key=genome name; val=tuple: contig, pcr product size, bin pair (contig, num1, num2)
    """
    # open the file
    with open(params.queryFn, 'w') as fh:
        for pair in pairs.keys():
            # recast the pair as a strings
            fwd,rev = map(str, pair)
            
            # the "key" column will be the pair separated by a pipe
            key = __KEY_SEP.join((fwd,rev))
            
            # write the data for each pair to file; one per line
            fh.write(__ROW_SEP.join((key,fwd,rev)) + "\n")


def __runPcr(params:Parameters) -> dict[tuple[str,str],set[tuple[str,int]]]:
    """runs isPcr on the primer pairs

    Args:
        params (Parameters): a Parameters object

    Returns:
        dict[tuple[str,str],set[tuple[str,int]]]: key=primer pair (as strings); val=pcr product size
    """
    # constants
    CMD = 'isPcr'
    ARGS = ['stdout', '-out=bed']
    EMPTY = ['']
    
    # initialize output
    out = dict()
    
    # build command
    cmd = [CMD, params.allContigsFna, params.queryFn]
    cmd.extend(ARGS)
    
    # run the command
    results = subprocess.run(cmd, capture_output=True, check=True)
    
    # convert the results to a list of values; do not keep empty lines
    results = list(map(lambda x: x.split(__ROW_SEP), results.stdout.decode().split('\n')))
    results = [x for x in results if x != EMPTY]
    
    # for each result
    for res in results:
        # get the contig id, pcr size, and pair
        contig = res[0]
        size = abs(int(res[1]) - int(res[2]))
        pair = tuple(res[3].split(__KEY_SEP))
        
        # save the results
        out[pair] = out.get(pair, set())
        out[pair].add((contig, size))
    
    return out


def __filterPairs(pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]], pcrs:dict[tuple[str,str],list[tuple[str,int]]]) -> None:
    """filters primer pairs based on isPcr results; does not return

    Args:
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]): key=Primer pair; dict:key=genome name; val=tuple: contig, pcr product size, bin pair (contig, num1, num2)
        pcrs (dict[tuple[str,str],list[tuple[str,int]]]): the dictionary produced by __runPcr
    """
    # for each pair
    for pair,observed in pcrs.items():
        # get a set of the expected product sizes (contig, size)
        expected = {x[:2] for x in pairs[pair].values() if x[0] != 'NA'}

        # remove any pairs that do not match the observed product sizes
        if expected != observed:
            del pairs[pair]


def _validatePrimerPairs(params:Parameters, pairs:dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]) -> None:
    """uses in silico PCR to validate primer pairs

    Args:
        params (Parameters): a Parameters object
        pairs (dict[tuple[Primer,Primer],dict[str,tuple[str,int,tuple[str,int,int]]]]): key=Primer pair; dict:key=genome name; val=tuple: contig, pcr product size, bin pair (contig, num1, num2)

    Raises:
        RuntimeError: no valid primer pairs
    """
    # messages
    MSG_1 = 'validating primer pairs with isPcr'
    ERR_MSG = 'failed to find valid primer pairs'
    
    # initialize clock
    clock = Clock()
    
    # log data
    params.log.rename(_validatePrimerPairs.__name__)
    
    # print and log status
    clock.printStart(MSG_1)
    params.log.info(MSG_1)
    
    # create the fasta and query files and run isPcr
    __makeFastaFile(params)
    __makeQueryFile(params, pairs)
    pcrs = __runPcr(params)
    
    # filter out bad pairs
    __filterPairs(pairs, pcrs)
    
    # print and log status
    clock.printDone()
    params.log.info(f"done {clock.getTimeString()}")
    
    # make sure that some primer pairs still remain
    if len(pairs) == 0:
        params.log.error(ERR_MSG)
        raise RuntimeError(ERR_MSG)
