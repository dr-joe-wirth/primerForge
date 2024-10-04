from Bio import SeqIO
from bin.Clock import Clock
from typing import Generator
from bin.Primer import Primer
from bin.Product import Product
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters
import multiprocessing, os, subprocess


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


def __makeQueryString(pairs:list[tuple[Primer,Primer]]) -> str:
    """makes the query string from a list of primer pairs

    Args:
        pairs (list[tuple[Primer,Primer]]): a list of primer pairs
    
    Returns:
        str: the query string
    """
    # initialize the output
    out = ""
    
    # for each pair
    for pair in pairs:
        # recast the pair as a strings
        fwd,rev = map(str, pair)
        
        # the "key" column will be the pair separated by a pipe
        key = __KEY_SEP.join((fwd,rev))
        
        # save the data for each pair; one per line
        out += __ROW_SEP.join((key,fwd,rev)) + "\n"
    
    return out


def __runOnePcr(allContigsFna:str, query:str, minGood:int, minPerfect:int, tileSize:int) -> dict[tuple[str,str],set[tuple[str,int]]]:
    """runs isPcr on a properly formatted query string

    Args:
        allContigsFna (str): the filename containing all the contigs
        query (str): the query to search
        minGood (int): the isPcr '-minGood' parameter
        minPerfect (int): the isPcr '-minPerfect' parameter
        tileSize (int): the isPcr '-tileSize' parameter

    Returns:
        dict[tuple[str,str],set[tuple[str,int]]]: key=primer pair (as strings); val=(contig, pcr product size)
    """
    # constants
    CMD = 'isPcr'
    ARGS = ['stdin', 'stdout', '-out=bed']#, '-minGood=6', '-minPerfect=8']
    GOOD = '-minGood='
    PERF = '-minPerfect='
    TILE = '-tileSize='
    EMPTY = ['']
    
    # initialize output
    out = dict()
    
    # build command
    cmd = [CMD, allContigsFna]
    cmd.extend(ARGS)
    cmd.append(GOOD + str(minGood))
    cmd.append(PERF + str(minPerfect))
    cmd.append(TILE + str(tileSize))
    
    # run the command; inject the query string using a pipe
    results = subprocess.run(cmd, input=query.encode(), capture_output=True, check=True)
    
    # convert the results to a collection of rows (lists of strings)
    results = map(lambda x: x.split(__ROW_SEP), results.stdout.decode().split('\n'))
    
    # do not keep empty rows
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


def __runAllPcrs(params:Parameters, pairs:list[tuple[Primer,Primer]]) -> dict[tuple[str,str],set[tuple[str,int]]]:
    """runs all pcrs in parallel

    Args:
        params (Parameters): a Parameters object
        pairs (list[tuple[Primer,Primer]]): a list of Primer pairs

    Returns:
        dict[tuple[str,str],set[tuple[str,int]]]: key=Primer pair; val=(contig, pcr product size)
    """
    # helper function
    def genArgs() -> Generator[tuple[str,str],None,None]:
        """generates arguments for __runOnePcr
        """
        # constant
        MAX_CHUNK = 2000
        
        # initialize values for iteration
        numPairs = len(pairs)
        chunk = min((MAX_CHUNK, numPairs // params.numThreads))
        start = 0
        end = chunk
        
        # keep iterating until all the pairs have been evaluated
        while start <= numPairs:
            # make the query for this chunk
            query = __makeQueryString(pairs[start:end])
            
            # update the start and end for the next chunk
            start = end
            end += chunk
            
            yield (params.allContigsFna, query, params.isPcr_minGood, params.isPcr_minPerfect, params.isPcr_tileSize)
    
    # run pcrs in parallel
    pool = multiprocessing.Pool(params.numThreads)
    results = pool.starmap(__runOnePcr, genArgs())
    pool.close()
    pool.join()
    
    # combine the results and return them
    return {k:v for r in results for k,v in r.items()}


def __filterPairs(pairs:dict[tuple[Primer,Primer],dict[str,Product]], pcrs:dict[tuple[str,str],list[tuple[str,int]]]) -> None:
    """filters primer pairs based on isPcr results; does not return

    Args:
        pairs (dict[tuple[Primer,Primer],dict[str,Product]]): key=Primer pair; dict:key=genome name; val=Product
        pcrs (dict[tuple[str,str],list[tuple[str,int]]]): the dictionary produced by __runAllPcrs
    """
    # for each pair
    for pair,observed in pcrs.items():
        # get a set of the expected product sizes (contig, size)
        expected = {(x.contig, x.size) for x in pairs[pair].values() if x.contig != 'NA'}

        # remove any pairs that do not match the observed product sizes
        if expected != observed:
            del pairs[pair]


def _validatePrimerPairs(params:Parameters, pairs:dict[tuple[Primer,Primer],dict[str,Product]]) -> None:
    """uses in silico PCR to validate primer pairs

    Args:
        params (Parameters): a Parameters object
        pairs (dict[tuple[Primer,Primer],dict[str,Product]]): key=Primer pair; dict:key=genome name; val=Product

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
    
    # create the fasta file
    __makeFastaFile(params)
    
    # run isPcr in parallel
    pcrs = __runAllPcrs(params, list(pairs.keys()))
    
    # filter out bad pairs
    __filterPairs(pairs, pcrs)
    
    # print and log status
    clock.printDone()
    params.log.info(f"done {clock.getTimeString()}")
    
    # make sure that some primer pairs still remain
    if len(pairs) == 0:
        params.log.error(ERR_MSG)
        raise RuntimeError(ERR_MSG)
