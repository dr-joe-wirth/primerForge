from Bio.Seq import Seq
from bin.Clock import Clock
from bin.Primer import Primer
import multiprocessing, primer3
from Bio.SeqRecord import SeqRecord
from bin.Parameters import Parameters

# functions
def __getUniqueKmers(seqs:list[SeqRecord], minLen:int, maxLen:int, name:str) -> dict[str,dict[Seq,dict[str,tuple[str,int,int,str]]]]:
    """gets a collection of kmers from a genome that are:
        * not repeated anywhere in the genome
        * at least one end is a G or C

    Args:
        seqs (list[SeqRecord]): a list of contigs as SeqRecord objects
        minLen (int): the minimum kmer length
        maxLen (int): the maximum kmer length
        name (str): the name of the genome
    
    Returns:
        dict[str,dict[Seq,dict[str,tuple[str,int,int,str]]]]: key=strand; val=dict: key=kmer seq; val=dict: key=name; val=tuple: contig, start, klen, strand
    """
    # constant
    GC = {"G", "C", 'g', 'c'}
    
    # helper functions to clarify boolean expressions
    def isOneEndGc(seq:Seq) -> bool:
        """ evaluates if one end of a sequence is a GC
        """
        return seq[-1] in GC or seq[0] in GC

    def isDuplicated(seq:Seq) -> bool:
        return seq in kmers[Primer.PLUS].keys() or seq in kmers[Primer.MINUS].keys()

    # initialize variables
    kmers = dict()
    kmers[Primer.PLUS] = dict()
    kmers[Primer.MINUS] = dict()
    bad = set()
    krange = range(minLen, maxLen + 1)
    smallest = min(krange)
    
    # go through each contig
    for contig in seqs:
        # extract the contig sequences
        fwdSeq:Seq = contig.seq
        revSeq:Seq = fwdSeq.reverse_complement()
        
        # get the length of the contig
        contigLen = len(contig)
        done = False
        
        # get every possible kmer start position
        for start in range(contigLen):
            # for each allowed kmer length
            for klen in krange:
                # stop looping through the contig once we're past the smallest kmer
                if start+smallest > contigLen:
                    done = True
                    break
                
                # stop looping through the kmers once the length is too long
                elif start+klen > contigLen:
                    break
            
                # proceed if the extracted kmer length is good
                else:
                    # extract the kmer sequences
                    fwdKmer = fwdSeq[start:start+klen]
                    
                    # handle the first base differently for the (-) strand (special slicing condition)
                    if start == 0:
                        revKmer = revSeq[-(start+klen):]
                    else:
                        revKmer = revSeq[-(start+klen):-start]
                    
                    # mark duplicate kmers for removal
                    if isDuplicated(fwdKmer):
                        bad.add(fwdKmer)
                        bad.add(revKmer)
                    
                    # only save kmers that have GC at one end
                    elif isOneEndGc(fwdKmer):
                        kmers[Primer.PLUS][fwdKmer]  = {name: (contig.id, start, klen, Primer.PLUS)}
                        kmers[Primer.MINUS][revKmer] = {name: (contig.id, start, klen, Primer.MINUS)}
            
            # stop iterating through the contig when we're done with it
            if done:
                break

    # discard any sequences that were duplicated
    for seq in bad:
        try: kmers[Primer.PLUS].pop(seq)
        except KeyError: pass
        
        try: kmers[Primer.MINUS].pop(seq)
        except KeyError: pass
    
    return kmers
            

def __getSharedKmers(seqs:dict[str,list[SeqRecord]], params:Parameters) -> dict[Seq,dict[str,tuple[str,int,int,str]]]:
    """retrieves all the kmers that are shared between the input genomes

    Args:
        seqs (dict[str,list[SeqRecord]]): key=genome name; val=list of contigs
        params (Parameters): a Parameters object

    Raises:
        RuntimeError: could not extract kmers from a single genome
        RuntimeError: a single genome has no shared kmers with other genomes

    Returns:
        dict[Seq,dict[str,tuple[str,int,int,str]]]: key=kmer; val=dict: key=genome name: val=tuple: contig, start, klen, strand
    """
    # messages
    ERR_MSG_1 = 'failed to extract kmers from '
    ERR_MSG_2 = 'no shared kmers after processing '
    DBG_MSG_1 = ' shared kmers after processing '
    
    # intialize variables
    sharedKmers = dict()
    firstGenome = True

    # set the debugger if required
    params.log.rename(__getSharedKmers.__name__)

    # for each genome
    for name in seqs.keys():
        # get the unique kmers for this genome/kmer lengths
        kmers = __getUniqueKmers(seqs[name], params.minLen, params.maxLen, name)
        
        # make sure there are kmers for this genome
        if kmers == dict():
            params.log.error(f"{ERR_MSG_1}{name}")
            raise RuntimeError(f"{ERR_MSG_1}{name}")
        
        # keep only the (+) strand kmers if this is the first genome
        if firstGenome:
            sharedKmers.update(kmers[Primer.PLUS])
            firstGenome = False
        
        # otherwise keep the shared kmers (both + and -)
        else:
            # get all the kmers for this genome (both + and - strands)
            thisGenome = set(kmers[Primer.PLUS].keys()).union(kmers[Primer.MINUS].keys())
            
            # go through each kmer that is currently shared by all processed genomes
            shared = set(sharedKmers.keys())
            for kmer in shared:
                # save data for kmers shared with this genome
                if kmer in thisGenome:
                    # strand information is now lost; try both
                    try:
                        sharedKmers[kmer].update(kmers[Primer.PLUS][kmer])
                    except KeyError:
                        sharedKmers[kmer].update(kmers[Primer.MINUS][kmer])
                
                # remove any kmers that are not shared with this genome
                else:
                    sharedKmers.pop(kmer)
        
        # make sure the kmers did not depopulate
        if sharedKmers == dict():
            params.log.error(f"{ERR_MSG_2}{name}")
            raise RuntimeError(f"{ERR_MSG_2}{name}")
        
        # log the number of shared kmers if debugging
        params.log.debug(f"{' '*8}{len(sharedKmers)}{DBG_MSG_1}{name}")
    
    return sharedKmers


def __reorganizeDataByPosition(name:str, kmers:dict[Seq,dict[str,tuple[str,int,int,str]]]) -> dict[str,dict[int,list[tuple[Seq,int,str]]]]:
    """reorganizes data from __getSharedKmers by its genomic position

    Args:
        name (str): the name of the genome to be processed
        kmers (dict[Seq,dict[str,tuple[str,int,int]]]): the dictionary produced by __getSharedKmers

    Returns:
        dict[str,dict[int,list[tuple[Seq,int,str]]]]: key=contig; val=dict: key=start position; val=list of tuples: kmer, klen, strand
    """
    # initialize output
    out = dict()
    
    # for each kmer
    for kmer in kmers.keys():
        # extract data from the dictionary
        contig, start, length, strand = kmers[kmer][name]
        
        # contig = top level key; start position = second level key; val = list
        out[contig] = out.get(contig, dict())
        out[contig][start] = out[contig].get(start, list())
        
        # add the sequence and its length to the list
        out[contig][start].append((kmer, length, strand))
    
    return out


def __evaluateKmersAtOnePosition(contig:str, start:int, posL:list[tuple[Seq,int,str]], minGc:float, maxGc:float, minTm:float, maxTm:float) -> Primer:
    """evaluates all the primers at a single position in the genome; designed for parallel calls

    Args:
        contig (str): the name of the contig
        start (int): the start position in the sequence
        posL (list[tuple[Seq,int]]): a list of tuples: kmer, klen, strand
        minGc (float): the minimum percent GC allowed
        maxGc (float): the maximum percent GC allowed
        minTm (float): the minimum melting temperature allowed
        maxTm (float): the maximum melting temperature allowed
    
    Returns:
        Primer: a suitable primer at the given position
    """
    # constant
    FIVE_DEGREES = 5
    
    # define helper functions to make booleans below more readable
    def isGcWithinRange(primer:Primer) -> bool:
        """is the percent GC within the acceptable range?"""
        return primer.gcPer >= minGc and primer.gcPer <= maxGc

    def isTmWithinRange(primer:Primer) -> bool:
        """is the Tm within the acceptable range?"""
        return primer.Tm >= minTm and primer.Tm <= maxTm
    
    def noLongRepeats(primer:Primer) -> bool:
        """verifies that a primer does not have long repeats
        """
        # constants
        MAX_LEN = 4
        REPEATS = ("A"*MAX_LEN, "T"*MAX_LEN, "C"*MAX_LEN, "G"*MAX_LEN)

        # check for each repeat in the primer
        for repeat in REPEATS:
            if repeat in primer.seq:
                return False
        return True

    def noHairpins(primer:Primer) -> bool:
        """verifies that the primer does not form hairpins
        """
        # hairpin tm should be less than (minTm - 5°); need to check both strands
        fwdOk = primer3.calc_hairpin_tm(str(primer)) < (minTm - FIVE_DEGREES)
        revOk = primer3.calc_hairpin_tm(str(primer.reverseComplement())) < (minTm - FIVE_DEGREES)
        
        return fwdOk and revOk

    def noHomodimers(primer:Primer) -> bool:
        """verifies that the primer does not form homodimers
        """
        # homodimer tm should be less than (minTm - 5°); need to check both strands
        fwdOk = primer3.calc_homodimer_tm(str(primer)) < (minTm - FIVE_DEGREES)
        revOk = primer3.calc_homodimer_tm(str(primer.reverseComplement())) < (minTm - FIVE_DEGREES)
        
        return fwdOk and revOk
    
    # initialize values for the while loop
    idx = 0
    
    # continue to iterate through each primer in the list until a primer is found 
    for idx in range(len(posL)):
        # extract data from the list
        seq,length,strand = posL[idx]
        
        # create a Primer object
        primer = Primer(seq, contig, start, length, strand)
        
        # evaluate the primer's percent GC, Tm, hairpin potential, and homodimer potential; save if passes
        if isGcWithinRange(primer) and isTmWithinRange(primer): # O(1)
            if noLongRepeats(primer): # O(1)
                if noHairpins(primer):
                    if noHomodimers(primer):
                        return primer


def __evaluateAllKmers(kmers:dict[str,dict[int,list[tuple[Seq,int,str]]]], minGc:float, maxGc:float, minTm:float, maxTm:float, numThreads:int) -> list[Primer]:
    """evaluates kmers at each position for their suitability as primers

    Args:
        kmers (dict[str,dict[int,list[tuple[Seq,int]]]]): the dictionary produced by __reorganizeDataByPosition
        minGc (float): the minimum percent G+C
        maxGc (float): the maximum percent G+C
        minTm (float): the minimum melting temperature
        maxTm (float): the maximum melting temperature
        numThreads (int): number of threads for parallel processing

    Returns:
        list[Primer]: a list of suitable primers as Primer objects
    """
    # initialize a list of arguments
    args = list()

    # each contig needs to be evalutated
    for contig in kmers.keys():
        # each start position within the contig needs to be evaluated
        for start in kmers[contig].keys():
            # save arguments to pass in parallel
            args.append((contig, start, kmers[contig][start], minGc, maxGc, minTm, maxTm))

    # parallelize primer evaluations
    pool = multiprocessing.Pool(processes=numThreads)
    results = pool.starmap(__evaluateKmersAtOnePosition, args)
    pool.close()
    pool.join()

    # remove failed searches before returning
    return [x for x in results if x is not None]
 

def __buildOutput(kmers:dict[str,dict[Seq,tuple[str,int,int,str]]], candidates:list[Primer]) -> dict[str,dict[str,list[Primer]]]:
    """builds the candidate primer output

    Args:
        kmers (dict[str,dict[Seq,tuple[str,int,int,str]]]): the dictionary produced by __getSharedKmers
        candidates (list[Primer]): the list produced by __evaluateAllKmers

    Returns:
        dict[str,dict[str,list[Primer]]]: key=genome name; val=dict: key=contig; val=list of Primers
    """
    # initialize output
    out = dict()
    
    # for each candidate primer
    for cand in candidates:
        # identify which sequence is present in the genome
        try:
            entry = kmers[cand.seq]
        except:
            entry = kmers[cand.seq.reverse_complement()]
        
        # for each genome
        for name in entry.keys():
            # extract the data from the entry
            contig, start, length, strand = entry[name]
        
            # save a Primer object in this contig's list
            out[name] = out.get(name, dict())
            out[name][contig] = out[name].get(contig, list())
            out[name][contig].append(Primer(cand.seq, contig, start, length, strand))
    
    return out


def __getCandidatesForOneGenome(name:str, kmers:dict[Seq,dict[str,tuple[str,int,int]]], params:Parameters) -> dict[str,dict[str,list[Primer]]]:
    """gets candidate primers for a single genome

    Args:
        name (str): the name of the genome to get primers for
        kmers (dict[Seq,dict[str,tuple[str,int,int]]]): the dictionary produced by __getSharedKmers
        params (Parameters): a Parameters object

    Returns:
        dict[str.dict[str,list[Primer]]]: key=genome name; val=dict: key=contig; val=list of Primers
    """
    # reorganize data by each unique start positions for one genome
    positions = __reorganizeDataByPosition(name, kmers)
    
    # get a list of the kmers that pass the evaulation
    candidates = __evaluateAllKmers(positions, params.minGc, params.maxGc, params.minTm, params.maxTm, params.numThreads)
    
    # create a dictionary whose keys are contigs and values are lists of candidate primers
    return __buildOutput(kmers, candidates)


def _getAllCandidateKmers(ingroup:dict[str,list[SeqRecord]], params:Parameters, sharedExists:bool) -> dict[str,dict[str,list[Primer]]]:
    """gets all the candidate kmer sequences for a given ingroup

    Args:
        ingroup (dict[str,list[SeqRecord]]): key=genome name; val=contigs as a list of SeqRecord (generator actually)
        params (Parameters): a Parameters object
        sharedExists (bool): indicates if the shared kmers are already pickled

    Raises:
        RuntimeError: no shared ingroup kmers
        RuntimeError: shared ingroup kmers are not suitable for use as primers

    Returns:
        dict[str,dict[str,list[Primer]]]: key=genome name; val=dict: key=contig; val=list of Primers
    """
    # constants
    SHARED_KMER_NUM = 0
    
    # messages
    GAP = " "*4
    MSG_1 = f"{GAP}getting shared ingroup kmers that appear once in each genome"
    MSG_2 = f"{GAP}evaluating kmers"
    MSG_3A = f"{GAP}identified "
    MSG_3B = " candidate kmers"
    ERR_MSG_1 = "failed to identify a set of kmers shared between the ingroup genomes"
    ERR_MSG_2 = "none of the ingroup kmers are suitable for use as a primer"
    
    # initialize variables
    clock = Clock()
    numCand = 0
    out = dict()
    
    # setup debugger
    params.log.rename(_getAllCandidateKmers.__name__)
    
    if sharedExists:
        # close the genome files (not actually a list of seqRecords; they are generators)
        for handle in ingroup.values():
            handle.stream.close()
        
        # load existing shared kmers from file
        kmers = params.loadObj(params.pickles[SHARED_KMER_NUM])

    else:
        # get all non-duplicated kmers that are shared in the ingroup
        params.log.info(MSG_1)
        clock.printStart(MSG_1)
        kmers = __getSharedKmers(ingroup, params)
        clock.printDone()
        
        # move log back to this function
        params.log.rename(_getAllCandidateKmers.__name__)
        params.log.info(f'{GAP}done {clock.getTimeString()}')
    
        # make sure that ingroup kmers were identified
        if kmers == dict():
            params.log.error(ERR_MSG_1)
            raise RuntimeError(ERR_MSG_1)
        
        # dump the shared kmers to file
        params.dumpObj(kmers, params.pickles[SHARED_KMER_NUM], "shared kmers")
    
    # print status
    clock.printStart(MSG_2)
    params.log.info(MSG_2)
    
    # go through each genome name
    names = list(next(iter(kmers.values())).keys())
    for name in names:
        # get the candidate kmers for the genome
        candidates = __getCandidatesForOneGenome(name, kmers, params)
        
        # for each genome in the candidates
        for genome in candidates.keys():
            # create a sub dictionary if one does not already exist
            out[genome] = out.get(genome, dict())
            
            # for each contig in the genome
            for contig in candidates[genome].keys():
                # store a set of all the candidates for this contig
                out[genome][contig] = out[genome].get(contig, set())
                out[genome][contig].update(candidates[genome][contig])
    
    # for each genome
    for name in out.keys():
        # for each contig
        for contig in out[name].keys():
            # sort kmers by their start position on the (+) strand
            out[name][contig] = sorted(out[name][contig], key=lambda x: min(x.start, x.end))
            
            # count the number of candidate kmers for the first genome only
            if name == names[0]:
                numCand += len(out[name][contig])

    # make sure candidates were found
    if numCand == 0:
        params.log.error(ERR_MSG_2)
        raise RuntimeError(ERR_MSG_2)
    
    # print status
    clock.printDone()
    print(f"{MSG_3A}{numCand}{MSG_3B}")
    
    # log status
    params.log.info(f'{GAP}done {clock.getTimeString()}')
    params.log.info(f'{MSG_3A}{numCand}{MSG_3B}')

    return out
