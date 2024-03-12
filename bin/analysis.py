import numpy, os
from Bio import SeqIO
from Bio.Seq import Seq
from bin.Clock import Clock
from matplotlib import pyplot
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from scipy.stats import gaussian_kde
from bin.Parameters import Parameters
from bin.AnalysisData import AnalysisData, Level
from matplotlib.backends.backend_pdf import PdfPages


def _initializeAnalysisData(candKmers:dict[str,dict[str,list[Primer]]]) -> dict[tuple[Seq,str],AnalysisData]:
    """initializes a dictionary of analysis objects

    Args:
        candKmers (dict[str,dict[str,list[Primer]]]): key = genome name; val=dict: key=contig name; val=list of candidate kmers

    Returns:
        dict[tuple[Seq,str],AnalysisData]: key=(kmer sequence, contig name); val=AnalysisData object
    """
    # initialize variables
    index = 0
    out = dict()
    
    # for each genome
    for name in candKmers.keys():
        # for each contig
        for contig in candKmers[name].keys():
            # for each kmer
            for kmer in candKmers[name][contig]:
                # save the kmer analysis data
                out[(kmer.seq,contig)] = AnalysisData(kmer, index, name)
                
                # update the index
                index += 1
    
    return out


def _updateAnalysisData(analysis:dict[tuple[Seq,str],AnalysisData], pairs:list[tuple[Primer,Primer]]) -> None:
    """updates a dictionary of AnalysisData objects with pair and data level info

    Args:
        analysis (dict[tuple[Seq,str],AnalysisData]): the dictionary produced by _initializeAnalysisData
        pairs (list[tuple[Primer,Primer]]): a list of primer pairs
    """
    # initialize variables
    seen = set()
    allKmers = {k for k,c in analysis.keys()}
    allContigs = {c for k,c in analysis.keys()}
    
    # for each pair of primers
    for fwd,rev in pairs:
        # make sure all primers match the existing entries
        if not fwd.seq in allKmers:
            fwd = fwd.reverseComplement()
        if not rev.seq in allKmers:
            rev = rev.reverseComplement()
        
        # store the contigs for this pair of primers
        contigsThisPair = set()
        
        # for each contig
        for contig in allContigs:
            try:
                # extract the indices of the pair
                fIdx = analysis[(fwd.seq, contig)].getIndex()
                rIdx = analysis[(rev.seq, contig)].getIndex()

                # link the records by their indices
                analysis[(fwd.seq, contig)].updatePairs(rIdx)
                analysis[(rev.seq, contig)].updatePairs(fIdx)
                
                # keep the contigs for this pair
                contigsThisPair.add(contig)
            
            # only one contig per genome will be represented
            except KeyError:
                pass
        
        # only increment forward once
        if fwd not in seen:
            seen.add(fwd)
            for contig in contigsThisPair:
                analysis[(fwd.seq, contig)].incrementLevel() 
            
        # only increment reverse once
        if rev not in seen:
            seen.add(rev)
            for contig in contigsThisPair:
                analysis[(rev.seq, contig)].incrementLevel()


def __initializeCounts(files:list[str], frmt:str) -> dict[str,dict[str,dict[int,dict[Level,int]]]]:
    """initializes a dictionary of counts with 0 values

    Args:
        files (list[str]): a list of ingroup files
        frmt (str): the format of the ingroup files

    Returns:
        dict[str,dict[str,dict[int,dict[Level,int]]]]: key=genome name; val=dict: key=contig name; val=dict: key=position; val=dict: key=level; val=count
    """
    # initialize variables
    out = dict()
    contig:SeqRecord

    # for each filename
    for fn in files:
        # get the genome name
        name = os.path.basename(fn)
        
        # create a subdictionary for each genome
        out[name] = dict()
        
        # for each contig
        for contig in SeqIO.parse(fn, frmt):
            # save a dictionary at each position with 0 counts for each level
            out[name][contig.id] = {pos:{level:0 for level in AnalysisData.LEVELS} for pos in range(len(contig))}
    
    return out
    
    
def __countPositions(entries:dict[tuple[Seq,str],AnalysisData], params:Parameters) -> dict[str,dict[str,dict[int,dict[Level,int]]]]:
    """counts the number of kmers that appear at each genomic position

    Args:
        entries (dict[tuple[Seq,str],AnalysisData]): key=(kmer sequence, contig name); val=AnalysisData object
        params (Parameters): a Parameters object

    Returns:
        dict[str,dict[str,dict[int,dict[str,int]]]]: key=name; val=dict: key=contig name; val=dict: key=position; val=dict: key=level; val=count
    """
    # intitialize output
    out = __initializeCounts(params.ingroupFns, params.format)
          
    # start counting
    for kmer,contig in entries.keys():
        entry = entries[(kmer,contig)]
        
        # for each level
        for level in AnalysisData.LEVELS:
            # only count if the kmer is within this level
            if entry.getLevel() >= level:
                # get the start and end in order
                start = min(entry.primer.start, entry.primer.end)
                end   = max(entry.primer.start, entry.primer.end)
                
                # count each position in the kmer
                for pos in range(start, end+1):
                    out[entry.name][contig][pos][level] += 1
    
    return out


def __concatenateContigCounts(counts:dict[str,dict[str,dict[int,dict[Level,int]]]]) -> tuple[dict[str,dict[Level,dict[int,int]]], dict[str,dict[str,tuple[int,int]]]]:
    """concatenates the counts for each contig for plotting

    Args:
        counts (dict[str,dict[str,dict[int,dict[Level,int]]]]): the dictionary produced by __coutnPositions

    Returns:
        tuple[dict[str,dict[Level,dict[int,int]]], dict[str,dict[str,int]]]: (key=genome name; val=dict: key=level; val=dict: key=position; val=count), (key=genome name; val=dict: key=contig name; val=(shift,end in plot)
    """
    # initialize outputs
    concatCounts = dict()
    contigBreaks = dict()
    
    # for each genome
    for name in counts.keys():
        # reset the shift
        shift = 0
        
        # initialize new subdictionaries
        concatCounts[name] = dict()
        contigBreaks[name] = dict()
        
        # for each contig
        for contig in counts[name].keys():
            # for each position
            for position in counts[name][contig].keys():
                # shift the position
                newPos = shift + position
                
                # for each level
                for level in AnalysisData.LEVELS:
                    # initialize a new subdictionary if one does not exist
                    concatCounts[name][level] = concatCounts[name].get(level, dict())
                    
                    # save the counts under the new position
                    concatCounts[name][level][newPos] = counts[name][contig][position][level]
        
            # keep track of where the contig ended
            contigBreaks[name][contig] = shift,newPos
            
            # update the shift
            shift = newPos + 1
    
    return concatCounts, contigBreaks


def __makeOnePlot(plotData:dict[int,int], boundaries:dict[str,tuple[int,int]], title:str, pdf:PdfPages) -> None:
    """makes a single plot and writes it to file

    Args:
        plotData (dict[int,int]): key=position; val=count
        boundaries (dict[str,int]): key=contig name; val=plot position
        title (str): plot title
        pdf (PdfPages): a filehandle to the pdf
    """
    # constants
    PLOT_WIDTH  = 11
    PLOT_HEIGHT = 8
    LINE_WIDTH  = 0.4
    LINE_STYLE  = '-'
    NUM_SPACE = 120
    RAW = " raw"
    DENSITY = " density"

    # get arrays of the positions and their counts
    positions = numpy.array(list(plotData.keys()))
    counts = numpy.array(list(plotData.values()))
    
    # plot the raw values
    pyplot.figure(figsize=(PLOT_WIDTH,PLOT_HEIGHT))
    pyplot.plot(positions, counts, linestyle=LINE_STYLE, linewidth=LINE_WIDTH)
    
    # plot the edge of each contig
    for shift,edge in boundaries.values():
        pyplot.axvline(x=edge, color='red', linestyle=LINE_STYLE, linewidth=LINE_WIDTH)
   
    # add labels 
    pyplot.xlabel('Position')
    pyplot.ylabel('Count')
    pyplot.title(f"{title}{RAW}")
    
    # write the plot to file
    pdf.savefig()
    pyplot.close()
    
    # calculate the density (can fail with empty data)
    try:
        x = numpy.linspace(min(positions), max(positions), NUM_SPACE)
        kde = gaussian_kde(numpy.repeat(positions, counts))
        y = kde(x)

        # plot the density
        pyplot.figure(figsize=(PLOT_WIDTH,PLOT_HEIGHT))
        pyplot.plot(x, y)
        
        # plot the edge of each contig
        for shift,edge in boundaries.values():
            pyplot.axvline(x=edge, color='red', linestyle=LINE_STYLE, linewidth=LINE_WIDTH)
        
        # add labels
        pyplot.xlabel('Position')
        pyplot.ylabel('Density')
        pyplot.title(f'{title}{DENSITY}')
        
        # write the plot to file
        pdf.savefig()
        pyplot.close()
    
    except:
        pass


def __restructureAnalysisDataForWriting(analysisData:dict[tuple[Seq,str],AnalysisData]) -> dict[tuple[Seq,Level], dict[str,tuple[str,AnalysisData]]]:
    """restructures analysis data in preparation for writing it to file

    Args:
        analysisData (dict[tuple[Seq,str],AnalysisData]): the dictionary produced by __initializeAnalysisData and modified by __updateAnalysisData

    Returns:
        dict[tuple[Seq,Level], dict[str,tuple[str,AnalysisData]]]: key=(kmer, Level); val=dict: key = genome name; val=(contig name, AnalysisData)
    """
    # initialize output
    out = dict()
    
    # for each kmer
    for kmer,contig in analysisData.keys():
        # extract the entry and the name
        entry = analysisData[(kmer,contig)]
        name = entry.name
        
        # create the new output key
        key = (kmer,entry.getLevel())
        
        # build data structure
        out[key] = out.get(key, dict())
        out[key][name] = (contig, entry)
    
    return out


def __writeAnalysisData(analysisData:dict[tuple[Seq,str],AnalysisData], contigBreaks:dict[str,dict[str,tuple[int,int]]], fn:str) -> None:
    """writes analysis data to file

    Args:
        analysisData (dict[tuple[Seq,str],AnalysisData]): the dictionary made and edited by __initializeAnalysisData and __updateAnalysisData, respectively
        contigBreaks (dict[str,dict[str,tuple[int,int]]]): the second dictionary produced by __concatenateContigCounts
        fn (str): the filename to write
    """
    # constants
    EOL = "\n"
    SEP = "\t"
    GENOME_COORD = "genome coord"
    PLOT_COORD = "plot coord"
    NA = "NA"
    
    # build the header
    header = ["kmer", "dataset"]
    names = list(contigBreaks.keys())
    for name in names:
        header.append(f"{name} {GENOME_COORD}")
        header.append(f"{name} {PLOT_COORD}")
    
    # restructure the data for writing
    outData = __restructureAnalysisDataForWriting(analysisData)
    
    with open(fn, 'w') as fh:
        # write the header to file
        fh.write(SEP.join(header) + EOL)
            
        # for each kmer
        for kmer,level in outData.keys():
            # initialize the row
            row = [kmer, level]
            
            # for each genome
            for name in names:
                # extract the contig name, entry, and shift
                try:
                    contig, entry = outData[(kmer,level)][name]
                    shift,edge = contigBreaks[name][contig]
                    
                    # get the genomic and plot coordinates
                    gCoord = f"{contig} ({entry.primer.start}, {entry.primer.end})"
                    pCoord = f"({entry.primer.start + shift}, {entry.primer.end + shift})"
                
                # if the genome doesn't have this kmer (rev-comp) then add filler data
                except KeyError:
                    gCoord = NA
                    pCoord = NA
            
                # save them to the row
                row.append(gCoord)
                row.append(pCoord)
            
            # write the row to file
            fh.write(SEP.join(map(str,row)) + EOL)


def _plotAnalysisData(analysisData:dict[tuple[Seq,str],AnalysisData], params:Parameters) -> None:
    """plots and writes the analysis data

    Args:
        analysisData (dict[tuple[Seq,str],AnalysisData]): the dictionary produced by _initializeAnalysisData and modified by __updateAnalysisData
        params (Parameters): a Parameters object
    """
    # messages
    GAP = 4 * ' '
    DONE = "done"
    MSG_1 = GAP + 'preparing data for plotting'
    MSG_2 = GAP + 'plotting data for '
    MSG_3 = GAP + 'writing raw plot data to file'
    
    # initialize the clock
    clock = Clock()
    
    # print the status and log
    params.log.rename(_plotAnalysisData.__name__)
    params.log.info(MSG_1)
    clock.printStart(MSG_1)
    
    # count the positions
    counts = __countPositions(analysisData, params)
    
    # concatenate contigs
    catCounts, contigBreaks = __concatenateContigCounts(counts)
    
    # print status
    clock.printDone()
    params.log.info(f"{DONE} {clock.getTimeString()}")
    
    # open the pdf
    with PdfPages(params.plotsFn) as pdf:
        # for each genome
        for name in catCounts.keys():
            # print the status
            params.log.info(f"{MSG_2}{name}")
            clock.printStart(f"{MSG_2}{name}")
            
            # for each level
            for level in catCounts[name].keys():
                # make the plots
                __makeOnePlot(catCounts[name][level], contigBreaks[name], f"{name} ({level})", pdf)
            
            # print the status
            clock.printDone()
            params.log.info(f"{DONE} {clock.getTimeString()}")
    
    # print status
    clock.printStart(MSG_3)
    params.log.info(MSG_3)
    
    # write the analysis data
    __writeAnalysisData(analysisData, contigBreaks, params.plotDataFn)

    # print status
    clock.printDone()
    params.log.info(f"{DONE} {clock.getTimeString()}")
