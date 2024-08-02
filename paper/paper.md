---
title: 'primerForge: a Python program for identifying primer pairs capable of distinguishing groups of genomes from each other'
tags:
  - Python
  - microbial ecology
  - microbiology
  - molecular biology
  - molecular epidemiology
  - PCR
  - polymerase chain reaction
  - primers
authors:
  - name: Joseph S. Wirth
    orcid: 0000-0002-9750-2845
    affiliation: "1,2"
  - name: Lee S. Katz
    orcid: 0000-0002-2533-9161
    affiliation: "1"
  - name: Grant M. Williams
    orcid: 0000-0002-6033-485X
    affiliation: "1"
  - name: Jessica C. Chen
    orcid: 0000-0002-9320-6774
    affiliation: "1"
affiliations:
 - name: Centers for Disease Control and Prevention, Atlanta, GA, USA
   index: "1"
 - name: Oak Ridge Institute for Science and Education, Oak Ridge, TN, USA
   index: "2"
date: 15 April 2024
bibliography: paper.bib
---
# Summary
In both molecular epidemiology and microbial ecology, it is useful to be able
to categorize specific strains of microorganisms in either an ingroup or an
outgroup in a given population. While whole genome sequencing and downstream
phylogenetic analyses can be employed to do this, these techniques are often
slow and can be resource intensive. Additionally, the laboratory would have to
sequence the whole genome to use these tools to determine whether or not a new
sample is part of the ingroup or outgroup. Alternatively, polymerase chain
reaction (PCR) can be used to amplify regions of genetic material that are
specific to the strain(s) of  interest. PCR is faster, less expensive, and more
accessible than whole genome sequencing, so having a PCR-based approach can
accelerate the detection of specific strain(s) of microbes and facilitate
diagnoses and/or population studies.

# Statement of need
In order to perform PCR, a pair of DNA primers capable of amplifying a region
of interest is required. Traditional primer design involves the selection of a
target region of DNA to amplify, followed by primer pair selection and
subsequent validation of the primer pair. Identifying a good pair of primers
and a suitable target region often requires several iterations of the primer
design process, which can be tedious and time consuming.

`primerForge` seeks to assist biologists with the process of primer design.
Instead of requiring the identification of specific target sequences (as is
required with existing tools), `primerForge` identifies all suitable pairs of
primers capable of producing PCR products of a specific size in a set of whole
genome sequences. Optionally, it can also filter those primer pairs and limit
the output to primer pairs that can be used to distinguish one set of genomes
from another set of genomes via PCR amplification. `primerForge` relies on the
`khmer` package to extract k-mers from genomic sequences and the `primer3-py`
package to evaluate specific characteristics of primer pairs including melting
temperature, hairpin potential, and dimer formation
[@10.12688/f1000research.6924.1; @10.1093/nar/gks596]. It also uses _in silico_
PCR via the `isPcr` program to validate and filter the primer pairs
[@10.1093/bib/bbs038].

There are many use cases for what `primerForge` offers. One use case would be
surveillance of an outbreak clone of a particular pathogen. A laboratory could
develop a set of PCR reactions to track the population of this outbreak clone
which could help inform if the population were to grow, shrink, or migrate.

# Comparing `primerForge` to `swga2`
Another software package, `swga2`, was developed to choose sets of primers that
selectively amplify one set of genomes but not in another set of genomes
[@10.1371/journal.pcbi.1010137]. Given this similarity, it was important to
compare the performance of `primerForge` to that of `swga2`. To do this, the
sequence files listed in \autoref{tab:sequences} were used as inputs. The
default parameters were used except that the melting temperature range was set
at 55°C to 68°C and only a single processor was used. The primer pairs
identified by the `swga2` and `primerForge` were evaluated using `isPcr` with
the following additional parameters: `-tileSize=8`, `-minGood=8`, and
`-minPerfect=8`. These parameters were necessary because the primers identified
by `swga2` were too short to use the default value of 15. The results of these
comparisons are shown in \autoref{tab:comparisons}.

Table: \label{tab:sequences} Datasets used to compare `primerForge` to `swga2`.

|Dataset ^[1]|Name|NCBI Accession|Group|
|:------|:--:|:------------:|:---:|
|plasmid|pcDNA|not provided|ingroup|
|plasmid|pLTR|not provided|outgroup|
|_M. mycoides_|_Mycoplasma mycoides_ subsp. mycoides str. KH3J|GCF_003034305.1|ingroup|
|_M. mycoides_|_Mycoplasma mycoides_ subsp. mycoides str. B345/93|GCF_003034275.1|ingroup|
|_M. mycoides_|_Mycoplasma mycoides_ subsp. mycoides str. Gemu Goffa|GCF_003034345.1|ingroup|
|_M. mycoides_|_Mycoplasma mycoides_ subsp. capri str. GM12|GCF_900489555.1|outgroup|
|_M. mycoides_|_Mycoplasma mycoides_ subsp. capri str. 80/93|GCF_018389745.1|outgroup|
|_E. coli_|_Escherichia coli_ O157 str. 644-PT8|GCF_001650295.1|ingroup|
|_E. coli_|_Escherichia coli_ O157 str. AR-0428|GCF_008727175.1|ingroup|
|_E. coli_|_Escherichia coli_ O157 str. FDAARGOS_293|GCF_002208865.2|ingroup|
|_E. coli_|_Escherichia coli_ K12 str. MG1655|GCF_000005845|outgroup|
|_E. coli_|_Salmonella enterica_ subsp. enterica serovar Typhimurium str. LT2|GCF_000006945|outgroup|
|SARS-CoV-2|SARS-CoV-2 isolate human/USA/MA_MGH_00230/2020|MT520374|ingroup|
|SARS-CoV-2|SARS-CoV-2 isolate human/USA/MA_MGH_00229/2020|MT520263|ingroup|
|SARS-CoV-2|SARS-CoV-2 isolate human/USA/MA_MGH_00257/2020|MT520479|ingroup|
|SARS-CoV-2|SARS-CoV-2 isolate Wuhan-Hu-1|NC_045512|outgroup|

[^1]: The plasmid dataset is provided as an example in the `swga2` repository. The _M. mycoides_ is provided as an example in the `primerForge` repository.

Table: \label{tab:comparisons} Comparing `swga2` to `primerForge`

|program|dataset|runtime (mm:ss)|RAM (Gb)|primer pairs ^[2]|`isPcr`-compatible pairs ^[3]|validated pairs ^[4]|optimized pairs ^[5]|
|:-----:|:-----:|:-------------:|:------:|:-----------:|:-----------------------:|:--------------:|:--------------:|
|`swga2`|plasmid|23:21|0.136|94|22|22|11|
|`primerForge`|plasmid|00:10|0.051|3,210|3,210|3,168|2,934|
|`swga2`|_M. mycoides_|05:13|0.221|run failed|NA|NA|NA|
|`primerForge`|_M. mycoides_|02:52|1.478|1446|1,446|989|884|
|`swga2`|_E. coli_|21:10|4.452|run failed|NA|NA|NA|
|`primerForge`|_E. coli_|83:33|10.329|1,451,164|1,451,164|318,927|125,932|
|`swga2`|SARS-CoV-2|10:38|0.141|63|7|0|0|
|`primerForge`|SARS-CoV-2|00:19|0.122|39|39|15|15|

^[2]: the number of primer pairs identified by the program

^[3]: the number of primer pairs that generated PCR products with
`isPcr`

^[4]: the number of primer pairs that produced a PCR product in every
ingroup genome and no products in any of the outgroup genomes

[^5]: the number of valid primer pairs that produced exactly one PCR
product in each ingroup genome

Although many of the primer pairs predicted by `primerForge` were not validated
by `isPcr`, this can be attributed to the fact that `primerForge` allows primer
pairs to produce PCR products in the outgroup provided they are outside of the
user-specified range. For example, 1,118,742 of the 1,132,237 primer pairs
(98.8 %) identified by `primerForge` in the _E. coli_ dataset that were not
validated by `isPcr` were predicted to produce a PCR product in one or more of
the outgroup genomes. Similarly, 453 of the 457 primer pairs identified in the
_M. mycoides_ dataset were not validated with `isPcr` for the same reason. The
remaining invalidated primer pairs can be attributed to the fact that the
`isPcr` parameters `tileSize`, `minGood`, and `minPerfect` were set to very low
values in order to directly compare the results of `primerForge` with those
produced by `swga2`.

In addition to the improvements observed in \autoref{tab:comparisons},
`primerForge` also allows the user to specify PCR product size ranges, specify
primer size ranges, specify melting temperature ranges, specify the difference
in melting temperature between the forward and reverse primers, exclude PCR
product sizes from the outgroup genomes, and use sequences in genbank file
format. Both `primerForge` and `swga2` can take advantage of parallel
processing, which speeds up runtimes. For example, `primerForge` ran on the
"`primerForge` _Escherichia coli_ example" dataset in 37 minutes and 7 seconds
when using 8 cpus.

# Acknowledgements
This work was made possible by support from the Advanced Molecular Detection
initiative at the Centers for Disease Control and Prevention and is covered by
activities approved by the Centers for Disease Control and Prevention Internal
Review Board (approval no. 7172).

# References
