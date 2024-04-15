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
  - name: Jessica C. Chen
    orcid: 0000-0002-9320-6774
    affiliation: "1"
  - name: Lee S. Katz
    orcid: 0000-0002-2533-9161
    affiliation: "1"
  - name: Grant M. Williams
    orcid: 0000-0002-6033-485X
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
to categorize specific strains of microorganisms in either an ingroup or an outgroup in a
given population. While whole genome sequencing and pangenomics pipelines can be
employed to do this, these techniques are often slow and can be expensive.
Examples of these pangenomics pipelines include PIRATE and Rory [@pirate] [@rory].
<!-- I think it might be more appropriate to have pangenomics instead metagenomics as the contrast ~~~~Lee -->
The conclusion of these pangenomics pipelines would tell the laboratory that,
in order to find whether or not a new sample is part of the ingroup or outgroup,
they would have to sequence the whole genome.
Alternatively, polymerase chain reaction (PCR) can be used to amplify regions
of genetic material that are specific to the strain(s) of  interest. PCR is
faster, less expensive, and more accessible than whole genome sequencing, so
having a PCR-based approach can accelerate the detection of specific strain(s)
of microbes and faciliate diagnoses and/or population studies.

# Statement of need
<!-- I feel like you might want this explainer either here or in the previous section ~~~~Lee -->
There are many use cases for what `primerForge` offers.
One use case would be surveillance of an outbreak clone of a particular pathogen.
A laboratory could develop a set of PCR reactions to track the population
of this outbreak clone which could help inform whether the population
were to grow, shrink, or migrate.
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
the output to primer pairs tha can be used to distinguish one set of genomes
from another set of genomes via PCR amplification.

# Acknowledgements
This work was made possible by support from the Advanced Molecular Detection
initiative at the Centers for Disease Control and Prevention and is covered by
activities approved by the Centers for Disease Control and Prevention Internal
Review Board (approval no. 7172).

# References
