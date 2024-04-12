---
title: 'primerForge: a Python package for identifying primer pairs capable of distinguishing groups of genomes from each other'
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
    affiliation: 1
  - name: Lee S. Katz
    orcid: 0000-0002-2533-9161
    affiliation: 1
  - name: Grant Williams
    orcid: 0000-0002-6033-485X
    affiliation: 1
affiliations:
 - name: Centers for Disease Control and Prevention, Atlanta, GA, USA
   index: 1
 - name: Oak Ridge Institute for Science and Education, Oak Ridge, TN, USA
   index: 2
date: 15 April 2024
bibliography: paper.bib
---
# Summary
In both molecular epidemiology and microbial ecology, it is useful to be able to confirm
or deny the presence of specific strains of microorganisms. While whole genome sequencing
and metagenomics can be employed to do this, these techniques are often slow and expensive.
Alternatively, polymerase chain reaction (PCR) can be used to amplify regions of genetic
material that are specific to the strain(s) of  interest. Because PCR is both faster and
less expensive, having a PCR-based approach can accelerate the detection of specific strain(s)
of microbes which faciliates diagnoses and/or population studies.

In order to perform PCR, a pair of DNA primers capable of amplifying a region of interest is
required. Traditional primer design involves the selection of a target region of DNA to amplify, 
followed by primer pair selection and subsequent validation of the primer pair. Identifying a
good pair of primers often requires several iterations of the primer design process, and can be
tedious and time consuming.

`primerForge` seeks to assist biologists with the process of primer design. Instead of requiring
the identification of specific target sequences (as is required with existing tools), `primerForge`
identifies all suitable pairs of primers capable of producing PCR products in a set of genomes.
Optionally, it can also filter those primer pairs and limit the output to pairs tha can be used
to distinguish one set of genomes from another set of genomes.

# Statement of need
# Acknowledgements
# References