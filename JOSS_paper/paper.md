---
title: 'osTIR: An entirely open source package for estimating translation initiation.'
tags:
  - Python
  - open source
  - ribosome
  - rbs
  - translation
  - calculator
authors:
  - name: Cameron Roots
    affiliation: "1"
  - name: Alexandra Lukasiewicz
    affiliation: "2"
  - name: Jeffrey Barrick
    affiliation: "1"
affiliations:
  - Department of Molecular Biosciences, University of Texas at Austin
    index: 1
  - McKetta Department of Chemical Engineering, University of Texas at Austin
    index: 2
date: 3 December 2020
bibliography: paper.bib
---

# Summary

Transcriptional initiation is one of the fundamental steps of the central dogma,
defining the rate at which proteins are generated from an mRNA template. The
ability to accurately predict the relative activity of a ribosome binding site
is an important part of metabolic engineering studies, enabling researchers to
rationally design initiation sites and identify positions where unintentional
initiation may occur. While prior works have produced software packages with
this goal in mind, they either rely on closed source dependencies or are closed
source themselves.

# Statement of Need

`osTIR` (open source Translation Initiation of Ribosomes) is an open source
Python package that integrates RNA spacing considerations and the
ViennaRNA [@Lorenz:2011] software suite to perform free energy calculations
for binding events involved during protein translation initiation. This work is
derived from the last open source version of the calculator described in Salis
2009, and accuracy of osTIR was benchmarked against the same experimental data
with comparable or better accuracy compared to previous works
[@Salis:2009; @Na:2010; @Reis:2020; @Bonde:2016].

`osTIR` was written with a user-friendly interface and high throughput in mind.
The most recent open source version of the calculator this work derives from
requires NUPack [@Zadeh:2011], which has licensing restrictions, and any updates
beyond the initial 1.0 version has been restricted to use through a webserver,
restricting the ease of use within workflows. Other related software face
similar restrictions [@Na:2010]. This work is installable by the user through
Bioconda [@Zadeh:2018], enabling easy integration. `osTIR` supports multiFASTA
input with command line parameters or CSV input with support for defining
parameters on a per-sample basis. Additionally, `osTIR` supports asynchronous
processing, accelerating use cases that require screening of many mRNAs.  

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

# Acknowledgements

???

# References
