# A nonconvex optimization approach to IMRT planning with dose-volume constraints

[![DOI](https://zenodo.org/badge/224053659.svg)](https://zenodo.org/badge/latestdoi/224053659)

Fluence map optimization for intensity-modulated radiation therapy planning can be formulated as a large-scale inverse problem with competing objectives and constraints associated with the tumors and organs-at-risk.
Unfortunately, clinically relevant dose-volume constraints are nonconvex, so standard algorithms for convex problems cannot be directly applied.
While prior work focused on convex approximations for these constraints, we propose a novel relaxation approach to handle nonconvex dose-volume constraints.
We develop efficient, provably convergent algorithms based on partial minimization, and show how to adapt them to handle maximum-dose constraints and infeasible problems.
We demonstrate our approach using the CORT dataset, and show that it is easily adaptable to radiation treatment planning with dose-volume constraints for multiple tumors and organs-at-risk.

## Documents
* [Preprint](https://arxiv.org/abs/1907.10712) (August 2021)
* [Poster](https://github.com/kels271828/FluenceMapOpt/blob/master/poster.pdf) (April 2018)

## Code
Download data and solver from the links below, and unzip in the same directory as code.
* [Main](https://github.com/kels271828/FluenceMapOpt/blob/master/FluenceMapOpt.m): Functions for loading data, computing fluence maps, and plotting results
* [Script](https://github.com/kels271828/FluenceMapOpt/blob/master/run.m): Example script for using FluenceMapOpt
* [Examples](https://github.com/kels271828/FluenceMapOpt/tree/master/Examples): Code to reproduce the examples from our paper
* [Figures](https://github.com/kels271828/FluenceMapOpt/tree/master/Figures): Figures from our paper

## Links
* [minConf solver](https://www.cs.ubc.ca/~schmidtm/Software/minConf.zip)
* [CORT dataset](https://gigadb.org/dataset/100110)
  * Prostate case data: ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100110/PROSTATE.zip
  * Prostate DICOM data: ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100110/Prostate_Dicom.zip
