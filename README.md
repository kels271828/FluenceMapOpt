# A nonconvex optimization approach to IMRT planning with dose-volume constraints

Fluence map optimization for intensity-modulated radiation therapy planning can be formulated as a large-scale inverse problem with multi-objectives on the tumors and organs-at-risk. Unfortunately, clinically relevant dose-volume constraints are nonconvex, so convex formulations and algorithms cannot be directly applied to the problem. We propose a novel approach to handle dose-volume constraints while preserving their nonconvexity, as opposed to previous efforts which focused on iterative convexification. The proposed method is amenable to efficient algorithms based on partial minimization and naturally adapts to handle maximum and mean dose constraints, which are prevalent in current practice, and cases of infeasibility. We demonstrate our approach using the [CORT dataset](https://doi.org/10.1186/2047-217X-3-37), and show that it is easily adaptable to radiation treatment planning with dose-volume constraints for multiple tumors and organs-at-risk.

## Documents
* [Preprint](https://arxiv.org/abs/1907.10712) (July 2019)
* [Poster](https://github.com/kels271828/FluenceMapOpt/blob/master/poster.pdf) (April 2018)

## Code
Download data and solver from the links below, and unzip in the same directory as code.
* [Main](https://github.com/kels271828/FluenceMapOpt/blob/master/FluenceMapOpt.m): Functions for loading data, computing fluence maps, and plotting results
* [Script](https://github.com/kels271828/FluenceMapOpt/blob/master/run.m): Example script for using FluenceMapOpt
* [Examples](https://github.com/kels271828/FluenceMapOpt/tree/master/Examples): Code to reproduce the examples from our paper

## Links
* [CORT data](https://gigadb.org/dataset/100110)
  * Prostate case data: ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100110/PROSTATE.zip
  * Prostate DICOM data: ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100110/Prostate_Dicom.zip
* [minConf solver](https://www.cs.ubc.ca/~schmidtm/Software/minConf.zip)
