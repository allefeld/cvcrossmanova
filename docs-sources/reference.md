---
title: Reference documentation
author: Carsten Allefeld
date: 2023-3-23
keep-md: true
keep-ipynb: true
---



## [function]{.smallcaps} cvCrossManovaSearchlight

cross-validated (cross-) MANOVA on searchlight

```matlab
cvCrossManovaSearchlight(dirName, slRadius, Cs, permute = false, lambda = 0)
```

----------  ---------------------------------------------------------------------------
`dirName`   directory where the SPM.mat file referring to an estimated model is located
`slRadius`  radius of the searchlight sphere in voxels
`analyses`  cell array of analysis specifications
`permute`   whether to compute permutation values
`lambda`    regularization parameter (0&#x2013;1)
----------  ---------------------------------------------------------------------------

Output files are written to the same directory:

-----------------------  --------------------------------------------------------------------------------------------
`spmD_C####_P####.nii`   images of the pattern discriminability D, contrast and permutation are identified by numbers
`spmDs_C####_P####.nii`  images of standardized pattern discriminability D~s~
`VPSL.nii`               image of the number of voxels for each searchlight
`cmsParameters.mat`      record of the analysis parameters
-----------------------  --------------------------------------------------------------------------------------------



## [function]{.smallcaps} cvCrossManovaRegion

cross-validated MANOVA on region

```matlab
[D, p] = cvCrossManovaRegion(dirName, regions, analyses, permute = false, lambda = 0)
```

----------  ---------------------------------------------------------------------------
`dirName`   directory where the SPM.mat file referring to an estimated model is located
`regions`   region mask(s), (cell array of) logical 3D volume(s) or filename(s)
`analyses`  cell array of analysis specifications
`lambda`    regularization parameter (0&#x2013;1)
`permute`   whether to compute permutation values
`D`         pattern distinctness, contrasts &#xD7; permutations &#xD7; regions
`p`         number of voxels in the region(s)
----------  -----------------------------------------------------------------



## [class]{.smallcaps} Analysis

object representing an analysis

```matlab
Analysis(CA, CB, sessionsA, sessionsB)
```


## [class]{.smallcaps} CvCrossManova

```
Index exceeds the number of array elements. Index must not exceed 1.
Error in help2markdown (line 45)
    syntax = helpText{2};
```


