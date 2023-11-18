---
title: CvCrossManova Toolbox
subtitle: MVPA by Cross-Validated (Cross-) MANOVA
listing:
  id: subsections
  contents:
    - installation.md
    - tutorial-fmri-spm/index.md
    - tutorial-other/index.md
    - reference.qmd
  type: default
  fields: [title, subtitle]
order: 0
---

::: callout-warning
This project is under development.
:::

CvCrossManova is a toolbox for [Matlab](https://www.mathworks.com/products/matlab.html) which implements the Cross-validated MANOVA and Cross-validated Cross-MANOVA algorithms.

Cross-validated MANOVA was introduced by @allefeld2014 as a method for multivariate pattern analysis (MVPA) in fMRI. Unlike most MVPA methods, it is not based on the accuracy of a classification algorithm, but estimates the multivariate variance explained by a contrast applied to a general linear model (GLM), called pattern distinctness *D*. Using the GLM has the advantage of a common basis with univariate fMRI analyses, resulting in an interpretable measure of effect size, and being applicable to designs where there are no classes. Other than standard GLM-based methods however, it uses cross-validation to obtain an unbiased estimator of that explained multivariate variance.

A limitation of Cross-validated MANOVA was that it provided no analog of cross-classification, i.e. training a classifier on one pair of classes and testing it on another pair of classes, though the original paper pointed to the possibility to use interaction contrasts for a similar purpose. In recent papers, we used an adapted version of Cross-validated MANOVA which uses two different contrasts for 'training' and 'validation', called Cross-validated Cross-MANOVA, but it was published only in the methods sections of these papers. A forthcoming methodological paper will remedy that, and the CvCrossManova toolbox is the corresponding published implementation. Cross-validated Cross-MANOVA estimates the explained multivariate variance shared between two contrasts, called pattern stability *D*^×^.

An implementation of Cross-validated MANOVA was provided previously by the [CvManova toolbox](https://github.com/allefeld/cvmanova). Because Cross-validated MANOVA is a special case of Cross-validated Cross-MANOVA, the CvCrossManova toolbox implements both and is therefore an extended revision of the CvManova toolbox.

These pages document how to use the CvCrossManova toolbox to estimate both pattern distinctness and pattern stability, including examples of analysing fMRI data supported by SPM and other data (in this case EEG).


::: {#subsections}
:::

***

The CvCrossManova toolbox is copyrighted © 2023 by [Carsten Allefeld](https://allefeld.github.io/) and released under the terms of the [GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.en.html) or later. [Remi Gau](https://remi-gau.github.io/), [Polina Iamshchinina](https://www.timbuschman.com/LabMembers/Polina-Iamshchinina), and [Thomas Christophel](https://discolab.eu/team/thomas-christophel/) contributed to earlier versions of the code.

This documentation was created with [Quarto](https://quarto.org/) using a slightly modified [Cosmo Bootswatch theme](https://bootswatch.com/cosmo/).
It is set in [Source Sans 3](https://fonts.google.com/specimen/Source+Sans+3) and [Source Code Pro](https://fonts.google.com/specimen/Source+Code+Pro). Matlab code was executed using [MKernel](https://github.com/allefeld/mkernel).
