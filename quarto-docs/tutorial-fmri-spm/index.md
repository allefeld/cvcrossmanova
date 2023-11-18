---
title: "Tutorial: fMRI with SPM"
subtitle: Example analysis of data of Haxby et al. (2001)
listing:
  contents: "/*.qmd"
  type: default
  fields: [title, subtitle]
order: 10
---


As an example for how to use the toolbox to analyze fMRI data, here we apply Cross-validated MANOVA as well as Cross-validated Cross-MANOVA to the data of subject 1 from @haxby2001.

The details of the experiment and analysis are presented in footnotes. From Footnote 18:

>   Stimuli were gray-scale images of faces, houses, cats, bottles, scissors, shoes, chairs, and nonsense patterns. The categories were chosen so that all stimuli from a given category would have the same base level name. […]\
>   Twelve time series were obtained in each subject. Each time series began and ended with 12 s of rest and contained eight stimulus blocks of 24-s duration, one for each category, separated by 12-s intervals of rest. Stimuli were presented for 500 ms with an interstimulus interval of 1500 ms. […]\
>   To determine the patterns of response to each category on even-numbered and odd-numbered runs, we used a 16-regressor model—eight regressors to model the response to each category relative to rest on even runs and eight regressors to model the response to each category on odd runs with no regressor that contrasted all stimulus blocks to rest.

From Footnote 19:

>   Analysis of the accuracy with which the category being viewed could be identified focused on comparisons between patterns of response for pairs of categories […]\
>   If the within-category correlation (for example, response to category A on even and odd runs) was larger than the between-category correlation (correlation of the response to category A on even runs with the response to category B on odd runs), that comparison was counted as a correct identiﬁcation.

The main result of the paper is the identification accuracy for each category, determined in several different regions of interest and presented in Table 1.

Our analysis is broken down into several steps in the following subpages:


<!-- Copyright © 2023 Carsten Allefeld
SPDX-License-Identifier: GPL-3.0-or-later -->
