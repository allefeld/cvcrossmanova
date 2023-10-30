---
title: "Tutorial: Other data"
subtitle: Example analysis of data of Kappenman et al. (2021)
listing:
  contents: "/*.qmd"
  type: default
  fields: [title, subtitle]
order: 20
---


As an example for how to use the toolbox to analyze data other than fMRI, here we apply Cross-validated MANOVA as well as Cross-validated Cross-MANOVA to data from @kappenman2021. Their [ERP CORE resource](https://osf.io/thsqg/) provides experimental material and data from paradigms eliciting seven standard ERP components.

In particular, we use the data of subject 1 from their P3 paradigm. In the Supplementary Materials the authors write:

>   The P3 was elicited in an active visual oddball paradigm. […] On each trial, one of five uppercase letters (A, B, C, D, or E […]) was presented for 200 ms in the center of the screen over the fixation point. […] Participants completed a total of 200 trials, divided into five blocks of 40 trials each. In each block, one letter was designated the target stimulus and the other four letters were designated non-targets. Participants pressed one button for targets and another button for non-targets. Each of the five letters served as a target in one block of the experiment and as a non-target in the other four blocks, with the order of blocks randomized across participants. Each letter was presented with equal probability within a block of trials (p = .2), such that the target category was rare (p = .2) and the non-target category was frequent (p = .8).

The analysis is broken down into several steps in the following subpages: