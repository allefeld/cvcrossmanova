% Cross-validated Cross-MANOVA
% 
% converted from quarto-docs/tutorial-fmri-spm/cvcrossmanova.qmd

% ## Heading
% 
% This is an introductory text which explains all the very important things
% that are going to be done in this section.

% select subject
sub = 'subj1';
% load information
load(fullfile(sub, 'info.mat'))
% directory with estimated model
modelDir = fullfile(sub, 'model');

%%

% >   8 contrasts between animate ("face", "cat") and inanimate objects
% ("bottle", "scissors", "shoe", "chair"). Cross-MANOVA for all pairs of
% contrasts.
% 
% 
% <!-- Copyright © 2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later -->

