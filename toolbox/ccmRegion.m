function [Ds, ps] = ccmRegion(modelDir, regions, analyses, nvargs_md, nvargs_ccm)

% run Cross-validated (Cross-) MANOVA on regions
%
% [Ds, ps] = ccmRegion(modelDir, regions, analyses, wf = true, lambda = 1e-8)
%
% `modelDir` is the directory where the `SPM.mat` file referring to an
% estimated model is located.
%
% `regions` is a cell array of logical region masks specified as
% three-dimensional arrays or filenames.
% 
% `analyses` is a cell array of `Analysis` objects specifying analyses.
%
% The optional `wf` specifies whether to apply the whitening and high-pass
% filtering set up in SPM to data and design matrices. It should usually be
% kept at its default value.
%
% The optional `lambda` (from 0 to 1) controls the amount of shrinkage
% regularization applied to the estimate of the error covariance matrix. It
% should usually be kept at its default value.
%
% `Ds` is a two-dimensional cell array of analysis results. Each cell
% corresponds to the combination of an analysis (rows) and a region
% (columns) and contains either a scalar value or an array of permutation
% values. Whether a result is an estimate of pattern distinctness *D* or
% pattern stability *D*^×^ depends on the contrasts of the corresponding
% analysis and the regressors involved in them.
%
% `ps` is an array with the numbers of voxels contained in each region.


arguments
    modelDir           (1, :)  char
    regions            (1, :)  cell
    analyses           (1, :)  cell
    nvargs_md.wf       (1, 1)  logical
    nvargs_ccm.lambda  (1, 1)  double
end

fprintf('\nccmRegion\n\n')

% create ModeledData object encapsulating data and design matrices
nvargs_md = namedargs2cell(nvargs_md);      % pass on wf
[modeledData, misc] = ModeledData.fromSPM(modelDir, nvargs_md{:}, ...
    regions = regions);
nRegions = numel(misc.rmvi);

% create CvCrossManova object implementing
% Cross-validated (Cross-) MANOVA analysis
nvargs_ccm = namedargs2cell(nvargs_ccm);    % pass on lambda
ccm = CvCrossManova(modeledData, analyses, nvargs_ccm{:});

% compute on regions
fprintf('\ncomputing Cross-validated (Cross-) MANOVA\n')
disp(ccm)
Ds = cell(ccm.nAnalyses, nRegions);
for i = 1 : nRegions
    Ds(:, i) = ccm.runAnalyses(misc.rmvi{i});
end
clear ccm        % release memory

% determine number of voxels per region
ps = cellfun(@numel, misc.rmvi) .';


% Copyright © 2015–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
