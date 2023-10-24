function [Ds, ps] = cvCrossManovaRegion(modelDir, regions, analyses, lambda)

% cross-validated (Cross-) MANOVA on regions
%
% [Ds, ps] = cvCrossManovaRegion(modelDir, regions, analyses, lambda = 1e-8)
%
% ----------  -----------------------------------------------------------------------------
% `modelDir`  directory where the `SPM.mat` file referring to an estimated model is located
% `regions`   region mask(s), cell array of logical 3D volumes or filenames
% `analyses`  cell array of analysis specifications
% `lambda`    regularization parameter (0–1)
% `Ds`        pattern distinctness or pattern stability, cell array analyses × regions
% `ps`        number of voxels in the regions
% ----------  -----------------------------------------------------------------------------
%
% Regarding the parameters `analyses` and `lambda` and the contents of the cells of `Ds`,
% see also Analysis, CvCrossManova.CvCrossManova, CvCrossManova.runAnalyses.

% TODO varargin for any keyword parameters, to be passed to CvCrossManova
% TODO support `loadDataSPM` option `whitenfilter`?

fprintf('\ncvCrossManovaRegion\n\n')

if ~exist('lambda', 'var')
    lambda = 1e-8;
end

% load data, design matrix etc.
[Ys, Xs, ~, misc] = loadDataSPM(modelDir, regions);
nRegions = numel(misc.rmvi);

% compute on regions
fprintf('\ncomputing Cross-validated (Cross-) MANOVA\n')
cm = CvCrossManova(Ys, Xs, analyses, fs=misc.fs, lambda=lambda);
disp(cm)
clear Ys Xs     % release memory
Ds = cell(numel(cm.nResults), nRegions);
for i = 1 : nRegions
    Ds(:, i) = cm.runAnalyses(misc.rmvi{i});
end
clear cm        % release memory

% determine number of voxels per region
ps = cellfun(@numel, misc.rmvi) .';

% Copyright © 2015–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
