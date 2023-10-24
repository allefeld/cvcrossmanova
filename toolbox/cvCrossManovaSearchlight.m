function cvCrossManovaSearchlight(modelDir, slRadius, analyses, lambda)

% cross-validated (Cross-) MANOVA on searchlight
%
% cvCrossManovaSearchlight(modelDir, slRadius, analyses, lambda = 1e-8)
%
% ----------  -----------------------------------------------------------------------------
% `modelDir`  directory where the `SPM.mat` file referring to an estimated model is located
% `slRadius`  radius of the searchlight sphere in voxels
% `analyses`  cell array of analysis specifications
% `lambda`    regularization parameter (0–1)
% ----------  -----------------------------------------------------------------------------
%
% Results are written to files the same directory:
%
% -----------------------  -------------------------------------------------------------------------------------------------------
% `spmD_A####_P####.nii`   images of pattern distinctness or pattern stability; analysis and permutation are identified by numbers
% `VPSL.nii`               image of the number of voxels for each searchlight
% `ccmsParameters.mat`     record of the analysis parameters
% -----------------------  -------------------------------------------------------------------------------------------------------
%
% The searchlight procedure includes a checkpointing mechanism:
% Intermediate results are saved to a file `ccmsCheckpoint….mat`, where `…`
% stands for a 32-digit hexadecimal checksum. If an analysis is
% interrupted, running it again recovers partial results from the
% checkpoint file and continues from there.
%
% Note that this only works if the analysis parameters remain identical. In
% particular, if the analysis includes randomly selected permutations, make
% sure to initialize Matlab's random number generator before the selection. It
% is recommended to run `s = rng('shuffle')` once, note the values of `s.Seed`
% and `s.Type`, and then to include `rng(<seed>, <type>)` in your analysis
% pipeline.
%
% Regarding the parameters `slRadius`, `analyses` and `lambda`
% see also slSize, Analysis, CvCrossManova.CvCrossManova.

% TODO varargin for any keyword parameters, to be passed to CvCrossManova
% TODO support `loadDataSPM` option `whitenfilter`?

fprintf('\ncvCrossManovaSearchlight\n\n')

if ~exist('lambda', 'var')
    lambda = 1e-8;
end

% load data, design matrix etc.
[Ys, Xs, mask, misc] = loadDataSPM(modelDir);

% for checkpointing, compute unique filename which encodes parameters:
%   SPM.mat & referenced data -> <timestamp of SPM.mat> & `modelDir`,
%   `slRadius`, `analyses`, `lambda`
% encode parameters as string
if ~exist('gencode.m', 'file')
    addpath(fullfile(spm('dir'), 'matlabbatch'))
end
uid = gencode({getfield(dir(fullfile(modelDir, 'SPM.mat')), 'date'), ...
    modelDir, slRadius, analyses, lambda});
uid = sprintf('%s\n', uid{:});
% compute MD5 checksum as hexadecimal string
uid = java.security.MessageDigest.getInstance('MD5').digest(unicode2native(uid, 'UTF-8'));
uid = reshape(dec2hex(uid).', 1, 32);
% checkpoint filename
checkpoint = fullfile(modelDir, ['ccmsCheckpoint' uid '.mat']);

% run searchlight
fprintf('\ncomputing Cross-validated (Cross-) MANOVA\n')
cm = CvCrossManova(Ys, Xs, analyses, fs=misc.fs, lambda=lambda);
disp(cm)
clear Ys Xs     % release memory
fprintf('  running searchlight of radius %d (%d voxels)\n', slRadius, slSize(slRadius))
fprintf('  intermediate results are saved to\n')
fprintf('      %s\n', checkpoint)
[res, ps] = runSearchlight(checkpoint, slRadius, mask, cm);
nResults = cm.nResults;
nAnalyses = numel(nResults);
clear cm        % release memory

% save results
for i = 1 : nAnalyses
    for j = 1 : nResults(i)
        volume = reshape(res{i}(:, j), size(mask));
        fn = fullfile(modelDir, sprintf('spmD_A%04d_P%04d.nii', i, j));
        spmWriteImage(volume, fn, misc.mat, ...
            'descrip', 'pattern distinctness / stability')
    end
end

% save voxels per searchlight as image
spmWriteImage(reshape(ps, size(mask)), fullfile(modelDir, 'VPSL.nii'), ...
    misc.mat, 'descrip', 'voxels per searchlight')

% save analysis parameters
save(fullfile(modelDir, 'ccmsParameters.mat'), ...
    'slRadius', 'analyses', 'lambda', 'misc')
% parameter `modelDir` is implicit in the location of the file

% Copyright © 2013–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
