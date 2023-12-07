function ccmSearchlight(modelDir, radius, analyses, ...
    nvargs, nvargs_md, nvargs_ccm)

% run Cross-validated (Cross-) MANOVA on searchlight
%
% ccmSearchlight(modelDir, radius, analyses, ...
%     mmUnits = false, wf = true, lambda = 1e-8)
%
% `modelDir` is the directory where the `SPM.mat` file referring to an
% estimated model is located.
%
% `radius` is the radius of the searchlight sphere.
%
% `analyses` is a cell array of `Analysis` objects specifying analyses.
%
% The optional `mmUnits` specifies whether the searchlight radius is in
% voxels or millimeters. By default the value of `radius` is interpreted as
% voxels, while interpretation as mm has to be requested by `mmUnits =
% true`.
%
% The optional `wf` specifies whether to apply the whitening and high-pass
% filtering set up in SPM to data and design matrices. It should usually
% be kept at its default value.
%
% The optional `lambda` (from 0 to 1) controls the amount of shrinkage
% regularization applied to the estimate of the error covariance matrix. It
% should usually be kept at its default value.
%
% Analysis results are written to image files in the same directory
% `modelDir`, with names of the form `spmD_A####_P####.nii`. The number
% following `A` indicates the analysis and the number following `P`
% indicates the permutation. Whether a result is an estimate of pattern
% distinctness *D* or pattern stability *D*^×^ depends on the contrasts of
% the corresponding analysis and the regressors involved in them.
% 
% In addition, an image file `VPSL.nii` is written containing the numbers
% of voxels for each searchlight, and a Matlab file
% `ccmSearchlightParams.mat` is written containing a record of the analysis
% parameters.
%
% The searchlight procedure includes a checkpointing mechanism:
% Intermediate results are saved to a file `ccmsCheckpoint….mat`, where `…`
% stands for a 32-digit hexadecimal checksum of the analysis parameters. If
% a run of `ccmSearchlight` is interrupted, running it again recovers
% partial results from the checkpoint file and continues from there. Note
% that checkpointing only works if *all* analysis parameters remain
% identical. In particular, if the analysis includes randomly selected
% permutations, make sure to initialize Matlab's random number generator
% before the selection.

arguments
    modelDir           (1, :)  char
    radius             (1, 1)  double
    analyses           (1, :)  cell
    nvargs.mmUnits     (1, 1)  logical = false
    nvargs_md.wf       (1, 1)  logical
    nvargs_ccm.lambda  (1, 1)  double
end

fprintf('\nccmSearchlight\n\n')

% create ModeledData object encapsulating data and design matrices
nvargs_md = namedargs2cell(nvargs_md);      % pass on wf
[modeledData, misc] = ModeledData.fromSPM(modelDir, nvargs_md{:});

% create CvCrossManova object implementing
% Cross-validated (Cross-) MANOVA analysis
nvargs_ccm = namedargs2cell(nvargs_ccm);    % pass on lambda
ccm = CvCrossManova(modeledData, analyses, nvargs_ccm{:});

% for checkpointing, compute unique filename which encodes parameters:
%   `modelDir`, `radius`, `analyses`, `wf`, `lambda`
%   + SPM.mat & referenced data -> <timestamp of SPM.mat>
params = struct;
params.date = getfield(dir(fullfile(modelDir, 'SPM.mat')), 'date');
params.modelDir = modelDir;
params.radius = radius;
params.analyses = analyses;
params.mmUnits = nvargs.mmUnits;
params.wf = misc.wf;
params.lambda = ccm.lambda;
% encode parameters as string
if ~exist('gencode.m', 'file')
    addpath(fullfile(spm('dir'), 'matlabbatch'))
end
uid = gencode(params);
uid = sprintf('%s\n', uid{:});
% compute MD5 checksum as hexadecimal string
md5Instance = java.security.MessageDigest.getInstance('MD5');
uid =  md5Instance.digest(unicode2native(uid, 'UTF-8'));
uid = reshape(dec2hex(uid).', 1, 32);
% checkpoint filename
checkpoint = fullfile(modelDir, ['ccmSearchlightCheckpoint' uid '.mat']);

% get prototype searchlight
if nvargs.mmUnits
    PSL = searchlight(radius, misc.mat);
    unit = 'mm';
else
    PSL = searchlight(radius);
    unit = 'voxels';
end

% run searchlight
fprintf('\ncomputing Cross-validated (Cross-) MANOVA\n')
disp(ccm)
fprintf('  running searchlight of radius %d %s (%d voxels)\n', ...
    radius, unit, nnz(PSL))
fprintf('  intermediate results are saved to\n')
fprintf('      %s\n', checkpoint)
[res, ps] = runSearchlight(checkpoint, PSL, misc.mask, ccm);
nResults = ccm.nResults;
clear ccm        % release memory

% save results
for i = 1 : numel(nResults)
    for j = 1 : nResults(i)
        volume = nan(size(misc.mask));
        volume(misc.mask) = res{i}(:, j);
        fn = fullfile(modelDir, sprintf('spmD_A%04d_P%04d.nii', i, j));
        spmWriteImage(volume, fn, misc.mat, ...
            'descrip', 'pattern distinctness / stability')
    end
end

% save voxels per searchlight as image
volume = nan(size(misc.mask));
volume(misc.mask) = ps;
fn = fullfile(modelDir, 'VPSL.nii');
spmWriteImage(volume, fn, misc.mat, 'descrip', 'voxels per searchlight')

% save analysis parameters
save(fullfile(modelDir, 'ccmSearchlightParams.mat'), '-struct', 'params')


% Copyright © 2013–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
