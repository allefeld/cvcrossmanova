function cvCrossManovaSearchlight(dirName, slRadius, analyses, permute, lambda)

% cross-validated (cross-) MANOVA on searchlight
%
% cvCrossManovaSearchlight(dirName, slRadius, Cs, permute = false, lambda = 0)
%
% dirName:   directory where the SPM.mat file referring to an estimated
%            model is located
% slRadius:  radius of the searchlight sphere in voxels
% analyses:  cell array of analysis specifications
% permute:   whether to compute permutation values
% lambda:    regularization parameter (0–1)
%
% Output files are written to the same directory:
% spmD_C####_P####.nii:   images of the pattern discriminability D
%                         contrast and permutation are identified by numbers
% spmDs_C####_P####.nii:  images of standardized pattern discriminability D_s
% VPSL.nii:               image of the number of voxels for each searchlight
% cmsParameters.mat:      record of the analysis parameters


fprintf('\n\ncvCrossManovaSearchlight\n\n')

if nargin < 4
    permute = false;
end
if nargin < 5
    lambda = 0;
end

if dirName(end) ~= filesep
    dirName = [dirName filesep];
end

% simplify things by changing to directory
wd = cd(dirName);
% ensure change back on exit
cleanupObj = onCleanup(@() cd(wd));     

% load data, design matrix etc.
[Ys, Xs, mask, misc] = loadDataSPM(dirName);

% for checkpointing, compute unique ID that encodes parameters:
%   SPM.mat & referenced data -> <timestamp of SPM.mat>,
%   slRadius, analyses, permute, lambda
% encode as string
if ~exist('gencode.m', 'file')
    addpath([spm('dir') filesep 'matlabbatch'])
end
uid = gencode({['SPM.mat of ' getfield(dir('SPM.mat'), 'date')], ...
    slRadius, analyses, permute, lambda});
uid = sprintf('%s\n', uid{:});
% compute Fletcher-16 checksum
uid = dec2hex(fletcher16(uid), 4);

% run searchlight
fprintf('\ncomputing cross-validated MANOVA on searchlight\n')
fEMin = sum(misc.fE) - max(misc.fE);
pMax = slSize(slRadius);
if pMax > fEMin * 0.9       % ensures halfways decent numerical precision
    error('data insufficient for searchlight of size %d!', pMax);
end
fprintf(' running searchlight\n')
fprintf('  searchlight size: %d\n', pMax)
cm = CvCrossManova(Ys, Xs, analyses, fs=misc.fE, lambda=lambda);
[D, p] = runSearchlight(['cmsCheckpoint' uid '.mat'], slRadius, mask, cm);
clear Xs Ys cm    % clear memory
 
% separate contrast and permutation dimensions of result
nContrasts = numel(analyses);
nPerms = size(D, 2) / nContrasts;
D = reshape(D, [], nContrasts, nPerms);

% save results
for ci = 1 : nContrasts
    for pi = 1 : nPerms
        spmWriteImage(reshape(D(:, ci, pi), size(mask)), ...
            sprintf('spmD_C%04d_P%04d.nii', ci, pi), ...
            misc.mat, 'descrip', 'pattern discriminability')
        spmWriteImage(reshape(D(:, ci, pi) ./ sqrt(p), size(mask)), ...
            sprintf('spmDs_C%04d_P%04d.nii', ci, pi), ...
            misc.mat, 'descrip', 'standardized pattern discriminability')
    end
end

% save voxels per searchlight as image
spmWriteImage(reshape(p, size(mask)), 'VPSL.nii', misc.mat, ...
    'descrip', 'voxels per searchlight')

% save analysis parameters
save cmsParameters.mat slRadius analyses permute misc nPerms


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
