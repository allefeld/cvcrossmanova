function [D, p] = cvCrossManovaRegion(dirName, regions, analyses, permute, lambda)

% cross-validated MANOVA on region
%
% [D, p] = cvCrossManovaRegion(dirName, regions, analyses, permute = false, lambda = 0)
%
% ----------  ---------------------------------------------------------------------------
% `dirName`   directory where the SPM.mat file referring to an estimated model is located
% `regions`   region mask(s), (cell array of) logical 3D volume(s) or filename(s)
% `analyses`  cell array of analysis specifications
% `lambda`    regularization parameter (0–1)
% `permute`   whether to compute permutation values
% `D`         pattern distinctness, contrasts × permutations × regions
% `p`         number of voxels in the region(s)
% ----------  -----------------------------------------------------------------


fprintf('\n\ncvCrossManovaRegion\n\n')

if nargin < 4
    permute = false;
end
if nargin < 5
    lambda = 0;
end    
assert(~islogical(lambda), 'incorrect order of parameters?')
assert(~isempty(regions), 'no region mask specified!')

if dirName(end) ~= filesep
    dirName = [dirName filesep];
end

% load data, design matrix etc.
[Ys, Xs, ~, misc] = loadDataSPM(dirName, regions);
nRegions = numel(misc.rmvi);

% compute on regions
fprintf('\ncomputing cross-validated MANOVA on region\n')
cm = CvCrossManova(Ys, Xs, analyses, fs=misc.fE, lambda=lambda);
clear Ys Xs
fEMin = sum(misc.fE) - max(misc.fE);
p = cellfun(@numel, misc.rmvi);
D = [];
for i = 1 : nRegions
    s = warning('off', 'MATLAB:nearlySingularMatrix');
    D = cat(3, D, cm.runAnalyses(misc.rmvi{i}));
    warning(s.state, 'MATLAB:nearlySingularMatrix')
    if p(i) > fEMin * 0.9   % ensures decent numerical precision
        D(:, :, i) = nan;
        fprintf(2, 'data insufficient for the %d voxels of region %d!\n', p(i), i);
    end
end
clear cm          % clear memory

% separate contrast and permutation dimensions of result
D = reshape(D, numel(analyses), [], nRegions);


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.