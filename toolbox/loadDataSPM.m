function [Ys, Xs, fs, names, misc] = loadDataSPM(modelDir, regions, wf)

% load fMRI data via SPM.mat
%
% [Ys, Xs, fs, names, misc] = loadDataSPM(modelDir, regions = {}, wf = true)
%
% modelDir:  name of directory that contains SPM.mat
% regions:   optional additional region mask(s),
%            cell array of logical 3D volume(s) or filename(s)
% wf:        whether to whiten and filter data and design matrices
%
% Ys:        fMRI data matrices, 1 × sessions cell array of arrays scans × voxels
% Xs:        design matrices, 1 × sessions cell array of arrays scans × regressors
% fs:        residual degrees of freedom, 1 × sessions array
% names:     names of regressors, 1 × sessions cell array of string arrays 1 × regressors
% misc:      struct with additional information:
%   mask       analysis brain mask, logical 3D volume;
%              intersected with union of region masks if present
%   mat      3D voxel indices to mm transformation matrix
%   rmvi     1 × regions cell array of indices into columns of Y
%
% Y includes only voxels within the mask, in the linear order of the mask.

% default argument values
if nargin < 2
    regions = {};
end
if nargin < 3
    wf = true;
end

% load SPM.mat
SPMname = fullfile(modelDir, 'SPM.mat');
fprintf('loading data via %s\n', SPMname)
load(SPMname, 'SPM');

% get data volumes and check matching voxel grid
VY = SPM.xY.VY;
spm_check_orientations(VY);

% check whether data might have been moved
if ~exist(VY(1).fname, 'file')
    SPMold = fullfile(SPM.swd, 'SPM.mat');
    comLen = min(numel(SPMname), numel(SPMold));
    comPart = comLen - find(diff(...
        [SPMname(end - comLen + 1 : end) == SPMold(end - comLen + 1 : end) 1] ...
        ), 1, 'last');
    fprintf(2, '  if analysis and data folders were moved together, try\n');
    fprintf(2, '    spm_changepath(''%s'', ''%s'', ''%s'')\n', ...
        modelDir, SPMold(1 : end - comPart), SPMname(1 : end - comPart));
end

% read analysis brain mask image
assert(isfield(SPM, 'VM'), 'no analysis brain mask in SPM.VM!')
try
    mask = (spm_read_vols(SPM.VM) > 0);
catch
    % SPM8 stores the filename without the path
    VM = spm_vol(fullfile(modelDir, SPM.VM.fname));
    mask = (spm_read_vols(VM) > 0);
end
fprintf('  volume of %d × %d × %d = %d voxels\n', size(mask), numel(mask))
fprintf('  %d voxels within brain mask\n', sum(mask(:)));

% possibly apply region mask(s)
if isempty(regions)
    fprintf('  no region masks\n')
    rmvi = {};
else
    nRegions = numel(regions);
    regionNames = cell(1, nRegions);
    for i = 1 : nRegions
        if ~isnumeric(regions{i})
            regionNames{i} = regions{i};
            regions{i} = (spmReadVolMatched(char(regions{i}), VY(1)) > 0);
        else
            regionNames{i} = sprintf('region %d', i);
        end
    end
    try
        regions = (cat(4, regions{:}) > 0);
    catch
        error('region masks don''t match!')
    end
    assert(isequal(size(regions(:, :, :, 1)), size(mask)), ...
        'region masks don''t match!')
    % intersect brain mask with union of regions
    mask = mask & any(regions, 4);
    % determine mask voxel indices for each region
    regions = reshape(regions, [], nRegions);
    rmvi = cell(nRegions, 1);
    for i = 1 : nRegions
        rmvi{i} = find(regions(mask(:), i));
        fprintf('    %d voxels within brain mask in %s\n', numel(rmvi{i}), regionNames{i})
    end
    fprintf('  %d voxels within brain mask in all regions\n', sum(mask(:)));
end

% read and mask data
pattern = SPM.xY.P(1, :);
pattern(~all(diff(SPM.xY.P) == 0)) = '?';
if strfind(pattern, pwd) == 1
    pattern = pattern(numel(pwd) + 2: end);
end
fprintf('  reading images from %s\n', pattern)
[Y, mask] = spmReadVolsMasked(VY, mask);

% whitening and filtering
% get design matrix and nonsphericity
X = SPM.xX.X;
V = SPM.xVi.V;
if wf
    % whiten data matrix, design matrix, and nonsphericity
    if isfield(SPM.xX, 'W')
        fprintf('  whitening\n')
        W = SPM.xX.W;
        Y = W * Y;
        X = W * X;
        V = W * V * W';
    else
        fprintf('  * SPM.mat does not define whitening matrix!\n')
    end
    % high-pass filter data matrix, design matrix, and nonsphericity
    fprintf('  high-pass-filtering\n')
    Y = spm_filter(SPM.xX.K, Y);
    X = spm_filter(SPM.xX.K, X);
    V = spm_filter(SPM.xX.K, spm_filter(SPM.xX.K, V)')';
    % check consistency with SPM's results
    assert(norm(X - SPM.xX.xKXs.X) < SPM.xX.xKXs.tol)
    assert(sqrt(sum(sum((V - SPM.xX.V) .^ 2))) < SPM.xX.xKXs.tol)
end

% extract session-wise data matrix, design matrix, and regressor names
% and calculate degrees of freedom
m = numel(SPM.nscan);
Ys = cell(m, 1);
Xs = cell(m, 1);
fs = nan(m, 1);
names = cell(m, 1);
for si = 1 : m
    % identify rows belonging to session
    rows = SPM.Sess(si).row;
    % identify columns belonging to session from regressor names
    prefix = sprintf('Sn(%d) ', si);    % session prefix
    cols = find(cellfun(@(x) startsWith(x, prefix), SPM.xX.name));
    assert(isequal(cols, [SPM.Sess(si).col, SPM.xX.iB(si)]))
    % data matrix
    Ys{si} = Y(rows, :);
    % design matrix
    Xs{si} = X(rows, cols);
    % degrees of freedom
    [trRV, trRVRV] = spm_SpUtil('trRV', Xs{si}, V(rows, rows));
    fs(si) = trRV^2 / trRVRV;
    % regressor names, without prefix, as string array
    names{si} = string( ...
        cellfun(@(s) s(numel(prefix) + 1 : end), SPM.xX.name(cols), ...
        'UniformOutput', false));
end

% miscellaneous output
% analysis brain mask
misc.mask = mask;
% voxels to mm transformation
misc.mat = VY(1).mat;
% mask voxel indices for each region
misc.rmvi = rmvi;


% Copyright © 2013–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
