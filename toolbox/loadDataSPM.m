function [Ys, Xs, mask, misc] = loadDataSPM(modelDir, regions, whitenfilter)

% load fMRI data via SPM.mat
%
% [Ys, Xs, mask, misc] = loadDataSPM(dirName, regions = {})
%
% modelDir: name of directory that contains SPM.mat
% regions:  optional additional region mask(s),
%           cell array of logical 3D volume(s) or filename(s)
% Ys:       MR data (within mask), cell array with one element for each
%           session, containing an array of size scans × voxels
% Xs:       design matrix for each session, cell array with one element
%           for each session, containing an array of size scans × regressors
% mask:     analysis brain mask, logical 3D volume;
%           possibly combined with union of region masks
% misc:     struct with additional data:
%     mat   voxels to mm transformation matrix
%     fs    residual degrees of freedom for each session
%     rmvi  cell array of mask voxel indices for each region
%
% Y & X and are high-pass filtered and whitened.
% Y includes only those voxels selected through mask.


% default argument values
if nargin < 2
    regions = {};
end
if nargin < 3
    whitenfilter = true;
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
fprintf('  %d voxels within brain mask\n', sum(mask(:)));

% possibly apply region mask(s)
if isempty(regions)
    fprintf('  no region mask\n')
    rmvi = {};
else
    nRegions = numel(regions);
    regionNames = cell(1, nRegions);
    for i = 1 : nRegions
        if ~isnumeric(regions{i})
            regionNames{i} = regions{i};
            regions{i} = (spmReadVolMatched(regions{i}, VY(1)) > 0);
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
    % restrict brain mask to conjunction of regions
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

% get design matrix
X = SPM.xX.X;

if whitenfilter
    % whiten data and design matrix
    if isfield(SPM.xX, 'W')
        fprintf('  whitening\n')
        W = SPM.xX.W;
        Y = W * Y;
        X = W * X;
    else
        fprintf('  * SPM.mat does not define whitening matrix!\n')
    end

    % high-pass filter data and design matrix
    fprintf('  high-pass-filtering\n')
    Y = spm_filter(SPM.xX.K, Y);
    X = spm_filter(SPM.xX.K, X);
end

% separate Y and X into session blocks; also for Bcov, W, XK = K.X0
m = numel(SPM.nscan);
Xs = cell(m, 1);
Ys = cell(m, 1);
Bcovs = cell(m, 1);
Ws = cell(m, 1);
XKs = cell(m, 1);
if isfield(SPM.xX, 'W')
    W = SPM.xX.W;
else
    W = speye(size(Y, 1), size(Y, 1));
end
for si = 1 : m
    Ys{si} = Y(SPM.Sess(si).row, :);
    % SPM.Sess(:).col does not include constant regressors,
    % get those from SPM.xX.iB
    col = [SPM.Sess(si).col, SPM.xX.iB(si)];
    Xs{si} = X(SPM.Sess(si).row, col);
    Bcovs{si} = SPM.xX.Bcov(col, col);
    Ws{si} = W(SPM.Sess(si).row, SPM.Sess(si).row);
    XKs{si} = SPM.xX.K(si).X0;
end
clear Y X

% degrees of freedom for each session
Tdf = nan(m, 1);
Kdf = nan(m, 1);
Xdf = nan(m, 1);
Rdf = nan(m, 1);
for si = 1 : m
    Tdf(si) = SPM.nscan(si);                    % total
    Kdf(si) = rank(SPM.xX.K(si).X0);            % loss from filter
    Xdf(si) = rank(Xs{si});                     % loss from regressors
    Rdf(si) = Tdf(si) - Kdf(si) - Xdf(si);      % residual
end
fprintf('  sum(fs) = %d - %d - %d = %d', sum(Tdf), sum(Kdf), sum(Xdf), sum(Rdf));
% other than SPM, we assume that whitening is perfect; for comparison
fprintf('   [SPM: trRV = %g  erdf = %g]\n', SPM.xX.trRV, SPM.xX.erdf)

% miscellaneous output
% voxels to mm transformation
misc.mat = VY(1).mat;
% residual degrees of freedom for each session
misc.fs = Rdf;
% mask voxel indices for each region
misc.rmvi = rmvi;
% parameter estimation covariance
misc.Bcovs = Bcovs;
if ~whitenfilter
    % whitening matrix
    misc.Ws = Ws;
    % high-pass filter regressors
    misc.XKs = XKs;
end

% Copyright © 2013–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
