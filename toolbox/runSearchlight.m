function [res, ps] = runSearchlight(checkpoint, slRadius, mask, cm)

% searchlight: apply function to data contained in a sliding spherical window
%
% [res, p] = runSearchlight(checkpoint, slRadius, mask, cm)
%
% checkpoint:   name of checkpoint file ([] to disable checkpointing)
% slRadius:     radius of searchlight in voxels
% mask:         3-dimensional logical array indicating which voxels to use
% cm:           CrossManova object
% res:          cell array with results per analysis
% p:            number of voxels in each searchlight, column vector
%
% Run `slSize` for a table of meaningful values and resulting searchlight
% sizes.
%
% `res` is a cell array with one cell for each analysis included in
% `cm`, each containing an array of size voxels × results.
%
% Intermediate results are saved at regular intervals to the checkpoint
% file, if given. If the checkpoint file already exists, its contents are
% loaded and the computation continues from that point.
%
% See also slSize.


% normalize checkpoint file name, preserving emptiness
if ~isempty(checkpoint)
    [~, ~, ext] = fileparts(checkpoint);
    if ~strcmp(ext, '.mat')
        checkpoint = [checkpoint '.mat'];
    end
end

% volume dimensions
dim = size(mask);
nVolumeVoxels = prod(dim);
nMaskVoxels = sum(mask(:));

% determine searchlight voxel offsets relative to center voxel
% prototype searchlight
[dxi, dyi, dzi] = ndgrid(-ceil(slRadius) : ceil(slRadius));
PSL = (dxi .^ 2 + dyi .^ 2 + dzi .^ 2 <= slRadius .^ 2);
% spatial offsets
dxi = dxi(PSL);
dyi = dyi(PSL);
dzi = dzi(PSL);
% index offsets
PSL(dim(1), dim(2), dim(3)) = 0;        % zero-pad to full volume
di = find(PSL);
cInd = find((dxi == 0) & (dyi == 0) & (dzi == 0));
di = di - di(cInd);                                                         %#ok<FNDSB>
clear PSL cInd

% sort offsets by increasing distance from center
[~, ind] = sort(dxi .^ 2 + dyi .^ 2 + dzi .^ 2);
di = di(ind);
dxi = dxi(ind);
dyi = dyi(ind);
dzi = dzi(ind);

% mapping from volume to mask voxel indices
vvi2mvi = nan(nVolumeVoxels, 1);
vvi2mvi(mask) = 1 : nMaskVoxels;

% initialize result volume(s)
nResults = cm.nResults;
nAnalyses = numel(nResults);
res = cell(1, nAnalyses);
for i = 1 : nAnalyses
    res{i} = nan(nVolumeVoxels, nResults(i));
end
ps = nan(nVolumeVoxels, 1);

tic
t = 0;
line = '';
cvvi = 0;   % searchlight center *volume voxel index*
cmvi = 0;   % searchlight center *mask voxel index*
while cvvi < nVolumeVoxels
    cvvi = cvvi + 1;        % next volume voxel

    if (cvvi == 1) && ~isempty(checkpoint)
        % load checkpoint file, if existing
        if exist(checkpoint, 'file')
            load(checkpoint, 'res', 'ps', 'cvvi', 'cmvi')
            fprintf('    continuing from  %d voxels  %.1f %%\n', ...
                cmvi, cmvi / nMaskVoxels * 100)
        end
    end

    % is center within mask?
    if mask(cvvi)
        cmvi = cmvi + 1;    % next mask voxel

        % searchlight center coordinates
        [xi, yi, zi] = ind2sub(dim, cvvi);
        % searchlight voxel coordinates; limit to volume boundaries
        ind = (xi + dxi >= 1) & (xi + dxi <= dim(1)) & ...
            (yi + dyi >= 1) & (yi + dyi <= dim(2)) & ...
            (zi + dzi >= 1) & (zi + dzi <= dim(3));
        % searchlight voxel volume indices
        vvi = cvvi + di(ind);
        % discard out-of-mask voxels
        vvi = vvi(mask(vvi) == 1);
        % translate to mask voxel indices
        mvi = vvi2mvi(vvi);

        % call CrossManova object and store output
        Ds = cm.runAnalyses(mvi);
        for i = 1 : nAnalyses
            res{i}(cvvi, :) = Ds{i};
        end
        ps(cvvi) = numel(mvi);
    end

    nt = toc;
    if (nt - t > 30) || (cvvi == nVolumeVoxels)
        % progress report
        t = nt;
        fprintf(repmat('\b', 1, numel(line)))
        line = sprintf('    %.1f min  %d voxels  %.1f %%', ...
            t / 60, cmvi, cmvi / nMaskVoxels * 100);
        fprintf('%s', line)
        if (cvvi < nVolumeVoxels) && ~isempty(checkpoint)
            % save checkpoint file
            save(checkpoint, 'res', 'ps', 'cvvi', 'cmvi')
        end
    end
end
fprintf('\n')
% delete checkpoint file after completion
if ~isempty(checkpoint)
    spm_unlink(checkpoint)
end

% Copyright © 2013–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
