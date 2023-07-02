function [res, p] = runSearchlight(checkpoint, slRadius, mask, cm)

% general purpose searchlight:
% apply function to data contained in a sliding spherical window
%
% [res, p] = runSearchlight(checkpoint, slRadius, mask, cm)
%
% checkpoint:   name of checkpoint file ([] to disable checkpointing)
% slRadius:     radius of searchlight in voxels
% mask:         3-dimensional logical array indicating which voxels to use
% cm:           CrossManova object
% res:          results, array of size voxels × output values
% p:            number of voxels in each searchlight, column vector
%
% A voxel is included in the searchlight if its distance from the center is
% *smaller than or equal to* the radius. Note that fractional values are
% possible. Run slSize for a table of meaningful values and resulting
% searchlight sizes
%
% Intermediate results are saved at regular intervals to the checkpoint
% file, if given. On a subsequent run, if the checkpoint file exists, its
% contents are loaded and the computation continues from that point. 
%
% See also slSize
%
%
% This file is part of v3 of cvmanova, see
% https://github.com/allefeld/cvmanova/releases
%
% Copyright (C) 2013–2019 Carsten Allefeld


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
res = nan(nVolumeVoxels, cm.nResults);
p = nan(nVolumeVoxels, 1);

tic
t = 0;
cvvi = 0;   % searchlight center *volume voxel index*
cmvi = 0;   % searchlight center *mask voxel index*
while cvvi < nVolumeVoxels
    cvvi = cvvi + 1;        % next volume voxel

    if (cvvi == 1) && ~isempty(checkpoint)
        % load checkpoint file, if existing
        if exist(checkpoint, 'file')
            load(checkpoint, 'res', 'p', 'cvvi', 'cmvi')
            fprintf('  *restart*  %6d voxels  %5.1f %%\n', ...
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
        res(cvvi, :) = cm.runAnalyses(mvi);
        p(cvvi) = numel(mvi);
    end
    
    nt = toc;
    if (nt - t > 30) || (cvvi == nVolumeVoxels)
        % progress report
        t = nt;
        fprintf(' %6.1f min  %6d voxels  %5.1f %%\n', ...
            t / 60, cmvi, cmvi / nMaskVoxels * 100)
        
        if (cvvi < nVolumeVoxels) && ~isempty(checkpoint)
            % save checkpoint file
            save(checkpoint, 'res', 'p', 'cvvi', 'cmvi')
        end
    end
end
% delete checkpoint file after completion
if ~isempty(checkpoint)
    spm_unlink(checkpoint)
end
