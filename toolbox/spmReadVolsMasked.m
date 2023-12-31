function [Y, mask] = spmReadVolsMasked(V, mask)

% read in-mask MR image data from an array of SPM volume structs
%
% Y = spmReadVolsMasked(V, mask)
% [Y, mask] = spmReadVolsMasked(V)
%
% V:    array of volume structs (or filenames)
% mask: logical array with dimensions matching data volumes
%       if not specified, read all voxels
% Y:    data (volumes × in-mask voxels)
%
% A replacement for spm_read_vols that reads only in-mask voxels to avoid
% memory spikes and to gain speed. This is achieved by directly accessing
% memory-mapped data via the volume structs' file_array object V.private.raw.
%
% See also spmUnMask, spmReadVolMatched, spmCoords

% This function is about 2–3 times faster than spm_read_vols on a local
% file system and about 20 times faster on SMB/CIFS.
% tested on both SPM8 and SPM12 with SPM.xY.VY generated by both SPM8 and SPM12


% support (cell array) of filename(s)
if iscell(V)
    V = char(V);
end
if ischar(V)
    V = spm_vol(V);
end

% check orientations
spm_check_orientations(V);

% check mask
if (nargin < 2) || isempty(mask)
    mask = true(V(1).dim);
end
if ~islogical(mask)
    if ~isequal(mask, logical(mask))
        fprintf(2, 'mask converted to logical!\n');
    end
    mask = logical(mask);
end
assert(numel(mask) == prod(V(1).dim), 'mask does not match data!')

% determine linear indices of in-mask voxels & corresponding image planes
maskInd = find(mask);
[~, ~, plane] = ind2sub(V(1).dim, maskInd);

% read and mask data
nVols = numel(V);
Y = nan(nVols, numel(maskInd));
reported = false;
line = '';
for vi = 1 : nVols
    try
        % shift mask indices according to sub-volume (n) of image file
        ind = maskInd + (V(vi).n(1) - 1) * prod(V(vi).private.raw.dim(1 : 3));
        if V(vi).n(2) > 1
            % including 5th dimension, if present
            ind = ind + (V(vi).n(2) - 1) * prod(V(vi).private.raw.dim(1 : 4));
        end
        % prepare linear transformation parameters
        if size(V(vi).pinfo, 2) == 1
            % uniform
            slope = V(vi).pinfo(1);
            inter = V(vi).pinfo(2);
        else
            % plane-specific
            slope = V(vi).pinfo(1, plane)';
            inter = V(vi).pinfo(2, plane)';
        end
        % access raw memory-mapped data and apply linear transform
        Y(vi, :) = double(V(vi).private.raw(ind)) .* slope + inter;
            % The last line could simply be
            %   Y(i, :) = V(i).private.dat(ind);
            % so we wouldn't need slope and inter. But unfortunately SPM8
            % doesn't apply global scaling to dat.scl_slope and
            % dat.scl_inter, so we have to use the scaled V(i).pinfo with
            % raw. Which is basically the same thing that spm_read_vols does.
    catch
        % if something goes wrong with accessing memory-mapped data,
        % fall back on spm_read_vols
        y = spm_read_vols(V(vi));
        Y(vi, :) = y(maskInd);
        if ~reported  % not yet reported fallback message?
            fprintf(2, ['error accessing memory-mapped data, ' ...
                'falling back on spm_read_vols!\n']);
            reported = true;
        end
    end
    % progress report
    if (mod(vi, 100) == 0) || (vi == nVols)
        fprintf(repmat('\b', 1, numel(line)))
        line = sprintf('  %d of %d volumes loaded', vi, nVols);
        fprintf(line)
    end
end
fprintf('\n')


% Copyright © 2016–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later

