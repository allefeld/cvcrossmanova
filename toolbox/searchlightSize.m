function ret = searchlightSize(radius, nvargs)

% searchlight size as a function of searchlight radius
%
% p = searchlightSize(radius, mat = eye(3), file = "")
% tbl = searchlightSize(mat = eye(3), file = "")
% tbl = searchlightSize([radiusMin, radiusMax], mat = eye(3), file = "")
%
% With the first syntax, the number of voxels contained in a searchlight of
% radius `radius` is returned as `p`.
%
% With the other syntaxes, a table of corresponding values of `radius` and
% `p` is returned as `tbl`. By default, the table includes `radius` values
% from 0 to 5, but another range can be specified as `[radiusMin,
% radiusMax]`. The tabulated values of `radius` are chosen such that all
% possible searchlight sizes in that range are listed.
%
% Without further information, `radius` is in voxels. The optional `mat`
% can be used to specify a transformation matrix from voxel space to
% physical space, and then `radius` is in the same physical units as `mat`.
% Alternatively, the optional `file` can be used to specify the filename of
% an `SPM.mat` file or NIfTI file from which to read the transformation
% matrix, and then `radius` is in mm. If both are given, `file` overrides
% `mat`.

arguments
    radius      (1, 2) double  = [0 5]
    % no argument uses the default
    % 1-element argument is duplicated
    % 2-element argument is used as is
    % more than 2 elements is rejected
    nvargs.mat  (:, :) double  = eye(3)
    nvargs.file (1, :) char    = ''
end

% get range of radius
radiusMin = min(radius);
radiusMax = max(radius);

% get transformation from voxel space to physical space
mat = nvargs.mat;
if ~isempty(nvargs.file)
    % obtain volume from SPM.mat or NIfTI file
    if endsWith(nvargs.file, ".mat")
        load(nvargs.file, "SPM");
        V = SPM.xY.VY;
    else
        V = spm_vol(nvargs.file);
    end
    mat = V(1).mat;
end

% get distances from center voxel on grid
[~, dist] = searchlight(radiusMax, mat);

% first syntax
if radiusMin == radiusMax
    ret = nnz(dist <= radiusMax);
    return
end

% list radii and sizes
radius = unique(dist);
radius = radius(radius <= radiusMax);
p = nan(size(radius));
for i = 1 : numel(radius)
    % number of voxels within radius
    p(i) = nnz(dist <= radius(i));
end

% round up radii within necessary precision
for i = 1 : numel(radius) - 1
    % try increasing precisions
    for digits = 0 : 16
        prec = 10 ^ digits;
        % round up to precision
        rr = ceil(radius(i) * prec) / prec;
        % still smaller than next larger radius?
        sufficient = (rr < radius(i + 1));
        if sufficient
            break
        end
    end
    assert(sufficient)
    radius(i) = rr;
end

% limit to range
ind = ((radius >= radiusMin) & (radius <= radiusMax));
p = p(ind);
radius = radius(ind);

% create table
ret = table(radius, p);


% Copyright © 2016–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
