function ret = searchlightSize(radius)

% searchlight size as a function of searchlight radius
%
% pMax = searchlightSize(radius)
% tbl = searchlightSize()
% tbl = searchlightSize([radiusMin, radiusMax])
%
% With the first syntax, the number of voxels corresponding to `radius` is
% returned as `pMax`.
%
% With the second and third syntax, a table of corresponding values of
% `radius` and `pMax` is returned. By default, the table includes `radius`
% values from 0 to 5, but other limits can be specified. The values of
% `radius` are chosen such that all possible searchlight sizes in that
% range are listed. If the output argument is omitted, the table is printed
% in Markdown-compatible form.

arguments
    radius  (1, 2) double  = [0 5]
    % no argument uses the default
    % 1-element argument is duplicated
    % 2-element argument is used as is
    % more than 2 elements is rejected
end

radiusMin = min(radius);
radiusMax = max(radius);

% compute distances from center voxel on grid
[dxi, dyi, dzi] = ndgrid(-ceil(radiusMax + 1) : ceil(radiusMax + 1));
d = sqrt(dxi .^ 2 + dyi .^ 2 + dzi .^ 2);

% first syntax
if radiusMin == radiusMax
    ret = nnz(d <= radiusMax);
    return
end

% list radii and sizes
radius = unique(d);
radius = radius(radius <= radiusMax + 1);
pMax = nan(size(radius));
for i = 1 : numel(radius)
    % number of voxels within radius
    pMax(i) = nnz(d <= radius(i));
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
pMax = pMax(ind);
radius = radius(ind);

% create table
ret = table(radius, pMax);

if nargout == 0
    % print table
    fprintf('`radius`   `pMax`\n--------  ------\n')
    fprintf('%-8g  %6d\n', [radius, pMax] .')
    % column widths are sufficient up to radius = 62.04, pMax = 999665
end


% Copyright © 2016–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
