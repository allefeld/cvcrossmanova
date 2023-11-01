function pMax = slSize(slRadius)

% searchlight size as a function of searchlight radius
%
% pMax = slSize(slRadius)
% slSizes(slRadius = 5)
%
% With the first syntax, the number of voxels corresponding to `slRadius`
% is returned as `pMax`. Note that this is the maximum number of voxels for
% a searchlight completely within the brain mask, at the borders the actual
% number will be smaller.
%
% With the second syntax, a table of corresponding values of `slRadius` and
% `pMax` is printed. By default, the table includes `slRadius` values from
% 0 to 5, but another upper limit can be specified. The values of
% `slRadius` are chosen such that all possible searchlight sizes in that
% range are listed.

if nargin == 0
    slRadius = 5;
end

% distances from center voxel on grid
[dxi, dyi, dzi] = ndgrid(-ceil(slRadius) : ceil(slRadius));
d = sqrt(dxi .^ 2 + dyi .^ 2 + dzi .^ 2);

% single query
if nargout > 0
    pMax = nnz(d <= slRadius);
    return
end

% tabulate radii and sizes
r = unique(d(d <= slRadius));
pMax = nan(size(r));
prec = 1 - floor(log10(min(diff(r))));      % necessary precision
fprintf('slRadius  pMax\n--------  ----\n')
for i = 1 : numel(r)
    % number of voxels within radius
    pMax(i) = nnz(d <= r(i));
    rd = ceil(r(i) * 10^prec) / 10^prec;    % round up
    fprintf('  %-5g   % 4d\n', rd, pMax(i))
    assert(numel(find(d <= rd)) == pMax(i)) % sufficient precision?
end
clear pMax

% Copyright © 2016–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
