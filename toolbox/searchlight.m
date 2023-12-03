function [PSL, dist] = searchlight(radius)

% create prototype searchlight
%
% [PSL, d] = searchlight(radius)
%
% PSL:   prototype searchlight matrix, 3-dimensional matrix
%        positive values indicate membership, value 2 indicates center
% dist:  distances from center voxel

[dx, dy, dz] = ndgrid(-ceil(radius) : ceil(radius));
dist = sqrt(dx .^ 2 + dy .^ 2 + dz .^ 2);
dist(dist > radius) = nan;

PSL = (dist <= radius) + (dist == 0);
