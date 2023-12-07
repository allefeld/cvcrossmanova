function [PSL, dist] = searchlight(radius, mat)

% create prototype searchlight
%
% [PSL, d] = searchlight(radius)
% [PSL, d] = searchlight(radius, mat)
%
% radius:  radius of the searchlight sphere
% mat:     transformation from voxel space to physical space
%
% PSL:     prototype searchlight, 3-dimensional array
%          positive values indicate membership, value 2 indicates center
% dist:    distances from center voxel, 3-dimensional array
%
% The prototype searchlight approximates a sphere in physical space, which
% is derived from voxel space via the transformation `mat`. Due to the
% transformation, the prototype searchlight in general approximates an
% ellipsoid in voxel space, except when `mat` is a (scaled) identity
% matrix. If it is not specified, the default `eye(3)` is used.

if nargin < 2
    mat = eye(3);
end

% bounding box for searchlight in voxel space
%   The searchlight is a sphere in physical space,
% which in general is an ellipsoid in voxel space.
%   See https://math.stackexchange.com/questions/3926884
M = mat(1:3, 1:3) / radius;
B = ceil(sqrt(diag(inv(M' * M))));      % `floor` should be sufficient
% grid covering bounding box
% positions relative to center in voxel space
[di, dj, dk] = ndgrid(-B(1) : B(1), -B(2) : B(2), -B(3) : B(3));
dijk = [di(:), dj(:), dk(:)] .';
% positions relative to center in physical space
dxyz = mat(1:3, 1:3) * dijk;
% distances to center in physical space
dist = sqrt(sum(dxyz .^ 2));
dist = reshape(dist, size(di));

% prototype searchlight
PSL = (dist <= radius) + (dist == 0);


% [dx, dy, dz] = ndgrid(-ceil(radius) : ceil(radius));
% dist = sqrt(dx .^ 2 + dy .^ 2 + dz .^ 2);
% 
% PSL = (dist <= radius) + (dist == 0);
