function [estimable, error] = contrastEstimable(C, X)

% check whether a contrast is estimable w.r.t. a design matrix
%
% [estimable, error] = contrastEstimable(C, X)
% contrastEstimable(C, X)
%
% -----------  ------------------------------------------
% `C`          contrast matrix, regressors × subcontrasts
% `X`          design matrix, observations × regressors
% `estimable`  whether the contrast is estimable, logical
% `error`      degree of inestimability
% -----------  ------------------------------------------
%
% Estimability is determined by checking whether the contrast is contained
% in the row space of the design matrix (see Friston et al. 2007, p. 106).
% `error` is the underlying numerical measure of inestimability, which
% should be 0 for an estimable contrast, but may just be close to 0 due to
% numerical error. The criterion used is `estimable = (error < max(size(C))
% * eps)`.
%
% The second syntax reports the result instead of returning it.

if size(C, 1) < size(X, 2)
    if nargout == 0
        % interactive use
        fprintf("contrast does not extend across all regressors!\n")
    end
    C(size(X, 2), end) = 0;
end

% Friston et al. (2007), p. 106: "c is a contrast, if it is unchanged by
% post-multiplication with pinv(X' X) X' X. This test is used in SPM for
% user-specified contrasts." This function implements a slightly revised
% version of the procedure in SPM12 r7771, which is accessible via
% ```
% spm_SpUtil('allCon', X, C)
% ```
% On lines 279 and 286 it calls
% ```
% sX = spm_sp('Set', X)
% i = spm_sp('inspp', sX, C)
% ```
% which in turn is implemented in lines 1126–1135, 1218, and 947.

% calculate projection operator into the row space of X
[~, S, V] = svd(X, 0);
s = diag(S);
tol = max(size(X)) * eps(max(s));
ind = (s > tol);
opp = V(:, ind) * V(:, ind)';
% check whether contrast C is unchanged under this operator
error = norm(opp * C - C) / norm(C);
estimable = (error < max(size(C)) * 10 * eps);

% output for interactive use
if nargout == 0
    if estimable
        fprintf("contrast is estimable, ")
    else
        fprintf("contrast is not estimable, ")
    end
    fprintf("error: %g\n", error)
    clear estimable error
end

% Copyright © 2013–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
