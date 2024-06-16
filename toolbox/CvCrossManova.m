classdef CvCrossManova < handle & matlab.mixin.Scalar

    % A `CvCrossManova` object implements the Cross-validated (Cross-)
    % MANOVA algorithm. It combines the data and design matrices
    % encapsulated in a `ModeledData` object with the analysis
    % specifications of one or more `Analysis` objects. If its method
    % `runAnalyses` is called for a set of dependent variables, all
    % specified analyses are performed on these variables concurrently.
    % This enables a more efficient implementation of the algorithm,
    % because partial results shared between analyses only have to be
    % computed once.

    properties
        modeledData     % data and design information, `ModeledData` object
        analyses        % analysis specifications, cell array of `Analysis` objects
        lambda          % amount of shrinkage regularization
        nAnalyses       % number of analyses
    end

    methods

        function obj = CvCrossManova(modeledData, analyses, nvargs)
            % create `CvCrossManova` object
            %
            % ccm = CvCrossManova(modeledData, analyses, lambda = 1e-8)
            %
            % `modeledData` is a `ModeledData` object encapsulating data
            % and design matrices for multiple sessions.
            %
            % `analyses` is a cell array of `Analysis` objects specifying
            % analyses.
            %
            % The optional `lambda` (from 0 to 1) controls the amount of
            % shrinkage regularization applied to the estimate of the error
            % covariance matrix. It should usually be kept at its default
            % value.


            arguments
                modeledData    (1, 1)  ModeledData
                analyses       (:, :)  cell
                nvargs.lambda  (1, 1)  double  = 1e-8
            end

            % store arguments
            obj.modeledData = modeledData;
            obj.analyses = analyses(:) .';
            obj.lambda = nvargs.lambda;

            % determine number of analyses
            obj.nAnalyses = numel(analyses);
            assert(all(cellfun(@(x) isa(x, 'Analysis'), analyses)), ...
                "analyses must contain Analysis objects.")
            assert(all(cellfun(@(x) x.m, analyses) == modeledData.m), ...
                "Number of sessions of modeledData and analyses must match.");

            % check estimability of contrasts w.r.t. design matrices
            [estimability, problems] = obj.checkEstimability();
            if problems
                warning("Some contrasts are not estimable.")
                fprintf("Estimability:\n")
                disp(estimability)
            end
        end

        function Ds = runAnalyses(obj, vi)
            % run analyses on (a subset of) the dependent variables
            %
            % Ds = ccm.runAnalyses(vi)
            % Ds = ccm.runAnalyses()
            %
            % `vi` specifies the variables to be included in the analysis,
            % as column indices into the data matrices `Ys` of
            % `modeledData`. If omitted, all variables are included.
            %
            % `Ds` is a cell array of analysis results where each cell
            % corresponds to a cell of `analyses`, either a scalar value or
            % an array of permutation values. Whether a result is an
            % estimate of pattern distinctness *D* or pattern stability
            % *D*^×^ depends on the contrasts of the corresponding analysis
            % and the regressors involved in them.
            %
            % To determine how many values will be included in each of the
            % cells of `Ds` before actually running the analyses, e.g. for
            % preallocation, use the method `nResults`.

            if nargin < 2
                vi = 1 : obj.modeledData.p;
            else
                vi = vi(:).';
            end

            % short names to simplify code
            md = obj.modeledData;
            m = md.m;

            % estimate error covariance matrix Sigma
            hXiXis = 0;
            for k = 1 : m
                x = md.hXis{k}(:, vi);      % faster if in variable
                hXiXis = hXiXis + x' * x;
            end
            % We estimate such that the inverse is an unbiased estimator
            % of the inverse of the true error covariance matrix.
            hSigma = hXiXis / (sum(md.fs) - numel(vi) - 1);

            % regularization target: scaled identity matrix
            hSigmaTarget = eye(numel(vi)) * mean(diag(hXiXis)) / (sum(md.fs) - 2);
            % Because in this case the inverse is the element-wise and
            % therefore 1-dimensional inverse, the correction to the
            % degrees of freedom, − p − 1, becomes − 2.
            % The scaled identity matrix can be seen as implementing a
            % Euclidean metric, because it preserves angles and relative
            % distances in original data space. Apart from regularization,
            % this is useful if the user intends to test orthogonality in
            % that space.
            clear hXiXis                % release memory

            % regularize towards target
            hSigma = (1 - obj.lambda) * hSigma + obj.lambda * hSigmaTarget;
            % The regularization parameter lambda interpolates between the
            % proper (lambda = 0) and the scaled-identity estimate (lambda
            % = 1).
            clear hSigmaTarget          % release memory

            % check condition number (invertibility)
            cond_hSigma = cond(hSigma);
            fprintf("Condition number of error covariance matrix: %g\n", ...
                cond_hSigma)
            if cond_hSigma > 1000
                warning("Spatial whitening may be unreliable; " ...
                    + "choose a larger regularization parameter lambda.")
            end

            % extract parameter estimates for specified variables and
            % spatially pre-whiten them by dividing by the Cholesky factor
            % of Sigma
            chol_hSigma = chol(hSigma);
            hBs = cell(1, m);
            for k = 1 : m
                hBs{k} = md.hBetas{k}(:, vi) / chol_hSigma;
            end
            clear hSigma chol_hSigma    % release memory

            % run analyses
            Ds = cell(obj.nAnalyses, 1);
            for i = 1 : obj.nAnalyses
                % short name to simplify code
                an = obj.analyses{i};

                % restrict to involved independent variables
                % dim 1
                riA = find(any(an.CA ~= 0, 2));
                riB = find(any(an.CB ~= 0, 2));

                % pre-compute parameter-effect matrices
                % and design matrix inner products for each session
                hDBAs = cell(1, m);
                hDBs = cell(1, m);
                XXns = cell(1, m);
                for k = 1 : m
                    if any(an.sessionsA(:, k))
                        hDBAs{k} = an.CB(riB, :) * pinv(an.CA(riA, :)) ...
                            * hBs{k}(riA, :);
                    end
                    if any(an.sessionsB(:, k))
                        hDBs{k} = an.CB(riB, :) * pinv(an.CB(riB, :)) ...
                            * hBs{k}(riB, :);
                        x = md.Xs{k}(:, riB);       % faster if in variable
                        XXns{k} = x' * x / md.ns(k);
                    end
                end

                % for each fold and permutation
                nPerms = size(an.perms, 1);
                D = zeros(1, nPerms);
                for l = 1 : an.L
                    A = an.sessionsA(l, :);
                    B = an.sessionsB(l, :);
                    for j = 1 : nPerms
                        perm = double(permute(an.perms(j, :), [1, 3, 2]));

                        % this is "Strategy 1" in the draft
                        mhDBA = mean(cat(3, hDBAs{A}) .* perm(A), 3);
                        mXXn = mean(cat(3, XXns{B}), 3);
                        mhDB = mean(cat(3, hDBs{B}) .* perm(B), 3);
                        D(j) = D(j) + trace(mhDBA' * mXXn * mhDB);
                    end
                end
                Ds{i} = D / an.L;
            end
        end

        function [ol, MSE] = optimizeLambda(obj, vi)
            % optimize regularization parameter for (a subset of) the dependent variables
            %
            % [ol, MSE] = ccm.optimizeLambda(vi)
            % [ol, MSE] = ccm.optimizeLambda()
            %
            % `vi` specifies the variables to be included in the analysis,
            % as column indices into the data matrices `Ys` of
            % `modeledData`. If omitted, all variables are included.
            %
            % `ol` is the optimal shrinkage regularization parameter
            % `lambda`, a number between 0 and 1.
            %
            % `MSE` is the mean squared error of whitening at `ol`.
            %
            % The optimal value is determined via leave-one-session-out
            % cross-validation. The error covariance matrix is estimated
            % from 'training' and 'validation' sessions, and the latter is
            % whitened using a regularized version of the former. The
            % deviation from the optimal result, the identity matrix, is
            % quantified by the squared error, averaged across matrix
            % elements and cross-validation folds. The regularization
            % parameter is chosen such that the mean squared error is
            % minimal. Because the 'training' error covariance matrix used
            % for whitening is calculated from less than the full number of
            % samples, `ol` will likely overestimate the optimum for the
            % actual analysis.

            if nargin < 2
                vi = 1 : obj.modeledData.p;
            else
                vi = vi(:).';
            end

            % short names to simplify code
            md = obj.modeledData;
            m = md.m;

            % compute inner products
            hXiXis = cell(1, m);
            for k = 1 : m
                x = md.hXis{k}(:, vi);      % faster if in variable
                hXiXis{k} = x' * x;
            end

            % pre-compute 'training' error variance, regularization target,
            % and 'validation' error variance for each fold
            hSigmaA = cell(1, m);
            hSigmaATarget = cell(1, m);
            hSigmaB = cell(1, m);
            I = eye(numel(vi));
            for ll = 1 : m
                A = ((1 : m) ~= ll);    % 'training' sessions
                B = ((1 : m) == ll);    % 'validation' session
                % 'training' error variance
                % degrees of freedom for unbiased inverse
                hXiXisA = sum(cat(3, hXiXis{A}), 3);
                fA = sum(md.fs(A));
                hSigmaA{ll} = hXiXisA / (fA - numel(vi) - 1);
                % regularization target, scaled identity matrix
                hSigmaATarget{ll} = I * mean(diag(hXiXisA)) / (fA - 2);
                % 'validation' error variance
                % degrees of freedom for unbiased estimate
                hXiXisB = sum(cat(3, hXiXis{B}), 3);
                fB = sum(md.fs(B));
                hSigmaB{ll} = hXiXisB / fB;
            end

            % objective function
            function MSE = MSEfun(lambda)
                MSE = 0;
                % leave-one-session-out cross-validation
                for l = 1 : m
                    % regularize 'training' error variance by lambda
                    hSigmaAReg = (1 - lambda) * hSigmaA{l} + lambda * hSigmaATarget{l};
                    % use it to whiten 'validation' error variance,
                    cholhSigmaAReg = chol(hSigmaAReg);
                    % and quantify error
                    error = cholhSigmaAReg' \ hSigmaB{l} / cholhSigmaAReg - I;
                    MSE = MSE + mean(error(:) .^ 2);  % equivalent to Frobenius norm
                end
                MSE = MSE / m;
            end

            % determine lambda optimal for whitening
            % by minimizing the objective function
            [ol, MSE] = fminbnd(@MSEfun, 0, 1);
        end

        function nr = nResults(obj)
            % return the number of values returned by `runAnalyses`
            %
            % nr = ccm.nResults()
            %
            % `nr` is an array where each element corresponds to an
            % analysis, and the value specifies the number of results
            % returned for that analysis. For an analysis without
            % permutations, that value is 1; otherwise it is the number of
            % permutations.

            nr = cellfun(@(x) size(x.perms, 1), obj.analyses);
        end

        function [estimability, problems] = checkEstimability(obj)
            % check estimability of all contrasts  of all analyses in all sessions
            %
            % [estimability, problems] = ccm.checkEstimability()
            %
            % `estimability` is a table with one row per session and one
            % column per analysis & contrast. 'true' means that the
            % contrast is estimable, 'false' that it is not estimable, '–'
            % that the contrast does not apply to the session.
            %
            % `problems` indicates whether there are any inestimable
            % contrasts (logical).
            %
            % It is not usually necessary to use this method explicitly,
            % because estimability is checked upon creation of a
            % `CvCrossManova` object.

            m = obj.modeledData.m;
            estimability = table('RowNames', ...
                cellfun(@(k) sprintf("session %d", k), num2cell(1 : m)));
            problems = false;   % shared with nested function `est`
            for i = 1 : obj.nAnalyses
                analysis = obj.analyses{i};
                if ~isequal(analysis.CA, analysis.CB)
                    % different contrasts checked in their respective sessions
                    estimability.(sprintf("%d A", i)) = est(...
                        any(analysis.sessionsA, 1), ...
                        analysis.CA);
                    estimability.(sprintf("%d B", i)) = est(...
                        any(analysis.sessionsB, 1), ...
                        analysis.CB);
                else
                    % same contrast checked in all sessions
                    estimability.(sprintf("%d", i)) = est(...
                        any([analysis.sessionsA ; analysis.sessionsB], 1), ...
                        analysis.CA);
                end
            end
            % nested function to check estimability of one contrast over
            % multiple sessions
            function e = est(sessions, C)
                e = repmat("–", m, 1);
                for kk = find(sessions)
                    e(kk) = contrastEstimable(C, obj.modeledData.Xs{kk});
                end
                e = categorical(e);
                problems = problems || ismember("false", e);
            end
        end

        function disp(obj)
            % textually display information about the object
            %
            % ccm.disp()
            %
            % This method overrides Matlab's `disp`, so you can also use
            % `disp(ccm)` or simply `ccm` without semicolon to
            % get the same output.

            % TODO make more informative using regressor labels if available

            % prepare information string
            str = "  CvCrossManova:";
            str = str + sprintf("\n    ModeledData:");
            lines = splitlines(obj.modeledData.disp());
            for j = 2 : numel(lines)
                str = str + sprintf("\n  %s", lines(j));
            end
            for i = 1 : obj.nAnalyses
                str = str + sprintf("\n    Analysis %d:", i);
                lines = splitlines(obj.analyses{i}.disp());
                for j = 2 : numel(lines)
                    str = str + sprintf("\n  %s", lines(j));
                end
            end
            str = str + sprintf("\n    lambda: %g", obj.lambda);
            % display information string
            disp(str)
        end

    end

end

% Copyright © 2016–2024 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
