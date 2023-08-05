classdef CvCrossManova < handle

    % data type representing data and analyses

    properties
        Ys          % cell array of per-session design matrices, observations × variables
        Xs          % cell array of per-session data matrices, observations × regressors
        analyses    % cell array of Analysis objects
        lambda      % strength of regularization (0–1) towards Euclidean metric
        fs          % array of per-session residual degrees of freedom
        m           % number of sessions
        ns          % array of per-session numbers of observations (rows)
    end

    properties (Hidden = true)
        nVariables  % number of data variables (columns)
        hBetas      % cell array of per-session parameter estimates
        hXis        % cell array of per-session error estimates
        nAnalyses   % number of analyses
    end

    properties (Dependent)
        nResults    % array of per-analysis returned numbers of results
    end

    methods

        function self = CvCrossManova(Ys, Xs, analyses, kwargs)
            % create `CvCrossManova` object
            %
            % ccm = CvCrossManova(Ys, Xs, analyses, lambda = 1e-8, fs = ...)
            %
            % A `CvCrossManova` object stores data matrices, design
            % matrices, analysis definitions, and further parameters. Upon
            % creation, it also estimates GLM parameters and errors. Actual
            % analyses are then performed on subsets of variables by
            % calling the method `runAnalyses`.
            %
            % The parameter λ (from 0 to 1) controls the amount of
            % shrinkage regularization applied to the estimate of the error
            % covariance matrix. A small value can be used to improve the
            % numerical and statistical stability of the estimate; the
            % default is 10^−8^. A value of 1 can be used to disregard the
            % covariance structure, because it replaces the estimated error
            % covariance matrix by a diagonal matrix where every diagonal
            % element is the average variance across variables. This can be
            % useful for Cross-MANOVA if it is intended to quantify
            % orthogonality w.r.t the original data space (Euclidean
            % metric) instead of the whitened space.
            %
            % If the per-session residual degrees of freedom `fs` are not
            % specified, they are calculated under the assumption that the
            % data observations are uncorrelated, as `ns(k) - rank(Xs{k})`.
            % If the observations have been only approximately decorrelated
            % or have been filtered, correct values should be explicitly
            % specified.

            arguments
                Ys             (:, :)  cell
                Xs             (:, :)  cell
                analyses       (:, :)  cell
                kwargs.lambda  (1, 1)  double  = 1e-8
                kwargs.fs      (:, :)  double  = []
            end

            % store arguments
            self.Ys = Ys(:).';
            self.Xs = Xs(:).';
            self.analyses = analyses(:).';
            self.lambda = kwargs.lambda;
            self.fs = kwargs.fs(:).';

            % determine number of sessions
            self.m = numel(self.Ys);
            assert(self.m == numel(self.Xs), ...
                "Number of sessions of Ys and Xs must match.");

            % determine number of observations for each session
            self.ns = cellfun(@(x) size(x, 1), self.Ys);
            assert(isequal(self.ns, cellfun(@(x) size(x, 1), self.Xs)), ...
                "Number of rows of Ys and Xs must match in each session.");

            % if not specified, calculate residual degrees of freedom for each session
            if isempty(self.fs)
                self.fs = self.ns - cellfun(@rank, self.Xs);
            else
                assert(self.m == numel(self.fs), ...
                    "Number of sessions of Ys and fs must match.");
            end

            % determine number of variables
            self.nVariables = size(self.Ys{1}, 2);
            assert(all(cellfun(@(x) size(x, 2), self.Ys) == self.nVariables), ...
                "Number of columns must match between sessions of Ys.");

            % estimate GLM parameters and errors for each session
            self.hBetas = cell(1, self.m);
            self.hXis = cell(1, self.m);
            for k = 1 : self.m
                self.hBetas{k} = pinv(self.Xs{k}) * self.Ys{k};
                self.hXis{k} = self.Ys{k} - self.Xs{k} * self.hBetas{k};
            end

            % determine number of analyses
            self.nAnalyses = numel(self.analyses);
            assert(all(cellfun(@(x) isequal(class(x), 'Analysis'), analyses)), ...
                "analyses must contain Analysis objects.")
            assert(all(cellfun(@(x) x.m, self.analyses) == self.m), ...
                "Number of sessions of Ys and analyses must match.");

            % check estimability of contrasts w.r.t. design matrices
            [estimability, problems] = self.checkEstimability();
            if problems
                warning("Some contrasts are not estimable.")
                fprintf("Estimability:\n")
                disp(estimability)
            end
        end

        function Ds = runAnalyses(self, vi)
            % run the defined analyses on (a subset of) the variables
            %
            % Ds = ccm.runAnalyses(vi)
            % Ds = ccm.runAnalyses()
            %
            % `vi` specifies the variables to be included in the analysis,
            % as column indices into the data matrices `Ys`. If omitted,
            % all variables are included.
            %
            % `Ds` is a cell array of analysis results where each element
            % results from the corresponding element in `analyses`. Whether
            % this result is an estimate of pattern distinctness *D* or
            % pattern stability *D*^×^ depends on the contrasts and
            % regressors involved in them.
            % 
            % For an analysis which does not include permutations, the cell
            % array element is a scalar; otherwise it is an array of
            % permutation values, where the first element is the actual
            % estimate (corresponding to the neutral permutation).
            %
            % To determine how many values will be included in each of the
            % elements of `Ds` before actually running an analysis, e.g.
            % for preallocation, use the property `nResults`. It is a
            % vector where each element gives the number of values returned
            % for the corresponding analysis (i.e. the number of
            % permutations).

            if nargin < 2
                vi = 1 : self.nVariables;
            else
                vi = vi(:).';
            end

            % determine number of variables
            p = numel(vi);

            % estimate error covariance matrix Sigma for specified variables
            hXiXis = 0;
            for k = 1 : self.m
                x = self.hXis{k}(:, vi);
                hXiXis = hXiXis + x' * x;
            end
            % we estimate such that the inverse is an unbiased estimator
            hSigma = hXiXis / (sum(self.fs) - p - 1);

            % regularize Sigma towards Euclidean metric
            hSigmaEuc = eye(p) * mean(diag(hXiXis)) / (sum(self.fs) - 2);
            hSigma = (1 - self.lambda) * hSigma + self.lambda * hSigmaEuc;
            % Euclidean metric means that an estimate of the error covariance
            % matrix is used which is a scaled identity matrix. Because in this
            % case the inverse is the element-wise and therefore 1-dimensional
            % inverse, the correction -p-1 becomes -2.
            % The regularization parameter lambda interpolates between the
            % proper (lambda = 0) and the scaled-identity estimate (lambda = 1).
            % A small default value of is used to avoid numerical instability.

            % extract parameter estimates for specified variables,
            % and pre-whiten them by dividing by the Cholesky factor of Sigma
            hBs = cell(1, self.m);
            for k = 1 : self.m
                hBs{k} = self.hBetas{k}(:, vi) / chol(hSigma);
            end

            % run analyses
            Ds = cell(1, self.nAnalyses);
            for i = 1 : self.nAnalyses
                % extract analysis definition
                CA = self.analyses{i}.CA;
                CB = self.analyses{i}.CB;
                sessionsA = self.analyses{i}.sessionsA;
                sessionsB = self.analyses{i}.sessionsB;
                L = self.analyses{i}.L;
                perms = self.analyses{i}.perms;

                % restrict to involved regressors
                % dim 1
                riA = find(any(CA ~= 0, 2));
                riB = find(any(CB ~= 0, 2));

                % pre-compute parameter-effect matrices
                % and design matrix inner products for each session
                hDBAs = cell(1, self.m);
                hDBs = cell(1, self.m);
                XXns = cell(1, self.m);
                for k = 1 : self.m
                    if any(sessionsA(:, k))
                        hDBAs{k} = CB(riB, :) * pinv(CA(riA, :)) * hBs{k}(riA, :);
                    end
                    if any(sessionsB(:, k))
                        hDBs{k} = CB(riB, :) * pinv(CB(riB, :)) * hBs{k}(riB, :);
                        XXns{k} = self.Xs{k}(:, riB)' * self.Xs{k}(:, riB) / self.ns(k);
                    end
                end

                % for each fold
                nPerms = size(perms, 1);
                D = zeros(1, nPerms);
                for l = 1 : L
                    A = sessionsA(l, :);
                    B = sessionsB(l, :);
                    for j = 1 : nPerms
                        perm = double(permute(perms(j, :), [1, 3, 2]));

                        % this is "Strategy 1" in the paper
                        mhDBA = mean(cat(3, hDBAs{A}) .* perm(A), 3);
                        mXXn = mean(cat(3, XXns{B}), 3);
                        mhDB = mean(cat(3, hDBs{B}) .* perm(B), 3);
                        D(j) = D(j) + trace(mhDBA' * mXXn * mhDB);
                    end
                end
                Ds{i} = D / L;
            end
        end

        function nr = get.nResults(self)
            % return the number of values returned by `runAnalyses`

            nr = cellfun(@(x) size(x.perms, 1), self.analyses);
        end

        function dispAnalyses(self)
            % textually display information about all analyses
            %
            % ccm.dispAnalyses()

            for i = 1 : self.nAnalyses
                fprintf("Analysis %d\n", i)
                self.analyses{i}.disp();
            end
        end

        function showAnalyses(self)
            % graphically display information about all analyses
            %
            % ccm.showAnalyses()

            for i = 1 : self.nAnalyses
                fig = self.analyses{i}.show();
                fig.Name = sprintf("Analysis %d", i);
            end
        end

        function [estimability, problems] = checkEstimability(self)
            % check estimability of all analyses' contrasts in all sessions
            %
            % [estimability, problems] = ccm.checkEstimability()
            %
            % `estimability` is a table with one row per session and one
            % column per analysis & contrast. 'true' means that the
            % contrast is estimable, 'false' that it is not estimable, '–'
            % that the contrast does not apply to the session. The check is
            % performed via the function `contrastEstimable`.
            %
            % `problems` indicates whether there are any inestimable
            % contrasts (logical).
            %
            % It is not usually necessary to use this method explicitly,
            % because estimability is checked upon creation of a
            % `CvCrossManova` object.

            estimability = table('RowNames', ...
                cellfun(@(k) sprintf("session %d", k), num2cell(1 : self.m)));
            problems = false;   % shared with nested function `est`
            for i = 1 : self.nAnalyses
                analysis = self.analyses{i};
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
                e = repmat("–", self.m, 1);
                for kk = find(sessions)
                    e(kk) = contrastEstimable(C, self.Xs{kk});
                end
                e = categorical(e);
                problems = problems || ismember("false", e);
            end

        end

    end

end

% TODO give asserts and warnings ids