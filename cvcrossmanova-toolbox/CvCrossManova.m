classdef CvCrossManova < handle

    % data type representing data and several analyses

    properties
        Ys          % cell array of per-session design matrices
        Xs          % cell array of per-session data matrices
        analyses    % cell array of Analysis objects
        lambda      % strength of regularization (0–1) towards Euclidean metric
        m           % number of sessions
        ns          % array of per-session numbers of observations (rows)
        fs          % array of per-session residual degrees of freedom
        nVariables  % number of data variables (columns)
        hBetas      % cell array of per-session parameter estimates
        hXis        % cell array of per-session error estimates
        nAnalyses   % number of analyses
    end

    methods

        function self = CvCrossManova(Ys, Xs, analyses, kwargs)
            % create CvCrossManova object
            %
            % CvCrossManova(Ys, Xs, analyses, lambda=, fs=)
            %
            % ----------  ---------------------------------------------------------------------------------------------------------
            % `Xs`        cell array of per-session design matrices
            % `Ys`        cell array of per-session data matrices
            % `analyses`  cell array of Analysis objects
            % `lambda`    strength of regularization (0–1) towards Euclidean metric, default: 1e-8
            % `fs`        array of per-session residual degrees of freedom, default: computed from size and rank of design matrices
            % ----------  ---------------------------------------------------------------------------------------------------------
            %
            % Upon creation, a `CvCrossManova` object stores data matrices,
            % design matrices, analysis definitions, and further parameters, and
            % it estimates GLM parameters and errors. Actual analyses are then
            % performed on subsets of variables by calling the method
            % `runAnalyses`.

            arguments
                Ys             (:, :)  cell
                Xs             (:, :)  cell
                analyses       (:, :)  cell
                kwargs.lambda  (1, 1)  double  = 1e-8
                kwargs.fs      (:, :)  double  = []
            end

            % store arguments
            self.Ys = Ys(:)';
            self.Xs = Xs(:)';
            self.analyses = analyses(:)';
            self.lambda = kwargs.lambda;
            self.fs = kwargs.fs(:)';

            % determine number of sessions
            self.m = numel(self.Ys);
            assert(self.m == numel(self.Xs), ...
                'Number of sessions of Ys and Xs must match.');

            % determine number of observations for each session
            self.ns = cellfun(@(x) size(x, 1), self.Ys);
            assert(isequal(self.ns, cellfun(@(x) size(x, 1), self.Xs)), ...
                'Number of rows of Ys and Xs must match in each session.');

            % if not specified, calculate residual degrees of freedom for each session
            if isempty(self.fs)
                self.fs = nan(1, self.m);
                for k = 1 : self.m
                    self.fs(k) = self.ns(k) - rank(self.Xs{k});
                end
            else
                assert(self.m == numel(self.fs), ...
                    'Number of sessions of Ys and fs must match.');
            end

            % determine number of variables
            self.nVariables = size(self.Ys{1}, 2);
            assert(all(cellfun(@(x) size(x, 2), self.Ys) == self.nVariables), ...
                'Number of columns must match between sessions of Ys.');

            % estimate GLM parameters and errors for each session
            self.hBetas = cell(1, self.m);
            self.hXis = cell(1, self.m);
            for k = 1 : self.m
                self.hBetas{k} = pinv(self.Xs{k}) * self.Ys{k};
                self.hXis{k} = self.Ys{k} - self.Xs{k} * self.hBetas{k};
            end

            % process analysis definitions
            self.nAnalyses = numel(self.analyses);
            for i = 1 : self.nAnalyses
                if iscell(self.analyses{i})
                    % a cell array is interpreted as arguments for
                    % Analysis.leaveOneSessionOut
                    self.analyses{i} = Analysis.leaveOneSessionOut(self.m, ...
                        self.analyses{i}{:});
                end
            end
            assert(all(cellfun(@(x) x.m, self.analyses) == self.m), ...
                'Number of sessions of Ys and analyses must match.');

%             % check estimability of contrasts w.r.t. design matrices
%             for i = 1 : self.nAnalyses
%                 self.analyses{i}.checkEstimability(self.Xs);
%             end
        end

        function showAnalyses(self)
            % graphically displays all defined analyses
            for i = 1 : self.nAnalyses
                self.analyses{i}.show(sprintf('Analysis %d', i));
            end
        end

        function dispAnalyses(self)
            % summarizes all defined analyses
            for i = 1 : self.nAnalyses
                % sprintf('Analysis %d', i)
                self.analyses{i}.disp();
            end
        end

        function n = nResults(self)
            % returns the number of values returned by runAnalyses
            n = self.nAnalyses;           % * nPermutations
        end

        function Ds = runAnalyses(self, vi)
            % runs the defined analyses on a subset of variables
            %
            % D = runAnalyses(self, vi)
            % 
            % vi  variable indices, i.e. column indices into the data matrices
            %     default: all variables
            % Ds  cell array of analysis results,
            %     either pattern distinctness or pattern stability
            %     Each cell element is a scalar if no permutations have been
            %     applied, or a vector of permutation values where the first
            %     one corresponds to the neutral permutation.

            if nargin < 2
                vi = 1 : self.nVariables;
            else
                vi = vi(:)';
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
                    for k = 1 : nPerms
                        perm = double(permute(perms(k, :), [1, 3, 2]));

                        % this is "Strategy 1" in the paper
                        mhDBA = mean(cat(3, hDBAs{A}) .* perm(A), 3);
                        mXXn = mean(cat(3, XXns{B}), 3);
                        mhDB = mean(cat(3, hDBs{B}) .* perm(B), 3);
                        D(k) = D(k) + trace(mhDBA' * mXXn * mhDB);
                    end
                end
                Ds{i} = D / L;
            end

        end

    end

end
