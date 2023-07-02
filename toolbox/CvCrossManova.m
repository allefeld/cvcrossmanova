classdef CvCrossManova

    properties
        Ys
        Xs
        analyses
        % permute   % TODO
        lambda

        m
        ns
        fs
        nVariables
        hBetas
        hXis
        nAnalyses
    end

    methods

        function self = CrossManova(Ys, Xs, analyses, kwargs)
            arguments
                Ys cell
                Xs cell
                analyses cell
                % kwargs.permute logical = false
                kwargs.lambda double = 1e-8
                kwargs.fs = []
            end

            % store arguments
            self.Ys = Ys(:)';
            self.Xs = Xs(:)';
            self.analyses = analyses(:)';
            % self.permute = kwargs.permute;
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

            % check estimability of contrasts w.r.t. design matrices
            for i = 1 : self.nAnalyses
                self.analyses{i}.checkEstimability(self.Xs);
            end
        end

        function showAnalyses(self)
            for i = 1 : self.nAnalyses
                self.analyses{i}.show(sprintf('Analysis %d', i));
            end
        end

        function dispAnalyses(self)
            for i = 1 : self.nAnalyses
                % sprintf('Analysis %d', i)
                self.analyses{i}.disp();
            end
        end

        function n = nResults(self)
            n = self.nAnalyses;           % * nPermutations
        end

        function D = runAnalyses(self, vi)
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
            % A small default value of lambda = 1e-8 is used to avoid numerical
            % errors.

            % extract parameter estimates for specified variables,
            % and pre-whiten them by dividing by the Cholesky factor of Sigma
            hBs = cell(1, self.m);
            for k = 1 : self.m
                hBs{k} = self.hBetas{k}(:, vi) / chol(hSigma);
            end

            % run analyses
            D = nan(1, self.nAnalyses);
            for i = 1 : self.nAnalyses
                % extract analysis definition
                CA = self.analyses{i}.CA;
                CB = self.analyses{i}.CB;
                sessionsA = self.analyses{i}.sessionsA;
                sessionsB = self.analyses{i}.sessionsB;
                L = self.analyses{i}.L;

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
                D(i) = 0;
                for l = 1 : L
                    A = sessionsA(l, :);
                    B = sessionsB(l, :);
                    mhDBA = mean(cat(3, hDBAs{A}), 3);
                    mXXn = mean(cat(3, XXns{B}), 3);
                    mhDB = mean(cat(3, hDBs{B}), 3);
                    D(i) = D(i) + trace(mhDBA' * mXXn * mhDB);
                end
                D(i) = D(i) / L;
            end

        end

    end

end


% (//|#|%|<!--|;|/\\*|^|^[ \\t]*(-|\\d+.))\\s*($TAGS)
% (//|#|%|<!--|;|/\*|^|^\s*(-|\d+.))\s*($TAGS)