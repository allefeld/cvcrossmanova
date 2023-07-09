classdef Analysis < handle

    % data type representing an analysis

    properties
        CA         % 'training' contrast
        CB         % 'validation' contrast
        sessionsA  % 'training' sessions for each fold, L × m logical array
        sessionsB  % 'validation' sessions for each fold, L × m logical array
        same       % whether CA and CB are the same
        L          % number of folds
        m          % number of sessions
        perms      % sign permutations
    end

    methods

        function self = Analysis(CA, CB, sessionsA, sessionsB)
            % create Analysis object
            %
            % analysis = Analysis(CA, CB, sessionsA, sessionsB)

            arguments
                CA         (:, :)  double
                CB         (:, :)  double
                sessionsA  (:, :)  logical
                sessionsB  (:, :)  logical
            end

            % store arguments
            self.CA = CA;
            self.CB = CB;
            self.sessionsA = sessionsA;
            self.sessionsB = sessionsB;

            % TODO check whether sessionsA and sessionsB overlap in any fold

            % determine whether contrasts are the same
            self.same = isequaln(self.CA, self.CB);

            % determine number of folds and sessions
            [self.L, self.m] = size(self.sessionsA);
            assert(isequaln([self.L, self.m], size(self.sessionsB)), ...
                "sessionsA and sessionsB must have the same size")

            % initialize permutations to neutral only
            self.perms = ones(1, self.m);
        end

        function addPermutations(self, maxPerms)
            % add unique sign permutations of per-session parameter estimates
            %
            % analysis.addPermutations(maxPerm = 1000)
            %
            % This method adds information to `analysis` that
            % sign-permutations should be applied, so that different values
            % of *D* for the different permutations are computed.
            %
            % ::: {.callout-warning}
            % Sign permutations are used to test a per-subject null
            % hypothesis of no effect in standard MANOVA analyses. It does
            % not make sense to apply them to Cross-MANOVAs.
            % :::
            % 
            % For *m* sessions there are formally 2^*m*^ per-session sign
            % permutations, but not all of them lead to different
            % permutation values. The number of unique permutations depends
            % on the 'training' and 'validation' sessions used in different
            % folds (`sessionsA`, `sessionsB`). This method determines
            % which permutations lead to different outcomes and makes sure
            % only those unique permutations are performed.

            % Two sign permutations are equivalent in a fold if the relative signs
            % applied to all involved sessions are the same. Two sign permutations are
            % equivalent if they are equivalent in every fold.

            arguments
                self
                maxPerms  (1, 1)  double  = 1000
            end

            assert(self.m <= 21, "Too many possible sign permutations to process.")

            % number of sign permutations
            nPerms = 2 ^ self.m;
            fprintf("%d possible sign permutations, ", nPerms)
            
            % Add warning that this only works for non-cross. Based on .same?
            
            % generate all sign permutations of sessions (permutations × sessions)
            % with the neutral permutation first
            self.perms = 1 - 2 * (int8(dec2bin(0 : nPerms - 1, self.m)) - '0');
            
            pINA = [];
            % for each fold
            for l = 1 : self.L
                % reduce permutations to Involved sessions
                involved = self.sessionsA(l, :) + self.sessionsB(l, :) > 0;
                pI = self.perms(:, involved);
                % Normalize sign permutations by multiplying with the first value,
                % because only relative signs are important
                pIN = pI(:, 1) .* pI;
                % collect involved, normalized permutations across All folds
                pINA = [pINA, pIN];                                                     %#ok<AGROW> 
            end
            % determine unique (non-equivalent) permutations
            [~, ind, ~] = unique(pINA, "rows");
            % select these permutations in order
            self.perms = self.perms(sort(ind), :);
            nPerms = size(self.perms, 1);
            fprintf("%d of them unique\n", nPerms)
            
            % Monte Carlo test, following Dwass (1957) and Ernst (2004)
            if nPerms > maxPerms
                fprintf("randomly selecting a subset of %d permutations\n", maxPerms)
                % neutral permutation (1) + a random sample with replacement
                % from permutations 2 : nPerms, for a total of maxPerms
                ind = [1, randperm(nPerms - 1, maxPerms - 1) + 1];
                self.perms = self.perms(sort(ind), :);
            end
        end

        function checkEstimability(self, Xs)
            % TODO
        end

        % TODO add function to check whether cross or not based on
        % regressor ids?

        function disp(self)
            % TODO add permutations information
            % see https://uk.mathworks.com/help/matlab/ref/matlab.mixin.customdisplay-class.html
            % for adding to instead of replacing original `disp`
            str = sprintf("Analysis:\n");
            str = str + sprintf("  %d folds, %d sessions\n", self.L, self.m);
            str = str + sprintf("  CA:     %d x %d, %d-dimensional\n", ...
                size(self.CA), rank(self.CA));
            if ~self.same
                str = str + sprintf("  CB:     %d x %d, %d-dimensional\n", ...
                    size(self.CB), rank(self.CB));
                % check cross parameter-effect extracting matrix
                CBpCA = self.CB * pinv(self.CA);
                sf = sum(CBpCA(:) .^ 2);
                rk = rank(CBpCA);
                str = str + sprintf("  CA->CB: %d-dimensional, %g%% variance\n", ...
                    rk, sf / rk * 100);
                if sf / rk - 1 > 1e-14
                    str = str + "  Warning: Variance is not preserved!\n";
                end
            end
            disp(str)
        end

        function fig = show(self, name)
            % TODO add permutations information
            fig = figure;
            if nargin > 1
                fig.Name = name;
            end

            % colormap red → white → green
            cmap = interp1(0 : 2, [1 0 0 ; 1 1 1 ; 0 1 0], linspace(0, 2, 511));
            colormap(cmap)

            subplot(2, 2, 1)
            heatmap(self.CA, ColorMap=cmap, ColorbarVisible=false)
            clim([-1 1])
            title('Contrast A')
            xlabel('subcontrasts')
            ylabel('regressors')

            subplot(2, 2, 2)
            heatmap(self.CB, ColorMap=cmap, ColorbarVisible=false)
            clim([-1 1])
            title('Contrast B')
            xlabel('subcontrasts')
            ylabel('regressors')

            subplot(2, 2, 3)
            heatmap(double(self.sessionsA), ColorMap=cmap, ColorbarVisible=false)
            clim([-1 1])
            title('Sessions A')
            xlabel('sessions')
            ylabel('folds')

            subplot(2, 2, 4)
            heatmap(double(self.sessionsB), ColorMap=cmap, ColorbarVisible=false)
            clim([-1 1])
            title('Sessions B')
            xlabel('sessions')
            ylabel('folds')
        end

    end

    methods (Static)

        function analysis = leaveOneSessionOut(m, CA, CB)
            % create Analysis object for leave-one-session-out cross-validation
            %
            % Analysis.leaveOneSessionOut(m, CA, CB)
            %
            % This is a convenience method as an alternative to calling the
            % constructor with manually specified `sessionsA` and
            % `sessionsB`.
            arguments
                m   (1, 1)  double
                CA  (:, :)  double
                CB  (:, :)  double  = CA
            end

            sessionsB = logical(eye(m));
            sessionsA = not(sessionsB);
            analysis = Analysis(CA, CB, sessionsA, sessionsB);
        end

    end

end
