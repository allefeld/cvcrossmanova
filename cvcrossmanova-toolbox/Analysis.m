classdef Analysis < handle

    % data type representing an analysis

    properties
        CA           % 'training' contrast, regressors × subcontrasts
        CB           % 'validation' contrast, regressors × subcontrasts
        sessionsA    % 'training' sessions, folds × sessions logical
        sessionsB    % 'validation' sessions, folds × sessions logical
        L            % number of folds
        m            % number of sessions
        perms        % sign permutations, permutations × sessions
    end

    properties (Hidden = true)
        dimensionA   % rank of `CA`
        dimensionB   % rank of `CB`
        dimensionAB  % dimension (rank) of cross-contrast
        normAB       % norm of cross-contrast
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
            % TODO check whether CB = P CA; isequal(sortrows(CA), sortrows(CB))

            % determine number of folds and sessions
            [self.L, self.m] = size(self.sessionsA);
            % check
            assert(isequaln([self.L, self.m], size(self.sessionsB)), ...
                "cvcrossmanova:sessionssize", ...
                "sessionsA and sessionsB must have the same size.")
            
            % characterize contrasts & check
            self.dimensionA = rank(self.CA);
            self.dimensionB = rank(self.CB);
            % check
            assert(size(CA, 2) == size(CB, 2), ...
                "cvcrossmanova:subcontrasts", ...
                "CA and CB must have the same number of subcontrasts (columns).")
            if self.dimensionA ~= self.dimensionB
                warning("cvcrossmanova:dimensions", ...
                    "CA and CB should have the same dimension.")
            end

            % characterize cross-contrast
            CBpCA = self.CB * pinv(self.CA);    % (cross-) parameter-effect extracting matrix
            self.dimensionAB = rank(CBpCA);
            self.normAB = sum(CBpCA(:) .^ 2);
            % check
            if self.dimensionAB < min(self.dimensionA, self.dimensionB)
                warning("cvcrossmanova:preservedim", ...
                    "CA → CB should preserve dimensionality.")
            end
            if abs(self.normAB / self.dimensionAB - 1) > 1e-14
                warning("cvcrossmanova:preservevar", ...
                    "CA → CB should preserve variance.")
            end

            % initialize permutations to neutral only
            self.perms = ones(1, self.m);
        end

        function addPermutations(self, maxPerms)
            % add sign permutations of per-session parameter estimates
            %
            % analysis.addPermutations(maxPerms = 1000)
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
            %
            % If the number of unique permutations exceeds `maxPerms`, a
            % randomly chosen subset of `maxPerms` is chosen (including the
            % neutral permutation, which is always permutation 1).

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
%             fprintf("%d possible sign permutations, ", nPerms)
            
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
%             fprintf("%d of them unique\n", nPerms)
            
            % Monte Carlo test, following Dwass (1957) and Ernst (2004)
            if nPerms > maxPerms
%                 fprintf("randomly selecting a subset of %d permutations\n", maxPerms)
                % neutral permutation (1) + a random sample with replacement
                % from permutations 2 : nPerms, for a total of maxPerms
                ind = [1, randperm(nPerms - 1, maxPerms - 1) + 1];
                self.perms = self.perms(sort(ind), :);
            end
        end

        function disp(self)
            % textually display information about analysis
            %
            % analysis.disp()
            %
            % This method overrides Matlab's `disp`, so you can also use
            % `disp(analysis)` or simply `analysis` without semicolon to
            % get the same output.

            str = sprintf("  Analysis:");
            str = str + sprintf("\n    %d fold(s), %d session(s)", self.L, self.m);
            if ~isequaln(self.CA, self.CB)
                str = str + sprintf("\n    CA:       %d × %d, %d-dimensional", ...
                    size(self.CA), self.dimensionA);
                str = str + sprintf("\n    CB:       %d × %d, %d-dimensional", ...
                    size(self.CB), self.dimensionB);
                % check cross parameter-effect extracting matrix
                str = str + sprintf("\n    CA ↔ CB:  %d-dimensional, %g %% variance", ...
                    self.dimensionAB, self.normAB / self.dimensionAB * 100);
                % TODO variance check equivalent to permutation check?
            else
                str = str + sprintf("\n    CA = CB:  %d × %d, %d-dimensional", ...
                    size(self.CA), self.dimensionA);
            end
            if size(self.perms, 1) > 1
                str = str + sprintf("\n    %d permutations (including neutral)", size(self.perms, 1));
            else
                str = str + sprintf("\n    no permutations");
            end
            disp(str)
        end

        function fig = show(self)
            % graphically display information about analysis
            %
            % fig = analysis.show()
            %
            % The method creates a figure and returns the handle.

            % create figure
            fig = figure();

            % size and layout of figure
            showPerms = size(self.perms, 1) > 1;
            if ~showPerms
                nRows = 2;
                fig.Position(3 : 4) = [640, 640];
            else
                nRows = 3;
                fig.Position(3 : 4) = [640, 960];
            end

            % colormap red → white → green
            cmap = interp1(0 : 2, [1 0 0 ; 1 1 1 ; 0 1 0], linspace(0, 2, 511));
            colormap(cmap)

            % determine scaling of contrasts
            maxC = max(abs([self.CA(:) ; self.CB(:)]));

            % contrasts
            subplot(nRows, 2, 1)
            heatmap(self.CA, ColorMap=cmap, ColorbarVisible=false)
            clim([-maxC, maxC])
            title('Contrast A')
            xlabel('subcontrasts')
            ylabel('regressors')
            subplot(nRows, 2, 2)
            heatmap(self.CB, ColorMap=cmap, ColorbarVisible=false)
            clim([-maxC maxC])
            title('Contrast B')
            xlabel('subcontrasts')
            ylabel('regressors')

            % sessions
            subplot(nRows, 2, 3)
            heatmap(double(self.sessionsA), ColorMap=cmap, ColorbarVisible=false)
            clim([0 1])
            title('Sessions A')
            xlabel('sessions')
            ylabel('folds')
            subplot(nRows, 2, 4)
            heatmap(double(self.sessionsB), ColorMap=cmap, ColorbarVisible=false)
            clim([0 1])
            title('Sessions B')
            xlabel('sessions')
            ylabel('folds')

            % permutations
            if showPerms
                subplot(nRows, 2, 5:6)
                hm = heatmap(self.perms, ColorMap=cmap, ColorbarVisible=false);
                % no row labels
                hm.YDisplayLabels = repmat({''}, size(hm.YDisplayLabels, 1), 1);
                clim([-1 1])
                title('Sign Permutations')
                xlabel('sessions')
                ylabel('permutations')
            end

            if nargout == 0
                clear fig
            end
        end
    end

    methods (Static)

        function analysis = leaveOneSessionOut(m, CA, CB)
            % create Analysis object for leave-one-session-out cross-validation
            %
            % analysis = Analysis.leaveOneSessionOut(m, C)
            % analysis = Analysis.leaveOneSessionOut(m, CA, CB)
            %
            % This is a convenience method which calls the constructor with
            % `sessionsB = logical(eye(m))` and `sessionsA = not(sessionsB)`.
            %
            % If only one contrast `C` is specified, it is used for both
            % 'training' (`CA`) and 'validation' (`CB`).

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
