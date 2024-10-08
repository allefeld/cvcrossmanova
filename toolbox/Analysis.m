classdef Analysis < handle & matlab.mixin.Scalar

    % An `Analysis` object encapsulates 'training' and 'validation'
    % contrasts and sessions as well as permutations, which together
    % specify a Cross-validated (Cross-) MANOVA analysis.

    properties
        CA           % 'training' contrast matrix, regressors × subcontrasts
        CB           % 'validation' contrast matrix, regressors × subcontrasts
        sessionsA    % 'training' sessions logical matrix, folds × sessions
        sessionsB    % 'validation' sessions logical matrix, folds × sessions
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

        function obj = Analysis(CA, CB, sessionsA, sessionsB)
            % create `Analysis` object
            %
            % analysis = Analysis(CA, CB, sessionsA, sessionsB)
            %
            % `CA` and `CB` are 'training' and 'validation' contrast
            % matrices, respectively.
            %
            % `sessionsA` and `sessionsB` are logical matrices which
            % indicate which sessions (columns) are used for 'training' and
            % 'validation', respectively, in each fold (rows).

            arguments
                CA         (:, :)  double
                CB         (:, :)  double
                sessionsA  (:, :)  logical
                sessionsB  (:, :)  logical
            end

            % store arguments
            obj.CA = CA;
            obj.CB = CB;
            obj.sessionsA = sessionsA;
            obj.sessionsB = sessionsB;

            % TODO check whether sessionsA and sessionsB overlap in any fold
            % TODO check whether CB = P CA; isequal(sortrows(CA), sortrows(CB))

            % determine number of folds and sessions
            [obj.L, obj.m] = size(obj.sessionsA);
            % check
            assert(isequaln([obj.L, obj.m], size(obj.sessionsB)), ...
                "cvcrossmanova:sessionssize", ...
                "sessionsA and sessionsB must have the same size.")

            % characterize contrasts & check
            obj.dimensionA = rank(obj.CA);
            obj.dimensionB = rank(obj.CB);
            % check
            assert(size(CA, 2) == size(CB, 2), ...
                "cvcrossmanova:subcontrasts", ...
                "CA and CB must have the same number of subcontrasts (columns).")
            if obj.dimensionA ~= obj.dimensionB
                warning("cvcrossmanova:ranks", ...
                    "CA and CB should have the same rank.")
            end

            % characterize cross-contrast
            CBpCA = obj.CB * pinv(obj.CA);    % (cross-) parameter-effect extracting matrix
            obj.dimensionAB = rank(CBpCA);
            obj.normAB = sum(CBpCA(:) .^ 2);
            % check
            if obj.dimensionAB < min(obj.dimensionA, obj.dimensionB)
                warning("cvcrossmanova:preservedim", ...
                    "CA → CB should preserve dimensionality.")
            end
            if abs(obj.normAB / obj.dimensionAB - 1) > 1e-14
                warning("cvcrossmanova:preservevar", ...
                    "CA → CB should preserve variance.")
            end

            % initialize permutations to neutral only
            obj.perms = ones(1, obj.m);
        end

        function addPermutations(obj, nvargs)
            % add sign permutations of per-session parameter estimates
            %
            % analysis.addPermutations(maxPerms = 1000)
            %
            % This method adds information to `analysis` that sign
            % permutations should be applied, so that different values of
            % pattern distinctness *D* for the different permutations are
            % computed.
            %
            % The optional `maxPerms` specifies the maximum number of
            % permutations.

            % Two sign permutations are equivalent in a fold if the
            % relative signs applied to all involved sessions are the same.
            % Two sign permutations are equivalent if they are equivalent
            % in every fold.

            arguments
                obj
                nvargs.maxPerms  (1, 1)  double  = 1000
            end
            maxPerms = nvargs.maxPerms;

            assert(obj.m <= 21, "Too many possible sign permutations to process.")

            % number of sign permutations
            nPerms = 2 ^ obj.m;

            % generate all sign permutations of sessions (permutations × sessions)
            % with the neutral permutation first
            obj.perms = 1 - 2 * (int8(dec2bin(0 : nPerms - 1, obj.m)) - '0');

            pINA = [];
            % for each fold
            for l = 1 : obj.L
                % reduce permutations to Involved sessions
                involved = obj.sessionsA(l, :) + obj.sessionsB(l, :) > 0;
                pI = obj.perms(:, involved);
                % Normalize sign permutations by multiplying with the first value,
                % because only relative signs are important
                pIN = pI(:, 1) .* pI;
                % collect involved, normalized permutations across All folds
                pINA = [pINA, pIN];                                                     %#ok<AGROW>
            end
            % determine unique (non-equivalent) permutations
            [~, ind, ~] = unique(pINA, "rows");
            % select these permutations in order
            obj.perms = obj.perms(sort(ind), :);
            nPerms = size(obj.perms, 1);
            fprintf("%d permutations possible\n", nPerms)

            % Monte Carlo test
            if nPerms > maxPerms
                fprintf("randomly selecting a subset of %d permutations\n", maxPerms)
                % neutral permutation (1) + a random sample with replacement
                % from permutations 2 : nPerms, for a total of maxPerms
                ind = [1, randperm(nPerms - 1, maxPerms - 1) + 1];
                obj.perms = obj.perms(sort(ind), :);
            end
        end

        function str = disp(obj)
            % textually display information about the object
            %
            % analysis.disp()
            %
            % This method overrides Matlab's `disp`, so you can also use
            % `disp(analysis)` or simply `analysis` without semicolon to
            % get the same output.

            % prepare information string
            str = "  Analysis:";
            str = str + sprintf("\n    %d fold(s), %d session(s)", obj.L, obj.m);
            if ~isequal(obj.CA, obj.CB)
                str = str + sprintf("\n    CA:       %d × %d, %d-dimensional", ...
                    size(obj.CA), obj.dimensionA);
                str = str + sprintf("\n    CB:       %d × %d, %d-dimensional", ...
                    size(obj.CB), obj.dimensionB);
                % check cross parameter-effect extracting matrix
                str = str + sprintf("\n    CA ↔ CB:  %d-dimensional, %g %% variance", ...
                    obj.dimensionAB, obj.normAB / obj.dimensionAB * 100);
                % TODO variance check equivalent to permutation check?
            else
                str = str + sprintf("\n    CA = CB:  %d × %d, %d-dimensional", ...
                    size(obj.CA), obj.dimensionA);
            end
            if size(obj.perms, 1) > 1
                str = str + sprintf("\n    %d permutations (including neutral)", ...
                    size(obj.perms, 1));
            else
                str = str + sprintf("\n    no permutations (besides neutral)");
            end
            % display information string, if not requested as output
            if nargout == 0
                disp(str)
                clear str
            end
        end

        function fig = show(obj)
            % graphically display information about the object
            %
            % fig = analysis.show()
            %
            % This method creates a figure showing contrast matrices,
            % session matrices, and if present permutations.
            %
            % `fig` is the handle of the created figure.

            % create figure
            fig = figure(PaperPositionMode='auto');

            % size and layout of figure
            showPerms = size(obj.perms, 1) > 1;
            if ~showPerms
                nRows = 2;
                fig.Position(3 : 4) = [640, 640];
            else
                nRows = 3;
                fig.Position(3 : 4) = [640, 960];
            end

            % colormap red → white → green
            cmap = interp1(0 : 2, [0 0.686 1 ; 1 1 1 ; 1 0.686 0], linspace(0, 2, 511));

            % determine scaling of contrasts
            maxC = max(abs([obj.CA(:) ; obj.CB(:)]));

            % heatmap doesn't respect the default font
            fontname = get(0, 'defaultTextFontName');

            % contrasts
            subplot(nRows, 2, 1)
            heatmap(obj.CA, ...
                ColorMap=cmap, ColorbarVisible=false, FontName=fontname)
            clim([-maxC, maxC])
            title('Contrast A')
            xlabel('subcontrast')
            ylabel('regressor')
            subplot(nRows, 2, 2)
            heatmap(obj.CB, ...
                ColorMap=cmap, ColorbarVisible=false, FontName=fontname)
            clim([-maxC maxC])
            title('Contrast B')
            xlabel('subcontrast')
            ylabel('regressor')

            % sessions
            subplot(nRows, 2, 3)
            heatmap(double(obj.sessionsA), ...
                ColorMap=cmap, ColorbarVisible=false, FontName=fontname)
            clim([-1 1])
            title('Sessions A')
            xlabel('session')
            ylabel('fold')
            subplot(nRows, 2, 4)
            heatmap(double(obj.sessionsB), ...
                ColorMap=cmap, ColorbarVisible=false, FontName=fontname)
            clim([-1 1])
            title('Sessions B')
            xlabel('session')
            ylabel('fold')

            % permutations
            if showPerms
                subplot(nRows, 2, 5:6)
                imagesc(obj.perms .')
                colormap(cmap)
                clim([-1 1])
                title('Sign Permutations')
                xlabel('permutation')
                ylabel('session')
            end

            if nargout == 0
                clear fig
            end
        end
    end

    methods (Static)

        function analysis = leaveOneSessionOut(m, CA, CB)
            % create `Analysis` object for leave-one-session-out
            % cross-validation
            %
            % analysis = Analysis.leaveOneSessionOut(m, CA, CB)
            % analysis = Analysis.leaveOneSessionOut(m, C)
            %
            % `m` is the number of sessions. The object is created with
            % 'training' sessions `not(logical(eye(m)))` and 'validation'
            % sessions `logical(eye(m))`, the specification of standard
            % leave-one-session-out cross-validation.
            %
            % `CA` and `CB` are 'training' and 'validation' contrast
            % matrices, respectively. If only one contrast matrix `C` is
            % specified, it is used for both.

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

% Copyright © 2023-24 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
