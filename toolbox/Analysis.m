classdef Analysis

    properties
        CA         % 'training' contrast
        CB         % 'validation' contrast
        sessionsA  % 'training' sessions for each fold, L × m logical array
        sessionsB  % 'validation' sessions for each fold, L × m logical array

        same       % whether CA and CB are the same
        L          % number of folds
        m          % number of sessions
    end

    methods

        function self = Analysis(CA, CB, sessionsA, sessionsB)
            % object representing an analysis

            arguments
                CA (:,:) double
                CB (:,:) double
                sessionsA (:,:) logical
                sessionsB (:,:) logical
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
        end

        function checkEstimability(self, Xs)
            % TODO
        end

        % TODO add function to check whether cross or not based on
        % regressor ids

        function disp(self)
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
            arguments
                m  (1, 1) double
                CA (:, :) double
                CB (:, :) double = CA
            end

            sessionsB = logical(eye(m));
            sessionsA = not(sessionsB);
            analysis = Analysis(CA, CB, sessionsA, sessionsB);
        end

    end

end
