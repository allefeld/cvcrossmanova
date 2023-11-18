classdef ModeledData < handle & matlab.mixin.Scalar

    % A `ModeledData` object encapsulates data and design matrices for
    % multiple sessions along with supplementary information. Upon
    % creation, it estimates and stores GLM parameters and errors. Combined
    % with one or more `Analysis` objects, it is the basis for analyses
    % within a `CvCrossManova` object.

    properties
        Xs      % cell array of per-session data matrices, observations × independent variables
        Ys      % cell array of per-session design matrices, observations × dependent variables
        fs      % array of per-session residual degrees of freedom
        m       % number of sessions 
        ns      % array of per-session numbers of observations
        p       % number of dependent variables
        qs      % array of per-session numbers of independent variables
        names   % cell array of per-session string arrays of names of independent variables
    end

    properties (Hidden = true)
        hBetas  % cell array of per-session parameter estimates
        hXis    % cell array of per-session error estimates
    end

    methods

        function obj = ModeledData(Ys, Xs, nvargs)
            % create `ModeledData` object
            %
            % modeledData = ModeledData(Ys, Xs, fs = ..., names = ...)
            %
            % `Ys` and `Xs` are cell arrays of per-session data and design
            % matrices, respectively.
            %
            % The optional `fs` is an array of per-session residual degrees
            % of freedom. Their default values are `size(Xs{k}, 1) -
            % rank(Xs{k})`, corresponding to the assumption that
            % observations are uncorrelated. If observations have been only
            % approximately decorrelated or have been filtered, correct
            % values should be explicitly specified.
            %
            % The optional `names` is a cell array of per-session string
            % arrays of names of independent variables. Default names are
            % identical in all sessions and of the form `"reg(#)"`, where
            % `#` is the index of the independent variable. If independent
            % variables with the same index in different sessions have
            % different meanings, correct names should be explicitly
            % specified.

            arguments
                Ys             (:, :)  cell
                Xs             (:, :)  cell
                nvargs.fs      (:, :)  double  = []
                nvargs.names   (:, :)  cell = {}
            end

            % store arguments
            obj.Ys = Ys(:) .';
            obj.Xs = Xs(:) .';
            obj.fs = nvargs.fs(:) .';
            obj.names = nvargs.names(:) .';

            % determine number of sessions
            obj.m = numel(obj.Ys);
            assert(obj.m == numel(obj.Xs), ...
                "Number of sessions of Ys and Xs must match.");

            % determine number of observations for each session
            obj.ns = cellfun(@(x) size(x, 1), obj.Ys);
            assert(isequal(obj.ns, cellfun(@(x) size(x, 1), obj.Xs)), ...
                "Number of rows of Ys and Xs must match in each session.");

            % determine number of dependent variables
            obj.p = size(obj.Ys{1}, 2);
            assert(all(cellfun(@(x) size(x, 2), obj.Ys) == obj.p), ...
                "Number of columns must match between sessions of Ys.");

            % determine numbers of independent variables
            obj.qs = cellfun(@(x) size(x, 2), obj.Xs);

            % if not specified, calculate residual degrees of freedom for each session
            if isempty(obj.fs)
                obj.fs = obj.ns - cellfun(@rank, obj.Xs);
            end
            assert(obj.m == numel(obj.fs), ...
                "Number of sessions of Ys and fs must match.");

            % if not specified, assign names to independent variables
            if isempty(obj.names)
                obj.names = arrayfun(@(q) compose("reg(%d)", 1 : q), ...
                    obj.qs, 'UniformOutput', false);
            end
            assert(obj.m == numel(obj.names), ...
                "Number of sessions of Ys and names must match.");
            assert(isequal(cellfun(@numel, obj.names), obj.qs), ...
                "Numbers of independent variables and names must match in each session.")

            % estimate GLM parameters and errors for each session
            obj.hBetas = cell(1, obj.m);
            obj.hXis = cell(1, obj.m);
            for k = 1 : obj.m
                obj.hBetas{k} = pinv(obj.Xs{k}) * obj.Ys{k};
                obj.hXis{k} = obj.Ys{k} - obj.Xs{k} * obj.hBetas{k};
            end
        end

        function str = disp(obj)
            % textually display information about the object
            %
            % modeledData.disp()
            %
            % This method overrides Matlab's `disp`, so you can also use
            % `disp(modeledData)` or simply `modeledData` without semicolon
            % to get the same output.

            % TODO make more informative using regressor names if available

            % prepare information string
            str = sprintf("  ModeledData:");
            ivns = strings(1, obj.m);
            for k = 1 : obj.m
                n = compose("%s", obj.names{k});
                if numel(n) > 4
                    n = n(1 : 4);
                    n(4) = '…';
                end
                ivns(k) = join(n, ", ");
            end
            cn = ["session", "n", "f", "q", "names"];
            arr = [
                compose("%d", 1 : obj.m)
                compose("%d", obj.ns)
                compose("%.1f", obj.fs)
                compose("%d", obj.qs)
                ivns
                ] .';
            ws = max(strlength([cn ; arr]));
            fmt = "\n    " + join(compose("%%-%ds", ws), "  ");
            str = str + sprintf(fmt, cn .');
            fmt = "\n    " + join(compose("%%%ds", ws), "  ");
            str = str + sprintf(fmt, arr .');
            str = str + sprintf("\n    p = %d", obj.p);
            % display information string, if not requested as output
            if nargout == 0
                disp(str)
                clear str
            end
        end

        function fig = show(obj, rescale)
            % graphically display information about the object
            %
            % fig = analysis.show(rescale = true)
            %
            % This method creates a figure showing the design matrices of
            % all sessions.
            %
            % The optional `rescale` specifies whether independent
            % variables are individually rescaled to an absolute maximum of
            % 1, to aid visibility of details.
            %
            % `fig` is the handle of the created figure.

            if nargin < 2
                rescale = true;
            end

            % create figure
            fig = figure();
            colormap(gray)

            % size and layout of figure
            fig.Position(3 : 4) = [640, 960];
            nCols = ceil(sqrt(obj.m));
            nRows = ceil(obj.m / nCols);

            % plot design matrices
            for k = 1 : obj.m
                subplot(nRows, nCols, k)
                X = obj.Xs{k};
                if rescale
                    maX = max(abs(X));
                    maX(maX < eps) = eps;   % avoid division by 0
                    X = X * diag(1 ./ maX);
                end
                imagesc(X)
                set(gca, 'XTick', 1 : obj.qs(k))
                set(gca, 'XTickLabel', obj.names{k})
            end
        end

    end

    methods (Static)

        function [md, misc] = fromSPM(modelDir, nvargs)
            % create `ModeledData` object from SPM data
            %
            % [modeledData, misc] = fromSPM(modelDir, regions = {}, wf = true)
            %
            % `modelDir` is the directory where the `SPM.mat` file
            % referring to an estimated model is located.
            %
            % The optional `regions` is a cell array of logical 3D volumes
            % or filenames specifying region masks. Without it, only the
            % SPM brain mask is applied.
            % 
            % The optional `wf` specifies whether to apply whitening and
            % high-pass filtering (set up in SPM) to data and design
            % matrices. It should usually be kept at its default value.
            %
            % Dependent variables (columns of data matrices `Ys`) are only
            % read from voxels within the SPM analysis brain mask after it
            % was intersected with the union of region masks (if provided).
            %
            % `misc` is a structure of fields containing additional
            % information:
            %
            % ------  ---------------------------------------------------------------------------------
            % `mask`  logical 3D volume indicating the voxels which dependent variables correspond to
            % `mat`   3D voxel indices to mm transformation matrix
            % `rmvi`  indices of dependent variables corresponding to each region, cell array of arrays
            % ------  ---------------------------------------------------------------------------------

            arguments
                modelDir        (1, :)  char
                nvargs.regions  (:, :)  cell     = {}
                nvargs.wf       (1, 1)  logical  = true
            end
    
            [Ys, Xs, fs, names, misc] = loadDataSPM(modelDir, ...
                nvargs.regions, nvargs.wf);
            md = ModeledData(Ys, Xs, fs=fs, names=names);
        end

    end

end


% Copyright © 2016–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later%d
