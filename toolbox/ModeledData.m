classdef ModeledData < handle & matlab.mixin.Scalar

    % object type containing data and design information

    properties
        Ys          % cell array of per-session design matrices, observations × variables
        Xs          % cell array of per-session data matrices, observations × regressors
        fs          % array of per-session residual degrees of freedom
        m           % number of sessions 
        ns          % array of per-session numbers of observations (rows)
    end

    properties (Hidden = true)
        nVariables  % number of data variables (columns)
        hBetas      % cell array of per-session parameter estimates
        hXis        % cell array of per-session error estimates
    end

    methods

        function obj = ModeledData(Ys, Xs, kwargs)
            % create `ModeledData` object
            %
            % modeledData = ModeledData(Ys, Xs, fs = ...)
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
                kwargs.fs      (:, :)  double  = []
            end

            % store arguments
            obj.Ys = Ys(:).';
            obj.Xs = Xs(:).';
            obj.fs = kwargs.fs(:).';

            % determine number of sessions
            obj.m = numel(obj.Ys);
            assert(obj.m == numel(obj.Xs), ...
                "Number of sessions of Ys and Xs must match.");

            % determine number of observations for each session
            obj.ns = cellfun(@(x) size(x, 1), obj.Ys);
            assert(isequal(obj.ns, cellfun(@(x) size(x, 1), obj.Xs)), ...
                "Number of rows of Ys and Xs must match in each session.");

            % if not specified, calculate residual degrees of freedom for each session
            if isempty(obj.fs)
                obj.fs = obj.ns - cellfun(@rank, obj.Xs);
            else
                assert(obj.m == numel(obj.fs), ...
                    "Number of sessions of Ys and fs must match.");
            end

            % determine number of variables
            obj.nVariables = size(obj.Ys{1}, 2);
            assert(all(cellfun(@(x) size(x, 2), obj.Ys) == obj.nVariables), ...
                "Number of columns must match between sessions of Ys.");

            % estimate GLM parameters and errors for each session
            obj.hBetas = cell(1, obj.m);
            obj.hXis = cell(1, obj.m);
            for k = 1 : obj.m
                obj.hBetas{k} = pinv(obj.Xs{k}) * obj.Ys{k};
                obj.hXis{k} = obj.Ys{k} - obj.Xs{k} * obj.hBetas{k};
            end
        end

        function disp(obj)
            % textually display information about the object
            %
            % modeledData.disp()
            %
            % This method overrides Matlab's `disp`, so you can also use
            % `disp(modeledData)` or simply modeledData without semicolon
            % to get the same output.

            % TODO make more informative using regressor names if available

            % prepare information string
            str = "  ModeledData:";
            sess = (1 : obj.m);
            ps = cellfun(@(x) size(x, 2), obj.Ys);
            qs = cellfun(@(x) size(x, 2), obj.Xs);
            tbl = table(sess(:), obj.ns(:), ps(:), qs(:), obj.fs(:), ...
                VariableNames=["session", 'n', 'p', 'q', 'f']);
            lines = splitlines(formattedDisplayText(tbl, ...
                SuppressMarkup=true));
            str = str + sprintf("\n%s", lines(1));
            for k = 1 : obj.m
                str = str + sprintf("\n%s", lines(k + 2));
            end
            % display information string
            disp(str)
        end

        % TODO show()

    end

    methods (Static)

        function [md, misc] = fromSPM(modelDir, kwargs)
            % create `ModeledData` object from SPM data
            %
            % [modeledData, misc] = fromSPM(modelDir, regions = {}, wf = true)
            %
            % ----------  -----------------------------------------------------------------------------
            % `modelDir`  directory where the `SPM.mat` file referring to an estimated model is located
            % `regions`   region mask(s), cell array of logical 3D volumes or filenames
            % `wf`        whether to apply the whitening & high-pass filtering set up by SPM 
            % ----------  -----------------------------------------------------------------------------
            %
            % Data variables (columns of data matrices `Ys`) are only
            % read from voxels within the SPM analysis brain mask, after it
            % was intersected with the union of region masks if present.
            %
            % ::: callout-warning
            % It is recommended to keep `wf` on its default value `true`.
            % While the pattern distinctness and pattern stability
            % estimators should still be unbiased because the degrees of
            % freedom are adjusted, not whitening may lead to loss of
            % statistical power and not filtering will retain low-frequency
            % signal components which may act as a confound for
            % experimental effects.
            % :::
            %
            % `misc` is a structure of fields containing additional information not
            % stored in `modeledData`:
            %
            % ------  --------------------------------------------------------------------------------------- 
            % `mask`  logical 3D volume indicating which voxels variables correspond to
            % `mat`   3D voxel indices to mm transformation matrix
            % `rmvi`  variable indices corresponding to voxels contained in each region, cell array of arrays
            % ------  ---------------------------------------------------------------------------------------

            arguments
                modelDir        (1, :)  char
                kwargs.regions  (:, :)  cell  = {}
                kwargs.wf       (1, 1)  logical = true
            end
    
            [Ys, Xs, fs, names, misc] = loadDataSPM(modelDir, ...
                kwargs.regions, kwargs.wf);
            md = ModeledData(Ys, Xs, fs=fs);
            % TODO add names to ModeledData
        end

    end

end


% Copyright © 2016–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
