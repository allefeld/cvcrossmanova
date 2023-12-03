function [res, ps] = runSearchlight(checkpoint, PSL, mask, ccm)

% apply object methods to data contained in a sliding spherical window
%
% [res, p] = runSearchlight(checkpoint, PSL, mask, ccm)
%
% checkpoint:   name for checkpoint file ([] to disable checkpointing)
% PSL:          prototype searchlight, 3-dimensional array
%               positive values indicate membership, value 2 indicates center
% mask:         3-dimensional logical array indicating which voxels to use
% ccm:          object providing `nResults` and `runAnalyses` methods
% res:          cell array with one cell for each analysis included in `cm`,
%               each containing an array of size mask voxels × results.
% p:            number of voxels in each searchlight, column vector
%
% Intermediate results are saved at regular intervals to the checkpoint
% file, if its name is given. If the checkpoint file already exists, its
% contents are loaded and the computation continues from that point.


% normalize checkpoint file name, preserving emptiness
if ~isempty(checkpoint)
    [~, ~, ext] = fileparts(checkpoint);
    if ~strcmp(ext, '.mat')
        checkpoint = [checkpoint '.mat'];
    end
end

assert(islogical(mask))

% searchlight voxels' subscripts relative to center
[s1, s2, s3] = ind2sub(size(PSL), find(PSL > 0));
[c1, c2, c3] = ind2sub(size(PSL), find(PSL == 2));
d = [s1, s2, s3] - [c1, c2, c3];

% volume dimensions
vdim = size(mask);
nMaskVoxels = nnz(mask);
factors = [1, cumprod(vdim(1 : end - 1))] .';

% mapping from volume to mask index
% mask index: index into extracted within-mask voxels
% can be accessed by 3 subscripts or single linear volume index
vol2mi = nan(vdim);
vol2mi(mask) = 1 : nMaskVoxels;

% initialize result volume(s)
% results are collected in arrays of size volume voxels × results
res = arrayfun(@(x) nan(nMaskVoxels, x), ccm.nResults(), ...
    'UniformOutput', false);
nAnalyses = size(res, 2);
ps = nan(nMaskVoxels, 1);

% for progress report
tic
t = 0;
line = '';
cmi = 0;
progressReport;
% iterate over subscripts of center (faster than `ind2sub`)
% `while` instead of `for` to support checkpointing
c3 = 0;
while c3 < vdim(3)                  % third volume index (slowest)
    c3 = c3 + 1;
    c2 = 0;
    while c2 < vdim(2)              % second volume index
        c2 = c2 + 1;
        c1 = 0;
        while c1 < vdim(1)          % first volume index (fastest)
            c1 = c1 + 1;

            if mask(c1, c2, c3)
                cmi = vol2mi(c1, c2, c3);      % center *mask* index

                % load checkpoint file at (re-) start and progress report 
                if cmi == 1
                    % if it exists
                    if exist(checkpoint, 'file')
                        load(checkpoint, 'c1', 'c2', 'c3', 'cmi', 'res', 'ps')
                        progressReport('loaded checkpoint file')
                    end
                end
    
                % save checkpoint file and progress report 
                nt = toc;
                if (nt - t > 30)    % every 30 s
                    t = nt;
                    if ~isempty(checkpoint) && (cmi < nMaskVoxels)
                        save(checkpoint, 'c1', 'c2', 'c3', 'cmi', 'res', 'ps')
                    end
                    progressReport()
                end

                % compute and store results
                %   determine subscripts of searchlight voxels
                s = d + [c1, c2, c3];
                %   discard out-of-volume voxels
                s = s(all((s >= [1, 1, 1]) & (s <= vdim), 2), :);
                %   volume indices of within-volume searchlight voxels
                vi = (s - 1) * factors + 1;   % faster than `sub2ind`
                %   discard out-of-mask voxels
                vi = vi(mask(vi));
                %   mask indices of within-mask searchlight voxels
                mi = vol2mi(vi);
                %   call CvCrossManova object method to run analyses
                Ds = ccm.runAnalyses(mi);
                % store results at center voxel
                for i = 1 : nAnalyses
                    res{i}(cmi, :) = Ds{i};
                end
                % store number of voxels within searchlight
                ps(cmi) = numel(mi);
            end

        end
    end
end
% last progress report
t = toc;
progressReport()
fprintf('\n')
% delete checkpoint file after completion
if ~isempty(checkpoint)
    spm_unlink(checkpoint)
end

function progressReport(message)
    fprintf(repmat('\b', 1, numel(line)))
    if nargin > 0
        fprintf('    %s\n', message)
        line = '';
    end
    line = sprintf('    %.1f min  %d voxels  %.1f %%', ...
        t / 60, cmi, cmi / nMaskVoxels * 100);
    fprintf('%s', line)
end

end


% Copyright © 2013–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
