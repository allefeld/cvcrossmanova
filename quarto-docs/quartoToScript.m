function quartoToScript(filename)

% convert a Quarto `qmd` file to a script file,
% and add HTML for a download button

% get name of notebook corresponding to document
if nargin < 1
    notebooks = dir("*.ipynb");
    [~, ind] = max([notebooks.datenum]);
    filename = notebooks(ind).name;
end
[~, name, ~] = fileparts(filename);
nb = jsondecode(fileread(name + ".ipynb"));

% get relative directory
wd = pwd();
reldir = wd(strfind(wd, "quarto-docs") : end);

% transform notebook to m-script
fid = fopen(name + ".m", 'w');
for i = 1 : numel(nb.cells)
    cell_type = nb.cells{i}.cell_type;
    source = strip(string(nb.cells{i}.source));
    switch cell_type
    case 'raw'
        % raw: metadata, keep only title, add relative directory
        for j = 1 : numel(source)
            if startsWith(source(j), "title: ")
                fprintf(fid, "%% %s\n", source{j}(8 : end));
                fprintf(fid, "%% \n");
            end
        end
        fprintf(fid, "%% converted from %s\n", fullfile(reldir, name + ".qmd"));
        fprintf(fid, "\n");
    case 'markdown'
        % markdown: wrapped and commented
        for j = 1 : numel(source)
            lines = strip(string(textwrap(source(j), 78)), "right");
            fprintf(fid, "%% %s\n", lines);
        end
        fprintf(fid, "\n");
    case 'code'
        % code: as is, plus cell delimiter
        if ~any(startsWith(source, "%| echo: false"))
            fprintf(fid, "%s\n", source);
            fprintf(fid, "\n%%%%\n");
            fprintf(fid, "\n");
        end
    end
end
fclose(fid);

% create download link (needs `output: asis`)
lines = {
    '<div id="script">'
    '<div class="quarto-title-meta-heading">Matlab script</div>'
    '<div class="quarto-title-meta-contents">'
    sprintf('<p><a href="%s" download=><code>%s</code></a></p>', ...
        name + ".m", name + ".m")
    '</div>'
    '</div>'
    '<script>'
    'let script = document.getElementById("script");'
    'let meta = document.getElementsByClassName("quarto-title-meta");'
    'meta[0].append(script);'
    '</script>'
};
fprintf('%s\n', string(lines))


% <!-- Copyright © 2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later -->
