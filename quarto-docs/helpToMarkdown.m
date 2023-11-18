function helpToMarkdown(names)
% print Matlab help texts as Markdown for use in Quarto

for i = 1 : numel(names)
    name = string(names{i});
    if exist(name, "class")
        classHelpToMarkdown(name)
    else
        functionHelpToMarkdown(name)
    end
end


function functionHelpToMarkdown(name, cls, prefix)
% print Matlab function or method help texts as Markdown

% determine if method
isMethod = exist('cls', 'var');
% fix prefix
if ~exist('prefix', 'var')
    prefix = '';
end
% get help text
if ~isMethod
    helpText = getHelp(name);
else
    helpText = getHelp(sprintf("%s.%s", cls, name));
end
% split into paragraphs
helpText = strsplit(helpText, '\n\n');
% format second line as code
if numel(helpText) >= 2
    helpText{2} = sprintf('```matlab\n%s\n```', helpText{2});
end
% adjust class of warnings
for i = find(startsWith(helpText, ':::'))
    helpText{i} = strrep(helpText{i}, 'Warning', 'callout-warning');
end
% print markdown
if ~isMethod
    id = lower(name);
    fprintf("\n## [function]{.smallcaps} `%s` {#%s}\n\n", name, id)
else
    id = lower(cls) + "." + lower(name);
    if ~isequal(name, cls)
        fprintf("\n### [%s method]{.smallcaps} `%s` {#%s}\n\n", prefix, name, id)
    else
        fprintf("\n### [constructor]{.smallcaps} `%s` {#%s}\n\n", name, id)
    end
end
fprintf("%s\n\n", helpText{:})


function classHelpToMarkdown(name)
% print Matlab class help texts as Markdown

% get help text
helpText = getHelp(name);
% check whether it's just the constructor help
constructorHelp = getHelp(sprintf("%s.%s", name, name));
% if same, discard it
if isequal(helpText, constructorHelp)
    helpText = {};
else
    % split into paragraphs
    helpText = strsplit(helpText, '\n\n');
end
% print markdown
id = lower(name);
fprintf("\n## [class]{.smallcaps} `%s` {#%s}\n\n", name, id)
if numel(helpText) > 0
    fprintf("%s\n\n", helpText{:})
end

% get meta.class
% mc = eval(sprintf("?%s", name));
mc = meta.class.fromName(name);

% process properties
propertyVisible = ~[mc.PropertyList.Hidden];
propertyNames = {mc.PropertyList(propertyVisible).Name};
propertyDescriptions = {mc.PropertyList(propertyVisible).Description};
nameColumn = char(cellfun(@(x) ['`' x '`'], ...
    propertyNames, 'UniformOutput', false));
descriptionColumn = char(propertyDescriptions);
tableMarker = [repmat('-', 1, size(nameColumn, 2), 1), '  ', ...
    repmat('-', 1, size(descriptionColumn, 2), 1)];
table = cellstr([tableMarker ; ...
    nameColumn, repmat('  ', size(nameColumn, 1), 1), descriptionColumn ; ...
    tableMarker]);
% print Markdown
id = lower(name) + "-properties";
fprintf("\n### [properties]{.smallcaps} {#%s}\n\n", id)
fprintf("%s\n", table{:})
fprintf("\n\n")

% process methods
mcml = fliplr(reshape(mc.MethodList, 1, []));
% restrict to non-inherited methods
mcml = mcml(arrayfun(@(m) isequal(m.DefiningClass.Name, name), mcml));
% restrict to non-hidden methods
mcml = mcml(~[mcml.Hidden]);
% remove constructor
mcml = mcml(arrayfun(@(m) ~isequal(m.Name, name), mcml));
% process constructor
functionHelpToMarkdown(name, name)
% process static methods
for m = mcml([mcml.Static]);
    functionHelpToMarkdown(m.Name, name, 'static')
end
% process non-static methods
for m = mcml(~[mcml.Static]);
    functionHelpToMarkdown(m.Name, name)
end


function helpText = getHelp(name)
% cleaned-up version of help()

% get help text
helpText = help(name);
% split into lines & remove two spaces prepended to each line
helpText = strsplit(helpText, '\n');
for i = 1 : numel(helpText)
    helpText{i} = helpText{i}(3 : end);
end
% remove link to documentation browser
if numel(helpText) >= 3
    if strfind(helpText{end - 2}, "  Documentation for ")
        helpText = helpText(1 : end - 3);
    end
end
% % turn "see also" into links
% if numel(helpText) >= 2
%     line = helpText{end - 1};
%     ind = strfind(lower(line), 'see also ');
%     if ~isempty(ind)
%         terms = line(ind + 9 : end);
%         line = line(1 : ind + 8);
%         if terms(end) == '.'
%             terms = terms(1 : end - 1);
%         end
%         terms = split(terms, ', ');
%         for i = 1 : numel(terms)
%             terms{i} = sprintf('[%s](#%s)', terms{i}, lower(terms{i}));
%         end
%         terms = join(terms, ', ');
%         line = [line, terms{1}, '.'];
%         helpText{end - 1} = line;
%     end
% end
% rejoin lines
helpText = strjoin(helpText, '\n');

% Copyright © 2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later
