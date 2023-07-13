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
    helpText{2} = sprintf("```matlab\n%s\n```", helpText{2});
end
% print markdown
if ~isMethod
    fprintf("\n## [function]{.smallcaps} `%s`\n\n", name)
else
    if ~isequal(name, cls)
        fprintf("\n### [%s method]{.smallcaps} `%s`\n\n", prefix, name)
    else
        fprintf("\n### [constructor]{.smallcaps} `%s`\n\n", name)
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
fprintf("\n## [class]{.smallcaps} `%s`\n\n", name)
if numel(helpText) > 0
    fprintf("%s\n\n", helpText{:})
end

% get meta.class
mc = eval(sprintf("?%s", name));

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
table = cellfun(@unicodeEntities, table, 'UniformOutput', false);
% print Markdown
fprintf("\n### [properties]{.smallcaps}\n\n")
fprintf("%s\n", table{:})
fprintf("\n\n")

% process methods
mcml = flip(mc.MethodList);
for i = 1 : numel(mcml)
    notInherited = isequal(mcml(i).DefiningClass.Name, name);
    if ~mcml(i).Hidden && notInherited
        methodName = mcml(i).Name;
        if ~mcml(i).Static
            functionHelpToMarkdown(methodName, name)
        else
            functionHelpToMarkdown(methodName, name, 'static')
        end
    end
end


function helpText = getHelp(name)
% cleaned-up version of help()

% get help text
helpText = help(name);
helpText = unicodeEntities(helpText);
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
% rejoin lines
helpText = strjoin(helpText, '\n');

function helpText = unicodeEntities(helpText)
% encode Unicode characters into entities

chars = unique(helpText(helpText > 127));
ents = arrayfun(@(x) sprintf('&#x%s;', dec2hex(x)), chars, ...
    'UniformOutput', false);
helpText = replace(helpText, num2cell(chars), ents);
