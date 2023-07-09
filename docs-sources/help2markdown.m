function help2markdown(names)

% There is a difference between
%   help('CvCrossManova')                   % class help
% and
%   help('CvCrossManova.CvCrossManova')     % constructor help
% Constructor help is substituted for class help if the latter is empty.
%  methods('CvCrossManova')
% lists the constructor as one of the methods.



for i = 1 : numel(names)
    name = names{i};
    isclass = exist(name, "class");

    if isclass
        label = "[class]{.smallcaps}";
    else
        label = "[function]{.smallcaps}";
    end

%     fprintf("\n## %s\n\n", name)
%     fprintf("%s\n\n", label)
    fprintf("\n## %s %s\n\n", label, name)

    % get help text
    helpText = help(name);
    % encode Unicode characters into entities
    helpText = unicodeToEntities(helpText);
    % split into lines & remove two spaces prepended to each line
    helpText = strsplit(helpText, '\n');
    for j = 1 : numel(helpText)
        helpText{j} = helpText{j}(3 : end);
    end
    % remove link to documentation browser
    if strfind(helpText{end - 2}, '  Documentation for ')
        helpText = helpText(1 : end - 3);
    end
    % rejoin lines and split into paragraphs
    helpText = strjoin(helpText, '\n');
    helpText = strsplit(helpText, '\n\n');
    % separate brief description and syntax from rest
    description = helpText{1};
    syntax = helpText{2};
    helpText = helpText(3 : end);

    fprintf("%s\n\n", description)
    fprintf("```matlab\n%s\n```\n\n", syntax)
    fprintf("%s\n\n", helpText{:})

    if isclass
        % report on properties and methods
    end

end


function str = unicodeToEntities(str)
% encode Unicode characters beyond ASCII 127 into character entities
% workaround necessitated by bug in matlab_kernel / metakernel,
% see https://github.com/Calysto/matlab_kernel/issues/155#issuecomment-1483906288
chars = unique(str(str > 127));
ents = arrayfun(@(x) ['&#x' dec2hex(x) ';'], chars, 'UniformOutput', false);
str = replace(str, num2cell(chars), ents);
