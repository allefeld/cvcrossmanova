function mdTbl = tableToMarkdown(matTbl)

% format table object as Markdown table
%
% tableToMarkdown(matTbl)
% mdTbl = tableToMarkdown(matTbl)

% special space characters used for padding
puncSpace = char(0x2008);   % U+2008 PUNCTUATION SPACE
figSpace = char(0x2007);    % U+2007 FIGURE SPACE

% initialize Markdown table
pipes = string(repelem('|', size(matTbl, 1) + 2)');
mdTbl = pipes;

% for each variable
for i = 1 : size(matTbl, 2)
    % convert variable data to string array
    if isnumeric(matTbl.(i))
        % Numeric variables are formatted to strings, aligned on the
        % decimal point, and padded to the same length. Alignment and
        % padding use special space characters to ensure digits align no
        % matter the table alignment.
        num = matTbl.(i);
        % convert numbers to strings
        str = compose("%g", num);
        % do strings contain a decimal point?
        ind = contains(str, ".");
        if any(ind)     % if there are any
            % append punctuation space to strings without decimal point
            str(~ind) = str(~ind) + puncSpace;
            % determine position of decimal point or punctuation space
            decPos = cell2mat(strfind(strrep(str, puncSpace, "."), "."));
            % determine necessary amount of left padding to align them
            leftPad = max(decPos) - decPos;
        else            % if there are none
            % determine length of strings
            len = strlength(str);
            % determine necessary amount of left padding to right-align them
            leftPad = max(len) - len;
        end
        % left-pad with figure spaces
        str = arrayfun(@(lp) string(repelem(figSpace, lp)), leftPad) + str;
        % right-pad with figure spaces to same length
        len = strlength(str);
        rightPad = max(len) - len;
        str = str + arrayfun(@(rp) string(repelem(figSpace, rp)), rightPad);
    else
        % Other variables are just converted.
        str = string(matTbl.(i));
    end
    % format as column
    % prepend escaped column name
    name = "`" + matTbl.Properties.VariableNames(i) + "`";
    col = [name ; str];
    % right-pad with regular spaces to same length
    len = strlength(col);
    colWidth = max(len);
    rightPad = colWidth - len;
    col = col + arrayfun(@(rp) string(repelem(' ', rp)), rightPad);
    % additional spaces left and right
    col = " " + col + " ";
    % create and insert divider for horizontally centered column
    div = ":" + string(repelem('-', colWidth)) + ":";
    col = [col(1) ; div ; col(2 : end)];
    % append to table
    mdTbl = mdTbl + col + pipes;
end

% print if not requested as output
if nargout == 0
    fprintf("%s\n", mdTbl)
    clear mdTbl
end

% Possible improvements:
% - Make column alignments configurable.
% - Make `compose` formats configurable.
% - Make special space characters configurable?
% - Extra processing for scientific format, NaN, +-Inf?
% - replace HYPHEN-MINUS by proper MINUS


% <!-- Copyright © 2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later -->
