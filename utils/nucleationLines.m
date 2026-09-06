function nucleationLines(ax, tNuc, colors, labels)
%NUCLEATIONLINES  Draw thin vertical markers at the crystal nucleation times.
%
%   Every time-resolved figure in this project needs the nucleation instants
%   marked, because they separate the two regimes the thesis compares: the
%   flow before any crystal exists, and the flow driven by a growing crystal.
%   This function draws those markers in a deliberately unobtrusive way so
%   they never compete visually with the measured curves.
%
%   The lines are thin, translucent, drawn behind the data, and excluded from
%   the legend.
%
%   SYNTAX
%     nucleationLines(ax, tNuc)
%     nucleationLines(ax, tNuc, colors)
%     nucleationLines(ax, tNuc, colors, labels)
%
%   INPUTS
%     ax      Target axes. Pass [] to use the current axes.
%     tNuc    Vector of nucleation times in seconds. Non-finite entries are
%             ignored, so an unmeasured crystal can be passed as NaN.
%     colors  [] for a neutral dark grey, or an Nx3 RGB matrix with one row
%             per line. Rows are reused cyclically if fewer than numel(tNuc).
%     labels  [] or false  -> no text
%             true         -> automatic numbering N1, N2, ...
%             cellstr      -> custom labels, one per line
%
%   See also FIGURESTYLE, XLINE.

if nargin < 1 || isempty(ax), ax = gca; end
if nargin < 3, colors = []; end
if nargin < 4, labels = []; end

tNuc = tNuc(:)';
tNuc = tNuc(isfinite(tNuc));
if isempty(tNuc), return; end

% Resolve the label argument into either a cellstr or an empty value.
if islogical(labels) && isscalar(labels)
    if labels
        labels = arrayfun(@(k) sprintf('N%d', k), 1:numel(tNuc), ...
                          'UniformOutput', false);
    else
        labels = [];
    end
elseif ischar(labels) || isstring(labels)
    labels = cellstr(labels);
end

for k = 1:numel(tNuc)
    if isempty(colors)
        col = [0.25 0.25 0.25];
    else
        col = colors(min(k, size(colors,1)), :);
    end

    xl = xline(ax, tNuc(k), '-');
    xl.Color            = col;
    xl.LineWidth        = 0.4;      % deliberately hairline
    xl.Alpha            = 0.55;     % and translucent
    xl.HandleVisibility = 'off';    % keep it out of the legend

    if iscell(labels) && k <= numel(labels)
        xl.Label                    = labels{k};
        xl.Interpreter              = 'latex';
        xl.FontSize                 = 9;
        xl.FontName                 = 'Times New Roman';
        xl.LabelColor               = col;
        xl.LabelHorizontalAlignment = 'left';
        xl.LabelVerticalAlignment   = 'top';
        xl.LabelOrientation         = 'horizontal';
    end
end
end
