function figureStyle(ax, tAll, tStart)
%FIGURESTYLE  Apply the single common figure style used throughout the thesis.
%
%   Every figure in the document goes through this function, so that fonts,
%   tick direction, grid weight and LaTeX interpreters are identical across
%   chapters. Call it LAST, after plotting, after the axis labels and after
%   the legend, because it overwrites those properties.
%
%   SYNTAX
%     figureStyle(ax)                  style only
%     figureStyle(ax, tAll)            style + time axis anchored at tAll(1)
%     figureStyle(ax, tAll, tStart)    style + time axis anchored at tStart
%
%   INPUTS
%     ax      Axes handle, vector of axes handles, or a figure handle (in
%             which case every axes in the figure is styled, top to bottom).
%             [] or omitted uses the current axes.
%     tAll    Either the full time vector or just a window [tmin tmax].
%     tStart  Time of the first camera frame, in seconds.
%
%   TIME AXIS HANDLING
%     Recordings do not start at t = 0: the first frame is taken at the real
%     clock time t_start read from the specifications file. MATLAB places
%     ticks at round multiples, so an axis running from 108 s to 221 s gets
%     its first labelled tick at 120, leaving the 108-120 stretch unlabelled
%     against the frame. A reader then interprets "120" as the beginning of
%     the record and, with a curve starting near y = 0, misreads the corner
%     as an origin at (0,0).
%
%     To prevent that, a tick is forced exactly at t_start, and any automatic
%     tick too close to it is dropped so the labels do not collide. A small
%     margin is added on both sides so that a nucleation line falling exactly
%     on t_start (a crystal already present in frame 1) is not hidden by the
%     axes box.
%
%     With stacked axes only the LAST one carries tick labels; the others
%     share the tick positions without text, which is the correct convention
%     for a vertical stack sharing an x axis.
%
%   NOTE
%     The tick interpreter is LaTeX, so axis labels must use LaTeX syntax:
%     '$\mu$m/s', '$\mathrm{d}A/\mathrm{d}t$', and so on.
%
%   See also NUCLEATIONLINES.

if nargin < 1 || isempty(ax), ax = gca; end
if nargin < 2, tAll   = []; end
if nargin < 3, tStart = []; end

% A figure handle means "style every axes it contains", ordered top to bottom.
if isa(ax, 'matlab.ui.Figure')
    ax = flipud(findall(ax, 'Type', 'axes'));
end
ax = ax(isvalid(ax));
if isempty(ax), return; end

%% ---------- time axis ----------
if ~isempty(tAll)
    t1 = min(tAll(:));
    t2 = max(tAll(:));
    if isempty(tStart), tStart = t1; end

    if isfinite(t1) && isfinite(t2) && t2 > t1
        span = t2 - t1;
        pad  = 0.030 * span;                    % ~3 % margin on each side

        set(ax, 'XLim', [t1 - pad, t2 + pad]);

        tk = get(ax(end), 'XTick');             % automatic ticks
        tk = tk(tk > t1 & tk < t2);
        tk = tk(abs(tk - tStart) > 0.09*span);  % drop ticks crowding t_start
        tk = unique([tStart, tk]);              % and force one at t_start

        lbl = cell(size(tk));
        for k = 1:numel(tk)
            if abs(tk(k) - round(tk(k))) < 1e-6
                lbl{k} = sprintf('%d', round(tk(k)));
            else
                lbl{k} = sprintf('%.1f', tk(k));
            end
        end

        set(ax, 'XTick', tk);
        for k = 1:numel(ax)
            if k == numel(ax)
                set(ax(k), 'XTickLabel', lbl);   % only the bottom axes labels
            else
                set(ax(k), 'XTickLabel', []);
            end
        end
    end
end

%% ---------- style ----------
for k = 1:numel(ax)
    a = ax(k);
    if ~isvalid(a), continue; end

    set(a, 'FontName',             'Times New Roman', ...
           'FontSize',             12, ...
           'TickLabelInterpreter', 'latex', ...
           'Box',                  'on', ...
           'LineWidth',            0.9, ...
           'TickDir',              'in', ...
           'TickLength',           [0.012 0.012], ...
           'XMinorTick',           'on', ...
           'YMinorTick',           'on', ...
           'Layer',                'top', ...
           'Color',                'w');

    % Discreet grid, so it never competes with the data.
    grid(a, 'on');
    set(a, 'GridLineStyle', '-', 'GridColor', [0.15 0.15 0.15], ...
           'GridAlpha', 0.12, 'MinorGridAlpha', 0.05);

    % Title and labels in LaTeX, slightly larger than the ticks.
    set(a.Title,  'Interpreter','latex', 'FontSize',13, 'FontWeight','normal');
    set(a.XLabel, 'Interpreter','latex', 'FontSize',13);
    set(a.YLabel, 'Interpreter','latex', 'FontSize',13);
    if isprop(a, 'ZLabel')
        set(a.ZLabel, 'Interpreter','latex', 'FontSize',13);
    end

    % Double y axes carry two YLabel objects; the right-hand one needs it too.
    if numel(a.YAxis) > 1
        try
            activeSide = a.YAxisLocation;
            yyaxis(a, 'right'); set(a.YLabel, 'Interpreter','latex', 'FontSize',13);
            yyaxis(a, 'left');  set(a.YLabel, 'Interpreter','latex', 'FontSize',13);
            yyaxis(a, activeSide);
        catch
            % Not a double-axis plot after all; nothing to do.
        end
    end

    lg = get(a, 'Legend');
    if ~isempty(lg) && isvalid(lg)
        set(lg, 'Interpreter','latex', 'FontSize',11, ...
                'Box','off', 'Color','none');
    end
end
end
