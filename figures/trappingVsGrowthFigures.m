%% ========================================================================
%  TRAPPINGVSGROWTHFIGURES  Does a faster-growing crystal trap more particles?
%  ========================================================================
%
%  PRODUCES: FIGURES ONLY. No data file is written; both inputs are CSVs
%  produced by earlier stages.
%
%  THE QUESTION
%    The accumulation analysis showed that particles end up at the crystals.
%    This asks whether the effect scales: within a droplet, does the crystal
%    growing faster also trap a larger fraction of the particles that reach
%    it? If crystallisation-driven flow is the transport mechanism, it should.
%
%  WHY EVERYTHING IS CENTRED PER DROPLET
%    Droplets differ enormously in concentration, ambient conditions and
%    lifetime, and both growth rate and trapped fraction shift with all of
%    them. Comparing crystals across droplets would therefore measure the
%    difference between droplets, not between crystals. Subtracting each
%    droplet's own mean removes that and leaves the within-droplet
%    comparison, which is the one that bears on the mechanism.
%
%    The correlation is computed with the degrees of freedom corrected
%    accordingly: centring within groups costs one degree of freedom per
%    group, so df = n - nGroups - 1 rather than n - 2. Ignoring that would
%    overstate the significance.
%
%  INPUT FILES, MATCHED BY POSITION
%    Table2_percrystal_growth_vs_flow.csv   from growthVsFlowFigures
%    trap_percrystal_all.csv                from accum5_trappingPerCrystal
%
%    The two files walk the crystals in the same order, droplet by droplet,
%    crystal 1..N. They are cross-checked on the Crystal column and the row
%    count before anything is computed, because a silent misalignment would
%    pair each crystal with the wrong one and still produce a plausible plot.
%
%  FIGURES
%    trap_vs_growth_centred.png   global scatter, centred per droplet
%    trap_growth_opt1_rows.png    one row per droplet, one point per crystal
%    trap_growth_opt2_grid.png    a mini-panel per droplet
%    trap_growth_opt3_bars.png    Spearman rho per droplet
%
%    Four views of the same 58 crystals, from most aggregated to most
%    granular, so the reader can see whether the overall trend survives when
%    droplets are inspected one at a time.
%
%  STATISTICS
%    Pearson and Spearman are implemented locally, so this runs without the
%    Statistics Toolbox. Spearman is reported per droplet because with four
%    or fewer crystals a rank correlation is more honest than a fitted slope.
%
%  See also GROWTHVSFLOWFIGURES, ACCUM5_TRAPPINGPERCRYSTAL.

clear; close all; clc;

%% <<<  CONFIGURATION  >>>
DATA_ROOT = fullfile('C:', 'Data', 'FG');

growthFile = fullfile(DATA_ROOT,'growthflowcomparison', ...
                      'Table2_percrystal_growth_vs_flow.csv');
trapFile   = fullfile(DATA_ROOT,'Accumulation','trap_percrystal_all.csv');

outDir = fullfile(DATA_ROOT,'Interdroplet');
if ~exist(outDir,'dir'), mkdir(outDir); end

FONT_SZ = 12;
TICK_SZ = 10;

%% ---------- read and check ----------
G = readtable(growthFile);
G = G(~strcmpi(string(G.Droplet), 'MEAN (all crystals)'), :);   % drop the mean row

Tr = readtable(trapFile);

% Both assertions guard the same hazard: the files are joined by position,
% so any divergence in ordering must stop the analysis rather than silently
% pair the wrong crystals.
assert(height(G) == height(Tr), ...
       'The files do not have the same number of rows (%d vs %d).', ...
       height(G), height(Tr));
assert(all(G.Crystal == Tr.Crystal), ...
       'Crystal indices do not match: check the ordering of the files.');

droplet = string(G.Droplet);
growth  = G.MeanRadialGrowth;
pivLoc  = G.MeanPIV;
tmLoc   = G.MeanTM;
trapped = Tr.PercentTrapped;

%% ---------- centre within each droplet ----------
[uD, ~, gIdx] = unique(droplet, 'stable');
nG = numel(uD);

growth_c  = nan(size(growth));
trapped_c = nan(size(trapped));
piv_c     = nan(size(pivLoc));
tm_c      = nan(size(tmLoc));

for g = 1:nG
    s = (gIdx == g);
    growth_c(s)  = growth(s)  - mean(growth(s),  'omitnan');
    trapped_c(s) = trapped(s) - mean(trapped(s), 'omitnan');
    piv_c(s)     = pivLoc(s)  - mean(pivLoc(s),  'omitnan');
    tm_c(s)      = tmLoc(s)   - mean(tmLoc(s),   'omitnan');
end

%% ---------- correlations, with the corrected degrees of freedom ----------
fprintf('=== Correlations centred per droplet (n = %d crystals, %d droplets) ===\n', ...
        height(G), nG);
fprintf('%-28s %8s %10s\n', 'vs trapped fraction', 'r', 'p');
fprintf('%s\n', repmat('-',1,50));
pairs = {'growth', growth_c; 'local PIV', piv_c; 'local TrackMate', tm_c};
for k = 1:size(pairs,1)
    [r, p] = pearsonGrouped(pairs{k,2}, trapped_c, height(G), nG);
    fprintf('%-28s %8.3f %10.4f\n', pairs{k,1}, r, p);
end

%% ---------- rho per droplet, reused by the three views below ----------
rhoOrd = nan(nG,1);
for g = 1:nG
    s = (gIdx == g);
    if sum(s) >= 2
        rhoOrd(g) = spearmanRP(growth(s), trapped(s));
    end
end

fprintf('\n=== Spearman rho within each droplet ===\n');
for g = 1:nG
    s = (gIdx == g);
    if sum(s) < 2
        fprintf('  %-5s n=%d  rho = ---   (single crystal)\n', uD(g), sum(s));
    else
        [rho, p] = spearmanRP(growth(s), trapped(s));
        flag = '';
        if sum(s) < 4, flag = '  (n < 4, unreliable)'; end
        fprintf('  %-5s n=%d  rho = %6.3f  (p = %.3f)%s\n', ...
                uD(g), sum(s), rho, p, flag);
    end
end

% Common ordering for all three views: ascending rho, droplets without a rho
% at the end. Keeping the order identical lets the reader carry a droplet's
% position from one figure to the next.
key = rhoOrd;
key(isnan(key)) = -Inf;
[~, ordG] = sort(key);

%% ========================================================================
%  BASE FIGURE: global scatter, centred per droplet (58 crystals)
%  ========================================================================
fig0 = figure('Color','w','Position',[100 100 760 560]);
ax0  = gca; hold(ax0,'on');

PAL  = lines(nG);
hLeg = gobjects(nG,1);
for g = 1:nG
    s = (gIdx == g);
    hLeg(g) = scatter(ax0, growth_c(s), trapped_c(s), 110, PAL(g,:), 'filled', ...
                      'MarkerEdgeColor', PAL(g,:)*0.6, 'LineWidth', 0.8);
end

ok = isfinite(growth_c) & isfinite(trapped_c);
pf = polyfit(growth_c(ok), trapped_c(ok), 1);
xf = linspace(min(growth_c(ok)), max(growth_c(ok)), 50);
plot(ax0, xf, polyval(pf,xf), 'k--', 'LineWidth', 1.3);

[rAll, pAll] = pearsonGrouped(growth_c, trapped_c, height(G), nG);
text(ax0, 0.03, 0.96, ...
     sprintf('$r = %.2f$, $p = %.3f$, $n = %d$', rAll, pAll, sum(ok)), ...
     'Units','normalized', 'Interpreter','latex', 'FontSize', TICK_SZ, ...
     'VerticalAlignment','top', 'BackgroundColor',[1 1 1], 'Margin', 2);

xlabel(ax0, ['Mean radial growth rate, droplet mean removed ' ...
             '($\mu$m s$^{-1}$)'], 'Interpreter','latex', 'FontSize', FONT_SZ);
ylabel(ax0, 'Trapped fraction, droplet mean removed (\%)', ...
        'Interpreter','latex', 'FontSize', FONT_SZ);

grid(ax0,'on');
ax0.GridAlpha = 0.35;
set(ax0, 'TickLabelInterpreter','latex', 'FontSize', TICK_SZ, 'Box','on');

legend(ax0, hLeg, cellstr(uD), 'Interpreter','none', 'Location','eastoutside', ...
       'FontSize', TICK_SZ-1, 'Box','off', 'NumColumns',1);
hold(ax0,'off');

exportgraphics(fig0, fullfile(outDir,'trap_vs_growth_centred.png'), 'Resolution',300);
fprintf('\nSaved: trap_vs_growth_centred.png\n');

%% ========================================================================
%  VIEW 1: one row per droplet, one point per crystal
%   x = centred growth, colour = trapped fraction.
%   All 14 droplets and all 58 crystals are visible at once.
%  ========================================================================
figA = figure('Color','w','Position',[100 100 640 520]);
axA  = gca; hold(axA,'on');

yl = cell(nG,1);
for k = 1:nG
    g = ordG(k);
    s = (gIdx == g);
    scatter(axA, growth_c(s), repmat(k,sum(s),1), 120, trapped(s), 'filled', ...
            'MarkerEdgeColor',[0.3 0.3 0.3], 'LineWidth',0.6);
    yl{k} = sprintf('%s ($n=%d$)', char(uD(g)), sum(s));
end

plot(axA, [0 0], [0.4 nG+0.6], 'k--', 'LineWidth', 1);   % the droplet mean

colormap(axA, parula);
cb = colorbar(axA);
cb.Label.String       = 'Trapped fraction (\%)';
cb.Label.Interpreter  = 'latex';
cb.Label.FontSize     = FONT_SZ;
cb.TickLabelInterpreter = 'latex';

xlabel(axA, ['Mean radial growth rate, droplet mean removed ' ...
             '($\mu$m s$^{-1}$)'], 'Interpreter','latex', 'FontSize', FONT_SZ);
axA.XGrid = 'on';  axA.YGrid = 'off';  axA.GridAlpha = 0.35;
set(axA, 'YTick', 1:nG, 'YTickLabel', yl, ...
         'TickLabelInterpreter','latex', 'FontSize', TICK_SZ, ...
         'Box','on', 'YLim',[0.4 nG+0.6]);
hold(axA,'off');

exportgraphics(figA, fullfile(outDir,'trap_growth_opt1_rows.png'), 'Resolution',300);
fprintf('Saved: trap_growth_opt1_rows.png\n');

%% ========================================================================
%  VIEW 2: a mini-panel per droplet
%   Growth (uncentred) against trapping, on a common scale so the panels
%   are directly comparable.
%  ========================================================================
nCol = 5;
nRow = ceil(nG / nCol);

figB = figure('Color','w','Position',[100 100 1100 240*nRow]);
tl = tiledlayout(figB, nRow, nCol, 'TileSpacing','compact', 'Padding','compact');

xLim = [min(growth)-0.2, max(growth)+0.2];
yLim = [0, 100];

for k = 1:nG
    g = ordG(k);
    s = (gIdx == g);
    ax = nexttile(tl);
    hold(ax,'on');

    scatter(ax, growth(s), trapped(s), 90, [0.180 0.459 0.710], 'filled', ...
            'MarkerEdgeColor',[0.11 0.30 0.46], 'LineWidth',0.8);

    % A fitted line only where there are enough crystals for it to mean
    % something; with two points it would just join them.
    if sum(s) >= 3
        pfk = polyfit(growth(s), trapped(s), 1);
        xfk = linspace(min(growth(s)), max(growth(s)), 20);
        plot(ax, xfk, polyval(pfk,xfk), 'k--', 'LineWidth', 1.1);
    end

    if sum(s) >= 2
        txt = sprintf('%s: $\\rho = %.2f$', char(uD(g)), rhoOrd(g));
    else
        txt = sprintf('%s: $n = 1$', char(uD(g)));
    end
    title(ax, txt, 'Interpreter','latex', 'FontSize', TICK_SZ);

    grid(ax,'on'); ax.GridAlpha = 0.3;
    set(ax, 'TickLabelInterpreter','latex', 'FontSize', TICK_SZ-2, ...
            'Box','on', 'XLim',xLim, 'YLim',yLim);
    hold(ax,'off');
end

xlabel(tl, 'Mean radial growth rate ($\mu$m s$^{-1}$)', ...
       'Interpreter','latex', 'FontSize', FONT_SZ);
ylabel(tl, 'Trapped fraction (\%)', ...
       'Interpreter','latex', 'FontSize', FONT_SZ);

exportgraphics(figB, fullfile(outDir,'trap_growth_opt2_grid.png'), 'Resolution',300);
fprintf('Saved: trap_growth_opt2_grid.png\n');

%% ========================================================================
%  VIEW 3: rho bars for droplets with two or more crystals
%   Solid for n >= 4, hollow for n < 4. The single-crystal droplet is
%   omitted and mentioned in the figure caption instead.
%  ========================================================================
rhoV = []; nameV = {}; nV = []; okV = [];
for g = 1:nG
    s = (gIdx == g);
    if sum(s) < 2, continue; end
    rhoV(end+1)  = spearmanRP(growth(s), trapped(s)); %#ok<SAGROW>
    nameV{end+1} = char(uD(g));                       %#ok<SAGROW>
    nV(end+1)    = sum(s);                            %#ok<SAGROW>
    okV(end+1)   = (sum(s) >= 4);                     %#ok<SAGROW>
end
[rhoV, ord] = sort(rhoV);
nameV = nameV(ord);  nV = nV(ord);  okV = okV(ord);

ylabels = cell(1, numel(rhoV));
for k = 1:numel(rhoV)
    ylabels{k} = sprintf('%s ($n=%d$)', nameV{k}, nV(k));
end

figC = figure('Color','w','Position',[100 100 560 500]);
axC  = gca; hold(axC,'on');

POS = [0.180 0.459 0.710];
NEG = [0.839 0.153 0.157];
hSolid = []; hHollow = [];
for k = 1:numel(rhoV)
    c = POS;
    if rhoV(k) < 0, c = NEG; end
    if okV(k)
        h = barh(axC, k, rhoV(k), 0.62, 'FaceColor', c, ...
                 'EdgeColor', c*0.6, 'LineWidth', 0.9);
        if isempty(hSolid), hSolid = h; end
    else
        % Hollow: the value is shown but flagged as resting on too few points
        % to be taken at face value.
        h = barh(axC, k, rhoV(k), 0.62, 'FaceColor','none', ...
                 'EdgeColor', c, 'LineWidth', 1.2, 'LineStyle','--');
        if isempty(hHollow), hHollow = h; end
    end
end

plot(axC, [0 0], [0.4 numel(rhoV)+0.6], 'k-', 'LineWidth', 1);

xlabel(axC, ['Rank correlation between growth rate and ' ...
             'trapped fraction, $\rho$'], ...
       'Interpreter','latex', 'FontSize', FONT_SZ);
axC.XGrid = 'on';  axC.YGrid = 'off';  axC.GridAlpha = 0.35;
set(axC, 'YTick', 1:numel(rhoV), 'YTickLabel', ylabels, ...
         'TickLabelInterpreter','latex', 'FontSize', TICK_SZ, ...
         'Box','on', 'XLim',[-1.15 1.15], 'YLim',[0.4 numel(rhoV)+0.6]);

hh = [hSolid, hHollow];
ll = {'$n \geq 4$ crystals', '$n < 4$ crystals'};
hh = hh(isgraphics(hh));
legend(axC, hh, ll(1:numel(hh)), 'Interpreter','latex', ...
       'Location','southeast', 'FontSize', TICK_SZ-1, 'Box','off');
hold(axC,'off');

exportgraphics(figC, fullfile(outDir,'trap_growth_opt3_bars.png'), 'Resolution',300);
fprintf('Saved: trap_growth_opt3_bars.png\n');

fprintf('\nDone. Output in: %s\n', outDir);


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [r, p] = pearsonGrouped(x, y, nObs, nGroups)
%PEARSONGROUPED  Pearson r with the degrees of freedom corrected for grouping.
%
%   Centring within each group consumes one degree of freedom per group, so
%   df = nObs - nGroups - 1 rather than the usual nObs - 2. Using the
%   uncorrected value would inflate the significance of the correlation.

    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    C = corrcoef(x, y);
    r = C(1,2);
    df = nObs - nGroups - 1;
    if df < 1 || ~isfinite(r) || abs(r) >= 1
        p = NaN; return
    end
    t = r * sqrt(df / (1 - r^2));
    p = betainc(df / (df + t^2), df/2, 0.5);
end


function [rho, p] = spearmanRP(x, y)
%SPEARMANRP  Spearman rank correlation, without the Statistics Toolbox.
%   Ranks are computed with average ranks for ties, then correlated.
    ok = isfinite(x) & isfinite(y);
    x = rankTies(x(ok));  y = rankTies(y(ok));
    n = numel(x);
    if n < 2, rho = NaN; p = NaN; return; end
    C = corrcoef(x, y);
    rho = C(1,2);
    if n < 3 || ~isfinite(rho) || abs(rho) >= 1
        p = NaN; return
    end
    t = rho * sqrt((n-2) / (1 - rho^2));
    p = betainc((n-2) / ((n-2) + t^2), (n-2)/2, 0.5);
end


function rk = rankTies(x)
%RANKTIES  Ranks, averaging over ties.
    x = x(:);
    [xs, ord] = sort(x);
    n = numel(x);
    rk = zeros(n,1);
    i = 1;
    while i <= n
        j = i;
        while j < n && xs(j+1) == xs(i), j = j + 1; end
        rk(ord(i:j)) = (i + j) / 2;
        i = j + 1;
    end
end
