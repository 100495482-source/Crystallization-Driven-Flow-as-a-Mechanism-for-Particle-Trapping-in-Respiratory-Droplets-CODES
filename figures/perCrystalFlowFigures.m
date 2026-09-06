%% ========================================================================
%  PERCRYSTALFLOWFIGURES  Flow in each crystal's buffer zone.
%  ========================================================================
%
%  PRODUCES: FIGURES AND TABLES. Nothing is recomputed; the velocities come
%  from the .mat files that growthFlowComparison already wrote.
%
%  PREREQUISITE
%    growthFlowComparison must have been run for every droplet in the list,
%    so that compare_individual_<name>.mat exists.
%
%  WHAT THESE FIGURES SHOW
%    The flow measured inside each crystal's region of influence, plotted
%    against time since THAT crystal nucleated rather than against absolute
%    time. Aligning every crystal on its own nucleation is what makes them
%    comparable: crystals appear at different moments in different droplets,
%    so on an absolute axis a common post-nucleation behaviour would be
%    smeared out and invisible.
%
%    Both measurement methods are shown for every crystal, which doubles as
%    a validation: PIV and particle tracking are independent, and where they
%    agree the measurement is trustworthy.
%
%  READS  (per droplet)
%    compare_individual_<name>.mat  -> vFlowPIV_ind, vFlowTM_ind, tGrowth, tPIV
%    crystal_areas_<name>.mat       -> nucList, the nucleation frames
%
%  OUTPUTS  (in <BASE_FG>/flowmeasurements/crystals)
%    fig_crystal_PIV_all.png            all crystals, PIV, blue ramp
%    fig_crystal_TM_all.png             all crystals, TrackMate, orange ramp
%    fig_crystal_INDIV_<name>.png       per droplet, both methods together
%    table_crystal_comparison.csv       mean PIV, mean TM and their difference
%    table_crystal_comparison_latex.tex the same table, ready for the document
%
%  See also GROWTHFLOWCOMPARISON, CRYSTALGROWTHANALYSIS.

clear; close all; clc;

%% ========================================================================
%  <<< CONFIGURE HERE >>>
%  ========================================================================
DATA_ROOT = fullfile('C:', 'Data', 'FG');

sessions = { ...
    'TP2-B(2811)', 'S4' ; ...
    'TP3(1202)',   'S5' ; ...
    'TP3(1202)',   'S6' ; ...
    'TP5-E(0806)', 'S12'; ...
    'TP5-E(0806)', 'S13'; ...
    'TP5-E(0806)', 'S15'; ...
    'TP5-E(0806)', 'S18'; ...
    'TP5-E(0806)', 'S20'; ...
    'TP6-(2907)',  'S21'; ...
    'TP6-(2907)',  'S22'; ...
    'TP6-(2907)',  'S24'; ...
    'TP6-(2907)',  'S30'; ...
    'TP6-(2907)',  'S32'; ...
    'TP6-(2907)',  'S35'  ...
};

dropletDirs = cell(size(sessions,1),1);
for i = 1:size(sessions,1)
    dropletDirs{i} = fullfile(DATA_ROOT, sessions{i,1}, sessions{i,2});
end

outFolder = fullfile(DATA_ROOT, 'flowmeasurements', 'crystals');

% ---------- style ----------
FONT_SZ    = 15;
LW_PLOT    = 1.5;              % crystal lines
LW_NUC     = 0.9;              % nucleation marker
COLOR_PIV  = [0.15 0.35 0.80]; % reference blue   (PIV)
COLOR_TM   = [0.88 0.38 0.12]; % reference orange (TrackMate)
CLR_NUC    = [0.35 0.35 0.35]; % nucleation line
SMOOTH_WIN = 5;                % moving-mean window (frames)

% Line styles are a SECONDARY hint inside the per-droplet plots. The primary
% discriminator is the position in the colour ramp; the style cycle only
% helps when two adjacent shades are hard to tell apart.
LS_CYCLE   = {'-', '--', ':', '-.'};   % cycles every 4 crystals

if ~isfolder(outFolder), mkdir(outFolder); end

%% ========================================================================
%  LOAD EVERY CRYSTAL INTO ONE FLAT LIST
%  ========================================================================
%  Flattening across droplets is what allows the two global figures: every
%  crystal in the study on one axis, each aligned on its own nucleation.

ND = numel(dropletDirs);
allCrystals = struct([]);
Ntotal = 0;

for d = 1:ND
    dDir = dropletDirs{d};
    if ~isfolder(dDir)
        warning('Folder not found, skipping: %s', dDir); continue
    end

    specs  = readSpecs(dDir);
    name   = specs.nombre;
    resDir = fullfile(dDir, 'RESULTS_crystals');

    % ---------- rebuild frameTimes from the specifications ----------
    % The frame count comes from the thresholded image folder; the anchor
    % and the rate come from the specifications, as everywhere else.
    dd2 = dir(dDir); dd2 = dd2([dd2.isdir]);
    otsuDir = '';
    for i = 1:numel(dd2)
        if startsWith(dd2(i).name,'FotosOtsu','IgnoreCase',true)
            otsuDir = fullfile(dDir, dd2(i).name); break
        end
    end
    if ~isempty(otsuDir)
        ff = dir(fullfile(otsuDir,'*.tif'));
        if isempty(ff), ff = dir(fullfile(otsuDir,'*.png')); end
        nFr = numel(ff);
    else
        nFr = 600;   % fallback; only affects the axis length, not the data
    end
    t0 = 0;
    if isfield(specs,'t_start') && ~isnan(specs.t_start), t0 = specs.t_start; end
    frameTimes = t0 + (0:nFr-1) / specs.fps;

    % ---------- nucleation time of each crystal ----------
    cf = fullfile(resDir, ['crystal_areas_' name '.mat']);
    if ~isfile(cf)
        warning('No crystal_areas for %s, skipping.', name); continue
    end
    C_crys = load(cf, 'nucList');
    nucList = C_crys.nucList(:);
    NCrys   = numel(nucList);
    nucList = max(1, min(nucList, numel(frameTimes)));
    nucTimes = frameTimes(nucList)';

    % ---------- per-crystal buffer-zone velocities ----------
    indivFile = fullfile(resDir, ['compare_individual_' name '.mat']);
    if ~isfile(indivFile)
        warning(['No compare_individual for %s.\n' ...
                 '  -> Run growthFlowComparison first.'], name); continue
    end
    Cv = load(indivFile, 'vFlowPIV_ind', 'vFlowTM_ind', 'tGrowth', 'tPIV');

    vPIV_mat = Cv.vFlowPIV_ind;   % [nP  x NCrys]
    vTM_mat  = Cv.vFlowTM_ind;    % [nTM x NCrys]
    tPIV_abs = Cv.tPIV(:);        % absolute time, PIV
    tTM_abs  = Cv.tGrowth(:);     % absolute time, TrackMate and growth

    % A crystal with no usable flow data yields fewer columns than there are
    % crystals; padding keeps the indices aligned with nucList.
    if size(vPIV_mat,2) < NCrys
        vPIV_mat(:, end+1:NCrys) = NaN;
    end
    if size(vTM_mat,2) < NCrys
        vTM_mat(:, end+1:NCrys) = NaN;
    end

    % ---------- one struct entry per crystal ----------
    for ci = 1:NCrys
        Ntotal = Ntotal + 1;
        idx    = Ntotal;

        t_nuc_ci = nucTimes(ci);   % absolute nucleation time of this crystal

        % PIV: clip to the available frames, smooth, shift to relative time
        nP_avail = min(size(vPIV_mat,1), numel(tPIV_abs));
        vPIV_ci  = vPIV_mat(1:nP_avail, ci);
        vPIV_ci  = smoothdata(vPIV_ci,'movmean',SMOOTH_WIN,'omitnan');
        tPIV_rel = tPIV_abs(1:nP_avail) - t_nuc_ci;

        % TrackMate: identical treatment, so the two remain comparable
        nTM_avail = min(size(vTM_mat,1), numel(tTM_abs));
        vTM_ci    = vTM_mat(1:nTM_avail, ci);
        vTM_ci    = smoothdata(vTM_ci,'movmean',SMOOTH_WIN,'omitnan');
        tTM_rel   = tTM_abs(1:nTM_avail) - t_nuc_ci;

        allCrystals(idx).label     = sprintf('%sC%d', name, ci);
        allCrystals(idx).name      = name;
        allCrystals(idx).ci        = ci;
        allCrystals(idx).t_nuc     = t_nuc_ci;
        allCrystals(idx).tPIV_rel  = tPIV_rel(:);
        allCrystals(idx).vPIV      = vPIV_ci(:);
        allCrystals(idx).tTM_rel   = tTM_rel(:);
        allCrystals(idx).vTM       = vTM_ci(:);

        % Temporal means taken from nucleation onwards only. Including the
        % pre-nucleation flow would dilute exactly the effect being measured.
        mPIV = mean(vPIV_ci(tPIV_rel >= 0 & isfinite(vPIV_ci)), 'omitnan');
        mTM  = mean(vTM_ci( tTM_rel  >= 0 & isfinite(vTM_ci)),  'omitnan');
        allCrystals(idx).vPIV_mean = mPIV;
        allCrystals(idx).vTM_mean  = mTM;
    end

    fprintf('Loaded %-5s  %d crystal(s)\n', name, NCrys);
end

if Ntotal == 0
    error(['No crystal data loaded. Check the paths and whether ' ...
           'growthFlowComparison has been run.']);
end
fprintf('\nTotal crystals loaded: %d\n\n', Ntotal);

%% ========================================================================
%  COLOUR RAMPS
%  ========================================================================
%  Crystal k keeps the same position in both ramps, so the k-th blue and the
%  k-th orange always refer to the same crystal. Both ramps stay firmly
%  within their hue, so a curve is identifiable as PIV or TrackMate at a
%  glance regardless of which crystal it belongs to.

N = Ntotal;
pivCM = [linspace(0.05, 0.55, N)', ...   % deep blue -> sky blue
         linspace(0.20, 0.75, N)', ...
         linspace(0.65, 1.00, N)'];

tmCM  = [linspace(0.78, 1.00, N)', ...   % deep orange -> pale gold
         linspace(0.18, 0.65, N)', ...
         linspace(0.00, 0.18, N)'];

%% ========================================================================
%  FIGURE 1: PIV buffer velocity, every crystal
%  ========================================================================
fig1 = figure('Color','w','Position',[80 80 1300 680]);
ax1  = axes(fig1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
ax1.FontSize = FONT_SZ; ax1.TickLabelInterpreter = 'latex';

xline(ax1, 0, 'Color', CLR_NUC, 'LineWidth', LW_NUC);
text(ax1, 0.5, 0, 'Nucleation', 'FontSize',8,'Color',CLR_NUC, ...
     'Interpreter','latex','VerticalAlignment','bottom', ...
     'HorizontalAlignment','left');

for k = 1:Ntotal
    Ck = allCrystals(k);
    ok = isfinite(Ck.vPIV);
    if ~any(ok), continue; end

    plot(ax1, Ck.tPIV_rel, Ck.vPIV, '-', ...
         'Color', pivCM(k,:), 'LineWidth', LW_PLOT);

    % Label at the end of each line rather than in a legend: with 58 curves
    % a legend would be longer than the figure.
    idxE = find(ok, 1, 'last');
    text(ax1, Ck.tPIV_rel(idxE) + 1.5, Ck.vPIV(idxE), Ck.label, ...
         'FontSize', 7, 'Color', pivCM(k,:), ...
         'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
         'Interpreter','none');
end

xlabel(ax1, 'Time since crystal nucleation (s)', 'Interpreter','latex','FontSize',FONT_SZ);
ylabel(ax1, 'PIV velocity in buffer zone ($\mu$m\,s$^{-1}$)', ...
       'Interpreter','latex','FontSize',FONT_SZ);
title(ax1, 'PIV buffer-zone velocity --- all crystals, all droplets', ...
      'Interpreter','latex','FontSize',FONT_SZ+4,'FontWeight','normal');

exportgraphics(fig1, fullfile(outFolder,'fig_crystal_PIV_all.png'), 'Resolution',300);
fprintf('Saved: fig_crystal_PIV_all.png\n');

%% ========================================================================
%  FIGURE 2: TrackMate buffer speed, every crystal
%  ========================================================================
fig2 = figure('Color','w','Position',[80 80 1300 680]);
ax2  = axes(fig2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
ax2.FontSize = FONT_SZ; ax2.TickLabelInterpreter = 'latex';

xline(ax2, 0, 'Color', CLR_NUC, 'LineWidth', LW_NUC);
text(ax2, 0.5, 0, 'Nucleation', 'FontSize',8,'Color',CLR_NUC, ...
     'Interpreter','latex','VerticalAlignment','bottom', ...
     'HorizontalAlignment','left');

for k = 1:Ntotal
    Ck = allCrystals(k);
    ok = isfinite(Ck.vTM);
    if ~any(ok), continue; end

    plot(ax2, Ck.tTM_rel, Ck.vTM, '-', ...
         'Color', tmCM(k,:), 'LineWidth', LW_PLOT);

    idxE = find(ok, 1, 'last');
    text(ax2, Ck.tTM_rel(idxE) + 1.5, Ck.vTM(idxE), Ck.label, ...
         'FontSize', 7, 'Color', tmCM(k,:), ...
         'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
         'Interpreter','none');
end

xlabel(ax2, 'Time since crystal nucleation (s)', 'Interpreter','latex','FontSize',FONT_SZ);
ylabel(ax2, 'TrackMate speed in buffer zone ($\mu$m\,s$^{-1}$)', ...
       'Interpreter','latex','FontSize',FONT_SZ);
title(ax2, 'TrackMate buffer-zone speed --- all crystals, all droplets', ...
      'Interpreter','latex','FontSize',FONT_SZ+4,'FontWeight','normal');

exportgraphics(fig2, fullfile(outFolder,'fig_crystal_TM_all.png'), 'Resolution',300);
fprintf('Saved: fig_crystal_TM_all.png\n');

%% ========================================================================
%  FIGURE 3: one figure per droplet, both methods together
%  ========================================================================
%  CONVENTION
%    Line COLOUR  -> the method   (blue = PIV, orange = TrackMate)
%    Shade + style -> the crystal (C1 dark and solid, then lighter)
%
%  A crystal occupies the same position in both ramps and carries the same
%  line style in both methods, so a dashed blue and a dashed orange are the
%  same crystal measured two ways. That pairing is the point of the figure.

dropletNames = unique({allCrystals.name}, 'stable');

for di = 1:numel(dropletNames)
    dName   = dropletNames{di};
    idxList = find(strcmp({allCrystals.name}, dName));
    NCrys_d = numel(idxList);

    % Per-droplet ramps, so a droplet with two crystals gets the full
    % contrast rather than two nearly identical shades from a global ramp.
    nc = max(NCrys_d, 1);
    t_ramp = linspace(0, 1, nc)';
    pivRamp = [0.05 + 0.55*t_ramp, ...   % navy -> sky blue
               0.18 + 0.57*t_ramp, ...
               0.60 + 0.38*t_ramp];
    tmRamp  = [0.75 + 0.23*t_ramp, ...   % burnt orange -> pale gold
               0.20 + 0.45*t_ramp, ...
               0.00 + 0.18*t_ramp];

    figD = figure('Color','w','Position',[80 80 1200 650]);
    axD  = axes(figD); hold(axD,'on'); box(axD,'on'); grid(axD,'on');
    axD.FontSize = FONT_SZ; axD.TickLabelInterpreter = 'latex';

    xline(axD, 0, 'Color', CLR_NUC, 'LineWidth', LW_NUC);

    for ci = 1:NCrys_d
        k  = idxList(ci);
        Ck = allCrystals(k);

        pivClr = pivRamp(ci,:);
        tmClr  = tmRamp(ci,:);
        ls = LS_CYCLE{mod(ci-1, numel(LS_CYCLE)) + 1};

        if any(isfinite(Ck.vPIV))
            plot(axD, Ck.tPIV_rel, Ck.vPIV, ls, ...
                 'Color', pivClr, 'LineWidth', LW_PLOT + 0.2, ...
                 'DisplayName', sprintf('C%d PIV', ci));
        end

        if any(isfinite(Ck.vTM))
            plot(axD, Ck.tTM_rel, Ck.vTM, ls, ...
                 'Color', tmClr, 'LineWidth', LW_PLOT, ...
                 'DisplayName', sprintf('C%d TM', ci));
        end
    end

    nCols = min(3, ceil(NCrys_d / 3));
    legend(axD, 'Interpreter','none', 'Location','northeast', ...
           'FontSize', 9, 'NumColumns', nCols);

    % Key explaining the two-ramp convention, since it is not self-evident.
    text(axD, 0.01, 0.97, ...
         ['Shade position $\equiv$ crystal (C1 dark $\to$ C' num2str(NCrys_d) ...
          ' light)$\quad|\quad$Blue $=$ PIV$\quad|\quad$Orange $=$ TrackMate'], ...
         'Units','normalized','FontSize',8,'Interpreter','latex', ...
         'HorizontalAlignment','left','VerticalAlignment','top', ...
         'Color',[0.38 0.38 0.38]);

    xlabel(axD,'Time since crystal nucleation (s)', ...
           'Interpreter','latex','FontSize',FONT_SZ);
    ylabel(axD,'Buffer-zone velocity ($\mu$m\,s$^{-1}$)', ...
           'Interpreter','latex','FontSize',FONT_SZ);
    title(axD, sprintf('Droplet %s --- per-crystal buffer-zone velocity', dName), ...
          'Interpreter','none','FontSize',FONT_SZ+3,'FontWeight','normal');

    fname = sprintf('fig_crystal_INDIV_%s.png', dName);
    exportgraphics(figD, fullfile(outFolder, fname), 'Resolution',300);
    fprintf('Saved: %s\n', fname);
end

%% ========================================================================
%  TABLE: mean PIV vs TrackMate per crystal, with the difference
%  ========================================================================
%  Written twice: a CSV for further analysis, and a LaTeX table ready to be
%  included in the document. Writing the .tex here rather than transcribing
%  the numbers by hand removes a whole class of copying error.

% ---------- CSV ----------
fid = fopen(fullfile(outFolder,'table_crystal_comparison.csv'),'w');
fprintf(fid,'Crystal,Droplet,MeanPIV_um_s,MeanTM_um_s,Diff_PIV_minus_TM\n');
for k = 1:Ntotal
    Ck = allCrystals(k);
    fprintf(fid,'%s,%s,%.3f,%.3f,%.3f\n', ...
            Ck.label, Ck.name, Ck.vPIV_mean, Ck.vTM_mean, ...
            Ck.vPIV_mean - Ck.vTM_mean);
end
fclose(fid);

% ---------- LaTeX ----------
fid_tex = fopen(fullfile(outFolder,'table_crystal_comparison_latex.tex'),'w');
fprintf(fid_tex, ['\\begin{table}[htbp]\n\\centering\n' ...
    '\\caption{Per-crystal comparison of mean buffer-zone velocities from PIV ' ...
    'and TrackMate. $\\bar{v}$ is the temporal mean computed from each crystal''s ' ...
    'nucleation frame onwards. ' ...
    '$\\Delta = \\bar{v}_{\\mathrm{PIV}} - \\bar{v}_{\\mathrm{TM}}$; ' ...
    'positive values indicate PIV exceeds TrackMate. ' ...
    'Units: $\\mu$m\\,s$^{-1}$.}\n' ...
    '\\label{tab:crystal_comparison}\n']);
fprintf(fid_tex,'\\begin{tabular}{lccc}\\toprule\n');
fprintf(fid_tex,['Crystal & $\\bar{v}_{\\mathrm{PIV}}$ & ' ...
                 '$\\bar{v}_{\\mathrm{TM}}$ & $\\Delta$ \\\\\\midrule\n']);

prevName = '';
diffs_all = nan(Ntotal,1);
for k = 1:Ntotal
    Ck = allCrystals(k);
    % A rule between droplets, so the crystals of each one read as a block.
    if ~strcmp(Ck.name, prevName) && k > 1
        fprintf(fid_tex,'\\midrule\n');
    end
    prevName  = Ck.name;
    diff_k    = Ck.vPIV_mean - Ck.vTM_mean;
    diffs_all(k) = diff_k;
    fprintf(fid_tex,'%s & %.2f & %.2f & %+.2f \\\\\n', ...
            Ck.label, Ck.vPIV_mean, Ck.vTM_mean, diff_k);
end

fprintf(fid_tex,'\\midrule\n');
fprintf(fid_tex,['\\textbf{Mean (all)} & \\textbf{%.2f} & \\textbf{%.2f} & ' ...
                 '\\textbf{%+.2f} \\\\\n'], ...
        mean([allCrystals.vPIV_mean],'omitnan'), ...
        mean([allCrystals.vTM_mean], 'omitnan'), ...
        mean(diffs_all,'omitnan'));
fprintf(fid_tex,'\\bottomrule\n\\end{tabular}\n\\end{table}\n');
fclose(fid_tex);

fprintf('\nSaved: table_crystal_comparison.csv\n');
fprintf('Saved: table_crystal_comparison_latex.tex\n');
fprintf('\n>>> Overall mean PIV  (all crystals): %.3f um/s\n', ...
        mean([allCrystals.vPIV_mean],'omitnan'));
fprintf('>>> Overall mean TM   (all crystals): %.3f um/s\n', ...
        mean([allCrystals.vTM_mean], 'omitnan'));
fprintf('>>> Overall mean diff (PIV - TM):     %.3f um/s\n', ...
        mean(diffs_all,'omitnan'));
fprintf('\nAll outputs saved to: %s\n', outFolder);
