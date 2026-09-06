%% ========================================================================
%  INDIVIDUALGROWTHPANELS  Two-panel growth figure, one per droplet.
%  ========================================================================
%
%  PRODUCES: FIGURES ONLY. Nothing is recomputed and no data file is
%  written; everything is read from the .mat files that crystalGrowthAnalysis
%  already produced.
%
%  WHAT IT DRAWS
%    One figure per droplet, two panels side by side:
%      left   individual crystal area (um^2) vs time
%      right  individual radial growth velocity (um/s) vs time
%
%    Both panels put every crystal of the droplet on the same axes, in the
%    green colour range, so that crystals within a droplet can be compared
%    directly. The area panel shows how much material each crystal captured;
%    the velocity panel shows how fast its front advanced, which is the
%    quantity comparable with the flow measurements.
%
%  WHY EACH CURVE STARTS WHERE IT DOES
%    Every crystal is blanked to NaN before ITS OWN nucleation frame. A
%    crystal that nucleates late would otherwise appear to sit at zero area
%    from the beginning of the recording, implying it existed and was not
%    growing, which is not what was observed.
%
%  STYLE
%    Deliberately identical to crystalGrowthAnalysis: same green range, same
%    LW_PLOT, same time window, same nucleation markers. These figures sit
%    next to those in the document and any difference would read as meaning
%    something.
%
%  READS  (from each droplet's RESULTS_crystals)
%    growth_global_<n>.mat      -> frameTimes
%    growth_individual_<n>.mat  -> areaCrystal_um2
%    radial_velocity_<n>.mat    -> vRadial
%    plus the specifications, for the time axis and nucleation frames
%
%  OUTPUT
%    <outDir>/individual_growth_<droplet>.png
%
%  ONLY EDIT THE SECTIONS MARKED <<<
%
%  See also CRYSTALGROWTHANALYSIS, FIGURESTYLE, NUCLEATIONLINES.

clear; close all; clc;

%% <<<  DATA ROOT  >>>
DATA_ROOT = fullfile('C:', 'Data', 'FG');

%% <<<  DROPLETS  >>>
sessions = { ...
    'TP2-B(2811)', 'S4' ; ...
    'TP3(1202)',   'S5' ; ...
    'TP3(1202)',   'S6' ; ...
    'TP5-E(0806)', 'S12'; ...
    'TP5-E(0806)', 'S13'; ...
    'TP4-D(0106)', 'S15'; ...
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

%% <<<  OUTPUT FOLDER  >>>
outDir = fullfile(DATA_ROOT, 'individual_growth');

%% ------------------------------------------------------------------------
%  STYLE PARAMETERS, matching crystalGrowthAnalysis
%  ------------------------------------------------------------------------
LW_PLOT   = 1.6;   % line width, same as crystalGrowthAnalysis
FONT_SZ   = 13;    % axis labels
TITLE_SZ  = 17;    % panel titles
SUPTIT_SZ = 19;    % overall title (droplet name)
TICK_SZ   = 9;     % tick labels
LEG_SZ    = 13;    % legend

if ~isfolder(outDir), mkdir(outDir); end

%% ------------------------------------------------------------------------
%  MAIN LOOP: one figure per droplet
%  ------------------------------------------------------------------------
for iD = 1:numel(dropletDirs)
    dDir = dropletDirs{iD};
    [~, folderName] = fileparts(dDir);

    if ~isfolder(dDir)
        warning('Folder not found: %s', dDir);
        continue
    end

    resF = fullfile(dDir, 'RESULTS_crystals');

    % ---------- recover the droplet name from the saved results ----------
    % The folder name and the name used in the .mat filenames do not always
    % match (some droplets carry a suffix), so the name is taken from the
    % file itself rather than assumed.
    name = folderName;
    hit  = dir(fullfile(resF, 'growth_global_*.mat'));
    if isempty(hit)
        warning('No growth_global_*.mat in %s - skipping.', resF);
        continue
    end
    tk = extractBetween(hit(1).name, 'growth_global_', '.mat');
    if ~isempty(tk), name = tk{1}; end
    dropletName = name;
    fprintf('Processing %s ...\n', dropletName);

    % ---------- frameTimes, first from the .mat ----------
    ggFile = fullfile(resF, ['growth_global_' name '.mat']);
    if ~isfile(ggFile)
        warning('Not found: %s - skipping.', ggFile);
        continue
    end
    gg         = load(ggFile, 'frameTimes');
    frameTimes = gg.frameTimes(:)';
    NF         = numel(frameTimes);

    % ---------- then overwritten from the specifications ----------
    % The specifications are authoritative for the time axis. The stored
    % frameTimes is only a fallback for droplets whose specs cannot be read.
    nucList = 1;
    try
        [~, sp] = evalc('readSpecs(dDir)');
        tsR = 0;
        if ~isnan(sp.t_start), tsR = sp.t_start; end
        frameTimes = tsR + (0:NF-1) / sp.fps;

        if isfield(sp,'nucleation_frames') && ~isempty(sp.nucleation_frames) ...
                                           && any(~isnan(sp.nucleation_frames))
            nfR     = sort(sp.nucleation_frames(~isnan(sp.nucleation_frames)));
            nucList = min(max(round(nfR(:)), 1), NF);
        elseif ~isnan(sp.nucleation_frame)
            nucList = min(max(round(sp.nucleation_frame), 1), NF);
        end
    catch MEr
        warning('%s: readSpecs failed (%s). Using frameTimes from the .mat.', ...
                name, MEr.message);
    end

    NCrys  = numel(nucList);
    nucT   = frameTimes(nucList);
    tPlot0 = nucT(1);
    tWin   = [tPlot0, frameTimes(end)];

    % ---------- green colour range, as everywhere else ----------
    crystalColors = [linspace(0.05, 0.55, NCrys)', ...
                     linspace(0.45, 0.80, NCrys)', ...
                     linspace(0.22, 0.45, NCrys)'];

    % ---------- individual areas ----------
    giFile = fullfile(resF, ['growth_individual_' name '.mat']);
    if ~isfile(giFile)
        warning('Not found: %s - skipping.', giFile);
        continue
    end
    gi = load(giFile, 'areaCrystal_um2');
    if ~isfield(gi, 'areaCrystal_um2')
        warning('%s: growth_individual has no areaCrystal_um2.', name);
        continue
    end
    areaCrystal_um2 = double(gi.areaCrystal_um2);
    if size(areaCrystal_um2,1) ~= NF, areaCrystal_um2 = areaCrystal_um2'; end
    NCrys_a = min(NCrys, size(areaCrystal_um2,2));
    for ci = 1:NCrys_a
        f0 = nucList(ci);
        if f0 > 1, areaCrystal_um2(1:f0-1, ci) = NaN; end   % before it existed
    end

    % ---------- radial velocity ----------
    vrFile = fullfile(resF, ['radial_velocity_' name '.mat']);
    if ~isfile(vrFile)
        warning('Not found: %s - skipping.', vrFile);
        continue
    end
    vr      = load(vrFile, 'vRadial');
    vRadial = double(vr.vRadial);
    if size(vRadial,1) ~= NF, vRadial = vRadial'; end
    NCrys_v = min(NCrys, size(vRadial,2));
    for ci = 1:NCrys_v
        f0 = nucList(ci);
        if f0 > 1, vRadial(1:f0-1, ci) = NaN; end
    end

    % ====================================================================
    %  FIGURE: flat two-panel layout, one shared legend underneath
    % ====================================================================
    figIG = figure('Name',  ['Individual growth - ' dropletName], ...
                   'Color', 'w', ...
                   'Position', [100 100 1900 430]);

    tl = tiledlayout(figIG, 1, 2, 'TileSpacing','compact','Padding','compact');

    title(tl, ['Droplet ' regexp(dropletName, '\d+', 'match', 'once')], ...
          'Interpreter','latex', ...
          'FontSize',   SUPTIT_SZ, ...
          'FontWeight', 'bold');

    % ---------- panel 1: individual area ----------
    ax1 = nexttile(tl, 1);
    hold(ax1, 'on');
    for c = 1:min(NCrys_a, NCrys)
        plot(ax1, frameTimes, areaCrystal_um2(:,c), ...
            'LineWidth', LW_PLOT, ...
            'Color',     crystalColors(c,:), ...
            'DisplayName', sprintf('Crystal %d', c));
    end
    grid(ax1, 'on');
    nucleationLines(ax1, nucT, crystalColors, true);
    figureStyle(ax1, tWin, tPlot0);
    % Labels AFTER figureStyle, which would otherwise overwrite them.
    xlabel(ax1, 'Time (s)',                  'FontSize', FONT_SZ);
    ylabel(ax1, 'Crystal area ($\mu$m$^2$)', 'FontSize', FONT_SZ);
    title(ax1,  'Individual crystal area',   'FontSize', TITLE_SZ);

    % ---------- panel 2: radial growth velocity ----------
    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    for c = 1:min(NCrys_v, NCrys)
        plot(ax2, frameTimes, vRadial(:,c), ...
            'LineWidth', LW_PLOT, ...
            'Color',     crystalColors(c,:), ...
            'DisplayName', sprintf('Crystal %d', c));
    end
    grid(ax2, 'on');
    nucleationLines(ax2, nucT, crystalColors, true);
    figureStyle(ax2, tWin, tPlot0);
    xlabel(ax2, 'Time (s)',                          'FontSize', FONT_SZ);
    ylabel(ax2, 'Radial growth velocity ($\mu$m/s)', 'FontSize', FONT_SZ);
    title(ax2,  'Radial growth velocity per crystal', 'FontSize', TITLE_SZ);

    % ---------- one horizontal legend under both panels ----------
    lgd = legend(ax1, 'FontSize', LEG_SZ, ...
                 'Orientation', 'horizontal', ...
                 'Location',    'southoutside');
    lgd.Layout.Tile = 'south';

    outFile = fullfile(outDir, ['individual_growth_' dropletName '.png']);
    exportgraphics(figIG, outFile, 'Resolution', 300);
    fprintf('  Saved: %s\n', outFile);
    close(figIG);
end

fprintf('\nDone.\n');
