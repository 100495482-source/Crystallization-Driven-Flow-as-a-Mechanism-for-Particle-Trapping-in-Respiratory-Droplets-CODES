%% ========================================================================
%  GROWTHVSFLOWFIGURES  Crystal growth against flow, across all droplets.
%  ========================================================================
%
%  PRODUCES: FIGURES AND TABLES. Nothing is recomputed; every quantity is
%  read from the RESULTS_crystals folder of each droplet.
%
%  Table2_percrystal_growth_vs_flow.csv, written here, is consumed by
%  trappingVsGrowthFigures. Run this first.
%
%  WHAT IT DOES
%    Takes the per-droplet comparison that growthFlowComparison produced and
%    repeats it across the whole data set, at two levels:
%
%      GLOBAL      one figure per droplet, and one summary table with the
%                  mean flow before and after nucleation for both methods
%      PER CRYSTAL one panel per crystal, and one table row per crystal
%
%    The before/after split in the global table is the headline comparison:
%    if crystallisation drives the flow, the post-nucleation means should
%    exceed the pre-nucleation ones systematically, across droplets that
%    otherwise differ in every respect.
%
%  WHERE THE DATA COME FROM
%    crystals  : growth_global_<S>.mat   -> frontVel_perim
%                radial_velocity_<S>.mat -> vRadial, per crystal
%    PIV       : piv_fields_<S>.mat, compare_individual_<S>.mat
%    TrackMate : trackmate_meanspeed_*.mat, compare_individual_<S>.mat
%
%  TWO DROPLET LISTS, ON PURPOSE
%    The global analysis uses all 26 droplets. The per-crystal analysis uses
%    only the 14 whose crystals were individually segmented, since it needs
%    the per-crystal masks and buffers that only those have.
%
%  FIXED COLOURS: PIV blue, TrackMate orange, crystals green, as everywhere.
%
%  See also GROWTHFLOWCOMPARISON, TRAPPINGVSGROWTHFIGURES.

clear; close all; clc;

%% ---------- 1. DROPLETS, global analysis ----------
DATA_ROOT = fullfile('C:', 'Data', 'FG');

sessionsAll = { ...
    'TP2-B(2811)','S4' ; 'TP3(1202)','S5'  ; 'TP3(1202)','S6'  ; ...
    'TP4-D(0106)','S7' ; 'TP4-D(0106)','S8'; 'TP4-D(0106)','S9' ; ...
    'TP4-D(0106)','S10'; 'TP5-E(0806)','S12'; 'TP5-E(0806)','S13'; ...
    'TP5-E(0806)','S14'; 'TP5-E(0806)','S15'; 'TP5-E(0806)','S18'; ...
    'TP5-E(0806)','S19'; 'TP5-E(0806)','S20'; 'TP6-(2907)','S21' ; ...
    'TP6-(2907)','S22' ; 'TP6-(2907)','S23' ; 'TP6-(2907)','S24' ; ...
    'TP6-(2907)','S25' ; 'TP6-(2907)','S29' ; 'TP6-(2907)','S30' ; ...
    'TP6-(2907)','S31' ; 'TP6-(2907)','S32' ; 'TP6-(2907)','S33' ; ...
    'TP6-(2907)','S34' ; 'TP6-(2907)','S35'   ...
};

droplets = cell(size(sessionsAll,1),1);
for i = 1:size(sessionsAll,1)
    droplets{i} = fullfile(DATA_ROOT, sessionsAll{i,1}, sessionsAll{i,2});
end

%% ---------- 1b. DROPLETS, per-crystal analysis ----------
%  Only those with individual crystal segmentation available.
sessionsCrys = { ...
    'TP2-B(2811)','S4' ; 'TP3(1202)','S5'  ; 'TP3(1202)','S6'  ; ...
    'TP5-E(0806)','S12'; 'TP5-E(0806)','S13'; 'TP5-E(0806)','S15'; ...
    'TP5-E(0806)','S18'; 'TP5-E(0806)','S20'; 'TP6-(2907)','S21' ; ...
    'TP6-(2907)','S22' ; 'TP6-(2907)','S24' ; 'TP6-(2907)','S30' ; ...
    'TP6-(2907)','S32' ; 'TP6-(2907)','S35'   ...
};

dropletsCrys = cell(size(sessionsCrys,1),1);
for i = 1:size(sessionsCrys,1)
    dropletsCrys{i} = fullfile(DATA_ROOT, sessionsCrys{i,1}, sessionsCrys{i,2});
end

%% ---------- 2. OPTIONS ----------
USE_MASKED_TM = true;   % prefer the TrackMate data with crystal spots removed
SMOOTH_PIV    = 9;      % display only, for the global PIV curve; 1 = none
SMOOTH_TM     = 5;      % display only, for the global TrackMate curve
GROWTH_SCALE  = 1;      % 1 = true scale. Raise it only to make the front
                        % visible, and say so in the caption if you do.
TITLE_FS      = 20;
CLOSE_FIGS    = true;   % close each figure after saving; 26 droplets

colorPIV     = [0.20 0.40 0.80];   % blue
colorTM      = [0.85 0.40 0.20];   % orange
colorCrystal = [0.10 0.60 0.30];   % green

%% ---------- 3. OUTPUT FOLDER ----------
outDir = fullfile(DATA_ROOT, 'growthflowcomparison');
if ~isfolder(outDir), mkdir(outDir); end
fprintf('Results in: %s\n\n', outDir);

%% ---------- 4. LOAD EVERY DROPLET ----------
G = struct([]);
for i = 1:numel(droplets)
    G(i).d = loadGF(droplets{i}, USE_MASKED_TM);
    fprintf('  %-5s  PIV:%d  TM:%d  front:%d  crystals:%d\n', G(i).d.name, ...
        ~isempty(G(i).d.vPIV), ~isempty(G(i).d.vTM), ...
        ~isempty(G(i).d.vFront), size(G(i).d.vRad,2));
end

%% ---------- 5. GLOBAL FIGURE PER DROPLET: PIV + TM + FRONT ----------
for i = 1:numel(G)
    g = G(i).d;
    if isempty(g.vPIV) && isempty(g.vTM) && isempty(g.vFront)
        warning('%s: no global data', g.name); continue
    end

    fig = figure('Color','w','Position',[80 80 1300 560]); hold on
    h = gobjects(0); lab = {};

    if ~isempty(g.vPIV)
        h(end+1) = plot(g.tPIV, smoothSig(g.vPIV, SMOOTH_PIV), '-', ...
                        'LineWidth', 1.6, 'Color', colorPIV);
        lab{end+1} = 'PIV';
    end
    if ~isempty(g.vTM)
        h(end+1) = plot(g.tTM, smoothSig(g.vTM, SMOOTH_TM), '-', ...
                        'LineWidth', 1.6, 'Color', colorTM);
        lab{end+1} = 'TrackMate';
    end
    if ~isempty(g.vFront)
        h(end+1) = plot(g.tFront, GROWTH_SCALE*g.vFront, '-', ...
                        'LineWidth', 1.8, 'Color', colorCrystal);
        if GROWTH_SCALE == 1
            lab{end+1} = 'Crystal front';
        else
            % If the curve is scaled, the legend says so: an unlabelled
            % scaling would misrepresent the comparison.
            lab{end+1} = sprintf('Crystal front ($\\times$%g)', GROWTH_SCALE);
        end
    end

    % Nucleations, all drawn alike and labelled N1, N2, ...
    for k = 1:numel(g.tNuc)
        xline(g.tNuc(k), '-', sprintf('N%d', k), ...
              'Color',[0.55 0.55 0.55], 'LineWidth', 1.2, ...
              'Interpreter','latex', 'FontSize', 11, ...
              'LabelVerticalAlignment','top', ...
              'LabelHorizontalAlignment','center', 'HandleVisibility','off');
    end

    grid on; box on
    set(gca, 'TickLabelInterpreter','latex', 'FontSize',12);
    xlabel('Time (s)', 'Interpreter','latex', 'FontSize',14);
    ylabel('Velocity ($\mu$m s$^{-1}$)', 'Interpreter','latex', 'FontSize',14);
    title(sprintf('Droplet %s: crystal growth and flow', g.name), ...
          'Interpreter','latex', 'FontSize', TITLE_FS);
    legend(h, lab, 'Interpreter','latex', 'Location','northeast', 'FontSize',11);
    annotation(fig, 'textbox', [0.012 0.012 0.24 0.035], ...
        'String', 'N$k$ = nucleation of crystal $k$', 'Interpreter','latex', ...
        'FontSize', 11, 'EdgeColor',[0.75 0.75 0.75], 'BackgroundColor','w', ...
        'FitBoxToText','on', 'VerticalAlignment','middle');

    exportgraphics(fig, fullfile(outDir, ...
        sprintf('GLOBAL_%s_growth_vs_flow.png', g.name)), 'Resolution',300);
    if CLOSE_FIGS, close(fig); end
end

%% ---------- 6. GLOBAL TABLE, ONE ROW PER DROPLET ----------
%  The before/after columns are split at the FIRST nucleation, which is the
%  moment the droplet stops being a plain evaporating drop.
nm = {}; pAll=[]; pBef=[]; pAft=[]; tAll=[]; tBef=[]; tAft=[]; fAft=[];
for i = 1:numel(G)
    g = G(i).d;
    nm{end+1,1} = g.name;
    tN = g.tNuc; if isempty(tN), tN = NaN; end, tN = tN(1);
    pAll(end+1,1) = mn(g.vPIV);
    pBef(end+1,1) = mn(g.vPIV(g.tPIV <  tN));
    pAft(end+1,1) = mn(g.vPIV(g.tPIV >= tN));
    tAll(end+1,1) = mn(g.vTM);
    tBef(end+1,1) = mn(g.vTM(g.tTM <  tN));
    tAft(end+1,1) = mn(g.vTM(g.tTM >= tN));
    fAft(end+1,1) = mn(g.vFront);          % already NaN before nucleation
end

T1 = table(nm, pAll, pBef, pAft, tAll, tBef, tAft, fAft, ...
           pAft./fAft, tAft./fAft, pAft - tAft, ...
    'VariableNames', {'Droplet', ...
      'PIV_all','PIV_before','PIV_after', ...
      'TM_all','TM_before','TM_after', ...
      'FrontVel_after','Ratio_PIV_Front','Ratio_TM_Front','Diff_PIV_TM_after'});
T1 = [T1; [{'MEAN (all droplets)'}, ...
      num2cell(varfun(@(x) mean(x(isfinite(x))), T1(:,2:end), 'OutputFormat','uniform'))]];

writetable(T1, fullfile(outDir,'Table1_global_growth_vs_flow.csv'));
xlsPath = fullfile(outDir,'GrowthFlowComparison.xlsx');
if isfile(xlsPath), delete(xlsPath); end
writetable(T1, xlsPath, 'Sheet','Global');
disp(T1);

%% ---------- 7. PER CRYSTAL: THE THREE MEASUREMENTS ----------
GC = struct([]);
for i = 1:numel(dropletsCrys)
    GC(i).d = loadGF(dropletsCrys{i}, USE_MASKED_TM);
end

for i = 1:numel(GC)
    g = GC(i).d;
    nC = max([size(g.vRad,2) size(g.pInd,2) size(g.tInd,2)]);
    if nC == 0, warning('%s: no per-crystal data', g.name); continue; end

    nCol = min(3, nC); nRow = ceil(nC/nCol);
    fig = figure('Color','w','Position',[80 80 460*nCol 380*nRow]);
    tl = tiledlayout(fig, nRow, nCol, 'TileSpacing','compact', 'Padding','compact');

    for c = 1:nC
        ax = nexttile(tl); hold(ax,'on')
        % Time relative to THIS crystal's nucleation, so panels from
        % different crystals and droplets are directly comparable.
        tN = 0; if c <= numel(g.tNuc), tN = g.tNuc(c); end
        h = gobjects(0); lab = {};
        if c <= size(g.pInd,2)
            h(end+1) = plot(ax, g.tPind - tN, g.pInd(:,c), '-', ...
                            'LineWidth', 1.5, 'Color', colorPIV);
            lab{end+1} = 'PIV';
        end
        if c <= size(g.tInd,2)
            h(end+1) = plot(ax, g.tTind - tN, g.tInd(:,c), '-', ...
                            'LineWidth', 1.5, 'Color', colorTM);
            lab{end+1} = 'TrackMate';
        end
        if c <= size(g.vRad,2)
            h(end+1) = plot(ax, g.tRad - tN, GROWTH_SCALE*g.vRad(:,c), '-', ...
                            'LineWidth', 1.7, 'Color', colorCrystal);
            lab{end+1} = 'Radial growth';
        end
        xline(ax, 0, '-', 'Color',[0.55 0.55 0.55], 'LineWidth', 1.2);
        grid(ax,'on'); box(ax,'on')
        set(ax, 'TickLabelInterpreter','latex', 'FontSize',11);
        title(ax, sprintf('Crystal %d', c), 'Interpreter','latex', 'FontSize',16);
        if c == 1
            legend(ax, h, lab, 'Interpreter','latex', 'Location','northeast', ...
                   'FontSize',10);
        end
    end
    xlabel(tl, 'Time from crystal nucleation (s)', 'Interpreter','latex', 'FontSize',14);
    ylabel(tl, 'Velocity ($\mu$m s$^{-1}$)', 'Interpreter','latex', 'FontSize',14);
    title(tl, sprintf('Droplet %s: growth and flow per crystal', g.name), ...
          'Interpreter','latex', 'FontSize', TITLE_FS);
    exportgraphics(fig, fullfile(outDir, ...
        sprintf('PERCRYSTAL_%s_growth_vs_flow.png', g.name)), 'Resolution',300);
    if CLOSE_FIGS, close(fig); end
end

%% ---------- 8. PER-CRYSTAL TABLE ----------
%  This is the file trappingVsGrowthFigures reads. Its row order, droplet by
%  droplet and crystal 1..N, is what that script matches against.
dn={}; cn=[]; cP=[]; cT=[]; cR=[];
for i = 1:numel(GC)
    g = GC(i).d;
    nC = max([size(g.vRad,2) size(g.pInd,2) size(g.tInd,2)]);
    for c = 1:nC
        dn{end+1,1} = g.name;  cn(end+1,1) = c;
        if c <= size(g.pInd,2), cP(end+1,1) = mn(g.pInd(:,c)); else, cP(end+1,1) = NaN; end
        if c <= size(g.tInd,2), cT(end+1,1) = mn(g.tInd(:,c)); else, cT(end+1,1) = NaN; end
        if c <= size(g.vRad,2), cR(end+1,1) = mn(g.vRad(:,c)); else, cR(end+1,1) = NaN; end
    end
end
T2 = table(dn, cn, cP, cT, cR, cP-cT, cP./cR, cT./cR, ...
    'VariableNames', {'Droplet','Crystal','MeanPIV','MeanTM','MeanRadialGrowth', ...
                      'Diff_PIV_TM','Ratio_PIV_Growth','Ratio_TM_Growth'});
T2 = [T2; {'MEAN (all crystals)', NaN, mn(cP), mn(cT), mn(cR), ...
           mn(cP-cT), mn(cP./cR), mn(cT./cR)}];

writetable(T2, fullfile(outDir,'Table2_percrystal_growth_vs_flow.csv'));
writetable(T2, xlsPath, 'Sheet','PerCrystal');
disp(T2);

save(fullfile(outDir,'growthflow_all.mat'), 'G', 'GC', 'droplets', 'dropletsCrys');
fprintf('\nDone. Results in %s\n', outDir);


%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================

function y = smoothSig(v, w)
%SMOOTHSIG  Moving mean, or a pass-through when w <= 1.
if isempty(v), y = v; return; end
if w > 1, y = smoothdata(v, 'movmean', w, 'omitnan'); else, y = v; end
end


function m = mn(v)
%MN  Mean over the finite entries, NaN if there are none.
%   Used everywhere below so that a droplet missing one measurement yields a
%   blank cell rather than breaking the whole table.
v = v(isfinite(v));
if isempty(v), m = NaN; else, m = mean(v); end
end


function g = loadGF(dDir, useMasked)
%LOADGF  Collect every growth and flow series of one droplet.
%
%   Returns a struct with empty fields for whatever is not on disk, so the
%   caller can plot what exists instead of failing on the first gap. Files
%   are located by pattern rather than by exact name, because the droplet
%   name inside the filenames does not always equal the folder name.

[~, name] = fileparts(dDir);
g = struct('name',name, 'fps',NaN, 't0',0, 'tNuc',[], ...
           'tPIV',[], 'vPIV',[], 'tTM',[], 'vTM',[], ...
           'tFront',[], 'vFront',[], 'tRad',[], 'vRad',[], ...
           'tPind',[], 'pInd',[], 'tTind',[], 'tInd',[]);
resF = fullfile(dDir, 'RESULTS_crystals');
if ~isfolder(resF), warning('%s: no RESULTS_crystals', name); return; end

% evalc keeps the specifications printout out of the console: it would
% otherwise repeat 26 times and bury this function's own output.
evalc('specs = readSpecs(dDir);');
g.fps = specs.fps;
if isfield(specs,'t_start') && ~isnan(specs.t_start), g.t0 = specs.t_start; end

% ---------- nucleations, in crystal order ----------
% The specifications win; the stored nucList is only a fallback.
nucAll = [];
if isfield(specs,'nucleation_frames') && ~isempty(specs.nucleation_frames) ...
        && any(~isnan(specs.nucleation_frames))
    nucAll = sort(specs.nucleation_frames(~isnan(specs.nucleation_frames)));
else
    fca = dir(fullfile(resF,'crystal_areas_*.mat'));
    if ~isempty(fca)
        Sa = load(fullfile(resF, fca(1).name));
        if isfield(Sa,'nucList'), nucAll = sort(Sa.nucList); end
    end
    if isempty(nucAll) && ~isnan(specs.nucleation_frame)
        nucAll = specs.nucleation_frame;
    end
end
g.tNuc = g.t0 + (nucAll(:)' - 1)/g.fps;

% ---------- global PIV ----------
pf = dir(fullfile(resF,'piv_fields_*.mat'));
if ~isempty(pf)
    P = load(fullfile(resF, pf(1).name), 'velocity_magnitude_from_u_v_avg');
    Vc = P.velocity_magnitude_from_u_v_avg; nP = numel(Vc);
    vm = nan(nP,1);
    for k = 1:nP
        V = Vc{k}(:); V = V(isfinite(V));
        if ~isempty(V), vm(k) = mean(V); end
    end
    g.vPIV = vm;  g.tPIV = g.t0 + ((1:nP)' - 1)/g.fps;
end

% ---------- global TrackMate ----------
tmAll = dir(fullfile(resF,'trackmate_meanspeed_*.mat'));
if ~isempty(tmAll)
    isM = contains({tmAll.name},'masked');
    if useMasked && any(isM), idx = find(isM,1);
    elseif any(~isM),         idx = find(~isM,1);
    else,                     idx = 1;
    end
    S = load(fullfile(resF, tmAll(idx).name), 'framesU','meanSpeed');
    g.vTM = S.meanSpeed(:);
    % TrackMate counts frames from 0, hence framesU/fps without the -1.
    g.tTM = g.t0 + double(S.framesU(:))/g.fps;
end

% ---------- global front velocity ----------
gf = dir(fullfile(resF,'growth_global_*.mat'));
if ~isempty(gf)
    S = load(fullfile(resF, gf(1).name), 'frontVel_perim','frameTimes');
    if isfield(S,'frontVel_perim')
        g.vFront = S.frontVel_perim(:);
        g.tFront = S.frameTimes(:);
        n = min(numel(g.vFront), numel(g.tFront));
        g.vFront = g.vFront(1:n); g.tFront = g.tFront(1:n);
    end
end

% ---------- radial velocity per crystal ----------
rf = dir(fullfile(resF,'radial_velocity_*.mat'));
if ~isempty(rf)
    S = load(fullfile(resF, rf(1).name), 'vRadial','frameTimes');
    if isfield(S,'vRadial')
        g.vRad = S.vRadial;
        g.tRad = S.frameTimes(:);
        n = min(size(g.vRad,1), numel(g.tRad));
        g.vRad = g.vRad(1:n,:); g.tRad = g.tRad(1:n);
    end
end

% ---------- per-crystal flow inside the buffer ----------
cf = dir(fullfile(resF,'compare_individual_*.mat'));
if ~isempty(cf)
    S = load(fullfile(resF, cf(1).name), ...
             'vFlowPIV_ind','vFlowTM_ind','tGrowth','tPIV');
    if isfield(S,'vFlowPIV_ind') && ~isempty(S.vFlowPIV_ind)
        g.pInd  = S.vFlowPIV_ind;
        g.tPind = S.tPIV(:);
    end
    if isfield(S,'vFlowTM_ind') && ~isempty(S.vFlowTM_ind)
        g.tInd  = S.vFlowTM_ind;
        g.tTind = S.tGrowth(1:size(S.vFlowTM_ind,1)); g.tTind = g.tTind(:);
    end
end
end
