%% ========================================================================
%  GLOBALFLOWFIGURES  Flow across every droplet, aligned at nucleation.
%  ========================================================================
%
%  PRODUCES: FIGURES AND TABLES. Nothing is recomputed; every quantity is
%  read from the RESULTS_crystals folder of each droplet.
%
%  WHAT IT SHOWS
%    The central result of the flow chapter: every droplet's mean flow speed
%    on one axis, with time measured from ITS OWN first nucleation. Droplets
%    nucleate at wildly different moments, so on an absolute axis a common
%    response to crystallisation would be invisible. Aligning them at t = 0
%    is what makes the pattern legible.
%
%    The stretch before nucleation is deliberately NOT trimmed. Seeing the
%    flow that exists beforehand is what distinguishes an increase caused by
%    the crystal from a flow that was already there.
%
%  SMOOTHING AND TABLES
%    Smoothing is applied ONLY to the plotted curves. Every table is computed
%    from the raw data, so no reported number depends on a display choice.
%
%  SECTIONS
%    5-6   global mean speed, PIV and TrackMate, all droplets aligned
%    7     tables: PIV, TrackMate, and the two compared
%    8     scatter of PIV against TrackMate, one point per droplet
%    9     per droplet: both methods on that droplet's own time axis
%    10-11 per crystal: buffer-zone flow, in panels and in detail
%    12    per-crystal table
%
%  WHY THE PIV vs TRACKMATE SCATTER MATTERS
%    The two methods are independent: PIV correlates image texture on a fixed
%    grid, tracking follows individual particles. Points near the 1:1 line
%    mean the measurement is a property of the flow rather than of the
%    method. The systematic offset that remains is discussed in the text.
%
%  READS  (per droplet)
%    piv_fields_<S>.mat                       spatial mean and max per frame,
%                                             plus mean vorticity and divergence
%    trackmate_meanspeed_[masked_]<S>.mat     meanSpeed, maxSpeed
%    compare_individual_<S>.mat               per-crystal buffer velocities
%
%  OUTPUT  (in <FG>/Flow measurements results)
%    PIV_globalmeanspeed.png, TRACKMATE_globalmeanspeed.png
%    SCATTER_PIV_vs_TrackMate_mean.png / _max.png
%    DROPLET_<S>_PIV_vs_TM.png, DROPLET_<S>_percrystal_PIV_TM.png
%    PIV_percrystal_set<k>.png, TRACKMATE_percrystal_set<k>.png
%    Table1_PIV.csv, Table2_TrackMate.csv, Table3_PIV_vs_TM.csv,
%    Table4_perCrystal.csv, FlowGlobal_PIV_TrackMate.xlsx
%    flow_all_droplets.mat
%
%  See also PIVFLOWANALYSIS, TRACKFLOWANALYSIS, GROWTHFLOWCOMPARISON.

clear; close all; clc;

%% ---------- 1. DROPLETS ----------
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

% Subset with individual crystal segmentation, used by sections 10 to 12.
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

% Sections 9 and 11 used to carry their own copies of these lists; they are
% the same sets, so they reuse them here.
dropletsSel    = droplets;
dropletsDetail = dropletsCrys;

%% ---------- 2. OPTIONS ----------
USE_MASKED_TM = true;    % TrackMate without the spots inside crystals
SMOOTH_PIV    = 9;       % display window for PIV; 1 = no smoothing
SMOOTH_TM     = 5;       % display window for TrackMate; 1 = no smoothing
XLIM_REL      = [];      % e.g. [-100 250]; [] = automatic
SHOW_LABELS   = true;    % droplet number at the end of each curve
TITLE_FS      = 18;
PANELS_PER_FIG = 4;      % droplets per multi-panel figure (2x2)

colorPIV = [0.20 0.40 0.80];   % blue
colorTM  = [0.85 0.40 0.20];   % orange

%% ---------- 3. OUTPUT FOLDER ----------
outDir = fullfile(DATA_ROOT, 'Flow measurements results');
if ~isfolder(outDir), mkdir(outDir); end
fprintf('Results in: %s\n\n', outDir);

%% ---------- 4. READ EVERY DROPLET ----------
nD = numel(droplets);
D  = struct([]);

for i = 1:nD
    dDir = droplets{i};
    [~, name] = fileparts(dDir);
    resF = fullfile(dDir, 'RESULTS_crystals');

    D(i).name = name;
    D(i).num  = regexprep(name, '\D', '');     % '4', '21', ...
    [D(i).tPIV, D(i).vPIV, D(i).mPIV, D(i).wPIV, D(i).dPIV] = deal([]);
    [D(i).tTM,  D(i).vTM,  D(i).mTM] = deal([]);
    D(i).fps = NaN; D(i).nucFrame = NaN;

    if ~isfolder(dDir), warning('%s does not exist', dDir); continue; end

    try
        evalc('specs = readSpecs(dDir);');
    catch ME
        warning('%s: no specifications (%s)', name, ME.message); continue
    end
    fps = specs.fps;

    % --- first nucleation (frame, 1-based) ---
    % Specifications first, stored nucList as a fallback. Without it the
    % droplet cannot be aligned and is skipped rather than plotted wrongly.
    nucFrame = NaN;
    if isfield(specs,'nucleation_frames') && ~isempty(specs.nucleation_frames) ...
            && any(~isnan(specs.nucleation_frames))
        nucFrame = min(specs.nucleation_frames(~isnan(specs.nucleation_frames)));
    elseif ~isnan(specs.nucleation_frame)
        nucFrame = specs.nucleation_frame;
    else
        fca = dir(fullfile(resF,'crystal_areas_*.mat'));
        if ~isempty(fca)
            S = load(fullfile(resF, fca(1).name));
            if isfield(S,'nucList') && ~isempty(S.nucList), nucFrame = min(S.nucList); end
        end
    end
    if isnan(nucFrame)
        warning('%s: no nucleation frame -> skipped', name); continue
    end
    D(i).fps = fps; D(i).nucFrame = nucFrame;

    % --- PIV: spatial mean and max per frame, plus vorticity and divergence ---
    pf = dir(fullfile(resF, 'piv_fields_*.mat'));
    if isempty(pf)
        fprintf('  %-5s  PIV: NO', name);
    else
        P = load(fullfile(resF, pf(1).name), ...
                 'velocity_magnitude_from_u_v_avg', ...
                 'vorticity_from_u_v_avg', 'divergence_from_u_v_avg');
        Vc = P.velocity_magnitude_from_u_v_avg;
        Wc = P.vorticity_from_u_v_avg;
        Gc = P.divergence_from_u_v_avg;
        nP = numel(Vc);
        [vm, mx, wm, gm] = deal(nan(nP,1));
        for k = 1:nP
            V = Vc{k}(:); V = V(isfinite(V));
            if ~isempty(V), vm(k) = mean(V); mx(k) = max(V); end
            if k <= numel(Wc) && ~isempty(Wc{k})
                W = Wc{k}(:); W = W(isfinite(W));
                if ~isempty(W), wm(k) = mean(W); end
            end
            if k <= numel(Gc) && ~isempty(Gc{k})
                G = Gc{k}(:); G = G(isfinite(G));
                if ~isempty(G), gm(k) = mean(G); end
            end
        end
        D(i).vPIV = vm;  D(i).mPIV = mx;  D(i).wPIV = wm;  D(i).dPIV = gm;
        D(i).tPIV = ((1:nP)' - nucFrame) / fps;   % nucleation at t = 0
        fprintf('  %-5s  PIV: %d frames', name, nP);
        clear P Vc Wc Gc     % these are large; free them before the next droplet
    end

    % --- TrackMate: the stored meanSpeed and maxSpeed ---
    tmAll = dir(fullfile(resF, 'trackmate_meanspeed_*.mat'));
    tmPath = '';
    if ~isempty(tmAll)
        isMasked = contains({tmAll.name}, 'masked');
        if USE_MASKED_TM && any(isMasked), idx = find(isMasked,1);
        elseif any(~isMasked),             idx = find(~isMasked,1);
        else,                              idx = 1;
        end
        tmPath = fullfile(resF, tmAll(idx).name);
    end
    if isempty(tmPath)
        fprintf('   |  TM: NO\n');
    else
        S  = load(tmPath, 'framesU', 'meanSpeed', 'maxSpeed');
        fr = double(S.framesU(:));               % TrackMate frame 0 = image 1
        D(i).vTM = S.meanSpeed(:);
        if isfield(S,'maxSpeed'), D(i).mTM = S.maxSpeed(:); else, D(i).mTM = nan(size(fr)); end
        D(i).tTM = (fr + 1 - nucFrame) / fps;    % hence the +1 before aligning
        fprintf('   |  TM: %d frames\n', numel(fr));
    end
end

%% ---------- 5. PIV, ALL DROPLETS ALIGNED ----------
fig1 = figure('Color','w','Position',[100 100 1500 560]); hold on
nPIV = 0;
for i = 1:numel(D)
    if isempty(D(i).vPIV), continue; end
    v = smoothSig(D(i).vPIV, SMOOTH_PIV);
    plot(D(i).tPIV, v, '-', 'LineWidth', 1.3, 'Color', colorPIV);
    % Every curve the same colour: this figure is about the shared shape,
    % not about telling droplets apart. The end label identifies them.
    if SHOW_LABELS, endLabel(D(i).tPIV, v, D(i).num, colorPIV); end
    nPIV = nPIV + 1;
end
xline(0, 'k--', 'LineWidth', 1.3);
grid on; box on
if ~isempty(XLIM_REL), xlim(XLIM_REL); end
set(gca, 'TickLabelInterpreter','latex', 'FontSize',12);
xlabel('Time relative to nucleation (s)', 'Interpreter','latex');
ylabel('Mean flow velocity $\langle |u| \rangle$ ($\mu$m/s)', 'Interpreter','latex');
title('Global mean speed from PIV, aligned at nucleation', ...
      'Interpreter','latex', 'FontSize', TITLE_FS);
noteLabels(nPIV);
exportgraphics(fig1, fullfile(outDir,'PIV_globalmeanspeed.png'), 'Resolution',300);

%% ---------- 6. TRACKMATE, ALL DROPLETS ALIGNED ----------
fig2 = figure('Color','w','Position',[100 100 1500 560]); hold on
nTM = 0;
for i = 1:numel(D)
    if isempty(D(i).vTM), continue; end
    v = smoothSig(D(i).vTM, SMOOTH_TM);
    plot(D(i).tTM, v, '-', 'LineWidth', 1.3, 'Color', colorTM);
    if SHOW_LABELS, endLabel(D(i).tTM, v, D(i).num, colorTM); end
    nTM = nTM + 1;
end
xline(0, 'k--', 'LineWidth', 1.3);
grid on; box on
if ~isempty(XLIM_REL), xlim(XLIM_REL); end
set(gca, 'TickLabelInterpreter','latex', 'FontSize',12);
xlabel('Time relative to nucleation (s)', 'Interpreter','latex');
ylabel('Mean particle speed $\bar{v}$ ($\mu$m/s)', 'Interpreter','latex');
title('Global mean speed from particle tracking, aligned at nucleation', ...
      'Interpreter','latex', 'FontSize', TITLE_FS);
noteLabels(nTM);
exportgraphics(fig2, fullfile(outDir,'TRACKMATE_globalmeanspeed.png'), 'Resolution',300);

%% ---------- 7. TABLES (raw data, unsmoothed) ----------
nm = {D.name}';

% --- SHEET 1: PIV ---
T1 = table(nm, ...
    col(D,'vPIV','tPIV','all'), col(D,'vPIV','tPIV','bef'), col(D,'vPIV','tPIV','aft'), ...
    col(D,'mPIV','tPIV','all'), col(D,'mPIV','tPIV','bef'), col(D,'mPIV','tPIV','aft'), ...
    col(D,'wPIV','tPIV','all'), col(D,'wPIV','tPIV','bef'), col(D,'wPIV','tPIV','aft'), ...
    col(D,'dPIV','tPIV','all'), col(D,'dPIV','tPIV','bef'), col(D,'dPIV','tPIV','aft'), ...
    'VariableNames', {'Droplet', ...
      'MeanSpeed','MeanSpeed_before','MeanSpeed_after', ...
      'MeanMaxSpeed','MeanMaxSpeed_before','MeanMaxSpeed_after', ...
      'MeanVorticity','MeanVorticity_before','MeanVorticity_after', ...
      'MeanDivergence','MeanDivergence_before','MeanDivergence_after'});

% --- SHEET 2: TrackMate ---
T2 = table(nm, ...
    col(D,'vTM','tTM','all'), col(D,'vTM','tTM','bef'), col(D,'vTM','tTM','aft'), ...
    col(D,'mTM','tTM','all'), col(D,'mTM','tTM','bef'), col(D,'mTM','tTM','aft'), ...
    'VariableNames', {'Droplet', ...
      'MeanSpeed','MeanSpeed_before','MeanSpeed_after', ...
      'MeanMaxSpeed','MeanMaxSpeed_before','MeanMaxSpeed_after'});

% --- SHEET 3: PIV vs TrackMate (difference = PIV - TM) ---
sP = [T1.MeanSpeed T1.MeanSpeed_before T1.MeanSpeed_after];
sT = [T2.MeanSpeed T2.MeanSpeed_before T2.MeanSpeed_after];
xP = [T1.MeanMaxSpeed T1.MeanMaxSpeed_before T1.MeanMaxSpeed_after];
xT = [T2.MeanMaxSpeed T2.MeanMaxSpeed_before T2.MeanMaxSpeed_after];
dS = sP - sT;   dX = xP - xT;

T3 = table(nm, sP(:,1), sT(:,1), dS(:,1), sP(:,2), sT(:,2), dS(:,2), ...
                sP(:,3), sT(:,3), dS(:,3), ...
                xP(:,1), xT(:,1), dX(:,1), xP(:,2), xT(:,2), dX(:,2), ...
                xP(:,3), xT(:,3), dX(:,3), ...
    'VariableNames', {'Droplet', ...
      'MeanSpeed_PIV','MeanSpeed_TM','Diff', ...
      'MeanSpeed_PIV_before','MeanSpeed_TM_before','Diff_before', ...
      'MeanSpeed_PIV_after','MeanSpeed_TM_after','Diff_after', ...
      'MaxSpeed_PIV','MaxSpeed_TM','DiffMax', ...
      'MaxSpeed_PIV_before','MaxSpeed_TM_before','DiffMax_before', ...
      'MaxSpeed_PIV_after','MaxSpeed_TM_after','DiffMax_after'});

% Summary rows: mean and standard deviation of every column across droplets.
agg = @(f) varfun(@(x) f(x(isfinite(x))), T3(:,2:end), 'OutputFormat','uniform');
T3 = [T3; [{'MEAN (all droplets)'}, num2cell(agg(@mean))]];
T3 = [T3; [{'STD  (all droplets)'},  num2cell(agg(@std))]];

writetable(T1, fullfile(outDir,'Table1_PIV.csv'));
writetable(T2, fullfile(outDir,'Table2_TrackMate.csv'));
writetable(T3, fullfile(outDir,'Table3_PIV_vs_TM.csv'));

xlsPath = fullfile(outDir,'FlowGlobal_PIV_TrackMate.xlsx');
if isfile(xlsPath), delete(xlsPath); end
writetable(T1, xlsPath, 'Sheet','PIV');
writetable(T2, xlsPath, 'Sheet','TrackMate');
writetable(T3, xlsPath, 'Sheet','PIV_vs_TM');

save(fullfile(outDir,'flow_all_droplets.mat'), 'D', 'droplets', 'USE_MASKED_TM');

disp(T3(:,{'Droplet','MeanSpeed_PIV','MeanSpeed_TM','Diff','Diff_before','Diff_after'}));
fprintf('\n%d droplets with PIV, %d with TrackMate.\nTables in %s\n', nPIV, nTM, outDir);

%% ---------- 8. SCATTER: PIV vs TRACKMATE, ONE POINT PER DROPLET ----------
scatterPIVTM(sP(:,1), sT(:,1), D, colorPIV, TITLE_FS, ...
    'Mean PIV velocity ($\mu$m s$^{-1}$)', ...
    'Mean TrackMate velocity ($\mu$m s$^{-1}$)', ...
    'PIV vs TrackMate, mean velocity per droplet', ...
    fullfile(outDir,'SCATTER_PIV_vs_TrackMate_mean.png'));

scatterPIVTM(xP(:,1), xT(:,1), D, colorPIV, TITLE_FS, ...
    'Mean maximum PIV velocity ($\mu$m s$^{-1}$)', ...
    'Mean maximum TrackMate velocity ($\mu$m s$^{-1}$)', ...
    'PIV vs TrackMate, maximum velocity per droplet', ...
    fullfile(outDir,'SCATTER_PIV_vs_TrackMate_max.png'));

%% ---------- 9. PER DROPLET, ON ITS OWN TIME AXIS ----------
%  Absolute time here, not aligned, with every nucleation marked. This is
%  where a droplet with several crystals can be inspected one at a time.
for i = 1:numel(dropletsSel)
    G = loadOne(dropletsSel{i}, USE_MASKED_TM);
    if isempty(G.vPIV) && isempty(G.vTM)
        warning('%s: no data', G.name); continue
    end

    figS = figure('Color','w','Position',[100 100 1300 560]); hold on

    h = gobjects(0); lab = {};
    if ~isempty(G.vPIV)
        h(end+1) = plot(G.tPIV, smoothSig(G.vPIV, SMOOTH_PIV), '-', ...
                        'LineWidth', 1.6, 'Color', colorPIV);
        lab{end+1} = 'PIV';
    end
    if ~isempty(G.vTM)
        h(end+1) = plot(G.tTM, smoothSig(G.vTM, SMOOTH_TM), '-', ...
                        'LineWidth', 1.6, 'Color', colorTM);
        lab{end+1} = 'TrackMate';
    end

    % All nucleations drawn alike and labelled N1, N2, ...
    for k = 1:numel(G.nucAll)
        xline(G.t0 + (G.nucAll(k)-1)/G.fps, '-', sprintf('N%d', k), ...
              'Color', [0.55 0.55 0.55], 'LineWidth', 1.2, ...
              'Interpreter','latex', 'FontSize', 11, ...
              'LabelVerticalAlignment','top', ...
              'LabelHorizontalAlignment','center', ...
              'HandleVisibility','off');
    end

    grid on; box on
    set(gca, 'TickLabelInterpreter','latex', 'FontSize',12);
    xlabel('Time (s)', 'Interpreter','latex');
    ylabel('Mean velocity ($\mu$m s$^{-1}$)', 'Interpreter','latex');
    title(sprintf('Droplet %s: PIV and TrackMate', G.name), ...
          'Interpreter','latex', 'FontSize', TITLE_FS);
    legend(h, lab, 'Interpreter','latex', 'Location','northeast');
    text(0.015, 0.94, 'N$k$ = nucleation of crystal $k$', ...
         'Units','normalized', 'Interpreter','latex', 'FontSize', 10);

    exportgraphics(figS, fullfile(outDir, ...
        sprintf('DROPLET_%s_PIV_vs_TM.png', G.name)), 'Resolution',300);
    fprintf('  %s: %d nucleations\n', G.name, numel(G.nucAll));
end

%% ---------- 10. PER-CRYSTAL BUFFER FLOW, IN PANELS ----------
C = struct([]);
for i = 1:numel(dropletsCrys)
    C(i).d = loadCrys(dropletsCrys{i});
    if isempty(C(i).d.vPIV) && isempty(C(i).d.vTM)
        warning('%s: compare_individual_*.mat is missing', C(i).d.name);
    end
end

crysPanels(C, 'PIV', colorPIV, PANELS_PER_FIG, outDir, ...
           'Mean PIV velocity in crystal buffer', 'PIV_percrystal');
crysPanels(C, 'TM',  colorTM,  PANELS_PER_FIG, outDir, ...
           'Mean TrackMate velocity in crystal buffer', 'TRACKMATE_percrystal');

%% ---------- 11. PER-CRYSTAL DETAIL: ONE PANEL PER CRYSTAL ----------
for q = 1:numel(dropletsDetail)
    G  = loadCrys(dropletsDetail{q});
    nC = max([size(G.vPIV,2) size(G.vTM,2)]);
    if nC == 0, warning('%s: no per-crystal data', G.name); continue; end

    nCol = min(3, nC); nRow = ceil(nC/nCol);
    figD = figure('Color','w','Position',[80 80 460*nCol 380*nRow]);
    tl = tiledlayout(figD, nRow, nCol, 'TileSpacing','compact', 'Padding','compact');

    for c = 1:nC
        ax = nexttile(tl); hold(ax,'on')
        h = gobjects(0); lab = {};
        % Each panel shifted to ITS crystal's nucleation, so panels can be
        % compared across crystals and droplets.
        if c <= size(G.vPIV,2)
            h(end+1) = plot(ax, G.tPIV - G.tNuc(c), G.vPIV(:,c), '-', ...
                            'LineWidth', 1.6, 'Color', colorPIV);
            lab{end+1} = 'PIV';
        end
        if c <= size(G.vTM,2)
            h(end+1) = plot(ax, G.tTM - G.tNuc(c), G.vTM(:,c), '-', ...
                            'LineWidth', 1.6, 'Color', colorTM);
            lab{end+1} = 'TrackMate';
        end
        xline(ax, 0, '-', 'Color',[0.55 0.55 0.55], 'LineWidth', 1.2);
        box(ax,'on'); grid(ax,'on')
        set(ax, 'TickLabelInterpreter','latex', 'FontSize',11);
        title(ax, sprintf('Crystal %d', c), 'Interpreter','latex', 'FontSize',16);
        if c == 1, legend(ax, h, lab, 'Interpreter','latex', ...
                          'Location','northeast', 'FontSize',10); end
    end
    xlabel(tl, 'Time from crystal nucleation (s)', 'Interpreter','latex', 'FontSize',14);
    ylabel(tl, 'Mean velocity in buffer ($\mu$m s$^{-1}$)', ...
           'Interpreter','latex', 'FontSize',14);
    title(tl, sprintf('Droplet %s: PIV and TrackMate per crystal', G.name), ...
          'Interpreter','latex', 'FontSize', 20);
    exportgraphics(figD, fullfile(outDir, ...
        sprintf('DROPLET_%s_percrystal_PIV_TM.png', G.name)), 'Resolution',300);
end

%% ---------- 12. PER-CRYSTAL TABLE ----------
dn = {}; cn = []; mP = []; mT = [];
for i = 1:numel(C)
    G = C(i).d;
    nC = max([size(G.vPIV,2) size(G.vTM,2)]);
    for c = 1:nC
        dn{end+1,1} = G.name;  cn(end+1,1) = c;
        if c <= size(G.vPIV,2), mP(end+1,1) = mean(G.vPIV(isfinite(G.vPIV(:,c)),c));
        else,                   mP(end+1,1) = NaN; end
        if c <= size(G.vTM,2),  mT(end+1,1) = mean(G.vTM(isfinite(G.vTM(:,c)),c));
        else,                   mT(end+1,1) = NaN; end
    end
end
T4 = table(dn, cn, mP, mT, mP - mT, 'VariableNames', ...
    {'Droplet','Crystal','MeanSpeed_PIV','MeanSpeed_TM','Diff'});
T4 = [T4; {'MEAN (all crystals)', NaN, mean(mP,'omitnan'), ...
           mean(mT,'omitnan'), mean(mP-mT,'omitnan')}];

writetable(T4, fullfile(outDir,'Table4_perCrystal.csv'));
if exist('xlsPath','var'), writetable(T4, xlsPath, 'Sheet','PerCrystal'); end
disp(T4);


%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================

function y = smoothSig(v, w)
%SMOOTHSIG  Moving mean for display, or a pass-through when w <= 1.
if w > 1, y = smoothdata(v, 'movmean', w, 'omitnan'); else, y = v; end
end


function endLabel(t, v, txt, c)
%ENDLABEL  Droplet number at the end of its curve.
%   With 26 same-coloured curves a legend is useless; a label at the tip of
%   each one is the only way to identify a droplet.
k = find(isfinite(v), 1, 'last');
if isempty(k), return; end
tf = t(isfinite(t));
if isempty(tf), dx = 0; else, dx = 0.012 * (max(tf) - min(tf)); end
text(t(k) + dx, v(k), txt, 'FontSize', 11, 'FontWeight', 'bold', ...
     'Color', 'k', 'BackgroundColor', [1 1 1 0.75], 'Margin', 0.8, ...
     'Interpreter','latex', 'VerticalAlignment','middle');
end


function noteLabels(n)
%NOTELABELS  Small box explaining what the end labels are.
annotation('textbox', [0.135 0.80 0.20 0.08], ...
    'String', sprintf('Numbers = droplet ID (S$n$), %d droplets', n), ...
    'Interpreter','latex', 'FontSize', 11, ...
    'EdgeColor', [0.75 0.75 0.75], 'BackgroundColor', 'w', ...
    'FitBoxToText', 'on', 'VerticalAlignment','middle');
end


function out = col(D, fv, ft, when)
%COL  Temporal mean of a signal, restricted to a stretch of the record.
%   'all' = everything, 'bef' = before nucleation (t < 0),
%   'aft' = after (t >= 0). The before/after split is what the whole
%   comparison rests on, so it is done once here and reused for every column.
out = nan(numel(D),1);
for i = 1:numel(D)
    v = D(i).(fv); t = D(i).(ft);
    if isempty(v) || isempty(t), continue; end
    n = min(numel(v), numel(t)); v = v(1:n); t = t(1:n);
    switch when
        case 'bef', v = v(t <  0);
        case 'aft', v = v(t >= 0);
    end
    v = v(isfinite(v));
    if ~isempty(v), out(i) = mean(v); end
end
end


function G = loadOne(dDir, useMasked)
%LOADONE  Global PIV and TrackMate series of one droplet, in absolute time.
[~, name] = fileparts(dDir);
G = struct('name',name, 'num',regexprep(name,'\D',''), 'fps',NaN, 't0',0, ...
           'nucAll',[], 'tPIV',[], 'vPIV',[], 'tTM',[], 'vTM',[]);
resF = fullfile(dDir, 'RESULTS_crystals');
evalc('specs = readSpecs(dDir);');
G.fps = specs.fps;
if isfield(specs,'t_start') && ~isnan(specs.t_start), G.t0 = specs.t_start; end

% --- every nucleation, not just the first ---
nucAll = [];
if isfield(specs,'nucleation_frames') && ~isempty(specs.nucleation_frames) ...
        && any(~isnan(specs.nucleation_frames))
    nucAll = specs.nucleation_frames(~isnan(specs.nucleation_frames));
else
    fca = dir(fullfile(resF,'crystal_areas_*.mat'));
    if ~isempty(fca)
        S = load(fullfile(resF, fca(1).name));
        if isfield(S,'nucList'), nucAll = S.nucList; end
    end
    if isempty(nucAll) && ~isnan(specs.nucleation_frame)
        nucAll = specs.nucleation_frame;
    end
end
G.nucAll = sort(nucAll(:))';

% --- PIV ---
pf = dir(fullfile(resF,'piv_fields_*.mat'));
if ~isempty(pf)
    P = load(fullfile(resF, pf(1).name), 'velocity_magnitude_from_u_v_avg');
    Vc = P.velocity_magnitude_from_u_v_avg; nP = numel(Vc);
    vm = nan(nP,1);
    for k = 1:nP
        V = Vc{k}(:); V = V(isfinite(V));
        if ~isempty(V), vm(k) = mean(V); end
    end
    G.vPIV = vm;
    G.tPIV = G.t0 + ((1:nP)' - 1)/G.fps;
end

% --- TrackMate ---
tmAll = dir(fullfile(resF,'trackmate_meanspeed_*.mat'));
if ~isempty(tmAll)
    isMasked = contains({tmAll.name},'masked');
    if useMasked && any(isMasked), idx = find(isMasked,1);
    elseif any(~isMasked),         idx = find(~isMasked,1);
    else,                          idx = 1;
    end
    S = load(fullfile(resF, tmAll(idx).name), 'framesU','meanSpeed');
    G.vTM = S.meanSpeed(:);
    G.tTM = G.t0 + double(S.framesU(:))/G.fps;
end
end


function scatterPIVTM(x, y, D, c, tfs, xlab, ylab, ttl, outPath)
%SCATTERPIVTM  One point per droplet, PIV against TrackMate, with a 1:1 line.
%   Square axes with equal limits, so the eye can judge the departure from
%   1:1 directly. The droplet number is written inside its own marker.
ok = isfinite(x) & isfinite(y);
if ~any(ok), warning('Scatter with no data: %s', ttl); return; end
fig = figure('Color','w','Position',[100 100 780 720]); hold on
lim = [0 1.08*max([x(ok); y(ok)])];
plot(lim, lim, 'k--', 'LineWidth', 1.2);
scatter(x(ok), y(ok), 260, 'filled', ...
        'MarkerFaceColor', c, 'MarkerEdgeColor','w', 'LineWidth', 1);
idx = find(ok);
for j = 1:numel(idx)
    text(x(idx(j)), y(idx(j)), D(idx(j)).num, 'Color','w', ...
         'FontSize', 9, 'FontWeight','bold', 'Interpreter','latex', ...
         'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
axis square; xlim(lim); ylim(lim); grid on; box on
set(gca, 'TickLabelInterpreter','latex', 'FontSize',12);
xlabel(xlab, 'Interpreter','latex');
ylabel(ylab, 'Interpreter','latex');
title(ttl, 'Interpreter','latex', 'FontSize', tfs);
text(0.04*lim(2), 0.96*lim(2), 'Dashed line: 1:1 agreement', ...
     'Interpreter','latex', 'FontSize', 11, 'VerticalAlignment','top');
exportgraphics(fig, outPath, 'Resolution', 300);
end


function G = loadCrys(dDir)
%LOADCRYS  Per-crystal buffer velocities of one droplet.
%   Returns empty fields if compare_individual_*.mat does not exist, so the
%   caller can warn and continue rather than stop the whole run.
[~, name] = fileparts(dDir);
G = struct('name',name, 'num',regexprep(name,'\D',''), 'fps',NaN, ...
           'tPIV',[], 'tTM',[], 'vPIV',[], 'vTM',[], 'tNuc',[]);
resF = fullfile(dDir, 'RESULTS_crystals');
cf = dir(fullfile(resF,'compare_individual_*.mat'));
if isempty(cf), return; end
S = load(fullfile(resF, cf(1).name), 'vFlowPIV_ind','vFlowTM_ind','tGrowth','tPIV');

evalc('specs = readSpecs(dDir);');
G.fps = specs.fps;
t0 = 0;
if isfield(specs,'t_start') && ~isnan(specs.t_start), t0 = specs.t_start; end

% Nucleation of each crystal. The crystals are already ordered by
% nucleation upstream, so sorting the frames keeps column c matched to
% crystal c.
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
end
G.tNuc = t0 + (nucAll(:)' - 1)/G.fps;

if isfield(S,'vFlowPIV_ind') && ~isempty(S.vFlowPIV_ind)
    G.vPIV = S.vFlowPIV_ind;
    G.tPIV = S.tPIV(:);
end
if isfield(S,'vFlowTM_ind') && ~isempty(S.vFlowTM_ind)
    G.vTM = S.vFlowTM_ind;
    G.tTM = S.tGrowth(1:size(S.vFlowTM_ind,1));
    G.tTM = G.tTM(:);
end
end


function crysPanels(C, which, col, perFig, outDir, ttl, tag)
%CRYSPANELS  Multi-panel figures: one panel per droplet, all its crystals.
%   One method per figure so the panels stay readable; the crystal index is
%   written at the end of each curve instead of using a legend.
nD = numel(C);
nFig = ceil(nD/perFig);
for f = 1:nFig
    idx = (f-1)*perFig + 1 : min(f*perFig, nD);
    fig = figure('Color','w','Position',[80 80 1150 760]);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');
    for j = idx
        G = C(j).d;
        ax = nexttile(tl); hold(ax,'on')
        if strcmp(which,'PIV'), V = G.vPIV; t = G.tPIV;
        else,                   V = G.vTM;  t = G.tTM;  end
        for c = 1:size(V,2)
            if c > numel(G.tNuc), break; end
            tr = t - G.tNuc(c);
            plot(ax, tr, V(:,c), '-', 'LineWidth', 1.3, 'Color', col);
            k = find(isfinite(V(:,c)), 1, 'last');
            if ~isempty(k)
                text(ax, tr(k), V(k,c), sprintf('%d',c), 'Color','k', ...
                     'FontSize', 10, 'FontWeight','bold', 'Interpreter','latex', ...
                     'BackgroundColor',[1 1 1 0.75], 'Margin', 0.5, ...
                     'VerticalAlignment','middle');
            end
        end
        xline(ax, 0, '-', 'Color',[0.55 0.55 0.55], 'LineWidth', 1.2);
        box(ax,'on'); grid(ax,'on')
        set(ax, 'TickLabelInterpreter','latex', 'FontSize',11);
        title(ax, sprintf('Droplet %s', G.name), 'Interpreter','latex', 'FontSize',16);
    end
    xlabel(tl, 'Time from crystal nucleation (s)', 'Interpreter','latex', 'FontSize',14);
    ylabel(tl, 'Mean velocity in buffer ($\mu$m s$^{-1}$)', ...
           'Interpreter','latex', 'FontSize',14);
    title(tl, ttl, 'Interpreter','latex', 'FontSize', 20);
    annotation(fig, 'textbox', [0.012 0.012 0.24 0.035], ...
        'String', 'Numbers = crystal index', 'Interpreter','latex', ...
        'FontSize', 11, 'EdgeColor', [0.75 0.75 0.75], ...
        'BackgroundColor','w', 'FitBoxToText','on', 'VerticalAlignment','middle');
    exportgraphics(fig, fullfile(outDir, sprintf('%s_set%d.png', tag, f)), ...
                   'Resolution', 300);
end
end
