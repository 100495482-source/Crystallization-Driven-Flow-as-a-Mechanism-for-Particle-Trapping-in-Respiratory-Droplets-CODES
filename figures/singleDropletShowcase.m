%% ========================================================================
%  SINGLEDROPLETSHOWCASE  Full flow summary of one droplet.
%  ========================================================================
%
%  PRODUCES: FIGURES ONLY. Nothing is recomputed that already exists on
%  disk, and no data file is written.
%
%  WHAT IT DRAWS
%    Four square figures for a single droplet, each combining a strip of
%    representative frames with the corresponding time series below:
%      figures 1-3 (PIV, blues)        velocity, vorticity, divergence
%      figure 4    (TrackMate, oranges) tracks coloured by speed
%
%    Putting the images and the curves in the same figure is the point. A
%    time series alone says the flow accelerated; the frames above it show
%    WHERE, and that the acceleration is localised around the crystals rather
%    than spread over the droplet.
%
%    A final section produces the same frame strip alone, without the curves,
%    for every droplet in two separate lists. Those rows are the compact
%    version used when many droplets have to fit on one page.
%
%  WHERE THE DATA COME FROM (nothing recomputed that already exists)
%    PIV  : OneDPlots_<quantity>_<droplet>.mat -> Q_mean, Q_max, Q_std
%    TM   : trackmate_meanspeed_[masked_]<droplet>.mat -> meanSpeed, maxSpeed
%           The standard deviation is NOT stored, so it is computed here and
%           smoothed with the same window, to stay comparable with the other
%           two curves.
%    Time : always from the specifications, never from the .mat files.
%
%  WHY THE TRACKS COME FROM THE RAW CSV
%    The track panels read *allspots*.csv directly and draw EVERY spot,
%    without the length and quality filters used for the measurements. The
%    filters exist to keep unreliable tracks out of the statistics; for a
%    visual impression of the flow field, showing everything is more honest
%    and denser. The two spatial cleanings are still applied, since the
%    stopwatch digits and the crystal texture are not particles at all.
%
%  ONLY EDIT SECTIONS 1, 2 AND 7.
%
%  OUTPUT
%    <FG>/dropletsummary/<droplet>/SUMMARY_PIV_<quantity>_<droplet>.png
%    <FG>/dropletsummary/<droplet>/SUMMARY_TRACKMATE_<droplet>.png
%    <FG>/dropletsummary/dropletrows/ROW_PIV_velocity_<droplet>.png
%    <FG>/dropletsummary/dropletrows/ROW_TRACKMATE_<droplet>.png
%
%  See also PIVFLOWANALYSIS, TRACKFLOWANALYSIS, GETCMAP.

clear; close all; clc;

%% ---------- 1. DROPLET ----------
DATA_ROOT  = fullfile('C:', 'Data', 'FG');
dropletDir = fullfile(DATA_ROOT, 'TP3(1202)', 'S5');
IMG_DIR_OVERRIDE = '';   % set this if the FPic<n> autodetection fails

%% ---------- 2. OPTIONS ----------
NFRAMES_PIV      = 4;      % PIV frames in the strip. 3 makes them larger
NFRAMES_TM       = 4;      % TrackMate frames, same layout as PIV
FRAMES_MANUAL    = [];     % e.g. [50 120 200 300]; [] = evenly spaced
FIRST_FRAME_MIN  = 25;     % never start before this index: the opening
                           % frames show nothing happening yet

% --- square figures ---
FIGSIDE_PIV      = 3000;   % side in px (width = height)
FIGSIDE_TM       = 3000;
ROWS_PIV         = 21;  ROWS_PIV_IMG = 9;    % 12 rows left -> 4+4+4
ROWS_TM          = 21;  ROWS_TM_IMG  = 9;    % same split as PIV

% --- PIV arrows ---
QUIVER_STEP      = 2;      % grid subsampling (1 = every node)
QUIVER_MAX_CELLS = 1.0;    % the LONGEST arrow spans this many grid cells,
                           % so arrows never overlap into neighbouring cells
QUIVER_LINEW     = 0.6;
QUIVER_COLOR     = [0.90 0 0];

% --- colour map ---
FIELD_ALPHA      = 0.30;
CMAP_NAME        = 'turbo';
CMAP_TRIM        = 0.78;   % trim the top: no dark red, which reads badly
CLIM_LO_PRCTILE  = 5;
CLIM_HI_PRCTILE  = 90;     % lower = more contrast
CLIM_MANUAL      = [];

% --- TrackMate ---
TM_TAIL          = 60;
TM_LINEW         = 0.25;   % very thin tracks: there are thousands
TM_VELOCITY_N    = 3;      % as in the pipeline, colouring only
TM_FLIP_Y        = false;  % as flipTrackMateY in the pipeline
SMOOTH_WINDOW_TM = 10;     % same window as SMOOTH_WINDOW in the pipeline
STAT_LINEW       = 2.2;

% Blues for PIV, oranges for TrackMate, as everywhere in the project.
cMeanP = [0.13 0.33 0.72];  cMaxP = [0.45 0.65 0.90];  cStdP = [0.65 0.80 0.95];
cMeanT = [0.85 0.40 0.20];  cMaxT = [0.95 0.63 0.42];  cStdT = [0.98 0.78 0.60];

CMAP = getCmap(CMAP_NAME, CMAP_TRIM);

%% ---------- 3. PATHS AND SPECIFICATIONS ----------
[~, dropletName] = fileparts(dropletDir);
dispName = dropletDisp(dropletName);   % just the number, without the S
resF = fullfile(dropletDir, 'RESULTS_crystals');
assert(isfolder(resF), '%s does not exist', resF);

fgDir  = fileparts(fileparts(dropletDir));
outDir = fullfile(fgDir, 'dropletsummary', dropletName);
if ~isfolder(outDir), mkdir(outDir); end

specs = readSpecs(dropletDir);
fps = specs.fps;
t0  = 0;
if isfield(specs,'t_start') && ~isnan(specs.t_start), t0 = specs.t_start; end

% --- the FPic<n> folder: starts with FPic and does NOT contain 'track' ---
% The trackmate folder holds the same frames with the detections burnt in,
% which is not what should sit under a figure.
dd = dir(dropletDir); dd = dd([dd.isdir]);
cand = {};
for i = 1:numel(dd)
    if ~isempty(regexpi(dd(i).name,'^FPIC','once')) && ...
       isempty(regexpi(dd(i).name,'track','once'))
        cand{end+1} = dd(i).name; %#ok<SAGROW>
    end
end
fprintf('Candidate FPic folders: %s\n', strjoin(cand, ', '));
if ~isempty(IMG_DIR_OVERRIDE)
    imgDir = IMG_DIR_OVERRIDE;
else
    assert(~isempty(cand), 'No FPic<n> folder found in %s', dropletDir);
    imgDir = fullfile(dropletDir, cand{1});
end
fprintf('Images used: %s\n', imgDir);

imgF = dir(fullfile(imgDir,'*.tif'));
if isempty(imgF), imgF = dir(fullfile(imgDir,'*.png')); end
assert(~isempty(imgF), 'No images in %s', imgDir);
imgF = sortNumeric(imgF);
fprintf('  %d images found.\n', numel(imgF));

%% ---------- 4. LOAD THE PIV FIELDS ----------
pf = dir(fullfile(resF,'piv_fields_*.mat'));
assert(~isempty(pf), 'piv_fields_*.mat is missing in %s', resF);
P = load(fullfile(resF, pf(1).name), ...
    'velocity_magnitude_from_u_v_avg','vorticity_from_u_v_avg', ...
    'divergence_from_u_v_avg','u_avg_um_s','v_avg_um_s','x_um','y_um','um_per_px');
um_per_px = P.um_per_px;

nPIV = numel(P.velocity_magnitude_from_u_v_avg);
nUse = min(nPIV, numel(imgF));
frameTimes = t0 + (0:nUse-1)/fps;

% Frames without a temporal average (the first and last K_AVG) are excluded
% from the candidates, so the strip never shows an empty panel.
okF = find(~cellfun(@isempty, P.velocity_magnitude_from_u_v_avg(1:nUse)));
if isempty(okF), okF = 1:nUse; end

framesPIV = pickFrames(FRAMES_MANUAL, okF, NFRAMES_PIV, nUse, FIRST_FRAME_MIN);
framesTM  = pickFrames(FRAMES_MANUAL, okF, NFRAMES_TM,  nUse, FIRST_FRAME_MIN);
fprintf('PIV frames: %s\nTM frames : %s\n', mat2str(framesPIV), mat2str(framesTM));

%% ---------- 5. PIV FIGURES ----------
% isSym marks the fields that are signed and must use a colour scale centred
% on zero; velocity magnitude is not one of them.
fieldSpecs = { ...
 'velocity',   P.velocity_magnitude_from_u_v_avg, 'Velocity magnitude', '$|v|$ ($\mu$m/s)', false
 'vorticity',  P.vorticity_from_u_v_avg,          'Vorticity',          '$\omega$ (1/s)',   true
 'divergence', P.divergence_from_u_v_avg,         'Divergence',         '$\nabla\cdot v$ (1/s)', true};

for q = 1:size(fieldSpecs,1)
    kind = fieldSpecs{q,1};  Fcell = fieldSpecs{q,2};
    ttl  = fieldSpecs{q,3};  ylab  = fieldSpecs{q,4};  isSym = fieldSpecs{q,5};

    if ~isempty(CLIM_MANUAL)
        climVal = CLIM_MANUAL;
    else
        climVal = climFromCells(Fcell, isSym, CLIM_LO_PRCTILE, CLIM_HI_PRCTILE);
    end
    fprintf('  %s: colour scale [%.3g  %.3g]\n', kind, climVal(1), climVal(2));

    % --- 1D series: prefer the stored ones, which are already smoothed ---
    % Recomputing here would give unsmoothed curves that do not match the
    % ones reported elsewhere in the document.
    mf = fullfile(resF, sprintf('OneDPlots_%s_%s.mat', kind, dropletName));
    if isfile(mf)
        S1 = load(mf, 'Q_mean','Q_max','Q_std');
        Qm = S1.Q_mean(:); Qx = S1.Q_max(:); Qs = S1.Q_std(:);
        tStats = t0 + (0:numel(Qm)-1)'/fps;      % time from the specifications
        fprintf('  %s: 1D series loaded from %s\n', kind, mf);
    else
        [Qm, Qx, Qs] = fieldStats(Fcell, nUse);
        tStats = frameTimes(:);
        warning('%s: .mat not found, recomputing from the fields.', kind);
    end

    fig = figure('Color','w','Position',[20 20 FIGSIDE_PIV FIGSIDE_PIV]);
    tl = tiledlayout(fig, ROWS_PIV, 1, 'TileSpacing','compact', 'Padding','compact');

    % --- frame strip, panels touching each other ---
    tlTop = tiledlayout(tl, 1, numel(framesPIV), ...
                        'TileSpacing','none', 'Padding','none');
    tlTop.Layout.Tile = 1;  tlTop.Layout.TileSpan = [ROWS_PIV_IMG 1];
    axLast = [];
    for s = 1:numel(framesPIV)
        k  = framesPIV(s);
        ax = nexttile(tlTop); axLast = ax;
        drawFieldPanel(ax, fullfile(imgDir, imgF(k).name), um_per_px, ...
            P.x_um{k}, P.y_um{k}, Fcell{k}, P.u_avg_um_s{k}, P.v_avg_um_s{k}, ...
            climVal, QUIVER_STEP, QUIVER_MAX_CELLS, ...
            QUIVER_LINEW, QUIVER_COLOR, FIELD_ALPHA, CMAP);
        stampTime(ax, sprintf('$t = %.1f$ s', frameTimes(k)), 16);
    end
    if ~isempty(axLast)
        cb = colorbar(axLast); cb.Layout.Tile = 'east';
        cb.Label.String = ttl; cb.TickLabelInterpreter = 'latex';
        cb.Label.FontSize = 15;
    end

    % --- the three curves, all given the same height ---
    % Equal heights matter: a taller panel makes its curve look more
    % significant than the others.
    nRest = ROWS_PIV - ROWS_PIV_IMG;  blk = floor(nRest/3);
    r1 = ROWS_PIV_IMG + 1;  r2 = r1 + blk;  r3 = r2 + blk;
    statRow(nexttile(tl,r1,[blk 1]), tStats, Qx, cMaxP,  ['Max ' ylab],  STAT_LINEW);
    statRow(nexttile(tl,r2,[blk 1]), tStats, Qm, cMeanP, ['Mean ' ylab], STAT_LINEW);
    ax3 = nexttile(tl,r3,[blk 1]);
    statRow(ax3, tStats, Qs, cStdP, ['Std ' ylab], STAT_LINEW);
    xlabel(ax3, 'Time (s)', 'Interpreter','latex', 'FontSize', 16);

    title(tl, sprintf('Droplet %s: %s (PIV)', dispName, ttl), ...
          'Interpreter','latex', 'FontSize', 24);

    exportgraphics(fig, fullfile(outDir, ...
        sprintf('SUMMARY_PIV_%s_%s.png', kind, dropletName)), 'Resolution',300);
    fprintf('Saved: SUMMARY_PIV_%s_%s.png\n', kind, dropletName);
end

%% ---------- 6. TRACKMATE FIGURE ----------
msAll = dir(fullfile(resF,'trackmate_meanspeed_*.mat'));
msPath = '';
if ~isempty(msAll)
    isM = contains({msAll.name},'masked');
    if any(isM), msPath = fullfile(resF, msAll(find(isM,1)).name);
    else,        msPath = fullfile(resF, msAll(1).name); end
end

fMv = dir(fullfile(resF,'trackmate_velocities_masked_*.mat'));
fUv = dir(fullfile(resF,'trackmate_velocities_*.mat'));
fUv = fUv(~contains({fUv.name},'masked'));

if isempty(msPath)
    warning('trackmate_meanspeed_*.mat is missing: skipping that figure.');
else
    Sm = load(msPath, 'framesU','meanSpeed','maxSpeed');

    if contains(msPath,'masked') && ~isempty(fMv)
        Sst = load(fullfile(resF, fMv(1).name), 'trackVelocityTable_reduced');
    elseif ~isempty(fUv)
        Sst = load(fullfile(resF, fUv(1).name), 'trackVelocityTable_reduced');
    else
        Sst = [];
    end

    tTM = t0 + double(Sm.framesU(:))/fps;   % time from the specifications
    Tm  = Sm.meanSpeed(:);                  % already computed and smoothed
    Tx  = Sm.maxSpeed(:);                   % already computed and smoothed

    % --- the standard deviation is not stored, so it is computed here ---
    % Smoothed with the same window as the mean and max, otherwise the three
    % panels would not be comparable.
    Ts = nan(numel(Sm.framesU),1);
    if ~isempty(Sst)
        Tstd = Sst.trackVelocityTable_reduced;
        for i = 1:numel(Sm.framesU)
            s = Tstd.Speed_um_s(Tstd.Frame == Sm.framesU(i));
            s = s(isfinite(s));
            if ~isempty(s), Ts(i) = std(s); end
        end
        validS = ~isnan(Ts);
        Ts = smoothdata(Ts, 'movmean', SMOOTH_WINDOW_TM, 'omitnan');
        Ts(~validS) = NaN;
    else
        warning('No trackmate_velocities_*.mat: cannot compute the std.');
    end

    % --- the table used for DRAWING: raw CSV, no length or quality filter,
    %     but with the stopwatch and crystal-interior spots removed ---
    Tv = loadAllSpotsVel(dropletDir, fps, um_per_px, TM_VELOCITY_N, TM_FLIP_Y);
    fprintf('Tracks drawn: %d, spots: %d\n', ...
            numel(unique(Tv.TrackID)), height(Tv));

    % Sanity check: coordinates outside the image mean the TrackMate run used
    % a differently cropped sequence, and the overlay would be meaningless.
    imgTest = imread(fullfile(imgDir, imgF(1).name));
    if ndims(imgTest)==3, imgTest = im2gray(imgTest); end
    [Nyi, Nxi] = size(imgTest);
    if max(Tv.X_px) > Nxi || max(Tv.Y_px) > Nyi
        warning('TrackMate coordinates fall outside the image (%dx%d).', Nxi, Nyi);
    end

    sp = Tv.Speed_um_s(isfinite(Tv.Speed_um_s));
    if ~isempty(CLIM_MANUAL)
        climV = CLIM_MANUAL;
    else
        climV = [prctile(sp, CLIM_LO_PRCTILE), prctile(sp, CLIM_HI_PRCTILE)];
        if climV(2) <= climV(1), climV = [0, max(sp)]; end
    end
    fprintf('  TrackMate: colour scale [%.3g  %.3g]\n', climV(1), climV(2));

    fig = figure('Color','w','Position',[20 20 FIGSIDE_TM FIGSIDE_TM]);
    tl  = tiledlayout(fig, ROWS_TM, 1, 'TileSpacing','compact', 'Padding','compact');

    tlTop = tiledlayout(tl, 1, numel(framesTM), 'TileSpacing','none', 'Padding','none');
    tlTop.Layout.Tile = 1;  tlTop.Layout.TileSpan = [ROWS_TM_IMG 1];
    axLast = [];
    for s = 1:numel(framesTM)
        k  = framesTM(s);
        ax = nexttile(tlTop); axLast = ax;
        % k-1 because TrackMate counts frames from 0 while image 1 is index 1.
        drawTrackPanel(ax, fullfile(imgDir, imgF(k).name), Tv, ...
                       k-1, TM_TAIL, climV, TM_LINEW, CMAP);
        stampTime(ax, sprintf('$t = %.1f$ s', frameTimes(k)), 17);
    end
    if ~isempty(axLast)
        cb = colorbar(axLast); cb.Layout.Tile = 'east';
        cb.Label.String = 'Track speed (\mum/s)';
        cb.TickLabelInterpreter = 'latex'; cb.Label.FontSize = 15;
    end

    ylabTM = '$v$ ($\mu$m/s)';
    nRest = ROWS_TM - ROWS_TM_IMG;  blk = floor(nRest/3);
    r1 = ROWS_TM_IMG + 1;  r2 = r1 + blk;  r3 = r2 + blk;
    statRow(nexttile(tl,r1,[blk 1]), tTM, Tx, cMaxT,  ['Max ' ylabTM],  STAT_LINEW);
    statRow(nexttile(tl,r2,[blk 1]), tTM, Tm, cMeanT, ['Mean ' ylabTM], STAT_LINEW);
    ax3 = nexttile(tl,r3,[blk 1]);
    statRow(ax3, tTM, Ts, cStdT, ['Std ' ylabTM], STAT_LINEW);
    xlabel(ax3, 'Time (s)', 'Interpreter','latex', 'FontSize', 16);

    title(tl, sprintf('Droplet %s: particle speed (TrackMate)', dispName), ...
          'Interpreter','latex', 'FontSize', 24);

    exportgraphics(fig, fullfile(outDir, ...
        sprintf('SUMMARY_TRACKMATE_%s.png', dropletName)), 'Resolution',300);
    fprintf('Saved: SUMMARY_TRACKMATE_%s.png\n', dropletName);
end

fprintf('\nDone. Results in %s\n', outDir);

%% ---------- 7. FRAME ROWS, ONE PER DROPLET ----------
%  The frame strip alone, without the curves, in the same layout and style
%  as above. Two independent lists so the PIV and TrackMate rows can be
%  generated for different subsets.

sessionsRow = { ...
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

dropletsRowPIV = cell(size(sessionsRow,1),1);
for i = 1:size(sessionsRow,1)
    dropletsRowPIV{i} = fullfile(DATA_ROOT, sessionsRow{i,1}, sessionsRow{i,2});
end
dropletsRowTM = dropletsRowPIV;   % edit here to use a different subset

NFRAMES_ROW  = 4;              % frames per row
FIGSIZE_ROW  = [2600 800];     % width x height of the row

rowDir = fullfile(fgDir, 'dropletsummary', 'dropletrows');
if ~isfolder(rowDir), mkdir(rowDir); end

% ===== PIV (velocity magnitude) =====
for d = 1:numel(dropletsRowPIV)
    dD = dropletsRowPIV{d};
    [~, dName] = fileparts(dD);
    rF = fullfile(dD, 'RESULTS_crystals');
    if ~isfolder(rF), warning('%s: no RESULTS_crystals', dName); continue; end

    evalc('sp2 = readSpecs(dD);');
    fps2 = sp2.fps;  tA = 0;
    if isfield(sp2,'t_start') && ~isnan(sp2.t_start), tA = sp2.t_start; end

    iDir = findImgDir(dD);
    if isempty(iDir), warning('%s: no FPic folder', dName); continue; end
    iF = dir(fullfile(iDir,'*.tif'));
    if isempty(iF), iF = dir(fullfile(iDir,'*.png')); end
    if isempty(iF), warning('%s: no images', dName); continue; end
    iF = sortNumeric(iF);

    pF = dir(fullfile(rF,'piv_fields_*.mat'));
    if isempty(pF), warning('%s: no piv_fields', dName); continue; end
    P2 = load(fullfile(rF, pF(1).name), ...
        'velocity_magnitude_from_u_v_avg','u_avg_um_s','v_avg_um_s', ...
        'x_um','y_um','um_per_px');
    Fc = P2.velocity_magnitude_from_u_v_avg;

    nU = min(numel(Fc), numel(iF));
    ft = tA + (0:nU-1)/fps2;
    ok2 = find(~cellfun(@isempty, Fc(1:nU)));
    if isempty(ok2), ok2 = 1:nU; end
    frRow = pickFrames([], ok2, NFRAMES_ROW, nU, FIRST_FRAME_MIN);

    % The colour scale is computed per droplet, not shared: droplets differ
    % by more than an order of magnitude in flow speed and a common scale
    % would flatten most of them.
    if ~isempty(CLIM_MANUAL)
        cV = CLIM_MANUAL;
    else
        cV = climFromCells(Fc, false, CLIM_LO_PRCTILE, CLIM_HI_PRCTILE);
    end

    fig = figure('Color','w','Position',[20 20 FIGSIZE_ROW]);
    tlR = tiledlayout(fig, 1, numel(frRow), ...
                      'TileSpacing','none', 'Padding','compact');
    axL = [];
    for s = 1:numel(frRow)
        k = frRow(s);
        ax = nexttile(tlR); axL = ax;
        drawFieldPanel(ax, fullfile(iDir, iF(k).name), P2.um_per_px, ...
            P2.x_um{k}, P2.y_um{k}, Fc{k}, P2.u_avg_um_s{k}, P2.v_avg_um_s{k}, ...
            cV, QUIVER_STEP, QUIVER_MAX_CELLS, ...
            QUIVER_LINEW, QUIVER_COLOR, FIELD_ALPHA, CMAP);
        stampTime(ax, sprintf('$t = %.1f$ s', ft(k)), 16);
    end
    if ~isempty(axL)
        cb = colorbar(axL); cb.Layout.Tile = 'east';
        cb.Label.String = 'Velocity magnitude';
        cb.TickLabelInterpreter = 'latex'; cb.Label.FontSize = 15;
    end
    title(tlR, sprintf('Droplet %s: velocity magnitude (PIV)', dropletDisp(dName)), ...
          'Interpreter','latex', 'FontSize', 26);

    exportgraphics(fig, fullfile(rowDir, ...
        sprintf('ROW_PIV_velocity_%s.png', dName)), 'Resolution',300);
    close(fig);
    fprintf('Saved: ROW_PIV_velocity_%s.png\n', dName);
end

% ===== TrackMate =====
for d = 1:numel(dropletsRowTM)
    dD = dropletsRowTM{d};
    [~, dName] = fileparts(dD);
    rF = fullfile(dD, 'RESULTS_crystals');
    if ~isfolder(rF), warning('%s: no RESULTS_crystals', dName); continue; end

    evalc('sp2 = readSpecs(dD);');
    fps2 = sp2.fps;  tA = 0;
    if isfield(sp2,'t_start') && ~isnan(sp2.t_start), tA = sp2.t_start; end

    iDir = findImgDir(dD);
    if isempty(iDir), warning('%s: no FPic folder', dName); continue; end
    iF = dir(fullfile(iDir,'*.tif'));
    if isempty(iF), iF = dir(fullfile(iDir,'*.png')); end
    if isempty(iF), warning('%s: no images', dName); continue; end
    iF = sortNumeric(iF);

    % Reference frame count: the PIV one if it exists, otherwise the images,
    % so the TrackMate rows show the same frames as the PIV rows.
    pF = dir(fullfile(rF,'piv_fields_*.mat'));
    if ~isempty(pF)
        M = matfile(fullfile(rF, pF(1).name));
        nU = min(numel(M.velocity_magnitude_from_u_v_avg), numel(iF));
    else
        nU = numel(iF);
    end
    ft = tA + (0:nU-1)/fps2;
    frRow = pickFrames([], 1:nU, NFRAMES_ROW, nU, FIRST_FRAME_MIN);

    try
        Tv2 = loadAllSpotsVel(dD, fps2, sp2.ppf, TM_VELOCITY_N, TM_FLIP_Y);
    catch ME
        warning('%s: %s', dName, ME.message); continue
    end
    sp3 = Tv2.Speed_um_s(isfinite(Tv2.Speed_um_s));
    if isempty(sp3), warning('%s: no velocities', dName); continue; end
    if ~isempty(CLIM_MANUAL)
        cV = CLIM_MANUAL;
    else
        cV = [prctile(sp3, CLIM_LO_PRCTILE), prctile(sp3, CLIM_HI_PRCTILE)];
        if cV(2) <= cV(1), cV = [0, max(sp3)]; end
    end

    fig = figure('Color','w','Position',[20 20 FIGSIZE_ROW]);
    tlR = tiledlayout(fig, 1, numel(frRow), ...
                      'TileSpacing','none', 'Padding','compact');
    axL = [];
    for s = 1:numel(frRow)
        k = frRow(s);
        ax = nexttile(tlR); axL = ax;
        drawTrackPanel(ax, fullfile(iDir, iF(k).name), Tv2, ...
                       k-1, TM_TAIL, cV, TM_LINEW, CMAP);
        stampTime(ax, sprintf('$t = %.1f$ s', ft(k)), 16);
    end
    if ~isempty(axL)
        cb = colorbar(axL); cb.Layout.Tile = 'east';
        cb.Label.String = 'Track speed (\mum/s)';
        cb.TickLabelInterpreter = 'latex'; cb.Label.FontSize = 15;
    end
    title(tlR, sprintf('Droplet %s: particle speed (TrackMate)', dropletDisp(dName)), ...
          'Interpreter','latex', 'FontSize', 26);

    exportgraphics(fig, fullfile(rowDir, ...
        sprintf('ROW_TRACKMATE_%s.png', dName)), 'Resolution',300);
    close(fig);
    fprintf('Saved: ROW_TRACKMATE_%s.png\n', dName);
end


%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================
%  getCmap lives in utils/ and is shared with the accumulation analysis.

function stampTime(ax, str, fs)
%STAMPTIME  Time label INSIDE the image, top left.
%   Placed inside rather than as a title so the frames can sit flush against
%   each other with no gap between panels.
text(ax, 0.035, 0.96, str, 'Units','normalized', ...
     'Interpreter','latex', 'FontSize', fs, 'Color','w', ...
     'FontWeight','bold', 'BackgroundColor',[0 0 0 0.45], 'Margin', 3, ...
     'VerticalAlignment','top', 'HorizontalAlignment','left');
end


function f = sortNumeric(f)
%SORTNUMERIC  Sort a dir() struct by the LAST number in each filename.
%   The last one, because the folder name is often repeated in the filename
%   before the frame index.
nm = {f.name}; num = nan(size(nm));
for i = 1:numel(nm)
    tk = regexp(nm{i}, '\d+', 'match');
    if ~isempty(tk), num(i) = str2double(tk{end}); end
end
if all(~isnan(num)), [~,ix] = sort(num); else, [~,ix] = sort(nm); end
f = f(ix);
end


function fr = pickFrames(manual, okF, n, nUse, minStart)
%PICKFRAMES  Evenly spaced frames, never starting before minStart.
%   The opening frames of a recording show a droplet in which nothing has
%   happened yet; spending a panel on one of them wastes a quarter of the
%   strip.
if nargin < 5 || isempty(minStart), minStart = 1; end
if ~isempty(manual)
    fr = manual(:)';
else
    f1 = max(okF(1), minStart);
    fEnd = okF(end);
    if f1 >= fEnd, f1 = okF(1); end          % sequence too short for the margin
    fr = round(linspace(f1, fEnd, n));
end
fr = unique(min(max(fr,1), nUse), 'stable');
end


function d = dropletDisp(name)
%DROPLETDISP  'S5' -> '5', for figure titles.
d = regexprep(name, '\D', '');
if isempty(d), d = name; end
end


function climVal = climFromCells(Fcell, isSym, prLo, prHi)
%CLIMFROMCELLS  Colour limits from percentiles over every frame.
%   Percentiles rather than min/max, so a single spurious cell does not
%   compress the whole sequence into one colour. Signed fields get a scale
%   centred on zero, so that positive and negative read symmetrically.
allVals = [];
for k = 1:numel(Fcell)
    S = Fcell{k};
    if ~isempty(S), allVals = [allVals; S(isfinite(S))]; end %#ok<AGROW>
end
if isempty(allVals), climVal = [0 1]; return; end
if isSym
    m = prctile(abs(allVals), prHi); if m <= 0, m = 1; end
    climVal = [-m m];
else
    lo = prctile(allVals, prLo);  hi = prctile(allVals, prHi);
    if hi <= lo, lo = min(allVals); hi = max(allVals); end
    if hi <= lo, hi = lo + 1; end
    climVal = [lo hi];
end
end


function [Qm, Qx, Qs] = fieldStats(Fcell, n)
%FIELDSTATS  Fallback only: same as pivFlowAnalysis part B, unsmoothed.
%   Used when the stored 1D series are missing. The curves it produces are
%   noisier than the ones in the rest of the document.
Qm = nan(n,1); Qx = nan(n,1); Qs = nan(n,1);
for k = 1:n
    Q = Fcell{k};
    if isempty(Q), continue; end
    Q = Q(isfinite(Q));
    if isempty(Q), continue; end
    Qm(k) = mean(Q); Qx(k) = max(Q); Qs(k) = std(Q);
end
end


function statRow(ax, t, y, col, ylab, lw)
%STATROW  One time-series panel, styled to match its neighbours.
n = min(numel(t), numel(y)); t = t(1:n); y = y(1:n);
plot(ax, t, y, '-', 'LineWidth', lw, 'Color', col);
grid(ax,'on'); grid(ax,'minor'); box(ax,'on')
set(ax, 'TickLabelInterpreter','latex', 'FontSize', 14, ...
        'GridAlpha', 0.25, 'MinorGridAlpha', 0.10);
ylabel(ax, ylab, 'Interpreter','latex', 'FontSize', 15);
tf = t(isfinite(t));
if numel(tf) > 1, xlim(ax, [min(tf) max(tf)]); end
end


function drawFieldPanel(ax, imgPath, um_per_px, X, Y, S, U, V, climVal, ...
                        step, maxCells, qlw, qcol, alph, cmap)
%DRAWFIELDPANEL  One frame with a translucent scalar field and velocity arrows.
img = imread(imgPath);
if ndims(img)==3, img = im2gray(img); end
img = double(img);
[Ny, Nx] = size(img);
x_img = (0:Nx-1)*um_per_px;  y_img = (0:Ny-1)*um_per_px;

% Background converted to explicit RGB, so the colormap applied to the
% scalar field afterwards cannot repaint it.
cmin = min(img(:)); cmax = max(img(:));
if cmax > cmin, idx = round((img-cmin)/(cmax-cmin)*255)+1; else, idx = ones(size(img)); end
rgb = ind2rgb(idx, gray(256));

hold(ax,'on');
imagesc(ax, x_img, y_img, rgb);

if ~isempty(S)
    hS = imagesc(ax, X(1,:), Y(:,1), S);
    set(hS, 'AlphaData', alph * ~isnan(S));
    colormap(ax, cmap); clim(ax, climVal);
end

if ~isempty(U)
    % --- arrow scale BOUNDED to the grid cell ---
    % MATLAB's automatic quiver scaling varies between frames, which would
    % make an unchanged flow appear to speed up. Here the longest arrow is
    % pinned to maxCells grid cells, so arrow length is comparable across
    % every panel and every droplet.
    Us = U(1:step:end, 1:step:end);
    Vs = V(1:step:end, 1:step:end);
    Xs = X(1:step:end, 1:step:end);
    Ys = Y(1:step:end, 1:step:end);
    mag = hypot(Us, Vs);  mmax = max(mag(isfinite(mag)));
    dxg = mean(diff(X(1,:)), 'omitnan') * step;
    if isempty(mmax) || mmax <= 0 || ~isfinite(dxg) || dxg <= 0
        scUse = 1;
    else
        scUse = maxCells * dxg / mmax;
    end
    qv = quiver(ax, Xs, Ys, scUse*Us, scUse*Vs, 0, ...
                'Color', qcol, 'LineWidth', qlw, 'MaxHeadSize', 0.9);
    qv.Clipping = 'on';
end

set(ax,'YDir','normal'); axis(ax,'image');
xlim(ax,[min(x_img) max(x_img)]); ylim(ax,[min(y_img) max(y_img)]);
set(ax,'XTick',[],'YTick',[]); box(ax,'on');
end


function drawTrackPanel(ax, imgPath, Tv, k, tail, climV, lw, cmap)
%DRAWTRACKPANEL  One frame with the recent track tails coloured by speed.
%   Built as a single patch: thousands of individual plot calls would make
%   the figure impossible to render or export.
img = imread(imgPath);
if ndims(img)==3, img = im2gray(img); end
img = double(img);
[Ny, Nx] = size(img);
cmin = min(img(:)); cmax = max(img(:));
if cmax > cmin, idx = round((img-cmin)/(cmax-cmin)*255)+1; else, idx = ones(size(img)); end
rgb = ind2rgb(idx, gray(256));

hold(ax,'on');
imagesc(ax, rgb);
set(ax,'YDir','normal'); axis(ax,'image');
xlim(ax,[0 Nx]); ylim(ax,[0 Ny]);

sub = Tv(Tv.Frame <= k & Tv.Frame >= k-tail, :);
if ~isempty(sub)
    sub = sortrows(sub, {'TrackID','Frame'});
    n = height(sub);
    if n > 1
        i1 = (1:n-1)';
        ok = sub.TrackID(i1) == sub.TrackID(i1+1);   % no segment across tracks
        i1 = i1(ok);
        if ~isempty(i1)
            Vtx = [sub.X_px, sub.Y_px];
            cd  = sub.Speed_um_s;  cd(~isfinite(cd)) = climV(1);
            patch(ax, 'Faces', [i1 i1+1], 'Vertices', Vtx, ...
                  'FaceVertexCData', cd, 'EdgeColor','interp', ...
                  'FaceColor','none', 'LineWidth', lw);
        end
    end
end
colormap(ax, cmap); clim(ax, climV);
set(ax,'XTick',[],'YTick',[]); box(ax,'on');
end


function Tv = loadAllSpotsVel(dropletDir, fps, um_per_px, Nv, flipY)
%LOADALLSPOTSVEL  Read the raw TrackMate CSV and return spots with velocity.
%
%   Does NOT filter by track length or detection quality, because these
%   tracks are for drawing and a denser field reads better. It DOES apply the
%   same two spatial cleanings as trackFlowAnalysis, since those remove
%   things that are not particles at all:
%     1) the stopwatch region (timer_box_*.mat), only if already stored
%     2) spots inside the crystals (masks_*.mat)
%
%   Nothing is saved and no .mat is modified.

if nargin < 5 || isempty(flipY), flipY = false; end

csvHit = dir(fullfile(dropletDir, '*allspots*.csv'));
assert(~isempty(csvHit), 'No *allspots*.csv found in %s', dropletDir);
T = readtable(fullfile(dropletDir, csvHit(1).name));
T(1,:) = [];  T = T(:,2:9);          % same as the main reading stage

fr = double(T.FRAME);  xx = double(T.POSITION_X);
yy = double(T.POSITION_Y);  id = double(T.TRACK_ID);
ok = ~isnan(fr) & ~isnan(xx) & ~isnan(yy) & ~isnan(id);
fr = fr(ok); xx = xx(ok); yy = yy(ok); id = id(ok);
nRaw = numel(fr);

resF = fullfile(dropletDir, 'RESULTS_crystals');

% --- 1) stopwatch region, only if the box has already been drawn ---
% This script never opens a drawing dialog: it is meant to run unattended
% over a whole list of droplets.
tb = dir(fullfile(resF, 'timer_box_*.mat'));
if ~isempty(tb)
    Sb = load(fullfile(resF, tb(1).name), 'timerBox');
    if isfield(Sb,'timerBox') && numel(Sb.timerBox) == 4
        b = Sb.timerBox;
        keep = ~(xx >= b(1) & xx <= b(3) & yy >= b(2) & yy <= b(4));
        fprintf('  stopwatch: %d spots removed\n', sum(~keep));
        fr = fr(keep); xx = xx(keep); yy = yy(keep); id = id(keep);
    end
else
    fprintf('  no timer_box_*.mat: the stopwatch is not filtered.\n');
end

% --- 2) spots inside the crystals ---
mk = dir(fullfile(resF, 'masks_*.mat'));
if ~isempty(mk)
    Smk = load(fullfile(resF, mk(1).name), 'masks');
    masks = Smk.masks;
    keep = true(size(fr));
    for s2 = 1:numel(fr)
        f1 = fr(s2) + 1;                       % TrackMate frame 0 -> index 1
        if f1 < 1 || f1 > numel(masks), continue; end
        M = masks{f1};
        if isempty(M), continue; end
        col = round(xx(s2)); row = round(yy(s2));
        if flipY, row = size(M,1) - row + 1; end
        if row >= 1 && row <= size(M,1) && col >= 1 && col <= size(M,2)
            if M(row,col), keep(s2) = false; end
        end
    end
    fprintf('  crystals: %d spots removed\n', sum(~keep));
    fr = fr(keep); xx = xx(keep); yy = yy(keep); id = id(keep);
else
    warning('No masks_*.mat found: tracks will be drawn over the crystals.');
end
fprintf('  spots: %d of %d after cleaning\n', numel(fr), nRaw);

% --- velocities on the cleaned spots ---
% The window is clamped to each track's own ends rather than skipped, so
% short tracks still get a velocity and can be drawn.
uid = unique(id);
aID=[]; aF=[]; aX=[]; aY=[]; aSp=[];
for i = 1:numel(uid)
    m = (id == uid(i));
    [f, o] = sort(fr(m));
    x = xx(m); y = yy(m); x = x(o); y = y(o);
    n = numel(f);
    v = nan(n,1);
    for j = 1:n
        ip = max(1, j-Nv);  in = min(n, j+Nv);
        df = f(in) - f(ip);
        if df <= 0, continue; end
        dt = df / fps;
        v(j) = hypot(x(in)-x(ip), y(in)-y(ip)) * um_per_px / dt;
    end
    aID = [aID; repmat(uid(i), n, 1)]; %#ok<AGROW>
    aF  = [aF;  f];  aX = [aX; x];  aY = [aY; y];  aSp = [aSp; v]; %#ok<AGROW>
end
Tv = table(aID, aF, aX, aY, aSp, 'VariableNames', ...
    {'TrackID','Frame','X_px','Y_px','Speed_um_s'});
end


function iDir = findImgDir(dD)
%FINDIMGDIR  The FPic<n> folder of a droplet: starts with FPic, no 'track'.
iDir = '';
dd = dir(dD); dd = dd([dd.isdir]);
for i = 1:numel(dd)
    if ~isempty(regexpi(dd(i).name,'^FPIC','once')) && ...
       isempty(regexpi(dd(i).name,'track','once'))
        iDir = fullfile(dD, dd(i).name); return
    end
end
end
