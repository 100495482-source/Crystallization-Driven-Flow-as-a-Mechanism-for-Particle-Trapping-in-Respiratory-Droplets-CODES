%% ========================================================================
%  ACCUM1_PREPARETRACKS  Clean the raw TrackMate exports of every droplet.
%  ========================================================================
%
%  STAGE 1 OF 7 of the particle accumulation pipeline. Run the stages in
%  order; each consumes what the previous one wrote.
%
%  WHAT THIS DOES
%    For each droplet:
%      1. read the raw TrackMate export (*allspots*.csv)
%      2. discard the spots falling inside the stopwatch box
%      3. rebuild the tracks from the surviving spots
%      4. discard tracks with fewer than MIN_SPOTS points
%      5. save the result to the shared working folder
%
%  WHY THE STOPWATCH IS REMOVED FIRST
%    The recordings carry a burnt-in timer whose digits change every frame,
%    so TrackMate detects them as spots and links them into tracks with large
%    apparent velocities. The box was drawn once per droplet during the flow
%    analysis and stored as timer_box_*.mat; it is reused here. If it is
%    missing the script warns and continues unfiltered, which is visible in
%    the summary table so it cannot pass unnoticed.
%
%  WHY TRACKS ARE REBUILT AFTER FILTERING
%    Removing spots can split or empty a track, so the grouping has to be
%    redone rather than kept from the CSV.
%
%  OUTPUT  (in cfg.workDir)
%    tmclean_<droplet>.mat      Tv   table [TrackID Frame X_px Y_px]
%                               info struct with the count at each step
%    preprocessing_summary.csv  one row per droplet
%
%  See also ACCUMULATIONCONFIG, ACCUM2_REMOVECRYSTALTEXTURE.

clear; close all; clc;

cfg = accumulationConfig;
droplets = cfg.droplets;

% ---------- options ----------
MIN_SPOTS = 5;      % minimum number of spots per track

outDir = cfg.workDir;
if ~isfolder(outDir), mkdir(outDir); end
fprintf('Output: %s\n\n', outDir);

nm = {}; s0 = []; s1 = []; t0 = []; t1 = []; tb = [];

for d = 1:numel(droplets)
    dDir = droplets{d};
    [~, dName] = fileparts(dDir);
    resF = fullfile(dDir, 'RESULTS_crystals');
    fprintf('=== %s ===\n', dName);
    tG = tic;

    % ---------- 1. raw CSV ----------
    csvHit = dir(fullfile(dDir, '*allspots*.csv'));
    if isempty(csvHit)
        warning('No *allspots*.csv found in %s', dDir); continue
    end
    T = readtable(fullfile(dDir, csvHit(1).name));
    T(1,:) = [];      % TrackMate writes a units row under the header
    T = T(:,2:9);
    fr = double(T.FRAME);       x  = double(T.POSITION_X);
    y  = double(T.POSITION_Y);  id = double(T.TRACK_ID);
    % Spots detected but never linked into a track come out with NaN IDs.
    ok = ~isnan(fr) & ~isnan(x) & ~isnan(y) & ~isnan(id);
    fr = fr(ok); x = x(ok); y = y(ok); id = id(ok);
    nSpots0  = numel(fr);
    nTracks0 = numel(unique(id));
    fprintf('  CSV: %d spots, %d tracks\n', nSpots0, nTracks0);

    % ---------- 2. stopwatch ----------
    hasBox = false;
    tbf = dir(fullfile(resF, 'timer_box_*.mat'));
    if isempty(tbf)
        fprintf(2,'  WARNING: no timer_box_*.mat found, not filtering\n');
    else
        S = load(fullfile(resF, tbf(1).name), 'timerBox');
        if ~isfield(S,'timerBox') || numel(S.timerBox) ~= 4
            fprintf(2,'  WARNING: timer_box_*.mat has no valid timerBox field\n');
        else
            b = S.timerBox;
            inBox = x >= b(1) & x <= b(3) & y >= b(2) & y <= b(4);
            fprintf('  stopwatch: %d spots removed\n', sum(inBox));
            fr = fr(~inBox); x = x(~inBox); y = y(~inBox); id = id(~inBox);
            hasBox = true;
        end
    end

    % ---------- 3 and 4. regroup into tracks, enforce minimum length ------
    [uid, ~, g] = unique(id);
    cnt  = accumarray(g, 1);
    keep = ismember(g, find(cnt >= MIN_SPOTS));
    fr = fr(keep); x = x(keep); y = y(keep); id = id(keep);
    nTracks1 = numel(unique(id));
    fprintf('  tracks with >= %d spots: %d (of %d)\n', ...
            MIN_SPOTS, nTracks1, numel(uid));

    % ---------- 5. save, sorted by track then frame ----------
    [~, o] = sortrows([id fr]);
    Tv = table(id(o), fr(o), x(o), y(o), ...
               'VariableNames', {'TrackID','Frame','X_px','Y_px'});
    info = struct('droplet', dName, 'spots_raw', nSpots0, ...
                  'spots_kept', height(Tv), 'tracks_raw', nTracks0, ...
                  'tracks_kept', nTracks1, 'min_spots', MIN_SPOTS, ...
                  'timer_box_found', hasBox);
    save(fullfile(outDir, ['tmclean_' dName '.mat']), 'Tv', 'info');

    nm{end+1,1} = dName;                                    %#ok<SAGROW>
    s0(end+1,1) = nSpots0;      s1(end+1,1) = height(Tv);   %#ok<SAGROW>
    t0(end+1,1) = nTracks0;     t1(end+1,1) = nTracks1;     %#ok<SAGROW>
    tb(end+1,1) = hasBox;                                   %#ok<SAGROW>

    fprintf('  saved (%.1f s)\n\n', toc(tG));
end

% ---------- summary ----------
if isempty(nm), error('No droplet was processed.'); end
T = table(nm, s0, s1, t0, t1, logical(tb), ...
    'VariableNames', {'Droplet','SpotsRaw','SpotsKept', ...
                      'TracksRaw','TracksKept','TimerBoxFound'});
disp(T);
writetable(T, fullfile(outDir, 'preprocessing_summary.csv'));
fprintf('Done. %d droplets in %s\n', numel(nm), outDir);
