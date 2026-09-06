%% ========================================================================
%  ACCUM6_REMOVETRACKSDYINGINSIDE  Drop tracks that end inside a crystal.
%  ========================================================================
%
%  STAGE 6 OF 7 of the particle accumulation pipeline.
%
%  WHAT THIS DOES
%    Starts from tmclean2_<droplet>.mat and removes the tracks whose LAST spot
%    falls inside a crystal mask.
%
%  WHY, GIVEN STAGE 5 COUNTED EXACTLY THOSE
%    The two stages ask different questions of the same data. Stage 5 counts
%    trapped particles, so a track ending on the crystal is the signal. Stage
%    7 measures how particle endpoints are distributed with DISTANCE from the
%    crystal edge, and a point inside the crystal has no meaningful distance
%    to that edge: it would sit in the zero bin and inflate it.
%
%    Those endpoints are also not reliable observations. A particle that
%    disappears over the crystal may have been buried, or may simply have
%    become undetectable against the bright textured surface. Excluding them
%    makes the radial profile a statement about the fluid outside the crystal,
%    which is what it is meant to measure.
%
%  IMPORTANT
%    The mask only, with NO epsilon collar. A track dying in the collar is
%    KEPT: that is the region the radial profile is most interested in.
%
%  OUTPUT  (in cfg.workDir)
%    tmclean3_<droplet>.mat     Tv3, info3
%    diesinside_summary.csv
%
%  See also ACCUM5_TRAPPINGPERCRYSTAL, ACCUM7_RADIALDENSITY.

clear; close all; clc;

cfg = accumulationConfig;
droplets = cfg.droplets;

FLIP_Y = false;

inDir = cfg.workDir;
assert(isfolder(inDir), '%s not found', inDir);
fprintf('Data in: %s\n\n', inDir);

nm = {}; n0 = []; nDead = []; n1 = [];

for d = 1:numel(droplets)
    dDir = droplets{d};
    [~, dName] = fileparts(dDir);
    resF = fullfile(dDir, 'RESULTS_crystals');
    fprintf('=== %s ===\n', dName);
    tG = tic;

    f2 = fullfile(inDir, ['tmclean2_' dName '.mat']);
    if ~isfile(f2), warning('%s is missing', f2); continue; end
    S = load(f2, 'Tv2');  Tv = S.Tv2;

    kf = dir(fullfile(resF, 'individual_crystals_*.mat'));
    if isempty(kf), warning('individual_crystals_*.mat is missing'); continue; end
    K = load(fullfile(resF, kf(1).name), 'crystalMasks');
    cmask = K.crystalMasks;
    [nMF, NC] = size(cmask);
    iF = find(cellfun(@(m) ~isempty(m), cmask(:,1)), 1);
    [Ny, Nx] = size(cmask{iF,1});

    % ---------- last spot of every track ----------
    [uid, ~, tidx] = unique(Tv.TrackID);
    nT = numel(uid);
    [~, o] = sortrows([tidx Tv.Frame]);
    fr = Tv.Frame(o);  xx = Tv.X_px(o);  yy = Tv.Y_px(o);  ti = tidx(o);
    ng = [true; diff(ti) ~= 0];
    lastRow = [find(ng(2:end)); numel(ti)];   % last row of each track

    endF = fr(lastRow);  endX = xx(lastRow);  endY = yy(lastRow);

    % ---------- test only those points, frame by frame ----------
    % Each crystal mask is fetched once per frame and queried for every
    % endpoint in that frame at the same time.
    diesInside = false(nT,1);
    diesCrys   = zeros(nT,1);
    for k = unique(endF)'
        kk = k + 1;                       % TrackMate frame 0 -> index 1
        if kk < 1 || kk > nMF, continue; end
        sel = find(endF == k);
        rr = round(endY(sel));  cc = round(endX(sel));
        if FLIP_Y, rr = Ny - rr + 1; end
        ok = rr>=1 & rr<=Ny & cc>=1 & cc<=Nx;
        sel = sel(ok);  rr = rr(ok);  cc = cc(ok);
        if isempty(sel), continue; end
        lin = sub2ind([Ny Nx], rr, cc);
        for c = 1:NC
            M = cmask{kk,c};
            if isempty(M) || ~any(M(:)), continue; end
            hit = M(lin);
            diesInside(sel(hit)) = true;
            diesCrys(sel(hit))   = c;
        end
    end

    keepT = ~diesInside;
    Tv3 = Tv(ismember(Tv.TrackID, uid(keepT)), :);

    fprintf('  input tracks        : %d\n', nT);
    fprintf('  die inside a crystal: %d  (REMOVED)\n', sum(diesInside));
    for c = 1:NC
        nc = sum(diesCrys == c);
        if nc > 0, fprintf('     crystal %d: %d\n', c, nc); end
    end
    fprintf('  tracks kept         : %d\n', sum(keepT));

    info3 = struct('droplet', dName, 'tracks_in', nT, ...
                   'dies_inside_mask', sum(diesInside), ...
                   'tracks_out', sum(keepT), 'eps_used', 0);
    save(fullfile(inDir, ['tmclean3_' dName '.mat']), 'Tv3', 'info3');

    nm{end+1,1}    = dName;            %#ok<SAGROW>
    n0(end+1,1)    = nT;               %#ok<SAGROW>
    nDead(end+1,1) = sum(diesInside);  %#ok<SAGROW>
    n1(end+1,1)    = sum(keepT);       %#ok<SAGROW>
    fprintf('  saved (%.1f s)\n\n', toc(tG));
end

% ---------- summary ----------
if isempty(nm), error('No droplet was processed.'); end
T = table(nm, n0, nDead, n1, 100*nDead./max(n0,1), ...
    'VariableNames', {'Droplet','TracksIn','DieInsideMask','TracksOut', ...
                      'PercentRemoved'});
disp(T);
writetable(T, fullfile(inDir,'diesinside_summary.csv'));
fprintf('Done. Results in %s\n', inDir);
