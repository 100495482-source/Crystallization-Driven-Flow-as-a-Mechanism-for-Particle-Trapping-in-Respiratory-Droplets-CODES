%% ========================================================================
%  ACCUM2_REMOVECRYSTALTEXTURE  Discard detections that are crystal texture.
%  ========================================================================
%
%  STAGE 2 OF 7 of the particle accumulation pipeline.
%
%  WHAT THIS DOES
%    For each droplet already processed by stage 1:
%      1. open the individual crystal masks
%      2. remove the tracks that are ALWAYS inside some mask
%      3. report how many tracks are PARTIAL, with spots both inside and
%         outside a mask; those are kept, only counted
%
%  THE DISTINCTION THAT MATTERS HERE
%    The crystal surface has texture, and TrackMate detects it as spots. Those
%    are not particles, and leaving them in would make the crystal's own
%    growth register as particle motion, which is exactly the effect this
%    analysis is trying to measure.
%
%    But a track cannot simply be dropped for touching a crystal: a particle
%    that swims up to a crystal and stops there is the phenomenon under study.
%    The test is therefore whether the track is EVER seen outside a mask. One
%    spot outside in any frame is enough to keep the whole track. Only tracks
%    that never leave a crystal are treated as texture.
%
%    Partial tracks are counted separately because their number is a check on
%    the segmentation: an implausibly large count would suggest the masks are
%    too generous.
%
%  IMPLEMENTATION NOTE
%    The masks are only queried at the pixels where spots actually are, frame
%    by frame, rather than being scanned in full. With thousands of tracks
%    across hundreds of frames the difference is between minutes and hours.
%
%  OUTPUT  (in cfg.workDir)
%    tmclean2_<droplet>.mat     Tv2, info2
%    crystalfilter_summary.csv
%
%  See also ACCUM1_PREPARETRACKS, ACCUM3_BUILDBUFFERZONES.

clear; close all; clc;

cfg = accumulationConfig;
droplets = cfg.droplets;

FLIP_Y = false;   % true if the TrackMate y axis is inverted

inDir = cfg.workDir;
assert(isfolder(inDir), '%s not found', inDir);
fprintf('Data in: %s\n\n', inDir);

nm = {}; nIn = []; nOut = []; nPart = []; nRem = [];

for d = 1:numel(droplets)
    dDir = droplets{d};
    [~, dName] = fileparts(dDir);
    resF = fullfile(dDir, 'RESULTS_crystals');
    fprintf('=== %s ===\n', dName);
    tG = tic;

    % ---------- cleaned tracks from stage 1 ----------
    f1 = fullfile(inDir, ['tmclean_' dName '.mat']);
    if ~isfile(f1), warning('%s is missing', f1); continue; end
    S = load(f1, 'Tv');  Tv = S.Tv;

    % ---------- individual crystal masks ----------
    kf = dir(fullfile(resF, 'individual_crystals_*.mat'));
    if isempty(kf), warning('individual_crystals_*.mat is missing'); continue; end
    K = load(fullfile(resF, kf(1).name), 'crystalMasks');
    cmask = K.crystalMasks;
    [nMF, NC] = size(cmask);
    iF = find(cellfun(@(m) ~isempty(m), cmask(:,1)), 1);
    [Ny, Nx] = size(cmask{iF,1});

    % ---------- index the spots by frame ----------
    % Sorting by frame and recording where each frame's block starts and ends
    % means each mask is loaded once and queried for all its spots together.
    [uid, ~, tidx] = unique(Tv.TrackID);
    nT = numel(uid);
    fr = Tv.Frame;  xx = Tv.X_px;  yy = Tv.Y_px;
    [fr, o] = sort(fr);  xx = xx(o); yy = yy(o); tidx = tidx(o);
    dfr = [true; diff(fr) ~= 0];
    fs = find(dfr);  fe = [fs(2:end)-1; numel(fr)];
    frList = fr(fs);

    % ---------- inside / outside, queried only at the spot pixels --------
    hasIn  = false(nT,1);
    hasOut = false(nT,1);

    for q = 1:numel(frList)
        k = frList(q) + 1;                 % TrackMate frame 0 -> index 1
        idx = fs(q):fe(q);
        rr = round(yy(idx));  cc = round(xx(idx));  tt = tidx(idx);
        if FLIP_Y, rr = Ny - rr + 1; end
        ok = rr>=1 & rr<=Ny & cc>=1 & cc<=Nx;
        if ~any(ok), continue; end
        rr = rr(ok); cc = cc(ok); tt = tt(ok);
        lin = sub2ind([Ny Nx], rr, cc);

        inAny = false(numel(lin),1);
        if k >= 1 && k <= nMF
            for c = 1:NC
                M = cmask{k,c};
                if isempty(M) || ~any(M(:)), continue; end
                inAny = inAny | M(lin);
            end
        end
        hasIn(tt(inAny))   = true;
        hasOut(tt(~inAny)) = true;
    end

    % ---------- classify ----------
    alwaysIn = hasIn & ~hasOut;      % texture: removed
    partial  = hasIn &  hasOut;      % enters and leaves: kept
    neverIn  = ~hasIn;               % never touches a crystal

    keepT = ~alwaysIn;
    Tv2 = Tv(ismember(Tv.TrackID, uid(keepT)), :);

    fprintf('  total tracks       : %d\n', nT);
    fprintf('  always inside      : %d  (REMOVED)\n', sum(alwaysIn));
    fprintf('  partial (in+out)   : %d  (kept)\n', sum(partial));
    fprintf('  never in a crystal : %d\n', sum(neverIn));
    fprintf('  tracks kept        : %d\n', sum(keepT));

    info2 = struct('droplet', dName, 'tracks_in', nT, ...
                   'always_inside', sum(alwaysIn), 'partial', sum(partial), ...
                   'never_inside', sum(neverIn), 'tracks_out', sum(keepT));
    save(fullfile(inDir, ['tmclean2_' dName '.mat']), 'Tv2', 'info2');

    nm{end+1,1}    = dName;         %#ok<SAGROW>
    nIn(end+1,1)   = nT;            %#ok<SAGROW>
    nRem(end+1,1)  = sum(alwaysIn); %#ok<SAGROW>
    nPart(end+1,1) = sum(partial);  %#ok<SAGROW>
    nOut(end+1,1)  = sum(keepT);    %#ok<SAGROW>

    fprintf('  saved (%.1f s)\n\n', toc(tG));
end

% ---------- summary ----------
if isempty(nm), error('No droplet was processed.'); end
T = table(nm, nIn, nRem, nPart, nOut, ...
    'VariableNames', {'Droplet','TracksIn','AlwaysInside_removed', ...
                      'Partial_kept','TracksOut'});
disp(T);
writetable(T, fullfile(inDir, 'crystalfilter_summary.csv'));
fprintf('Done. Results in %s\n', inDir);
