%% ========================================================================
%  ACCUM5_TRAPPINGPERCRYSTAL  How many particles end up caught at a crystal.
%  ========================================================================
%
%  STAGE 5 OF 7 of the particle accumulation pipeline. This is the stage that
%  answers the question the thesis is asking.
%
%  WHAT THIS DOES
%    Uses the data prepared by the earlier stages:
%      tmclean2_<droplet>.mat  clean tracks (no stopwatch, >= MIN_SPOTS
%                              points, no crystal texture)
%      buffers_<droplet>.mat   non-overlapping buffer zones, clipped to the
%                              droplet outline
%
%    For each crystal it counts:
%      - how many tracks END inside its buffer, i.e. their last spot falls
%        within its region of influence
%      - of those, how many are TRAPPED: the last spot lies inside the
%        crystal mask, or within an EPS_PX collar around it
%
%    And for the droplet as a whole, the totals and percentages.
%
%  WHY THE LAST SPOT
%    A particle carried to a crystal and held there stops being detected as
%    moving. Its final position is therefore the observable that distinguishes
%    trapping from a particle merely passing nearby.
%
%  WHY THE COLLAR FOLLOWS THE CRYSTAL SHAPE
%    The collar is a dilation of the mask itself, not a circle, so it hugs a
%    faceted crystal instead of covering fluid far from its edges. It is
%    recomputed for each frame, so it always matches the crystal as it was
%    when that particular particle stopped, not as it ended up.
%
%    A track that both starts and ends in the collar is included, since its
%    final state is what defines trapping.
%
%  FIGURE
%    The last frame, with every track in translucent grey and the trapped
%    ones coloured by speed, with a red dot at each trapped endpoint.
%
%  OUTPUT  (in cfg.outDir)
%    trap_percrystal_<droplet>.csv, trap_map_<droplet>.png
%    trap_summary.csv, trap_percrystal_all.csv
%
%  See also ACCUM4_CLIPBUFFERSTODROPLET, BUILDBUFFERFRAME, DRAWTRACKSET.

clear; close all; clc;

cfg = accumulationConfig;
droplets = cfg.droplets;

% ---------- options ----------
EPS_PX   = 20;      % collar around the crystal mask (px)
FLIP_Y   = false;
DRAW_MAP = true;

TM_VEL_N  = 3;                      % half-width for the speed used in colour
CMAP_NAME = 'turbo';  CMAP_TRIM = 0.78;
CLIM_LO = 5;  CLIM_HI = 90;         % colour limits, as percentiles
LW_CAP = 1.3;  LW_FREE = 0.35;  ALPHA_FREE = 0.22;

colFree = [0.72 0.72 0.72];
colCrys = [0.10 0.60 0.30];         % green = crystal, as everywhere else
colDot  = [0.90 0.05 0.05];
CMAP    = getCmap(CMAP_NAME, CMAP_TRIM);

inDir  = cfg.workDir;
assert(isfolder(inDir), '%s not found', inDir);
outDir = cfg.outDir;
if ~isfolder(outDir), mkdir(outDir); end
fprintf('Data: %s\nOutput: %s\n\n', inDir, outDir);

nm = {}; nEnd = []; nTrap = [];  allC = {};

for d = 1:numel(droplets)
    dDir = droplets{d};
    [~, dName] = fileparts(dDir);
    resF = fullfile(dDir, 'RESULTS_crystals');
    fprintf('=== %s ===\n', dName);
    tG = tic;

    % ---------- inputs ----------
    f2 = fullfile(inDir, ['tmclean2_' dName '.mat']);
    if ~isfile(f2), warning('%s is missing', f2); continue; end
    S = load(f2, 'Tv2');  Tv = S.Tv2;

    bf = fullfile(inDir, ['buffers_' dName '.mat']);
    if ~isfile(bf), warning('%s is missing', bf); continue; end
    B = load(bf, 'buf');  buf = B.buf;
    Ny = buf.imgSize(1);  Nx = buf.imgSize(2);
    NC = buf.NCrys;
    dropMask = [];
    if isfield(buf,'dropPoly') && ~isempty(buf.dropPoly)
        dropMask = poly2mask(buf.dropPoly(:,1), buf.dropPoly(:,2), Ny, Nx);
    end

    kf = dir(fullfile(resF, 'individual_crystals_*.mat'));
    if isempty(kf), warning('individual_crystals_*.mat is missing'); continue; end
    K = load(fullfile(resF, kf(1).name), 'crystalMasks');
    cmask = K.crystalMasks;
    nMF = size(cmask,1);

    evalc('specs = readSpecs(dDir);');
    fps = specs.fps;  upp = specs.ppf;

    % ---------- index the spots and compute speeds ----------
    % The centred difference is clamped to each track's own first and last
    % row, so it never reaches across into a neighbouring track.
    [uid, ~, tidx] = unique(Tv.TrackID);
    nT = numel(uid);
    fr = Tv.Frame;  xx = Tv.X_px;  yy = Tv.Y_px;
    [~, o] = sortrows([tidx fr]);
    fr = fr(o); xx = xx(o); yy = yy(o); tidx = tidx(o);
    ng = [true; diff(tidx) ~= 0];
    rs = find(ng);  re = [rs(2:end)-1; numel(fr)];
    S0 = rs(tidx);  E0 = re(tidx);  j = (1:numel(fr))';
    ip = max(S0, j-TM_VEL_N);  in = min(E0, j+TM_VEL_N);
    df = fr(in) - fr(ip);
    vv = hypot(xx(in)-xx(ip), yy(in)-yy(ip)) * upp ./ (df/fps);
    vv(df <= 0) = NaN;

    lastRow  = re;                    % row of each track's last spot
    firstRow = rs;                    % row of each track's first spot
    fprintf('  %d tracks | %d crystals\n', nT, NC);

    % ---------- state of the first and last spot of every track ----------
    endOwn  = zeros(nT,1);            % crystal whose buffer holds the endpoint
    endStat = zeros(nT,1,'uint8');    % 0 fluid | 1 collar | 2 mask
    fstStat = zeros(nT,1,'uint8');

    % Only two rows per track are needed, so the loop runs over the frames
    % those rows fall in rather than over every frame of the sequence.
    rowsNeeded = [lastRow; firstRow];
    isLastRow  = [true(nT,1); false(nT,1)];
    tOfRow     = [(1:nT)'; (1:nT)'];
    frNeeded   = fr(rowsNeeded);

    for k = unique(frNeeded)'
        kk = k + 1;                                  % frame 0 -> index 1
        sel = find(frNeeded == k);
        if isempty(sel), continue; end
        rr = round(yy(rowsNeeded(sel)));  cc = round(xx(rowsNeeded(sel)));
        if FLIP_Y, rr = Ny - rr + 1; end
        ok = rr>=1 & rr<=Ny & cc>=1 & cc<=Nx;
        sel = sel(ok);  rr = rr(ok);  cc = cc(ok);
        if isempty(sel), continue; end
        lin = sub2ind([Ny Nx], rr, cc);

        % --- which buffer, in this frame ---
        if kk >= 1 && kk <= nMF
            O = buildBufferFrame(buf.centroids, buf.Rbuf_px, kk, Ny, Nx);
            if ~isempty(dropMask), O(~dropMask) = 0; end
        else
            O = zeros(Ny,Nx,'uint8');
        end
        ow = double(O(lin));

        % --- mask and collar, in this frame's crystal shape ---
        % Mask wins over collar: a point inside the crystal is state 2 and is
        % not downgraded by a neighbouring crystal's collar.
        st = zeros(numel(sel),1,'uint8');
        if kk >= 1 && kk <= nMF
            for c = 1:NC
                M = cmask{kk,c};
                if isempty(M) || ~any(M(:)), continue; end
                inM = M(lin);
                st(inM) = 2;
                Me = imdilate(M, strel('disk', EPS_PX));
                inE = Me(lin) & ~inM;
                st(inE & st == 0) = 1;
            end
        end

        isL = isLastRow(sel);  tt = tOfRow(sel);
        endOwn(tt(isL))   = ow(isL);
        endStat(tt(isL))  = st(isL);
        fstStat(tt(~isL)) = st(~isL);
    end

    % ---------- classification ----------
    endsInBuffer = endOwn > 0;
    trapped = endsInBuffer & (endStat >= 1);

    cnt = zeros(NC,4);
    for c = 1:NC
        cnt(c,1) = sum(endsInBuffer & endOwn == c);
        cnt(c,2) = sum(trapped      & endOwn == c);
        cnt(c,3) = sum(trapped      & endOwn == c & endStat == 2);
        cnt(c,4) = sum(trapped      & endOwn == c & endStat == 1);
    end
    pct = 100 * cnt(:,2) ./ max(cnt(:,1),1);

    Tc = table((1:NC)', cnt(:,1), cnt(:,2), pct, cnt(:,3), cnt(:,4), ...
        'VariableNames', {'Crystal','EndsInBuffer','Trapped','PercentTrapped', ...
                          'InsideMask','WithinEps'});
    disp(Tc);
    writetable(Tc, fullfile(outDir, sprintf('trap_percrystal_%s.csv', dName)));

    nE = sum(endsInBuffer);  nTr = sum(trapped);
    fprintf('  DROPLET: %d tracks end in a buffer, %d trapped (%.1f%%)\n', ...
            nE, nTr, 100*nTr/max(nE,1));
    fprintf('           of all tracks: %.1f%% in a buffer, %.1f%% trapped\n', ...
            100*nE/nT, 100*nTr/nT);

    % ---------- figure ----------
    if DRAW_MAP
        kLast = 0;
        for c = 1:NC
            lc = find(cellfun(@(m) ~isempty(m) && any(m(:)), cmask(:,c)),1,'last');
            if ~isempty(lc), kLast = max(kLast,lc); end
        end
        % Percentile colour limits, so one fast track does not compress the
        % rest into a single colour.
        sp = vv(isfinite(vv));
        climV = [prctile(sp,CLIM_LO) prctile(sp,CLIM_HI)];
        if climV(2) <= climV(1), climV = [0 max(sp)]; end

        fig = figure('Color','w','Position',[60 60 1000 1000]);
        ax = axes(fig); hold(ax,'on');
        ipth = findFrameImage(dDir, kLast);
        if ~isempty(ipth)
            % Converted to explicit RGB so the later colormap does not
            % repaint the greyscale background.
            I = imread(ipth);  if ndims(I) == 3, I = im2gray(I); end
            I = double(I);  a = min(I(:)); b = max(I(:));
            if b > a, ii = round((I-a)/(b-a)*255)+1; else, ii = ones(size(I)); end
            imagesc(ax, ind2rgb(ii, gray(256)));
        end
        set(ax,'YDir','normal'); axis(ax,'image');
        xlim(ax,[0 Nx]); ylim(ax,[0 Ny]);

        for c = 1:NC
            lc = find(cellfun(@(m) ~isempty(m) && any(m(:)), cmask(:,c)),1,'last');
            if isempty(lc), continue; end
            Bb = bwboundaries(cmask{lc,c});
            for q = 1:numel(Bb)
                plot(ax, Bb{q}(:,2), Bb{q}(:,1), '-', 'Color', colCrys, 'LineWidth', 1.8);
            end
        end

        % Free tracks first in flat grey, trapped ones on top coloured by
        % speed, so the trapped set reads clearly against the background.
        drawTrackSet(ax, xx, yy, tidx, vv, find(~trapped), [], colFree, ...
                     LW_FREE, ALPHA_FREE, []);
        drawTrackSet(ax, xx, yy, tidx, vv, find(trapped), climV, [], ...
                     LW_CAP, 1, CMAP);
        sel = lastRow(trapped);
        plot(ax, xx(sel), yy(sel), 'o', 'MarkerSize', 4.5, ...
             'MarkerFaceColor', colDot, 'MarkerEdgeColor','k', 'LineWidth', 0.4);

        cb = colorbar(ax); cb.Label.String = 'Track speed (\mum/s)';
        cb.TickLabelInterpreter = 'latex'; cb.Label.FontSize = 13;
        set(ax,'XTick',[],'YTick',[]); box(ax,'on');
        xlabel(ax, sprintf(['%d of %d tracks ending in a buffer are ' ...
               'trapped (%.1f\\%%)'], nTr, nE, 100*nTr/max(nE,1)), ...
               'Interpreter','latex','FontSize',14);
        exportgraphics(fig, fullfile(outDir, sprintf('trap_map_%s.png', dName)), ...
                       'Resolution',300);
        close(fig);
    end

    nm{end+1,1}    = dName;  %#ok<SAGROW>
    nEnd(end+1,1)  = nE;     %#ok<SAGROW>
    nTrap(end+1,1) = nTr;    %#ok<SAGROW>
    allC{end+1,1}  = Tc;     %#ok<SAGROW>
    fprintf('  time: %.1f s\n\n', toc(tG));
    clear Tv cmask K buf
end

% ---------- summary ----------
if isempty(nm), error('No droplet was processed.'); end
T = table(nm, nEnd, nTrap, 100*nTrap./max(nEnd,1), ...
    'VariableNames', {'Droplet','EndsInBuffer','Trapped','PercentTrapped'});
disp(T);
writetable(T, fullfile(outDir,'trap_summary.csv'));
Tall = vertcat(allC{:});
writetable(Tall, fullfile(outDir,'trap_percrystal_all.csv'));
fprintf('Done. Results in %s\n', outDir);
