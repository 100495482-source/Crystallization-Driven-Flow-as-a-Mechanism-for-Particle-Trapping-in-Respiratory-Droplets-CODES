%% ========================================================================
%  ACCUM7_RADIALDENSITY  Endpoint density against distance from the crystal.
%  ========================================================================
%
%  STAGE 7 OF 7 of the particle accumulation pipeline, and the final result.
%
%  WHAT THIS MEASURES
%    Whether particles pile up near the crystals, expressed as the density of
%    track endpoints in bands of increasing distance from the crystal edge,
%    normalised by what a uniform distribution would give. A value of 1 means
%    no preference; values above 1 near the edge mean accumulation.
%
%  ONE POINT PER PARTICLE
%    Each particle contributes exactly ONE point: its last known position, in
%    the frame where it stops being visible. Counting every detection instead
%    would weight slow particles far more heavily than fast ones and produce
%    points that are not independent, since consecutive detections of the same
%    particle are highly correlated.
%
%  THE GEOMETRY CHANGES WITH TIME, AND IS TREATED AS SUCH
%    - the distance of each endpoint is measured against the mask of ITS OWN
%      frame, that is, against the crystal as it was when that particle
%      stopped, not against the final crystal
%    - the area of each band is likewise recomputed per frame, and the uniform
%      model is accumulated weighted by how many endpoints that frame
%      contributed
%
%    Measuring everything against the final crystal would push early endpoints
%    into bands that did not exist yet when those particles stopped.
%
%  PERFORMANCE
%    Crystals grow slowly, so the band areas change little between consecutive
%    frames. They are recomputed every AREA_STRIDE frames and reused in
%    between, which is what makes this tractable.
%
%  OUTPUT  (in cfg.outDir)
%    bandsEP_<droplet>.png     the bands, illustrative
%    histEP_<droplet>.png      raw endpoint counts per crystal
%    densityEP_<droplet>.png   normalised density per crystal
%    densityEP_ALL.png         all crystals, with the mean and its error band
%    densityEP_all.csv
%
%  See also ACCUM6_REMOVETRACKSDYINGINSIDE, BUILDBUFFERFRAME.

clear; close all; clc;

cfg = accumulationConfig;
droplets = cfg.droplets;

% ---------- options ----------
BIN_UM        = 20;      % band width (um)
AREA_STRIDE   = 5;       % recompute band areas every N frames
MIN_AREA_FRAC = 0.35;    % merge the last band if its area is below this
MIN_EXPEC     = 4;       % skip bands expecting fewer than N endpoints.
                         %  Avoids meaningless spikes: a band with a tiny
                         %  area containing a single endpoint yields a
                         %  density of 8 or 24 that means nothing.
MIN_CRYS_MEAN = 3;       % in the global figure, minimum crystals per band
FLIP_Y        = false;
DRAW_ZONES    = true;

colD    = [0.85 0.40 0.20];
colCrys = [0.10 0.60 0.30];

inDir = cfg.workDir;
assert(isfolder(inDir), '%s not found', inDir);
outDir = cfg.outDir;
if ~isfolder(outDir), mkdir(outDir); end
fprintf('Output: %s\n\n', outDir);

ALL = struct([]);

for d = 1:numel(droplets)
    dDir = droplets{d};
    [~, dName] = fileparts(dDir);
    resF = fullfile(dDir, 'RESULTS_crystals');
    fprintf('=== %s ===\n', dName);
    tG = tic;

    f3 = fullfile(inDir, ['tmclean3_' dName '.mat']);
    if ~isfile(f3), warning('%s is missing', f3); continue; end
    S = load(f3, 'Tv3');  Tv = S.Tv3;

    bf = fullfile(inDir, ['buffers_' dName '.mat']);
    if ~isfile(bf), warning('%s is missing', bf); continue; end
    B = load(bf, 'buf');  buf = B.buf;
    Ny = buf.imgSize(1);  Nx = buf.imgSize(2);  NC = buf.NCrys;
    % Without the outline the reference areas would include dry substrate,
    % so this is a hard requirement rather than a warning.
    assert(isfield(buf,'dropPoly') && ~isempty(buf.dropPoly), ...
           '%s has no droplet outline drawn. Run stage 4 first.', dName);
    dropMask = poly2mask(buf.dropPoly(:,1), buf.dropPoly(:,2), Ny, Nx);

    kf = dir(fullfile(resF,'individual_crystals_*.mat'));
    if isempty(kf), warning('individual_crystals_*.mat is missing'); continue; end
    K = load(fullfile(resF, kf(1).name), 'crystalMasks');
    cmask = K.crystalMasks;
    nMF = size(cmask,1);

    evalc('specs = readSpecs(dDir);');
    upp = specs.ppf;

    cen  = buf.centroids;
    Rbuf = buf.Rbuf_px;

    % ---------- band grid, common to the whole droplet ----------
    hiD = max(Rbuf) * upp;
    ed  = 0:BIN_UM:(ceil(hiD/BIN_UM)*BIN_UM);
    ctr = ed(1:end-1) + BIN_UM/2;
    nB  = numel(ctr);

    obs    = zeros(nB, NC);   % observed endpoints per band
    expec  = zeros(nB, NC);   % endpoints expected if uniform
    aCache = cell(1, NC);     % cached band areas, refreshed every AREA_STRIDE
    nSpotsTot = 0;

    % ---------- one point per track: its last position ----------
    [uid, ~, tix] = unique(Tv.TrackID);
    nTk = numel(uid);
    [~, o] = sortrows([tix Tv.Frame]);
    frS = Tv.Frame(o);  xS = Tv.X_px(o);  yS = Tv.Y_px(o);  tS = tix(o);
    ng  = [true; diff(tS) ~= 0];
    lastRow = [find(ng(2:end)); numel(tS)];
    Tv = table(uid(tS(lastRow)), frS(lastRow), xS(lastRow), yS(lastRow), ...
               'VariableNames', {'TrackID','Frame','X_px','Y_px'});

    [ufr, ~, fidx] = unique(Tv.Frame);
    fprintf('  %d tracks -> %d endpoints in %d frames\n', ...
            nTk, height(Tv), numel(ufr));

    cnt = 0;
    for q = 1:numel(ufr)
        k  = ufr(q);
        kk = k + 1;
        if kk < 1 || kk > nMF, continue; end
        rows = find(fidx == q);
        rr = round(Tv.Y_px(rows));  cc = round(Tv.X_px(rows));
        if FLIP_Y, rr = Ny - rr + 1; end
        ok = rr>=1 & rr<=Ny & cc>=1 & cc<=Nx;
        rows = rows(ok);  rr = rr(ok);  cc = cc(ok);
        if isempty(rows), continue; end
        px = Tv.X_px(rows);  py = Tv.Y_px(rows);

        % --- assign each endpoint to a buffer, analytically ---
        % Same power criterion as buildBufferFrame, but evaluated at the
        % handful of endpoint coordinates rather than over the whole image.
        cxk = reshape(cen(kk,:,1),1,[]);  cyk = reshape(cen(kk,:,2),1,[]);
        pres = find(isfinite(Rbuf(:)') & isfinite(cxk) & isfinite(cyk));
        if isempty(pres), continue; end
        best = inf(numel(rows),1);  own = zeros(numel(rows),1);
        for c = pres
            d2 = (px-cxk(c)).^2 + (py-cyk(c)).^2;
            pw = d2 - Rbuf(c)^2;
            m  = d2 <= Rbuf(c)^2 & pw < best;
            best(m) = pw(m);  own(m) = c;
        end
        inDrop = dropMask(sub2ind([Ny Nx], rr, cc));
        own(~inDrop) = 0;

        doArea = (mod(cnt, AREA_STRIDE) == 0);
        cnt = cnt + 1;

        for c = pres
            mine = own == c;
            M = cmask{kk,c};
            if isempty(M) || ~any(M(:)), continue; end
            if ~any(mine) && ~doArea, continue; end

            % --- bounding box of this crystal's buffer ---
            c0 = max(1, floor(cxk(c)-Rbuf(c)));  c1 = min(Nx, ceil(cxk(c)+Rbuf(c)));
            r0 = max(1, floor(cyk(c)-Rbuf(c)));  r1 = min(Ny, ceil(cyk(c)+Rbuf(c)));
            if c1 < c0 || r1 < r0, continue; end
            Mc = M(r0:r1, c0:c1);
            if ~any(Mc(:)), continue; end
            % Distance transform from the crystal outline, in um. This is why
            % the bands follow the crystal shape rather than being annuli.
            Dc = bwdist(Mc) * upp;

            % --- bin this crystal's endpoints ---
            if any(mine)
                sr = rr(mine) - r0 + 1;  sc = cc(mine) - c0 + 1;
                inb = sr>=1 & sr<=size(Dc,1) & sc>=1 & sc<=size(Dc,2);
                if any(inb)
                    ds = Dc(sub2ind(size(Dc), sr(inb), sc(inb)));
                    onMask = Mc(sub2ind(size(Mc), sr(inb), sc(inb)));
                    ds = ds(~onMask);    % points on the mask have no distance
                    bI = discretize(ds, ed);
                    bI = bI(~isnan(bI));
                    if ~isempty(bI)
                        obs(:,c) = obs(:,c) + accumarray(bI, 1, [nB 1]);
                    end
                end
            end

            % --- band areas, refreshed every AREA_STRIDE frames ---
            % The zone is this crystal's buffer, minus its own mask, minus the
            % part of the droplet outline that excludes it, minus whatever a
            % neighbouring crystal wins under the same power criterion.
            if doArea || isempty(aCache{c})
                [XX,YY] = meshgrid(c0:c1, r0:r1);
                dcen = (XX-cxk(c)).^2 + (YY-cyk(c)).^2;
                zone = dcen <= Rbuf(c)^2 & ~Mc & dropMask(r0:r1, c0:c1);
                for c2 = pres
                    if c2 == c, continue; end
                    d22 = (XX-cxk(c2)).^2 + (YY-cyk(c2)).^2;
                    zone = zone & ~(d22 <= Rbuf(c2)^2 & ...
                                    (d22 - Rbuf(c2)^2) < (dcen - Rbuf(c)^2));
                end
                if any(zone(:))
                    aCache{c} = histcounts(Dc(zone), ed)';
                end
            end
            % The uniform expectation for this frame: this frame's endpoints
            % spread over this frame's band areas.
            if ~isempty(aCache{c}) && any(mine)
                aB = aCache{c};
                if sum(aB) > 0
                    expec(:,c) = expec(:,c) + sum(mine) * (aB / sum(aB));
                end
            end
        end
        nSpotsTot = nSpotsTot + sum(own > 0);
    end
    fprintf('  endpoints assigned to a buffer: %d\n', nSpotsTot);

    % ---------- merge the outermost band if it is small ----------
    % The last band is usually a thin sliver where the disc is cut by the
    % droplet outline, and on its own it is too noisy to report.
    edU = ed;  ctrU = ctr;
    typ = median(expec(expec > 0));
    for c = 1:NC
        if expec(end,c) > 0 && expec(end,c) < MIN_AREA_FRAC*typ
            expec(end-1,c) = expec(end-1,c) + expec(end,c);
            obs(end-1,c)   = obs(end-1,c)   + obs(end,c);
            expec(end,c) = 0;  obs(end,c) = 0;
        end
    end

    dens = obs ./ max(expec, eps);
    dens(expec == 0) = NaN;
    nDrop = sum(expec(:) > 0 & expec(:) < MIN_EXPEC);
    dens(expec < MIN_EXPEC) = NaN;        % bands without enough area
    if nDrop > 0
        fprintf('  %d bands discarded for insufficient area\n', nDrop);
    end

    % ---------- FIGURE 1: the bands, illustrative ----------
    kLast = 0;
    for c = 1:NC
        lc = find(cellfun(@(m) ~isempty(m) && any(m(:)), cmask(:,c)),1,'last');
        if ~isempty(lc), kLast = max(kLast,lc); end
    end
    if DRAW_ZONES && kLast > 0
        O = buildBufferFrame(cen, Rbuf, kLast, Ny, Nx);
        O(~dropMask) = 0;
        zoneImg = zeros(Ny,Nx,'uint8');
        for c = 1:NC
            M = cmask{kLast,c};
            if isempty(M) || ~any(M(:)), continue; end
            Dum = bwdist(M)*upp;
            z = (O==c) & ~M;
            for b = 1:nB
                band = z & Dum >= ed(b) & Dum < ed(b+1);
                zoneImg(band) = b;
            end
        end
        nb = double(max(zoneImg(:)));
        cmapB = parula(max(nb,2));
        RGB = ones(Ny,Nx,3);
        for b = 1:nb
            m = zoneImg == b;
            if ~any(m(:)), continue; end
            for ch = 1:3
                t = RGB(:,:,ch); t(m) = cmapB(b,ch); RGB(:,:,ch) = t;
            end
        end
        fig = figure('Color','w','Position',[60 60 950 950]);
        ax = axes(fig); hold(ax,'on');
        imagesc(ax, RGB);
        set(ax,'YDir','normal'); axis(ax,'image');
        xlim(ax,[0 Nx]); ylim(ax,[0 Ny]);
        for c = 1:NC
            M = cmask{kLast,c};
            if isempty(M) || ~any(M(:)), continue; end
            Bb = bwboundaries(M);
            for qq = 1:numel(Bb)
                plot(ax, Bb{qq}(:,2), Bb{qq}(:,1), '-', 'Color', colCrys, 'LineWidth',1.6);
            end
        end
        set(ax,'XTick',[],'YTick',[]); box(ax,'on');
        xlabel(ax, sprintf('%d bands of %g $\\mu$m from the crystal edge', nb, BIN_UM), ...
               'Interpreter','latex','FontSize',14);
        exportgraphics(fig, fullfile(outDir, sprintf('bandsEP_%s.png', dName)), ...
                       'Resolution',300);
        close(fig);
    end

    % ---------- FIGURE 2: raw counts per crystal ----------
    % Shown alongside the density so the reader can see how many endpoints
    % each normalised value rests on.
    have = find(sum(obs,1) > 0);
    nCr = numel(have);
    if nCr > 0
        nCol = min(4,nCr);  nRow = ceil(nCr/nCol);
        fig = figure('Color','w','Position',[60 60 340*nCol 300*nRow]);
        tl = tiledlayout(fig,nRow,nCol,'TileSpacing','compact','Padding','compact');
        for i = 1:nCr
            ax = nexttile(tl);
            bar(ax, ctr, obs(:,have(i)), 1, 'FaceColor', colD, 'EdgeColor','none');
            grid(ax,'on'); box(ax,'on')
            set(ax,'TickLabelInterpreter','latex','FontSize',10);
        end
        xlabel(tl,'Distance to crystal edge ($\mu$m)','Interpreter','latex','FontSize',14);
        ylabel(tl,'Track endpoints','Interpreter','latex','FontSize',14);
        exportgraphics(fig, fullfile(outDir, sprintf('histEP_%s.png', dName)), ...
                       'Resolution',300);
        close(fig);
    end

    % ---------- FIGURE 3: normalised density per crystal ----------
    fig = figure('Color','w','Position',[60 60 900 580]); hold on
    lg = cell(0,1);
    for i = 1:nCr
        plot(ctr, dens(:,have(i)), '-o', 'LineWidth',1.6, 'MarkerSize',4);
        lg{end+1,1} = sprintf('Crystal %d', have(i)); %#ok<SAGROW>
    end
    yline(1,'k--','LineWidth',1.4);      % uniform reference
    grid on; box on
    set(gca,'TickLabelInterpreter','latex','FontSize',13);
    xlabel('Distance to crystal edge ($\mu$m)','Interpreter','latex','FontSize',15);
    ylabel('Endpoint density relative to uniform','Interpreter','latex','FontSize',15);
    if ~isempty(lg), legend(lg,'Interpreter','latex','Location','northeast'); end
    exportgraphics(fig, fullfile(outDir, sprintf('densityEP_%s.png', dName)), ...
                   'Resolution',300);
    close(fig);

    for i = 1:nCr
        ALL(end+1).droplet = dName; %#ok<SAGROW>
        ALL(end).c    = have(i);
        ALL(end).ctr  = ctr(:);
        ALL(end).dens = dens(:,have(i));
        ALL(end).n    = sum(obs(:,have(i)));
    end
    fprintf('  time: %.1f s\n\n', toc(tG));
    clear Tv cmask K buf
end

%% ========================================================================
%  GLOBAL RESULT: every crystal of every droplet on one axis
%  ========================================================================
if isempty(ALL), error('No droplet was processed.'); end

% Droplets have different buffer radii, so their band grids differ in length.
% They are resampled onto a common grid with 'nearest' rather than a linear
% interpolation, so no band value is invented between measured points.
hiA = 0;
for a = 1:numel(ALL), hiA = max(hiA, max(ALL(a).ctr)); end
ctrG = (BIN_UM/2):BIN_UM:hiA;
MG = nan(numel(ALL), numel(ctrG));
for a = 1:numel(ALL)
    MG(a,:) = interp1(ALL(a).ctr, ALL(a).dens, ctrG, 'nearest', NaN);
end
writetable([table({ALL.droplet}', [ALL.c]', [ALL.n]', ...
                  'VariableNames',{'Droplet','Crystal','Endpoints'}) ...
            array2table(MG,'VariableNames', ...
            matlab.lang.makeValidName(compose('d_%g',ctrG)))], ...
           fullfile(outDir,'densityEP_all.csv'));

fig = figure('Color','w','Position',[60 60 1000 640]); hold on
% Individual crystals in light grey behind the mean, so the spread is visible
% rather than hidden by the average.
for r = 1:size(MG,1)
    plot(ctrG, MG(r,:), '-', 'LineWidth',0.8, 'Color',[0.78 0.78 0.78]);
end
mu  = mean(MG,1,'omitnan');
n   = sum(isfinite(MG),1);
sem = std(MG,0,1,'omitnan') ./ max(sqrt(n),1);
% Only bands backed by at least MIN_CRYS_MEAN crystals are drawn: the far
% bands exist for only a few crystals and their mean would be unreliable.
ok  = isfinite(mu) & n >= MIN_CRYS_MEAN;
fill([ctrG(ok) fliplr(ctrG(ok))], [mu(ok)+sem(ok) fliplr(mu(ok)-sem(ok))], ...
     colD,'FaceAlpha',0.25,'EdgeColor','none');
plot(ctrG(ok), mu(ok), '-o','LineWidth',2.4,'Color',colD,'MarkerSize',5);
yline(1,'k--','LineWidth',1.4);
xmax = max(ctrG(ok));
if ~isempty(xmax), xlim([0 xmax + BIN_UM]); end
grid on; box on
set(gca,'TickLabelInterpreter','latex','FontSize',13);
xlabel('Distance to crystal edge ($\mu$m)','Interpreter','latex','FontSize',15);
ylabel('Endpoint density relative to uniform','Interpreter','latex','FontSize',15);
exportgraphics(fig, fullfile(outDir,'densityEP_ALL.png'),'Resolution',300);

fprintf('%d crystals from %d droplets\nDone. Results in %s\n', ...
        numel(ALL), numel(unique({ALL.droplet})), outDir);
