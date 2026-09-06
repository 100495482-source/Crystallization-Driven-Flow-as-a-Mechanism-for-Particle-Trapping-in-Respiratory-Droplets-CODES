%% ========================================================================
%  ACCUM3_BUILDBUFFERZONES  Define a region of influence around each crystal.
%  ========================================================================
%
%  STAGE 3 OF 7 of the particle accumulation pipeline.
%
%  WHAT THIS DOES
%    For each droplet:
%      1. read centroids and radius_mean from individual_crystals_<S>.mat
%         (already measured by the crystal analysis; nothing is recomputed)
%      2. give each crystal a FIXED buffer radius: BUFFER_FACTOR times its
%         mean radius in the last frame in which it exists
%      3. in every frame, place that disc on the centroid of THAT frame, so
%         the zone follows the crystal without changing size
%      4. where two or more discs touch, split the shared region along the
%         common chord of each pair (see BUILDBUFFERFRAME)
%
%  WHY THE RADIUS IS FIXED AND THE CENTRE IS NOT
%    A buffer that grew with the crystal would change the reference area from
%    frame to frame, so a constant particle density would show up as a
%    changing one. Fixing the radius at its final value keeps the reference
%    region constant, and letting the centre follow the crystal keeps it in
%    the right place as the crystal shifts.
%
%  STORAGE
%    Storing the mask for every frame would run to gigabytes. Only the
%    PARAMETERS are saved (centroid per frame, radius per crystal), and
%    BUILDBUFFERFRAME reconstructs any frame instantly. The masks of a few
%    example frames are saved as well, for the figures.
%
%  OUTPUT  (in cfg.workDir)
%    buffers_<droplet>.mat   buf struct with the parameters and examples
%    buffers_<droplet>.png   example frames
%    buffers_summary.csv
%
%  See also BUILDBUFFERFRAME, ACCUM4_CLIPBUFFERSTODROPLET.

clear; close all; clc;

cfg = accumulationConfig;
droplets = cfg.droplets;

% ---------- options ----------
BUFFER_FACTOR = 2.5;    % buffer radius = factor * final mean crystal radius
DRAW_FRAMES   = 3;      % example frames in the figure
SAVE_EXAMPLES = true;   % store the masks of those frames

colC = lines(12);

outDir = cfg.workDir;
if ~isfolder(outDir), mkdir(outDir); end
fprintf('Output: %s\n\n', outDir);

nm = {}; nc = []; rr = [];

for d = 1:numel(droplets)
    dDir = droplets{d};
    [~, dName] = fileparts(dDir);
    resF = fullfile(dDir, 'RESULTS_crystals');
    fprintf('=== %s ===\n', dName);
    tG = tic;

    kf = dir(fullfile(resF, 'individual_crystals_*.mat'));
    if isempty(kf), warning('individual_crystals_*.mat is missing'); continue; end
    K = load(fullfile(resF, kf(1).name), 'centroids','radius_mean','crystalMasks');
    cen   = K.centroids;      % [frames x NCrys x 2], px
    rmn   = K.radius_mean;    % [frames x NCrys],     px
    cmask = K.crystalMasks;
    [nMF, NC] = size(rmn);
    iF = find(cellfun(@(m) ~isempty(m), cmask(:,1)), 1);
    [Ny, Nx] = size(cmask{iF,1});

    % evalc suppresses the specifications printout, which would otherwise
    % repeat for every droplet and bury the output of this stage.
    evalc('specs = readSpecs(dDir);');
    upp = specs.ppf;

    % ---------- fixed radius per crystal: from its last frame ----------
    Rbuf = nan(1, NC);  lastF = nan(1, NC);
    for c = 1:NC
        lc = find(isfinite(rmn(:,c)) & rmn(:,c) > 0, 1, 'last');
        if isempty(lc), continue; end
        lastF(c) = lc;
        Rbuf(c)  = BUFFER_FACTOR * rmn(lc,c);
    end
    fprintf('  %d crystals | buffer radii (um): %s\n', NC, ...
            strjoin(compose('%.0f', Rbuf*upp), ', '));

    % ---------- frames containing at least one crystal ----------
    okFrame = any(isfinite(rmn) & rmn > 0, 2);
    frIdxList = find(okFrame)';

    exFrames = round(linspace(frIdxList(1), frIdxList(end), DRAW_FRAMES));
    exMasks  = cell(numel(exFrames),1);

    % ---------- how often the buffers overlap, for information only ------
    % A high count is not an error, since overlaps are resolved by the power
    % diagram, but it indicates how much of the partition is contested.
    nOverlapFrames = 0;
    for k = frIdxList
        pres = find(isfinite(Rbuf(:)') & isfinite(reshape(cen(k,:,1),1,[])));
        if numel(pres) < 2, continue; end
        hit = false;
        for a = 1:numel(pres)-1
            for b = a+1:numel(pres)
                ca = pres(a); cb = pres(b);
                dd = hypot(cen(k,ca,1)-cen(k,cb,1), cen(k,ca,2)-cen(k,cb,2));
                if dd < Rbuf(ca) + Rbuf(cb), hit = true; break; end
            end
            if hit, break; end
        end
        if hit, nOverlapFrames = nOverlapFrames + 1; end
    end
    fprintf('  frames with overlapping buffers: %d of %d\n', ...
            nOverlapFrames, numel(frIdxList));

    % ---------- example masks ----------
    if SAVE_EXAMPLES
        for e = 1:numel(exFrames)
            exMasks{e} = buildBufferFrame(cen, Rbuf, exFrames(e), Ny, Nx);
        end
    end

    % ---------- save the parameters ----------
    buf.droplet   = dName;
    buf.NCrys     = NC;
    buf.Rbuf_px   = Rbuf;
    buf.Rbuf_um   = Rbuf * upp;
    buf.centroids = cen;          % [frames x NCrys x 2]
    buf.lastFrame = lastF;
    buf.factor    = BUFFER_FACTOR;
    buf.imgSize   = [Ny Nx];
    buf.exFrames  = exFrames;
    buf.exMasks   = exMasks;      % uint8: 0 outside, c = crystal c
    save(fullfile(outDir, ['buffers_' dName '.mat']), 'buf', '-v7.3');

    % ---------- figure ----------
    fig = figure('Color','w','Position',[60 60 520*numel(exFrames) 560]);
    tl = tiledlayout(fig, 1, numel(exFrames), ...
                     'TileSpacing','compact','Padding','compact');
    th = linspace(0, 2*pi, 300);
    for e = 1:numel(exFrames)
        k = exFrames(e);
        ax = nexttile(tl); hold(ax,'on')
        O = buildBufferFrame(cen, Rbuf, k, Ny, Nx);

        % assigned zones, one colour per crystal
        RGB = ones(Ny,Nx,3);
        for c = 1:NC
            m = O == c;
            if ~any(m(:)), continue; end
            for ch = 1:3
                tmp = RGB(:,:,ch);  tmp(m) = colC(mod(c-1,12)+1, ch);
                RGB(:,:,ch) = tmp;
            end
        end
        imagesc(ax, RGB);
        set(ax,'YDir','normal'); axis(ax,'image');
        xlim(ax,[0 Nx]); ylim(ax,[0 Ny]);

        % Full circles dashed on top, so it is visible where each disc was
        % cut back by a neighbour, plus the crystal outline itself.
        for c = 1:NC
            if ~isfinite(Rbuf(c)) || ~isfinite(cen(k,c,1)), continue; end
            plot(ax, cen(k,c,1)+Rbuf(c)*cos(th), cen(k,c,2)+Rbuf(c)*sin(th), ...
                 'k--', 'LineWidth', 0.9);
            plot(ax, cen(k,c,1), cen(k,c,2), 'k+', 'MarkerSize', 8, 'LineWidth', 1.3);
            M = cmask{k,c};
            if ~isempty(M) && any(M(:))
                B = bwboundaries(M);
                for q = 1:numel(B)
                    plot(ax, B{q}(:,2), B{q}(:,1), 'k-', 'LineWidth', 1.4);
                end
            end
        end
        set(ax,'XTick',[],'YTick',[]); box(ax,'on');
        title(ax, sprintf('Frame %d', k), 'Interpreter','latex','FontSize',14);
    end
    exportgraphics(fig, fullfile(outDir, sprintf('buffers_%s.png', dName)), ...
                   'Resolution',300);
    close(fig);

    nm{end+1,1} = dName;                    %#ok<SAGROW>
    nc(end+1,1) = NC;                       %#ok<SAGROW>
    rr(end+1,1) = mean(Rbuf*upp,'omitnan'); %#ok<SAGROW>
    fprintf('  saved (%.1f s)\n\n', toc(tG));
end

% ---------- summary ----------
if isempty(nm), error('No droplet was processed.'); end
T = table(nm, nc, rr, ...
    'VariableNames', {'Droplet','NCrystals','MeanBufferRadius_um'});
disp(T);
writetable(T, fullfile(outDir,'buffers_summary.csv'));
fprintf('Done. Results in %s\n', outDir);
