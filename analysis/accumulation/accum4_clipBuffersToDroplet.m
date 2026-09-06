%% ========================================================================
%  ACCUM4_CLIPBUFFERSTODROPLET  Cut the buffer zones to the droplet outline.
%  ========================================================================
%
%  STAGE 4 OF 7 of the particle accumulation pipeline.
%
%  WHAT THIS DOES
%    For each droplet:
%      1. show the FIRST frame
%      2. you trace the droplet outline by hand
%      3. that region is intersected with the buffer zones
%      4. the clipped buffers are shown again for checking
%
%  WHY THIS IS NECESSARY
%    A buffer is a disc, and a crystal near the contact line has part of its
%    disc lying outside the droplet, over dry substrate where no particle can
%    ever be. Counting that area as available space would make the reference
%    density too low and manufacture an apparent accumulation. The intersection
%    removes it.
%
%    The outline is drawn on the first frame because the contact line is
%    pinned: the droplet's footprint does not move as it dries, only its
%    height falls. One outline is therefore valid for the whole sequence.
%
%    The outline is stored inside buffers_<droplet>.mat (field dropPoly), so
%    it only has to be drawn ONCE per droplet. Set REDRAW = true to redo it.
%
%  OUTPUT  (in cfg.workDir)
%    buffers_<droplet>.mat updated with dropPoly, clipped exMasks and
%                          the clipped flag
%    buffers_clipped_<droplet>.png
%
%  See also ACCUM3_BUILDBUFFERZONES, BUILDBUFFERFRAME, FINDFRAMEIMAGE.

clear; close all; clc;

cfg = accumulationConfig;
droplets = cfg.droplets;

REDRAW = false;     % true = redraw even if an outline is already stored
colC   = lines(12);

outDir = cfg.workDir;
assert(isfolder(outDir), '%s not found', outDir);

for d = 1:numel(droplets)
    dDir = droplets{d};
    [~, dName] = fileparts(dDir);
    fprintf('=== %s ===\n', dName);

    bf = fullfile(outDir, ['buffers_' dName '.mat']);
    if ~isfile(bf), warning('%s is missing', bf); continue; end
    S = load(bf, 'buf');  buf = S.buf;
    Ny = buf.imgSize(1);  Nx = buf.imgSize(2);

    % ---------- 1. droplet outline ----------
    if isfield(buf,'dropPoly') && ~isempty(buf.dropPoly) && ~REDRAW
        fprintf('  outline already stored, reusing it\n');
    else
        ip = findFrameImage(dDir, 1);
        assert(~isempty(ip), 'First image of %s not found', dDir);
        I = imread(ip);
        if ndims(I) == 3, I = im2gray(I); end

        fig = figure('Color','w','Position',[80 60 950 950]);
        ax = axes(fig);
        imshow(I, [], 'Parent', ax); hold(ax,'on');
        set(ax,'YDir','normal');
        title(ax, sprintf('%s: trace the droplet outline and press Enter', dName), ...
              'FontSize', 14);
        % Interactive modes are switched off so a stray click does not zoom
        % the axes in the middle of drawing the polygon.
        zoom(fig,'off'); pan(fig,'off'); datacursormode(fig,'off');
        fprintf(['  Click around the droplet edge (12-20 points is enough).\n' ...
                 '  Double-click, or click the first point, to close it.\n' ...
                 '  Vertices can be dragged to correct. Then press ENTER.\n']);
        h = drawpolygon(ax, 'Color', [1 1 0], 'LineWidth', 1.5, 'FaceAlpha', 0.10);
        wait(h);
        buf.dropPoly = h.Position;      % [N x 2] in (x, y)
        fprintf('  outline: %d vertices\n', size(buf.dropPoly,1));
        close(fig);
        save(bf, 'buf', '-v7.3');
        fprintf('  outline saved (%d points)\n', size(buf.dropPoly,1));
    end

    dropMask = poly2mask(buf.dropPoly(:,1), buf.dropPoly(:,2), Ny, Nx);
    fprintf('  droplet: %.1f%% of the image\n', 100*nnz(dropMask)/(Ny*Nx));

    % ---------- 2. clipped buffers ----------
    exF = buf.exFrames;
    exM = cell(numel(exF),1);
    for e = 1:numel(exF)
        O = buildBufferFrame(buf.centroids, buf.Rbuf_px, exF(e), Ny, Nx);
        O(~dropMask) = 0;               % intersection with the droplet
        exM{e} = O;
    end
    buf.exMasks = exM;
    buf.clipped = true;
    save(bf, 'buf', '-v7.3');

    % ---------- how much buffer area the clip removed, for information ----
    O0 = buildBufferFrame(buf.centroids, buf.Rbuf_px, exF(end), Ny, Nx);
    lost = 100 * (nnz(O0>0) - nnz(exM{end}>0)) / max(nnz(O0>0),1);
    fprintf('  buffer area clipped away in the last example: %.1f%%\n', lost);

    % ---------- 3. figure ----------
    fig = figure('Color','w','Position',[60 60 520*numel(exF) 580]);
    tl = tiledlayout(fig, 1, numel(exF), 'TileSpacing','compact','Padding','compact');
    th = linspace(0, 2*pi, 300);
    Bd = bwboundaries(dropMask);
    for e = 1:numel(exF)
        k = exF(e);
        ax = nexttile(tl); hold(ax,'on')
        O = exM{e};
        RGB = ones(Ny,Nx,3);
        for c = 1:buf.NCrys
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
        for q = 1:numel(Bd)
            plot(ax, Bd{q}(:,2), Bd{q}(:,1), 'k-', 'LineWidth', 1.6);
        end
        cxk = reshape(buf.centroids(k,:,1),1,[]);
        cyk = reshape(buf.centroids(k,:,2),1,[]);
        for c = 1:buf.NCrys
            if ~isfinite(buf.Rbuf_px(c)) || ~isfinite(cxk(c)), continue; end
            plot(ax, cxk(c)+buf.Rbuf_px(c)*cos(th), ...
                     cyk(c)+buf.Rbuf_px(c)*sin(th), 'k--', 'LineWidth', 0.9);
            plot(ax, cxk(c), cyk(c), 'k+', 'MarkerSize', 8, 'LineWidth', 1.3);
        end
        set(ax,'XTick',[],'YTick',[]); box(ax,'on');
        title(ax, sprintf('Frame %d', k), 'Interpreter','latex','FontSize',14);
    end
    exportgraphics(fig, fullfile(outDir, sprintf('buffers_clipped_%s.png', dName)), ...
                   'Resolution',300);
    fprintf('  figure saved\n\n');
end

fprintf('Done. Results in %s\n', outDir);
