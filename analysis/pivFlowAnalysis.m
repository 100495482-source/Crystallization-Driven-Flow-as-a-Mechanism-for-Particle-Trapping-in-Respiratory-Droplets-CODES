%% ========================================================================
%  PIVFLOWANALYSIS  Eulerian flow measurement from PIVlab exports.
%  ========================================================================
%
%  WHAT THIS SCRIPT DOES
%    Turns a raw PIVlab export into calibrated velocity fields and into the
%    time series the thesis actually reports. It runs in three parts, each
%    behind its own flag in runDropletPipeline:
%
%      DO_PIV_FIELDS       read the PIVlab .mat, convert px/frame to um/s,
%                          blank out the velocities that fall inside the
%                          crystals, average in time, and derive vorticity,
%                          divergence and speed
%      DO_PIV_1DPLOTS      collapse each field to one number per frame
%                          (spatial mean, maximum and standard deviation)
%                          and plot it against real time
%      DO_VIDEO_*          render any of the fields over the raw photographs
%
%    This is a SCRIPT, not a function: it is called with run() from
%    runDropletPipeline and shares its workspace.
%
%  WHY THE CRYSTAL IS MASKED OUT
%    PIV cross-correlates image texture. The crystal surface has plenty of
%    texture, so PIVlab happily reports "velocities" inside it that are not
%    fluid motion at all. Those nodes are set to NaN before any statistic is
%    computed.
%
%  INPUTS FROM THE PIPELINE
%    fps, um_per_px, dropletDir, dropletName, resultsFolder, photoDir,
%    frameTimes, nucTimes, frame_min, frame_max, SMOOTH_WINDOW,
%    USE_CRYSTAL_MASK, K_AVG, and 'masks' when USE_CRYSTAL_MASK is true.
%
%  OUTPUT FILES  (in resultsFolder)
%    piv_fields_<droplet>.mat        all calibrated fields
%    OneDPlots_<quantity>_<droplet>.mat / .png
%    <quantity>_field_<droplet>.mp4  when the video flags are on
%
%  See also RUNDROPLETPIPELINE, TRACKFLOWANALYSIS, FIGURESTYLE.

%% ========================================================================
%  PART A : FIELD COMPUTATION
%  ========================================================================
if DO_PIV_FIELDS

    % ---------- locate the PIVlab export in the droplet folder ----------
    sessionDir = dropletDir;
    pivMat = dir(fullfile(sessionDir, 'PIVlab*.mat'));
    assert(~isempty(pivMat), 'No PIVlab*.mat found in %s', sessionDir);
    pivMatPath = fullfile(sessionDir, pivMat(1).name);
    fprintf('PIVlab: %s\n', pivMatPath);
    load(pivMatPath, 'u_original', 'v_original', 'x', 'y');

    % ---------- calibration, always from the specifications ----------
    % PIVlab returns displacements in px/frame; both factors come from the
    % specifications file rather than being hard-coded, because the optical
    % scale changed between sessions.
    scale_um_s = um_per_px * fps;    % px/frame -> um/s
    fprintf('PIV calibration: um_per_px=%.4f  fps=%g  scale=%.4f um/s\n', ...
            um_per_px, fps, scale_um_s);

    % ---------- blank the nodes lying inside a crystal ----------
    u_masked = cell(size(u_original));
    v_masked = cell(size(v_original));

    for k = 1:numel(u_original)
        if isempty(u_original{k})
            u_masked{k} = u_original{k};
            v_masked{k} = v_original{k};
            continue
        end
        u = u_original{k};
        v = v_original{k};

        if USE_CRYSTAL_MASK && exist('masks','var') && k <= numel(masks) ...
                            && ~isempty(masks{k}) && any(masks{k}(:)) ...
                            && k <= numel(x) && ~isempty(x{k})
            % The PIV grid does NOT cover the whole image: it starts half an
            % interrogation window inside the frame and advances in steps.
            % Resizing the mask onto the grid therefore stretched it over the
            % full field and left the crystal offset from the vectors it was
            % meant to blank. The correct approach is to sample the mask at
            % the real pixel coordinates (x{k}, y{k}) of each PIV node.
            M   = masks{k};
            col = round(x{k});   row = round(y{k});
            inb = col >= 1 & col <= size(M,2) & row >= 1 & row <= size(M,1);
            mask_piv = false(size(u));
            mask_piv(inb) = M(sub2ind(size(M), row(inb), col(inb)));
            u(mask_piv) = NaN;
            v(mask_piv) = NaN;
        end

        u_masked{k} = u;
        v_masked{k} = v;
    end

    % ---------- velocities in um/s ----------
    u_original_um_s = cell(size(u_original));
    v_original_um_s = cell(size(v_original));
    for k = 1:numel(u_original)
        u_original_um_s{k} = u_masked{k} * scale_um_s;
        v_original_um_s{k} = v_masked{k} * scale_um_s;
    end

    % ---------- node coordinates in um ----------
    x_um = cell(size(x));
    y_um = cell(size(y));
    for k = 1:numel(x)
        x_um{k} = x{k} * um_per_px;
        y_um{k} = y{k} * um_per_px;
    end

    % ---------- temporal averaging ----------
    % Single-frame PIV of a slow flow is dominated by correlation noise, so
    % each frame is replaced by the mean over [n-K_AVG, n+K_AVG]. The first
    % and last K_AVG frames have no complete window and are left empty rather
    % than averaged over a truncated one.
    Nframes    = numel(u_original_um_s);
    u_avg_um_s = cell(size(u_original_um_s));
    v_avg_um_s = cell(size(v_original_um_s));

    for n = 1:Nframes
        if n > K_AVG && n <= Nframes - K_AVG
            u_stack = cat(3, u_original_um_s{n-K_AVG : n+K_AVG});
            v_stack = cat(3, v_original_um_s{n-K_AVG : n+K_AVG});
            u_avg_um_s{n} = mean(u_stack, 3, 'omitnan');
            v_avg_um_s{n} = mean(v_stack, 3, 'omitnan');
        else
            u_avg_um_s{n} = [];
            v_avg_um_s{n} = [];
        end
    end

    % ---------- grid spacing in um, needed for the gradients ----------
    dx = mean(diff(x_um{1}(1,:)));
    dy = mean(diff(y_um{1}(:,1)));

    % ---------- derived fields, instantaneous ----------
    vorticity_from_u_v          = cell(size(u_original_um_s));
    divergence_from_u_v         = cell(size(u_original_um_s));
    velocity_magnitude_from_u_v = cell(size(u_original_um_s));

    for k = 1:Nframes
        u = u_original_um_s{k};
        v = v_original_um_s{k};
        % gradient(F, hx, hy): hx is the spacing along columns, hy along rows
        [du_dx, du_dy] = gradient(u, dx, dy);
        [dv_dx, dv_dy] = gradient(v, dx, dy);
        vorticity_from_u_v{k}          = dv_dx - du_dy;   % 1/s
        divergence_from_u_v{k}         = du_dx + dv_dy;   % 1/s
        velocity_magnitude_from_u_v{k} = hypot(u, v);     % um/s
    end

    % ---------- derived fields, temporally averaged ----------
    vorticity_from_u_v_avg          = cell(size(u_avg_um_s));
    divergence_from_u_v_avg         = cell(size(u_avg_um_s));
    velocity_magnitude_from_u_v_avg = cell(size(u_avg_um_s));

    for k = 1:Nframes
        if ~isempty(u_avg_um_s{k})
            u = u_avg_um_s{k};
            v = v_avg_um_s{k};
            [du_dx, du_dy] = gradient(u, dx, dy);
            [dv_dx, dv_dy] = gradient(v, dx, dy);
            vorticity_from_u_v_avg{k}          = dv_dx - du_dy;
            divergence_from_u_v_avg{k}         = du_dx + dv_dy;
            velocity_magnitude_from_u_v_avg{k} = hypot(u, v);
        else
            vorticity_from_u_v_avg{k}          = [];
            divergence_from_u_v_avg{k}         = [];
            velocity_magnitude_from_u_v_avg{k} = [];
        end
    end

    % ---------- save everything the later parts need ----------
    save(fullfile(resultsFolder, ['piv_fields_' dropletName '.mat']), ...
         'u_original_um_s','v_original_um_s','u_avg_um_s','v_avg_um_s', ...
         'x_um','y_um','dx','dy','fps','um_per_px','scale_um_s','Nframes', ...
         'vorticity_from_u_v','divergence_from_u_v','velocity_magnitude_from_u_v', ...
         'vorticity_from_u_v_avg','divergence_from_u_v_avg', ...
         'velocity_magnitude_from_u_v_avg', '-v7.3');
    fprintf('PIV fields computed and saved.\n');
end

%% ========================================================================
%  PART B : ONE-DIMENSIONAL TEMPORAL STATISTICS
%  ========================================================================
%  Each field is a map over space at every frame. What the thesis compares
%  against the crystal growth curve is a single number per frame, so each
%  field is collapsed to its spatial mean, spatial maximum and spatial
%  standard deviation, and plotted against real time.
%
%  The mean measures the bulk flow; the maximum captures the localised jets
%  near the crystal, which the mean dilutes over the whole droplet.

if DO_PIV_1DPLOTS

    % ---------- load the fields if part A did not run this pass ----------
    if ~exist('velocity_magnitude_from_u_v_avg','var')
        load(fullfile(resultsFolder, ['piv_fields_' dropletName '.mat']), ...
             'velocity_magnitude_from_u_v_avg','vorticity_from_u_v_avg', ...
             'divergence_from_u_v_avg','fps');
        fprintf('PIV fields loaded from disk.\n');
    end

    quantities = {'velocity', 'vorticity', 'divergence'};

    for q = 1:numel(quantities)
        quantity = quantities{q};

        switch lower(quantity)
            case 'velocity'
                data_cell = velocity_magnitude_from_u_v_avg;
                y_label   = '$|\mathbf{v}|$ ($\mu$m/s)';
                title_str = 'Velocity magnitude (temporal average)';
            case 'vorticity'
                data_cell = vorticity_from_u_v_avg;
                y_label   = '$\omega$ (1/s)';
                title_str = 'Vorticity (temporal average)';
            case 'divergence'
                data_cell = divergence_from_u_v_avg;
                y_label   = '$\nabla\!\cdot\!\mathbf{v}$ (1/s)';
                title_str = 'Divergence (temporal average)';
        end

        Nf = numel(data_cell);

        % ---------- real time axis, anchored at t_start ----------
        %  PIV starts at t_start (frame 1), never at 0 unless t_start is 0.
        %  If the PIVlab export has more frames than the axis, it is extended
        %  at 1/fps while keeping the anchor.
        if exist('frameTimes','var') && ~isempty(frameTimes)
            if numel(frameTimes) >= Nf
                t = frameTimes(1:Nf);
            else
                t = frameTimes(1) + (0:Nf-1) / fps;
            end
        elseif exist('t_start','var') && ~isnan(t_start)
            t = t_start + (0:Nf-1) / fps;
        else
            t = (0:Nf-1) / fps;
        end
        t = t(:).';

        % ---------- spatial statistics per frame ----------
        % Non-finite nodes are dropped rather than propagated: they are the
        % masked crystal interior and the PIVlab outlier rejections.
        Q_mean = nan(Nf,1);
        Q_max  = nan(Nf,1);
        Q_std  = nan(Nf,1);

        for k = 1:Nf
            Q = data_cell{k};
            if isempty(Q), continue; end
            Q = Q(isfinite(Q));
            if isempty(Q), continue; end
            Q_mean(k) = mean(Q);
            Q_max(k)  = max(Q);
            Q_std(k)  = std(Q);
        end

        % ---------- smoothing, with the window shared across the project ----
        %  'omitnan' filled the empty edges by copying from the neighbours
        %  (the first and last K_AVG frames have no temporal average), which
        %  invented data. The NaNs are restored afterwards.
        validQ = ~isnan(Q_mean);
        Q_mean = smoothdata(Q_mean, 'movmean', SMOOTH_WINDOW, 'omitnan');
        Q_max  = smoothdata(Q_max,  'movmean', SMOOTH_WINDOW, 'omitnan');
        Q_std  = smoothdata(Q_std,  'movmean', SMOOTH_WINDOW, 'omitnan');
        Q_mean(~validQ) = NaN;  Q_max(~validQ) = NaN;  Q_std(~validQ) = NaN;

        kv = find(validQ);
        if ~isempty(kv)
            fprintf(['%-10s: valid frames %d-%d  ->  t = %.2f to %.2f s ' ...
                     '(SMOOTH_WINDOW=%d)\n'], quantity, kv(1), kv(end), ...
                     t(kv(1)), t(kv(end)), SMOOTH_WINDOW);
        end

        matName = fullfile(resultsFolder, ...
                  sprintf('OneDPlots_%s_%s.mat', quantity, dropletName));
        save(matName, 'Q_mean', 'Q_max', 'Q_std', 't');
        fprintf('1D data saved: %s\n', matName);

        % ---------- console summary ----------
        fprintf('\n===== %s GLOBAL STATS =====\n', upper(quantity));
        fprintf('Mean over all frames : %.4f\n', mean(Q_mean,'omitnan'));
        fprintf('Max  over all frames : %.4f\n', max(Q_max,[],'omitnan'));
        fprintf('Mean std over frames : %.4f\n', mean(Q_std,'omitnan'));

        % ---------- figure: three stacked panels sharing the x axis --------
        figure('Color','w','Position',[200 200 900 650])

        ax1 = subplot(3,1,1);
        plot(t, Q_mean, '-', 'LineWidth',1.6, 'Color',[0.20 0.40 0.80]);
        ylabel(['Mean ', y_label]);
        title(['Temporal evolution of ', title_str]);

        ax2 = subplot(3,1,2);
        plot(t, Q_max, '-', 'LineWidth',1.6, 'Color',[0.85 0.40 0.20]);
        ylabel(['Max ', y_label]);

        ax3 = subplot(3,1,3);
        plot(t, Q_std, '-', 'LineWidth',1.6, 'Color',[0.20 0.60 0.35]);
        xlabel('Time (s)'); ylabel(['Std ', y_label]);

        linkaxes([ax1 ax2 ax3],'x')

        % Nucleation markers in all three panels, labelled only in the top one
        if exist('nucTimes','var') && ~isempty(nucTimes)
            for axk = [ax1 ax2 ax3]
                nucleationLines(axk, nucTimes, [], axk == ax1);
            end
        end
        figureStyle([ax1 ax2 ax3], t);   % only ax3 carries tick labels

        figName = fullfile(resultsFolder, ...
                  ['OneDPlots_' quantity '_' dropletName '.png']);
        if exist(figName, 'file'), delete(figName); end
        exportgraphics(gcf, figName, 'Resolution', 300);
        fprintf('Saved: %s\n', figName);
    end
end

%% ========================================================================
%  PART C : FIELD VIDEOS OVER THE RAW PHOTOGRAPHS
%  ========================================================================
%  Overlays a PIV field on the original images (FotosS<n>), restricted to the
%  frames between frame_min and frame_max. These are diagnostic and
%  presentation outputs; none of the reported numbers come from them.

needPIV = DO_VIDEO_VELOCITY || DO_VIDEO_DIVERGENCE || ...
          DO_VIDEO_VORTICITY || DO_VIDEO_MAGNITUDE;

if needPIV && ~exist('u_avg_um_s','var')
    load(fullfile(resultsFolder, ['piv_fields_' dropletName '.mat']), ...
         'x_um','y_um','u_avg_um_s','v_avg_um_s', ...
         'divergence_from_u_v_avg','vorticity_from_u_v_avg', ...
         'velocity_magnitude_from_u_v_avg','fps','um_per_px');
    fprintf('PIV fields loaded from disk.\n');
end

% ---------- list the raw photographs, in numerical order ----------
if needPIV
    imgFiles = dir(fullfile(photoDir, '*.tif'));
    if isempty(imgFiles), imgFiles = dir(fullfile(photoDir, '*.png')); end
    nmI = {imgFiles.name};
    numI = nan(size(nmI));
    for i = 1:numel(nmI)
        tk = regexp(nmI{i}, '\d+', 'match');
        if ~isempty(tk), numI(i) = str2double(tk{1}); end
    end
    [~, ixI] = sort(numI);
    imgFiles = imgFiles(ixI);

    fmin = 1; fmax = numel(imgFiles);
    if exist('frame_min','var') && ~isnan(frame_min), fmin = frame_min; end
    if exist('frame_max','var') && ~isnan(frame_max)
        fmax = min(frame_max, numel(imgFiles));
    end
    imgFiles = imgFiles(fmin:fmax);
    fprintf('Videos: %d photographs (frames %d-%d).\n', numel(imgFiles), fmin, fmax);
end

if DO_VIDEO_VELOCITY
    makeFieldVideo('velocity', photoDir, imgFiles, x_um, y_um, ...
        u_avg_um_s, v_avg_um_s, [], um_per_px, fps, resultsFolder, dropletName);
end

if DO_VIDEO_DIVERGENCE
    makeFieldVideo('divergence', photoDir, imgFiles, x_um, y_um, ...
        u_avg_um_s, v_avg_um_s, divergence_from_u_v_avg, um_per_px, fps, ...
        resultsFolder, dropletName);
end

if DO_VIDEO_VORTICITY
    makeFieldVideo('vorticity', photoDir, imgFiles, x_um, y_um, ...
        u_avg_um_s, v_avg_um_s, vorticity_from_u_v_avg, um_per_px, fps, ...
        resultsFolder, dropletName);
end

if DO_VIDEO_MAGNITUDE
    makeFieldVideo('magnitude', photoDir, imgFiles, x_um, y_um, ...
        u_avg_um_s, v_avg_um_s, velocity_magnitude_from_u_v_avg, um_per_px, fps, ...
        resultsFolder, dropletName);
end


%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================

function makeFieldVideo(kind, photoDir, imgFiles, x_um, y_um, ...
                        u_avg, v_avg, scalarField, um_per_px, fps, ...
                        resultsFolder, dropletName)
%MAKEFIELDVIDEO  Render one PIV field over the photographs as an MP4.
%
%   'velocity'  draws only the vectors on the photograph.
%   anything else draws a colour map of the scalar field (one colour per PIV
%   cell) with the velocity vectors on top.

switch kind
    case 'velocity',   base = 'velocity_field';   ttl = 'Velocity field';     cbl = '';
    case 'divergence', base = 'divergence_field'; ttl = 'Divergence';         cbl = '1/s';
    case 'vorticity',  base = 'vorticity_field';  ttl = 'Vorticity';          cbl = '1/s';
    case 'magnitude',  base = 'magnitude_field';  ttl = 'Velocity magnitude'; cbl = '$\mu$m/s';
end

videoName = fullfile(resultsFolder, [base '_' dropletName '.mp4']);
if exist(videoName, 'file'), delete(videoName); end
vw = VideoWriter(videoName, 'MPEG-4');
vw.FrameRate = fps; vw.Quality = 100;
open(vw);

step        = 1;   % vector decimation; 1 draws every node
scaleVisual = 3;   % arrow length multiplier, cosmetic only

% ---------- fix a common colour scale across all frames ----------
% Percentiles rather than min/max, so a single outlier cell does not flatten
% the whole colour range for the entire video.
climVal = [];
if ~strcmp(kind,'velocity')
    allVals = [];
    for k = 1:numel(scalarField)
        S = scalarField{k};
        if ~isempty(S), allVals = [allVals; S(isfinite(S))]; end %#ok<AGROW>
    end
    if ~isempty(allVals)
        if strcmp(kind,'magnitude')
            climVal = [0, prctile(allVals, 98)];
        else
            m = prctile(abs(allVals), 98);   % symmetric about zero
            climVal = [-m, m];
        end
    end
end

fig = figure('Color','w','Position',[100 100 1200 900]);
nF = min(numel(imgFiles), numel(u_avg));

for k = 1:nF
    img = double(imread(fullfile(photoDir, imgFiles(k).name)));
    [Ny, Nx] = size(img);
    x_img_um = (0:Nx-1) * um_per_px;
    y_img_um = (0:Ny-1) * um_per_px;

    clf(fig)
    ax = axes(fig); hold(ax,'on'); ax.Clipping = 'on';

    if strcmp(kind,'velocity')
        imagesc(ax, x_img_um, y_img_um, img);
        colormap(ax, gray);
    else
        imagesc(ax, x_img_um, y_img_um, img);
        colormap(ax, gray);
        freezeColors_local(ax);   % lock the greyscale before applying parula

        S = scalarField{k};
        if ~isempty(S)
            X = x_um{k}; Y = y_um{k};
            hS = imagesc(ax, X(1,:), Y(:,1), S);
            set(hS, 'AlphaData', 0.6 * ~isnan(S));  % translucent, NaN invisible
            colormap(ax, parula);
            if ~isempty(climVal), caxis(ax, climVal); end
            cb = colorbar(ax); cb.Label.String = cbl;
        end
    end

    set(ax,'YDir','normal');
    xlim(ax,[min(x_img_um) max(x_img_um)]); ylim(ax,[min(y_img_um) max(y_img_um)]);
    axis(ax,'equal'); set(ax,'XLimMode','manual','YLimMode','manual');

    % ---------- velocity vectors, always on top ----------
    U = u_avg{k}; V = v_avg{k};
    if ~isempty(U)
        X = x_um{k}; Y = y_um{k};
        quiver(ax, X(1:step:end,1:step:end), Y(1:step:end,1:step:end), ...
               scaleVisual*U(1:step:end,1:step:end), ...
               scaleVisual*V(1:step:end,1:step:end), ...
               0, 'Color', [0 0 0], 'LineWidth', 0.8, 'MaxHeadSize', 2, ...
               'Clipping','on');
    end

    xlabel(ax,'$x$ ($\mu$m)'); ylabel(ax,'$y$ ($\mu$m)');
    title(ax, sprintf('%s - frame %d', ttl, k)); set(ax,'FontSize',14)
    drawnow
    writeVideo(vw, getframe(fig));
end

close(vw);
fprintf('Saved video: %s\n', videoName);
end


function freezeColors_local(ax)
%FREEZECOLORS_LOCAL  Lock the greyscale of the background photograph.
%
%   An axes has a single colormap. Applying parula to the scalar field would
%   otherwise repaint the greyscale photograph underneath it. Converting the
%   background image to explicit RGB makes it immune to later colormap
%   changes, which is enough for this use.

im = findobj(ax,'Type','image');
if isempty(im), return; end
im = im(end);                     % the background, i.e. the first one drawn
cdata = im.CData;
cmin = min(cdata(:)); cmax = max(cdata(:));
if cmax > cmin
    idx = round((cdata - cmin)/(cmax - cmin) * 255) + 1;
else
    idx = ones(size(cdata));
end
im.CData = ind2rgb(idx, gray(256));
end
