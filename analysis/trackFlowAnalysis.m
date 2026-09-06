%% ========================================================================
%  TRACKFLOWANALYSIS  Lagrangian flow measurement from TrackMate tracks.
%  ========================================================================
%
%  WHAT THIS SCRIPT DOES
%    Reads the trajectories exported by TrackMate (Fiji), cleans them, and
%    turns them into velocities. Where the PIV analysis measures the flow on
%    a fixed grid, this one follows the individual tracer particles, so the
%    two are complementary: PIV averages over an interrogation window and
%    smooths out the fastest excursions, while tracking retains them.
%
%    It runs in four parts, each behind its own flag in runDropletPipeline:
%
%      DO_TRACKMATE_READ      read the CSV, filter tracks, remove the
%                             stopwatch region and the spots inside crystals
%      DO_TRACKMATE_VELOCITY  per-point velocities by centred differences
%      DO_TRACKMATE_MEANVEL   mean speed per frame against real time
%      DO_TRACKMATE_VELPLOT   tracks coloured by speed (figure or video)
%
%    This is a SCRIPT, not a function: it is called with run() from
%    runDropletPipeline and shares its workspace.
%
%  THE TWO CLEANING STEPS, AND WHY THEY MATTER
%    1. STOPWATCH. The recordings carry a burnt-in timer. Its digits change
%       every frame, so TrackMate detects them as spots and links them into
%       spurious tracks with large apparent velocities. The timer region is
%       drawn once per droplet and every spot inside it is discarded.
%    2. CRYSTAL INTERIOR. The crystal surface has texture, which is detected
%       as spots too. Those are not tracer particles, and leaving them in
%       would let the crystal's own growth register as fluid motion, which is
%       precisely the effect this thesis is trying to measure. This filter
%       must stay on.
%
%  INPUTS FROM THE PIPELINE
%    fps, um_per_px, dropletName, resultsFolder, frameTimes, nucTimes,
%    colorTM, SMOOTH_WINDOW, TM_* parameters, flipTrackMateY,
%    tm_remove_timer, tm_timer_box, tm_redraw_box, and 'masks' when
%    TM_REMOVE_INSIDE_MASK is true.
%
%  OUTPUT FILES  (in resultsFolder)
%    trackmate_reduced_<tag>_<droplet>.mat     cleaned tracks
%    trackmate_velocities_<tag>_<droplet>.mat  per-point velocity table
%    trackmate_meanspeed_<tag>_<droplet>.mat   mean speed per frame
%    matching .png figures, and .mp4 when the video flags are on
%
%    <tag> is 'masked' when the crystal-interior filter was applied, so that
%    filtered and unfiltered results can never be confused on disk.
%
%  See also RUNDROPLETPIPELINE, PIVFLOWANALYSIS, FIGURESTYLE.

%% ========================================================================
%  LOCATE THE INPUT FILES
%  ========================================================================

% ---------- the TrackMate spot export (*allspots*.csv) ----------
sessionDir = fileparts(resultsFolder);
csvHit = dir(fullfile(sessionDir, '*allspots*.csv'));
assert(~isempty(csvHit), 'No *allspots*.csv found in %s', sessionDir);
trackmateCsvPath = fullfile(sessionDir, csvHit(1).name);
fprintf('TrackMate CSV: %s\n', trackmateCsvPath);

% ---------- the image folder TrackMate ran on (FPICTrackmate) ----------
dd = dir(sessionDir); dd = dd([dd.isdir]);
trackmateImageFolder = '';
for i = 1:numel(dd)
    if startsWith(dd(i).name, 'FPICTrack', 'IgnoreCase', true) || ...
       (contains(dd(i).name,'FPIC','IgnoreCase',true) && ...
        contains(dd(i).name,'rack','IgnoreCase',true))
        trackmateImageFolder = fullfile(sessionDir, dd(i).name); break
    end
end
assert(~isempty(trackmateImageFolder), ...
       'No FPICTrackmate folder found in %s', sessionDir);
fprintf('TrackMate images: %s\n', trackmateImageFolder);

%% ========================================================================
%  PART A : READ AND CLEAN THE TRACKS
%  ========================================================================
if DO_TRACKMATE_READ

    % ---------- 1. read the CSV ----------
    T = readtable(trackmateCsvPath);
    T(1,:) = [];      % TrackMate writes a units row under the header
    T = T(:,2:9);     % keep the identity, position, frame and quality columns

    frame_r   = double(T.FRAME);
    x_r       = double(T.POSITION_X);
    y_r       = double(T.POSITION_Y);
    trackID_r = double(T.TRACK_ID);
    quality_r = double(T.QUALITY);

    % Spots that were detected but never linked into a track come out with
    % NaN identifiers and are of no use here.
    ok = ~isnan(frame_r) & ~isnan(x_r) & ~isnan(y_r) & ...
         ~isnan(trackID_r) & ~isnan(quality_r);
    frame_r = frame_r(ok); x_r = x_r(ok); y_r = y_r(ok);
    trackID_r = trackID_r(ok); quality_r = quality_r(ok);

    fprintf('\n=== TRACKMATE ===\nValid spots before filtering: %d\n', numel(frame_r));

    % ---------- 2. group the spots into trajectories ----------
    uniqueTracks = unique(trackID_r);
    tracks_tmp = struct();
    for i = 1:numel(uniqueTracks)
        idx = (trackID_r == uniqueTracks(i));
        [fs, order] = sort(frame_r(idx));
        xi = x_r(idx); yi = y_r(idx); qi = quality_r(idx);
        tracks_tmp(i).id    = uniqueTracks(i);
        tracks_tmp(i).frame = fs;
        tracks_tmp(i).x     = xi(order);
        tracks_tmp(i).y     = yi(order);
        tracks_tmp(i).q     = qi(order);
    end

    % ---------- 3. filter by track length and detection quality ----------
    % Short tracks cannot give a centred velocity difference and are usually
    % spurious links between unrelated detections.
    lens = arrayfun(@(s) numel(s.frame), tracks_tmp);
    minQ = arrayfun(@(s) min(s.q), tracks_tmp);
    keep = find(lens >= TM_MIN_LENGTH & minQ >= TM_MIN_QUALITY);

    tracks_reduced = tracks_tmp(keep);
    validIDs       = uniqueTracks(keep);

    idxSpots   = ismember(trackID_r, validIDs);
    frame_rf   = frame_r(idxSpots);
    x_rf       = x_r(idxSpots);
    y_rf       = y_r(idxSpots);
    trackID_rf = trackID_r(idxSpots);

    % ---------- 3a. remove the spots on the burnt-in stopwatch ----------
    if tm_remove_timer
        boxFile  = fullfile(resultsFolder, ['timer_box_' dropletName '.mat']);
        timerBox = [];
        if ~isempty(tm_timer_box)
            timerBox = tm_timer_box;
        elseif isfile(boxFile) && ~tm_redraw_box
            S = load(boxFile, 'timerBox'); timerBox = S.timerBox;
            fprintf('Timer region reused: [%.0f %.0f %.0f %.0f]\n', timerBox);
        end

        if isempty(timerBox) || tm_redraw_box
            it = dir(fullfile(trackmateImageFolder, '*.tif'));
            if isempty(it), it = dir(fullfile(trackmateImageFolder, '*.png')); end
            nmt = {it.name}; nt = nan(size(nmt));
            for ii = 1:numel(nmt)
                tk = regexp(nmt{ii}, '\d+', 'match');
                if ~isempty(tk), nt(ii) = str2double(tk{1}); end
            end
            if all(~isnan(nt)), [~,ix] = sort(nt); else, [~,ix] = sort(nmt); end
            it = it(ix);
            img1 = imread(fullfile(trackmateImageFolder, it(1).name));
            if ndims(img1) == 3, img1 = im2gray(img1); end
            hfB = figure('Name','Draw over the stopwatch','Color','w');
            imshow(imadjust(img1));
            title('Enclose the stopwatch digits. Double-click when done.');
            rec = drawrectangle; wait(rec);
            pos = rec.Position;
            timerBox = [pos(1), pos(2), pos(1)+pos(3), pos(2)+pos(4)];
            close(hfB);
            save(boxFile, 'timerBox');
            fprintf('Timer region saved: [%.0f %.0f %.0f %.0f]\n', timerBox);
        end

        keepSpot = ~(x_rf >= timerBox(1) & x_rf <= timerBox(3) & ...
                     y_rf >= timerBox(2) & y_rf <= timerBox(4));
        nRem = sum(~keepSpot);
        frame_rf   = frame_rf(keepSpot);
        x_rf       = x_rf(keepSpot);
        y_rf       = y_rf(keepSpot);
        trackID_rf = trackID_rf(keepSpot);
        fprintf('Spots removed by the stopwatch region: %d\n', nRem);

        % Rebuild the trajectories from the surviving spots.
        for t = 1:numel(tracks_reduced)
            id = tracks_reduced(t).id;
            sel = (trackID_rf == id);
            fsel = frame_rf(sel); xsel = x_rf(sel); ysel = y_rf(sel);
            [fs, o] = sort(fsel);
            tracks_reduced(t).frame = fs;
            tracks_reduced(t).x = xsel(o);
            tracks_reduced(t).y = ysel(o);
        end
        emptyTr = arrayfun(@(s) isempty(s.frame), tracks_reduced);
        tracks_reduced(emptyTr) = [];
    end

    % ---------- 3b. remove the spots falling inside a crystal ----------
    maskedTag = '';
    if TM_REMOVE_INSIDE_MASK
        if ~exist('masks','var')
            load(fullfile(resultsFolder, ['masks_' dropletName '.mat']), 'masks');
            fprintf('masks loaded from disk for spot filtering.\n');
        end

        keepSpot = true(size(frame_rf));
        for s = 1:numel(frame_rf)
            fr = frame_rf(s) + 1;             % TrackMate frame 0 -> index 1
            if fr < 1 || fr > numel(masks), continue; end
            M = masks{fr};
            if isempty(M), continue; end
            col = round(x_rf(s)); row = round(y_rf(s));
            if flipTrackMateY, row = size(M,1) - row + 1; end
            if row >= 1 && row <= size(M,1) && col >= 1 && col <= size(M,2)
                if M(row,col), keepSpot(s) = false; end
            end
        end

        nRemoved = sum(~keepSpot);
        frame_rf   = frame_rf(keepSpot);
        x_rf       = x_rf(keepSpot);
        y_rf       = y_rf(keepSpot);
        trackID_rf = trackID_rf(keepSpot);
        fprintf('Spots removed for falling inside crystals: %d\n', nRemoved);

        for t = 1:numel(tracks_reduced)
            id = tracks_reduced(t).id;
            sel = (trackID_rf == id);
            fsel = frame_rf(sel); xsel = x_rf(sel); ysel = y_rf(sel);
            [fs, o] = sort(fsel);
            tracks_reduced(t).frame = fs;
            tracks_reduced(t).x = xsel(o);
            tracks_reduced(t).y = ysel(o);
        end
        emptyTr = arrayfun(@(s) isempty(s.frame), tracks_reduced);
        tracks_reduced(emptyTr) = [];

        maskedTag = 'masked';   % goes into every output filename
    end

    nTracks_reduced = numel(tracks_reduced);
    fprintf('Trajectories after filtering: %d\n', nTracks_reduced);
    assert(nTracks_reduced > 0, ...
           'No trajectory survived. Check the filters or the mask.');

    % ---------- 4. list the images in numerical order ----------
    imgF = dir(fullfile(trackmateImageFolder, '*.tif'));
    if isempty(imgF), imgF = dir(fullfile(trackmateImageFolder, '*.png')); end
    assert(~isempty(imgF), 'No images in %s', trackmateImageFolder);
    nmT = {imgF.name}; numT = nan(size(nmT));
    for i = 1:numel(nmT)
        tk = regexp(nmT{i}, '\d+', 'match');
        if ~isempty(tk), numT(i) = str2double(tk{1}); end
    end
    if all(~isnan(numT)), [~,ixT] = sort(numT); else, [~,ixT] = sort(nmT); end
    imgF = imgF(ixT);
    nImg = numel(imgF);
    nRender = min(nImg, max(frame_rf) + 1);

    % ---------- 5. optional video of the cleaned tracks ----------
    if TM_MAKE_VIDEO
        colorsT = hsv(nTracks_reduced);
        vName = fullfile(resultsFolder, ['trackmate_' maskedTag '_' dropletName '.mp4']);
        vName = strrep(vName, '__', '_');

        if exist(vName,'file'), delete(vName); end
        vw = VideoWriter(vName,'MPEG-4'); vw.FrameRate = fps; vw.Quality = 100;
        open(vw);

        figR = figure('Color','w','Position',[100 100 1280 960]);
        axR = axes(figR); axR.Clipping = 'on';

        for k = 0:(nRender-1)
            img = imread(fullfile(trackmateImageFolder, imgF(k+1).name));
            if ndims(img) == 3, img = im2gray(img); end
            [Ny,Nx] = size(img);
            cla(axR); imagesc(axR, img); colormap(axR, gray);
            set(axR,'YDir','normal'); hold(axR,'on');
            axis(axR,'equal'); xlim(axR,[0 Nx]); ylim(axR,[0 Ny]);
            set(axR,'XLimMode','manual','YLimMode','manual');

            for t = 1:nTracks_reduced
                tf = tracks_reduced(t).frame;
                valid = (tf <= k) & (tf >= k - TM_TAIL);
                if sum(valid) < TM_MIN_POINTS, continue; end
                xs = tracks_reduced(t).x(valid); ys = tracks_reduced(t).y(valid);
                if numel(xs) > 4
                    xs = smoothdata(xs,'movmean',TM_SMOOTH_WIN);
                    ys = smoothdata(ys,'movmean',TM_SMOOTH_WIN);
                end
                fsub  = tf(valid);
                age   = (k - fsub(end)) / TM_TAIL;   % fade the older tails
                alpha = max(0.2, 1 - age);
                plot(axR, xs, ys, 'Color',[colorsT(t,:) alpha], 'LineWidth',1.6);
            end

            idxk = (frame_rf == k);
            if any(idxk), scatter(axR, x_rf(idxk), y_rf(idxk), 14, 'r', 'filled'); end

            title(axR, sprintf('TrackMate - frame %d', k), ...
                  'FontSize',16,'FontWeight','bold');
            drawnow limitrate; writeVideo(vw, getframe(figR));
            if mod(k,50) == 0, fprintf('  video frame %d/%d\n', k+1, nRender); end
        end
        close(vw);
        fprintf('Video saved: %s\n', vName);
    end

    % ---------- 6. summary figure with TM_NSHOW representative frames ------
    framesShow = round(linspace(0, nRender-1, TM_NSHOW));
    figS = figure('Color','w','Position',[100 100 1400 1050]);
    colorsT = hsv(nTracks_reduced);
    for s = 1:TM_NSHOW
        k = framesShow(s);
        img = imread(fullfile(trackmateImageFolder, imgF(k+1).name));
        if ndims(img) == 3, img = im2gray(img); end
        [Ny,Nx] = size(img);
        subplot(2,2,s); imagesc(img); colormap(gray); set(gca,'YDir','normal'); hold on
        axis equal; xlim([0 Nx]); ylim([0 Ny]);
        for t = 1:nTracks_reduced
            tf = tracks_reduced(t).frame;
            valid = (tf <= k) & (tf >= k - TM_TAIL);
            if sum(valid) < TM_MIN_POINTS, continue; end
            xs = tracks_reduced(t).x(valid); ys = tracks_reduced(t).y(valid);
            if numel(xs) > 4
                xs = smoothdata(xs,'movmean',TM_SMOOTH_WIN);
                ys = smoothdata(ys,'movmean',TM_SMOOTH_WIN);
            end
            plot(xs, ys, 'Color', colorsT(t,:), 'LineWidth',1.5);
        end
        idxk = (frame_rf == k);
        if any(idxk), scatter(x_rf(idxk), y_rf(idxk), 1.5, 'r', 'filled'); end
        title(sprintf('Frame %d', k), 'FontSize',14,'FontWeight','bold');
    end
    sgtitle('TrackMate reduced - representative frames', ...
            'FontSize',16,'FontWeight','bold');
    sName = fullfile(resultsFolder, ['trackmate_' maskedTag '_' dropletName '.png']);
    sName = strrep(sName, '__', '_');
    if exist(sName,'file'), delete(sName); end
    exportgraphics(figS, sName, 'Resolution', 300);
    fprintf('Summary saved: %s\n', sName);

    % ---------- 7. save for the velocity stage ----------
    outName = ['trackmate_reduced_' maskedTag '_' dropletName '.mat'];
    outName = strrep(outName, '__', '_');   % avoid a double _ if the tag is empty
    save(fullfile(resultsFolder, outName), ...
         'tracks_reduced','nTracks_reduced','frame_rf','x_rf','y_rf','maskedTag');
    fprintf('TrackMate read and filtered. Saved: %s\n', outName);
end

%% ========================================================================
%  PART B : PER-POINT VELOCITIES
%  ========================================================================
%  Velocity at each point of each track, from a centred difference spanning
%  TM_VELOCITY_N points either side. A centred difference over several points
%  rather than between consecutive frames is what makes these velocities
%  usable: the particle displacement per frame is comparable to the
%  localisation uncertainty, so a one-frame difference is mostly noise.
%
%  Points closer than TM_VELOCITY_N to either end of a track have no
%  symmetric window and are left as NaN rather than estimated one-sidedly.

if DO_TRACKMATE_VELOCITY

    % Load if part A did not run this pass, preferring the masked version.
    if ~exist('tracks_reduced','var')
        fM = fullfile(resultsFolder, ['trackmate_reduced_masked_' dropletName '.mat']);
        fP = fullfile(resultsFolder, ['trackmate_reduced_' dropletName '.mat']);
        if isfile(fM)
            load(fM, 'tracks_reduced','nTracks_reduced','maskedTag');
            fprintf('tracks_reduced (masked) loaded from disk.\n');
        else
            load(fP, 'tracks_reduced','nTracks_reduced');
            maskedTag = '';
            fprintf('tracks_reduced (unmasked) loaded from disk.\n');
        end
    end
    if ~exist('maskedTag','var'), maskedTag = ''; end

    Nv = TM_VELOCITY_N;
    fprintf('\n=== VELOCITIES (centred difference, N=%d) ===\n', Nv);

    tracks_vel = tracks_reduced;

    for t = 1:nTracks_reduced
        f = tracks_reduced(t).frame(:);
        x = tracks_reduced(t).x(:);
        y = tracks_reduced(t).y(:);
        n = numel(f);
        vx = nan(n,1); vy = nan(n,1); v = nan(n,1);

        for i = 1:n
            ip = i - Nv; in = i + Nv;
            if ip < 1 || in > n, continue; end
            df = f(in) - f(ip);
            if df <= 0, continue; end
            % The frame gap is used rather than the index gap, because tracks
            % can have missing frames where the detection dropped out.
            dt = df / fps;
            vx(i) = (x(in) - x(ip)) * um_per_px / dt;   % um/s
            vy(i) = (y(in) - y(ip)) * um_per_px / dt;
            v(i)  = hypot(vx(i), vy(i));
        end

        tracks_vel(t).vx = vx;
        tracks_vel(t).vy = vy;
        tracks_vel(t).v  = v;
    end

    % ---------- flat table with every point of every track ----------
    aID = []; aF = []; aX = []; aY = []; aVx = []; aVy = []; aSp = [];
    for t = 1:nTracks_reduced
        f = tracks_vel(t).frame(:);
        aID = [aID; repmat(tracks_vel(t).id, numel(f), 1)]; %#ok<AGROW>
        aF  = [aF;  f];                                     %#ok<AGROW>
        aX  = [aX;  tracks_vel(t).x(:)];                    %#ok<AGROW>
        aY  = [aY;  tracks_vel(t).y(:)];                    %#ok<AGROW>
        aVx = [aVx; tracks_vel(t).vx(:)];                   %#ok<AGROW>
        aVy = [aVy; tracks_vel(t).vy(:)];                   %#ok<AGROW>
        aSp = [aSp; tracks_vel(t).v(:)];                    %#ok<AGROW>
    end

    trackVelocityTable_reduced = table(aID, aF, aX, aY, aVx, aVy, aSp, ...
        'VariableNames', {'TrackID','Frame','X_px','Y_px', ...
                          'Vx_um_s','Vy_um_s','Speed_um_s'});

    outName = ['trackmate_velocities_' maskedTag '_' dropletName '.mat'];
    outName = strrep(outName, '__', '_');
    save(fullfile(resultsFolder, outName), ...
         'tracks_vel','trackVelocityTable_reduced','maskedTag');
    fprintf('Velocities saved: %s\n', outName);
    fprintf('Velocities computed for %d trajectories.\n', nTracks_reduced);
end

%% ========================================================================
%  PART C : MEAN SPEED PER FRAME AGAINST TIME
%  ========================================================================
%  The single curve compared against the crystal growth rate: the mean speed
%  of all particles visible in each frame. The maximum is kept alongside it,
%  because the flow near a crystal is localised and the spatial mean over the
%  whole droplet dilutes it.

if DO_TRACKMATE_MEANVEL

    if ~exist('trackVelocityTable_reduced','var')
        fM = fullfile(resultsFolder, ['trackmate_velocities_masked_' dropletName '.mat']);
        fP = fullfile(resultsFolder, ['trackmate_velocities_' dropletName '.mat']);
        if isfile(fM), load(fM,'trackVelocityTable_reduced','maskedTag');
        else,          load(fP,'trackVelocityTable_reduced'); maskedTag = ''; end
    end
    if ~exist('maskedTag','var'), maskedTag = ''; end

    Tv = trackVelocityTable_reduced;

    % ---------- mean and maximum speed per frame ----------
    framesU   = unique(Tv.Frame);
    meanSpeed = nan(numel(framesU),1);
    maxSpeed  = nan(numel(framesU),1);
    for i = 1:numel(framesU)
        s = Tv.Speed_um_s(Tv.Frame == framesU(i));
        s = s(isfinite(s));
        if ~isempty(s)
            meanSpeed(i) = mean(s);
            maxSpeed(i)  = max(s);
        end
    end

    % ---------- smoothing, same window as crystals and PIV ----------
    % Frames with no data keep their NaN afterwards, so nothing is invented.
    validS    = ~isnan(meanSpeed);
    meanSpeed = smoothdata(meanSpeed, 'movmean', SMOOTH_WINDOW, 'omitnan');
    maxSpeed  = smoothdata(maxSpeed,  'movmean', SMOOTH_WINDOW, 'omitnan');
    meanSpeed(~validS) = NaN;
    maxSpeed(~validS)  = NaN;
    fprintf('TrackMate meanSpeed: SMOOTH_WINDOW = %d frames\n', SMOOTH_WINDOW);

    % ---------- real time of each frame ----------
    %  TrackMate counts frames from 0, so frame k maps to index k+1 in
    %  frameTimes. Anchored at t_start; a frame beyond the axis is extended
    %  at 1/fps rather than left as NaN.
    framesU = framesU(:);
    if exist('frameTimes','var') && ~isempty(frameTimes)
        tSpeed = frameTimes(1) + framesU / fps;
        idxT   = framesU + 1;
        good   = idxT >= 1 & idxT <= numel(frameTimes);
        tSpeed(good) = frameTimes(idxT(good));
    elseif exist('t_start','var') && ~isnan(t_start)
        tSpeed = t_start + framesU / fps;
    else
        tSpeed = framesU / fps;
    end
    tSpeed = tSpeed(:);

    % ---------- plot ----------
    figMV = figure('Name','Mean track speed vs time','Color','w', ...
                   'Position',[100 100 1000 620]);
    plot(tSpeed, meanSpeed, '-', 'LineWidth', 1.8, 'Color', colorTM);
    xlabel('Time (s)'); ylabel('Mean particle speed ($\mu$m/s)');
    title('Mean track speed per frame'); grid on

    % The axis spans the full recording, not just the frames with tracks, so
    % that this figure is directly comparable with the PIV and crystal ones.
    tRef = tSpeed;
    if exist('frameTimes','var') && ~isempty(frameTimes), tRef = frameTimes; end
    if exist('nucTimes','var') && ~isempty(nucTimes)
        nucleationLines(gca, nucTimes, [], true);
    end
    figureStyle(gca, tRef);

    sName = fullfile(resultsFolder, ...
            ['trackmate_meanspeed_' maskedTag '_' dropletName '.png']);
    sName = strrep(sName, '__','_');
    if exist(sName,'file'), delete(sName); end
    exportgraphics(figMV, sName, 'Resolution', 300);
    fprintf('Mean speed plot saved: %s\n', sName);

    mvName = ['trackmate_meanspeed_' maskedTag '_' dropletName '.mat'];
    mvName = strrep(mvName, '__', '_');
    save(fullfile(resultsFolder, mvName), ...
         'framesU','meanSpeed','maxSpeed','tSpeed','maskedTag');
    fprintf('Mean speed saved: %s\n', mvName);
end

%% ========================================================================
%  PART D : TRACKS COLOURED BY SPEED
%  ========================================================================
%  Where the mean-speed curve shows when the flow is fast, this shows where.
%  It is what makes the flow towards the crystals visible rather than merely
%  measurable.

if DO_TRACKMATE_VELPLOT

    if ~exist('trackVelocityTable_reduced','var')
        fM = fullfile(resultsFolder, ['trackmate_velocities_masked_' dropletName '.mat']);
        fP = fullfile(resultsFolder, ['trackmate_velocities_' dropletName '.mat']);
        if isfile(fM), load(fM,'trackVelocityTable_reduced','maskedTag');
        else,          load(fP,'trackVelocityTable_reduced'); maskedTag = ''; end
    end
    if ~exist('maskedTag','var'), maskedTag = ''; end

    Tv = trackVelocityTable_reduced;

    imgF = dir(fullfile(trackmateImageFolder, '*.tif'));
    if isempty(imgF), imgF = dir(fullfile(trackmateImageFolder, '*.png')); end
    nmT = {imgF.name}; numT = nan(size(nmT));
    for i = 1:numel(nmT)
        tk = regexp(nmT{i}, '\d+', 'match');
        if ~isempty(tk), numT(i) = str2double(tk{1}); end
    end
    if all(~isnan(numT)), [~,ixT] = sort(numT); else, [~,ixT] = sort(nmT); end
    imgF = imgF(ixT);
    nImg = numel(imgF);
    maxF = max(Tv.Frame);
    nRender = min(nImg, maxF + 1);

    % Colour scale from the 98th percentile, so a single fast outlier does
    % not compress every other track into the bottom of the colormap.
    if isempty(TM_VELPLOT_CLIM)
        sp = Tv.Speed_um_s(isfinite(Tv.Speed_um_s));
        climV = [0, prctile(sp, 98)];
    else
        climV = TM_VELPLOT_CLIM;
    end

    tail     = TM_VELPLOT_TAIL;
    trackIDs = unique(Tv.TrackID);

    if strcmpi(TM_VELPLOT_MODE, 'video')
        % ---------- video ----------
        vName = fullfile(resultsFolder, ...
                ['trackmate_speed_' maskedTag '_' dropletName '.mp4']);
        vName = strrep(vName, '__','_');
        if exist(vName,'file'), delete(vName); end
        vw = VideoWriter(vName,'MPEG-4'); vw.FrameRate = fps; vw.Quality = 100;
        open(vw);

        fig = figure('Color','w','Position',[100 100 1280 960]);
        for k = 0:(nRender-1)
            clf(fig); ax = axes(fig);
            drawSpeedFrameTM(ax, k, imgF, trackmateImageFolder, Tv, ...
                             trackIDs, tail, climV);
            cb = colorbar(ax); cb.Label.String = 'Speed ($\mu$m/s)';
            title(ax, sprintf('Track speed - frame %d', k), 'FontSize',15);
            drawnow; writeVideo(vw, getframe(fig));
            if mod(k,50) == 0, fprintf('  speed video %d/%d\n', k+1, nRender); end
        end
        close(vw);
        fprintf('Speed video: %s\n', vName);
    else
        % ---------- static figure with N representative frames ----------
        nShow = TM_VELPLOT_NFRAMES;
        framesShow = round(linspace(0, nRender-1, nShow));
        nc = ceil(sqrt(nShow)); nr = ceil(nShow/nc);
        fig = figure('Color','w','Position',[100 100 1500 1000]);
        for s = 1:nShow
            ax = subplot(nr, nc, s);
            drawSpeedFrameTM(ax, framesShow(s), imgF, trackmateImageFolder, ...
                             Tv, trackIDs, tail, climV);
            if exist('frameTimes','var') && framesShow(s)+1 <= numel(frameTimes)
                title(ax, sprintf('t = %.1f s', frameTimes(framesShow(s)+1)));
            else
                title(ax, sprintf('Frame %d', framesShow(s)));
            end
        end
        cb = colorbar('Position',[0.93 0.15 0.02 0.7]);
        cb.Label.String = 'Speed ($\mu$m/s)';
        colormap(parula); caxis(climV);
        sgtitle('Track speed - representative frames', ...
                'FontSize',15,'FontWeight','bold');

        sName = fullfile(resultsFolder, ...
                ['trackmate_speed_frames_' maskedTag '_' dropletName '.png']);
        sName = strrep(sName, '__','_');
        if exist(sName,'file'), delete(sName); end
        exportgraphics(fig, sName, 'Resolution', 300);
        fprintf('Speed figure: %s\n', sName);
    end
end


%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================

function drawSpeedFrameTM(ax, k, imgF, imgFolder, Tv, trackIDs, tail, climV)
%DRAWSPEEDFRAMETM  Draw one frame with the recent track tails coloured by speed.
%
%   Each tail segment is drawn separately with its own colour, because the
%   speed varies along a trajectory and a single colour per track would hide
%   exactly the acceleration this analysis is looking for.

img = imread(fullfile(imgFolder, imgF(k+1).name));
if ndims(img) == 3, img = im2gray(img); end
[Ny,Nx] = size(img);
imagesc(ax, img); colormap(ax, gray); set(ax,'YDir','normal'); hold(ax,'on');
axis(ax,'equal'); xlim(ax,[0 Nx]); ylim(ax,[0 Ny]);
freezeColors_localTM(ax);

for id = trackIDs'
    rows = Tv(Tv.TrackID == id & Tv.Frame <= k & Tv.Frame >= k-tail, :);
    if height(rows) < 2, continue; end
    xx = rows.X_px; yy = rows.Y_px; ss = rows.Speed_um_s;
    for j = 1:height(rows)-1
        cval = mean([ss(j) ss(j+1)],'omitnan');
        if isnan(cval), cval = 0; end
        col = colorFromClim(cval, climV);
        plot(ax, xx(j:j+1), yy(j:j+1), '-', 'Color', col, 'LineWidth', 2);
    end
end
colormap(ax, parula); caxis(ax, climV);
end


function col = colorFromClim(val, clim)
%COLORFROMCLIM  Map a scalar onto the parula colormap within given limits.
cmap = parula(256);
t = (val - clim(1)) / (clim(2) - clim(1));
t = min(max(t,0),1);
col = cmap(round(t*255) + 1, :);
end


function freezeColors_localTM(ax)
%FREEZECOLORS_LOCALTM  Lock the greyscale of the background image.
%   See the equivalent helper in pivFlowAnalysis for the reasoning.
im = findobj(ax,'Type','image');
if isempty(im), return; end
im = im(end);
cdata = im.CData;
cmin = min(cdata(:)); cmax = max(cdata(:));
if cmax > cmin
    idx = round((cdata-cmin)/(cmax-cmin)*255)+1;
else
    idx = ones(size(cdata));
end
im.CData = ind2rgb(idx, gray(256));
end
