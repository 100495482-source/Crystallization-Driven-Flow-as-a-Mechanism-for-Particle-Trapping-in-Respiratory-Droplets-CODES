%% ========================================================================
%  CRYSTALGROWTHANALYSIS  Crystal segmentation and growth measurement.
%  ========================================================================
%
%  WHAT THIS SCRIPT DOES
%    Everything concerning the crystals themselves: finding them in the image
%    sequence, measuring how their area and radius evolve, and converting
%    that into the growth velocities that are later compared against the
%    measured flow.
%
%    It runs in four parts, each behind its own flag in runDropletPipeline:
%
%      DO_GLOBAL_MASKS          segment the crystals in every frame
%      DO_GLOBAL_VELOCITIES     total area, global growth rate and the two
%                               definitions of front velocity
%      DO_INDIVIDUAL_MASKS      split the global mask into one mask per
%                               crystal, plus the buffer ring around each
%      DO_INDIVIDUAL_VELOCITIES per-crystal area, growth rate and radial
%                               growth velocity
%
%    This is a SCRIPT, not a function: it is called with run() from
%    runDropletPipeline and shares its workspace.
%
%  THE SEGMENTATION RUNS BACKWARDS, AND THAT IS DELIBERATE
%    Crystals are easiest to see at the end of the sequence, when they are
%    large and dense, and hardest at the beginning, when they are a few
%    faint pixels. Segmenting forwards therefore fails exactly where the
%    growth rate matters most.
%
%    Instead the last frame is used as a maximum envelope (drawn by hand, or
%    thresholded), and the sequence is walked from the last frame to the
%    first. Each frame is intersected with the already-processed frame that
%    follows it. Going backwards in time a crystal can only stay the same
%    size or shrink, so the mask is anchored to a region that is known to be
%    correct and cannot jump to an unrelated bright patch.
%
%  INPUTS FROM THE PIPELINE
%    otsuDir, resultsFolder, dropletName, fps, um_per_px, NCrys, frame_max,
%    frameTimes, nucList, crystalColors, colorCrystal, LW_PLOT,
%    INVERT_IMAGE, OTSU_SCALE, GAUSS_SIGMA, OTSU_MIN_SIZE, MORPH_OPEN,
%    MORPH_CLOSE, CONVEX_STRENGTH, USE_POLY_SMOOTH, POLY_TOL, MANUAL_LAST,
%    DISPLAY_N, SHOW_MASK_MOVIE, N_FINITE_DIFF, SMOOTH_WINDOW, bufferFactor,
%    and 'specs' so the nucleation frames can be written back.
%
%  OUTPUT FILES  (in resultsFolder)
%    masks_<droplet>.mat                global mask per frame
%    envelope_<droplet>.mat             hand-drawn maximum envelope
%    crystal_areas_<droplet>.mat        per-object areas and nucleation list
%    growth_global_<droplet>.mat        total area, growth rate, front velocities
%    individual_crystals_<droplet>.mat  per-crystal masks, centroids, buffers
%    growth_individual_<droplet>.mat    per-crystal area and growth rate
%    radial_velocity_<droplet>.mat      per-crystal and global radial velocity
%    matching .png figures
%
%  See also RUNDROPLETPIPELINE, PIVFLOWANALYSIS, TRACKFLOWANALYSIS.

%% ========================================================================
%  PART 0 : TIME AXIS CHECK  (runs unconditionally)
%  ========================================================================
%  The specifications are authoritative. Any .mat already on disk carrying a
%  stale time vector is rewritten in place; areas, masks and velocities are
%  left untouched. This is what allows t_start to be corrected without
%  re-running any analysis.

if exist('t_start','var') && ~isnan(t_start)
    t0 = t_start;                                   % from the pipeline (specs)
elseif exist('specs','var') && isfield(specs,'t_start') && ~isnan(specs.t_start)
    t0 = specs.t_start;
else
    t0 = 0;
    warning('No t_start in the specifications: the axis will start at 0 s.');
end

% Nucleation frames, also from the specifications.
if exist('specs','var') && isfield(specs,'nucleation_frames') ...
                        && ~isempty(specs.nucleation_frames) ...
                        && all(~isnan(specs.nucleation_frames))
    nucSpecs = sort(specs.nucleation_frames(:)');
elseif exist('nucFrame','var') && ~isnan(nucFrame)
    nucSpecs = nucFrame;
else
    nucSpecs = [];
end

% runDropletPipeline already performs the migration. This call only covers
% the case of running this script on its own.
if ~exist('t_start','var')
    syncTimeAxis(resultsFolder, dropletName, t0, fps, nucSpecs);
end

if ~exist('LW_PLOT','var'), LW_PLOT = 1.6; end
area_scale_um2 = um_per_px^2;

% ---------- list the frames in numerical order ----------
% Sorting by filename would put frame 10 before frame 2, so the leading
% integer in each name is extracted and sorted numerically.
f = dir(fullfile(otsuDir, '*.tif'));
if isempty(f), f = dir(fullfile(otsuDir, '*.png')); end
assert(~isempty(f), 'No images in %s', otsuDir);

nm  = {f.name};
num = nan(size(nm));
for i = 1:numel(nm)
    t = regexp(nm{i}, '\d+', 'match');
    if ~isempty(t), num(i) = str2double(t{1}); end
end
if all(~isnan(num)), [~,ix] = sort(num); else, [~,ix] = sort(nm); end
names = nm(ix);

totalFrames = numel(names);
if ~isnan(frame_max) && frame_max <= totalFrames
    totalFrames = frame_max;
end
names = names(1:totalFrames);

% The time axis is always rebuilt here and never inherited from a .mat file.
frameTimes = t0 + (0:totalFrames-1) / fps;

%% ========================================================================
%  PART 1 : GLOBAL CRYSTAL MASKS
%  ========================================================================
if DO_GLOBAL_MASKS

    % ---------- reuse existing masks if the user wants to ----------
    % Segmentation involves manual drawing, so it is worth not repeating.
    masksFile = fullfile(resultsFolder, ['masks_' dropletName '.mat']);
    SKIP_ALL = false;

    if isfile(masksFile)
        resp = questdlg(sprintf(['Masks already exist for %s.\n' ...
                        'Keep the existing masks?'], dropletName), ...
                        'Existing masks', ...
                        'Yes, keep them', 'No, recompute', 'Yes, keep them');
        if strcmp(resp, 'Yes, keep them')
            fprintf('Loading existing masks from %s\n', masksFile);
            load(masksFile, 'masks');
            af = fullfile(resultsFolder, ['crystal_areas_' dropletName '.mat']);
            if isfile(af)
                % Deliberately NOT loading frameTimes: the axis comes from
                % the specifications, never from a stored result.
                Sld = load(af, 'areasPx','areasUm','nucList');
                areasPx = Sld.areasPx; areasUm = Sld.areasUm; nucList = Sld.nucList;
                clear Sld
            end
            SKIP_ALL = true;
        else
            fprintf('Masks will be recomputed and will replace the old ones.\n');
        end
    end

    if ~SKIP_ALL

        % ---------- LAST FRAME = maximum envelope ----------
        % Every earlier mask is constrained to lie inside this region.
        if MANUAL_LAST
            gLast = im2double(im2gray(imread(fullfile(otsuDir, names{end}))));
            if INVERT_IMAGE, gLast = 1 - gLast; end

            envFile = fullfile(resultsFolder, ['envelope_' dropletName '.mat']);
            drawNew = true;

            if isfile(envFile)
                r = questdlg('A saved envelope already exists. Reuse it?', ...
                             'Envelope', 'Yes, reuse', 'No, redraw', 'Yes, reuse');
                if strcmp(r, 'Yes, reuse')
                    load(envFile, 'lastMask');
                    drawNew = false;
                    fprintf('Envelope reused from %s\n', envFile);
                end
            end

            if drawNew
                lastMask = false(size(gLast));

                hMan = figure('Name','Draw the crystals (last frame)','Color','w');
                imshow(imadjust(gLast)); hold on
                title(sprintf(['Draw up to %d crystals.\n' ...
                               'Double-click to close each polygon. ' ...
                               'Close the figure when finished.'], NCrys));

                nDrawn = 0;
                while nDrawn < NCrys
                    if ~ishandle(hMan), break; end
                    try
                        poly = drawpolygon('Color','r');
                    catch
                        break
                    end
                    if isempty(poly.Position) || size(poly.Position,1) < 3
                        break
                    end
                    bw = poly2mask(poly.Position(:,1), poly.Position(:,2), ...
                                   size(gLast,1), size(gLast,2));
                    lastMask = lastMask | bw;
                    nDrawn = nDrawn + 1;
                    fprintf('Crystal %d drawn.\n', nDrawn);

                    if nDrawn < NCrys
                        more = questdlg('Draw another crystal?', 'Crystals', ...
                                        'Yes','No, finish','Yes');
                        if ~strcmp(more,'Yes'), break; end
                    end
                end
                if ishandle(hMan), close(hMan); end
                fprintf('%d crystals defined manually.\n', nDrawn);

                if USE_POLY_SMOOTH && any(lastMask(:))
                    lastMask = polygonizeMask(lastMask, POLY_TOL);
                end

                save(envFile, 'lastMask');
                fprintf('Envelope saved to %s\n', envFile);
            end
        else
            % Automatic envelope, thresholded rather than drawn.
            lastMask = segFrame(fullfile(otsuDir, names{end}), INVERT_IMAGE, ...
                                GAUSS_SIGMA, OTSU_SCALE, OTSU_MIN_SIZE, ...
                                MORPH_OPEN, MORPH_CLOSE, CONVEX_STRENGTH);
            lastMask = keepNbiggest(lastMask, NCrys);
            if USE_POLY_SMOOTH && any(lastMask(:))
                lastMask = polygonizeMask(lastMask, POLY_TOL);
            end
        end

        % ---------- NUCLEATION FRAMES ----------
        reuse = false;
        if exist('specs','var') && isfield(specs,'nucleation_frames') ...
                                && ~isempty(specs.nucleation_frames) ...
                                && all(~isnan(specs.nucleation_frames))
            prev = specs.nucleation_frames;
            r = questdlg(sprintf(['Nucleation frames are already stored:\n%s' ...
                         '\n\nReuse them?'], mat2str(prev)), 'Nucleation', ...
                         'Yes, reuse', 'No, mark again', 'Yes, reuse');
            if strcmp(r, 'Yes, reuse')
                nucList = sort(prev(:)');
                reuse = true;
            end
        end

        if ~reuse
            nucList = markNucFrames(otsuDir, names, frameTimes, NCrys);
            nucList = sort(nucList(:)');
        end
        fprintf('Nucleation frames: %s\n', mat2str(nucList));

        % ---------- SEGMENT BACKWARDS, LAST FRAME TO FIRST ----------
        masks     = cell(totalFrames, 1);
        areasPx   = cell(totalFrames, 1);
        areasUm   = cell(totalFrames, 1);
        emptyMask = false(size(lastMask));

        prevMask = lastMask;            % start from the drawn envelope
        prevArea = nnz(lastMask);

        for i = totalFrames:-1:1

            % How many crystals have nucleated by this frame.
            nActive = sum(i >= nucList);

            if nActive == 0
                % Before the first nucleation there is nothing to find.
                masks{i}   = emptyMask;
                areasPx{i} = [];
                areasUm{i} = [];
                prevMask   = emptyMask;
                prevArea   = 0;
                continue
            end

            m = segFrame(fullfile(otsuDir, names{i}), INVERT_IMAGE, ...
                         GAUSS_SIGMA, OTSU_SCALE, OTSU_MIN_SIZE, ...
                         MORPH_OPEN, MORPH_CLOSE, CONVEX_STRENGTH);

            % Intersect with the already-processed following frame. The
            % dilation gives the crystal room to have moved slightly between
            % frames without the intersection cutting it away.
            tol = 35;
            fm  = m & imdilate(prevMask, strel('disk', tol));

            % Keep only as many objects as have nucleated so far.
            fm = keepNbiggest(fm, nActive);

            % Going backwards, the total area cannot increase. If the
            % threshold picked up more than the following frame had, erode
            % until it does not. The iteration cap stops a runaway erosion
            % from deleting the crystal entirely on a bad frame.
            nIter = 0;
            while nnz(fm) > prevArea && any(fm(:)) && nIter < 15
                fm = imerode(fm, strel('disk', 1));
                nIter = nIter + 1;
            end

            if USE_POLY_SMOOTH && any(fm(:))
                fm = polygonizeMask(fm, POLY_TOL);
            end

            masks{i} = fm;
            prevMask = fm;              % this frame anchors the previous one
            prevArea = nnz(fm);

            if any(fm(:))
                st = regionprops(bwlabel(fm), 'Area');
                areasPx{i} = [st.Area];
                areasUm{i} = [st.Area] * area_scale_um2;
            else
                areasPx{i} = []; areasUm{i} = [];
            end

            if mod(i,50) == 0, fprintf('  frame %d/%d\n', i, totalFrames); end
        end

        % ---------- save (replaces any previous masks) ----------
        save(fullfile(resultsFolder, ['masks_' dropletName '.mat']), ...
             'masks', '-v7.3');
        save(fullfile(resultsFolder, ['crystal_areas_' dropletName '.mat']), ...
             'areasPx','areasUm','frameTimes','nucList');
        fprintf('Saved in %s\n', resultsFolder);

        % ---------- write the nucleation frames back to the specs ----------
        % So that they only ever have to be marked by hand once.
        if exist('specs','var') && isfield(specs,'file')
            nucStr = strjoin(string(nucList), ',');
            updateSpecs(specs.file, struct( ...
                'nucleation_frames', char(nucStr), ...
                'nucleation_frame',  nucList(1)));   % first crystal = global
            fprintf('Nucleation saved: frames=%s, first frame=%d\n', ...
                    nucStr, nucList(1));
        else
            warning(['Nucleation frames were not written to the ' ...
                     'specifications (the specs variable is missing).']);
        end

    end   % end of if ~SKIP_ALL

    % ---------- summary figure, new or loaded masks alike ----------
    dispIdx = round(linspace(1, totalFrames, DISPLAY_N));
    figS = figure('Name','Crystal detection','Color','w', ...
                  'Position',[100 100 1400 800]);
    for i = 1:DISPLAY_N
        fi = dispIdx(i);
        g  = im2double(im2gray(imread(fullfile(otsuDir, names{fi}))));
        mk = masks{fi};
        ov = repmat(g, [1 1 3]);
        ov(:,:,1) = min(ov(:,:,1) + 0.6*double(mk), 1);   % mask tinted red
        subplot(3, DISPLAY_N, i);             imshow(g);  title(sprintf('Frame %d', fi));
        subplot(3, DISPLAY_N, i+DISPLAY_N);   imshow(mk); title(sprintf('%d obj', numel(areasPx{fi})));
        subplot(3, DISPLAY_N, i+2*DISPLAY_N); imshow(ov); title('Overlay');
    end
    exportgraphics(figS, ...
        fullfile(resultsFolder, ['maskssummary_' dropletName '.png']), ...
        'Resolution', 200);
    fprintf('Detection finished.\n');

    % ---------- optional playback of the whole mask sequence ----------
    % Visual check only, nothing is saved.
    if exist('SHOW_MASK_MOVIE','var') && SHOW_MASK_MOVIE
        hMov = figure('Name','Masks (start -> end)','Color','w');
        for i = 1:totalFrames
            if ~ishandle(hMov), break; end
            imshow(masks{i});
            title(sprintf('Frame %d / %d   t = %.2f s', i, totalFrames, frameTimes(i)));
            drawnow;
        end
    end

else
    % ---------- segmentation disabled: load the masks from disk ----------
    load(fullfile(resultsFolder, ['masks_' dropletName '.mat']), 'masks');
    af = fullfile(resultsFolder, ['crystal_areas_' dropletName '.mat']);
    if isfile(af)
        Sld = load(af, 'areasPx','areasUm','nucList');   % again, no frameTimes
        areasPx = Sld.areasPx; areasUm = Sld.areasUm; nucList = Sld.nucList;
        clear Sld
    end
    fprintf('Segmentation disabled: masks loaded from disk.\n');
end

% ---------- reconcile lengths and take nucList from the specifications ----
totalFrames = min([totalFrames, numel(masks), numel(names)]);
names       = names(1:totalFrames);
masks       = masks(1:totalFrames);
frameTimes  = t0 + (0:totalFrames-1) / fps;
if ~isempty(nucSpecs), nucList = nucSpecs; end
nucList = sort(double(nucList(:)'));
nucList(nucList < 1 | nucList > totalFrames) = [];
assert(~isempty(nucList), 'No valid nucleation frames.');
fprintf('Axis: %d frames, t = %.2f -> %.2f s | nucleation at t = %s s\n', ...
        totalFrames, frameTimes(1), frameTimes(end), ...
        mat2str(round(frameTimes(nucList),2)));

%% ========================================================================
%  PART 2 : GLOBAL GROWTH RATE AND FRONT VELOCITY
%  ========================================================================
nCryst = numel(nucList);
colors = crystalColors;   % green range

if DO_GLOBAL_VELOCITIES

    if ~exist('masks','var')
        load(fullfile(resultsFolder, ['masks_' dropletName '.mat']), 'masks');
    end

    % ---------- plotting window for the crystal-only figures ----------
    %  There is nothing to show before the first crystal exists, so these
    %  plots start at N1 rather than at t_start. The figures that compare
    %  against PIV or TrackMate DO start at t_start, because there the point
    %  is to see the flow before anything nucleates.
    tPlot0 = frameTimes(min(nucList(1), numel(frameTimes)));
    tWin   = [tPlot0, frameTimes(end)];

    % ---------- total crystal area and perimeter per frame ----------
    totalArea_um2 = zeros(numel(masks),1);
    totalPerim_um = zeros(numel(masks),1);
    for i = 1:numel(masks)
        totalArea_um2(i) = nnz(masks{i}) * area_scale_um2;
        if any(masks{i}(:))
            stP = regionprops(bwlabel(masks{i}), 'Perimeter');
            totalPerim_um(i) = sum([stP.Perimeter]) * um_per_px;
        end
    end

    figArea = figure('Name','Total crystal area','Color','w', ...
                     'Position',[100 100 1000 650]);
    plot(frameTimes, totalArea_um2, '-', 'LineWidth', LW_PLOT, 'Color', colorCrystal);
    xlabel('Time (s)'); ylabel('Total crystal area ($\mu$m$^2$)');
    title('Total crystal area'); grid on
    nucleationLines(gca, frameTimes(nucList), crystalColors, true);
    figureStyle(gca, tWin, tPlot0);
    exportgraphics(figArea, ...
        fullfile(resultsFolder, ['total_area_' dropletName '.png']), 'Resolution', 300);
    fprintf('Saved: total_area_%s.png\n', dropletName);

    % ---------- global growth rate, dA/dt ----------
    % Centred difference with a gap of N_FINITE_DIFF frames:
    %   dA/dt|i = ( A(i+N) - A(i-N) ) / ( t(i+N) - t(i-N) )
    % The gap is what makes the derivative usable: over one frame the area
    % change is comparable to the segmentation noise. One-sided differences
    % are used at the two ends.
    Nfd = N_FINITE_DIFF;
    Np  = numel(totalArea_um2);
    growthGlobal = nan(Np,1);

    for i = 1:Np
        if i <= Nfd
            if i + Nfd <= Np
                growthGlobal(i) = (totalArea_um2(i+Nfd) - totalArea_um2(i)) ...
                                / (frameTimes(i+Nfd) - frameTimes(i));
            end
        elseif i > Np - Nfd
            if i - Nfd >= 1
                growthGlobal(i) = (totalArea_um2(i) - totalArea_um2(i-Nfd)) ...
                                / (frameTimes(i) - frameTimes(i-Nfd));
            end
        else
            growthGlobal(i) = (totalArea_um2(i+Nfd) - totalArea_um2(i-Nfd)) ...
                            / (frameTimes(i+Nfd) - frameTimes(i-Nfd));
        end
    end

    growthGlobal = smoothdata(growthGlobal, 'movmean', SMOOTH_WINDOW);
    growthGlobal(growthGlobal < 0) = 0;   % crystals do not shrink
    % Before the first nucleation there is no crystal, so the curve starts
    % there rather than sitting at zero.
    preNuc = (1:Np)' < nucList(1);
    growthGlobal(preNuc) = NaN;
    nucT_glob = frameTimes(nucList);

    figGrowth = figure('Name','Global growth rate','Color','w', ...
                       'Position',[100 100 1000 650]);
    plot(frameTimes, growthGlobal, '-', 'LineWidth', LW_PLOT, 'Color', colorCrystal);
    xlabel('Time (s)'); ylabel('Growth rate ($\mu$m$^2$/s)');
    title('Global crystal growth rate'); grid on
    nucleationLines(gca, nucT_glob, crystalColors, true);
    figureStyle(gca, tWin, tPlot0);
    exportgraphics(figGrowth, ...
        fullfile(resultsFolder, ['global_growth_rate_' dropletName '.png']), ...
        'Resolution', 300);
    fprintf('Saved: global_growth_rate_%s.png\n', dropletName);

    % ====================================================================
    %  GLOBAL FRONT VELOCITY (um/s), two alternative definitions
    % ====================================================================
    %  The growth rate is an area per unit time; the flow is a length per
    %  unit time. Comparing them requires converting the first into a front
    %  velocity, and there is more than one way to do that.

    % ---------- (A) dA/dt normalised by the real perimeter ----------
    %  Mean advance of the interface per unit interface length, since
    %  dA/dt = integral of v_n along the perimeter. This is the definition
    %  used in the thesis.
    %
    %  growthGlobal is already smoothed, but dividing by a RAW perimeter
    %  reintroduces the noise, because the perimeter is the quantity most
    %  sensitive to segmentation. Measured on one droplet:
    %     dA/dt smoothed  1.59%   P raw  1.95%   ->  quotient 2.63%
    %     with P smoothed                        ->  quotient 1.82%
    %  So the perimeter is smoothed before dividing. totalPerim_um is still
    %  saved unsmoothed.
    SMOOTH_PERIMETER = 20;
    perim_sm = smoothdata(totalPerim_um, 'movmean', SMOOTH_PERIMETER, 'omitnan');

    frontVel_perim = nan(Np,1);
    valid = perim_sm > 0;
    frontVel_perim(valid) = growthGlobal(valid) ./ perim_sm(valid);
    frontVel_perim(frontVel_perim < 0) = 0;
    frontVel_perim(preNuc) = NaN;

    % ---------- (B) derivative of the equivalent radius ----------
    %  The rate at which the radius would grow if all the material formed a
    %  single circular crystal. Insensitive to shape, depending only on area,
    %  which is the more robust measurement. Its drawback is that it collapses
    %  every crystal into one circle: with N crystals of radius r,
    %  A = N*pi*r^2 gives R_eq = r*sqrt(N), so dR_eq/dt = sqrt(N)*dr/dt. The
    %  velocity is inflated by sqrt(N) without any crystal growing faster,
    %  and jumps artificially at each nucleation. Kept here for comparison.
    Req_um = sqrt(totalArea_um2 / pi);
    frontVel_Req = nan(Np,1);
    for i = 1:Np
        if i <= Nfd
            if i + Nfd <= Np
                frontVel_Req(i) = (Req_um(i+Nfd) - Req_um(i)) ...
                                / (frameTimes(i+Nfd) - frameTimes(i));
            end
        elseif i > Np - Nfd
            if i - Nfd >= 1
                frontVel_Req(i) = (Req_um(i) - Req_um(i-Nfd)) ...
                                / (frameTimes(i) - frameTimes(i-Nfd));
            end
        else
            frontVel_Req(i) = (Req_um(i+Nfd) - Req_um(i-Nfd)) ...
                            / (frameTimes(i+Nfd) - frameTimes(i-Nfd));
        end
    end
    frontVel_Req = smoothdata(frontVel_Req, 'movmean', SMOOTH_WINDOW);
    frontVel_Req(frontVel_Req < 0) = 0;
    frontVel_Req(preNuc) = NaN;

    % ---------- the two definitions side by side ----------
    figFront = figure('Name','Global front velocity','Color','w', ...
                      'Position',[100 100 1100 650]);
    plot(frameTimes, frontVel_perim, '-', 'LineWidth', LW_PLOT, ...
         'Color', [0.20 0.40 0.80]); hold on
    plot(frameTimes, frontVel_Req,   '-', 'LineWidth', LW_PLOT, ...
         'Color', [0.80 0.25 0.20])
    xlabel('Time (s)'); ylabel('Front velocity ($\mu$m/s)');
    title('Global front velocity: two definitions'); grid on
    legend({'(A) $\dot{A}/P_{\mathrm{tot}}$', ...
            '(B) $\mathrm{d}R_{\mathrm{eq}}/\mathrm{d}t$'}, 'Location','best');
    nucleationLines(gca, nucT_glob, crystalColors, true);
    figureStyle(gca, tWin, tPlot0);
    exportgraphics(figFront, ...
        fullfile(resultsFolder, ['front_velocity_' dropletName '.png']), ...
        'Resolution', 300);
    fprintf('Saved: front_velocity_%s.png\n', dropletName);

    save(fullfile(resultsFolder, ['growth_global_' dropletName '.mat']), ...
         'totalArea_um2', 'totalPerim_um', 'growthGlobal', 'frameTimes', ...
         'Req_um', 'frontVel_perim', 'frontVel_Req');
end

%% ========================================================================
%  PART 3 : PER-CRYSTAL MASKS AND BUFFER RINGS
%  ========================================================================
%  The global mask lumps all crystals together. To relate the growth of one
%  crystal to the flow around THAT crystal, each has to be isolated, and a
%  region of influence has to be defined around it.
%
%  The buffer ring is that region: a disc of radius bufferFactor times the
%  crystal's own mean radius, centred on its centroid. The flow measured
%  inside it is what gets compared against that crystal's growth.

if DO_INDIVIDUAL_MASKS

    finalMask = masks{end};

    % ---------- draw one region per crystal, on the final mask ----------
    figure('Name','Draw regions','Color','w')
    imshow(finalMask)
    title('Draw one region per crystal')

    regionMasks = cell(NCrys,1);
    for c = 1:NCrys
        fprintf('Draw region for crystal %d\n', c)
        h = drawpolygon;
        regionMasks{c} = createMask(h);
    end

    manualRegionsFig = fullfile(resultsFolder, ['manual_regions_' dropletName '.png']);
    if exist(manualRegionsFig, 'file'), delete(manualRegionsFig); end
    exportgraphics(gcf, manualRegionsFig, 'Resolution', 300);
    fprintf('Saved figure: %s\n', manualRegionsFig);

    % ---------- order the regions by nucleation ----------
    %  The regions are drawn in whatever order is convenient, but every plot
    %  in the thesis assumes "Crystal 1" is the one that nucleated first, and
    %  that it keeps the same colour everywhere. Sorting by first appearance
    %  enforces that.
    firstSeen = inf(1, NCrys);
    for c = 1:NCrys
        for i = 1:totalFrames
            if any(masks{i}(:) & regionMasks{c}(:))
                firstSeen(c) = i;
                break
            end
        end
    end
    [~, order]  = sort(firstSeen);
    regionMasks = regionMasks(order);
    fprintf('Regions reordered by nucleation: %s\n', mat2str(order));

    % ---------- split the global mask by region ----------
    crystalMasks = cell(totalFrames, NCrys);
    for frameIdx = 1:totalFrames
        mask = masks{frameIdx};
        for c = 1:NCrys
            crystalMasks{frameIdx, c} = mask & regionMasks{c};
        end
    end

    % ---------- visual check: six frames per crystal ----------
    Nshow = 6;
    frames_to_show = round(linspace(1, totalFrames, Nshow));
    for c = 1:NCrys
        figure('Name', sprintf('Crystal %d', c), 'Position', [100 100 1200 500])
        for i = 1:Nshow
            frameIdx = frames_to_show(i);
            subplot(2,3,i)
            imshow(crystalMasks{frameIdx,c})
            title(sprintf('Frame %d', frameIdx))
        end
        sgtitle(sprintf('Crystal %d (manual segmentation)', c))
        figName = fullfile(resultsFolder, ...
                  sprintf('crystal_%d_manual_segmentation_%s.png', c, dropletName));
        if exist(figName, 'file'), delete(figName); end
        exportgraphics(gcf, figName, 'Resolution', 300);
        fprintf('Saved figure: %s\n', figName);
    end

    % ---------- centroid of each crystal in each frame ----------
    % The largest object is taken when a region briefly contains more than
    % one, which happens while a crystal is still fragmenting into view.
    centroids = nan(totalFrames, NCrys, 2);
    for frameIdx = 1:totalFrames
        for c = 1:NCrys
            mask = crystalMasks{frameIdx, c};
            if any(mask(:))
                stats = regionprops(mask, 'Area', 'Centroid');
                if ~isempty(stats)
                    [~, idxMax] = max([stats.Area]);
                    centroids(frameIdx, c, :) = stats(idxMax).Centroid;
                end
            end
        end
    end

    % ---------- mean radius from the outline ----------
    %  Mean distance from the centroid to the boundary pixels. This is a real
    %  measurement of the crystal outline rather than an equivalent radius
    %  derived from the area, so it stays meaningful for the non-circular
    %  shapes the crystals actually have.
    radius_mean = nan(totalFrames, NCrys);
    for frameIdx = 1:totalFrames
        for c = 1:NCrys
            mask = crystalMasks{frameIdx, c};
            if any(mask(:))
                cx = centroids(frameIdx, c, 1);
                cy = centroids(frameIdx, c, 2);
                edge = bwperim(mask);
                [y, x] = find(edge);
                dist = sqrt((x - cx).^2 + (y - cy).^2);
                radius_mean(frameIdx, c) = mean(dist);
            end
        end
    end

    % ---------- buffer ring: the region of influence ----------
    maskPlusBuffer = cell(totalFrames, NCrys);
    for frameIdx = 1:totalFrames
        [Ny, Nx] = size(masks{frameIdx});
        [X, Y] = meshgrid(1:Nx, 1:Ny);
        for c = 1:NCrys
            cx = centroids(frameIdx, c, 1);
            cy = centroids(frameIdx, c, 2);
            R  = radius_mean(frameIdx, c);
            if ~isnan(cx) && ~isnan(R)
                dist = sqrt((X - cx).^2 + (Y - cy).^2);
                maskPlusBuffer{frameIdx, c} = dist <= bufferFactor * R;
            else
                maskPlusBuffer{frameIdx, c} = false(Ny, Nx);
            end
        end
    end

    % ---------- visual check: mask, centroid and buffer together ----------
    Nshow = 6;
    framesToShow = round(linspace(1, totalFrames, Nshow));
    colorsK = crystalColors;   % green range: "green = crystal" everywhere
    for c = 1:NCrys
        figure('Name', sprintf('Crystal %d - Mask + Centroid + Buffer', c), ...
               'Position', [100 100 1200 400])
        for i = 1:Nshow
            frameIdx = framesToShow(i);
            mask   = crystalMasks{frameIdx, c};
            buffer = maskPlusBuffer{frameIdx, c};
            subplot(1, Nshow, i)
            imshow(mask); hold on
            cx = centroids(frameIdx, c, 1);
            cy = centroids(frameIdx, c, 2);
            if ~isnan(cx)
                plot(cx, cy, 'ro', 'MarkerSize', 2, 'LineWidth', 2)
            end
            if any(mask(:))
                visboundaries(mask, 'Color', colorsK(c,:), 'LineWidth', 2)
            end
            if ~isempty(buffer) && any(buffer(:))
                visboundaries(buffer, 'Color', colorsK(c,:), ...
                    'LineStyle','--', 'LineWidth', 1.5)
            end
            title(sprintf('F%d', frameIdx))
        end
        sgtitle(sprintf('Crystal %d (mask + centroid + buffer)', c), 'FontSize', 14)
        figName = fullfile(resultsFolder, ...
                  sprintf('crystal_%d_mask_centroid_buffer_%s.png', c, dropletName));
        if exist(figName, 'file'), delete(figName); end
        exportgraphics(gcf, figName, 'Resolution', 300);
        fprintf('Saved figure: %s\n', figName);
    end

    save(fullfile(resultsFolder, ['individual_crystals_' dropletName '.mat']), ...
         'regionMasks', 'crystalMasks', 'centroids', 'radius_mean', ...
         'maskPlusBuffer', 'bufferFactor', '-v7.3');
    fprintf('Individual crystal masks saved.\n');
end

%% ========================================================================
%  PART 4 : PER-CRYSTAL GROWTH RATE AND RADIAL VELOCITY
%  ========================================================================
if DO_INDIVIDUAL_VELOCITIES

    % ---------- discard anything in memory from another droplet ----------
    % Without this, a stale workspace silently supplies the previous
    % droplet's masks to the "load only if missing" blocks below.
    if ~exist('dataDroplet','var') || ~strcmp(dataDroplet, dropletName)
        clear crystalMasks radius_mean centroids maskPlusBuffer
        dataDroplet = dropletName;
    end

    % ---------- load everything part 3 produces, in one go ----------
    % All four variables are loaded together: loading only some of them left
    % the others undefined and failed later.
    indFile = fullfile(resultsFolder, ['individual_crystals_' dropletName '.mat']);
    needInd = {'crystalMasks','radius_mean','centroids','maskPlusBuffer'};
    missing = {};
    for vi = 1:numel(needInd)
        if ~exist(needInd{vi}, 'var'), missing{end+1} = needInd{vi}; end %#ok<AGROW>
    end
    if ~isempty(missing)
        assert(isfile(indFile), ['%s is missing.\n' ...
               'Set DO_INDIVIDUAL_MASKS = true ONCE to generate it.'], indFile);
        Sld = load(indFile, missing{:});
        for vi = 1:numel(missing)
            eval([missing{vi} ' = Sld.' missing{vi} ';']);   %#ok<EVLEQ>
        end
        clear Sld
        fprintf('individual_crystals loaded from disk: %s\n', strjoin(missing, ', '));
    end

    % ---------- reconcile lengths with the pipeline time axis ----------
    nRowsInd    = size(crystalMasks, 1);
    totalFrames = min([nRowsInd, numel(frameTimes)]);
    if nRowsInd ~= numel(frameTimes)
        warning(['individual_crystals has %d frames and the axis has %d. ' ...
                 'Using %d. If this does not add up, redo ' ...
                 'DO_INDIVIDUAL_MASKS.'], nRowsInd, numel(frameTimes), totalFrames);
    end
    crystalMasks = crystalMasks(1:totalFrames, :);
    frameTimes   = frameTimes(1:totalFrames);
    if size(radius_mean,1)    >= totalFrames, radius_mean    = radius_mean(1:totalFrames,:);    end
    if size(centroids,1)      >= totalFrames, centroids      = centroids(1:totalFrames,:,:);    end
    if size(maskPlusBuffer,1) >= totalFrames, maskPlusBuffer = maskPlusBuffer(1:totalFrames,:); end

    % ---------- Ncrys: the specs win, but cannot exceed what is on disk ----
    NCrysDisk = size(crystalMasks, 2);
    if NCrys > NCrysDisk
        warning(['Ncrys=%d in the specifications but individual_crystals ' ...
                 'only has %d crystals. Using %d.'], NCrys, NCrysDisk, NCrysDisk);
        NCrys = NCrysDisk;
        crystalColors = [linspace(0.05,0.55,NCrys)', ...
                         linspace(0.45,0.80,NCrys)', ...
                         linspace(0.22,0.45,NCrys)'];
    end

    % ---------- nucleation frames clipped to this totalFrames ----------
    nucInd = sort(double(nucList(:)'));
    nucInd(nucInd < 1) = 1;
    nucInd(nucInd > totalFrames) = totalFrames;
    nucT = frameTimes(nucInd);
    fprintf('Individual: %d frames, %d crystals, nucleation at t = %s s\n', ...
            totalFrames, NCrys, mat2str(round(nucT,2)));

    % ---------- area of each crystal in each frame ----------
    areaCrystal_um2 = zeros(totalFrames, NCrys);
    for frameIdx = 1:totalFrames
        for c = 1:NCrys
            areaCrystal_um2(frameIdx, c) = ...
                nnz(crystalMasks{frameIdx, c}) * area_scale_um2;
        end
    end

    % ---------- growth rate of each crystal ----------
    % Same method as the global one: centred difference with a gap of
    % N_FINITE_DIFF, one-sided at the ends, then a moving average.
    Nfd = N_FINITE_DIFF;
    growthCrystal = nan(totalFrames, NCrys);

    for c = 1:NCrys
        A = areaCrystal_um2(:, c);
        g = nan(totalFrames, 1);
        for i = 1:totalFrames
            if i <= Nfd
                if i + Nfd <= totalFrames
                    g(i) = (A(i+Nfd) - A(i)) / (frameTimes(i+Nfd) - frameTimes(i));
                end
            elseif i > totalFrames - Nfd
                if i - Nfd >= 1
                    g(i) = (A(i) - A(i-Nfd)) / (frameTimes(i) - frameTimes(i-Nfd));
                end
            else
                g(i) = (A(i+Nfd) - A(i-Nfd)) / (frameTimes(i+Nfd) - frameTimes(i-Nfd));
            end
        end
        gs = smoothdata(g, 'movmean', SMOOTH_WINDOW);
        gs(gs < 0) = 0;                          % no negative growth
        if c <= numel(nucInd) && nucInd(c) > 1   % NaN before ITS own nucleation
            gs(1:min(nucInd(c)-1, numel(gs))) = NaN;
        end
        growthCrystal(:, c) = gs;
    end

    colorsG = crystalColors;   % same colour per crystal across the whole project

    % ---------- plotting window ----------
    tPlot0 = nucT(1);                      % first nucleation of this droplet
    tWin   = [tPlot0, frameTimes(end)];

    % ---------- area of each crystal over time ----------
    figAreaInd = figure('Name','Individual crystal area','Color','w', ...
                        'Position',[100 100 1000 650]);
    hold on
    for c = 1:NCrys
        plot(frameTimes, areaCrystal_um2(:, c), 'LineWidth', LW_PLOT, ...
             'Color', colorsG(c,:), 'DisplayName', sprintf('Crystal %d', c));
    end
    xlabel('Time (s)'); ylabel('Crystal area ($\mu$m$^2$)');
    title('Individual crystal area'); grid on
    legend('Location','northwest');
    nucleationLines(gca, nucT, crystalColors, true);
    figureStyle(gca, tWin, tPlot0);
    exportgraphics(figAreaInd, ...
        fullfile(resultsFolder, ['individual_area_' dropletName '.png']), ...
        'Resolution', 300);
    fprintf('Saved: individual_area_%s.png\n', dropletName);

    % ---------- growth rate of each crystal ----------
    figGrowInd = figure('Name','Individual growth rate','Color','w', ...
                        'Position',[100 100 1000 650]);
    hold on
    for c = 1:NCrys
        plot(frameTimes, growthCrystal(:, c), 'LineWidth', LW_PLOT, ...
             'Color', colorsG(c,:), 'DisplayName', sprintf('Crystal %d', c));
    end
    xlabel('Time (s)'); ylabel('Growth rate ($\mu$m$^2$/s)');
    title('Individual crystal growth rate'); grid on
    legend('Location','northwest');
    nucleationLines(gca, nucT, crystalColors, true);
    figureStyle(gca, tWin, tPlot0);
    exportgraphics(figGrowInd, ...
        fullfile(resultsFolder, ['individual_growth_rate_' dropletName '.png']), ...
        'Resolution', 300);
    fprintf('Saved: individual_growth_rate_%s.png\n', dropletName);

    save(fullfile(resultsFolder, ['growth_individual_' dropletName '.mat']), ...
         'areaCrystal_um2', 'growthCrystal', 'frameTimes');
    fprintf('Individual growth analysis done.\n');

    % ====================================================================
    %  RADIAL GROWTH VELOCITY (um/s)
    % ====================================================================
    %  dR/dt of each crystal's mean radius. Unlike the growth rate, which is
    %  an area per unit time, this is a length per unit time and therefore
    %  directly comparable with the PIV and TrackMate velocities.
    %
    %  The global value is the mean over the crystals that EXIST in each
    %  frame, not over all of them: averaging in the crystals that have not
    %  nucleated yet would drag the curve down for no physical reason.

    radius_um = radius_mean * um_per_px;      % [totalFrames x NCrys]

    Nfd = N_FINITE_DIFF;
    vRadial = nan(totalFrames, NCrys);

    for c = 1:NCrys
        R = radius_um(:, c);
        g = nan(totalFrames, 1);
        for i = 1:totalFrames
            if i <= Nfd
                if i + Nfd <= totalFrames && ~isnan(R(i)) && ~isnan(R(i+Nfd))
                    g(i) = (R(i+Nfd) - R(i)) / (frameTimes(i+Nfd) - frameTimes(i));
                end
            elseif i > totalFrames - Nfd
                if i - Nfd >= 1 && ~isnan(R(i)) && ~isnan(R(i-Nfd))
                    g(i) = (R(i) - R(i-Nfd)) / (frameTimes(i) - frameTimes(i-Nfd));
                end
            else
                if ~isnan(R(i+Nfd)) && ~isnan(R(i-Nfd))
                    g(i) = (R(i+Nfd) - R(i-Nfd)) / (frameTimes(i+Nfd) - frameTimes(i-Nfd));
                end
            end
        end
        gs = smoothdata(g, 'movmean', SMOOTH_WINDOW, 'omitnan');
        gs(gs < 0)   = 0;                        % no negative growth
        gs(isnan(R)) = NaN;                      % no mask means no radius
        if c <= numel(nucInd) && nucInd(c) > 1   % NaN before ITS own nucleation
            gs(1:min(nucInd(c)-1, numel(gs))) = NaN;
        end
        vRadial(:, c) = gs;
    end

    % ---------- global: mean over the crystals that exist ----------
    % A crystal "exists" in a frame if it has a radius greater than zero.
    vRadialGlobal = nan(totalFrames, 1);
    for i = 1:totalFrames
        exists = radius_um(i, :) > 0 & ~isnan(vRadial(i, :));
        if any(exists)
            vRadialGlobal(i) = mean(vRadial(i, exists), 'omitnan');
        end
    end
    vRadialGlobal(vRadialGlobal < 0) = 0;
    vRadialGlobal((1:totalFrames)' < nucInd(1)) = NaN;   % before the first one

    % ---------- radial velocity of each crystal ----------
    figVR = figure('Name','Radial growth velocity (per crystal)','Color','w', ...
                   'Position',[100 100 1000 650]);
    hold on
    for c = 1:NCrys
        plot(frameTimes, vRadial(:, c), 'LineWidth', LW_PLOT, ...
             'Color', crystalColors(c,:), 'DisplayName', sprintf('Crystal %d', c));
    end
    xlabel('Time (s)'); ylabel('Radial growth velocity ($\mu$m/s)');
    title('Radial growth velocity per crystal'); grid on
    legend('Location','northwest');
    nucleationLines(gca, nucT, crystalColors, true);
    figureStyle(gca, tWin, tPlot0);
    exportgraphics(figVR, ...
        fullfile(resultsFolder, ['radial_velocity_individual_' dropletName '.png']), ...
        'Resolution', 300);
    fprintf('Saved: radial_velocity_individual_%s.png\n', dropletName);

    % ---------- global radial velocity ----------
    figVRG = figure('Name','Global radial growth velocity','Color','w', ...
                    'Position',[100 100 1000 650]);
    plot(frameTimes, vRadialGlobal, '-', 'LineWidth', LW_PLOT, 'Color', colorCrystal);
    xlabel('Time (s)'); ylabel('Mean radial growth velocity ($\mu$m/s)');
    title('Global radial growth velocity'); grid on
    nucleationLines(gca, nucT, crystalColors, true);
    figureStyle(gca, tWin, tPlot0);
    exportgraphics(figVRG, ...
        fullfile(resultsFolder, ['radial_velocity_global_' dropletName '.png']), ...
        'Resolution', 300);
    fprintf('Saved: radial_velocity_global_%s.png\n', dropletName);

    save(fullfile(resultsFolder, ['radial_velocity_' dropletName '.mat']), ...
         'radius_um', 'vRadial', 'vRadialGlobal', 'frameTimes');
    fprintf('Radial growth velocity saved.\n');
end


%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================

function nucList = markNucFrames(otsuDir, names, frameTimes, NCrys)
%MARKNUCFRAMES  Mark the nucleation frame of each crystal by hand.
%
%   Nucleation is not reliably detectable automatically: the first frames of
%   a crystal are a handful of faint pixels indistinguishable from noise by
%   any threshold. Scrubbing through the sequence and marking the frame where
%   a new crystal appears is both faster and more accurate.
%
%   Returns up to NCrys frame indices.

N = numel(names);
nucList = [];
frameNow = 1;

hf = figure('Name','Mark the nucleation frames','Color','w', ...
            'Position',[80 80 950 780]);
ax = axes('Parent',hf,'Position',[0.05 0.20 0.9 0.72]);
drawNow();

sld = uicontrol('Parent',hf,'Style','slider','Units','normalized', ...
    'Position',[0.05 0.11 0.9 0.05],'Min',1,'Max',N,'Value',1, ...
    'SliderStep',[1/(N-1) 20/(N-1)]);
addlistener(sld,'Value','PostSet',@(~,~) setFrame(round(get(sld,'Value'))));

uicontrol('Parent',hf,'Style','pushbutton','String','A crystal nucleates here', ...
    'Units','normalized','Position',[0.06 0.03 0.40 0.06],'FontWeight','bold', ...
    'Callback',@(~,~) addFrame());
uicontrol('Parent',hf,'Style','pushbutton','String','Finish', ...
    'Units','normalized','Position',[0.54 0.03 0.40 0.06],'FontWeight','bold', ...
    'Callback',@(~,~) uiresume(hf));

uiwait(hf);
if ishandle(hf), close(hf); end

    function setFrame(k), frameNow = k; drawNow(); end

    function addFrame()
        nucList(end+1) = frameNow; %#ok<AGROW>
        fprintf('  Crystal %d nucleates at frame %d\n', numel(nucList), frameNow);
        if numel(nucList) >= NCrys, uiresume(hf); end
    end

    function drawNow()
        g = im2double(im2gray(imread(fullfile(otsuDir, names{frameNow}))));
        imshow(imadjust(g),'Parent',ax);
        title(ax, sprintf(['Scrub through and click when a new crystal ' ...
              'appears.\nFrame %d / %d   t = %.2f s   (marked: %d)'], ...
              frameNow, N, frameTimes(frameNow), numel(nucList)));
    end
end


function mask = segFrame(imgPath, doInvert, sigma, otsuScale, minSize, ...
                         morphOpen, morphClose, convexStrength)
%SEGFRAME  Segment the crystals in a single frame.
%
%   Gaussian blur, then Otsu thresholding scaled by otsuScale, then
%   morphological cleanup. The scale factor on the Otsu threshold is what
%   makes this usable across sessions: the automatic threshold alone is
%   biased by how much of the frame the crystal occupies, which changes
%   completely between the start and the end of a sequence.

g = im2double(im2gray(imread(imgPath)));
if doInvert, g = 1 - g; end
g = imgaussfilt(g, sigma);
mask = imbinarize(g, graythresh(g)*otsuScale);
mask = imfill(mask, 'holes');
mask = bwareaopen(mask, minSize);
if morphOpen  > 0, mask = imopen(mask,  strel('disk', morphOpen));  end
if morphClose > 0, mask = imclose(mask, strel('disk', morphClose)); end
mask = imfill(mask, 'holes');

% Light convex filling. Crystals are faceted and largely convex, but
% thresholding leaves bites out of their edges where the contrast is poor.
% Only the missing pixels close to the existing mask are filled in, so a
% genuinely concave shape is not turned into its convex hull.
if convexStrength > 0 && any(mask(:))
    hull    = bwconvhull(mask, 'objects');
    missing = hull & ~mask;
    nearby  = missing & imdilate(mask, strel('disk', convexStrength));
    mask = mask | nearby;
end
end


function mask = keepNbiggest(mask, N)
%KEEPNBIGGEST  Keep only the N largest connected components.
%   Used to discard the residual specks that survive the morphological
%   cleanup, keeping exactly as many objects as crystals have nucleated.

[labeled, numObj] = bwlabel(mask);
if numObj == 0, return; end
stats = regionprops(labeled, 'Area');
[~, sIdx] = sort([stats.Area], 'descend');
keep = min(N, numObj);
newMask = false(size(mask));
for k = 1:keep, newMask = newMask | (labeled == sIdx(k)); end
mask = newMask;
end


function out = polygonizeMask(mask, tol)
%POLYGONIZEMASK  Replace each object outline by a simplified polygon.
%   Optional. Straightens the pixellated boundary, which mainly matters for
%   the perimeter, the quantity most sensitive to segmentation noise.

out = false(size(mask));
B = bwboundaries(mask, 'noholes');
for k = 1:numel(B)
    b = B{k};
    p = reducepoly(b, tol);
    if size(p,1) < 3, continue; end
    bw = poly2mask(p(:,2), p(:,1), size(mask,1), size(mask,2));
    out = out | bw;
end
end


function updateSpecs(specPath, newFields)
%UPDATESPECS  Update or add keys in the specifications file, in place.
%
%   Existing keys are overwritten and new ones appended; every other line is
%   left exactly as it was, including comments and the fields this pipeline
%   does not use.

fid = fopen(specPath, 'r', 'n', 'UTF-8');
raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);
lines = raw{1};
keys = fieldnames(newFields);

for k = 1:numel(keys)
    key = keys{k};
    val = newFields.(key);
    if isnumeric(val), valStr = num2str(val); else, valStr = char(val); end
    newLine = sprintf('%s=%s', key, valStr);

    found = false;
    for i = 1:numel(lines)
        L = strtrim(lines{i});
        if isempty(L) || ~contains(L,'='), continue; end
        p = strsplit(L,'=');
        if strcmpi(strtrim(p{1}), key)
            lines{i} = newLine; found = true; break
        end
    end
    if ~found, lines{end+1} = newLine; end %#ok<AGROW>
end

fid = fopen(specPath, 'w', 'n', 'UTF-8');
for i = 1:numel(lines), fprintf(fid, '%s\n', lines{i}); end
fclose(fid);
end
