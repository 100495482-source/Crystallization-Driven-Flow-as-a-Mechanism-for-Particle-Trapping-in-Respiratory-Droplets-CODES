%% ========================================================================
%  BUFFERZONEFIGURE  Illustrate the per-crystal region of influence.
%  ========================================================================
%
%  PRODUCES: FIGURES ONLY. Nothing is recomputed and no data file is
%  written; the masks, centroids and radii are read from the .mat that
%  crystalGrowthAnalysis already produced.
%
%  WHAT IT DRAWS
%    Two or more panels of the SAME droplet at different times. In each:
%      - the outline of every crystal, one colour per crystal
%      - its centroid marked
%      - its region of influence, a disc of radius beta*R, dashed
%
%    This is the figure that makes the buffer definition legible. The
%    per-crystal flow analysis averages the velocity inside those dashed
%    circles, and showing them at two times makes clear that the disc grows
%    with the crystal rather than being a fixed region of the image.
%
%  READS PARTIALLY, ON PURPose
%    individual_crystals_*.mat holds one mask per frame per crystal and runs
%    to gigabytes. It is opened with matfile() so only the requested frames
%    are pulled from disk; loading it whole would exhaust memory for the sake
%    of three panels.
%
%  ONLY EDIT THE BLOCK MARKED BELOW.
%
%  OUTPUT
%    <outDir>/buffer_example_<droplet>.pdf  (vector)
%    <outDir>/buffer_example_<droplet>.png
%
%  See also CRYSTALGROWTHANALYSIS, ACCUM3_BUILDBUFFERZONES.

clear; close all; clc;

% ---------- WHICH DROPLET ----------
DATA_ROOT  = fullfile('C:', 'Data', 'FG');
dropletDir = fullfile(DATA_ROOT, 'TP3(1202)', 'S6');
outDir     = dropletDir;
IMGDIR     = fullfile(dropletDir, 'FPICTrackmate');   % '' = autodetect FPIC*
FRAMES     = [200 300 400];   % [] = automatic (2 frames); or e.g. [80 150 220]
CRYSTALS   = [];              % [] = all crystals; or [1 3] for a subset
MAXPX      = 800;             % downsample, so the vector PDF stays workable
ZOOM       = true;            % crop around the crystals
% -----------------------------------

if ~isfolder(outDir), mkdir(outDir); end
resF = fullfile(dropletDir,'RESULTS_crystals');
[~, folderName] = fileparts(dropletDir);

%% ---------- data, read partially ----------
hit = dir(fullfile(resF,'individual_crystals_*.mat'));
assert(~isempty(hit), ...
       'individual_crystals_*.mat not found (run the per-crystal mask stage).');
mfK = matfile(fullfile(resF, hit(1).name));      % partial read, does NOT load all
centroids   = mfK.centroids;
radius_mean = mfK.radius_mean;
beta        = mfK.bufferFactor;
sz = size(mfK,'crystalMasks');
nF = sz(1); NCrys = sz(2);
getM = @(fi,c) mfK.crystalMasks(fi,c);           % returns a 1x1 cell
fprintf('%d frames, %d crystals, beta = %g\n', nF, NCrys, beta);

%% ---------- background images ----------
if ~isempty(IMGDIR)
    imgDir = IMGDIR;
else
    dd = dir(dropletDir); dd = dd([dd.isdir]);
    imgDir = '';
    for i = 1:numel(dd)
        if startsWith(dd(i).name,'FPIC','IgnoreCase',true)
            imgDir = fullfile(dropletDir, dd(i).name); break
        end
    end
end
assert(isfolder(imgDir), 'Image folder not found: %s', imgDir);
fprintf('Background: %s\n', imgDir);

% ---------- list the images, discarding the auxiliary files ----------
% The processing folders also contain background.tif, preview.mp4 and mask
% exports, which are not frames and would corrupt the frame numbering.
ext = {'*.tif','*.tiff','*.png','*.jpg'};
f = [];
for e = 1:numel(ext)
    f = [f; dir(fullfile(imgDir, ext{e}))];   %#ok<AGROW>
end
f = f(~[f.isdir]);
bad = false(numel(f),1);
for i = 1:numel(f)
    n = lower(f(i).name);
    if contains(n,'background') || contains(n,'preview') || contains(n,'mask')
        bad(i) = true;                       % auxiliary file, not a frame
    elseif isempty(regexp(n,'\d','once'))
        bad(i) = true;                       % no number: not a frame
    end
end
f = f(~bad);
assert(~isempty(f), 'No usable images in %s', imgDir);

nm = {f.name}; num = nan(size(nm));
for i = 1:numel(nm)
    t = regexp(nm{i},'\d+','match');
    if ~isempty(t), num(i) = str2double(t{end}); end
end
[~,ix] = sort(num); f = f(ix);
fprintf('%d frames in %s\n', numel(f), imgDir);

%% ---------- time axis ----------
fps = NaN; t0 = 0;
try
    specs = readSpecs(dropletDir);
    fps = specs.fps;
    if isfield(specs,'t_start') && ~isnan(specs.t_start), t0 = specs.t_start; end
catch
    warning('No specifications: titles will use frame numbers.');
end

%% ---------- which frames ----------
% radius_mean is cheap to test and does not touch the masks, so it is used
% to find the first frame in which every crystal already exists.
exists   = radius_mean > 0 & ~isnan(radius_mean);
allThere = all(exists, 2);
if isempty(FRAMES)
    first = find(allThere,1,'first');
    if isempty(first), first = find(any(exists,2),1,'first'); end
    if isempty(first), first = 1; end
    last  = min(nF, numel(f));
    FRAMES = round([first + 0.15*(last-first), last]);
end
FRAMES = min(max(FRAMES,1), min(nF,numel(f)));
fprintf('Frames: %s\n', mat2str(FRAMES));

%% ---------- common crop ----------
% Computed once from the last frame, so all panels share the same field of
% view and the crystals can be seen to grow rather than the camera to zoom.
tmp = getM(FRAMES(end),1); acc = false(size(tmp{1}));
for c = 1:NCrys
    tmp = getM(FRAMES(end),c); acc = acc | tmp{1};
end
if ZOOM && any(acc(:))
    st = regionprops(acc,'BoundingBox');
    bb = cat(1, st.BoundingBox);
    x1 = min(bb(:,1)); y1 = min(bb(:,2));
    x2 = max(bb(:,1)+bb(:,3)); y2 = max(bb(:,2)+bb(:,4));
    Rmax = max(radius_mean(FRAMES(end),:), [], 'omitnan');
    pad  = 1.6*beta*Rmax;          % leave room for the buffer circles
    cropX = [max(1,x1-pad), x2+pad];
    cropY = [max(1,y1-pad), y2+pad];
else
    cropX = []; cropY = [];
end

%% ---------- figure ----------
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
cols = lines(NCrys);

if isempty(CRYSTALS), CRYSTALS = 1:NCrys; end
fig = figure('Color','w','Position',[60 80 430*numel(FRAMES) 660]);
tl = tiledlayout(fig,1,numel(FRAMES),'TileSpacing','compact','Padding','compact');

for k = 1:numel(FRAMES)
    fi = FRAMES(k);
    ax = nexttile(tl); hold(ax,'on')

    try
        I = im2double(im2gray(imread(fullfile(imgDir, f(fi).name))));
    catch ME
        error('Cannot read %s (%s)', fullfile(imgDir, f(fi).name), ME.message);
    end
    sc = MAXPX / max(size(I));
    if sc > 1, sc = 1; end
    Is = imresize(I, sc);
    imshow(Is,'Parent',ax);

    for c = CRYSTALS
        tmp = getM(fi,c); M = tmp{1};
        if ~any(M(:)), continue; end
        Ms = imresize(M, size(Is), 'nearest');
        B = bwboundaries(Ms,'noholes');
        for b = 1:numel(B)
            % Simplified outline: the pixel staircase of a raw boundary
            % renders badly in a vector PDF and bloats the file.
            p = reducepoly(B{b}, 0.002);
            plot(ax, p(:,2), p(:,1), '-', 'Color', cols(c,:), 'LineWidth', 1.8);
        end
        cx = centroids(fi,c,1)*sc; cy = centroids(fi,c,2)*sc;
        R  = radius_mean(fi,c)*sc;
        plot(ax, cx, cy, 'o', 'MarkerSize',6, 'MarkerFaceColor',cols(c,:), ...
             'MarkerEdgeColor','k','LineWidth',0.8);
        th = linspace(0,2*pi,600);
        xc = cx + beta*R*cos(th);
        yc = cy + beta*R*sin(th);
        % Clip the circle at the image edge rather than letting it wrap.
        outside = xc < 1 | xc > size(Is,2) | yc < 1 | yc > size(Is,1);
        xc(outside) = NaN; yc(outside) = NaN;
        plot(ax, xc, yc, '--', 'Color', cols(c,:), 'LineWidth', 1.3);
    end

    if ~isempty(cropX)
        xlim(ax, [max(1,cropX(1)*sc), min(size(Is,2),cropX(2)*sc)]);
        ylim(ax, [max(1,cropY(1)*sc), min(size(Is,1),cropY(2)*sc)]);
    end
    axis(ax,'image'); set(ax,'XTick',[],'YTick',[]);

    if ~isnan(fps)
        title(ax, sprintf('$t = %.0f$ s \\quad (frame %d)', t0 + (fi-1)/fps, fi), ...
              'FontSize',15);
    else
        title(ax, sprintf('Frame %d', fi), 'FontSize',15);
    end
end

% ---------- crystal legend, built from dummy lines ----------
lg = gobjects(numel(CRYSTALS),1);
for k = 1:numel(CRYSTALS)
    lg(k) = plot(ax, NaN, NaN, '-', 'Color', cols(CRYSTALS(k),:), 'LineWidth', 2);
end
lgd = legend(lg, arrayfun(@(c) sprintf('Crystal %d',c), CRYSTALS, ...
             'UniformOutput',false), 'Interpreter','latex','FontSize',14, ...
             'Orientation','horizontal','Box','off');
lgd.Layout.Tile = 'south';

% The pause lets the layout settle before export; without it the tiled
% legend is occasionally cut off in the PDF.
drawnow; pause(0.5);
base = fullfile(outDir, ['buffer_example_' folderName]);
exportgraphics(fig, [base '.pdf'], 'ContentType','vector');
exportgraphics(fig, [base '.png'], 'Resolution', 300);
fprintf('Saved:\n  %s.pdf\n  %s.png\n', base, base);

% Restore the global interpreters, so this script does not leak its settings
% into whatever runs next in the same MATLAB session.
set(groot,'defaultTextInterpreter','remove');
set(groot,'defaultAxesTickLabelInterpreter','remove');
