%% ========================================================================
%  PIPELINESUMMARYFIGURE  The image processing chain, for several droplets.
%  ========================================================================
%
%  PRODUCES: FIGURES ONLY. No data file is written.
%
%  WHAT IT DRAWS
%    A 4-column x N-row grid, one row per droplet, showing the same frame at
%    each stage of the processing chain:
%      col 1  FotosS<n>            raw frame, uncropped
%      col 2  FPICS<n>             processed
%      col 3  FPICS<n>trackmate    processed and background subtracted
%      col 4  the same, with the crystal mask overlaid in red
%
%    This is the methodological figure: it shows in one image why each
%    processing step exists. The tracer particles are barely visible in the
%    raw frame and clearly resolved by column three, which is what makes
%    particle tracking possible at all.
%
%    Showing three droplets rather than one demonstrates that the same chain
%    works across different concentrations and ambient conditions, rather
%    than having been tuned to a single favourable recording.
%
%  WHICH FRAME IS USED
%    The last usable frame of each droplet, taken from frame_max in the
%    specifications. In FotosS<n> the frame is picked BY POSITION, because
%    the filenames in that folder do not carry the frame number.
%
%  ONLY EDIT THE BLOCK MARKED BELOW.
%
%  OUTPUT
%    <outDir>/<outName>.pdf  (vector)
%    <outDir>/<outName>.png
%
%  See also SEGMENTATIONSUMMARYFIGURE, PREPROCESSSESSION, SUBTRACTBACKGROUND.

clear; close all; clc;

% ---------- WHICH DROPLETS ----------
DATA_ROOT = fullfile('C:', 'Data', 'FG');

dropletDirs = { ...
    fullfile(DATA_ROOT,'TP3(1202)','S6'), ...
    fullfile(DATA_ROOT,'TP5-E(0806)','S18'), ...
    fullfile(DATA_ROOT,'TP6-(2907)','S22')};

outDir      = fullfile(DATA_ROOT, 'FotosResults');
outName     = 'pipeline_summary';
MASK_ALPHA  = 0.40;          % transparency of the red mask tint
LEFT_MARGIN = 0.115;         % space reserved for the droplet labels
% ------------------------------------

if ~isfolder(outDir), mkdir(outDir); end

nD = numel(dropletDirs);
colLabel = {'Raw frame', 'Processed', ...
            'Processed and background subtraction', 'Crystal mask overlay'};

IMG      = cell(nD, 4);
MASKS    = cell(nD, 1);
rowLines = cell(nD, 1);

%% ---------- gather everything before drawing ----------
% Collecting first and drawing afterwards keeps the layout code free of file
% handling, and lets a missing folder be reported before any figure is opened.
for d = 1:nD
    dDir = dropletDirs{d};
    assert(isfolder(dDir), 'Folder %s does not exist', dDir);

    [~, folderName] = fileparts(dDir);
    t = regexp(folderName, '\d+', 'match');
    assert(~isempty(t), 'Cannot extract the droplet number from %s', folderName);
    dNo = str2double(t{1});

    fprintf('\n=== Droplet %d (%s) ===\n', dNo, folderName);

    % ---- specifications: the last usable frame ----
    T = NaN; RH = NaN; fMax = NaN;
    try
        evalc('specs = readSpecs(dDir);');
        T    = specs.temperatura;
        RH   = specs.humedad;
        fMax = specs.frame_max;
    catch
        warning('No specifications in %s', dDir);
    end

    % ---- locate the three folders ----
    % Several naming conventions were used across sessions, so each is tried
    % in turn rather than assuming one.
    fotosDir = findDirLike(dDir, {sprintf('FotosS%d', dNo), ...
                                  sprintf('fotosS%d', dNo), ...
                                  sprintf('Fotos%d',  dNo), ...
                                  sprintf('FOTOSS%d', dNo)});

    fpicDir  = findDirLike(dDir, {sprintf('FPICS%d', dNo), ...
                                  sprintf('FPIC%d',  dNo)});

    tmDir    = findDirLike(dDir, {sprintf('FPICS%dtrackmate', dNo), ...
                                  sprintf('FPIC%dtrackmate',  dNo)});

    assert(~isempty(fpicDir), 'FPICS%d not found in %s', dNo, dDir);
    assert(~isempty(tmDir),   'FPICS%dtrackmate not found in %s', dNo, dDir);

    % ---- columns 2 and 3: last frame of each folder ----
    fF = listImagesOrdered(fpicDir);
    fT = listImagesOrdered(tmDir);
    IMG{d,2} = im2double(im2gray(imread(fullfile(fpicDir, fF(end).name))));
    IMG{d,3} = im2double(im2gray(imread(fullfile(tmDir,   fT(end).name))));
    fprintf('  FPICS      : %s\n', fF(end).name);
    fprintf('  trackmate  : %s\n', fT(end).name);

    % ---- column 1: FotosS BY POSITION (its names carry no frame number) ----
    if ~isempty(fotosDir)
        fR = listImagesOrdered(fotosDir);
        if ~isnan(fMax) && fMax >= 1 && fMax <= numel(fR)
            k = fMax;
        else
            k = numel(fR);
            warning(['Droplet %d: frame_max=%s is out of range (%d images), ' ...
                     'using the last one.'], dNo, mat2str(fMax), numel(fR));
        end
        IMG{d,1} = im2double(im2gray(imread(fullfile(fotosDir, fR(k).name))));
        fprintf('  FotosS     : %s   (image %d of %d)\n', fR(k).name, k, numel(fR));
    else
        warning('Droplet %d: FotosS%d not found, leaving the cell empty.', dNo, dNo);
        IMG{d,1} = [];
    end

    % ---- column 4: trackmate frame with the mask in red ----
    resF = fullfile(dDir, 'RESULTS_crystals');
    hit  = dir(fullfile(resF, 'masks_*.mat'));
    if isempty(hit)
        warning('Droplet %d: no masks_*.mat, column 4 will have no mask.', dNo);
        IMG{d,4}  = repmat(IMG{d,3}, [1 1 3]);
        MASKS{d}  = [];
    else
        S = load(fullfile(resF, hit(1).name), 'masks');
        M = S.masks{end};
        base = IMG{d,3};
        if ~isequal(size(M), size(base))
            M = imresize(M, size(base), 'nearest');
        end
        ov = repmat(base, [1 1 3]);
        ov(:,:,1) = min(ov(:,:,1) + MASK_ALPHA*double(M), 1);
        ov(:,:,2) = ov(:,:,2) .* (1 - 0.6*MASK_ALPHA*double(M));
        ov(:,:,3) = ov(:,:,3) .* (1 - 0.6*MASK_ALPHA*double(M));
        IMG{d,4}  = ov;
        MASKS{d}  = M;
        fprintf('  mask       : %s\n', hit(1).name);
    end

    % ---- row label: several lines, horizontal text ----
    if ~isnan(T) && ~isnan(RH)
        rowLines{d} = { sprintf('\\textbf{Droplet %d}', dNo), ...
                        sprintf('$T = %.1f^\\circ$C', T), ...
                        sprintf('RH $= %.0f$\\%%', RH) };
    else
        rowLines{d} = { sprintf('\\textbf{Droplet %d}', dNo) };
    end
end

%% ---------- figure ----------
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

fig = figure('Color','w', 'Position', [40 40 360*4 330*nD]);
tl  = tiledlayout(fig, nD, 4, 'TileSpacing','tight', 'Padding','compact');
tl.OuterPosition = [LEFT_MARGIN 0 1-LEFT_MARGIN 1];   % room on the left

AX = gobjects(nD, 4);

for d = 1:nD
    for c = 1:4
        ax = nexttile(tl, (d-1)*4 + c);
        AX(d,c) = ax;

        if isempty(IMG{d,c})
            imshow(zeros(200,200,'uint8'), 'Parent', ax);
        else
            imshow(IMG{d,c}, [], 'Parent', ax);
        end

        % Crisp red outline on top of the tinted mask.
        if c == 4 && ~isempty(MASKS{d}) && any(MASKS{d}(:))
            hold(ax, 'on');
            B = bwboundaries(MASKS{d}, 'noholes');
            for b = 1:numel(B)
                plot(ax, B{b}(:,2), B{b}(:,1), '-', ...
                     'Color', [0.85 0.15 0.10], 'LineWidth', 1.2);
            end
            hold(ax, 'off');
        end

        % Column titles only on the first row.
        if d == 1
            title(ax, colLabel{c}, 'Interpreter','latex', 'FontSize', 12);
        end

        % Droplet label: HORIZONTAL text to the left of the row. A rotated
        % ylabel would be unreadable at this aspect ratio.
        if c == 1
            text(ax, -0.05, 0.5, rowLines{d}, ...
                 'Units','normalized', ...
                 'HorizontalAlignment','right', ...
                 'VerticalAlignment','middle', ...
                 'Interpreter','latex', ...
                 'FontSize', 13, ...
                 'Clipping','off');
        end
    end
end

%% ---------- thin separators between rows ----------
% Drawn in figure coordinates from the settled axes positions, which is why
% drawnow is needed first.
drawnow;
fpix = getpixelposition(fig);
Wf = fpix(3);  Hf = fpix(4);

for d = 1:nD-1
    pTop = getpixelposition(AX(d,1),   true);   % row above
    pBot = getpixelposition(AX(d+1,1), true);   % row below
    pRig = getpixelposition(AX(d,4),   true);   % last column

    yMid = (pTop(2) + (pBot(2) + pBot(4))) / 2;

    x1 = 0.015;
    x2 = (pRig(1) + pRig(3)) / Wf;

    annotation(fig, 'line', [x1 x2], [yMid/Hf yMid/Hf], ...
               'Color', [0.78 0.78 0.78], 'LineWidth', 0.4);
end

%% ---------- export ----------
base = fullfile(outDir, outName);
exportgraphics(fig, [base '.pdf'], 'ContentType','vector');
exportgraphics(fig, [base '.png'], 'Resolution', 300);
fprintf('\nSaved:\n  %s.pdf\n  %s.png\n', base, base);

set(groot, 'defaultTextInterpreter', 'remove');
set(groot, 'defaultAxesTickLabelInterpreter', 'remove');


%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================

function p = findDirExact(parentDir, targetName)
%FINDDIREXACT  Path of the subfolder whose name matches EXACTLY.
%   Exact matching matters here: a prefix match on 'FPICS7' would also
%   return 'FPICS7trackmate', which is a different processing stage.
p  = '';
dd = dir(parentDir);
dd = dd([dd.isdir]);
for i = 1:numel(dd)
    if strcmpi(dd(i).name, targetName)
        p = fullfile(parentDir, dd(i).name);
        return
    end
end
end


function p = findDirLike(parentDir, candidates)
%FINDDIRLIKE  Try several possible names and return the first that exists.
p = '';
for k = 1:numel(candidates)
    p = findDirExact(parentDir, candidates{k});
    if ~isempty(p), return; end
end
end


function f = listImagesOrdered(imgDir)
%LISTIMAGESORDERED  List images in NATURAL order, as the file explorer shows.
%
%   Does not assume the filename carries a frame number, because the raw
%   FotosS<n> folders do not. Numeric runs inside each name are zero-padded
%   before sorting, so frame 2 comes before frame 10.

exts = {'*.tif','*.tiff','*.png','*.jpg','*.jpeg','*.bmp'};
f = [];
for k = 1:numel(exts)
    f = dir(fullfile(imgDir, exts{k}));
    if ~isempty(f), break; end
end
assert(~isempty(f), 'No images in %s', imgDir);

nm  = {f.name};
key = cell(size(nm));
for i = 1:numel(nm)
    parts = regexp(nm{i}, '\d+|\D+', 'match');
    s = '';
    for p = 1:numel(parts)
        if all(isstrprop(parts{p}, 'digit'))
            s = [s sprintf('%012d', str2double(parts{p}))]; %#ok<AGROW>
        else
            s = [s parts{p}]; %#ok<AGROW>
        end
    end
    key{i} = s;
end
[~, ix] = sort(key);
f = f(ix);
end
