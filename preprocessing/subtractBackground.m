function subtractBackground(procDir, varargin)
%SUBTRACTBACKGROUND  Remove the static background from an image sequence.
%
%   WHY THIS IS NEEDED
%     The tracer particles are small and faint against a background that is
%     not uniform: dust on the optics, illumination gradients and marks on
%     the substrate all sit in the frame and never move. Those features are
%     as bright as the particles and are picked up by the particle detector,
%     so they have to be removed before tracking.
%
%     Because they are static, they can be estimated from the sequence itself
%     and subtracted. The reference image is built from the frames BEFORE
%     nucleation, where the field of view is at its cleanest.
%
%   This is equivalent to Fiji's Image Calculator > Subtract with a fixed
%   reference image. It is NOT Fiji's rolling-ball 'Subtract Background'.
%
%   It expects images that have already been inverted and contrast-stretched
%   by PREPROCESSSESSION.
%
%   SYNTAX
%     subtractBackground(procDir)
%     subtractBackground(procDir, 'method', 'min', 'bgRange', [1 120])
%
%   INPUTS
%     procDir  Folder of the already-preprocessed sequence.
%
%   NAME-VALUE PAIRS
%     'method'       'min' (default), 'median' or 'mean'. See below.
%     'bgRange'      [first last] frames used to build the background. If
%                    omitted, the function asks interactively and defaults to
%                    the frames before nucleation.
%     'maxBgFrames'  Cap on how many frames enter the estimate (default 80).
%     'medFilt'      Spatial median filter size applied after subtraction
%                    (default 0, off).
%     'outDir'       Output folder. Asked interactively if omitted.
%     'previewScale' Downscale factor for the preview video (default 0.4).
%     'previewFps'   Preview frame rate (default 25).
%     'showLive'     Show progress on screen (default true).
%     'liveEvery'    Refresh the live view every N frames (default 5).
%
%   BACKGROUND ESTIMATION METHODS
%     'min'     Temporal minimum. On INVERTED images, where the particles are
%               bright, this keeps only what is dark in EVERY frame, which is
%               exactly the static background. Recommended here.
%     'median'  Temporal median. More tolerant of illumination drift.
%     'mean'    Temporal mean. Contaminated by the particles themselves,
%               since a bright particle passing through pulls the average up.
%
%   The subtraction itself is D = frame - background, clipped to [0, 255].
%   No offset and no rescaling, so the intensities stay comparable between
%   frames and the result can be thresholded consistently.
%
%   OUTPUTS  (written to outDir)
%     frame_%05d.tif    the subtracted sequence
%     background.tif    the reference image that was subtracted
%     preview.mp4       downscaled preview of the result
%     bgsub_meta.mat    every parameter used, for reproducibility
%
%   See also PREPROCESSSESSION.

p = inputParser;
p.addParameter('method', 'min');
p.addParameter('bgRange', []);
p.addParameter('maxBgFrames', 80);
p.addParameter('medFilt', 0);
p.addParameter('outDir', '');
p.addParameter('previewScale', 0.4);
p.addParameter('previewFps', 25);
p.addParameter('showLive', true);
p.addParameter('liveEvery', 5);
p.parse(varargin{:});
opt = p.Results;

assert(isfolder(procDir), 'Folder %s does not exist', procDir);

%% ---------- 1. load the sequence ----------
f = dir(fullfile(procDir, '*.tif'));
if isempty(f), f = dir(fullfile(procDir, '*.png')); end
assert(~isempty(f), 'No images in %s', procDir);

names = sortNumeric({f.name});
N = numel(names);

fprintf('Loading %d images from %s\n', N, procDir);
frames = cell(N,1);
for i = 1:N
    A = imread(fullfile(procDir, names{i}));
    if size(A,3) == 3, A = rgb2gray(A); end
    frames{i} = A;
end

%% ---------- 2. read the preprocessing metadata, if present ----------
% Tells us whether the sequence was inverted, which decides whether 'min' is
% the right estimator, and where nucleation happens, which bounds the range
% of clean frames.
metaFile = fullfile(procDir, 'preproc_meta.mat');
prevMeta = [];
nucLocal = NaN;
if isfile(metaFile)
    S = load(metaFile);
    prevMeta = S.meta;
    fprintf('\n--- from preproc_meta.mat ---\n');
    fprintf('  inverted     : %s\n', ternary(prevMeta.invert,'yes','no'));
    fprintf('  contrast     : %s\n', ternary(prevMeta.contrast,'yes','no'));
    fprintf('  source frames: %d - %d\n', prevMeta.firstFrame, prevMeta.lastFrame);
    if isfield(prevMeta,'nucFrame') && ~isnan(prevMeta.nucFrame)
        nucLocal = prevMeta.nucFrame - prevMeta.firstFrame + 1;
        fprintf('  nucleation   : frame %d in the source -> %d here\n', ...
                prevMeta.nucFrame, nucLocal);
    end
    if ~prevMeta.invert
        fprintf(['  WARNING: the sequence is NOT inverted. With dark\n' ...
                 '  particles, ''min'' keeps the particles instead of the\n' ...
                 '  background. Use ''median'' here.\n']);
    end
end

%% ---------- 3. choose the frame range for the background ----------
if isempty(opt.bgRange)
    if askYesNo('Browse the sequence before choosing the range?')
        browseSequence(frames);
    end
    defEnd = N;
    if ~isnan(nucLocal) && nucLocal >= 1 && nucLocal <= N, defEnd = nucLocal; end
    fprintf('\nFrame range used to build the background image\n');
    fprintf('(from the first frame to just before the first nucleation)\n');
    bgIni = askNumber('First frame', 1, 1, N);
    bgFin = askNumber('Last frame ', defEnd, bgIni, N);
else
    bgIni = opt.bgRange(1);
    bgFin = min(opt.bgRange(2), N);
end
fprintf('Background built from frames %d - %d (%d frames)\n', ...
        bgIni, bgFin, bgFin-bgIni+1);

%% ---------- 4. estimate the background ----------
validMethods = {'median','min','mean'};
method = lower(opt.method);
assert(ismember(method, validMethods), 'method must be min, median or mean');

% Subsample the range evenly rather than using every frame: the estimate
% converges quickly and the full stack would not fit comfortably in memory.
bgIdx = bgIni:bgFin;
k   = unique(round(linspace(1, numel(bgIdx), min(opt.maxBgFrames, numel(bgIdx)))));
idx = bgIdx(k);

A1 = double(frames{idx(1)});
stack = zeros([size(A1) numel(idx)]);
for i = 1:numel(idx)
    stack(:,:,i) = double(frames{idx(i)});
end

bg = computeBg(stack, method);
fprintf('Background = %s of %d frames.\n', upper(method), numel(idx));

%% ---------- 5. compare the three methods visually (optional) ----------
if askYesNo('Compare min / median / mean visually?')
    hcmp = figure('Name','Background methods','NumberTitle','off');
    subplot(1,3,1); imshow(uint8(min(stack,[],3))); title('MIN');
    subplot(1,3,2); imshow(uint8(median(stack,3))); title('MEDIAN');
    subplot(1,3,3); imshow(uint8(mean(stack,3)));   title('MEAN');
    m = strtrim(input('Which one? [min/median/mean, ENTER = keep]: ','s'));
    if ~isempty(m) && ismember(lower(m), validMethods)
        method = lower(m);
        bg = computeBg(stack, method);
        fprintf('Method changed to %s.\n', upper(method));
    end
    if isvalid(hcmp), close(hcmp); end
end
clear stack

%% ---------- 6. check the result on one frame before committing ----------
% The whole sequence takes minutes to write, so the settings are validated on
% a single frame first and can be adjusted in the loop.
medK = opt.medFilt;
while true
    iTest = askNumber('Test frame (0 = skip)', round(N/2), 0, N);
    if iTest == 0, break; end

    Araw = frames{iTest};
    D = applySub(Araw, bg, medK);

    hc = figure('Name','Before / Background / After','NumberTitle','off');
    subplot(1,3,1); imshow(Araw);      title(sprintf('INPUT (frame %d)', iTest));
    subplot(1,3,2); imshow(uint8(bg)); title(sprintf('BACKGROUND (%s)', upper(method)));
    subplot(1,3,3); imshow(D);         title('SUBTRACTED');

    if askYesNo('Happy with this?')
        if isvalid(hc), close(hc); end
        break;
    end
    if isvalid(hc), close(hc); end

    medK = askNumber('Median filter after subtraction (0 = none)', medK, 0, 15);
    if askYesNo('Change the background method?')
        m = strtrim(input('  [min/median/mean]: ','s'));
        if ismember(lower(m), validMethods)
            method = lower(m);
            stack = zeros([size(A1) numel(idx)]);
            for i = 1:numel(idx), stack(:,:,i) = double(frames{idx(i)}); end
            bg = computeBg(stack, method);
            clear stack
        end
    end
end

%% ---------- 7. output folder ----------
[parentDir, srcName] = fileparts(procDir);
if isempty(opt.outDir)
    sug = [srcName 'BgSub'];
    fprintf('\nOutput folder (will be created in %s)\n', parentDir);
    nameIn = strtrim(input(sprintf('Name [%s]: ', sug), 's'));
    if isempty(nameIn), nameIn = sug; end
    opt.outDir = fullfile(parentDir, nameIn);
end

if isfolder(opt.outDir)
    fprintf('  Folder %s ALREADY EXISTS.\n', opt.outDir);
    if askYesNo('  Delete its contents and overwrite?')
        rmdir(opt.outDir, 's');
    else
        fprintf('Cancelled.\n'); return
    end
end
mkdir(opt.outDir);
imwrite(uint8(bg), fullfile(opt.outDir,'background.tif'));

%% ---------- 8. process the whole sequence ----------
vw = VideoWriter(fullfile(opt.outDir,'preview.mp4'),'MPEG-4');
vw.FrameRate = opt.previewFps; vw.Quality = 50; open(vw);

fprintf('Subtracting');
hIm = [];
if opt.showLive
    hLive = figure('Name','Subtracting background - live','NumberTitle','off');
    axL = axes('Parent', hLive);
end

for i = 1:N
    out = applySub(frames{i}, bg, medK);

    imwrite(out, fullfile(opt.outDir, sprintf('frame_%05d.tif', i)));
    writeVideo(vw, imresize(out, opt.previewScale));

    if opt.showLive && (mod(i, opt.liveEvery) == 0 || i == 1 || i == N)
        if isempty(hIm) || ~isvalid(hIm)
            hIm = imshow(out, 'Parent', axL);
        else
            set(hIm, 'CData', out);   % updating CData is far faster than imshow
        end
        title(axL, sprintf('%d / %d', i, N));
        drawnow limitrate;
    end
    if mod(i,50) == 0, fprintf('.'); end
end
close(vw);
fprintf(' done.\n');

%% ---------- 9. metadata ----------
bgMeta = struct('inputDir',  procDir, ...
                'outDir',    opt.outDir, ...
                'nFrames',   N, ...
                'method',    method, ...
                'bgRange',   [bgIni bgFin], ...
                'nBgFrames', numel(idx), ...
                'medFilt',   medK, ...
                'operation', 'frame - background, clipped to [0,255]', ...
                'prevMeta',  prevMeta, ...
                'date',      datestr(now,'yyyy-mm-dd HH:MM:SS'));
save(fullfile(opt.outDir,'bgsub_meta.mat'), 'bgMeta');

fprintf('\nSaved in %s\n', opt.outDir);
fprintf('  %d images + preview.mp4 + background.tif + bgsub_meta.mat\n', N);
fprintf('  method = %s   background from frames %d-%d\n', ...
        upper(method), bgIni, bgFin);
end


%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================

function out = applySub(Araw, bg, medK)
%APPLYSUB  Plain subtraction, clipped to [0,255].
%   Deliberately no offset and no rescaling: keeping the absolute intensity
%   scale means the same detection threshold applies to every frame.
D = double(Araw) - bg;
D = min(max(D, 0), 255);
if medK > 0, D = medfilt2(D, [medK medK]); end
out = uint8(D);
end


function bg = computeBg(stack, method)
%COMPUTEBG  Collapse a stack of frames into one background image.
switch method
    case 'min',    bg = min(stack, [], 3);
    case 'median', bg = median(stack, 3);
    case 'mean',   bg = mean(stack, 3);
end
end


function names = sortNumeric(names)
%SORTNUMERIC  Sort filenames by the numbers they contain, not alphabetically.
%   Alphabetical order would put frame 10 before frame 2. A second number, if
%   present, is used as a fractional tiebreak.
n1 = nan(size(names)); n2 = zeros(size(names));
for i = 1:numel(names)
    t = regexp(names{i}, '\d+', 'match');
    if ~isempty(t)
        n1(i) = str2double(t{1});
        if numel(t) >= 2, n2(i) = str2double(t{2}) / 10^numel(t{2}); end
    end
end
if all(~isnan(n1)), [~, idx] = sortrows([n1(:), n2(:)]);
else,               [~, idx] = sort(names); end
names = names(idx);
end


function browseSequence(frames)
%BROWSESEQUENCE  Scrub through the sequence with a slider.
N = numel(frames);
hf = figure('Name','Browse the sequence','NumberTitle','off');
ax = axes('Parent',hf,'Position',[0.05 0.15 0.9 0.78]);
step = [1/max(N-1,1), 10/max(N-1,1)];
sld = uicontrol('Parent',hf,'Style','slider','Units','normalized', ...
    'Position',[0.05 0.04 0.9 0.05],'Min',1,'Max',N,'Value',1,'SliderStep',step);
drawIt = @(i) drawFrame(ax, frames{round(i)}, round(i), N);
addlistener(sld,'Value','PostSet',@(~,~) drawIt(get(sld,'Value')));
sld.Callback = @(s,~) drawIt(s.Value);
drawIt(1);
fprintf('\nMove the slider, then return to the console.\n');
end


function drawFrame(ax, A, i, N)
imshow(A, 'Parent', ax);
title(ax, sprintf('Frame %d / %d', i, N));
drawnow limitrate;
end


function tf = askYesNo(prompt)
while true
    r = lower(strtrim(input([prompt ' (y/n): '], 's')));
    if ismember(r, {'y','yes','s','si'}), tf = true;  return; end
    if ismember(r, {'n','no'}),           tf = false; return; end
    fprintf('  Answer y or n.\n');
end
end


function x = askNumber(prompt, defaultVal, lo, hi)
while true
    r = strtrim(input(sprintf('%s [%g]: ', prompt, defaultVal), 's'));
    if isempty(r), x = defaultVal; return; end
    x = str2double(r);
    if ~isnan(x) && x >= lo && x <= hi, x = round(x); return; end
    fprintf('  Enter a number between %g and %g.\n', lo, hi);
end
end


function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
