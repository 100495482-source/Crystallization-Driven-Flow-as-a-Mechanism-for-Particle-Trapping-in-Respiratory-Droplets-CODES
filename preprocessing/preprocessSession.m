function preprocessSession(sessionDir, varargin)
%PREPROCESSSESSION  Prepare one acquisition session for analysis.
%
%   WHAT THIS DOES
%     The first stage of the pipeline. It takes a raw acquisition, whether
%     recorded as a folder of photographs or as a video file, and produces
%     the cleaned, calibrated image sequence that PIV and particle tracking
%     both run on.
%
%     It also fixes the two pieces of timing information the rest of the
%     project depends on, and writes them back into the session's
%     specifications file: the real clock time of frame 1, and the frame at
%     which the first crystal nucleates. Everything downstream reads those
%     from the specifications rather than re-deriving them, so they are
%     decided once, here.
%
%   SYNTAX
%     preprocessSession(sessionDir)
%     preprocessSession(sessionDir, 'liveEvery', 10, 'previewScale', 0.3)
%
%   NAME-VALUE PAIRS
%     'outDir'        Output folder. Asked interactively if omitted.
%     'previewScale'  Downscale factor for the preview video (default 0.4).
%     'previewFps'    Preview frame rate (default 25).
%     'showLive'      Show progress on screen (default true).
%     'liveEvery'     Refresh the live view every N frames (default 5).
%
%   STEPS
%     1.  Detect whether the session holds photographs or a video file
%     2.  Read the specifications (name, temperature, humidity, fps, ppf)
%     3.  Load the frames
%     4.  Browse the sequence and choose the first and last usable frame
%     5.  Set the time reference and the nucleation frame
%     5b. Choose the processing options
%     5c. Compute the global intensity limits
%     5d. Check the result on one frame before committing
%     6.  Name the output folder
%     7.  Process and save
%     8.  Write the metadata
%
%   WHY THE LIMITS ARE GLOBAL AND NOT PER FRAME
%     Normalising and contrast-stretching each frame independently would make
%     the overall brightness fluctuate from frame to frame. PIV works by
%     cross-correlating image texture between consecutive frames, so an
%     artificial brightness flicker degrades the correlation and corrupts the
%     velocities. Both the intensity range and the contrast limits are
%     therefore computed ONCE from a sample of the whole sequence and applied
%     identically to every frame.
%
%   OUTPUTS  (written to outDir)
%     frame_%05d.tif    the processed sequence
%     preview.mp4       downscaled preview
%     background.tif    the subtracted background, if background subtraction
%                       was used
%     preproc_meta.mat  every decision made here, for reproducibility
%
%   See also READSPECS, SUBTRACTBACKGROUND.

p = inputParser;
p.addParameter('outDir', '');
p.addParameter('previewScale', 0.4);
p.addParameter('previewFps', 25);
p.addParameter('showLive', true);
p.addParameter('liveEvery', 5);
p.parse(varargin{:});
opt = p.Results;

assert(isfolder(sessionDir), 'Folder %s does not exist', sessionDir);

%% ---------- 1. photographs or video? ----------
d = dir(sessionDir);
d = d(~ismember({d.name}, {'.','..'}));

vidExt = {'.avi','.mp4','.mov','.mj2'};
imgFolder = ''; vidFile = '';

for i = 1:numel(d)
    if d(i).isdir
        inner = dir(fullfile(sessionDir, d(i).name, '*.tif'));
        if isempty(inner), inner = dir(fullfile(sessionDir, d(i).name, '*.png')); end
        if isempty(inner), inner = dir(fullfile(sessionDir, d(i).name, '*.jpg')); end
        if ~isempty(inner), imgFolder = fullfile(sessionDir, d(i).name); end
    else
        [~,~,e] = fileparts(d(i).name);
        if ismember(lower(e), vidExt), vidFile = fullfile(sessionDir, d(i).name); end
    end
end

if ~isempty(imgFolder)
    mode = 'images';
    fprintf('Mode: IMAGES  (%s)\n', imgFolder);
elseif ~isempty(vidFile)
    mode = 'video';
    fprintf('Mode: VIDEO  (%s)\n', vidFile);
else
    error('No image folder and no video file found in %s', sessionDir);
end

if strcmp(mode,'images')
    [parentDir, srcName] = fileparts(imgFolder);
else
    [parentDir, srcName] = fileparts(vidFile);
end

%% ---------- 2. specifications: the only source of fps and scale ----------
specs = readSpecs(sessionDir);
fps = specs.fps;

%% ---------- 3. load the frames ----------
if strcmp(mode, 'video')
    v = VideoReader(vidFile);
    nNative = round(v.NumFrames);
    fprintf('Video: %d frames, %.2f native fps, %.1f s\n', ...
            nNative, v.FrameRate, v.Duration);
    if abs(v.FrameRate - fps) > 0.01
        % The container frame rate is not always what the camera actually
        % used, so the specifications value is the one that governs timing.
        fprintf(['  WARNING: video fps (%.2f) differs from the ' ...
                 'specifications (%g).\n  The specifications value is used ' ...
                 'for all timing.\n'], v.FrameRate, fps);
    end
    frames = cell(nNative,1);
    k = 0;
    while hasFrame(v)
        k = k + 1;
        A = readFrame(v);
        if size(A,3) == 3, A = rgb2gray(A); end
        frames{k} = A;
    end
    frames = frames(1:k);
else
    f = dir(fullfile(imgFolder, '*.tif'));
    if isempty(f), f = dir(fullfile(imgFolder, '*.png')); end
    if isempty(f), f = dir(fullfile(imgFolder, '*.jpg')); end
    names = sortNumeric({f.name});
    frames = cell(numel(names),1);
    for i = 1:numel(names)
        A = imread(fullfile(imgFolder, names{i}));
        if size(A,3) == 3, A = rgb2gray(A); end
        frames{i} = A;
    end
end

N = numel(frames);
fprintf('%d frames loaded.  Duration = %.1f s\n', N, (N-1)/fps);

%% ---------- 4. inspect and choose the frame range ----------
% Recordings usually start before anything is in focus and continue past the
% point where the droplet has fully dried.
browseSequence(frames, fps);

defFirst = 1;  defLast = N;
if ~isnan(specs.frame_min), defFirst = min(specs.frame_min, N); end
if ~isnan(specs.frame_max), defLast  = min(specs.frame_max, N); end

firstF = askNumber('First frame to use', defFirst, 1, N);
lastF  = askNumber('Last frame to use ', defLast, firstF, N);
sel = firstF:lastF;
fprintf('Range %d-%d  (%d frames, %.1f s)\n', ...
        firstF, lastF, numel(sel), (numel(sel)-1)/fps);

updateSpecs(specs.file, struct('frame_min', firstF, 'frame_max', lastF));
specs.frame_min = firstF;
specs.frame_max = lastF;

%% ---------- 5. time reference and nucleation ----------
%  The on-screen stopwatch and the recording do not start together, so frame
%  1 is not second 0. Anchoring the sequence to the stopwatch is what makes
%  droplets comparable with each other and with the evaporation timeline.
fprintf('\n--- Time reference ---\n');
fprintf('The stopwatch and the recording do not start at the same moment.\n');
fprintf('Enter the stopwatch time corresponding to frame 1.\n');

defT0 = 0;
if ~isnan(specs.t_start), defT0 = specs.t_start; end
tStart = askNumber('Real time (s) of the FIRST recorded frame', defT0, 0, 1e6);

fprintf('\n--- Nucleation ---\n');
fprintf('Enter the frame JUST BEFORE the first crystal appears.\n');
if askYesNo('Open the slider again to find it?')
    browseSequence(frames, fps);
end

defNuc = round(mean([firstF lastF]));
if ~isnan(specs.nucleation_frame) && specs.nucleation_frame >= firstF ...
                                 && specs.nucleation_frame <= lastF
    defNuc = specs.nucleation_frame;
end
nucF = askNumber('Frame before the first nucleation', defNuc, firstF, lastF);

tNuc = tStart + (nucF - 1)/fps;
fprintf('Nucleation: frame %d  ->  t = %.2f s (stopwatch)\n', nucF, tNuc);
fprintf('Pre-nucleation : frames %d - %d  (%d frames, %.1f s)\n', ...
        firstF, nucF, nucF-firstF+1, (nucF - firstF)/fps);
fprintf('Post-nucleation: frames %d - %d  (%d frames, %.1f s)\n', ...
        nucF+1, lastF, lastF-nucF, (lastF - nucF - 1)/fps);

updateSpecs(specs.file, struct('t_start', tStart, ...
                               'nucleation_frame', nucF, ...
                               'nucleation_time', tNuc));
specs.t_start          = tStart;
specs.nucleation_frame = nucF;
specs.nucleation_time  = tNuc;

%% ---------- 5b. processing options ----------
fprintf('\n--- Processing options ---\n');
doInvert   = askYesNo('Invert?');
doContrast = askYesNo('Stretch contrast (global, whole sequence)?');
medK       = askNumber('Median filter (0 = none, typically 3)', 0, 0, 15);

fprintf('\n--- Background subtraction ---\n');
fprintf('The background is estimated as the MEDIAN of the PRE-nucleation\n');
fprintf('frames (%d - %d), where no crystals exist yet, and subtracted from\n', ...
        firstF, nucF);
fprintf('the whole sequence. Only moving particles and crystals remain.\n');
bgSub = askYesNo('Subtract background?');

bg = [];
if bgSub
    preRange = firstF:nucF;
    if numel(preRange) < 3
        warning(['Only %d pre-nucleation frames. The background estimate ' ...
                 'will be poor; consider moving frame_min earlier or ' ...
                 'rechecking the nucleation frame.'], numel(preRange));
    end
    % Subsample rather than using every frame: the median converges fast and
    % the full stack would be unnecessarily large in memory.
    k = unique(round(linspace(1, numel(preRange), min(60, numel(preRange)))));
    A1 = double(frames{preRange(k(1))});
    tmp = zeros([size(A1) numel(k)]);
    for i = 1:numel(k)
        tmp(:,:,i) = double(frames{preRange(k(i))});
    end
    bg = median(tmp, 3);
    clear tmp
    fprintf('Background = median of %d pre-nucleation frames (of %d available).\n', ...
            numel(k), numel(preRange));

    if askYesNo('View the estimated background?')
        hb = figure('Name','Estimated background','NumberTitle','off');
        imshow(mat2gray(bg));
        title(sprintf('Median of frames %d - %d', firstF, nucF));
        input('Press ENTER to continue... ');
        if isvalid(hb), close(hb); end
    end
end

%% ---------- 5c. global intensity limits ----------
% Computed ONCE from a sample of the whole sequence and applied identically
% to every frame. See the note in the header on why this must not be per
% frame.

kg = unique(round(linspace(1, numel(sel), min(40, numel(sel)))));

gmin =  inf; gmax = -inf;
for i = 1:numel(kg)
    A = double(frames{sel(kg(i))});
    if ~isempty(bg), A = A - bg + mean(bg(:)); end
    gmin = min(gmin, min(A(:)));
    gmax = max(gmax, max(A(:)));
end
grange = [gmin gmax];
fprintf('Global intensity range: [%.1f %.1f]\n', gmin, gmax);

climits = [];
if doContrast
    % Percentile limits rather than the full range, so that a few saturated
    % pixels do not compress everything else into the middle of the scale.
    pool = [];
    for i = 1:numel(kg)
        A = double(frames{sel(kg(i))});
        if ~isempty(bg), A = A - bg + mean(bg(:)); end
        A = (A - grange(1)) / (grange(2) - grange(1));
        pool = [pool; A(:)];                                      %#ok<AGROW>
    end
    pool = min(max(pool, 0), 1);
    climits = [prctile(pool, 0.5), prctile(pool, 99.5)];
    clear pool
    fprintf('Global contrast over %d frames: [%.4f %.4f]\n', ...
            numel(kg), climits(1), climits(2));
end

%% ---------- 5d. check on one frame, nothing is written ----------
while true
    iTest = askNumber('Test frame (0 = skip the check)', ...
                      round(mean([firstF lastF])), 0, N);
    if iTest == 0, break; end

    Araw = frames{iTest};
    A = applyPipeline(Araw, bg, grange, medK, climits, doInvert);

    hc = figure('Name','Before / After  (preview only, nothing saved)', ...
                'NumberTitle','off');
    subplot(1,2,1); imshow(Araw);         title(sprintf('RAW  (frame %d)', iTest));
    subplot(1,2,2); imshow(uint8(A*255)); title('PROCESSED');

    if askYesNo('Happy with this?')
        if isvalid(hc), close(hc); end
        break;
    end
    if isvalid(hc), close(hc); end

    doInvert   = askYesNo('Invert?');
    medK       = askNumber('Median filter (0 = none)', medK, 0, 15);
    doContrast = askYesNo('Stretch contrast (global)?');
    if ~doContrast
        climits = [];
    elseif isempty(climits)
        pool = [];
        for i = 1:numel(kg)
            A = double(frames{sel(kg(i))});
            if ~isempty(bg), A = A - bg + mean(bg(:)); end
            A = (A - grange(1)) / (grange(2) - grange(1));
            pool = [pool; A(:)];                                  %#ok<AGROW>
        end
        pool = min(max(pool, 0), 1);
        climits = [prctile(pool, 0.5), prctile(pool, 99.5)];
        clear pool
    end
end

%% ---------- 6. output folder ----------
if isempty(opt.outDir)
    sug = suggestName(srcName, doInvert, doContrast, bgSub);
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
        fprintf('Cancelled.\n');
        return
    end
end
mkdir(opt.outDir);
fprintf('Output: %s\n', opt.outDir);

if ~isempty(bg)
    imwrite(uint8(mat2gray(bg)*255), fullfile(opt.outDir,'background.tif'));
end

%% ---------- 7. process and save ----------
vw = VideoWriter(fullfile(opt.outDir,'preview.mp4'),'MPEG-4');
vw.FrameRate = opt.previewFps;
vw.Quality   = 50;
open(vw);

fprintf('Processing');
hIm = [];
if opt.showLive
    hLive = figure('Name','Processing - live result','NumberTitle','off');
    axL = axes('Parent', hLive);
end

for i = 1:numel(sel)
    A = applyPipeline(frames{sel(i)}, bg, grange, medK, climits, doInvert);
    out = uint8(A*255);

    imwrite(out, fullfile(opt.outDir, sprintf('frame_%05d.tif', i)));
    writeVideo(vw, imresize(out, opt.previewScale));

    if opt.showLive && (mod(i, opt.liveEvery) == 0 || i == 1 || i == numel(sel))
        if isempty(hIm) || ~isvalid(hIm)
            hIm = imshow(out, 'Parent', axL);
        else
            set(hIm, 'CData', out);   % far faster than calling imshow again
        end
        tAbs = tStart + (sel(i)-1)/fps;
        title(axL, sprintf(['Processed  %d / %d    t = %.2f s    ' ...
              '(t-t_{nuc} = %+.2f s)'], i, numel(sel), tAbs, tAbs - tNuc));
        drawnow limitrate;
    end
    if mod(i,50) == 0, fprintf('.'); end
end
close(vw);
fprintf(' done.\n');

%% ---------- 8. metadata ----------
meta = struct('sessionDir', sessionDir, ...
              'source',     ternary(strcmp(mode,'video'), vidFile, imgFolder), ...
              'mode',       mode, ...
              'outDir',     opt.outDir, ...
              'specs',      specs, ...
              'nFramesOrig',N, ...
              'firstFrame', firstF, ...
              'lastFrame',  lastF, ...
              'tStart',     tStart, ...
              'nucFrame',   nucF, ...
              'nucTime',    tNuc, ...
              'invert',     doInvert, ...
              'contrast',   doContrast, ...
              'contrastLimits', climits, ...
              'globalRange',grange, ...
              'medFilt',    medK, ...
              'bgSubtract', bgSub, ...
              'bgFromPreNucleation', bgSub, ...
              'bgFrameRange', [firstF nucF], ...
              'date',       datestr(now, 'yyyy-mm-dd HH:MM:SS'));
save(fullfile(opt.outDir,'preproc_meta.mat'), 'meta');

fprintf('\nSaved in %s\n', opt.outDir);
fprintf('  %d images  +  preview.mp4  +  preproc_meta.mat\n', numel(sel));
if bgSub, fprintf('  + background.tif\n'); end
fprintf('  fps = %g   ppf = %.4f um/px\n', specs.fps, specs.ppf);
fprintf('  frames %d-%d   t_start = %.2f s   nucleation at frame %d (t = %.2f s)\n', ...
        firstF, lastF, tStart, nucF, tNuc);
end


%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================

function A = applyPipeline(Araw, bg, grange, medK, climits, doInvert)
%APPLYPIPELINE  The processing chain, used for both the check and the output.
%
%   Defined once and called from both places, so that what is previewed in
%   step 5d is exactly what gets written in step 7.
%
%     bg      background image (pre-nucleation median), or []
%     grange  [min max] global -> normalisation consistent across frames
%     climits [lo hi]   global -> contrast consistent across frames

A = double(Araw);

% 1. Subtract the static background. The mean level is added back so the
%    sequence does not simply go dark.
if ~isempty(bg), A = A - bg + mean(bg(:)); end

% 2. Normalise with the GLOBAL range, not this frame's range.
A = (A - grange(1)) / (grange(2) - grange(1));
A = min(max(A, 0), 1);

% 3. Median filter, for salt-and-pepper noise.
if medK > 0, A = medfilt2(A, [medK medK]); end

% 4. Contrast stretch with GLOBAL limits.
if ~isempty(climits), A = imadjust(A, climits, [0 1]); end

% 5. Inversion.
if doInvert, A = imcomplement(A); end
end


function name = suggestName(srcName, doInvert, doContrast, bgSub)
%SUGGESTNAME  Build an output folder name recording what was applied.
%   FotosS5 -> FotosProcessedInvertedBgSubS5
%   The name carries the processing history, so a folder is never ambiguous
%   about how it was produced.

tag = 'Processed';
if doInvert,   tag = [tag 'Inverted'];  end
if doContrast, tag = [tag 'Contrast'];  end
if bgSub,      tag = [tag 'BgSub'];     end

pat = '^(Fotos|Photos|Imagenes|Images|Video|Video)(.*)$';
t = regexpi(srcName, pat, 'tokens', 'once');
if ~isempty(t)
    name = [t{1} tag t{2}];
else
    name = [srcName tag];
end
end


function updateSpecs(specPath, newFields)
%UPDATESPECS  Update or add keys in the specifications file, in place.
%   Existing keys are overwritten and new ones appended; every other line is
%   left exactly as it was.

lines = readLinesUTF8(specPath);
keys  = fieldnames(newFields);

for k = 1:numel(keys)
    key = keys{k};
    val = newFields.(key);
    if isnumeric(val)
        if mod(val,1) == 0, valStr = sprintf('%d', val);
        else,               valStr = sprintf('%.4f', val); end
    else
        valStr = val;
    end
    newLine = sprintf('%s=%s', key, valStr);

    found = false;
    for i = 1:numel(lines)
        L = strtrim(lines{i});
        if isempty(L) || ~contains(L,'='), continue; end
        pp = strsplit(L,'=');
        if strcmpi(strtrim(pp{1}), key)
            lines{i} = newLine;
            found = true;
            break
        end
    end
    if ~found, lines{end+1} = newLine; end                        %#ok<AGROW>
end

fid = fopen(specPath, 'w', 'n', 'UTF-8');
for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines{i});
end
fclose(fid);
fprintf('Specifications updated: %s\n', strjoin(keys', ', '));
end


function lines = readLinesUTF8(p)
%READLINESUTF8  Read a text file as a cell array of lines, UTF-8.
fid = fopen(p, 'r', 'n', 'UTF-8');
raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);
lines = raw{1};
end


function names = sortNumeric(names)
%SORTNUMERIC  Sort filenames by the numbers they contain, not alphabetically.
%   Alphabetical order would put frame 10 before frame 2. A second number, if
%   present, is used as a fractional tiebreak.
n1 = nan(size(names));
n2 = zeros(size(names));
for i = 1:numel(names)
    t = regexp(names{i}, '\d+', 'match');
    if ~isempty(t)
        n1(i) = str2double(t{1});
        if numel(t) >= 2
            n2(i) = str2double(t{2}) / 10^numel(t{2});
        end
    end
end
if all(~isnan(n1))
    [~, idx] = sortrows([n1(:), n2(:)]);
else
    [~, idx] = sort(names);
end
names = names(idx);
fprintf('Order: %s ... %s\n', names{1}, names{end});
end


function browseSequence(frames, fps)
%BROWSESEQUENCE  Scrub through the sequence with a slider.
%   The figure is left open so the frame numbers can be noted down while
%   answering the prompts in the console.
N = numel(frames);
hf = figure('Name','Browse the sequence - note the frames you need', ...
            'NumberTitle','off');
ax = axes('Parent',hf,'Position',[0.05 0.15 0.9 0.78]);

step = [1/max(N-1,1), 10/max(N-1,1)];
sld = uicontrol('Parent',hf,'Style','slider','Units','normalized', ...
    'Position',[0.05 0.04 0.9 0.05], 'Min',1,'Max',N,'Value',1, ...
    'SliderStep',step);

drawIt = @(i) drawFrame(ax, frames{round(i)}, round(i), N, fps);
addlistener(sld, 'Value', 'PostSet', @(~,~) drawIt(get(sld,'Value')));
sld.Callback = @(s,~) drawIt(s.Value);
drawIt(1);

fprintf('\nMove the slider to inspect, then return to the console.\n');
end


function drawFrame(ax, A, i, N, fps)
% imadjust is applied for display only: faint particles are otherwise
% invisible on screen. The stored data are untouched.
imshow(imadjust(A), 'Parent', ax);
title(ax, sprintf('Frame %d / %d    t = %.2f s (from the start of the recording)', ...
      i, N, (i-1)/fps));
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
    if ~isnan(x) && x >= lo && x <= hi
        if hi > 1e5, return; end     % times: do not round
        x = round(x); return;
    end
    fprintf('  Enter a number between %g and %g.\n', lo, hi);
end
end


function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
