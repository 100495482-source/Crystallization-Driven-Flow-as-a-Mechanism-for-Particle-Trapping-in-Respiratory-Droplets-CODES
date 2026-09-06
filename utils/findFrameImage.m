function p = findFrameImage(dDir, frIdx)
%FINDFRAMEIMAGE  Path of one raw frame of a droplet, for figure backgrounds.
%
%   Locates the droplet's raw image folder and returns the path of frame
%   frIdx, sorted numerically. Returns '' if nothing is found, so the caller
%   can draw the figure without a background rather than failing.
%
%   SYNTAX
%     p = findFrameImage(dDir, frIdx)
%
%   The folder must start with FPIC and must NOT contain 'track': the
%   TrackMate working folder (FPICTrackmate) holds the same frames with the
%   detection overlays burnt in, which is not what a figure background wants.
%
%   frIdx is clamped to the available range, so asking for the last frame of
%   a shorter sequence returns its actual last frame instead of erroring.

p = '';
dd = dir(dDir); dd = dd([dd.isdir]);
imgDir = '';
for i = 1:numel(dd)
    if ~isempty(regexpi(dd(i).name, '^FPIC', 'once')) && ...
       isempty(regexpi(dd(i).name, 'track', 'once'))
        imgDir = fullfile(dDir, dd(i).name); break
    end
end
if isempty(imgDir), return; end

f = dir(fullfile(imgDir, '*.tif'));
if isempty(f), f = dir(fullfile(imgDir, '*.png')); end
if isempty(f), return; end

nmf = {f.name}; num = nan(size(nmf));
for i = 1:numel(nmf)
    tk = regexp(nmf{i}, '\d+', 'match');
    if ~isempty(tk), num(i) = str2double(tk{end}); end
end
if all(~isnan(num)), [~,ix] = sort(num); else, [~,ix] = sort(nmf); end
f = f(ix);
frIdx = min(max(frIdx, 1), numel(f));
p = fullfile(imgDir, f(frIdx).name);
end
