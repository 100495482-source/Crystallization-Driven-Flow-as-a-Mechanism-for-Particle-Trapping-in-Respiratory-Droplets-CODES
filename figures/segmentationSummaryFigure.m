%% ========================================================================
%  SEGMENTATIONSUMMARYFIGURE  Show that the crystal masks are correct.
%  ========================================================================
%
%  PRODUCES: FIGURES ONLY. No data file is written; the masks are read from
%  the .mat that crystalGrowthAnalysis already produced.
%
%  WHAT IT DRAWS
%    A 3 x N grid for one droplet:
%      row 1  the processed frame, as the segmentation saw it
%      row 2  the crystal mask derived from it
%      row 3  the two superimposed, with the mask outlined in red
%
%    Every measured quantity in this thesis, from area to front velocity to
%    the buffer zones, descends from these masks. This figure is the evidence
%    that they follow the crystal rather than the noise, and the third row is
%    the one that carries that: a mask can look plausible on its own and
%    still be visibly offset once laid over the image.
%
%    The columns span the sequence from the first frame containing a crystal
%    to the last, so the reader sees the segmentation holding up while the
%    crystal grows by orders of magnitude in area.
%
%  ONLY EDIT THE BLOCK MARKED BELOW.
%
%  OUTPUT
%    <outDir>/masks_summary_<droplet>.pdf  (vector)
%    <outDir>/masks_summary_<droplet>.png
%
%  See also CRYSTALGROWTHANALYSIS, PIPELINESUMMARYFIGURE.

clear; close all; clc;

% ---------- WHICH DROPLET ----------
DATA_ROOT  = fullfile('C:', 'Data', 'FG');
dropletDir = fullfile(DATA_ROOT, 'TP3(1202)', 'S5');
outDir     = fullfile(DATA_ROOT, 'FotosResuls');
FRAMES     = [];      % [] = 4 automatic frames; or e.g. [40 90 150 210]
NSHOW      = 4;
% -----------------------------------

if ~isfolder(outDir), mkdir(outDir); end
resF = fullfile(dropletDir, 'RESULTS_crystals');
[~, folderName] = fileparts(dropletDir);

%% ---------- masks ----------
hit = dir(fullfile(resF, 'masks_*.mat'));
assert(~isempty(hit), 'masks_*.mat not found in %s', resF);
name = extractBetween(hit(1).name, 'masks_', '.mat'); name = name{1};
load(fullfile(resF, hit(1).name), 'masks');
nF = numel(masks);
fprintf('Masks: %d frames (%s)\n', nF, name);

%% ---------- processed image folder (FPIC*) ----------
dd = dir(dropletDir); dd = dd([dd.isdir]);
imgDir = '';
for i = 1:numel(dd)
    if startsWith(dd(i).name,'FPIC','IgnoreCase',true)
        imgDir = fullfile(dropletDir, dd(i).name); break
    end
end
assert(~isempty(imgDir), 'No FPIC* folder found in %s', dropletDir);
fprintf('Images: %s\n', imgDir);

f = dir(fullfile(imgDir,'*.tif'));
if isempty(f), f = dir(fullfile(imgDir,'*.png')); end
nm = {f.name}; num = nan(size(nm));
for i = 1:numel(nm)
    t = regexp(nm{i},'\d+','match');
    if ~isempty(t), num(i) = str2double(t{1}); end
end
[~,ix] = sort(num); f = f(ix);
nImg = numel(f);

%% ---------- time axis ----------
fps = NaN; t0 = 0;
try
    specs = readSpecs(dropletDir);
    fps = specs.fps;
    if isfield(specs,'t_start') && ~isnan(specs.t_start), t0 = specs.t_start; end
catch
    warning('No specifications: titles will use frame numbers.');
end

%% ---------- which frames are shown ----------
% Starting at the first frame that actually contains a crystal: the frames
% before nucleation are empty masks and would waste three of the columns.
if isempty(FRAMES)
    first = find(cellfun(@(m) any(m(:)), masks), 1, 'first');
    if isempty(first), first = 1; end
    last  = min(nF, nImg);
    FRAMES = round(linspace(first, last, NSHOW));
end
FRAMES = FRAMES(FRAMES >= 1 & FRAMES <= min(nF, nImg));
NSHOW  = numel(FRAMES);
fprintf('Frames shown: %s\n', mat2str(FRAMES));

%% ---------- figure ----------
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

fig = figure('Color','w','Position',[80 80 320*NSHOW 900]);
tl = tiledlayout(fig, 3, NSHOW, 'TileSpacing','compact', 'Padding','compact');

rowLabel = {'Processed frame', 'Crystal mask', 'Overlay'};

for r = 1:3
    for k = 1:NSHOW
        fi = FRAMES(k);
        ax = nexttile(tl, (r-1)*NSHOW + k);

        I = im2double(im2gray(imread(fullfile(imgDir, f(fi).name))));
        M = masks{fi};
        if ~isequal(size(M), size(I)), M = imresize(M, size(I), 'nearest'); end

        switch r
            case 1
                imshow(I, 'Parent', ax);
            case 2
                imshow(M, 'Parent', ax);
            case 3
                % Tint the mask red while keeping the underlying texture
                % visible, so a misaligned edge is obvious.
                ov = repmat(I, [1 1 3]);
                ov(:,:,1) = min(ov(:,:,1) + 0.45*double(M), 1);
                ov(:,:,2) = ov(:,:,2) .* (1 - 0.25*double(M));
                ov(:,:,3) = ov(:,:,3) .* (1 - 0.25*double(M));
                imshow(ov, 'Parent', ax); hold(ax,'on')
                if any(M(:))
                    B = bwboundaries(M,'noholes');
                    for b = 1:numel(B)
                        plot(ax, B{b}(:,2), B{b}(:,1), '-', ...
                             'Color',[0.85 0.15 0.10], 'LineWidth',1.1);
                    end
                end
        end

        if r == 1
            if ~isnan(fps)
                title(ax, sprintf('$t = %.0f$ s', t0 + (fi-1)/fps), ...
                      'Interpreter','latex', 'FontSize',15);
            else
                title(ax, sprintf('Frame %d', fi), 'Interpreter','latex','FontSize',15);
            end
        end
        if k == 1
            % Row label via the y axis: imshow hides the axes, so the axis
            % has to be made visible again with only the y colour showing.
            ylabel(ax, rowLabel{r}, 'Interpreter','latex', 'FontSize',15, ...
                   'Visible','on');
            set(ax,'YTick',[],'XTick',[],'Visible','on', ...
                   'XColor','none','YColor','k','Box','off');
            ax.YAxis.Color = 'k';
            ax.YLabel.Color = 'k';
        end
    end
end

%% ---------- export ----------
base = fullfile(outDir, ['masks_summary_' folderName]);
exportgraphics(fig, [base '.pdf'], 'ContentType','vector');
exportgraphics(fig, [base '.png'], 'Resolution', 300);
fprintf('Saved:\n  %s.pdf\n  %s.png\n', base, base);

% Restore the global interpreters so this script does not leak its settings
% into whatever runs next in the same session.
set(groot,'defaultTextInterpreter','remove');
set(groot,'defaultAxesTickLabelInterpreter','remove');
