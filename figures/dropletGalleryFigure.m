%% ========================================================================
%  DROPLETGALLERYFIGURE  Gallery of final frames across droplets.
%  ========================================================================
%
%  PRODUCES: FIGURES ONLY. No data file is written.
%
%  WHAT IT DRAWS
%    The last frame of each selected droplet, twelve per figure in a 4x3
%    grid, each labelled with its salt concentration, temperature and
%    relative humidity.
%
%    This is the figure that shows what the crystals actually look like at
%    the end of drying, across the range of conditions tested. Placed side
%    by side, the differences in crystal number and habit between
%    concentrations become visible in a way no plot conveys.
%
%  WHERE THE LABELS COME FROM
%    Temperature and humidity are read from each droplet's specifications
%    file, so they cannot drift out of step with the rest of the analysis.
%    The salt concentration is the one quantity kept in a table here, since
%    it is a property of the solution prepared rather than of the recording
%    and is not stored in the specifications.
%
%  LAYOUT NOTE
%    Unused tiles are filled with black placeholders rather than left empty.
%    Without them, a partly filled last figure stretches its remaining images
%    to fill the space and they no longer match the others in scale.
%
%  OUTPUT
%    <saveDir>/galleryred_fig%02d.pdf
%
%  See also READSPECS.

clear; clc; close all;

%% ------------------------------------------------------------------------
%  1.  SALT CONCENTRATION TABLE  (droplet number | c_s in g/L)
%      The only thing that stays hard-coded here.
%  ------------------------------------------------------------------------
cs_data = [
  1,200;  2,200;  3,150;  4,150;  5,280;  6,280;
  7,280;  8,280;  9,280; 10,280; 11,280;
 12,280; 13,280; 14,280; 15,280; 16,280; 17,280; 18,280; 19,280; 20,280;
 21,200; 22,150; 23,280; 24,280; 25,280; 26,150; 27,200; 28,200;
 29,200; 30,200; 31,280; 32,150; 33,280; 34,280; 35,280;
];

%% ------------------------------------------------------------------------
%  2.  FOLDER MAP
%      Indexed by droplet number, so dropletDirs{22} is droplet 22. Gaps are
%      left empty on purpose: those droplets were discarded during the
%      experiments and the numbering is not reused.
%  ------------------------------------------------------------------------
DATA_ROOT = fullfile('C:', 'Data', 'FG');

dropletDirs = repmat({''}, 35, 1);

dropletDirs{ 4} = fullfile(DATA_ROOT,'TP2-B(2811)','S4');
dropletDirs{ 5} = fullfile(DATA_ROOT,'TP3(1202)','S5');
dropletDirs{ 6} = fullfile(DATA_ROOT,'TP3(1202)','S6');
dropletDirs{ 7} = fullfile(DATA_ROOT,'TP4-D(0106)','S7');
dropletDirs{ 8} = fullfile(DATA_ROOT,'TP4-D(0106)','S8');
dropletDirs{ 9} = fullfile(DATA_ROOT,'TP4-D(0106)','S9');
dropletDirs{10} = fullfile(DATA_ROOT,'TP4-D(0106)','S10');
dropletDirs{12} = fullfile(DATA_ROOT,'TP5-E(0806)','S12');
dropletDirs{13} = fullfile(DATA_ROOT,'TP5-E(0806)','S13');
dropletDirs{14} = fullfile(DATA_ROOT,'TP5-E(0806)','S14');
dropletDirs{15} = fullfile(DATA_ROOT,'TP5-E(0806)','S15');
dropletDirs{18} = fullfile(DATA_ROOT,'TP5-E(0806)','S18');
dropletDirs{19} = fullfile(DATA_ROOT,'TP5-E(0806)','S19');
dropletDirs{20} = fullfile(DATA_ROOT,'TP5-E(0806)','S20');
dropletDirs{21} = fullfile(DATA_ROOT,'TP6-(2907)','S21');
dropletDirs{22} = fullfile(DATA_ROOT,'TP6-(2907)','S22');
dropletDirs{23} = fullfile(DATA_ROOT,'TP6-(2907)','S23');
dropletDirs{24} = fullfile(DATA_ROOT,'TP6-(2907)','S24');
dropletDirs{25} = fullfile(DATA_ROOT,'TP6-(2907)','S25');
dropletDirs{29} = fullfile(DATA_ROOT,'TP6-(2907)','S29');
dropletDirs{30} = fullfile(DATA_ROOT,'TP6-(2907)','S30');
dropletDirs{31} = fullfile(DATA_ROOT,'TP6-(2907)','S31');
dropletDirs{32} = fullfile(DATA_ROOT,'TP6-(2907)','S32');
dropletDirs{33} = fullfile(DATA_ROOT,'TP6-(2907)','S33');
dropletDirs{34} = fullfile(DATA_ROOT,'TP6-(2907)','S34');
dropletDirs{35} = fullfile(DATA_ROOT,'TP6-(2907)','S35');

%% ------------------------------------------------------------------------
%  3.  WHICH DROPLETS TO SHOW
%  ------------------------------------------------------------------------
dropletList = [5 6 7 9 13 18 19 22 23 25 29 30];

%% ------------------------------------------------------------------------
%  4.  LAYOUT AND EXPORT
%  ------------------------------------------------------------------------
N_PER_FIG = 12;
NCOLS     = 4;
NROWS     = 3;

FIG_W_CM  = 30;
FIG_H_CM  = 30;

saveDir   = fullfile(DATA_ROOT, 'FotosResults');

%% ------------------------------------------------------------------------
%  5.  MAIN LOOP
%  ------------------------------------------------------------------------
nDrop = numel(dropletList);
nFigs = ceil(nDrop / N_PER_FIG);

for figIdx = 1:nFigs

    iStart = (figIdx-1)*N_PER_FIG + 1;
    iEnd   = min(figIdx*N_PER_FIG, nDrop);
    batch  = dropletList(iStart:iEnd);
    nBatch = numel(batch);

    fig = figure('Units','centimeters', ...
                 'Position',[2 2 FIG_W_CM FIG_H_CM], ...
                 'Color','w');

    tl = tiledlayout(NROWS, NCOLS, ...
                     'TileSpacing', 'compact', ...
                     'Padding',     'compact');

    refImg = [];   % size of the first valid image, reused for placeholders

    for k = 1:nBatch
        dNo  = batch(k);
        dDir = dropletDirs{dNo};

        ax = nexttile;

        % ----------------------------------------------------------------
        %  A. Find the FPIC<n> folder.
        %     Matched exactly, so FPIC7trackmate is not picked up instead:
        %     that folder holds the same frames with detection overlays.
        % ----------------------------------------------------------------
        img        = [];
        imgDir     = '';
        targetName = sprintf('FPIC%d', dNo);

        if ~isempty(dDir) && isfolder(dDir)
            dd = dir(dDir);
            dd = dd([dd.isdir]);
            for i = 1:numel(dd)
                if strcmpi(dd(i).name, targetName)
                    imgDir = fullfile(dDir, dd(i).name);
                    break
                end
            end

            % ------------------------------------------------------------
            %  B. Last image, by name
            % ------------------------------------------------------------
            if ~isempty(imgDir)
                ff = [dir(fullfile(imgDir,'*.tif'));  ...
                      dir(fullfile(imgDir,'*.tiff')); ...
                      dir(fullfile(imgDir,'*.png'))];
                if ~isempty(ff)
                    [~, ord] = sort({ff.name});
                    ff  = ff(ord);
                    img = imread(fullfile(imgDir, ff(end).name));
                    if isempty(refImg), refImg = img; end
                end
            end
        end

        % ----------------------------------------------------------------
        %  C. Display
        % ----------------------------------------------------------------
        if ~isempty(img)
            imshow(img, [], 'Parent', ax);
        else
            placeholder = zeros(size(refImg,1), size(refImg,2), 'uint8');
            imshow(placeholder, [], 'Parent', ax);
            text(ax, 0.5, 0.5, sprintf('D%d not found', dNo), ...
                 'Color','r','FontSize',8, ...
                 'HorizontalAlignment','center', ...
                 'Units','normalized','Interpreter','none');
        end
        axis(ax, 'image', 'off');

        % ----------------------------------------------------------------
        %  D. Label: c_s from the table, T and RH from the specifications
        % ----------------------------------------------------------------
        cs_row = cs_data(cs_data(:,1) == dNo, :);
        cs_val = cs_row(2);

        specs = [];
        if ~isempty(dDir) && isfolder(dDir)
            try
                % evalc keeps the specifications printout out of the console:
                % it would otherwise repeat twelve times per figure.
                evalc('specs = readSpecs(dDir)');
            catch
            end
        end

        if ~isempty(specs)
            lbl = sprintf( ...
                'D%d $|$ $c_s{=}%g$\\,g/L $|$ $T{=}%.1f^\\circ$C $|$ RH${=}%.0f\\%%$', ...
                dNo, cs_val, specs.temperatura, specs.humedad);
        else
            lbl = sprintf('D%d $|$ $c_s{=}%g$\\,g/L $|$ specs not found', dNo, cs_val);
        end

        xlabel(ax, lbl, 'Interpreter','latex', 'FontSize', 8);
    end

    % --------------------------------------------------------------------
    %  Empty tiles: black placeholders keep the grid geometry, so the last
    %  images are not stretched to fill the figure.
    % --------------------------------------------------------------------
    for k = nBatch+1 : N_PER_FIG
        ax_empty = nexttile;
        if ~isempty(refImg)
            placeholder = zeros(size(refImg,1), size(refImg,2), 'uint8');
            imshow(placeholder, [], 'Parent', ax_empty);
        else
            imshow(zeros(512,512,'uint8'), 'Parent', ax_empty);
        end
        axis(ax_empty, 'off');
    end

    % --------------------------------------------------------------------
    %  Save
    % --------------------------------------------------------------------
    if ~isempty(saveDir)
        if ~isfolder(saveDir), mkdir(saveDir); end
        outFile = fullfile(saveDir, sprintf('galleryred_fig%02d.pdf', figIdx));
        exportgraphics(fig, outFile, 'Resolution', 300);
        fprintf('Saved: %s\n', outFile);
    end

end
