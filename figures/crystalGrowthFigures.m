%% ========================================================================
%  CRYSTALGROWTHFIGURES  Crystal growth results, global and per crystal.
%  ========================================================================
%
%  PRODUCES: FIGURES AND TABLES. Nothing is recomputed from images; every
%  quantity comes from the .mat files crystalGrowthAnalysis already wrote.
%
%  summary_crystal_growth.csv, written here, is the starting point of
%  interDropletFigures. Run this first.
%
%  WHAT IT DRAWS
%    FIGURE 1  (PLOT_GLOBAL_GROWTH) three separate figures:
%       1a  total crystallised area vs time
%       1b  global growth rate dA/dt vs time
%       1c  equivalent radial front velocity dR_eq/dt vs time
%      Each nucleation is marked with a numbered circle, so the jumps in the
%      curves can be attributed to a specific new crystal rather than read as
%      noise.
%
%    FIGURE 2  (PLOT_INDIVIDUAL) one figure per droplet, two rows:
%       row 1  N_FRAMES_SHOW frames with the mask in red, plus one CLEAN
%              frame with each crystal tinted in its own curve colour
%       row 2  individual area  |  radial growth velocity
%      The clean panel is the key to reading the figure: it tells you which
%      crystal in the image corresponds to which curve below.
%
%    LOG-LOG      crystal area vs time since ITS nucleation, with reference
%                 slopes. The exponent n in A ~ (t-t_nuc)^n distinguishes
%                 diffusive growth (n around 1/2) from linear (n around 1).
%
%    NORMALISED   every crystal aligned at its own nucleation and scaled to
%                 its own maximum, which compares the SHAPE of the growth
%                 curves independently of how large each crystal became.
%
%    RVEQ_FRAC    front velocity against crystallised fraction, i.e. how the
%                 front slows as the available salt is consumed.
%
%    SUMMARY      per-droplet statistics: console, CSV and a PNG table.
%
%  TIME AXIS
%    frameTimes and nucList stored inside the .mat files may be stale if
%    t_start or the nucleation frames were corrected afterwards, so both are
%    always rebuilt from the specifications. The .mat is only a fallback.
%
%    This script READS the specifications and never writes to them.
%
%  READS  (output of crystalGrowthAnalysis)
%    crystal_areas_<n>.mat       -> areasUm, nucList, frameTimes
%    growth_individual_<n>.mat   -> areaCrystal_um2
%    radial_velocity_<n>.mat     -> vRadial
%    individual_crystals_<n>.mat -> centroids
%    masks_<n>.mat               -> masks
%    images from the first FPIC* subfolder of the session directory
%
%  ONLY EDIT THE SECTIONS MARKED <<<
%
%  See also CRYSTALGROWTHANALYSIS, INTERDROPLETFIGURES.

clear; close all; clc;

%% <<<  0. DATA ROOT  >>>
DATA_ROOT = fullfile('C:', 'Data', 'FG');

% Helper to keep the lists below short and free of absolute paths.
D = @(session, droplet) fullfile(DATA_ROOT, session, droplet);

%% <<<  1. DROPLETS - GLOBAL FIGURE  >>>
dropletDirs_global = {
    D('TP3(1202)','S5')
    D('TP3(1202)','S6')
    D('TP5-E(0806)','S19')
    D('TP6-(2907)','S22')
    D('TP6-(2907)','S23')
    D('TP6-(2907)','S29')
    D('TP6-(2907)','S30')
    D('TP6-(2907)','S32')
    };

%% <<<  2. DROPLETS - INDIVIDUAL FIGURE  >>>
dropletDirs_indiv = {
    D('TP6-(2907)','S30')
    };

%% <<<  3. ADDITIONAL ANALYSES - FLAGS AND DROPLETS  >>>
PLOT_GLOBAL_GROWTH = false;  % figures 1a, 1b, 1c
PLOT_INDIVIDUAL    = false;  % one figure per droplet
PLOT_LOGLOG        = true;   % log-log: crystal area vs (t - t_nuc)
PLOT_NORM_CURVES   = true;   % normalised curves, aligned at nucleation
PLOT_RVEQ_FRAC     = false;  % dR_eq/dt vs crystallised fraction
MAKE_SUMMARY_TABLE = true;   % summary table: console + CSV + PNG

%  Droplets for each analysis. Leave a list as {} to reuse the droplets of
%  dropletDirs_global instead.
dropletDirs_loglog = {
    D('TP2-B(2811)','S4')
    D('TP3(1202)','S5')
    D('TP3(1202)','S6')
    D('TP5-E(0806)','S12')
    D('TP5-E(0806)','S13')
    D('TP5-E(0806)','S18')
    D('TP5-E(0806)','S20')
    D('TP6-(2907)','S21')
    D('TP6-(2907)','S22')
    D('TP6-(2907)','S24')
    D('TP6-(2907)','S30')
    };

dropletDirs_norm = dropletDirs_loglog;   % same set

dropletDirs_rveq = {
    D('TP2-B(2811)','S4')
    D('TP3(1202)','S5')
    D('TP3(1202)','S6')
    D('TP4-D(0106)','S7')
    D('TP4-D(0106)','S8')
    D('TP4-D(0106)','S9')
    D('TP4-D(0106)','S10')
    };

% The summary table covers ALL 26 droplets: it feeds interDropletFigures,
% which needs every one of them.
dropletDirs_table = {
    D('TP2-B(2811)','S4')
    D('TP3(1202)','S5');    D('TP3(1202)','S6')
    D('TP4-D(0106)','S7');  D('TP4-D(0106)','S8')
    D('TP4-D(0106)','S9');  D('TP4-D(0106)','S10')
    D('TP5-E(0806)','S12'); D('TP5-E(0806)','S13')
    D('TP5-E(0806)','S14'); D('TP5-E(0806)','S15')
    D('TP5-E(0806)','S18'); D('TP5-E(0806)','S19')
    D('TP5-E(0806)','S20')
    D('TP6-(2907)','S21');  D('TP6-(2907)','S22')
    D('TP6-(2907)','S23');  D('TP6-(2907)','S24')
    D('TP6-(2907)','S25');  D('TP6-(2907)','S29')
    D('TP6-(2907)','S30');  D('TP6-(2907)','S31')
    D('TP6-(2907)','S32');  D('TP6-(2907)','S33')
    D('TP6-(2907)','S34');  D('TP6-(2907)','S35')
    };

%% <<<  4. OUTPUT FOLDER  >>>
outDir = fullfile(DATA_ROOT, 'RESULTS_CRYSTALGROWTH');

%% ------------------------------------------------------------------------
%  PARAMETERS  (style only, apart from the two that must match the pipeline)
%  ------------------------------------------------------------------------
FONT_SZ       = 14;   % axis labels and titles
TICK_SZ       = 13;   % tick labels
LINE_W        = 0.6;  % data line width
N_FRAMES_SHOW = 4;    % masked frames (+1 clean frame = N_FRAMES_SHOW+1 panels)
ALPHA_MASK    = 0.40; % transparency of the red mask overlay
N_FINITE_DIFF = 4;    % finite-difference gap; MUST match crystalGrowthAnalysis
SMOOTH_WINDOW = 10;   % smoothing window;      MUST match crystalGrowthAnalysis
TITLE_FS      = 20;   % panel titles in the individual figure

% One colour per droplet. Extended with lines() if there are more droplets
% than rows here.
CMAP = [
    0.122  0.471  0.706;   0.890  0.102  0.110
    0.173  0.627  0.173;   0.580  0.404  0.741
    1.000  0.498  0.055;   0.694  0.349  0.157
    0.969  0.506  0.749;   0.502  0.502  0.502
    0.000  0.588  0.588;   0.420  0.557  0.137
    0.027  0.212  0.490;   0.502  0.000  0.000
    1.000  0.702  0.000;   0.388  0.604  0.804
    0.000  0.392  0.392;   0.698  0.000  0.698
    0.545  0.000  0.275;   0.914  0.588  0.478
    0.251  0.251  0.251;   0.855  0.647  0.125
    0.000  0.482  0.655;   0.500  0.784  0.196
    0.294  0.000  0.510;   1.000  0.388  0.278
    0.180  0.545  0.341;   0.600  0.196  0.800
    0.800  0.800  0.200;   0.400  0.200  0.600
];

if ~isfolder(outDir), mkdir(outDir); end

%% ========================================================================
%  LOAD THE GLOBAL SET
%  ========================================================================
%  Loaded once here through loadGStruct, the same function the other analyses
%  use, so there is a single definition of how a droplet is read. Only loaded
%  when something actually needs it.

NEED_G = PLOT_GLOBAL_GROWTH || ...
         (PLOT_LOGLOG        && isempty(dropletDirs_loglog))  || ...
         (PLOT_NORM_CURVES   && isempty(dropletDirs_norm))    || ...
         (PLOT_RVEQ_FRAC     && isempty(dropletDirs_rveq))    || ...
         (MAKE_SUMMARY_TABLE && isempty(dropletDirs_table));

G = [];
if NEED_G
    G = loadGStruct(dropletDirs_global, N_FINITE_DIFF, SMOOTH_WINDOW);
    if isempty(G)
        warning('No droplet loaded - the global figures will be skipped.');
    end
end

%% ========================================================================
%  FIGURE 1 - GLOBAL GROWTH
%  ========================================================================
%  Only t >= 0 is shown: before the first nucleation there is no crystal, so
%  plotting that stretch would only add a flat run of zeros.

if PLOT_GLOBAL_GROWTH && ~isempty(G)

    nG   = numel(G);
    cols = CMAP(1:min(nG, size(CMAP,1)), :);
    if nG > size(CMAP,1), cols = [cols; lines(nG - size(CMAP,1))]; end

    t_end = max(arrayfun(@(k) G(k).t(end), 1:nG));

    % ---------- 1a: total crystallised area ----------
    figA = figure('Color','w','Units','inches','Position',[1 1 14 5]);
    axA  = axes(figA); hold(axA,'on');
    plotGlobalCurve(axA, G, cols, 'A', 'nuc_A', LINE_W);
    xlim(axA, [0, t_end]);
    xlabel(axA, 'Time from the first nucleation (s)', ...
           'Interpreter','latex','FontSize',FONT_SZ);
    ylabel(axA, 'Total area ($\mu$m$^2$)', ...
           'Interpreter','latex','FontSize',FONT_SZ);
    title(axA, 'Total crystallized area', ...
          'Interpreter','latex','FontSize',FONT_SZ+1);
    legend(axA, 'Interpreter','latex','Location','northwest', ...
           'FontSize',9,'NumColumns',2);
    set(axA,'TickLabelInterpreter','latex','FontSize',TICK_SZ); box(axA,'on');
    exportgraphics(figA, fullfile(outDir,'global_crystal_area.png'), 'Resolution',300);
    fprintf('Saved: global_crystal_area.png\n');

    % ---------- 1b: global growth rate ----------
    figB = figure('Color','w','Units','inches','Position',[7.5 1 14 5]);
    axB  = axes(figB); hold(axB,'on');
    plotGlobalCurve(axB, G, cols, 'dA', 'nuc_dA', LINE_W);
    yline(axB, 0, '--', 'Color',[0.5 0.5 0.5], ...
          'LineWidth',0.9,'HandleVisibility','off');
    xlim(axB, [0, t_end]);
    xlabel(axB, 'Time from the first nucleation (s)', ...
           'Interpreter','latex','FontSize',FONT_SZ);
    ylabel(axB, '$\mathrm{d}A/\mathrm{d}t$ ($\mu$m$^2$ s$^{-1}$)', ...
           'Interpreter','latex','FontSize',FONT_SZ);
    title(axB, 'Global growth rate', ...
          'Interpreter','latex','FontSize',FONT_SZ+1);
    legend(axB, 'Interpreter','latex','Location','northeast', ...
           'FontSize',9,'NumColumns',2);
    set(axB,'TickLabelInterpreter','latex','FontSize',TICK_SZ); box(axB,'on');
    exportgraphics(figB, fullfile(outDir,'global_crystal_rate.png'), 'Resolution',300);
    fprintf('Saved: global_crystal_rate.png\n');

    % ---------- 1c: equivalent radial front velocity ----------
    %  R_eq = sqrt(A_total/pi), so this is a length per unit time and is
    %  therefore comparable with the measured flow velocities.
    figC = figure('Color','w','Units','inches','Position',[14 1 14 5]);
    axC2 = axes(figC); hold(axC2,'on');
    plotGlobalCurve(axC2, G, cols, 'dR_eq', 'nuc_dReq', LINE_W);
    yline(axC2, 0, '--', 'Color',[0.5 0.5 0.5], ...
          'LineWidth',0.9,'HandleVisibility','off');
    xlim(axC2, [0, t_end]);
    xlabel(axC2, 'Time from the first nucleation (s)', ...
           'Interpreter','latex','FontSize',FONT_SZ);
    ylabel(axC2, '$\dot{R}_\mathrm{eq}$ ($\mu$m s$^{-1}$)', ...
           'Interpreter','latex','FontSize',FONT_SZ);
    title(axC2, 'Equivalent radial growth velocity', ...
          'Interpreter','latex','FontSize',FONT_SZ+1);
    text(axC2, 0.97, 0.97, ...
         ['$\dot{R}_\mathrm{eq} = \mathrm{d}R_\mathrm{eq}/\mathrm{d}t,' ...
          '\quad R_\mathrm{eq}=\sqrt{A_\mathrm{tot}/\pi}$'], ...
         'Units','normalized','HorizontalAlignment','right', ...
         'VerticalAlignment','top','Interpreter','latex', ...
         'FontSize',8,'Color',[0.45 0.45 0.45]);
    legend(axC2, 'Interpreter','latex','Location','northeast', ...
           'FontSize',9,'NumColumns',2);
    set(axC2,'TickLabelInterpreter','latex','FontSize',TICK_SZ); box(axC2,'on');
    exportgraphics(figC, fullfile(outDir,'global_radial_velocity.png'), 'Resolution',300);
    fprintf('Saved: global_radial_velocity.png\n');
end

%% ========================================================================
%  FIGURE 2 - INDIVIDUAL GROWTH, one figure per droplet
%  ========================================================================
%  row 1 : [masked frame 1] ... [masked frame N] [clean frame, crystals tinted]
%  row 2 : [ area per crystal (left) ]           [ radial velocity (right) ]

if PLOT_INDIVIDUAL

    for i = 1:numel(dropletDirs_indiv)
        dDir = dropletDirs_indiv{i};
        [~, folderName] = fileparts(dDir);
        if ~isfolder(dDir)
            warning('Folder not found: %s', dDir); continue
        end

        resF = fullfile(dDir, 'RESULTS_crystals');

        % The droplet name inside the filenames is not always the folder name.
        name = folderName;
        hit  = dir(fullfile(resF, 'crystal_areas_*.mat'));
        if isempty(hit)
            warning('No crystal_areas_*.mat in %s', resF); continue
        end
        tk = extractBetween(hit(1).name, 'crystal_areas_', '.mat');
        if ~isempty(tk), name = tk{1}; end

        fprintf('Individual: loading %s ...\n', name);

        try
            ca = load(fullfile(resF, ['crystal_areas_'       name '.mat']), ...
                      'nucList','frameTimes');
            gi = load(fullfile(resF, ['growth_individual_'   name '.mat']), ...
                      'areaCrystal_um2');
            rv = load(fullfile(resF, ['radial_velocity_'     name '.mat']), ...
                      'vRadial');
            ic = load(fullfile(resF, ['individual_crystals_' name '.mat']), ...
                      'centroids');
            Sm = load(fullfile(resF, ['masks_'               name '.mat']), ...
                      'masks');
        catch ME
            warning('Error loading %s:\n  %s', name, ME.message);
            continue
        end

        nucList  = ca.nucList(:);
        frameT   = ca.frameTimes(:)';
        areaCrys = double(gi.areaCrystal_um2);   % [NF x NCrys]
        vRad     = double(rv.vRadial);           % [NF x NCrys]
        cents    = ic.centroids;                 % [NF x NCrys x 2]
        masks    = Sm.masks;                     % {NF x 1}

        NF      = numel(frameT);
        nucList = min(max(round(nucList), 1), NF);

        % ---- TIME AND NUCLEATION ALWAYS FROM THE SPECIFICATIONS ----
        [frameT, nucList] = timeAxisFromSpecs(dDir, NF, frameT, nucList, name);

        NCrys = numel(nucList);
        t0    = frameT(nucList(1));
        t_rel = frameT - t0;

        % NaN before each crystal's own nucleation
        areaNaN = areaCrys;
        vRadNaN = vRad;
        for ci = 1:NCrys
            f0 = nucList(ci);
            if f0 > 1
                areaNaN(1:f0-1, ci) = NaN;
                vRadNaN(1:f0-1, ci) = NaN;
            end
        end

        % ---- image folder (first FPIC*) ----
        imgFolders = dir(fullfile(dDir, 'FPIC*'));
        imgFolders = imgFolders([imgFolders.isdir]);
        hasImgs = false;  flist = [];  nImgs = 0;  imgDir = '';

        if ~isempty(imgFolders)
            imgDir = fullfile(dDir, imgFolders(1).name);
            fl = dir(fullfile(imgDir,'*.tif'));
            if isempty(fl), fl = dir(fullfile(imgDir,'*.TIF'));  end
            if isempty(fl), fl = dir(fullfile(imgDir,'*.tiff')); end
            if isempty(fl), fl = dir(fullfile(imgDir,'*.png'));  end

            if ~isempty(fl)
                nums = nan(numel(fl),1);
                for k = 1:numel(fl)
                    tok = regexp(fl(k).name,'\d+','match','once');
                    if ~isempty(tok), nums(k) = str2double(tok); end
                end
                if ~all(isnan(nums)), [~,ix] = sort(nums);
                else,                 [~,ix] = sort({fl.name}); end
                flist   = fl(ix);
                nImgs   = numel(flist);
                hasImgs = true;
            end
        end
        if ~hasImgs
            warning('No FPIC* images found in %s', dDir);
        end

        % ---- which frames are shown ----
        % Starting at the first nucleation: the earlier frames have no mask.
        nFrameMax   = max(1, min(NF, nImgs));
        f_start     = max(1, nucList(1));
        disp_f_mask = unique(round(linspace(f_start, nFrameMax, N_FRAMES_SHOW)));
        disp_f_mask(disp_f_mask < 1)         = [];
        disp_f_mask(disp_f_mask > nFrameMax) = [];
        nMask = numel(disp_f_mask);
        if nMask < 1, nMask = 1; disp_f_mask = f_start; end
        f_clean = nFrameMax;
        nCols   = max(nMask + 1, 2);

        cry_clr = lines(max(NCrys, 1));

        % ---- figure with manual axes positioning ----
        % tiledlayout cannot mix an image strip and two wide plots at these
        % proportions, so the axes are placed by hand.
        figW = max(960, nCols * 230);
        figH = 920;
        fig2 = figure('Color','w','Position',[50 50 figW figH], ...
                      'Name',['Individual: ' name]);

        lm      = 0.08;                          % left margin, for the ylabels
        rm      = 0.025;                         % right margin
        tw      = 1 - lm - rm;                   % usable width
        imgGap  = 0.008;                         % gap between images
        eiw     = (tw - (nCols-1)*imgGap) / nCols;
        imgBot  = 0.64;                          % bottom of the image strip
        imgH_n  = 0.23;                          % height of the image strip
        pBot    = 0.08;                          % bottom of the plots
        pTop    = imgBot - 0.06;
        pH_n    = pTop - pBot;
        pMidGap = 0.1;
        epw     = (tw - pMidGap) / 2;

        annotation(fig2,'textbox',[0, 0.935, 1, 0.055], ...
                   'String',[strrep(name,'_','\_') ' --- Individual crystals analysis'], ...
                   'Interpreter','latex','FontSize',FONT_SZ+2,'FontWeight','bold', ...
                   'EdgeColor','none','HorizontalAlignment','center', ...
                   'VerticalAlignment','middle');
        annotation(fig2,'textbox',[lm, imgBot+imgH_n+0.020, tw, 0.026], ...
                   'String','Crystal masks overlay', ...
                   'Interpreter','latex','FontSize',14,'EdgeColor','none', ...
                   'HorizontalAlignment','center','VerticalAlignment','middle');

        % ---- top strip: frames with the mask in red ----
        for k = 1:nMask
            fi   = disp_f_mask(k);
            xpos = lm + (k-1)*(eiw + imgGap);
            axF  = axes('Parent',fig2, ...                          %#ok<LAXES>
                        'Position',[xpos, imgBot, eiw, imgH_n]);

            im = zeros(100,100);
            if hasImgs && fi <= nImgs
                try
                    raw = imread(fullfile(imgDir, flist(fi).name));
                    if size(raw,3)==3, im = im2double(rgb2gray(raw));
                    else,              im = im2double(raw(:,:,1)); end
                    im = imadjust(im);   % display only, the data are untouched
                catch
                    warning('Could not read frame %d of %s', fi, name);
                end
            end
            imshow(im,'Parent',axF);
            hold(axF,'on');

            if fi <= numel(masks) && any(masks{fi}(:))
                mk    = single(masks{fi});
                [H,W] = size(mk);
                ov    = cat(3, ones(H,W,'single'), ...
                               zeros(H,W,'single'), zeros(H,W,'single'));
                hov   = image(axF, ov);
                set(hov,'AlphaData', mk * ALPHA_MASK);
            end
            text(axF, 0.5, 0.97, sprintf('$t = %.1f$ s', frameT(fi)), ...
                 'Units','normalized','HorizontalAlignment','center', ...
                 'VerticalAlignment','top','Interpreter','latex','FontSize',9, ...
                 'Color','w','BackgroundColor',[0 0 0 0.45]);
            set(axF,'XTick',[],'YTick',[]);
        end

        % ---- last panel: clean frame with each crystal in its curve colour ----
        %  This is what links the images to the plots below. Each crystal is
        %  identified by finding the connected component that contains its
        %  centroid, then tinting that component.
        xpos_cl = lm + (nCols-1)*(eiw + imgGap);
        axCl = axes('Parent',fig2, ...                              %#ok<LAXES>
                    'Position',[xpos_cl, imgBot, eiw, imgH_n]);

        im = zeros(100,100);
        if hasImgs && f_clean <= nImgs
            try
                raw = imread(fullfile(imgDir, flist(f_clean).name));
                if size(raw,3)==3, im = im2double(rgb2gray(raw));
                else,              im = im2double(raw(:,:,1)); end
                im = imadjust(im);
            catch
                warning('Could not read the clean frame of %s', name);
            end
        end
        imshow(im,'Parent',axCl);
        hold(axCl,'on');

        if f_clean <= numel(masks) && any(masks{f_clean}(:))
            combined_mask = masks{f_clean};
            labeled_mask  = bwlabel(combined_mask);
            [H, W]        = size(combined_mask);

            for ci = 1:NCrys
                if f_clean <= size(cents,1)
                    cx = cents(f_clean, ci, 1);
                    cy = cents(f_clean, ci, 2);
                    if ~isnan(cx) && ~isnan(cy)
                        rx = max(1, min(W, round(cx)));
                        ry = max(1, min(H, round(cy)));
                        lbl_id = labeled_mask(ry, rx);
                        if lbl_id > 0
                            crys_mk = single(labeled_mask == lbl_id);
                            cc = cry_clr(ci,:);
                            ov = cat(3, ones(H,W,'single')*cc(1), ...
                                        ones(H,W,'single')*cc(2), ...
                                        ones(H,W,'single')*cc(3));
                            hov = image(axCl, ov);
                            set(hov,'AlphaData', crys_mk * 0.65);
                        end
                    end
                end
            end
        end

        text(axCl, 0.5, 0.97, 'Crystal identification', ...
             'Units','normalized','HorizontalAlignment','center', ...
             'VerticalAlignment','top','Interpreter','latex','FontSize',9, ...
             'Color','w','BackgroundColor',[0 0 0 0.6]);
        set(axCl,'XTick',[],'YTick',[]);

        % ---- left plot: individual crystal area ----
        axC = axes('Parent',fig2,'Position',[lm, pBot, epw, pH_n]); %#ok<LAXES>
        hold(axC,'on');  box(axC,'on');
        for ci = 1:NCrys
            plot(axC, t_rel, areaNaN(:,ci)', '-', ...
                 'Color', cry_clr(ci,:), 'LineWidth', 1.3, ...
                 'DisplayName', sprintf('Crystal $%d$', ci));
        end
        xlabel(axC,'Time (s)','Interpreter','latex','FontSize',FONT_SZ);
        ylabel(axC,'Crystal area ($\mu$m$^2$)','Interpreter','latex','FontSize',FONT_SZ);
        title(axC,'(a) Individual crystal area','Interpreter','latex','FontSize',TITLE_FS);
        legend(axC,'Interpreter','latex','NumColumns',min(NCrys,4), ...
               'FontSize',8,'Location','northwest');
        set(axC,'TickLabelInterpreter','latex','FontSize',TICK_SZ);

        % ---- right plot: radial growth velocity ----
        axD = axes('Parent',fig2, ...                               %#ok<LAXES>
                   'Position',[lm+epw+pMidGap, pBot, epw, pH_n]);
        hold(axD,'on');  box(axD,'on');
        for ci = 1:NCrys
            plot(axD, t_rel, vRadNaN(:,ci)', '-', ...
                 'Color', cry_clr(ci,:), 'LineWidth', 1.3, ...
                 'DisplayName', sprintf('Crystal $%d$', ci));
        end
        yline(axD,0,'--','Color',[0.5 0.5 0.5],'LineWidth',0.9,'HandleVisibility','off');
        xlabel(axD,'Time (s)','Interpreter','latex','FontSize',FONT_SZ);
        ylabel(axD,'Radial growth velocity $\dot{R}_i$ ($\mu$m s$^{-1}$)', ...
               'Interpreter','latex','FontSize',FONT_SZ);
        title(axD,'(b) Radial growth velocity','Interpreter','latex','FontSize',TITLE_FS);
        legend(axD,'Interpreter','latex','NumColumns',min(NCrys,3), ...
               'FontSize',8,'Location','northeast');
        % Formula bottom right, away from the legend.
        text(axD, 0.98, 0.04, '$\dot{R}_i = \mathrm{d}R_i/\mathrm{d}t$', ...
             'Units','normalized','HorizontalAlignment','right', ...
             'VerticalAlignment','bottom','Interpreter','latex', ...
             'FontSize',8,'Color',[0.45 0.45 0.45]);
        set(axD,'TickLabelInterpreter','latex','FontSize',TICK_SZ);

        linkaxes([axC, axD], 'x');

        figPath = fullfile(outDir,['individual_' name '.png']);
        exportgraphics(fig2, figPath, 'Resolution', 300);
        fprintf('Individual figure saved: %s\n', figPath);
    end
end

%% ========================================================================
%  LOG-LOG: crystal area vs time since ITS OWN nucleation
%  ========================================================================
%  The slope of the fit gives the exponent n in A ~ (t - t_nuc)^n.
%    n around 1/2  -> diffusion-limited growth
%    n around 1    -> linear, interface-limited growth
%  Reference slopes are drawn dashed so the data can be read against them
%  without fitting anything.

if PLOT_LOGLOG
    G_ll = pickSet(dropletDirs_loglog, G, N_FINITE_DIFF, SMOOTH_WINDOW);

    if ~isempty(G_ll)
    nG_ll   = numel(G_ll);
    cols_ll = CMAP(1:min(nG_ll, size(CMAP,1)), :);
    if nG_ll > size(CMAP,1), cols_ll = [cols_ll; lines(nG_ll - size(CMAP,1))]; end

    figLL = figure('Color','w','Units','inches','Position',[1 6 8 5]);
    axLL  = axes(figLL); hold(axLL,'on');

    for k = 1:nG_ll
        clr = cols_ll(k,:);
        if isempty(G_ll(k).areaCrystal), continue; end
        leg_done = false;   % one legend entry per droplet, not per crystal
        for ci = 1:G_ll(k).NCrys
            if size(G_ll(k).areaCrystal,2) < ci, continue; end
            t_since = G_ll(k).t - G_ll(k).nuc_t(ci);
            A_ci    = G_ll(k).areaCrystal(:, ci)';
            m = t_since > 0 & ~isnan(A_ci) & A_ci > 0;   % log scale: strictly > 0
            if sum(m) < 4, continue; end
            dn = 'off';
            if ~leg_done, dn = 'on'; leg_done = true; end
            plot(axLL, t_since(m), A_ci(m), '-', ...
                 'Color', [clr 0.65], 'LineWidth', 0.9, ...
                 'DisplayName', strrep(G_ll(k).name,'_','\_'), ...
                 'HandleVisibility', dn);
        end
    end

    set(axLL,'XScale','log','YScale','log');
    drawnow; xl = xlim(axLL); yl = ylim(axLL);

    % Common left edge: every curve visible from 1 s after nucleation. Data
    % before that are loaded but fall outside the window, where a crystal of
    % a few pixels is dominated by segmentation noise anyway.
    xl(1) = max(xl(1), 1);

    t_ref = exp(linspace(log(max(xl(1),0.2)), log(xl(2)), 80));
    t_anc = exp(log(xl(1)) + 0.35*(log(xl(2))-log(xl(1))));
    A_anc = exp(log(yl(1)) + 0.35*(log(yl(2))-log(yl(1))));
    ref_n   = [0.5,   1,    2];
    ref_lbl = {'$n=1/2$','$n=1$','$n=2$'};
    ref_clr = [0.65 0.65 0.65; 0.35 0.35 0.35; 0.10 0.10 0.10];
    for si = 1:3
        n  = ref_n(si);
        C0 = A_anc / t_anc^n;
        yr = C0 * t_ref.^n;
        in = yr >= yl(1)*0.3 & yr <= yl(2)*3;
        if any(in)
            plot(axLL, t_ref(in), yr(in), '--', ...
                 'Color', ref_clr(si,:), 'LineWidth', 1.0, ...
                 'DisplayName', ref_lbl{si});
        end
    end

    set(axLL,'XScale','log','YScale','log');
    xlim(axLL, xl); ylim(axLL, yl);
    xlabel(axLL,'$t - t_\mathrm{nuc}$ (s)', 'Interpreter','latex','FontSize',FONT_SZ);
    ylabel(axLL,'Crystal area ($\mu$m$^2$)', 'Interpreter','latex','FontSize',TITLE_FS);
    title(axLL,'Crystal area growth --- log-log scaling', ...
          'Interpreter','latex','FontSize',FONT_SZ+6,'FontWeight','normal');
    text(axLL, 0.97, 0.97, ...
         'Dashed lines: power-law slopes $A \propto (t - t_\mathrm{nuc})^n$', ...
         'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
         'Interpreter','latex','FontSize',8,'Color',[0.45 0.45 0.45]);
    legend(axLL,'Interpreter','latex','Location','southeast','FontSize',9,'NumColumns',2);
    set(axLL,'TickLabelInterpreter','latex','FontSize',TICK_SZ); box(axLL,'on');

    exportgraphics(figLL, fullfile(outDir,'loglog_crystal_area.png'),'Resolution',300);
    fprintf('Saved: loglog_crystal_area.png\n');
    end
end

%% ========================================================================
%  NORMALISED CURVES: aligned at nucleation, area scaled to its own maximum
%  ========================================================================
%  x: t - t_nuc,i        y: A_i / A_i,max
%  Scaling each crystal by its own maximum removes the size difference and
%  leaves only the SHAPE of the growth, which is what should be compared
%  between crystals that ended up very different in size.

if PLOT_NORM_CURVES
    G_nc = pickSet(dropletDirs_norm, G, N_FINITE_DIFF, SMOOTH_WINDOW);

    if ~isempty(G_nc)
    nG_nc   = numel(G_nc);
    cols_nc = CMAP(1:min(nG_nc, size(CMAP,1)), :);
    if nG_nc > size(CMAP,1), cols_nc = [cols_nc; lines(nG_nc - size(CMAP,1))]; end

    figNC = figure('Color','w','Units','inches','Position',[7.5 6 8 5]);
    axNC  = axes(figNC); hold(axNC,'on');

    for k = 1:nG_nc
        clr = cols_nc(k,:);
        if isempty(G_nc(k).areaCrystal), continue; end
        leg_done = false;
        for ci = 1:G_nc(k).NCrys
            if size(G_nc(k).areaCrystal,2) < ci, continue; end
            t_since = G_nc(k).t - G_nc(k).nuc_t(ci);
            A_ci    = G_nc(k).areaCrystal(:, ci)';
            A_max   = max(A_ci(~isnan(A_ci)));
            if isnan(A_max) || A_max == 0, continue; end
            A_norm  = A_ci / A_max;
            m = t_since >= 0 & ~isnan(A_norm);
            if sum(m) < 3, continue; end
            dn = 'off';
            if ~leg_done, dn = 'on'; leg_done = true; end
            plot(axNC, t_since(m), A_norm(m), '-', ...
                 'Color', [clr 0.65], 'LineWidth', 0.9, ...
                 'DisplayName', strrep(G_nc(k).name,'_','\_'), ...
                 'HandleVisibility', dn);
        end
    end

    ylim(axNC, [0, 1.05]);
    xlabel(axNC,'$t - t_\mathrm{nuc,i}$ (s)', 'Interpreter','latex','FontSize',FONT_SZ);
    ylabel(axNC,'Normalized area $A_i / A_{i,\mathrm{max}}$', ...
           'Interpreter','latex','FontSize',TITLE_FS);
    title(axNC,'Normalized crystal growth curves (aligned at nucleation)', ...
          'Interpreter','latex','FontSize',FONT_SZ+1);
    text(axNC, 0.97, 0.05, ...
         ['Each curve: one crystal. Time shifted to $t_{\mathrm{nuc},i}$; ' ...
          'area normalized to $[0,1]$.'], ...
         'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom', ...
         'Interpreter','latex','FontSize',8,'Color',[0.45 0.45 0.45]);
    legend(axNC,'Interpreter','latex','Location','southeast','FontSize',9,'NumColumns',2);
    set(axNC,'TickLabelInterpreter','latex','FontSize',TICK_SZ); box(axNC,'on');

    exportgraphics(figNC, fullfile(outDir,'normalized_growth_curves.png'),'Resolution',300);
    fprintf('Saved: normalized_growth_curves.png\n');
    end
end

%% ========================================================================
%  RVEQ vs CRYSTALLISED FRACTION
%  ========================================================================
%  x: A_tot / A_tot,max (0 to 1)    y: dR_eq/dt (um/s)
%  Replacing time by the crystallised fraction removes the difference in
%  drying speed between droplets, so the shape of the curve reflects how the
%  front slows as the available salt is used up rather than how long the
%  droplet happened to take.

if PLOT_RVEQ_FRAC
    G_rv = pickSet(dropletDirs_rveq, G, N_FINITE_DIFF, SMOOTH_WINDOW);

    if ~isempty(G_rv)
    nG_rv   = numel(G_rv);
    cols_rv = CMAP(1:min(nG_rv, size(CMAP,1)), :);
    if nG_rv > size(CMAP,1), cols_rv = [cols_rv; lines(nG_rv - size(CMAP,1))]; end

    figRV = figure('Color','w','Units','inches','Position',[14 6 8 5]);
    axRV  = axes(figRV); hold(axRV,'on');

    for k = 1:nG_rv
        clr = cols_rv(k,:);
        m   = G_rv(k).t >= 0;
        A_tot = G_rv(k).A(m);
        dR_v  = G_rv(k).dR_eq(m);
        t_m   = G_rv(k).t(m);

        A_max_k = max(A_tot(~isnan(A_tot)));
        if isnan(A_max_k) || A_max_k == 0, continue; end
        frac = A_tot / A_max_k;

        plot(axRV, frac, dR_v, '-', ...
             'Color', clr, 'LineWidth', LINE_W, ...
             'DisplayName', strrep(G_rv(k).name,'_','\_'));

        % Nucleation markers, same style as figures 1a-1c. The x position is
        % found by locating the nucleation time within the plotted stretch.
        for ci = 1:G_rv(k).NCrys
            if G_rv(k).nuc_t(ci) < 0, continue; end
            [~, fi_rv] = min(abs(t_m - G_rv(k).nuc_t(ci)));
            plot(axRV, frac(fi_rv), dR_v(fi_rv), 'o', ...
                 'MarkerSize', 8, 'MarkerFaceColor','w', ...
                 'MarkerEdgeColor', clr, 'LineWidth', 1.2, ...
                 'HandleVisibility','off');
            text(axRV, frac(fi_rv), dR_v(fi_rv), num2str(ci), ...
                 'Color', clr, 'FontSize', 5, 'FontWeight','bold', ...
                 'HorizontalAlignment','center','VerticalAlignment','middle', ...
                 'Interpreter','none');
        end
    end

    plot(axRV, NaN, NaN, 'o', 'MarkerSize', 8, ...
         'MarkerFaceColor','w','MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',1.2, ...
         'DisplayName','New crystal nucleation');
    yline(axRV, 0, '--', 'Color',[0.5 0.5 0.5],'LineWidth',0.9,'HandleVisibility','off');
    xlim(axRV,[0, 1]);
    xlabel(axRV,'Crystallized fraction $A_\mathrm{tot}/A_\mathrm{tot,max}$', ...
           'Interpreter','latex','FontSize',FONT_SZ);
    ylabel(axRV,'$\dot{R}_\mathrm{eq}$ ($\mu$m s$^{-1}$)', ...
           'Interpreter','latex','FontSize',FONT_SZ);
    title(axRV,'Equivalent radial velocity vs.\ crystallized fraction', ...
          'Interpreter','latex','FontSize',FONT_SZ+1);
    legend(axRV,'Interpreter','latex','Location','northeast','FontSize',9,'NumColumns',2);
    set(axRV,'TickLabelInterpreter','latex','FontSize',TICK_SZ); box(axRV,'on');

    exportgraphics(figRV, fullfile(outDir,'rveq_vs_fraction.png'),'Resolution',300);
    fprintf('Saved: rveq_vs_fraction.png\n');
    end
end

%% ========================================================================
%  SUMMARY TABLE: growth statistics per droplet
%  ========================================================================
%  Three outputs: console, CSV and a rendered PNG. The CSV is the one that
%  matters downstream, since interDropletFigures reads it.
%
%  The means are taken over the POSITIVE values only: a negative dA/dt is a
%  segmentation artefact, not shrinkage, and averaging it in would bias the
%  mean growth rate downwards.

if MAKE_SUMMARY_TABLE
    G_tb = pickSet(dropletDirs_table, G, N_FINITE_DIFF, SMOOTH_WINDOW);

    if ~isempty(G_tb)
    nT_tb = numel(G_tb);

    tbl_name  = cell(nT_tb,1);
    tbl_ncrys = zeros(nT_tb,1);
    tbl_tnuc1 = zeros(nT_tb,1);
    tbl_tend  = zeros(nT_tb,1);
    tbl_dA_mn = zeros(nT_tb,1);
    tbl_dA_mx = zeros(nT_tb,1);
    tbl_dR_mn = zeros(nT_tb,1);
    tbl_dR_mx = zeros(nT_tb,1);

    for ki = 1:nT_tb
        tbl_name{ki}  = G_tb(ki).name;
        tbl_ncrys(ki) = G_tb(ki).NCrys;
        tbl_tnuc1(ki) = G_tb(ki).t_nuc1;   % absolute time of first nucleation
        tbl_tend(ki)  = G_tb(ki).t_end;    % absolute time of the last frame
        m     = G_tb(ki).t >= 0;
        dA_v  = G_tb(ki).dA(m);
        dR_v  = G_tb(ki).dR_eq(m);
        tbl_dA_mn(ki) = mean(dA_v(~isnan(dA_v) & dA_v > 0));
        tbl_dA_mx(ki) = max(dA_v(~isnan(dA_v)));
        tbl_dR_mn(ki) = mean(dR_v(~isnan(dR_v) & dR_v > 0));
        tbl_dR_mx(ki) = max(dR_v(~isnan(dR_v)));
    end

    % ---- console ----
    W = 100;
    fprintf('\n%s\n', repmat('=',1,W));
    fprintf('  SUMMARY TABLE -- Crystal Growth\n');
    fprintf('%s\n', repmat('=',1,W));
    fprintf('  %-14s  %5s  %9s  %9s  %11s  %11s  %10s  %10s\n', ...
            'Session','NCrys','t_nuc1(s)','t_end(s)','dA mean','dA max','dReq mean','dReq max');
    fprintf('  %-14s  %5s  %9s  %9s  %11s  %11s  %10s  %10s\n', ...
            '','','','','(um2/s)','(um2/s)','(um/s)','(um/s)');
    fprintf('  %s\n', repmat('-',1,W-2));
    for ki = 1:nT_tb
        fprintf('  %-14s  %5d  %9.1f  %9.1f  %11.1f  %11.1f  %10.2f  %10.2f\n', ...
                tbl_name{ki}, tbl_ncrys(ki), tbl_tnuc1(ki), tbl_tend(ki), ...
                tbl_dA_mn(ki), tbl_dA_mx(ki), tbl_dR_mn(ki), tbl_dR_mx(ki));
    end
    fprintf('%s\n\n', repmat('=',1,W));

    % ---- CSV: the file interDropletFigures reads ----
    T_csv = table(tbl_name, tbl_ncrys, tbl_tnuc1, tbl_tend, ...
                  tbl_dA_mn, tbl_dA_mx, tbl_dR_mn, tbl_dR_mx, ...
                  'VariableNames', {'Session','N_crystals', ...
                                    't_nuc1_s','t_end_s', ...
                                    'dAdt_mean_um2ps','dAdt_max_um2ps', ...
                                    'dReq_mean_umps','dReq_max_umps'});
    writetable(T_csv, fullfile(outDir,'summary_crystal_growth.csv'));
    fprintf('CSV saved: summary_crystal_growth.csv\n');

    % ---- rendered PNG table ----
    % Drawn by hand rather than with uitable, which does not export cleanly.
    figH_tb = max(2.5, 1.3 + nT_tb * 0.44);
    figTB   = figure('Color','w','Units','inches', ...
                     'Position',[1 12 13 figH_tb], ...
                     'MenuBar','none','ToolBar','none');
    axTB = axes(figTB,'Position',[0 0.04 1 0.86]);
    set(axTB,'Visible','off','XLim',[0 1],'YLim',[0 1]);

    cX   = [0.01, 0.19, 0.27, 0.36, 0.45, 0.57, 0.69, 0.82, 1.00];
    cHd  = {'Session', '$N_\mathrm{crys}$', ...
            '$t_{\mathrm{nuc},1}$ (s)', '$t_\mathrm{end}$ (s)', ...
            '$\dot{A}$ mean', '$\dot{A}$ max', ...
            '$\dot{R}_\mathrm{eq}$ mean', '$\dot{R}_\mathrm{eq}$ max'};
    cUn  = {'', '', '', '', ...
            '($\mu$m$^2$/s)', '($\mu$m$^2$/s)', ...
            '($\mu$m/s)', '($\mu$m/s)'};
    nC_tb  = 8;
    nR_tb  = nT_tb + 2;          % header + units + data rows
    rh_tb  = 0.78 / nR_tb;
    top_tb = 0.88;

    % Horizontal rules: heavy at the top and bottom, medium under the header.
    for ri = 0:nR_tb
        y_l  = top_tb - ri*rh_tb;
        lw_l = 0.4;
        if ri == 0 || ri == nR_tb, lw_l = 1.5;
        elseif ri == 2,            lw_l = 0.9; end
        line(axTB,[0.01 0.99],[y_l y_l],'Color',[0.25 0.25 0.25],'LineWidth',lw_l);
    end
    for ci2 = 2:nC_tb
        line(axTB,[cX(ci2) cX(ci2)],[top_tb-nR_tb*rh_tb, top_tb], ...
             'Color',[0.75 0.75 0.75],'LineWidth',0.4);
    end

    for ci2 = 1:nC_tb
        xc = (cX(ci2)+cX(ci2+1))/2;
        ha_h = 'center';
        if ci2 == 1, ha_h = 'left'; xc = cX(ci2)+0.006; end
        text(axTB, xc, top_tb-0.5*rh_tb, cHd{ci2}, ...
             'Interpreter','latex','FontSize',9,'FontWeight','bold', ...
             'HorizontalAlignment',ha_h,'VerticalAlignment','middle');
        if ~isempty(cUn{ci2})
            text(axTB, xc, top_tb-1.5*rh_tb, cUn{ci2}, ...
                 'Interpreter','latex','FontSize',8,'Color',[0.4 0.4 0.4], ...
                 'HorizontalAlignment',ha_h,'VerticalAlignment','middle');
        end
    end

    for ki = 1:nT_tb
        y_row = top_tb - (ki+1)*rh_tb;
        if mod(ki,2) == 0        % alternating band, easier to read across
            fill(axTB, [0.01 0.99 0.99 0.01], ...
                 [y_row y_row y_row+rh_tb y_row+rh_tb], ...
                 [0.93 0.93 0.97], 'EdgeColor','none');
        end
        vals = {tbl_name{ki}, tbl_ncrys(ki), tbl_tnuc1(ki), tbl_tend(ki), ...
                tbl_dA_mn(ki), tbl_dA_mx(ki), tbl_dR_mn(ki), tbl_dR_mx(ki)};
        for ci2 = 1:nC_tb
            xc = (cX(ci2)+cX(ci2+1))/2;
            ha_d = 'center';
            if ci2 == 1
                txt  = strrep(vals{ci2},'_','\_');
                ha_d = 'left'; xc = cX(ci2)+0.006;
            elseif ci2 == 2
                txt = sprintf('%d', vals{ci2});
            elseif ci2 == 3 || ci2 == 4   % times in seconds, one decimal
                txt = sprintf('%.1f', vals{ci2});
            else
                txt = sprintf('%.2f', vals{ci2});
            end
            text(axTB, xc, y_row+0.5*rh_tb, txt, ...
                 'Interpreter','latex','FontSize',9, ...
                 'HorizontalAlignment',ha_d,'VerticalAlignment','middle');
        end
    end

    annotation(figTB,'textbox',[0, 0.93, 1, 0.07], ...
               'String','Crystal Growth --- Summary Statistics', ...
               'Interpreter','latex','FontSize',FONT_SZ+1,'FontWeight','bold', ...
               'EdgeColor','none','HorizontalAlignment','center', ...
               'VerticalAlignment','middle');

    exportgraphics(figTB, fullfile(outDir,'summary_table.png'),'Resolution',300);
    fprintf('Saved: summary_table.png\n');
    end
end

fprintf('\nDone.\n');


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function Gout = pickSet(dirList, Gfallback, N_FD, SW)
%PICKSET  Load the droplets for one analysis, or reuse the global set.
%   An empty list means "use the same droplets as the global figure", which
%   is what keeps the configuration section short.
    if ~isempty(dirList)
        Gout = loadGStruct(dirList, N_FD, SW);
    elseif ~isempty(Gfallback)
        Gout = Gfallback;
    else
        Gout = [];
    end
end


function plotGlobalCurve(ax, G, cols, fieldName, nucFieldName, lineW)
%PLOTGLOBALCURVE  One curve per droplet with numbered nucleation markers.
%
%   Shared by figures 1a, 1b and 1c so the three are guaranteed to look the
%   same and differ only in the quantity plotted. Only t >= 0 is drawn.
%
%   The numbered circles are what let a step in the curve be attributed to a
%   specific new crystal instead of read as noise.

    nG = numel(G);
    for k = 1:nG
        clr = cols(k,:);
        m   = G(k).t >= 0;
        plot(ax, G(k).t(m), G(k).(fieldName)(m), '-', ...
             'Color', clr, 'LineWidth', lineW, ...
             'DisplayName', strrep(G(k).name,'_','\_'));
        for ci = 1:G(k).NCrys
            if G(k).nuc_t(ci) >= 0
                yv = G(k).(nucFieldName)(ci);
                plot(ax, G(k).nuc_t(ci), yv, 'o', ...
                     'MarkerSize', 8, 'MarkerFaceColor', 'w', ...
                     'MarkerEdgeColor', clr, 'LineWidth', 1.2, ...
                     'HandleVisibility','off');
                text(ax, G(k).nuc_t(ci), yv, num2str(ci), ...
                     'Color', clr, 'FontSize', 5, 'FontWeight', 'bold', ...
                     'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', 'Interpreter', 'none');
            end
        end
    end
    % A dummy point so the marker gets one legend entry rather than none.
    plot(ax, NaN, NaN, 'o', ...
         'MarkerSize', 8, 'MarkerFaceColor', 'w', ...
         'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
         'DisplayName', 'New crystal nucleation');
end


function [frameT, nucList] = timeAxisFromSpecs(dDir, NF, frameT, nucList, name)
%TIMEAXISFROMSPECS  Rebuild the time axis and nucleation list from the specs.
%
%   The frameTimes and nucList stored inside a .mat may be stale if t_start or
%   the nucleation frames were corrected after that file was written. The
%   specifications are authoritative; the arguments passed in are only used
%   as a fallback when they cannot be read.
%
%   This function READS the specifications and never writes to them.

    try
        [~, spR] = evalc('readSpecs(dDir)');   % quiet console
        tsR = 0; if ~isnan(spR.t_start), tsR = spR.t_start; end
        frameT = tsR + (0:NF-1) / spR.fps;              % frame 1 = t_start
        if isfield(spR,'nucleation_frames') && ~isempty(spR.nucleation_frames) ...
                                            && any(~isnan(spR.nucleation_frames))
            nfR = sort(spR.nucleation_frames(~isnan(spR.nucleation_frames)));
            nucList = min(max(round(nfR(:)), 1), NF);
        elseif ~isnan(spR.nucleation_frame)
            nucList = min(max(round(spR.nucleation_frame), 1), NF);
        end
    catch MEr
        warning('%s: cannot read the specifications (%s). Using the .mat.', ...
                name, MEr.message);
    end
end


function G = loadGStruct(dirs, N_FD, SW)
%LOADGSTRUCT  Load the crystal growth data of every droplet in `dirs`.
%
%   The single place where a droplet is read, so every analysis in this file
%   works from identically derived quantities.
%
%   N_FD and SW must match N_FINITE_DIFF and SMOOTH_WINDOW in
%   crystalGrowthAnalysis, otherwise the curves here would be derived
%   differently from the ones the pipeline produced.

    G = struct('t',{},'A',{},'dA',{},'R_eq',{},'dR_eq',{},...
               'nuc_t',{},'nuc_A',{},'nuc_dA',{},'nuc_Req',{},'nuc_dReq',{},...
               'areaCrystal',{},'nucList',{},'name',{},'NCrys',{}, ...
               't_nuc1',{},'t_end',{});

    for i = 1:numel(dirs)
        dDir = dirs{i};
        [~, folderName] = fileparts(dDir);
        if ~isfolder(dDir)
            warning('loadGStruct: folder not found: %s', dDir); continue
        end
        resF = fullfile(dDir, 'RESULTS_crystals');
        name = folderName;
        hit  = dir(fullfile(resF, 'crystal_areas_*.mat'));
        if isempty(hit)
            warning('loadGStruct: no crystal_areas_*.mat in %s', resF); continue
        end
        tk = extractBetween(hit(1).name, 'crystal_areas_', '.mat');
        if ~isempty(tk), name = tk{1}; end

        caFile = fullfile(resF, ['crystal_areas_' name '.mat']);
        if ~isfile(caFile)
            warning('loadGStruct: not found: %s', caFile); continue
        end
        ca = load(caFile, 'areasUm', 'frameTimes', 'nucList');

        areasUm = ca.areasUm;
        frameT  = ca.frameTimes(:)';
        nucList = ca.nucList(:);
        NF      = numel(frameT);
        nucList = min(max(round(nucList), 1), NF);

        [frameT, nucList] = timeAxisFromSpecs(dDir, NF, frameT, nucList, name);

        totalArea = cellfun(@sum, areasUm);
        totalArea = totalArea(:)';

        % Growth rate and equivalent-radius velocity, both by the same
        % centred finite difference used in the pipeline.
        growthRate = centredDiff(totalArea, frameT, N_FD, SW);
        R_eq       = sqrt(totalArea / pi);
        dR_eq      = centredDiff(R_eq, frameT, N_FD, SW);

        t0    = frameT(nucList(1));
        t_rel = frameT - t0;

        G(end+1).t      = t_rel;                       %#ok<AGROW>
        G(end).A        = totalArea;
        G(end).dA       = growthRate;
        G(end).R_eq     = R_eq;
        G(end).dR_eq    = dR_eq;
        G(end).nuc_t    = frameT(nucList) - t0;
        G(end).nuc_A    = totalArea(nucList);
        G(end).nuc_dA   = growthRate(nucList);
        G(end).nuc_Req  = R_eq(nucList);
        G(end).nuc_dReq = dR_eq(nucList);
        G(end).name     = folderName;
        G(end).NCrys    = numel(nucList);

        % ---- t_nuc1 and t_end, derived from the specifications ----
        % Computed rather than stored, so they follow automatically if fps or
        % t_start are ever corrected. Nothing is written back to the
        % specifications file.
        try
            [~, sp_rc] = evalc('readSpecs(dDir)');
            fp_rc = sp_rc.fps;
            ts_rc = 0;
            if ~isnan(sp_rc.t_start), ts_rc = sp_rc.t_start; end

            if ~isnan(sp_rc.nucleation_time)
                G(end).t_nuc1 = sp_rc.nucleation_time;
            else
                % Derive it from the frame, which is the exact quantity.
                fn_rc = NaN;
                if ~isnan(sp_rc.nucleation_frame)
                    fn_rc = sp_rc.nucleation_frame;
                elseif isfield(sp_rc,'nucleation_frames') && ~isempty(sp_rc.nucleation_frames)
                    fn_arr = sp_rc.nucleation_frames(isfinite(sp_rc.nucleation_frames));
                    if ~isempty(fn_arr), fn_rc = min(fn_arr); end
                end
                if ~isnan(fn_rc)
                    G(end).t_nuc1 = ts_rc + (fn_rc - 1) / fp_rc;
                else
                    G(end).t_nuc1 = t0;
                end
            end

            if ~isnan(sp_rc.frame_max)
                G(end).t_end = ts_rc + (sp_rc.frame_max - 1) / fp_rc;
            else
                G(end).t_end = frameT(end);
            end
        catch
            G(end).t_nuc1 = t0;
            G(end).t_end  = frameT(end);
        end

        % ---- per-crystal areas, for the log-log and normalised figures ----
        giFile = fullfile(resF, ['growth_individual_' name '.mat']);
        if isfile(giFile)
            gi_tmp = load(giFile, 'areaCrystal_um2');
            ac_tmp = double(gi_tmp.areaCrystal_um2);
            if size(ac_tmp,1) ~= NF, ac_tmp = ac_tmp'; end   % force [NF x NCrys]
            NCrys_tmp = min(G(end).NCrys, size(ac_tmp,2));
            for ci_tmp = 1:NCrys_tmp
                f0_tmp = nucList(ci_tmp);
                if f0_tmp > 1 && f0_tmp <= NF
                    ac_tmp(1:f0_tmp-1, ci_tmp) = NaN;   % before it existed
                end
            end
            G(end).areaCrystal = ac_tmp;
            G(end).nucList     = nucList;
        else
            warning(['No growth_individual_%s.mat -> log-log and normalised ' ...
                     'figures unavailable for this droplet'], name);
            G(end).areaCrystal = [];
            G(end).nucList     = [];
        end
    end
end


function d = centredDiff(y, t, N_FD, SW)
%CENTREDDIFF  Centred finite difference over a gap of N_FD, then smoothed.
%
%   dy/dt|i = ( y(i+N) - y(i-N) ) / ( t(i+N) - t(i-N) )
%
%   One-sided at the two ends. The gap is what makes the derivative usable:
%   over a single frame the change is comparable to the segmentation noise.
%   Identical to the method in crystalGrowthAnalysis, which is why N_FD and
%   SW must match its settings.

    NF = numel(y);
    d  = nan(1, NF);
    for fi = 1:NF
        if fi <= N_FD
            if fi + N_FD <= NF
                d(fi) = (y(fi+N_FD) - y(fi)) / (t(fi+N_FD) - t(fi));
            end
        elseif fi > NF - N_FD
            if fi - N_FD >= 1
                d(fi) = (y(fi) - y(fi-N_FD)) / (t(fi) - t(fi-N_FD));
            end
        else
            d(fi) = (y(fi+N_FD) - y(fi-N_FD)) / (t(fi+N_FD) - t(fi-N_FD));
        end
    end
    d = smoothdata(d, 'movmean', SW);
end
