%% ========================================================================
%  DIMENSIONALANALYSISFIGURES  Predicted vs measured velocities.
%  ========================================================================
%
%  PRODUCES: FIGURES ONLY. Nothing here is recomputed and no data file is
%  written; this script draws the outcome of the dimensional analysis.
%
%  WHAT THESE FIGURES SHOW
%    Each candidate mechanism predicts a characteristic velocity. Plotting
%    those predictions against the velocities actually measured is what
%    discriminates between them: a mechanism that is off by five orders of
%    magnitude is not the one driving the flow, whatever else it might
%    explain.
%
%    FIGURE A : everything against the POST-nucleation flow.
%               Five groups: capillary vs mean, capillary vs max,
%               natural convection vs mean, natural convection vs max,
%               Marangoni vs mean. Marangoni vs max is omitted because it
%               is even further off and would only stretch the axis.
%
%    FIGURE B : capillary against the PRE-nucleation flow.
%               (a) ratios vs mean and vs max, (b) predicted vs measured.
%
%    The pre- and post-nucleation comparisons are kept apart on purpose.
%    Capillary replenishment is driven by evaporation and operates from the
%    start, so it belongs against the pre-nucleation flow. The other two are
%    driven by the gradient the crystal creates and only exist afterwards.
%
%  WHY THE MEDIAN AND NOT THE MEAN
%    The ratio distributions are strongly skewed across droplets, so the
%    median line drawn on each group is the representative statistic. Points
%    are jittered horizontally so overlapping droplets stay countable.
%
%  NOTE ON UNITS AND HEIGHTS
%    v_Ca is evaluated at the initial height h0; v_Ma and v_Ra at the height
%    at nucleation h_nuc. See dimensionalAnalysis.m for why.
%
%  IMPORTANT: THE NUMBERS BELOW ARE HARD-CODED
%    vCa, vMa and vRa are transcribed from the output of dimensionalAnalysis
%    rather than read from dimensional_per_droplet.csv. If that analysis is
%    ever rerun with different assumptions, this figure will NOT follow.
%    Re-transcribe the block below, or point it at the CSV.
%
%  OUTPUT  (in imagenes/dimenanalysis)
%    dimensional_postnucleation.pdf
%    dimensional_prenucleation.pdf
%
%  See also DIMENSIONALANALYSIS.

figureStyle;

outDir = fullfile('imagenes','dimenanalysis');
if ~exist(outDir,'dir'); mkdir(outDir); end

%% ---------------- data ----------------
id  = [4 5 6 7 8 9 10 12 13 14 15 18 19 20 21 22 23 24 25 29 30 31 32 33 34 35];

% Predicted characteristic velocities (um/s), from dimensionalAnalysis
vCa = [0.87 0.79 0.98 1.93 1.72 1.96 2.06 1.02 1.31 0.99 2.11 1.30 1.74 ...
       1.97 1.41 1.49 0.96 0.84 1.14 1.41 1.06 1.09 1.14 1.10 1.01 1.25];

vMa = [43.5 79.0 64.1 86.2 75.1 51.6 64.0 168.7 62.4 153.7 25.1 119.6 ...
       85.6 90.4 75.4 50.9 114.7 171.7 182.8 62.2 88.3 136.2 67.2 140.2 ...
       147.2 89.6]*1e3;

vRa = [1.16 7.08 3.43 13.79 8.33 2.39 4.78 31.70 3.58 45.87 0.21 25.93 ...
       5.66 6.60 7.00 3.78 12.73 49.79 74.84 4.60 17.15 36.65 5.39 50.40 ...
       80.23 11.62];

% Measured particle-tracking velocities (um/s). NaN where the droplet
% already had crystals when recording started, so it has no pre-nucleation
% interval.
tmPre    = [NaN 1.34 1.46 7.36 8.57 8.10 9.87 NaN 3.31 NaN 3.09 NaN 2.58 NaN ...
            2.43 2.83 6.47 9.26 4.27 2.76 1.54 4.16 0.68 2.84 NaN NaN];

tmMaxPre = [NaN 8.1 5.7 40.3 39.9 38.5 40.0 NaN 11.3 NaN 11.1 NaN 10.5 NaN ...
            10.6 10.1 22.7 26.7 15.4 10.0 10.5 19.5 4.1 26.7 NaN NaN];

tmPost = [2.58 4.19 3.44 11.24 12.07 11.77 13.67 3.51 6.92 3.34 1.83 3.31 ...
          6.67 7.56 2.24 2.16 8.62 9.25 5.68 1.81 1.70 3.87 0.83 2.85 2.83 2.02];

tmMax  = [15.0 19.8 23.5 53.7 52.2 55.0 47.9 16.0 37.0 19.5 9.3 17.9 30.0 ...
          28.6 12.7 16.0 28.6 31.1 32.1 11.8 12.5 32.0 6.4 22.4 29.4 19.5];

n = numel(id);

%% ---------------- colours, one per mechanism ----------------
cCa = [0.400 0.635 0.400];   % capillary
cRa = [0.204 0.427 0.616];   % natural convection
cMa = [0.839 0.376 0.302];   % Marangoni


%% ========================================================================
%  FIGURE A : everything against the POST-nucleation flow
%  ========================================================================

rCaPost = vCa ./ tmPost;
rCaMax  = vCa ./ tmMax;
rRa     = vRa ./ tmPost;
rRaMax  = vRa ./ tmMax;
rMa     = vMa ./ tmPost;

figA = figure('Units','centimeters','Position',[2 2 13 17],'Color','w');
tiledlayout(figA,2,1,'TileSpacing','compact','Padding','compact');

% ----------------- panel (a) : ratios -----------------
ax1 = nexttile; hold(ax1,'on');

xlimA = [0.4 5.9];
% Shaded band = within a factor of 10 either way, the usual tolerance for
% an order-of-magnitude scaling argument.
patch(ax1,[xlimA(1) xlimA(2) xlimA(2) xlimA(1)],[0.1 0.1 10 10], ...
      [0.92 0.95 0.92],'EdgeColor','none');
plot(ax1,xlimA,[1 1],'-','Color',[0.2 0.2 0.2],'LineWidth',1.2);

grupos = {rCaPost, rCaMax, rRa, rRaMax, rMa};
cols   = [cCa; cCa; cRa; cRa; cMa];
hollow = [false true false true false];   % open marker for the "vs max" pairs

rng(0);   % fixed seed: the jitter must be reproducible between runs
for j = 1:5
    r  = grupos{j};
    ok = ~isnan(r);
    xj = j + (rand(1,n)-0.5)*0.30;
    if hollow(j)
        plot(ax1,xj(ok), r(ok),'o','MarkerSize',6,'MarkerFaceColor','w', ...
             'MarkerEdgeColor',cols(j,:),'LineWidth',1.0);
    else
        plot(ax1,xj(ok), r(ok),'o','MarkerSize',6,'MarkerFaceColor',cols(j,:), ...
             'MarkerEdgeColor','w','LineWidth',0.5);
    end
    m = median(r,'omitnan');
    plot(ax1,[j-0.32 j+0.32],[m m],'-','Color',cols(j,:)*0.55,'LineWidth',2.5);

    % The Marangoni median is of order 1e4, so it needs its own format.
    if m >= 100
        e   = floor(log10(m));
        txt = sprintf('%.0f\\times10^{%d}\\times', m/10^e, e);
    elseif m < 0.1
        txt = sprintf('%.3f\\times', m);
    else
        txt = sprintf('%.2f\\times', m);
    end
    text(ax1,j+0.36, m, txt,'FontSize',8.5,'FontWeight','bold', ...
         'Color',cols(j,:)*0.55,'VerticalAlignment','middle','Interpreter','tex');
end

etiq = {sprintf('Capillary\\newlinevs mean'), ...
        sprintf('Capillary\\newlinevs max'), ...
        sprintf('Nat. conv.\\newlinevs mean'), ...
        sprintf('Nat. conv.\\newlinevs max'), ...
        sprintf('Marangoni\\newlinevs mean')};

set(ax1,'YScale','log','YLim',[1e-2 2e6],'XLim',xlimA, ...
        'XTick',1:5,'XTickLabel',etiq,'FontSize',8.5, ...
        'TickLabelInterpreter','tex', ...
        'YTick',10.^(-2:2:6),'TickDir','out','Box','on', ...
        'YGrid','on','GridAlpha',0.10,'YMinorGrid','off','Layer','top');
ylabel(ax1,'Predicted / measured','FontSize',9);

text(ax1,0.46, 1.8,'perfect agreement','FontSize',7.5,'Color',[0.2 0.2 0.2]);
text(ax1,0.46, 0.13,'within a factor of 10','FontSize',7.5,'Color',[0.35 0.5 0.35]);

% ----------------- panel (b) : predicted vs measured -----------------
ax2 = nexttile; hold(ax2,'on');

lim = [3e-1 5e5];
patch(ax2,[lim(1) lim(2) lim(2) lim(1)],[lim(1)*10 lim(2)*10 lim(2)/10 lim(1)/10], ...
      [0.92 0.95 0.92],'EdgeColor','none');
plot(ax2,lim,lim,'-','Color',[0.2 0.2 0.2],'LineWidth',1.2);

plot(ax2,tmPost, vCa,'s','MarkerSize',6,'MarkerFaceColor',cCa, ...
     'MarkerEdgeColor','w','LineWidth',0.5);
plot(ax2,tmMax,  vCa,'s','MarkerSize',6,'MarkerFaceColor','w', ...
     'MarkerEdgeColor',cCa,'LineWidth',1.0);
plot(ax2,tmPost, vRa,'o','MarkerSize',6,'MarkerFaceColor',cRa, ...
     'MarkerEdgeColor','w','LineWidth',0.5);
plot(ax2,tmMax,  vRa,'o','MarkerSize',6,'MarkerFaceColor','w', ...
     'MarkerEdgeColor',cRa,'LineWidth',1.0);
plot(ax2,tmPost, vMa,'v','MarkerSize',6,'MarkerFaceColor',cMa, ...
     'MarkerEdgeColor','w','LineWidth',0.5);

set(ax2,'XScale','log','YScale','log','XLim',lim,'YLim',lim, ...
        'XTick',10.^(0:1:5),'YTick',10.^(0:1:5),'FontSize',9, ...
        'TickDir','out','Box','on','XGrid','on','YGrid','on', ...
        'GridAlpha',0.10,'XMinorGrid','off','YMinorGrid','off','Layer','top');
pbaspect(ax2,[1 1 1]);
xlabel(ax2,'Measured post-nucleation velocity  (\mum s^{-1})','FontSize',9);
ylabel(ax2,'Predicted velocity  (\mum s^{-1})','FontSize',9);

% Labelled in place rather than with a legend: the three clouds sit far
% enough apart that a legend box would only cover data.
text(ax2,1.5e3, 2e5,'Marangoni','Color',cMa,'FontSize',9,'FontWeight','bold');
text(ax2,60, 4,'Nat. convection','Color',cRa,'FontSize',9,'FontWeight','bold');
text(ax2,1.0, 0.55,'Capillary','Color',cCa,'FontSize',9,'FontWeight','bold');
text(ax2,4e3, 1.5e3,'1:1','FontSize',8,'Color',[0.2 0.2 0.2],'Rotation',45);

annotation(figA,'textbox',[0.02 0.94 0.06 0.04],'String','(a)', ...
           'EdgeColor','none','FontSize',10,'FontWeight','bold');
annotation(figA,'textbox',[0.02 0.44 0.06 0.04],'String','(b)', ...
           'EdgeColor','none','FontSize',10,'FontWeight','bold');

exportgraphics(figA, fullfile(outDir,'dimensional_postnucleation.pdf'), ...
               'ContentType','vector','BackgroundColor','white');
fprintf('Figure A saved.\n');


%% ========================================================================
%  FIGURE B : capillary against the PRE-nucleation flow
%  ========================================================================

rCaPre    = vCa ./ tmPre;
rCaPreMax = vCa ./ tmMaxPre;

figB = figure('Units','centimeters','Position',[2 2 13 7],'Color','w');
tiledlayout(figB,1,2,'TileSpacing','compact','Padding','compact');

% ----------------- panel (a) : ratios -----------------
bx1 = nexttile; hold(bx1,'on');

xlimB = [0.4 2.9];
patch(bx1,[xlimB(1) xlimB(2) xlimB(2) xlimB(1)],[0.1 0.1 10 10], ...
      [0.92 0.95 0.92],'EdgeColor','none');
plot(bx1,xlimB,[1 1],'-','Color',[0.2 0.2 0.2],'LineWidth',1.2);

gruposB = {rCaPre, rCaPreMax};
hollowB = [false true];

rng(0);
for j = 1:2
    r  = gruposB{j};
    ok = ~isnan(r);
    xj = j + (rand(1,n)-0.5)*0.30;
    if hollowB(j)
        plot(bx1,xj(ok), r(ok),'s','MarkerSize',6,'MarkerFaceColor','w', ...
             'MarkerEdgeColor',cCa,'LineWidth',1.0);
    else
        plot(bx1,xj(ok), r(ok),'s','MarkerSize',6,'MarkerFaceColor',cCa, ...
             'MarkerEdgeColor','w','LineWidth',0.5);
    end
    m = median(r,'omitnan');
    plot(bx1,[j-0.32 j+0.32],[m m],'-','Color',cCa*0.55,'LineWidth',2.5);
    text(bx1,j+0.36, m, sprintf('%.2f\\times',m),'FontSize',8.5, ...
         'FontWeight','bold','Color',cCa*0.55, ...
         'VerticalAlignment','middle','Interpreter','tex');
end

set(bx1,'YScale','log','YLim',[1e-2 1e2],'XLim',xlimB, ...
        'XTick',1:2,'FontSize',8.5,'TickLabelInterpreter','tex', ...
        'XTickLabel',{sprintf('vs mean'),sprintf('vs max')}, ...
        'YTick',10.^(-2:1:2),'TickDir','out','Box','on', ...
        'YGrid','on','GridAlpha',0.10,'YMinorGrid','off','Layer','top');
ylabel(bx1,'Predicted / measured','FontSize',9);
text(bx1,0.46, 1.5,'perfect agreement','FontSize',7,'Color',[0.2 0.2 0.2]);

% ----------------- panel (b) : predicted vs measured -----------------
bx2 = nexttile; hold(bx2,'on');

limB = [3e-1 1e2];
patch(bx2,[limB(1) limB(2) limB(2) limB(1)], ...
      [limB(1)*10 limB(2)*10 limB(2)/10 limB(1)/10], ...
      [0.92 0.95 0.92],'EdgeColor','none');
plot(bx2,limB,limB,'-','Color',[0.2 0.2 0.2],'LineWidth',1.2);

plot(bx2,tmPre,    vCa,'s','MarkerSize',6,'MarkerFaceColor',cCa, ...
     'MarkerEdgeColor','w','LineWidth',0.5);
plot(bx2,tmMaxPre, vCa,'s','MarkerSize',6,'MarkerFaceColor','w', ...
     'MarkerEdgeColor',cCa,'LineWidth',1.0);

set(bx2,'XScale','log','YScale','log','XLim',limB,'YLim',limB, ...
        'XTick',10.^(0:1:2),'YTick',10.^(0:1:2),'FontSize',9, ...
        'TickDir','out','Box','on','XGrid','on','YGrid','on', ...
        'GridAlpha',0.10,'XMinorGrid','off','YMinorGrid','off','Layer','top');
pbaspect(bx2,[1 1 1]);
xlabel(bx2,'Measured pre-nucleation velocity  (\mum s^{-1})','FontSize',9);
ylabel(bx2,'Predicted velocity  (\mum s^{-1})','FontSize',9);
text(bx2,25, 15,'1:1','FontSize',8,'Color',[0.2 0.2 0.2],'Rotation',45);

annotation(figB,'textbox',[0.01 0.92 0.06 0.06],'String','(a)', ...
           'EdgeColor','none','FontSize',10,'FontWeight','bold');
annotation(figB,'textbox',[0.50 0.92 0.06 0.06],'String','(b)', ...
           'EdgeColor','none','FontSize',10,'FontWeight','bold');

exportgraphics(figB, fullfile(outDir,'dimensional_prenucleation.pdf'), ...
               'ContentType','vector','BackgroundColor','white');
fprintf('Figure B saved.\n');


%% ---------------- console check ----------------
fprintf('\nMedians of the ratios:\n');
fprintf('  POST-NUCLEATION\n');
fprintf('    Capillary  vs TM post : %7.3f\n', median(rCaPost));
fprintf('    Capillary  vs TM max  : %7.3f\n', median(rCaMax));
fprintf('    Nat. conv. vs TM post : %7.3f\n', median(rRa));
fprintf('    Nat. conv. vs TM max  : %7.3f\n', median(rRaMax));
fprintf('    Marangoni  vs TM post : %7.2e\n', median(rMa));
fprintf('  PRE-NUCLEATION (n = %d)\n', sum(~isnan(rCaPre)));
fprintf('    Capillary  vs TM pre  : %7.3f\n', median(rCaPre,'omitnan'));
fprintf('    Capillary  vs TM max  : %7.3f\n', median(rCaPreMax,'omitnan'));
