%% ========================================================================
%  GROWTHFLOWCOMPARISON  Crystal growth against the measured flow.
%  ========================================================================
%
%  PRODUCES: FIGURES AND DATA. This is a pipeline stage, not a figure
%  script. The .mat files it writes are consumed by the figure scripts in
%  figures/, so it must run before them.
%
%  WHAT THIS DOES
%    Puts the crystal growth rate and the two independent flow measurements
%    on the same axis, in the same units, with the same smoothing. That is
%    the whole point of the comparison: if the flow is driven by the growing
%    crystal, the two should rise together.
%
%    One global plot (all three curves) and one panel per crystal (all three
%    curves). Nothing else.
%
%  WHAT IS COMPARED
%    GLOBAL      crystal   -> frontVel_perim, that is dA/dt normalised by
%                             the total perimeter, which is the mean normal
%                             velocity of the interface in um/s
%                PIV       -> spatial mean of |u| per frame
%                TrackMate -> mean spot speed per frame
%
%    INDIVIDUAL  crystal c -> vRadial(:,c), dR/dt of its mean radius
%                PIV       -> mean inside maskPlusBuffer{frame,c}
%                TrackMate -> mean of the spots inside the same buffer
%
%    The per-crystal version is the stronger test: it asks whether the flow
%    around THIS crystal tracks the growth of THIS crystal, rather than
%    whether two droplet-wide averages happen to move together.
%
%  TIME AXES (the specifications are authoritative)
%    PIV and TrackMate start at t_start. The crystal curves start at their
%    nucleation: the first one for the global plot, its own for each panel.
%    The global x axis is forced to begin at t_start so that the flow is
%    visibly present BEFORE any crystal exists, which is what separates the
%    pre-nucleation regime from the crystallisation-driven one.
%
%  SMOOTHING
%    SMOOTH_WINDOW from the pipeline for everything, so the three curves
%    carry exactly the same filter and their comparison is not an artefact
%    of different processing.
%
%  This is a SCRIPT, not a function: it is called with run() from
%  runDropletPipeline and shares its workspace.
%
%  INPUTS FROM THE PIPELINE
%    resultsFolder, dropletName, frameTimes, t_start, fps, um_per_px,
%    nucFrames, nucTimes, t_nucleation, NCrys, SMOOTH_WINDOW, crystalColors,
%    colorPIV, colorTM, colorCrystal, SAME_AXIS, LW_PLOT, IND_PER_FIG,
%    GLOBAL_CRYSTAL_METRIC and the DO_COMPARE_* flags.
%
%  OUTPUT FILES  (in resultsFolder)
%    compare_global_<droplet>.mat / .png
%    compare_individual_<droplet>.mat
%    compare_ind_panel<k>_<droplet>.png
%
%  See also RUNDROPLETPIPELINE, CRYSTALGROWTHANALYSIS, PIVFLOWANALYSIS,
%  TRACKFLOWANALYSIS.

% ---- the time axis ALWAYS comes from the specifications ----
if exist('frameTimes_master','var'), frameTimes = frameTimes_master; end
tAll = frameTimes(:);
nT   = numel(tAll);

if ~exist('nucTimes','var') || isempty(nucTimes)
    nucTimes = tAll(min(max(nucFrames,1), nT));
end
if ~exist('t_nucleation','var') || isnan(t_nucleation)
    t_nucleation = nucTimes(1);
end

% Stored series can be one or two frames longer than the axis, depending on
% which stage wrote them; this clips any series to the canonical length.
trim = @(v) v(1:min(numel(v), nT));

%% ========================================================================
%  PART 1 : GLOBAL, crystal front velocity vs PIV vs TrackMate
%  ========================================================================
if DO_COMPARE_GLOBAL_ALL || DO_COMPARE_GLOBAL_PIV || DO_COMPARE_GLOBAL_TM

    % ---------- the crystal curve, per GLOBAL_CRYSTAL_METRIC ----------
    if ~exist('GLOBAL_CRYSTAL_METRIC','var'), GLOBAL_CRYSTAL_METRIC = 'perim'; end

    switch lower(GLOBAL_CRYSTAL_METRIC)
        case 'perim'
            rf = fullfile(resultsFolder, ['growth_global_' dropletName '.mat']);
            assert(isfile(rf), '%s is missing (global crystal stage).', rf);
            Rc = load(rf, 'frontVel_perim');
            vGrowth   = Rc.frontVel_perim(:);
            growthLbl = 'Crystal front velocity ($\mu$m/s)';
            growthLeg = 'Crystal front ($\dot{A}/P$)';

        case 'req'
            rf = fullfile(resultsFolder, ['growth_global_' dropletName '.mat']);
            assert(isfile(rf), '%s is missing (global crystal stage).', rf);
            Rc = load(rf, 'frontVel_Req');
            vGrowth   = Rc.frontVel_Req(:);
            growthLbl = 'Equivalent radius growth ($\mu$m/s)';
            growthLeg = 'Crystal $\mathrm{d}R_{\mathrm{eq}}/\mathrm{d}t$';
            if numel(nucFrames) > 1
                warning(['GLOBAL_CRYSTAL_METRIC = ''Req'' with %d crystals: ' ...
                         'R_eq lumps every area into ONE circle, so the ' ...
                         'velocity is inflated by about sqrt(%d) = %.2f and ' ...
                         'jumps at each nucleation. Consider ''perim'' or ' ...
                         '''radial''.'], numel(nucFrames), numel(nucFrames), ...
                         sqrt(numel(nucFrames)));
            end

        case 'radial'
            rf = fullfile(resultsFolder, ['radial_velocity_' dropletName '.mat']);
            assert(isfile(rf), '%s is missing (per-crystal radial stage).', rf);
            Rc = load(rf, 'vRadialGlobal');
            vGrowth   = Rc.vRadialGlobal(:);
            growthLbl = 'Mean radial growth velocity ($\mu$m/s)';
            growthLeg = 'Crystal $\langle\dot{R}\rangle$';

        otherwise
            error(['GLOBAL_CRYSTAL_METRIC = ''%s'' is not valid. ' ...
                   'Use ''perim'', ''Req'' or ''radial''.'], GLOBAL_CRYSTAL_METRIC);
    end

    vGrowth = trim(vGrowth);
    tGrowth = tAll(1:numel(vGrowth));
    vGrowth(tGrowth < t_nucleation) = NaN;      % no crystal before nucleation
    fprintf('Global: crystal curve = %s\n', GLOBAL_CRYSTAL_METRIC);

    % ---------- PIV: spatial mean per frame ----------
    vFlowPIV = []; tPIV = [];
    pf = fullfile(resultsFolder, ['piv_fields_' dropletName '.mat']);
    if isfile(pf)
        P  = load(pf, 'velocity_magnitude_from_u_v_avg');
        Vc = P.velocity_magnitude_from_u_v_avg;
        nP = numel(Vc);
        vFlowPIV = nan(nP,1);
        for k = 1:nP
            V = Vc{k};
            if isempty(V), continue; end
            V = V(isfinite(V));
            if ~isempty(V), vFlowPIV(k) = mean(V); end
        end
        validP   = ~isnan(vFlowPIV);
        vFlowPIV = smoothdata(vFlowPIV, 'movmean', SMOOTH_WINDOW, 'omitnan');
        vFlowPIV(~validP) = NaN;                % do not invent the edges
        vFlowPIV = trim(vFlowPIV);
        tPIV     = tAll(1:numel(vFlowPIV));
    else
        warning('%s not found: no PIV curve.', pf);
    end

    % ---------- TrackMate: mean speed per frame ----------
    % The masked version is preferred: it excludes the spots detected on the
    % crystal surface, which would otherwise contaminate this comparison with
    % the crystal's own growth.
    vFlowTM = []; tTM = [];
    fM = fullfile(resultsFolder, ['trackmate_meanspeed_masked_' dropletName '.mat']);
    fP = fullfile(resultsFolder, ['trackmate_meanspeed_' dropletName '.mat']);
    fTM = '';
    if isfile(fM), fTM = fM; elseif isfile(fP), fTM = fP; end
    if ~isempty(fTM)
        M  = load(fTM, 'meanSpeed', 'framesU');
        vFlowTM = M.meanSpeed(:);
        fu      = double(M.framesU(:));
        tTM     = tAll(1) + fu / fps;                    % anchored at t_start
        ok      = (fu + 1) >= 1 & (fu + 1) <= nT;
        tTM(ok) = tAll(fu(ok) + 1);
    else
        warning('No trackmate_meanspeed file: no TrackMate curve.');
    end

    % ---------- the single global plot ----------
    if ~exist('SAME_AXIS','var'), SAME_AXIS = true; end
    if ~exist('LW_PLOT','var'),   LW_PLOT   = 1.6;  end

    figG = figure('Name','Global: crystal vs PIV vs TrackMate','Color','w', ...
                  'Position',[100 100 1000 560]);
    ax = axes(figG);
    hold(ax,'on')
    hs = gobjects(0);  leg = {};

    if SAME_AXIS
        % All three on the SAME axis. They are all velocities in um/s, so
        % this compares real magnitudes. Two independent y scales could be
        % stretched to make the curves agree or disagree at will, which
        % would make the figure an argument rather than a measurement.
        hs(end+1) = plot(ax, tGrowth, vGrowth, '-', 'LineWidth',LW_PLOT, 'Color',colorCrystal);
        leg{end+1} = growthLeg;
        if ~isempty(vFlowPIV)
            hs(end+1) = plot(ax, tPIV, vFlowPIV, '-', 'LineWidth',LW_PLOT, 'Color',colorPIV);
            leg{end+1} = 'PIV (Eulerian)';
        end
        if ~isempty(vFlowTM)
            hs(end+1) = plot(ax, tTM, vFlowTM, '-', 'LineWidth',LW_PLOT, 'Color',colorTM);
            leg{end+1} = 'TrackMate (Lagrangian)';
        end
        ylabel(ax, 'Velocity ($\mu$m/s)');
        set(ax, 'YColor', [0.15 0.15 0.15]);
    else
        yyaxis(ax,'left')
        hs(end+1) = plot(ax, tGrowth, vGrowth, '-', 'LineWidth',LW_PLOT, 'Color',colorCrystal);
        leg{end+1} = growthLeg;
        ylabel(ax, growthLbl);
        set(ax,'YColor',[0.15 0.15 0.15]);
        yyaxis(ax,'right'); hold(ax,'on')
        if ~isempty(vFlowPIV)
            hs(end+1) = plot(ax, tPIV, vFlowPIV, '-', 'LineWidth',LW_PLOT, 'Color',colorPIV);
            leg{end+1} = 'PIV (Eulerian)';
        end
        if ~isempty(vFlowTM)
            hs(end+1) = plot(ax, tTM, vFlowTM, '-', 'LineWidth',LW_PLOT, 'Color',colorTM);
            leg{end+1} = 'TrackMate (Lagrangian)';
        end
        ylabel(ax, 'Flow velocity ($\mu$m/s)');
        set(ax,'YColor',[0.15 0.15 0.15]);
    end

    xlabel(ax, 'Time (s)');
    title(ax, sprintf('%s: crystal growth vs flow', strrep(dropletName,'_','\_')));
    legend(ax, hs, leg, 'Location','northwest');
    nucleationLines(ax, nucTimes, [], true);    % neutral marks, not competing
    % The global plot DOES start at t_start: the point of this figure is to
    % show that the flow is already there before the first crystal appears.
    figureStyle(ax, tAll, t_start);

    exportgraphics(figG, fullfile(resultsFolder, ...
        ['compare_global_' dropletName '.png']), 'Resolution', 300);
    fprintf('Saved: compare_global_%s.png\n', dropletName);

    % ---------- save the ALREADY SMOOTHED series ----------
    save(fullfile(resultsFolder, ['compare_global_' dropletName '.mat']), ...
         'tGrowth','vGrowth','tPIV','vFlowPIV','tTM','vFlowTM', ...
         'nucTimes','t_nucleation','t_start','fps','SMOOTH_WINDOW', ...
         'GLOBAL_CRYSTAL_METRIC');
    fprintf('Global comparison saved (smoothed series).\n');
end

%% ========================================================================
%  PART 2 : INDIVIDUAL, radial growth vs PIV vs TrackMate in the buffer
%  ========================================================================
if DO_COMPARE_IND_ALL || DO_COMPARE_IND_PIV || DO_COMPARE_IND_TM

    % ---------- radial growth per crystal ----------
    rf = fullfile(resultsFolder, ['radial_velocity_' dropletName '.mat']);
    assert(isfile(rf), '%s is missing (per-crystal radial stage).', rf);
    Rc = load(rf, 'vRadial');
    vRadial = Rc.vRadial;
    nInd    = min(size(vRadial,1), nT);
    vRadial = vRadial(1:nInd, :);
    tGrowth = tAll(1:nInd);

    nC = min(NCrys, size(vRadial,2));    % Ncrys from specs, capped by disk
    if nC < NCrys
        warning('Ncrys=%d but vRadial has %d columns. Using %d.', ...
                NCrys, size(vRadial,2), nC);
    end

    % ---------- buffer rings ----------
    kf = fullfile(resultsFolder, ['individual_crystals_' dropletName '.mat']);
    assert(isfile(kf), '%s is missing (crystal buffers).', kf);
    K = load(kf, 'maskPlusBuffer');
    maskPlusBuffer = K.maskPlusBuffer;

    % ---------- PIV inside the buffer ----------
    vFlowPIV_ind = []; tPIV = [];
    pf = fullfile(resultsFolder, ['piv_fields_' dropletName '.mat']);
    if isfile(pf)
        P = load(pf, 'u_avg_um_s', 'v_avg_um_s', 'x_um', 'y_um');
        u_avg = P.u_avg_um_s;  v_avg = P.v_avg_um_s;
        xg = P.x_um;           yg = P.y_um;
        nP = min(numel(u_avg), nT);
        vFlowPIV_ind = nan(nP, nC);

        for k = 1:nP
            if isempty(u_avg{k}) || k > numel(xg) || isempty(xg{k}), continue; end
            Vmag = hypot(u_avg{k}, v_avg{k});
            if k > size(maskPlusBuffer,1), continue; end
            % PIV node coordinates back in pixels
            colP = round(xg{k} / um_per_px);
            rowP = round(yg{k} / um_per_px);
            for c = 1:nC
                buf = maskPlusBuffer{k, c};
                if isempty(buf) || ~any(buf(:)), continue; end
                % IMPORTANT: the PIV grid does not cover the whole image, so
                % the buffer is SAMPLED at each node's real coordinates.
                % Using imresize left it offset by half an interrogation
                % window, which put the wrong nodes inside the ring.
                inb = colP >= 1 & colP <= size(buf,2) & ...
                      rowP >= 1 & rowP <= size(buf,1);
                bufP = false(size(Vmag));
                bufP(inb) = buf(sub2ind(size(buf), rowP(inb), colP(inb)));
                if ~any(bufP(:)), continue; end
                vals = Vmag(bufP);
                vals = vals(isfinite(vals));
                if ~isempty(vals), vFlowPIV_ind(k,c) = mean(vals); end
            end
        end
        for c = 1:nC                              % same window as everything
            vc = vFlowPIV_ind(:,c);  vv = ~isnan(vc);
            vc = smoothdata(vc, 'movmean', SMOOTH_WINDOW, 'omitnan');
            vc(~vv) = NaN;
            vFlowPIV_ind(:,c) = vc;
        end
        tPIV = tAll(1:nP);
    end

    % ---------- TrackMate inside the buffer ----------
    vFlowTM_ind = []; tTMi = [];
    fM = fullfile(resultsFolder, ['trackmate_velocities_masked_' dropletName '.mat']);
    fP = fullfile(resultsFolder, ['trackmate_velocities_' dropletName '.mat']);
    fTV = '';
    if isfile(fM), fTV = fM; elseif isfile(fP), fTV = fP; end
    if ~isempty(fTV)
        TM = load(fTV, 'trackVelocityTable_reduced');
        Tv = TM.trackVelocityTable_reduced;
        if ~exist('flipTrackMateY','var'), flipTrackMateY = false; end

        nfTM = min([nInd, max(Tv.Frame)+1, size(maskPlusBuffer,1)]);
        vFlowTM_ind = nan(nfTM, nC);

        for fr0 = 0:(nfTM-1)
            frIdx = fr0 + 1;
            rowsF = Tv(Tv.Frame == fr0, :);
            if isempty(rowsF), continue; end
            for c = 1:nC
                buf = maskPlusBuffer{frIdx, c};
                if isempty(buf) || ~any(buf(:)), continue; end
                [Nyb, Nxb] = size(buf);
                col = round(rowsF.X_px);
                row = round(rowsF.Y_px);
                if flipTrackMateY, row = Nyb - row + 1; end
                ok  = row>=1 & row<=Nyb & col>=1 & col<=Nxb;
                if ~any(ok), continue; end
                inBuf = false(size(ok));
                inBuf(ok) = buf(sub2ind([Nyb Nxb], row(ok), col(ok)));
                spd = rowsF.Speed_um_s(inBuf);
                spd = spd(isfinite(spd));
                if ~isempty(spd), vFlowTM_ind(frIdx,c) = mean(spd); end
            end
        end
        for c = 1:nC
            vc = vFlowTM_ind(:,c);  vv = ~isnan(vc);
            vc = smoothdata(vc, 'movmean', SMOOTH_WINDOW, 'omitnan');
            vc(~vv) = NaN;
            vFlowTM_ind(:,c) = vc;
        end
        tTMi = tAll(1:nfTM);
    end

    % ---------- multi-panel figure, one panel per crystal ----------
    if ~exist('SAME_AXIS','var'),   SAME_AXIS   = true; end
    if ~exist('LW_PLOT','var'),     LW_PLOT     = 1.6;  end
    if ~exist('IND_PER_FIG','var'), IND_PER_FIG = 4;    end

    nFig = ceil(nC / IND_PER_FIG);

    for fi = 1:nFig
        cs = ((fi-1)*IND_PER_FIG + 1) : min(fi*IND_PER_FIG, nC);
        nP = numel(cs);

        nCols = min(nP, 2);
        nRows = ceil(nP / nCols);

        figM = figure('Name', sprintf('Crystals %d-%d: growth vs flow', cs(1), cs(end)), ...
                      'Color','w', 'Position',[60 60 560*nCols 460*nRows]);
        tl = tiledlayout(figM, nRows, nCols, ...
                         'TileSpacing','compact', 'Padding','loose');

        for jj = 1:nP
            c     = cs(jj);
            tNucC = nucTimes(min(c, numel(nucTimes)));

            vG = vRadial(:, c);
            vG(tGrowth < tNucC) = NaN;          % starts at ITS OWN nucleation

            ax = nexttile(tl);
            hold(ax,'on')
            hs = gobjects(0);  leg = {};

            if SAME_AXIS
                hs(end+1) = plot(ax, tGrowth, vG, '-', ...
                                 'LineWidth',LW_PLOT, 'Color',colorCrystal);
                leg{end+1} = 'Crystal $\dot{R}$';
                if ~isempty(vFlowPIV_ind)
                    hs(end+1) = plot(ax, tPIV, vFlowPIV_ind(:,c), '-', ...
                                     'LineWidth',LW_PLOT, 'Color',colorPIV);
                    leg{end+1} = 'PIV in buffer';
                end
                if ~isempty(vFlowTM_ind)
                    hs(end+1) = plot(ax, tTMi, vFlowTM_ind(:,c), '-', ...
                                     'LineWidth',LW_PLOT, 'Color',colorTM);
                    leg{end+1} = 'TrackMate in buffer';
                end
                ylabel(ax, 'Velocity ($\mu$m/s)');
            else
                yyaxis(ax,'left')
                hs(end+1) = plot(ax, tGrowth, vG, '-', ...
                                 'LineWidth',LW_PLOT, 'Color',colorCrystal);
                leg{end+1} = 'Crystal $\dot{R}$';
                ylabel(ax, 'Radial growth ($\mu$m/s)');
                set(ax,'YColor',[0.15 0.15 0.15]);
                yyaxis(ax,'right'); hold(ax,'on')
                if ~isempty(vFlowPIV_ind)
                    hs(end+1) = plot(ax, tPIV, vFlowPIV_ind(:,c), '-', ...
                                     'LineWidth',LW_PLOT, 'Color',colorPIV);
                    leg{end+1} = 'PIV in buffer';
                end
                if ~isempty(vFlowTM_ind)
                    hs(end+1) = plot(ax, tTMi, vFlowTM_ind(:,c), '-', ...
                                     'LineWidth',LW_PLOT, 'Color',colorTM);
                    leg{end+1} = 'TrackMate in buffer';
                end
                ylabel(ax, 'Flow velocity ($\mu$m/s)');
            end
            set(ax, 'YColor', [0.15 0.15 0.15]);

            title(ax, sprintf('Crystal %d  ($t_{\\mathrm{nuc}} = %.1f$ s)', c, tNucC));
            xlabel(ax, 'Time (s)');

            % Legend only in the first panel, so it is not repeated four times
            if jj == 1
                legend(ax, hs, leg, 'Location','northeast');
            end

            % Unlike the global plot, each panel starts at ITS crystal's
            % nucleation rather than at t_start: there is nothing to show
            % before it exists, and including that stretch squashed the
            % curve against the right-hand edge.
            nucleationLines(ax, tNucC, [], {sprintf('N%d',c)});
            figureStyle(ax, [tNucC, tAll(end)], tNucC);
        end

        if nFig > 1
            ttl = sprintf('%s: individual crystals %d--%d of %d', ...
                          strrep(dropletName,'_','\_'), cs(1), cs(end), nC);
        else
            ttl = sprintf('%s: growth vs local flow, crystal by crystal', ...
                          strrep(dropletName,'_','\_'));
        end
        title(tl, ttl, 'Interpreter','latex', 'FontSize',15, ...
              'FontName','Times New Roman');

        outName = fullfile(resultsFolder, ...
                  sprintf('compare_ind_panel%d_%s.png', fi, dropletName));
        exportgraphics(figM, outName, 'Resolution', 300);
        fprintf('Saved: compare_ind_panel%d_%s.png  (crystals %d-%d)\n', ...
                fi, dropletName, cs(1), cs(end));
    end

    % ---------- save the ALREADY SMOOTHED series ----------
    save(fullfile(resultsFolder, ['compare_individual_' dropletName '.mat']), ...
         'vRadial','vFlowPIV_ind','vFlowTM_ind','tGrowth','tPIV','tTMi', ...
         'nucTimes','t_start','fps','SMOOTH_WINDOW');
    fprintf('Individual comparisons saved (smoothed series).\n');
end
