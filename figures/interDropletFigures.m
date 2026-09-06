%% ========================================================================
%  INTERDROPLETFIGURES  What varies between droplets, and what does not.
%  ========================================================================
%
%  PRODUCES: FIGURES AND TABLES. Nothing is recomputed; every quantity is
%  read either from summary_crystal_growth.csv or from each droplet's own
%  RESULTS_crystals folder.
%
%  PREREQUISITE
%    crystalGrowthFigures must have been run first: it writes
%    summary_crystal_growth.csv, which is the starting point here.
%
%  WHAT THIS ANSWERS
%    Each droplet is one experiment under slightly different conditions. This
%    script asks which of those conditions actually matter: does temperature
%    set the drying time, does concentration set the number of crystals, and
%    above all does a faster-growing crystal come with a faster flow.
%
%    The master table gathers every per-droplet quantity in one CSV, and the
%    scatter figures show the pairs worth looking at, each annotated with its
%    Pearson r, p-value and n so the reader can judge the strength directly
%    rather than by eye.
%
%  THE CONTROL FOR EXPERIMENTAL DAY
%    Correlations are computed on three subsets: all 26 droplets, only those
%    at 280 g/L, and the single day with the most droplets. The reason is
%    that droplets recorded on the same day share ambient conditions, so a
%    correlation across all of them could reflect a difference between days
%    rather than a relationship between droplets. If a trend survives within
%    one day and within one concentration, it is a real trend. Figure 17
%    shows the same thing visually, coloured by day.
%
%  TIME HANDLING
%    The specifications always win over the CSV for the first nucleation
%    time, and a mismatch beyond one frame raises a warning. fps is never
%    read from a .mat file: those can be stale, the specifications cannot.
%
%  ONLY EDIT THE CONFIGURATION SECTION.
%
%  OUTPUT  (in <FG>/Interdroplet)
%    summary_master.csv          every per-droplet quantity, one row each
%    summary_piv_velocity.csv, summary_salt_velocity.csv
%    summary_descriptive.csv, summary_by_day.csv, summary_by_salt.csv
%    summary_correlations.csv    r and rho for every pair, on three subsets
%    scatter_*.png               17 figures
%
%  See also CRYSTALGROWTHFIGURES, GLOBALFLOWFIGURES.

clear; close all; clc;

%% <<<  CONFIGURATION  >>>
DATA_ROOT = fullfile('C:', 'Data', 'FG');

csvFile = fullfile(DATA_ROOT, 'RESULTS_CRYSTALGROWTH', 'summary_crystal_growth.csv');
outDir  = fullfile(DATA_ROOT, 'Interdroplet');
if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('Created output folder: %s\n', outDir);
end

SMOOTH_WINDOW = 10;
SHOW_R = true;      % annotate Pearson r, p and n inside each panel

sessionsAll = { ...
    'TP2-B(2811)','S4' ; 'TP3(1202)','S5'  ; 'TP3(1202)','S6'  ; ...
    'TP4-D(0106)','S7' ; 'TP4-D(0106)','S8'; 'TP4-D(0106)','S9' ; ...
    'TP4-D(0106)','S10'; 'TP5-E(0806)','S12'; 'TP5-E(0806)','S13'; ...
    'TP5-E(0806)','S14'; 'TP5-E(0806)','S15'; 'TP5-E(0806)','S18'; ...
    'TP5-E(0806)','S19'; 'TP5-E(0806)','S20'; 'TP6-(2907)','S21' ; ...
    'TP6-(2907)','S22' ; 'TP6-(2907)','S23' ; 'TP6-(2907)','S24' ; ...
    'TP6-(2907)','S25' ; 'TP6-(2907)','S29' ; 'TP6-(2907)','S30' ; ...
    'TP6-(2907)','S31' ; 'TP6-(2907)','S32' ; 'TP6-(2907)','S33' ; ...
    'TP6-(2907)','S34' ; 'TP6-(2907)','S35'   ...
};

dropletDirs = cell(size(sessionsAll,1),1);
for i = 1:size(sessionsAll,1)
    dropletDirs{i} = fullfile(DATA_ROOT, sessionsAll{i,1}, sessionsAll{i,2});
end

% CSV column names, in case the file is ever edited by hand
COL_TNUC1     = 't_nuc1_s';        % absolute time of the first nucleation
COL_TEND      = 't_end_s';         % last frame, i.e. end of evaporation
COL_NCRYS     = 'N_crystals';
COL_DREQ_MEAN = 'dReq_mean_umps';
COL_DREQ_MAX  = 'dReq_max_umps';

FONT_SZ = 12;
TICK_SZ = 10;
PT_COLOR = [0.180 0.459 0.710];    % single colour for every point

% ---- salt concentration per droplet (g/L) ----
%  Column 1: droplet number. Column 2: concentration. This is a property of
%  the solution prepared, not of the recording, so it is not in the specs.
SALT_TABLE = [ ...
     1 200;  2 200;  3 150;  4 150;  5 280;  6 280;  7 280;
     8 280;  9 280; 10 280; 11 280; 12 280; 13 280; 14 280;
    15 280; 16 280; 17 280; 18 280; 19 280; 20 280; 21 200;
    22 150; 23 280; 24 280; 25 280; 26 150; 27 200; 28 200;
    29 200; 30 200; 31 280; 32 150; 33 280; 34 280; 35 280];

%% ------------------------------------------------------------------------
%  READ THE CSV TABLE
%  ------------------------------------------------------------------------
T     = readtable(csvFile);
nD    = height(T);

t_nuc1  = T.(COL_TNUC1);
t_end   = T.(COL_TEND);
N_crys  = T.(COL_NCRYS);
dR_mean = T.(COL_DREQ_MEAN);
dR_max  = T.(COL_DREQ_MAX);
names   = cellstr(T.Session);

% Numeric label from the session name: 'S12' -> '12'. Used inside the markers.
labels  = regexprep(names, '[^0-9]', '');

% ---- dA/dt from the CSV, if those columns exist ----
dAdt_mean = nan(nD,1);
dAdt_max  = nan(nD,1);
if ismember('dAdt_mean_um2ps', T.Properties.VariableNames)
    dAdt_mean = T.dAdt_mean_um2ps;
end
if ismember('dAdt_max_um2ps', T.Properties.VariableNames)
    dAdt_max  = T.dAdt_max_um2ps;
end

%% ------------------------------------------------------------------------
%  READ THE SPECIFICATIONS AND THE EXPERIMENTAL DAY
%  ------------------------------------------------------------------------
temp_v = nan(nD,1);
hum_v  = nan(nD,1);

% Timing parameters PER DROPLET: the specifications change from one to the next
fps_v    = nan(nD,1);
tstart_v = nan(nD,1);
tnuc_v   = nan(nD,1);      % first nucleation, in ABSOLUTE time

% Experimental day, taken from each droplet's parent folder
day_v = repmat({''}, nD, 1);

for i = 1:numel(dropletDirs)
    [~, sesName] = fileparts(dropletDirs{i});
    idx = find(strcmp(names, sesName), 1);
    if isempty(idx)
        warning('Session ''%s'' not found in the CSV.', sesName); continue
    end

    [parentPath, ~] = fileparts(char(dropletDirs{i}));
    [~, dayName]    = fileparts(parentPath);
    day_v{idx}      = dayName;

    try
        [~, sp] = evalc('readSpecs(dropletDirs{i})');   % quiet console
        temp_v(idx) = sp.temperatura;
        hum_v(idx)  = sp.humedad;

        fps_v(idx)    = sp.fps;
        tstart_v(idx) = 0;
        if ~isnan(sp.t_start), tstart_v(idx) = sp.t_start; end

        % ---- first nucleation: THE SPECIFICATIONS WIN over the CSV ----
        %  Priority: nucleation_frames(1) > nucleation_frame >
        %            nucleation_time > CSV. Frames are preferred over the
        %            recorded time because the time can be a rounded value
        %            typed by hand, while the frame is exact.
        fn = NaN;
        if isfield(sp,'nucleation_frames') && ~isempty(sp.nucleation_frames) ...
                                           && any(~isnan(sp.nucleation_frames))
            fn = min(sp.nucleation_frames(~isnan(sp.nucleation_frames)));
        elseif ~isnan(sp.nucleation_frame)
            fn = sp.nucleation_frame;
        end
        if ~isnan(fn)
            tnuc_v(idx) = tstart_v(idx) + (fn - 1) / sp.fps;   % frame 1 = t_start
        elseif ~isnan(sp.nucleation_time)
            tnuc_v(idx) = sp.nucleation_time;
        else
            tnuc_v(idx) = t_nuc1(idx);                         % last resort: CSV
        end

        if ~isnan(t_nuc1(idx)) && abs(t_nuc1(idx) - tnuc_v(idx)) > 1/sp.fps
            warning(['%s: the CSV says t_nuc1=%.2f s but the specifications ' ...
                     'give %.2f s. The specifications are used.'], ...
                     sesName, t_nuc1(idx), tnuc_v(idx));
        end
    catch ME
        warning('Could not read specs for %s: %s', sesName, ME.message);
    end
end

%% ------------------------------------------------------------------------
%  LOAD THE PIV VELOCITIES
%  ------------------------------------------------------------------------
piv_mean_total  = nan(nD,1);
piv_mean_before = nan(nD,1);
piv_mean_after  = nan(nD,1);

for i = 1:numel(dropletDirs)
    [~, sesName] = fileparts(dropletDirs{i});
    idx = find(strcmp(names, sesName), 1);
    if isempty(idx), continue; end

    sesNameC = char(sesName);
    fileName = ['piv_fields_' sesNameC '.mat'];
    pivFile  = '';
    % Several candidate locations: the results folder moved between early and
    % late sessions, and this keeps the older droplets readable.
    candidates = {
        fullfile(char(dropletDirs{i}), 'RESULTS_crystals', fileName)
        fullfile(char(dropletDirs{i}), 'RESULTS_PIV',      fileName)
        fullfile(char(dropletDirs{i}),                     fileName)
        fullfile(fileparts(char(dropletDirs{i})), 'RESULTS_crystals', fileName)
    };
    for ci = 1:numel(candidates)
        if exist(candidates{ci}, 'file')
            pivFile = candidates{ci};
            break;
        end
    end
    if isempty(pivFile)
        warning('No piv_fields_%s.mat found - skipping PIV for %s.', ...
                sesNameC, sesNameC);
        continue
    end

    fprintf('PIV data: %s  (t_start=%.2f s, fps=%g, t_nuc=%.2f s)\n', ...
            pivFile, tstart_v(idx), fps_v(idx), tnuc_v(idx));

    try
        % NOTE: 'fps' is deliberately NOT read from the .mat, which may hold
        % a stale value. It always comes from THIS droplet's specifications.
        piv  = load(pivFile, 'velocity_magnitude_from_u_v_avg');
        vmag = piv.velocity_magnitude_from_u_v_avg;   % cell array {1 x Nframes}
        fpsP = fps_v(idx);
        assert(~isnan(fpsP) && fpsP > 0, 'No fps in the specifications of %s', sesNameC);
        Nf   = numel(vmag);

        % Spatial mean velocity (um/s) per frame
        vmean_t = nan(Nf,1);
        for k = 1:Nf
            V = vmag{k};
            if isempty(V), continue; end
            V = V(isfinite(V));
            if ~isempty(V), vmean_t(k) = mean(V); end
        end
        % Smooth without inventing the edges that have no data
        validV  = ~isnan(vmean_t);
        vmean_t = smoothdata(vmean_t, 'movmean', SMOOTH_WINDOW, 'omitnan');
        vmean_t(~validV) = NaN;

        % ---- ABSOLUTE TIME AXIS ----
        %  Frame 1 is t_start, not second 0, and t_nuc1 is absolute. Getting
        %  this wrong would put the before/after split in the wrong place.
        t_piv = tstart_v(idx) + (0 : Nf-1)' / fpsP;
        tnuc  = tnuc_v(idx);

        ok = ~isnan(vmean_t);
        if any(ok)
            piv_mean_total(idx) = mean(vmean_t(ok));
        end

        if ~isnan(tnuc)
            before = ok & (t_piv <  tnuc);
            after  = ok & (t_piv >= tnuc);
            if any(before), piv_mean_before(idx) = mean(vmean_t(before)); end
            if any(after),  piv_mean_after(idx)  = mean(vmean_t(after));  end
        end

    catch ME
        warning('Could not load PIV for %s: %s', sesName, ME.message);
    end
end

%% ------------------------------------------------------------------------
%  PRINT AND SAVE THE PIV TABLE
%  ------------------------------------------------------------------------
fprintf('\n%-8s  %16s  %17s  %16s\n', ...
        'Session', 'mean_PIV_total', 'mean_PIV_before', 'mean_PIV_after');
fprintf('%s\n', repmat('-',1,62));
for i = 1:nD
    fprintf('%-8s  %16.3f  %17.3f  %16.3f\n', ...
            names{i}, piv_mean_total(i), piv_mean_before(i), piv_mean_after(i));
end
fprintf('\n(units: um/s)\n\n');

pivTbl = table(names, piv_mean_total, piv_mean_before, piv_mean_after, ...
               'VariableNames', {'Session','mean_PIV_total_umps', ...
                                 'mean_PIV_before_nuc_umps','mean_PIV_after_nuc_umps'});
pivCsvPath = fullfile(outDir, 'summary_piv_velocity.csv');
writetable(pivTbl, pivCsvPath);
fprintf('PIV table saved: %s\n', pivCsvPath);

%% ------------------------------------------------------------------------
%  LOAD THE FRONT VELOCITY AND THE SALT CONCENTRATION
%  ------------------------------------------------------------------------
%  frontVel_perim comes from growth_global_<session>.mat, already smoothed by
%  crystalGrowthAnalysis. It only exists after nucleation, since before that
%  there is no crystal perimeter, so its mean is inherently post-nucleation.

front_mean = nan(nD,1);
salt_v     = nan(nD,1);

for i = 1:nD
    nnum = str2double(labels{i});
    if isnan(nnum), continue; end
    k = find(SALT_TABLE(:,1) == nnum, 1);
    if ~isempty(k), salt_v(i) = SALT_TABLE(k,2); end
end

for i = 1:numel(dropletDirs)
    [~, sesName] = fileparts(dropletDirs{i});
    idx = find(strcmp(names, sesName), 1);
    if isempty(idx), continue; end
    sesNameC = char(sesName);

    fileName = ['growth_global_' sesNameC '.mat'];
    gFile = '';
    candidates = {
        fullfile(char(dropletDirs{i}), 'RESULTS_crystals', fileName)
        fullfile(char(dropletDirs{i}),                     fileName)
        fullfile(fileparts(char(dropletDirs{i})), 'RESULTS_crystals', fileName)
    };
    for ci = 1:numel(candidates)
        if exist(candidates{ci}, 'file'), gFile = candidates{ci}; break; end
    end
    if isempty(gFile)
        warning('No growth_global_%s.mat found - skipping front velocity.', sesNameC);
        continue
    end

    try
        G  = load(gFile, 'frontVel_perim');
        fv = G.frontVel_perim(:);
        fv = fv(isfinite(fv));
        if ~isempty(fv), front_mean(idx) = mean(fv); end
    catch ME
        warning('Could not load front velocity for %s: %s', sesNameC, ME.message);
    end
end

%% ------------------------------------------------------------------------
%  LOAD THE TRACKMATE MEAN SPEED
%  ------------------------------------------------------------------------
%  meanSpeed comes from trackmate_meanspeed_[masked_]<session>.mat, already
%  smoothed by the pipeline. The time axis is rebuilt from the specifications
%  (this droplet's t_start and fps), never from the .mat.

tm_mean_total  = nan(nD,1);
tm_mean_before = nan(nD,1);
tm_mean_after  = nan(nD,1);

for i = 1:numel(dropletDirs)
    [~, sesName] = fileparts(dropletDirs{i});
    idx = find(strcmp(names, sesName), 1);
    if isempty(idx), continue; end
    sesNameC = char(sesName);

    resF  = fullfile(char(dropletDirs{i}), 'RESULTS_crystals');
    tmAll = dir(fullfile(resF, 'trackmate_meanspeed_*.mat'));
    if isempty(tmAll)
        warning('No trackmate_meanspeed_*.mat for %s - skipping.', sesNameC);
        continue
    end
    isM = contains({tmAll.name}, 'masked');
    if any(isM), k = find(isM,1); else, k = 1; end
    tmFile = fullfile(resF, tmAll(k).name);

    try
        M    = load(tmFile, 'framesU', 'meanSpeed');
        fpsP = fps_v(idx);
        assert(~isnan(fpsP) && fpsP > 0, 'No fps in the specifications of %s', sesNameC);

        v_t  = M.meanSpeed(:);
        t_tm = tstart_v(idx) + double(M.framesU(:)) / fpsP;   % frame 0 = t_start
        tnuc = tnuc_v(idx);

        ok = isfinite(v_t);
        if any(ok), tm_mean_total(idx) = mean(v_t(ok)); end
        if ~isnan(tnuc)
            before = ok & (t_tm <  tnuc);
            after  = ok & (t_tm >= tnuc);
            if any(before), tm_mean_before(idx) = mean(v_t(before)); end
            if any(after),  tm_mean_after(idx)  = mean(v_t(after));  end
        end
    catch ME
        warning('Could not load TrackMate for %s: %s', sesNameC, ME.message);
    end
end

fprintf('\n%-8s  %9s  %13s  %13s  %15s\n', ...
        'Session', 'salt_g/L', 'mean_PIV', 'mean_TM', 'mean_front');
fprintf('%s\n', repmat('-',1,64));
for i = 1:nD
    fprintf('%-8s  %9.0f  %13.3f  %13.3f  %15.4f\n', ...
            names{i}, salt_v(i), piv_mean_total(i), tm_mean_total(i), front_mean(i));
end

saltTbl = table(names, salt_v, piv_mean_total, tm_mean_total, front_mean, ...
                'VariableNames', {'Session','salt_gL','mean_PIV_total_umps', ...
                                  'mean_TM_total_umps','mean_frontVel_umps'});
saltCsvPath = fullfile(outDir, 'summary_salt_velocity.csv');
writetable(saltTbl, saltCsvPath);
fprintf('Salt table saved: %s\n\n', saltCsvPath);

%% ========================================================================
%  MASTER TABLE: every per-droplet quantity in one CSV
%  ========================================================================
droplet_num  = str2double(labels);
nuc_frac     = t_nuc1 ./ t_end;                     % relative phase of nucleation
tm_minus_piv = tm_mean_total - piv_mean_total;      % method difference
piv_ratio    = piv_mean_after ./ piv_mean_before;   % amplification at nucleation

MASTER = table( ...
    names, droplet_num, day_v, salt_v, temp_v, hum_v, ...
    t_nuc1, t_end, nuc_frac, N_crys, ...
    front_mean, dR_mean, dR_max, dAdt_mean, dAdt_max, ...
    piv_mean_total, piv_mean_before, piv_mean_after, piv_ratio, ...
    tm_mean_total, tm_mean_before, tm_mean_after, tm_minus_piv, ...
    fps_v, tstart_v, tnuc_v, ...
    'VariableNames', { ...
      'Session','Droplet','Day','salt_gL','temp_C','RH_pct', ...
      't_nuc1_s','t_evap_s','nuc_frac','N_crystals', ...
      'front_mean_umps','dReq_mean_umps','dReq_max_umps', ...
      'dAdt_mean_um2ps','dAdt_max_um2ps', ...
      'PIV_total_umps','PIV_before_umps','PIV_after_umps','PIV_ratio', ...
      'TM_total_umps','TM_before_umps','TM_after_umps','TM_minus_PIV_umps', ...
      'fps','t_start_s','t_nuc_specs_s'});

MASTER = sortrows(MASTER, 'Droplet');

masterPath = fullfile(outDir, 'summary_master.csv');
writetable(MASTER, masterPath);
fprintf('Master table saved: %s  (%d droplets, %d columns)\n', ...
        masterPath, height(MASTER), width(MASTER));

fprintf('\n=== MASTER TABLE (key columns) ===\n');
disp(MASTER(:, {'Session','Day','salt_gL','temp_C','RH_pct','t_nuc1_s', ...
                't_evap_s','nuc_frac','N_crystals','front_mean_umps', ...
                'PIV_total_umps','PIV_after_umps','TM_total_umps'}));

% ---- descriptive statistics of every numeric variable ----
% The coefficient of variation is the useful column: it shows at a glance
% which quantities are reproducible across droplets and which are not.
numVars = MASTER.Properties.VariableNames( ...
              varfun(@isnumeric, MASTER, 'OutputFormat','uniform'));
statRows = {};
fprintf('\n%-22s %5s %10s %10s %10s %10s %8s\n', ...
        'Variable','n','mean','sd','min','max','CV(%)');
fprintf('%s\n', repmat('-',1,80));
for k = 1:numel(numVars)
    v = MASTER.(numVars{k});
    v = v(isfinite(v));
    if isempty(v), continue; end
    cv = 100*std(v)/abs(mean(v));
    fprintf('%-22s %5d %10.3f %10.3f %10.3f %10.3f %8.1f\n', ...
            numVars{k}, numel(v), mean(v), std(v), min(v), max(v), cv);
    statRows(end+1,:) = {numVars{k}, numel(v), mean(v), std(v), ...
                         min(v), max(v), cv}; %#ok<SAGROW>
end
STATTBL = cell2table(statRows, 'VariableNames', ...
    {'Variable','n','mean','sd','min','max','CV_pct'});
writetable(STATTBL, fullfile(outDir,'summary_descriptive.csv'));
fprintf('Descriptive table saved: summary_descriptive.csv\n');

% ---- group means, by experimental day and by concentration ----
gvars = {'temp_C','RH_pct','t_nuc1_s','t_evap_s','nuc_frac','N_crystals', ...
         'front_mean_umps','PIV_total_umps','PIV_after_umps','TM_total_umps'};

fprintf('\n=== BY EXPERIMENTAL DAY ===\n');
BYDAY = groupsummary(MASTER, 'Day', {'mean','std'}, gvars);
disp(BYDAY);
writetable(BYDAY, fullfile(outDir,'summary_by_day.csv'));

fprintf('\n=== BY SALT CONCENTRATION ===\n');
BYSALT = groupsummary(MASTER, 'salt_gL', {'mean','std'}, gvars);
disp(BYSALT);
writetable(BYSALT, fullfile(outDir,'summary_by_salt.csv'));

%% ========================================================================
%  CORRELATION TABLE: Pearson r and Spearman rho
%  ========================================================================
%  Computed on three subsets:
%    ALL      the 26 droplets, which is what goes in the main text
%    280 g/L  removes concentration as a variable
%    <day>    the day with the most droplets, which separates a real trend
%             between droplets from a mere difference between days
%
%  Spearman is reported alongside Pearson because several of these
%  relationships are monotonic but not linear, and n = 26 is small enough
%  that a single outlier can move r substantially.

PAIRS = { ...
  % --- 1. Ambient conditions ---
  'Ambient','temp_C','t_evap_s'
  'Ambient','RH_pct','t_evap_s'
  'Ambient','temp_C','t_nuc1_s'
  'Ambient','RH_pct','t_nuc1_s'
  'Ambient','temp_C','RH_pct'
  'Ambient','t_nuc1_s','t_evap_s'
  'Ambient','temp_C','PIV_total_umps'
  'Ambient','temp_C','PIV_after_umps'
  'Ambient','RH_pct','PIV_total_umps'
  'Ambient','RH_pct','PIV_after_umps'
  'Ambient','temp_C','front_mean_umps'
  'Ambient','RH_pct','front_mean_umps'
  'Ambient','temp_C','N_crystals'
  'Ambient','temp_C','nuc_frac'
  % --- 2. Composition ---
  'Composition','salt_gL','t_nuc1_s'
  'Composition','salt_gL','N_crystals'
  'Composition','salt_gL','front_mean_umps'
  'Composition','salt_gL','PIV_total_umps'
  'Composition','salt_gL','TM_total_umps'
  'Composition','N_crystals','front_mean_umps'
  'Composition','N_crystals','dReq_max_umps'
  'Composition','N_crystals','dAdt_max_um2ps'
  'Composition','N_crystals','t_evap_s'
  'Composition','N_crystals','PIV_after_umps'
  % --- 3. Flow-growth coupling ---
  'Coupling','front_mean_umps','PIV_total_umps'
  'Coupling','front_mean_umps','PIV_after_umps'
  'Coupling','front_mean_umps','TM_total_umps'
  'Coupling','dReq_max_umps','PIV_after_umps'
  'Coupling','dAdt_max_um2ps','PIV_after_umps'
  'Coupling','dAdt_mean_um2ps','PIV_after_umps'
  'Coupling','dReq_max_umps','TM_after_umps'
  'Coupling','t_nuc1_s','PIV_after_umps'
  'Coupling','front_mean_umps','dReq_max_umps'
  };

[dayNames,~,dayIdx] = unique(MASTER.Day);
dayCounts = accumarray(dayIdx,1);
[~,kBig]  = max(dayCounts);
bigDay    = dayNames{kBig};

SUBSETS = { ...
  'ALL',      true(height(MASTER),1)
  '280 g/L',  MASTER.salt_gL == 280
  bigDay,     strcmp(MASTER.Day, bigDay)
  };

corrRows = {};
for s = 1:size(SUBSETS,1)
    sel   = SUBSETS{s,2};
    sname = SUBSETS{s,1};
    fprintf('\n=== CORRELATIONS - subset: %s (n = %d) ===\n', sname, sum(sel));
    fprintf('%-12s %-20s %-20s %4s %8s %9s %8s %9s\n', ...
            'Block','X','Y','n','r','p(r)','rho','p(rho)');
    fprintf('%s\n', repmat('-',1,96));
    for k = 1:size(PAIRS,1)
        blk = PAIRS{k,1};  xn = PAIRS{k,2};  yn = PAIRS{k,3};
        x = MASTER.(xn)(sel);   y = MASTER.(yn)(sel);
        ok = isfinite(x) & isfinite(y);
        n  = sum(ok);
        if n < 4, continue; end     % below four points a correlation is noise
        [r,   pr  ] = pearsonRP(x(ok), y(ok));
        % Spearman is Pearson computed on the ranks.
        [rho, prho] = pearsonRP(rankTies(x(ok)), rankTies(y(ok)));
        fprintf('%-12s %-20s %-20s %4d %8.3f %9.4f %8.3f %9.4f\n', ...
                blk, xn, yn, n, r, pr, rho, prho);
        corrRows(end+1,:) = {sname, blk, xn, yn, n, r, pr, rho, prho}; %#ok<SAGROW>
    end
end

CORRTBL = cell2table(corrRows, 'VariableNames', ...
    {'Subset','Block','X','Y','n','pearson_r','pearson_p','spearman_rho','spearman_p'});
corrPath = fullfile(outDir, 'summary_correlations.csv');
writetable(CORRTBL, corrPath);
fprintf('\nCorrelation table saved: %s\n\n', corrPath);

%% ========================================================================
%  FIGURE 1: temperature and humidity vs evaporation time
%  ========================================================================
fig1 = figure('Color','w','Units','inches','Position',[1 1 12 4.6]);

subplot(1,2,1);
scatterLabeled(t_end, temp_v, labels, PT_COLOR, ...
    '$t_\mathrm{evap}$ (s)', 'Temperature ($^\circ$C)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

subplot(1,2,2);
scatterLabeled(t_end, hum_v, labels, PT_COLOR, ...
    '$t_\mathrm{evap}$ (s)', 'Relative humidity (\%)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

addPanelLabels(fig1);
exportgraphics(fig1, fullfile(outDir,'scatter_evap_time.png'), 'Resolution',300);
fprintf('Saved: scatter_evap_time.png\n');

%% ========================================================================
%  FIGURE 2: front velocity and peak growth rate vs number of crystals
%  ========================================================================
fig2 = figure('Color','w','Position',[100 100 1200 470]);

subplot(1,2,1);
scatterLabeled(N_crys, front_mean, labels, PT_COLOR, ...
    'Number of crystals $N_\mathrm{crys}$', ...
    'Mean front velocity ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

subplot(1,2,2);
scatterLabeled(N_crys, dR_max, labels, PT_COLOR, ...
    'Number of crystals $N_\mathrm{crys}$', ...
    'Max $\dot{R}_\mathrm{eq}$ ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

addPanelLabels(fig2);
exportgraphics(fig2, fullfile(outDir,'scatter_growthrate_vs_ncrys.png'), 'Resolution',300);
fprintf('Saved: scatter_growthrate_vs_ncrys.png\n');

%% ========================================================================
%  FIGURE 3: temperature and humidity vs first nucleation time
%  ========================================================================
fig3 = figure('Color','w','Units','inches','Position',[14 1 12 4.6]);

subplot(1,2,1);
scatterLabeled(t_nuc1, temp_v, labels, PT_COLOR, ...
    '$t_{\mathrm{nuc},1}$ (s)', 'Temperature ($^\circ$C)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

subplot(1,2,2);
scatterLabeled(t_nuc1, hum_v, labels, PT_COLOR, ...
    '$t_{\mathrm{nuc},1}$ (s)', 'Relative humidity (\%)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

addPanelLabels(fig3);
exportgraphics(fig3, fullfile(outDir,'scatter_nucleation_time.png'), 'Resolution',300);
fprintf('Saved: scatter_nucleation_time.png\n');

%% ========================================================================
%  FIGURE 4: crystals vs evaporation time, and t_evap vs t_nuc1
%  ========================================================================
fig4 = figure('Color','w','Units','inches','Position',[14 7 12 4.6]);

ax4a = subplot(1,2,1);
scatterLabeled(t_end, N_crys, labels, PT_COLOR, ...
    '$t_\mathrm{evap}$ (s)', 'Number of crystals $N_\mathrm{crys}$', ...
    FONT_SZ, TICK_SZ, SHOW_R);
set(ax4a, 'YTick', unique(N_crys));   % integer counts: no fractional ticks

subplot(1,2,2);
scatterLabeled(t_nuc1, t_end, labels, PT_COLOR, ...
    '$t_{\mathrm{nuc},1}$ (s)', '$t_\mathrm{evap}$ (s)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

addPanelLabels(fig4);
exportgraphics(fig4, fullfile(outDir,'scatter_bonus.png'), 'Resolution',300);
fprintf('Saved: scatter_bonus.png\n');

%% ========================================================================
%  FIGURE 5: total mean PIV velocity vs ambient conditions
%  ========================================================================
fig5 = figure('Color','w','Units','inches','Position',[1 1 12 4.6]);

subplot(1,2,1);
scatterLabeled(temp_v, piv_mean_total, labels, PT_COLOR, ...
    'Temperature ($^\circ$C)', 'Mean PIV velocity ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

subplot(1,2,2);
scatterLabeled(hum_v, piv_mean_total, labels, PT_COLOR, ...
    'Relative humidity (\%)', 'Mean PIV velocity ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

addPanelLabels(fig5);
exportgraphics(fig5, fullfile(outDir,'scatter_piv_total.png'), 'Resolution',300);
fprintf('Saved: scatter_piv_total.png\n');

%% ========================================================================
%  FIGURE 6: post-nucleation PIV velocity vs ambient conditions
%  ========================================================================
fig6 = figure('Color','w','Units','inches','Position',[14 1 12 4.6]);

subplot(1,2,1);
scatterLabeled(temp_v, piv_mean_after, labels, PT_COLOR, ...
    'Temperature ($^\circ$C)', ...
    'Mean PIV velocity after nuc.\ ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

subplot(1,2,2);
scatterLabeled(hum_v, piv_mean_after, labels, PT_COLOR, ...
    'Relative humidity (\%)', ...
    'Mean PIV velocity after nuc.\ ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

addPanelLabels(fig6);
exportgraphics(fig6, fullfile(outDir,'scatter_piv_after_nuc.png'), 'Resolution',300);
fprintf('Saved: scatter_piv_after_nuc.png\n');

%% ========================================================================
%  FIGURE 7: pre-nucleation PIV velocity vs ambient conditions
%  ========================================================================
%  Only 20 droplets have a pre-nucleation interval, so this is a candidate
%  for the appendix rather than the main text.
fig7 = figure('Color','w','Units','inches','Position',[1 7 12 4.6]);

subplot(1,2,1);
scatterLabeled(temp_v, piv_mean_before, labels, PT_COLOR, ...
    'Temperature ($^\circ$C)', ...
    'Mean PIV velocity before nuc.\ ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

subplot(1,2,2);
scatterLabeled(hum_v, piv_mean_before, labels, PT_COLOR, ...
    'Relative humidity (\%)', ...
    'Mean PIV velocity before nuc.\ ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

addPanelLabels(fig7);
exportgraphics(fig7, fullfile(outDir,'scatter_piv_before_nuc.png'), 'Resolution',300);
fprintf('Saved: scatter_piv_before_nuc.png\n');

%% ========================================================================
%  FIGURE 8: PIV velocity vs salt concentration
%  ========================================================================
fig8 = figure('Color','w','Position',[100 100 700 560]);
scatterLabeled(salt_v, piv_mean_total, labels, PT_COLOR, ...
    'Salt concentration (g/L)', ...
    'Mean PIV velocity ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);
exportgraphics(fig8, fullfile(outDir,'scatter_piv_vs_salt.png'), 'Resolution',300);
fprintf('Saved: scatter_piv_vs_salt.png\n');

%% ========================================================================
%  FIGURE 9: crystal front velocity vs salt concentration
%  ========================================================================
fig9 = figure('Color','w','Position',[100 100 700 560]);
scatterLabeled(salt_v, front_mean, labels, PT_COLOR, ...
    'Salt concentration (g/L)', ...
    'Mean front velocity ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);
exportgraphics(fig9, fullfile(outDir,'scatter_front_vs_salt.png'), 'Resolution',300);
fprintf('Saved: scatter_front_vs_salt.png\n');

%% ========================================================================
%  FIGURE 10: PIV velocity vs crystal front velocity
%  ========================================================================
%  The coupling figure, over the whole recording.
fig10 = figure('Color','w','Position',[100 100 700 560]);
scatterLabeled(front_mean, piv_mean_total, labels, PT_COLOR, ...
    'Mean front velocity ($\mu$m s$^{-1}$)', ...
    'Mean PIV velocity ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);
exportgraphics(fig10, fullfile(outDir,'scatter_piv_vs_front.png'), 'Resolution',300);
fprintf('Saved: scatter_piv_vs_front.png\n');

%% ========================================================================
%  FIGURE 11: number of crystals and nucleation time vs salt concentration
%  ========================================================================
fig11 = figure('Color','w','Position',[100 100 1200 470]);

subplot(1,2,1);
scatterLabeled(salt_v, N_crys, labels, PT_COLOR, ...
    'Salt concentration (g/L)', ...
    'Number of crystals $N_\mathrm{crys}$', ...
    FONT_SZ, TICK_SZ, SHOW_R);

subplot(1,2,2);
scatterLabeled(salt_v, t_nuc1, labels, PT_COLOR, ...
    'Salt concentration (g/L)', ...
    '$t_{\mathrm{nuc},1}$ (s)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

addPanelLabels(fig11);
exportgraphics(fig11, fullfile(outDir,'scatter_nucleation_vs_salt.png'), 'Resolution',300);
fprintf('Saved: scatter_nucleation_vs_salt.png\n');

%% ========================================================================
%  FIGURE 12: front velocity vs ambient conditions
%  ========================================================================
fig12 = figure('Color','w','Position',[100 100 1200 470]);

subplot(1,2,1);
scatterLabeled(temp_v, front_mean, labels, PT_COLOR, ...
    'Temperature ($^\circ$C)', ...
    'Mean front velocity ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

subplot(1,2,2);
scatterLabeled(hum_v, front_mean, labels, PT_COLOR, ...
    'Relative humidity (\%)', ...
    'Mean front velocity ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

addPanelLabels(fig12);
exportgraphics(fig12, fullfile(outDir,'scatter_front_vs_ambient.png'), 'Resolution',300);
fprintf('Saved: scatter_front_vs_ambient.png\n');

%% ========================================================================
%  FIGURE 13: flow vs growth, post-nucleation only
%  ========================================================================
%  Both quantities restricted to the same time window. Figure 10 averages the
%  flow over the whole recording, including the stretch before any crystal
%  exists, which dilutes the coupling. This is the fair comparison.
fig13 = figure('Color','w','Position',[100 100 700 560]);
scatterLabeled(front_mean, piv_mean_after, labels, PT_COLOR, ...
    'Mean front velocity ($\mu$m s$^{-1}$)', ...
    'Mean PIV velocity after nucleation ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);
exportgraphics(fig13, fullfile(outDir,'scatter_pivafter_vs_front.png'), 'Resolution',300);
fprintf('Saved: scatter_pivafter_vs_front.png\n');

%% ========================================================================
%  FIGURE 14: the same check with TrackMate
%  ========================================================================
%  Confirms the conclusion does not depend on the measurement method.
fig14 = figure('Color','w','Position',[100 100 1200 470]);

subplot(1,2,1);
scatterLabeled(salt_v, tm_mean_total, labels, PT_COLOR, ...
    'Salt concentration (g/L)', ...
    'Mean TrackMate velocity ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

subplot(1,2,2);
scatterLabeled(front_mean, tm_mean_total, labels, PT_COLOR, ...
    'Mean front velocity ($\mu$m s$^{-1}$)', ...
    'Mean TrackMate velocity ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

addPanelLabels(fig14);
exportgraphics(fig14, fullfile(outDir,'scatter_trackmate_check.png'), 'Resolution',300);
fprintf('Saved: scatter_trackmate_check.png\n');

%% ========================================================================
%  FIGURE 15: coupling with the growth RATE rather than the front velocity
%  ========================================================================
%  The flow scales better with dReq_max and dA/dt than with the front
%  velocity normalised by perimeter, which is worth showing explicitly: the
%  perimeter normalisation divides out part of the very signal being sought.
fig15 = figure('Color','w','Position',[100 100 1200 470]);

subplot(1,2,1);
scatterLabeled(dR_max, piv_mean_after, labels, PT_COLOR, ...
    'Max $\dot{R}_\mathrm{eq}$ ($\mu$m s$^{-1}$)', ...
    'Mean PIV velocity after nuc.\ ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

subplot(1,2,2);
scatterLabeled(dAdt_max, piv_mean_after, labels, PT_COLOR, ...
    'Max $\dot{A}$ ($\mu$m$^2$ s$^{-1}$)', ...
    'Mean PIV velocity after nuc.\ ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ, SHOW_R);

addPanelLabels(fig15);
exportgraphics(fig15, fullfile(outDir,'scatter_piv_vs_growthrate.png'), 'Resolution',300);
fprintf('Saved: scatter_piv_vs_growthrate.png\n');

%% ========================================================================
%  FIGURE 16: the relative phase of nucleation
%  ========================================================================
%  t_nuc1/t_evap is nearly constant across droplets: the ambient conditions
%  rescale the droplet's clock without changing WHEN it nucleates in relative
%  terms. That is why the aligned figures elsewhere work at all.
fig16 = figure('Color','w','Position',[100 100 1200 470]);

ax16a = subplot(1,2,1);
scatterLabeled(t_end, t_nuc1, labels, PT_COLOR, ...
    '$t_\mathrm{evap}$ (s)', '$t_{\mathrm{nuc},1}$ (s)', ...
    FONT_SZ, TICK_SZ, SHOW_R);
hold(ax16a,'on');
mf = mean(nuc_frac(isfinite(nuc_frac)));
xr = [min(t_end) max(t_end)];
plot(ax16a, xr, mf*xr, 'k--', 'LineWidth', 1.2);   % constant-fraction line
text(ax16a, 0.97, 0.06, sprintf('slope $= %.2f$', mf), ...
     'Units','normalized', 'HorizontalAlignment','right', ...
     'Interpreter','latex', 'FontSize', TICK_SZ);
hold(ax16a,'off');

ax16b = subplot(1,2,2);
scatterLabeled(t_end, nuc_frac, labels, PT_COLOR, ...
    '$t_\mathrm{evap}$ (s)', '$t_{\mathrm{nuc},1}/t_\mathrm{evap}$', ...
    FONT_SZ, TICK_SZ, SHOW_R);
hold(ax16b,'on');
yline(ax16b, mf, 'k--', 'LineWidth', 1.2);
ylim(ax16b, [0 1]);      % it is a fraction: the full range gives the scale
hold(ax16b,'off');

addPanelLabels(fig16);
exportgraphics(fig16, fullfile(outDir,'scatter_nucleation_phase.png'), 'Resolution',300);
fprintf('Saved: scatter_nucleation_phase.png\n');

%% ========================================================================
%  FIGURE 17: diagnostic, the coupling coloured by experimental day
%  ========================================================================
%  Shows whether the trend exists WITHIN each day or only between days. This
%  is a control figure meant for the appendix, and it deliberately does not
%  follow the PIV/TrackMate/crystal colour code: its colours mean days.
fig17 = figure('Color','w','Position',[100 100 1200 470]);

subplot(1,2,1);
scatterByGroup(front_mean, piv_mean_after, day_v, labels, ...
    'Mean front velocity ($\mu$m s$^{-1}$)', ...
    'Mean PIV velocity after nuc.\ ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ);

subplot(1,2,2);
scatterByGroup(temp_v, piv_mean_after, day_v, labels, ...
    'Temperature ($^\circ$C)', ...
    'Mean PIV velocity after nuc.\ ($\mu$m s$^{-1}$)', ...
    FONT_SZ, TICK_SZ);

addPanelLabels(fig17);
exportgraphics(fig17, fullfile(outDir,'scatter_grouped_by_day.png'), 'Resolution',300);
fprintf('Saved: scatter_grouped_by_day.png\n');

fprintf('\nDone. All output in: %s\n', outDir);


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function scatterLabeled(xv, yv, labels, color, xlbl, ylbl, ...
                        FONT_SZ, TICK_SZ, showR)
%SCATTERLABELED  Single-colour scatter with the droplet number inside each point.
%
%   No legend and no title: with 26 points a legend is unusable, and the
%   number written inside the marker identifies each droplet without adding
%   any clutter. The title is omitted because it would duplicate the axis
%   labels and the LaTeX caption; panels are identified by (a), (b), ...
%
%   showR (optional): annotate Pearson r, its p-value and n, top left.

    if nargin < 9, showR = false; end

    ax = gca;
    hold(ax, 'on');

    darkEdge = color * 0.65;   % slightly darker edge, so points separate

    for k = 1:numel(xv)
        if isnan(xv(k)) || isnan(yv(k)), continue; end

        scatter(ax, xv(k), yv(k), 140, color, 'filled', ...
                'MarkerEdgeColor', darkEdge, 'LineWidth', 0.8);

        text(ax, xv(k), yv(k), labels{k}, ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment',   'middle', ...
             'FontSize',   6, ...
             'FontWeight', 'bold', ...
             'Color',      'w', ...
             'Interpreter','none');
    end

    xlabel(ax, xlbl, 'Interpreter','latex', 'FontSize', FONT_SZ);
    ylabel(ax, ylbl, 'Interpreter','latex', 'FontSize', FONT_SZ);

    grid(ax, 'on');
    ax.GridAlpha      = 0.35;
    ax.GridLineStyle  = '-';
    ax.MinorGridAlpha = 0.18;

    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', TICK_SZ, 'Box', 'on');

    if showR
        ok = isfinite(xv) & isfinite(yv);
        if sum(ok) >= 4
            [r, p] = pearsonRP(xv(ok), yv(ok));
            text(ax, 0.03, 0.96, ...
                 sprintf('$r = %.2f$, $p = %.3f$, $n = %d$', r, p, sum(ok)), ...
                 'Units','normalized', 'Interpreter','latex', ...
                 'FontSize', TICK_SZ, 'VerticalAlignment','top', ...
                 'BackgroundColor',[1 1 1], 'Margin', 2);
        end
    end
end


function scatterByGroup(xv, yv, groups, labels, xlbl, ylbl, ...
                        FONT_SZ, TICK_SZ)
%SCATTERBYGROUP  Scatter coloured by group, here the experimental day.
%   Diagnostic figure: it does NOT follow the PIV/TrackMate/crystal code,
%   because here colour carries a different meaning entirely.

    ax = gca;
    hold(ax, 'on');

    PAL = [0.121 0.467 0.706
           0.890 0.467 0.098
           0.173 0.627 0.173
           0.839 0.153 0.157
           0.580 0.404 0.741
           0.549 0.337 0.294];

    [gNames, ~, gIdx] = unique(groups);
    hLeg = gobjects(numel(gNames),1);

    for g = 1:numel(gNames)
        c   = PAL(mod(g-1, size(PAL,1)) + 1, :);
        sel = (gIdx == g) & isfinite(xv) & isfinite(yv);
        if ~any(sel), continue; end
        hLeg(g) = scatter(ax, xv(sel), yv(sel), 140, c, 'filled', ...
                          'MarkerEdgeColor', c*0.65, 'LineWidth', 0.8);
        ii = find(sel);
        for k = 1:numel(ii)
            text(ax, xv(ii(k)), yv(ii(k)), labels{ii(k)}, ...
                 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                 'FontSize', 6, 'FontWeight','bold', 'Color','w', ...
                 'Interpreter','none');
        end
    end

    xlabel(ax, xlbl, 'Interpreter','latex', 'FontSize', FONT_SZ);
    ylabel(ax, ylbl, 'Interpreter','latex', 'FontSize', FONT_SZ);

    grid(ax, 'on');
    ax.GridAlpha     = 0.35;
    ax.GridLineStyle = '-';
    set(ax, 'TickLabelInterpreter','latex', 'FontSize', TICK_SZ, 'Box','on');

    valid = isgraphics(hLeg);
    legend(ax, hLeg(valid), gNames(valid), 'Interpreter','none', ...
           'Location','best', 'FontSize', TICK_SZ-1, 'Box','off');
end


function addPanelLabels(figH)
%ADDPANELLABELS  Label the axes of a multi-panel figure (a), (b), (c), ...
%   Ordered top to bottom and left to right. Does nothing on a single-axes
%   figure, so it can be called unconditionally. The legend axes is excluded,
%   otherwise it would be labelled as a panel.

    axAll = findall(figH, 'Type', 'axes');
    axAll = axAll(~strcmp(get(axAll,'Tag'), 'legend'));
    if numel(axAll) < 2, return; end

    pos = cell2mat(get(axAll, 'Position'));
    [~, ord] = sortrows([-pos(:,2), pos(:,1)]);   % top row first
    axAll = axAll(ord);

    letters = 'abcdefghijkl';
    for k = 1:min(numel(axAll), numel(letters))
        text(axAll(k), -0.14, 1.06, sprintf('(%s)', letters(k)), ...
             'Units','normalized', 'Interpreter','latex', ...
             'FontSize', 13, 'FontWeight','bold', ...
             'HorizontalAlignment','left', 'VerticalAlignment','bottom');
    end
end


function [r, p] = pearsonRP(x, y)
%PEARSONRP  Pearson r and its two-tailed p-value, base MATLAB only.
%   Written out here so the whole analysis runs without the Statistics
%   Toolbox, which not every machine in the lab has.
    x = x(:);  y = y(:);
    n = numel(x);
    if n < 3
        r = NaN; p = NaN; return
    end
    C = corrcoef(x, y);
    r = C(1,2);
    if ~isfinite(r) || abs(r) >= 1
        p = NaN; return
    end
    df    = n - 2;
    tstat = r * sqrt(df / (1 - r^2));
    p     = betainc(df / (df + tstat^2), df/2, 0.5);   % Student t, two-tailed
end


function rk = rankTies(x)
%RANKTIES  Ranks with ties averaged, for Spearman.
%   Equivalent to tiedrank without the Statistics Toolbox.
    x = x(:);
    [xs, ord] = sort(x);
    n  = numel(x);
    rk = zeros(n,1);
    i  = 1;
    while i <= n
        j = i;
        while j < n && xs(j+1) == xs(i), j = j + 1; end
        rk(ord(i:j)) = (i + j) / 2;
        i = j + 1;
    end
end
