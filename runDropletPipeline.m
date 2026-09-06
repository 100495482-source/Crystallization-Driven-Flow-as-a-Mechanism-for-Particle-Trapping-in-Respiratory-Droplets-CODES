%% ========================================================================
%  RUNDROPLETPIPELINE  Entry point of the per-droplet analysis.
%  ========================================================================
%
%  WHAT THIS SCRIPT DOES
%    This is the only file that needs editing to analyse a droplet. It sets
%    every parameter of the pipeline, reads the acquisition specifications,
%    builds the canonical time axis, and then calls the analysis stages that
%    are switched on:
%
%      crystalGrowthAnalysis   crystal segmentation, areas, growth rates and
%                              front velocities, globally and per crystal
%      pivFlowAnalysis         Eulerian velocity fields from PIVlab exports
%      trackFlowAnalysis       Lagrangian particle tracks from TrackMate
%      growthFlowComparison    growth curves against flow curves
%
%    Each stage is a script, not a function, so they share this workspace.
%    That is why the variable names below matter: the stages read them
%    directly.
%
%  HOW TO USE IT
%    1. Point DATA_ROOT at the folder holding the acquisition sessions.
%    2. Set dropletName to the droplet you want.
%    3. Switch on the stages you need with the flags below.
%    4. Run.
%
%  PARTICLE ACCUMULATION
%    The accumulation analysis is NOT called from here. It runs across all
%    droplets at once rather than one at a time, so it lives in its own
%    script: analysis/particleAccumulationAnalysis.m.
%
%  See also READSPECS, SYNCTIMEAXIS, CRYSTALGROWTHANALYSIS, PIVFLOWANALYSIS,
%  TRACKFLOWANALYSIS, GROWTHFLOWCOMPARISON.

close all; clc

%% ========================================================================
%  WORKSPACE CLEANUP  (do not remove)
%  ========================================================================
%  The specifications change from droplet to droplet. Without clearing, the
%  previous droplet's data (masks, crystalMasks, nucList, frameTimes...)
%  stays in the workspace, and every block written as
%      if ~exist('crystalMasks','var') ... load ...
%  finds them already defined and REUSES the wrong droplet's data without
%  ever loading the right one. The failure is silent and contaminates all
%  downstream results.
clearvars

addpath(genpath(fileparts(mfilename('fullpath'))));

%% ========================================================================
%  1. CONFIGURATION
%  ========================================================================

% ---------- data location ----------
% Root folder containing the acquisition sessions. Change this to wherever
% the image sequences live on the machine running the analysis.
DATA_ROOT   = fullfile('C:', 'Data', 'FG');
sessionName = 'TP6-(2907)';
dropletName = 'S35';

dropletDir  = fullfile(DATA_ROOT, sessionName, dropletName);
dataDroplet = dropletName;   % stamp: which droplet the workspace belongs to

% ---------- global LaTeX defaults for every figure ----------
set(groot, 'defaultTextInterpreter',             'latex');
set(groot, 'defaultAxesTickLabelInterpreter',    'latex');
set(groot, 'defaultLegendInterpreter',           'latex');
set(groot, 'defaultColorbarTickLabelInterpreter','latex');
set(groot, 'defaultAxesFontName',                'Times New Roman');
set(groot, 'defaultTextFontName',                'Times New Roman');

%% ========================================================================
%  2. STAGE SELECTION
%  ========================================================================
%  Flags marked REQUIRED must have been run at least once for a given
%  droplet before the later stages can work, because they produce the .mat
%  files those stages load.

% ---------- crystal segmentation and growth ----------
DO_CRYSTAL_ANALYSIS          = false;
    DO_GLOBAL_MASKS          = false;  % build the crystal masks   [REQUIRED]
    DO_GLOBAL_VELOCITIES     = true;   % global growth + front velocities
    DO_INDIVIDUAL_MASKS      = true;   % per-crystal masks + buffer ring
    DO_INDIVIDUAL_VELOCITIES = true;   % per-crystal growth and radial velocity

% ---------- PIV (Eulerian) ----------
DO_PIV = false;
    DO_PIV_FIELDS       = false;  % read PIVlab export, calibrate, mask, average
    DO_PIV_1DPLOTS      = true;   % spatial mean/max/std per frame vs time
    DO_VIDEO_VELOCITY   = false;  % vector field video
    DO_VIDEO_DIVERGENCE = false;  % divergence video
    DO_VIDEO_VORTICITY  = false;  % vorticity video
    DO_VIDEO_MAGNITUDE  = false;  % velocity magnitude video

% ---------- TrackMate (Lagrangian) ----------
DO_TRACKMATE = false;
    DO_TRACKMATE_READ     = true;  % read CSV, filter tracks, remove the timer
    DO_TRACKMATE_VELOCITY = true;  % per-point velocities
    DO_TRACKMATE_MEANVEL  = true;  % mean speed per frame vs time
    DO_TRACKMATE_VELPLOT  = false; % tracks coloured by speed (figure or video)
    TM_REMOVE_INSIDE_MASK = true;  % drop spots falling inside a crystal.
                                   % KEEP THIS TRUE: crystal texture is
                                   % otherwise tracked as if it were particles.

% ---------- growth vs flow comparison ----------
DO_COMPARISON = true;
    DO_COMPARE_GLOBAL_PIV = true;  % global crystal growth vs PIV
    DO_COMPARE_GLOBAL_TM  = true;  % global crystal growth vs TrackMate
    DO_COMPARE_GLOBAL_ALL = true;  % all three curves together
    DO_COMPARE_IND_PIV    = true;  % per crystal: growth vs PIV in its buffer
    DO_COMPARE_IND_TM     = true;  % per crystal: growth vs TrackMate
    DO_COMPARE_IND_ALL    = true;  % per crystal: all three together

%% ========================================================================
%  3. SEGMENTATION PARAMETERS
%  ========================================================================
INVERT_IMAGE    = false;  % the FotosOtsu sequence is already inverted
OTSU_SCALE      = 0.9;    % multiplies the automatic Otsu threshold
GAUSS_SIGMA     = 3;      % pre-smoothing before thresholding
OTSU_MIN_SIZE   = 1;      % minimum object area kept (px)
MORPH_OPEN      = 20;     % removes thin protrusions and noise (0 = off)
MORPH_CLOSE     = 40;     % fills gaps and joins the crystal (0 = off)
CONVEX_STRENGTH = 15;     % light convex filling; 0 = off, 2-3 subtle, 10+ strong
SKIP_PRENUC     = true;   % empty mask before nucleation
USE_POLY_SMOOTH = false;  % replace the outline by a simplified polygon
POLY_TOL        = 0.01;   % polygon tolerance: higher = fewer vertices
MANUAL_LAST     = true;   % draw the crystals by hand on the last frame
SHOW_MASK_MOVIE = true;   % play the mask sequence after segmenting
DISPLAY_N       = 8;      % frames shown in the segmentation summary figure

%% ========================================================================
%  4. GROWTH AND FLOW PARAMETERS
%  ========================================================================

% ---------- crystal growth ----------
N_FINITE_DIFF = 4;   % half-width in frames of the centred derivative
SMOOTH_WINDOW = 10;  % moving-average window, SHARED by crystals, PIV and
                     % TrackMate so the three curves stay comparable
bufferFactor  = 3;   % buffer ring radius, in multiples of the crystal radius

% ---------- PIV ----------
USE_CRYSTAL_MASK = true;  % set velocities inside crystals to NaN
K_AVG            = 4;     % temporal averaging radius, frames [n-K_AVG, n+K_AVG]

% ---------- TrackMate: track filtering ----------
TM_MIN_LENGTH  = 5;      % minimum track length (points)
TM_MIN_QUALITY = 0;      % minimum TrackMate quality score
TM_VELOCITY_N  = 4;      % half-width of the centred difference, in points
flipTrackMateY = false;  % true if the TrackMate y axis is inverted

% ---------- TrackMate: on-screen timer removal ----------
%  The recordings carry a burnt-in stopwatch whose digits are detected as
%  spots. The box is drawn once per droplet and stored.
tm_remove_timer = true;
tm_timer_box    = [];    % [x1 y1 x2 y2] in px; [] = draw it by hand
tm_redraw_box   = true;  % true = redraw even if a stored box exists

% ---------- TrackMate: cosmetic video settings ----------
TM_TAIL       = 30;     % track tail length in the video
TM_MIN_POINTS = 10;     % minimum points before a tail is drawn
TM_SMOOTH_WIN = 3;      % track smoothing, display only
TM_NSHOW      = 4;      % frames in the summary figure
TM_MAKE_VIDEO = false;  % false if only the data are needed

TM_VELPLOT_MODE    = 'frames';  % 'video' or 'frames'
TM_VELPLOT_NFRAMES = 6;         % frames shown if mode is 'frames'
TM_VELPLOT_TAIL    = 15;        % tail length, frames
TM_VELPLOT_CLIM    = [];        % [] = automatic (98th percentile), or [0 30]

%% ========================================================================
%  5. PLOTTING CONVENTIONS
%  ========================================================================
%  ONE COLOUR CODE ACROSS THE WHOLE PIPELINE
%    GREEN = crystal      BLUE = PIV      ORANGE = TrackMate
%  These three never change, in any plot, in any chapter.
colorPIV     = [0.20 0.40 0.80];   % blue   -> PIV (Eulerian)
colorTM      = [0.85 0.40 0.20];   % orange -> TrackMate (Lagrangian)
colorCrystal = [0.10 0.60 0.30];   % green  -> crystals

LW_PLOT     = 1.6;    % same line width for all three curves
SAME_AXIS   = true;   % true  -> one y axis (all curves are in um/s)
                      % false -> double axis (crystal left, flow right)
IND_PER_FIG = 4;      % crystals per multi-panel comparison figure

%% ========================================================================
%  6. GLOBAL CRYSTAL METRIC USED IN THE COMPARISONS
%  ========================================================================
%  'perim'   frontVel_perim = (dA/dt)/P_total
%            Mean normal velocity of the interface, since dA/dt = int(v_n dl).
%            Independent of the number of crystals and of their shape, which
%            is what should be compared against the flow. Its weakness is
%            that the perimeter is the quantity most sensitive to
%            segmentation (MORPH_CLOSE and CONVEX_STRENGTH smooth it).
%
%  'Req'     frontVel_Req = d/dt sqrt(A_total/pi)
%            Collapses ALL crystals into a single equivalent circle. With N
%            crystals of radius r, A = N*pi*r^2 gives R_eq = r*sqrt(N), so
%            dR_eq/dt = sqrt(N)*dr/dt. The velocity is inflated by sqrt(N)
%            without any crystal growing faster, and jumps artificially at
%            every nucleation. Its only advantage is depending on area alone,
%            which is more robust to measure.
%
%  'radial'  vRadialGlobal = mean of the individual dR/dt
%            Radial like Req but averaged PER CRYSTAL, without the sqrt(N).
%            Usually the best compromise.
GLOBAL_CRYSTAL_METRIC = 'perim';

%% ========================================================================
%  7. READ THE DROPLET
%  ========================================================================
specs     = readSpecs(dropletDir);
fps       = specs.fps;
um_per_px = specs.ppf;
NCrys     = specs.Ncrys;
frame_max = specs.frame_max;
nucFrame  = specs.nucleation_frame;

% One colour per crystal, all within the green range, so that
% "green = crystal" still holds when several crystals share a figure.
% (lines() produced blues and oranges, which clashed with PIV and TrackMate.)
if isnan(NCrys) || NCrys < 1, NCrys = 1; end
crystalColors = [linspace(0.05,0.55,NCrys)', ...
                 linspace(0.45,0.80,NCrys)', ...
                 linspace(0.22,0.45,NCrys)'];

% ---------- locate the thresholded image folder (FotosOtsu*) ----------
dd = dir(dropletDir); dd = dd([dd.isdir]);
otsuDir = '';
for i = 1:numel(dd)
    if startsWith(dd(i).name, 'FotosOtsu', 'IgnoreCase', true)
        otsuDir = fullfile(dropletDir, dd(i).name); break
    end
end
assert(~isempty(otsuDir), 'No FotosOtsu folder found in %s', dropletDir);

% ---------- locate the raw image folder (FotosS<n>) ----------
% Must match FotosS<n> exactly, not FotosProcessed and not FotosOtsu.
photoDir = '';
for i = 1:numel(dd)
    if ~isempty(regexpi(dd(i).name, '^FotosS?\d+$', 'once'))
        photoDir = fullfile(dropletDir, dd(i).name); break
    end
end
assert(~isempty(photoDir), 'No FotosS<n> folder found in %s', dropletDir);
fprintf('Raw images: %s\n', photoDir);

% ---------- results folder, inside the droplet folder ----------
resultsFolder = fullfile(dropletDir, 'RESULTS_crystals');
if ~isfolder(resultsFolder), mkdir(resultsFolder); end

fprintf('Droplet: %s   fps=%g  ppf=%.4f  Ncrys=%d\n', ...
        specs.nombre, fps, um_per_px, NCrys);

%% ========================================================================
%  8. CANONICAL TIME AXIS  -  the specifications are always authoritative
%  ========================================================================
%  t_start is the REAL clock time at which frame 1 was taken. If t_start is
%  124 s, frame 1 is not second 0: it is second 124.
%      frameTimes(i) = t_start + (i-1)/fps
%  Consistency check with a sample specifications file:
%      t_start=124, fps=2, nucleation_frame=152
%      -> 124 + 151/2 = 199.5 s = nucleation_time  [OK]
%
%  Nothing stored in a .mat file may override this. Whatever is wrong on
%  disk gets REWRITTEN (syncTimeAxis), never the other way round.

% Number of frames, counted from the thresholded image folder.
ff = dir(fullfile(otsuDir, '*.tif'));
if isempty(ff), ff = dir(fullfile(otsuDir, '*.png')); end
nTotal = numel(ff);
if ~isnan(frame_max) && frame_max <= nTotal, nTotal = frame_max; end

if isfield(specs, 't_start') && ~isnan(specs.t_start)
    t0 = specs.t_start;
else
    t0 = 0;
end

frameTimes        = t0 + (0:nTotal-1) / fps;
t_start           = t0;
frameTimes_master = frameTimes;

% ---------- nucleation frames, also from the specifications ----------
if isfield(specs,'nucleation_frames') && ~isempty(specs.nucleation_frames) ...
                                      && all(~isnan(specs.nucleation_frames))
    nucFrames = sort(specs.nucleation_frames(:)');
elseif ~isnan(nucFrame)
    nucFrames = nucFrame;
    warning(['Only nucleation_frame is present. Add ' ...
             'nucleation_frames=f1,f2,... for the per-crystal analysis.']);
else
    nucFrames = [];
end
nucFrames(nucFrames < 1 | nucFrames > nTotal) = [];
assert(~isempty(nucFrames), 'No valid nucleation frames in the specifications.');

nucTimes     = frameTimes(nucFrames);   % nucleation time of EACH crystal
t_nucleation = nucTimes(1);             % first nucleation starts the global curve

if ~isnan(specs.nucleation_time) && abs(specs.nucleation_time - t_nucleation) > 1/fps
    warning(['nucleation_time=%.2f s in the specifications, but frame %d ' ...
             'falls at t=%.2f s. THE FRAME WINS.'], ...
             specs.nucleation_time, nucFrames(1), t_nucleation);
end
if ~isnan(NCrys) && numel(nucFrames) ~= NCrys
    warning(['Ncrys=%d but there are %d nucleation frames. Ncrys from the ' ...
             'specifications wins.'], NCrys, numel(nucFrames));
end

fprintf('\n--- TIME AXIS (from the specifications) ---\n');
fprintf('  %d frames | t = %.2f -> %.2f s\n', nTotal, frameTimes(1), frameTimes(end));
fprintf('  t_start      : %.2f s   (frame 1)\n', t_start);
fprintf('  nucleation   : frames %s\n', mat2str(nucFrames));
fprintf('                 t      %s s\n', mat2str(round(nucTimes,2)));
fprintf('  t_nucleation : %.2f s   (start of the global growth curve)\n', t_nucleation);
fprintf('  SMOOTH_WINDOW: %d frames (shared by crystals, PIV and TrackMate)\n\n', ...
        SMOOTH_WINDOW);

% ---------- migrate any stale time variables already on disk ----------
%  Runs unconditionally, independent of the stage flags. It only touches the
%  time variables; areas, velocities and masks are left intact.
syncTimeAxis(resultsFolder, dropletName, t_start, fps, nucFrames);

%% ========================================================================
%  9. RUN THE SELECTED STAGES
%  ========================================================================
%  frameTimes is reset before every stage because the stages are scripts and
%  some of them load .mat files that may still redefine it.

if DO_CRYSTAL_ANALYSIS
    frameTimes = frameTimes_master;
    run(fullfile('analysis', 'crystalGrowthAnalysis.m'));
end

if DO_PIV
    frameTimes = frameTimes_master;
    run(fullfile('analysis', 'pivFlowAnalysis.m'));
end

if DO_TRACKMATE
    frameTimes = frameTimes_master;
    run(fullfile('analysis', 'trackFlowAnalysis.m'));
end

if DO_COMPARISON
    frameTimes = frameTimes_master;
    run(fullfile('analysis', 'growthFlowComparison.m'));
end
