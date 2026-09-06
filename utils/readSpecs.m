function specs = readSpecs(sessionDir)
%READSPECS  Read the per-droplet acquisition specifications file.
%
%   Each droplet folder contains a plain-text file named 'especificaciones'
%   holding the acquisition metadata that cannot be recovered from the images
%   themselves: frame rate, spatial calibration, ambient conditions, the real
%   clock time of the first frame, and the frames at which each crystal
%   nucleated.
%
%   This file is the single authoritative source for all of those quantities.
%   Nothing stored in a .mat result file is allowed to override it; the time
%   axis of the whole pipeline is rebuilt from t_start and fps every run (see
%   SYNCTIMEAXIS).
%
%   SYNTAX
%     specs = readSpecs(sessionDir)
%
%   FILE FORMAT
%     One 'key=value' pair per line. Blank lines and lines without '=' are
%     ignored. Keys are case-insensitive and several spellings are accepted
%     for each field, since the files were written by hand across sessions.
%     The spatial scale may be given as a fraction ('717/900'), which is how
%     it comes out of the calibration slide measurement.
%
%   OUTPUT FIELDS
%     nombre             Droplet identifier as written in the file
%     temperatura        Ambient temperature (C)
%     humedad            Relative humidity (%)
%     fps                Acquisition frame rate (frames/s)
%     ppf                Spatial calibration (um/px)
%     frame_min          First usable frame
%     frame_max          Last usable frame
%     nucleation_frame   Frame of the FIRST nucleation
%     nucleation_frames  Frames of every nucleation, one per crystal
%     nucleation_time    Clock time of the first nucleation (s)
%     t_start            Clock time of frame 1 (s)
%     Ncrys              Number of crystals in the droplet
%     file               Full path of the file that was read
%
%   The function errors if fps is missing or invalid, and warns if ppf or
%   t_start are missing, because a missing t_start silently shifts every time
%   axis to zero.
%
%   See also SYNCTIMEAXIS, RUNDROPLETPIPELINE.

f = dir(fullfile(sessionDir, 'especificaciones*'));
if isempty(f)
    error('No specifications file found in %s', sessionDir);
end
specPath = fullfile(sessionDir, f(1).name);

fid = fopen(specPath, 'r', 'n', 'UTF-8');
raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);
lines = raw{1};

specs = struct('nombre','', 'temperatura',NaN, 'humedad',NaN, ...
               'fps',NaN, 'ppf',NaN, 'frame_min',NaN, 'frame_max',NaN, ...
               'nucleation_frame',NaN, 'nucleation_time',NaN, ...
               't_start',NaN, 'Ncrys',NaN, 'file',specPath);

for i = 1:numel(lines)
    L = strtrim(lines{i});
    if isempty(L) || ~contains(L, '='), continue; end
    parts = strsplit(L, '=');
    key = lower(strtrim(parts{1}));
    val = strtrim(strjoin(parts(2:end), '='));

    switch key
        case 'nombre'
            specs.nombre = val;
        case {'temperatura','temperature'}
            specs.temperatura = str2double(strrep(val, ',', '.'));
        case {'humedad','humidity','rh'}
            specs.humedad = str2double(strrep(strrep(val,'%',''), ',', '.'));
        case 'fps'
            specs.fps = str2double(val);
        case {'ppf','um_per_px','umperpx'}
            specs.ppf = parseNumOrFraction(val);
        case 'frame_min'
            specs.frame_min = str2double(val);
        case 'frame_max'
            specs.frame_max = str2double(val);
        case {'nucleation_frame','frame_nucleacion'}
            specs.nucleation_frame = str2double(val);
        case {'nucleation_time','t_nucleacion'}
            specs.nucleation_time = str2double(val);
        case {'t_start','tstart','t_inicial','t_inicio','tiempo_inicial', ...
              'tframe1','t_frame1','tiempo_start','t0'}
            specs.t_start = str2double(val);
        case {'ncrys','n_crys','ncristales'}
            specs.Ncrys = str2double(val);
        case {'nucleation_frames','frames_nucleacion'}
            specs.nucleation_frames = str2double(strsplit(val, ','));
    end
end

%% ---------- validation ----------
if isnan(specs.fps) || specs.fps <= 0
    error('Invalid fps in %s', specPath);
end
if isnan(specs.ppf) || specs.ppf <= 0
    warning('Invalid ppf in %s', specPath);
end
if isnan(specs.t_start)
    warning(['No t_start in %s. The time axis WILL START AT 0 s. ' ...
             'Add a line  t_start=<seconds>  to the file.'], specPath);
end

%% ---------- console summary ----------
fprintf('--- specifications ---\n');
fprintf('  name        : %s\n', specs.nombre);
fprintf('  temperature : %.1f C\n', specs.temperatura);
fprintf('  humidity    : %.0f %%\n', specs.humedad);
fprintf('  fps         : %g  (%.3f s/frame)\n', specs.fps, 1/specs.fps);
fprintf('  ppf         : %.4f um/px\n', specs.ppf);
if isnan(specs.t_start)
    fprintf('  t_start     : (MISSING -> 0 s will be used)\n');
else
    fprintf('  t_start     : %.2f s\n', specs.t_start);
end
if ~isnan(specs.Ncrys)
    fprintf('  Ncrys       : %d\n', specs.Ncrys);
end
if ~isnan(specs.frame_min)
    fprintf('  frames      : %d - %d\n', specs.frame_min, specs.frame_max);
end
if ~isnan(specs.nucleation_frame)
    fprintf('  nucleation  : frame %d (t = %.2f s)\n', ...
            specs.nucleation_frame, specs.nucleation_time);
end
end


% =========================================================================
function x = parseNumOrFraction(s)
%PARSENUMORFRACTION  Accept either a decimal number or an 'a/b' fraction.
%   The spatial calibration is measured as (known length)/(measured pixels),
%   so it is convenient to record it in the specifications file as a literal
%   fraction rather than a rounded decimal.

s = strtrim(strrep(s, ' ', ''));
t = regexp(s, '^([\d\.]+)/([\d\.]+)$', 'tokens', 'once');
if ~isempty(t)
    x = str2double(t{1}) / str2double(t{2});
else
    x = str2double(strrep(s, ',', '.'));
end
end
