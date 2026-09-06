function syncTimeAxis(resultsFolder, dropletName, t_start, fps, nucFrames)
%SYNCTIMEAXIS  Rewrite stale time variables inside already-saved .mat files.
%
%   WHY THIS EXISTS
%     The specifications file is the single authoritative source for the time
%     axis. Results computed in earlier runs, however, carry whatever time
%     vector was in force when they were saved. If t_start is later corrected
%     in the specifications, every stored curve would silently keep plotting
%     against the old axis.
%
%     This function walks the .mat files of one droplet and, wherever a time
%     variable disagrees with the specifications, overwrites that variable in
%     place with save(...,'-append'). Nothing else is touched: areas, masks,
%     velocities and every other stored quantity are left exactly as they
%     were. A corrected t_start therefore repairs all existing results without
%     re-running any analysis.
%
%   SYNTAX
%     syncTimeAxis(resultsFolder, dropletName, t_start, fps, nucFrames)
%
%   TIME CONVENTION (identical everywhere in the pipeline)
%     frameTimes(i) = t_start + (i-1)/fps       frame 1 IS t_start
%     tSpeed(k)     = t_start + framesU(k)/fps  TrackMate counts frames from 0
%     nucList                                    always taken from the specs
%
%   VARIABLES MIGRATED
%     frameTimes   in crystal_areas, growth_global, growth_individual,
%                  radial_velocity, ...
%     t            in OneDPlots_*
%     tSpeed       in trackmate_meanspeed_*  (needs framesU in the same file)
%     nucList      in crystal_areas
%
%   Files containing none of those variables (masks_*, individual_crystals_*)
%   are never opened: only the header is inspected with whos('-file',...),
%   which matters because the mask files are large and stored as -v7.3.
%
%   See also READSPECS, RUNDROPLETPIPELINE.

if nargin < 5, nucFrames = []; end
TOL = 1e-9;

fprintf('\n=== TIME AXIS CHECK (specifications are authoritative) ===\n');
fprintf('  t_start = %.4f s | fps = %g  ->  frameTimes(i) = %.4f + (i-1)/%g\n', ...
        t_start, fps, t_start, fps);
if ~isempty(nucFrames)
    fprintf('  nucleation: frames %s\n', mat2str(nucFrames(:)'));
end

d = dir(fullfile(resultsFolder, ['*' dropletName '.mat']));
if isempty(d)
    fprintf('  No .mat files for %s in %s\n\n', dropletName, resultsFolder);
    return
end

timeVars = {'frameTimes','t','tSpeed'};
nChanged = 0;  nOk = 0;  nSkip = 0;

for k = 1:numel(d)
    f = fullfile(d(k).folder, d(k).name);

    % Read the header only: this does not load the data, which matters for
    % the large -v7.3 mask files.
    try
        info = whos('-file', f);
    catch ME
        fprintf('  [!] %-46s unreadable (%s)\n', d(k).name, ME.identifier);
        continue
    end
    present = {info.name};

    hasTime = any(ismember(timeVars, present));
    hasNuc  = ismember('nucList', present) && ~isempty(nucFrames);
    if ~hasTime && ~hasNuc
        nSkip = nSkip + 1;
        continue
    end

    changed = {};

    % ---------- frameTimes and t: same convention ----------
    for v = {'frameTimes','t'}
        vn = v{1};
        if ~ismember(vn, present), continue; end
        S = load(f, vn);
        old = S.(vn);
        if ~isnumeric(old) || isempty(old), continue; end
        N   = numel(old);
        new = t_start + (0:N-1) / fps;
        if size(old,1) > 1 && size(old,2) == 1   % column vector: keep shape
            new = new(:);
        else
            new = reshape(new, size(old));
        end
        if numel(old) == numel(new) && max(abs(old(:) - new(:))) <= TOL
            continue
        end
        eval([vn ' = new;']);                    %#ok<EVLEQ>  dynamic name
        save(f, vn, '-append');
        changed{end+1} = sprintf('%s: %.2f-%.2f -> %.2f-%.2f s (%d frames)', ...
                          vn, old(1), old(end), new(1), new(end), N);  %#ok<AGROW>
        clear(vn)
    end

    % ---------- tSpeed: indexed by framesU, which starts at 0 ----------
    if ismember('tSpeed', present)
        if ~ismember('framesU', present)
            fprintf('  [!] %-46s has tSpeed but no framesU: cannot migrate\n', ...
                    d(k).name);
        else
            S = load(f, 'tSpeed', 'framesU');
            old = S.tSpeed;  fu = double(S.framesU(:));
            new = t_start + fu / fps;
            new = reshape(new, size(old));
            if ~(numel(old) == numel(new) && max(abs(old(:) - new(:))) <= TOL)
                tSpeed = new;                                     %#ok<NASGU>
                save(f, 'tSpeed', '-append');
                changed{end+1} = sprintf(['tSpeed: %.2f-%.2f -> %.2f-%.2f s ' ...
                                  '(%d points)'], old(1), old(end), ...
                                  new(1), new(end), numel(new));   %#ok<AGROW>
                clear tSpeed
            end
        end
    end

    % ---------- nucList: always the one from the specifications ----------
    if hasNuc
        S = load(f, 'nucList');
        old = double(S.nucList(:)');
        new = sort(double(nucFrames(:)'));
        if numel(old) ~= numel(new) || any(old ~= new)
            nucList = new;                                        %#ok<NASGU>
            save(f, 'nucList', '-append');
            changed{end+1} = sprintf('nucList: %s -> %s', ...
                              mat2str(old), mat2str(new));         %#ok<AGROW>
            clear nucList
        end
    end

    if isempty(changed)
        nOk = nOk + 1;
        fprintf('  [ok]  %-46s already correct\n', d(k).name);
    else
        nChanged = nChanged + 1;
        fprintf('  [FIX] %s\n', d(k).name);
        for c = 1:numel(changed)
            fprintf('          %s\n', changed{c});
        end
    end
end

fprintf('  --- %d corrected | %d already correct | %d without time data ---\n\n', ...
        nChanged, nOk, nSkip);
end
