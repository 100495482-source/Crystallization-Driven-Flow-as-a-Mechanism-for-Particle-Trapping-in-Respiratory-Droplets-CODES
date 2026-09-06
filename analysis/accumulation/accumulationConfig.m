function cfg = accumulationConfig(dataRoot)
%ACCUMULATIONCONFIG  Shared settings for the particle accumulation pipeline.
%
%   The seven accumulation stages all operate on the same set of droplets and
%   write to the same folders. Declaring that list once here means a droplet
%   added or removed cannot end up applied to some stages and not others,
%   which is what happened while the list was copied into every script.
%
%   SYNTAX
%     cfg = accumulationConfig
%     cfg = accumulationConfig(dataRoot)
%
%   OUTPUT FIELDS
%     droplets  Cell array of droplet folder paths
%     dataRoot  Root folder holding the acquisition sessions
%     workDir   Intermediate results (tmclean*, buffers*)
%     outDir    Final figures and tables
%
%   WHICH DROPLETS ARE HERE
%     A subset of the full data set: only the droplets whose crystals were
%     individually segmented, since every stage below needs
%     individual_crystals_*.mat to exist.
%
%   See also ACCUM1_PREPARETRACKS.

if nargin < 1 || isempty(dataRoot)
    dataRoot = fullfile('C:', 'Data', 'FG');
end

sessions = { ...
    'TP2-B(2811)', 'S4' ; ...
    'TP3(1202)',   'S5' ; ...
    'TP3(1202)',   'S6' ; ...
    'TP5-E(0806)', 'S12'; ...
    'TP5-E(0806)', 'S13'; ...
    'TP5-E(0806)', 'S15'; ...
    'TP5-E(0806)', 'S18'; ...
    'TP5-E(0806)', 'S20'; ...
    'TP6-(2907)',  'S21'; ...
    'TP6-(2907)',  'S22'; ...
    'TP6-(2907)',  'S24'; ...
    'TP6-(2907)',  'S30'; ...
    'TP6-(2907)',  'S32'; ...
    'TP6-(2907)',  'S35'  ...
};

cfg.droplets = cell(size(sessions,1), 1);
for i = 1:size(sessions,1)
    cfg.droplets{i} = fullfile(dataRoot, sessions{i,1}, sessions{i,2});
end

cfg.dataRoot = dataRoot;
cfg.workDir  = fullfile(dataRoot, 'trackmateparticleaccumulation');
cfg.outDir   = fullfile(dataRoot, 'Accumulation');
end
