function dimensionalAnalysis(dataRoot)
%DIMENSIONALANALYSIS  Scaling comparison of the three candidate flow mechanisms.
%
%   WHAT THIS SCRIPT DOES
%     Measuring an internal flow does not identify what drives it. This
%     function estimates, for every droplet in the data set, the velocity
%     that each candidate mechanism would produce under the conditions of
%     that particular experiment, and compares those predictions against the
%     velocities actually measured by PIV and particle tracking.
%
%     The three candidates are:
%       capillary replenishment  flow feeding evaporative loss at the contact line
%       solutal Marangoni        surface-tension gradient from the salinity field
%       natural convection       buoyancy from the same salinity field
%
%     The comparison is carried out in four steps:
%       STEP 1  droplet geometry inferred from the evaporation law
%       STEP 2  the salinity gradient set up by the growing crystal
%       STEP 3  the characteristic velocity of each mechanism
%       STEP 4  the depleted layer and the governing dimensionless groups
%
%   HOW TO USE IT
%     dimensionalAnalysis                  uses the default data root below
%     dimensionalAnalysis(dataRoot)        uses a different root folder
%
%   INPUTS
%     Per-droplet ambient conditions, spatial calibration and nucleation time
%     are read from each droplet's specifications file via READSPECS. Every
%     other measured quantity (apparent diameter, drying time, front velocity
%     and the measured flow velocities) is tabulated in BLOCK 1 below,
%     transcribed from the results tables of the thesis.
%
%   OUTPUTS
%     Written to <dataRoot>/DimensionalAnalysisMaths:
%       dimensional_per_droplet.csv   one row per droplet, every quantity
%       dimensional.mat               the same table plus the property struct
%       dimensional_log.txt           full console transcript
%
%   See also READSPECS.

%% ========================================================================
%  BLOCK 0 : PATHS
%  ========================================================================
if nargin < 1 || isempty(dataRoot)
    dataRoot = fullfile('C:', 'Data', 'FG');
end

sessionDirs = { ...
    fullfile(dataRoot,'TP2-B(2811)','S4' ), ...
    fullfile(dataRoot,'TP3(1202)'  ,'S5' ), ...
    fullfile(dataRoot,'TP3(1202)'  ,'S6' ), ...
    fullfile(dataRoot,'TP4-D(0106)','S7' ), ...
    fullfile(dataRoot,'TP4-D(0106)','S8' ), ...
    fullfile(dataRoot,'TP4-D(0106)','S9' ), ...
    fullfile(dataRoot,'TP4-D(0106)','S10'), ...
    fullfile(dataRoot,'TP5-E(0806)','S12'), ...
    fullfile(dataRoot,'TP5-E(0806)','S13'), ...
    fullfile(dataRoot,'TP5-E(0806)','S14'), ...
    fullfile(dataRoot,'TP5-E(0806)','S15'), ...
    fullfile(dataRoot,'TP5-E(0806)','S18'), ...
    fullfile(dataRoot,'TP5-E(0806)','S19'), ...
    fullfile(dataRoot,'TP5-E(0806)','S20'), ...
    fullfile(dataRoot,'TP6-(2907)' ,'S21'), ...
    fullfile(dataRoot,'TP6-(2907)' ,'S22'), ...
    fullfile(dataRoot,'TP6-(2907)' ,'S23'), ...
    fullfile(dataRoot,'TP6-(2907)' ,'S24'), ...
    fullfile(dataRoot,'TP6-(2907)' ,'S25'), ...
    fullfile(dataRoot,'TP6-(2907)' ,'S29'), ...
    fullfile(dataRoot,'TP6-(2907)' ,'S30'), ...
    fullfile(dataRoot,'TP6-(2907)' ,'S31'), ...
    fullfile(dataRoot,'TP6-(2907)' ,'S32'), ...
    fullfile(dataRoot,'TP6-(2907)' ,'S33'), ...
    fullfile(dataRoot,'TP6-(2907)' ,'S34'), ...
    fullfile(dataRoot,'TP6-(2907)' ,'S35')  ...
};

outDir = fullfile(dataRoot, 'DimensionalAnalysisMaths');
if ~exist(outDir, 'dir'), mkdir(outDir); end

diaryFile = fullfile(outDir, 'dimensional_log.txt');
if exist(diaryFile, 'file'), delete(diaryFile); end
diary(diaryFile); diary on;

fprintf('\n');
fprintf('########################################################################\n');
fprintf('#  DIMENSIONAL ANALYSIS                                                #\n');
fprintf('#  Method and results: scaling comparison of the three mechanisms      #\n');
fprintf('########################################################################\n');

%% ========================================================================
%  BLOCK 1 : MEASURED DATA, PER DROPLET
%  ========================================================================
%  Everything that is NOT in the specifications files. Row order matches
%  sessionDirs above; the check at the end of this block enforces that.

% Identifiers of the 26 retained droplets
id = [4 5 6 7 8 9 10 12 13 14 15 18 19 20 21 22 23 24 25 29 30 31 32 33 34 35]';

% Mean apparent diameter, in pixels: average of the largest and smallest
% apparent diameters measured on the last frames. Converted to a contact
% radius with that droplet's own calibration:  R_c [um] = d_px * ppf / 2
d_px = [2110 2132 2030 1550 1483 1395 1428 1445 2159 1998 2070 2188 ...
        1690 1682 1485 1965 1561 1686 1881 2072 2064 1843 1932 1834 ...
        2009 1809]';

% Initial NaCl concentration [g/L]
conc = [150 280 280 280 280 280 280 280 280 280 280 280 ...
        280 280 200 150 280 280 280 200 200 280 150 280 ...
        280 280]';

% Total drying time [s]
t_f = [379.5 424.5 323.5 210.8 227.3 187.3 182.3 221.3 257.0 315.0 ...
       153.7 263.0 151.7 133.3 252.5 316.5 266.5 329.5 270.5 273.5 ...
       416.0 308.5 325.5 342.5 435.5 283.5]';

% Mean crystal front velocity [um/s]
v_front_um = [0.9 1.6 2.2 2.2 1.9 2.3 1.8 1.2 2.1 1.2 1.1 1.3 2.1 1.9 ...
              1.5 1.9 1.7 1.5 1.3 1.7 1.4 0.9 1.7 1.0 0.9 1.6]';

% Measured flow velocities [um/s], from PIV and from particle tracking.
%
%   NaN in the pre-nucleation columns of droplets 4, 12, 14, 18, 20 and 35:
%   these already had crystals when the recording started.
%
%   NaN also for droplet 34: it nucleates at 66 s, the lowest value in the
%   set and only 15 % of its lifetime against a mean of 52 %. Its
%   pre-nucleation interval is too short, and the flow there is already
%   conditioned by the imminent crystallisation.
%
%   That leaves 19 droplets for the pre-nucleation comparison and 26 for
%   everything else.

piv_pre  = [NaN  0.49 0.64 4.88 9.65 6.70 4.43 NaN  0.64 NaN  0.39 NaN  ...
            0.42 NaN  0.69 1.21 1.23 0.09 0.19 0.74 0.14 0.52 0.20 0.19 ...
            NaN  NaN]';

piv_post = [0.75 3.13 3.73 8.26 11.70 8.93 8.54 1.28 3.18 1.07 0.27 1.12 ...
            3.80 3.77 0.72 0.85 3.23 6.62 2.24 0.71 0.41 3.72 0.10 1.10 ...
            0.47 0.36]';

piv_max  = [12.3 15.7 24.6 44.1 48.3 43.0 34.5 54.9 30.5 15.5 5.8 22.1 ...
            25.4 22.7 5.1 6.8 28.4 32.0 29.0 11.0 11.1 62.9 5.3 48.5 ...
            17.2 12.7]';

tm_pre   = [NaN  1.34 1.46 7.36 8.57 8.10 9.87 NaN  3.31 NaN  3.09 NaN  ...
            2.58 NaN  2.43 2.83 6.47 9.26 4.27 2.76 1.54 4.16 0.68 2.84 ...
            NaN  NaN]';

tm_post  = [2.58 4.19 3.44 11.24 12.07 11.77 13.67 3.51 6.92 3.34 1.83 ...
            3.31 6.67 7.56 2.24 2.16 8.62 9.25 5.68 1.81 1.70 3.87 0.83 ...
            2.85 2.83 2.02]';

tm_max   = [15.0 19.8 23.5 53.7 52.2 55.0 47.9 16.0 37.0 19.5 9.3 17.9 ...
            30.0 28.6 12.7 16.0 28.6 31.1 32.1 11.8 12.5 32.0 6.4 22.4 ...
            29.4 19.5]';

n = numel(id);

% ---------- verify that the paths are in the right order ----------
fprintf('\n=== PATH CHECK ===\n');
ok = numel(sessionDirs) == n;
if ~ok
    error('There are %d paths and %d droplets.', numel(sessionDirs), n);
end
for i = 1:n
    [~, folder] = fileparts(sessionDirs{i});
    if ~strcmpi(folder, sprintf('S%d', id(i)))
        fprintf('  ERROR row %2d: folder %s, droplet S%d\n', i, folder, id(i));
        ok = false;
    end
    if ~exist(sessionDirs{i}, 'dir')
        fprintf('  ERROR row %2d: %s does not exist\n', i, sessionDirs{i});
        ok = false;
    end
end
if ~ok, error('Fix the paths before continuing.'); end
fprintf('  %d paths correct and in order.\n', n);

%% ========================================================================
%  BLOCK 2 : PHYSICAL PROPERTIES
%  ========================================================================
%  All in SI units. Properties depending on temperature or concentration are
%  defined as anonymous functions and evaluated droplet by droplet.

P.g      = 9.81;              % m/s2
P.kB     = 1.380649e-23;      % J/K
P.Mw     = 0.018015;          % kg/mol, molar mass of water
P.R      = 8.314;             % J/(mol K)
P.M_NaCl = 58.44;             % g/mol, molar mass of NaCl

% ---------- saturation and supersaturation ----------
P.w_sat = 0.265;              % Pinho & Macedo 2005, at 25 C
P.S     = 1.6;                % Desarnaud 2014; Shahidzadeh 2015
P.dw    = (P.S-1)*P.w_sat;    % jump imposed by the crystal = 0.6*w_sat

% ---------- density: Novotny & Sohnel, via Madani Sani & Nesic 2024 ----------
%  c in mol/L, t in C, rho in kg/m3.
P.rho_w = @(t) 999.65 + 0.20438*t - 0.06174*t.^1.5;
P.rho   = @(c,t) P.rho_w(t) ...
                 + (44.85 - 0.09634*t + 0.0006136*t.^2).*c ...
                 + (-2.712 + 0.01009*t).*c.^1.5;

% ---------- density at saturation, solved by iteration ----------
%  At saturation the mass fraction is known but the molarity is not, so the
%  pair (c, rho) has to be solved self-consistently:
%      c = w*rho/M   ->   rho = rho(c)   ->   repeat until it settles.
r = 1190;
for k = 1:50
    c_it = P.w_sat*r/P.M_NaCl;
    r    = P.rho(c_it, 25);
end
P.rho_sat = r;
P.c_sat   = c_it;
P.rho_s   = (P.rho_sat - P.rho_w(25))/P.w_sat;   % solutal expansion coefficient

% ---------- viscosity: Mao & Duan 2009, at saturation and 25 C ----------
P.mu = 1.80e-3;               % Pa s

% ---------- surface tension ----------
%  Pure water (CRC, 71.99 mN/m at 25 C) plus the slope of Pegram & Record
%  (1.73 mN/m per mol/kg) extrapolated to saturation (6.17 mol/kg).
P.gamma   = 82.7e-3;          % N/m
P.gamma_s = 0.040;            % N/m per unit mass fraction

% ---------- mutual diffusivity of the salt: Annunziata et al. 2000 ----------
P.D_s = 1.50e-9;              % m2/s

% ---------- vapour diffusivity: 25 C reference, Fuller scaling ----------
P.D_vap = @(t) 2.42e-5*((t+273.15)/298.15).^1.75;

% ---------- saturated vapour concentration: Magnus + ideal gas ----------
P.p_sat = @(t) 610.94*exp(17.625*t./(t+243.04));      % Pa,  t in C
P.c_s   = @(t) P.p_sat(t)*P.Mw./(P.R*(t+273.15));     % kg/m3

% ---------- water activity of saturated NaCl: Greenspan 1977 ----------
P.chi_w = 0.753;

% ---------- solid NaCl and tracer particle ----------
P.rho_c = 2165;               % kg/m3, solid NaCl
P.a     = 0.4595e-6;          % m, tracer particle radius

fprintf('\n=== PHYSICAL PROPERTIES ===\n');
fprintf('  rho water     (25 C) : %8.2f kg/m3\n', P.rho_w(25));
fprintf('  rho saturation(25 C) : %8.1f kg/m3   (c = %.3f mol/L)\n', P.rho_sat, P.c_sat);
fprintf('  rho_s  derived       : %8.0f kg/m3\n', P.rho_s);
fprintf('  mu     at saturation : %8.2e Pa s\n', P.mu);
fprintf('  gamma  at saturation : %8.4f N/m\n', P.gamma);
fprintf('  gamma_s              : %8.3f N/m\n', P.gamma_s);
fprintf('  D_s (mutual)         : %8.2e m2/s\n', P.D_s);
fprintf('  w_sat                : %8.3f\n', P.w_sat);
fprintf('  Delta w_s = 0.6*w_sat: %8.4f\n', P.dw);
fprintf('  D_vap  (25 C)        : %8.3e m2/s\n', P.D_vap(25));
fprintf('  c_s    (25 C)        : %8.5f kg/m3\n', P.c_s(25));
fprintf('  chi_w                : %8.3f\n', P.chi_w);

%% ========================================================================
%  BLOCK 3 : READ THE SPECIFICATIONS
%  ========================================================================
T = nan(n,1); Hr = nan(n,1); ppf = nan(n,1); t_nuc = nan(n,1);

fprintf('\n=== READING SPECIFICATIONS ===\n');
for i = 1:n
    s = readSpecs(sessionDirs{i});
    T(i)     = s.temperatura;
    Hr(i)    = s.humedad/100;
    ppf(i)   = s.ppf;
    t_nuc(i) = s.nucleation_time;
end

% Contact radius: apparent diameter in px times that droplet's scale, halved
Rc      = d_px.*ppf/2*1e-6;                   % m
v_front = v_front_um*1e-6;                    % m/s
c_molar = conc/P.M_NaCl;                      % mol/L

%% ========================================================================
%  BLOCK 4 : STEP 1, DROPLET GEOMETRY
%  ========================================================================
%  The droplets are too thin for their height to be measured directly from a
%  top view, so the contact angle is inferred instead from how long they take
%  to dry: the evaporation law relates total drying time to the initial
%  volume, and the contact radius is measured.
%
%  Density is evaluated at each droplet's INITIAL concentration, because the
%  evaporation law describes its whole lifetime rather than a single instant.
rho_i   = P.rho(c_molar, T);
D_vap_i = P.D_vap(T);
c_s_i   = P.c_s(T);

% Driving force for evaporation. chi_w accounts for the saline solution
% having a lower water activity than pure water.
dc      = c_s_i.*(P.chi_w - Hr);
dc_pure = c_s_i.*(1       - Hr);              % upper bound, chi_w = 1

% Contact angle from the drying time
theta0   = 16*D_vap_i.*dc     .*t_f./(pi*Rc.^2.*rho_i);
theta0_p = 16*D_vap_i.*dc_pure.*t_f./(pi*Rc.^2.*rho_i);

% Heights and volume. h_nuc is the height still remaining when the crystal
% nucleates, which is the relevant length scale for the mechanisms that act
% through the depth of the droplet.
h0      = Rc.*theta0/2;
h_nuc   = h0.*(1 - t_nuc./t_f);
V       = pi*Rc.^3.*theta0/4;

h0_p    = Rc.*theta0_p/2;
h_nuc_p = h0_p.*(1 - t_nuc./t_f);
V_p     = pi*Rc.^3.*theta0_p/4;

% Bond number and capillary length: both confirm the droplets are flat enough
% for gravity to be negligible in setting their shape.
Bo    = rho_i*P.g.*Rc.^2/P.gamma;
l_cap = sqrt(P.gamma./(rho_i*P.g));

%% ========================================================================
%  BLOCK 5 : STEP 2, SALINITY GRADIENT
%  ========================================================================
%  The growing crystal holds its own surface at saturation while the rest of
%  the droplet stays supersaturated, so the concentration jump falls across a
%  distance of order the droplet radius.
dwdr = P.dw./Rc;                              % 1/m

%% ========================================================================
%  BLOCK 6 : STEP 3, CHARACTERISTIC VELOCITIES
%  ========================================================================
%  Capillary replenishment is evaluated at the INITIAL height, because it
%  operates from the moment evaporation begins. The other two are evaluated
%  at the height at nucleation, since they are driven by the gradient the
%  crystal creates. Viscosity, gamma_s and rho_s are taken at saturation,
%  which is the state of the droplet when the crystal appears.
vCa = D_vap_i.*dc./(rho_i.*h0)             *1e6;   % um/s, capillary
vMa = P.gamma_s*h_nuc/(4*P.mu).*dwdr       *1e6;   % um/s, solutal Marangoni
vRa = P.rho_s*P.g*h_nuc.^3/(48*P.mu).*dwdr *1e6;   % um/s, natural convection

% The Marangoni-to-buoyancy ratio is independent of the gradient: both scale
% linearly with it, so it cancels. What remains depends only on the height,
% which is why the comparison is decided by how thin the droplet is.
ratio_MaRa = 12*P.gamma_s./(P.rho_s*P.g*h_nuc.^2);
h_equal    = sqrt(12*P.gamma_s/(P.rho_s*P.g));

%% ========================================================================
%  BLOCK 7 : STEP 4, DEPLETED LAYER AND DIMENSIONLESS GROUPS
%  ========================================================================
%  Post-nucleation quantities are evaluated with the saturation density,
%  because by that moment the droplet is at or above saturation regardless of
%  the concentration it started from.

% Depleted layer: the salt-poor region the advancing crystal leaves behind
delta = P.rho_sat*P.D_s*P.dw./(P.rho_c*v_front);     % m

% Tracer diffusivity, Stokes-Einstein
D_p = P.kB*(T+273.15)/(6*pi*P.mu*P.a);               % m2/s

U_piv  = piv_post*1e-6;
U_tm   = tm_post *1e-6;
L_band = 20e-6;                                      % observed band width

% Reynolds: confirms the flow is deeply viscous
Re_piv = P.rho_sat*U_piv.*h_nuc/P.mu;
Re_tm  = P.rho_sat*U_tm .*h_nuc/P.mu;

% Particle Peclet: advection against Brownian motion for the tracers
Pep_piv = U_piv*L_band./D_p;
Pep_tm  = U_tm *L_band./D_p;

% Solutal Peclet: advection against salt diffusion across the depleted layer
Pes_piv = U_piv.*delta/P.D_s;
Pes_tm  = U_tm .*delta/P.D_s;

% Vertical mixing criterion: whether the flow has time to redistribute salt
% over the droplet depth before advecting it away
eps_  = h_nuc./Rc;
epsPe = eps_.*U_tm.*h_nuc/P.D_s;

% Diffusion lengths over a representative observation window
t_ev = 100;
Ld_p = sqrt(2*mean(D_p)*t_ev);
Ld_s = sqrt(2*P.D_s   *t_ev);

%% ========================================================================
%  BLOCK 8 : PER-DROPLET OUTPUT
%  ========================================================================
fprintf('\n=== PER-DROPLET DETAIL ===\n');
for i = 1:n
    fprintf('\n--- Droplet %d  (%g g/L, %.1f C, %.0f%% RH) ---\n', ...
            id(i), conc(i), T(i), Hr(i)*100);
    fprintf('  INPUT      d = %d px, ppf = %.4f um/px, t_f = %.1f s, t_nuc = %.1f s\n', ...
            d_px(i), ppf(i), t_f(i), t_nuc(i));
    fprintf('  PROPERTY   rho = %.1f kg/m3 (c = %.2f mol/L)   c_s = %.5f   D_vap = %.3e\n', ...
            rho_i(i), c_molar(i), c_s_i(i), D_vap_i(i));
    fprintf('             Delta c = %.5f kg/m3\n', dc(i));
    fprintf('  GEOMETRY   R_c = %.0f um   theta_0 = %.1f deg   h_0 = %.0f um   h_nuc = %.0f um\n', ...
            Rc(i)*1e6, rad2deg(theta0(i)), h0(i)*1e6, h_nuc(i)*1e6);
    fprintf('             V = %.3f uL   Bo = %.3f   l_cap = %.2f mm   eps = %.3f   eps*Pe = %.3f\n', ...
            V(i)*1e9, Bo(i), l_cap(i)*1e3, eps_(i), epsPe(i));
    fprintf('  GRADIENT   dw/dr = %.1f 1/m\n', dwdr(i));
    fprintf('  VELOCITY   v_Ca = %.2f   v_Ma = %.3e   v_Ra = %.2f um/s   v_Ma/v_Ra = %.2e\n', ...
            vCa(i), vMa(i), vRa(i), ratio_MaRa(i));
    fprintf('  MEASURED   PIV  pre %5.2f  post %5.2f  max %5.1f um/s\n', ...
            piv_pre(i), piv_post(i), piv_max(i));
    fprintf('             TM   pre %5.2f  post %5.2f  max %5.1f um/s\n', ...
            tm_pre(i), tm_post(i), tm_max(i));
    fprintf('  LAYER      v_front = %.2f um/s   delta = %.0f um\n', ...
            v_front_um(i), delta(i)*1e6);
    fprintf('  DIMLESS    D_p = %.2e m2/s\n', D_p(i));
    fprintf('             Re   PIV %.2e   TM %.2e\n', Re_piv(i), Re_tm(i));
    fprintf('             Pe_p PIV %7.1f    TM %7.1f\n', Pep_piv(i), Pep_tm(i));
    fprintf('             Pe_s PIV %7.3f    TM %7.3f\n', Pes_piv(i), Pes_tm(i));
end

%% ========================================================================
%  BLOCK 9 : GLOBAL SUMMARY
%  ========================================================================
prt = @(lbl,x,u) fprintf('  %-24s %9.3f +/- %8.3f   [%8.3f , %8.3f] %s\n', ...
        lbl, mean(x,'omitnan'), std(x,'omitnan'), min(x), max(x), u);
% Scientific-notation variant, for quantities %f cannot display sensibly
% (Reynolds numbers here are of order 1e-4).
prte = @(lbl,x,u) fprintf('  %-24s %9.3e +/- %8.3e   [%8.3e , %8.3e] %s\n', ...
        lbl, mean(x,'omitnan'), std(x,'omitnan'), min(x), max(x), u);

fprintf('\n\n========================================================================\n');
fprintf('  GLOBAL SUMMARY  (mean +/- standard deviation  [minimum , maximum])\n');
fprintf('========================================================================\n');

fprintf('\n--- STEP 1: GEOMETRY ---\n');
prt('R_c',                 Rc*1e6,            'um');
prt('theta_0',             rad2deg(theta0),   'deg');
prt('h_0',                 h0*1e6,            'um');
prt('h at nucleation',     h_nuc*1e6,         'um');
prt('Implied volume',      V*1e9,             'uL');
prt('Bond',                Bo,                '-');
prt('Capillary length',    l_cap*1e3,         'mm');
prt('epsilon = h/R_c',     eps_,              '-');
fprintf('  epsilon*Pe : median %.4f, maximum %.4f\n', median(epsPe), max(epsPe));
fprintf('  Nucleation at %.3f +/- %.3f of the droplet lifetime\n', ...
        mean(t_nuc./t_f), std(t_nuc./t_f));

fprintf('\n  With chi_w = 1 (upper bound):\n');
prt('theta_0',             rad2deg(theta0_p), 'deg');
prt('h at nucleation',     h_nuc_p*1e6,       'um');
prt('Implied volume',      V_p*1e9,           'uL');

fprintf('\n--- STEP 2: GRADIENT ---\n');
prt('dw_s/dr',             dwdr,              '1/m');
fprintf('  Using delta instead of R_c would give %.0f 1/m, a factor %.1f larger\n', ...
        mean(P.dw./delta), mean(P.dw./delta)/mean(dwdr));

fprintf('\n--- STEP 3: CHARACTERISTIC VELOCITIES ---\n');
prt('v_Ca',                vCa,               'um/s');
prt('v_Ma',                vMa,               'um/s');
prt('v_Ra',                vRa,               'um/s');
prt('v_Ma/v_Ra',           ratio_MaRa,        '-');
fprintf('  They would be equal at h = %.1f mm (mean capillary length %.1f mm)\n', ...
        h_equal*1e3, mean(l_cap)*1e3);

fprintf('\n  Measured velocities:\n');
prt('PIV pre-nucleation',  piv_pre,  'um/s');
prt('PIV post-nucleation', piv_post, 'um/s');
prt('PIV spatial maximum', piv_max,  'um/s');
prt('TM  pre-nucleation',  tm_pre,   'um/s');
prt('TM  post-nucleation', tm_post,  'um/s');
prt('TM  spatial maximum', tm_max,   'um/s');
fprintf('  n pre-nucleation = %d, n post-nucleation = %d\n', ...
        sum(~isnan(tm_pre)), n);

% ---------- ratios: MEDIAN of the per-droplet ratios ----------
%  This is not the ratio of the means. The distributions are strongly skewed,
%  so the median is the appropriate statistic.
med = @(x) median(x,'omitnan');
fprintf('\n  Prediction/measurement ratios (median of the per-droplet ratios):\n');
fprintf('    %-12s %12s %12s %12s %12s\n','mechanism','vs PIV','vs TM','vs PIV max','vs TM max');
fprintf('    %-12s %12.2f %12.2f %12s %12s\n','Capillary', ...
        med(vCa./piv_pre), med(vCa./tm_pre), '-', '-');
fprintf('    %-12s %12.2e %12.2e %12.2e %12.2e\n','Marangoni', ...
        med(vMa./piv_post), med(vMa./tm_post), ...
        med(vMa./piv_max),  med(vMa./tm_max));
fprintf('    %-12s %12.2f %12.2f %12.2f %12.2f\n','Nat. conv.', ...
        med(vRa./piv_post), med(vRa./tm_post), ...
        med(vRa./piv_max),  med(vRa./tm_max));

in10 = @(p,m) sum(p./m > 0.1 & p./m < 10);
fprintf('\n  Within a factor of 10 of the measurement:\n');
fprintf('    Capillary  %2d/%d (PIV pre) , %2d/%d (TM pre)\n', ...
        in10(vCa,piv_pre), sum(~isnan(piv_pre)), ...
        in10(vCa,tm_pre),  sum(~isnan(tm_pre)));
fprintf('    Nat. conv. %2d/%d (PIV post), %2d/%d (TM post)\n', ...
        in10(vRa,piv_post), n, in10(vRa,tm_post), n);
fprintf('    Nat. conv. %2d/%d (PIV max) , %2d/%d (TM max)\n', ...
        in10(vRa,piv_max), n, in10(vRa,tm_max), n);

[rP,pP] = pearson(vRa, tm_post);
[rM,pM] = pearson(vRa, tm_max);
fprintf('\n  Between-droplet correlation:\n');
fprintf('    v_Ra vs TM post : r = %+.2f, p = %.3f\n', rP, pP);
fprintf('    v_Ra vs TM max  : r = %+.2f, p = %.3f\n', rM, pM);

fprintf('\n--- STEP 4: DEPLETED LAYER AND DIMENSIONLESS GROUPS ---\n');
prt('delta',                delta*1e6,  'um');
prt('D_p (Stokes-Einstein)',D_p*1e13,   'x1e-13 m2/s');
prte('Re  (PIV)',           Re_piv,     '-');
prte('Re  (TM)',            Re_tm,      '-');
prt('Pe_p (PIV)',           Pep_piv,    '-');
prt('Pe_p (TM)',            Pep_tm,     '-');
prt('Pe_s (PIV)',           Pes_piv,    '-');
prt('Pe_s (TM)',            Pes_tm,     '-');
fprintf('\n  Diffusion lengths in %d s:\n', t_ev);
fprintf('    particle : %5.1f um   (band width %.0f um)\n', Ld_p*1e6, L_band*1e6);
fprintf('    salt     : %5.0f um   (mean R_c %.0f um)\n', Ld_s*1e6, mean(Rc)*1e6);
fprintf('  Time for a particle to diffuse one band width: %.0f s\n', ...
        L_band^2/(2*mean(D_p)));

%% ========================================================================
%  BLOCK 10 : SAVE
%  ========================================================================
out = table(id, conc, T, Hr*100, ppf, d_px, Rc*1e6, rho_i, ...
            t_f, t_nuc, rad2deg(theta0), h0*1e6, h_nuc*1e6, V*1e9, ...
            Bo, l_cap*1e3, eps_, epsPe, dwdr, ...
            vCa, vMa, vRa, ratio_MaRa, ...
            piv_pre, piv_post, piv_max, tm_pre, tm_post, tm_max, ...
            v_front_um, delta*1e6, D_p, ...
            Re_piv, Re_tm, Pep_piv, Pep_tm, Pes_piv, Pes_tm, ...
    'VariableNames', {'id','conc_gL','T_C','RH_pct','ppf','d_px','Rc_um', ...
            'rho','t_f_s','t_nuc_s','theta0_deg','h0_um','hnuc_um','V_uL', ...
            'Bo','lcap_mm','eps','epsPe','dwdr', ...
            'vCa','vMa','vRa','vMa_vRa', ...
            'piv_pre','piv_post','piv_max','tm_pre','tm_post','tm_max', ...
            'v_front','delta_um','D_p', ...
            'Re_piv','Re_tm','Pep_piv','Pep_tm','Pes_piv','Pes_tm'});

writetable(out, fullfile(outDir, 'dimensional_per_droplet.csv'));
save(fullfile(outDir, 'dimensional.mat'), 'out', 'P');

fprintf('\n\nSaved in %s:\n', outDir);
fprintf('  dimensional_per_droplet.csv\n');
fprintf('  dimensional.mat\n');
fprintf('  dimensional_log.txt\n');

diary off;
end


% =========================================================================
function [r, p] = pearson(x, y)
%PEARSON  Pearson correlation coefficient and its p-value, toolbox-free.
%
%   Pairs containing NaN are dropped. The p-value comes from the statistic
%   t = r*sqrt((m-2)/(1-r^2)) with m-2 degrees of freedom, evaluated
%   two-tailed with the incomplete beta function (betainc ships with base
%   MATLAB, so the Statistics Toolbox is not required to run this analysis).

ok = ~isnan(x) & ~isnan(y);
x  = x(ok);  y = y(ok);
m  = numel(x);

if m < 3
    r = NaN; p = NaN; return
end

x = x - mean(x);
y = y - mean(y);
r = sum(x.*y)/sqrt(sum(x.^2)*sum(y.^2));

if abs(r) >= 1
    p = 0;
else
    df = m - 2;
    t  = r*sqrt(df/(1 - r^2));
    p  = betainc(df/(df + t^2), df/2, 0.5);
end
end
