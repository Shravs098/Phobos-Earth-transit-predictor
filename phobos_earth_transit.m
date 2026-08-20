%% ========================================================================
%  PHOBOS-EARTH TRANSIT PREDICTOR
%  As seen from Perseverance Rover, Jezero Crater, Mars
% =========================================================================
%
%  QUESTION: On what dates/times would Phobos (Mars's inner moon) pass
%  directly in front of Earth, as seen from Perseverance's exact landing
%  site? Phobos crosses the Martian sky ~3x/day; Earth drifts slowly
%  through that sky over weeks. This code searches for the rare moments
%  the two align.
%
%  METHOD:
%   1. Propagate Earth & Mars heliocentric orbits (Keplerian elements)
%   2. Propagate Phobos's orbit around Mars (Keplerian, Mars-equatorial frame)
%   3. Apply the real IAU Mars pole/rotation model to place Jezero Crater
%      precisely on the rotating Martian surface at each instant
%   4. Compute the topocentric (rover's-eye-view) angular separation
%      between Phobos and Earth, including PARALLAX (Phobos is close
%      enough that the rover's exact location matters)
%   5. Search for transits: separation < Phobos's angular radius
%
%  CAVEAT: Phobos's orbital PHASE (mean anomaly at epoch) is not pulled
%  from a live ephemeris here (no SPICE access in this environment) - a
%  representative epoch phase is used. The orbital mechanics, Mars
%  rotation model, and search method are all physically real; for
%  mission-grade exact calendar dates, replace the Phobos epoch state
%  with true values from JPL Horizons / SPICE kernels.
%
%  Requires: kepler_solve.m, kepler_orbit_to_frame.m,
%            kepler_orbit_to_frame_vec.m, heliocentric_position.m,
%            evaluate_geometry.m, mkdir_if_missing.m  (same folder)
% =========================================================================

clear; close all; clc;

% Anchor the output folder to THIS script's location, not MATLAB's current
% folder - avoids "invalid path" errors if the script is run while a
% different folder is active.
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end   % fallback if run as a pasted block
addpath(script_dir);   % make sure the helper .m files are found regardless of current folder
outdir = fullfile(script_dir, 'figures');
mkdir_if_missing(outdir);

%% ---------------- 1. Search window & rover location ----------------
date_start = datetime(2026,1,1,0,0,0);
search_days = 90;                          % length of search window
dt_coarse_min = 5;                         % coarse scan resolution (minutes)

% Perseverance rover, Jezero Crater landing site (areocentric, East longitude)
rover_lat = 18.4447;    % degrees N
rover_lon = 77.4508;    % degrees E
R_mars    = 3389.5;     % km, mean Mars radius

fprintf('=== Phobos-Earth Transit Search ===\n');
fprintf('Location: Perseverance, Jezero Crater (%.4f N, %.4f E)\n', rover_lat, rover_lon);
fprintf('Window: %s + %d days, %d-min resolution\n\n', datestr(date_start), search_days, dt_coarse_min);

%% ---------------- 2. Time array (Julian date / days since J2000) ----------------
JD0 = datenum(date_start) + 1721058.5;                 % Julian date at window start
JD_J2000 = 2451545.0;
n_coarse = round(search_days*24*60/dt_coarse_min);
t_days   = (0:n_coarse)' * (dt_coarse_min/1440);        % days since window start
JD       = JD0 + t_days;
d_since_J2000 = JD - JD_J2000;                          % days since J2000 (TT approx)
T_cent   = d_since_J2000/36525;                         % Julian centuries since J2000

%% ---------------- 3. Heliocentric orbits: Earth & Mars ----------------
% Low-precision Keplerian elements (JPL "Keplerian Elements for Approximate
% Positions of the Major Planets", valid 1800-2050), degrees & AU, with
% secular rates per Julian century.
% [a, e, I, L, long_peri, long_node]   and rates [.../cy]
earth_el0   = [1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0];
earth_rate  = [0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0];
mars_el0    = [1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891];
mars_rate   = [0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343];

R_earth_helio = heliocentric_position(earth_el0, earth_rate, T_cent);   % AU, equatorial J2000
R_mars_helio  = heliocentric_position(mars_el0,  mars_rate,  T_cent);   % AU, equatorial J2000

AU_km = 149597870.7;

%% ---------------- 4. Mars pole orientation & rotation (IAU model) ----------------
% Approximate IAU Mars rotational elements (secular terms only)
alpha0 = (317.269 - 0.108*T_cent) * pi/180;     % pole RA
delta0 = (54.432  - 0.061*T_cent) * pi/180;     % pole Dec
W0_deg = 176.630;  Wdot_deg = 350.89198226;     % prime meridian angle & rate (deg, deg/day)
W = mod(W0_deg + Wdot_deg*d_since_J2000, 360) * pi/180;

% Mars pole unit vector (J2000 equatorial) and equatorial-node direction
p_hat0 = [cos(delta0(1))*cos(alpha0(1)), cos(delta0(1))*sin(alpha0(1)), sin(delta0(1))];
N_hat0 = [-sin(alpha0(1)), cos(alpha0(1)), 0];             % ascending node dir (treat as const)
PxN0   = cross(p_hat0, N_hat0);

% Mars-equatorial INERTIAL frame basis (non-rotating; used for Phobos's orbit)
Xeq = N_hat0;  Yeq = PxN0;  Zeq = p_hat0;

% Mars BODY-FIXED frame basis (rotates with spin; used for the rover)
Xbf = cos(W).*N_hat0 + sin(W).*PxN0;     % Nx3
Ybf = cross(repmat(p_hat0,length(W),1), Xbf, 2);
Zbf = repmat(p_hat0, length(W), 1);

%% ---------------- 5. Rover position (Jezero Crater) in inertial frame ----------------
lat_r = rover_lat*pi/180; lon_r = rover_lon*pi/180;
r_bf = [cos(lat_r)*cos(lon_r), cos(lat_r)*sin(lon_r), sin(lat_r)];   % body-fixed unit dir

Rover_pos = R_mars*( r_bf(1).*Xbf + r_bf(2).*Ybf + r_bf(3).*Zbf );   % Nx3, km, Mars-centered
zenith_hat = Rover_pos / R_mars;                                     % unit vector, local "up"

% local topocentric East/North unit vectors (for sky-map azimuth), rotating with Mars
East_bf  = [-sin(lon_r), cos(lon_r), 0];
North_bf = [-sin(lat_r)*cos(lon_r), -sin(lat_r)*sin(lon_r), cos(lat_r)];
East_hat  = East_bf(1).*Xbf  + East_bf(2).*Ybf  + East_bf(3).*Zbf;    % Nx3
North_hat = North_bf(1).*Xbf + North_bf(2).*Ybf + North_bf(3).*Zbf;  % Nx3

%% ---------------- 6. Phobos orbit around Mars ----------------
a_phobos = 9376;          % km
e_phobos = 0.0151;
i_phobos = 1.093*pi/180;  % inclination to Mars equator
Om_phobos = 0*pi/180;     % node (arbitrary reference epoch phase - see caveat)
w_phobos  = 0*pi/180;     % argument of periapsis (arbitrary reference epoch phase)
M0_phobos = 0*pi/180;     % mean anomaly at window start (arbitrary reference epoch phase)
period_phobos = 0.31891023;    % days
n_phobos = 2*pi/period_phobos; % mean motion, rad/day
R_phobos_mean = 11.1;          % km, mean radius (irregular body ~13x11x9 km)

M_phobos = mod(M0_phobos + n_phobos*t_days, 2*pi);
[xr,yr,zr] = kepler_orbit_to_frame(a_phobos, e_phobos, i_phobos, Om_phobos, w_phobos, M_phobos);
Phobos_pos = xr.*Xeq + yr.*Yeq + zr.*Zeq;      % Nx3, km, Mars-centered inertial

%% ---------------- 7. Rover-relative vectors, angular separation ----------------
d_phobos = Phobos_pos - Rover_pos;                                    % km
d_earth  = (R_earth_helio - R_mars_helio)*AU_km - Rover_pos;          % km

dist_phobos = vecnorm(d_phobos,2,2);
dist_earth  = vecnorm(d_earth,2,2);

cos_sep = sum(d_phobos.*d_earth,2) ./ (dist_phobos.*dist_earth);
sep_deg = acosd(min(1,max(-1,cos_sep)));                              % angular separation, deg

ang_radius_phobos_deg = asind(R_phobos_mean./dist_phobos);            % Phobos's apparent radius

elev_phobos = asind( sum(zenith_hat.*d_phobos,2)./dist_phobos );
elev_earth  = asind( sum(zenith_hat.*d_earth,2)./dist_earth );
az_phobos = mod(atan2d(sum(East_hat.*d_phobos,2), sum(North_hat.*d_phobos,2)), 360);
az_earth  = mod(atan2d(sum(East_hat.*d_earth,2),  sum(North_hat.*d_earth,2)),  360);
visible = (elev_phobos > 0) & (elev_earth > 0);

fprintf('Phobos angular radius from rover: %.3f - %.3f arcmin (varies with orbit distance)\n', ...
    min(ang_radius_phobos_deg)*60, max(ang_radius_phobos_deg)*60);

%% ---------------- 8. Coarse search for candidate close approaches ----------------
candidate_thresh_deg = 1.0;      % flag anything within 1 deg for refinement
below = sep_deg < candidate_thresh_deg & visible;
idx_candidates = find(below(2:end-1) & sep_deg(2:end-1)<sep_deg(1:end-2) & sep_deg(2:end-1)<sep_deg(3:end)) + 1;

fprintf('\nFound %d candidate close-approach window(s) within %.1f deg (visible above horizon).\n', ...
    length(idx_candidates), candidate_thresh_deg);

%% ---------------- 9. Refine each candidate with fine time resolution ----------------
results = struct('t_days',{},'JD',{},'sep_min_deg',{},'phobos_radius_deg',{}, ...
                  'is_transit',{},'t_fine',{},'sep_f',{},'ang_rad_f',{});
for c = 1:length(idx_candidates)
    t_center = t_days(idx_candidates(c));
    t_fine = (t_center - 1/24 : (5/86400) : t_center + 1/24)';   % +/- 1 hr, 5-sec steps

    [sep_f, ang_rad_f, elev_ph_f, elev_ea_f, d_ph_f, d_ea_f, az_ph_f, az_ea_f] = evaluate_geometry( ...
        t_fine, JD0, JD_J2000, earth_el0, earth_rate, mars_el0, mars_rate, ...
        AU_km, a_phobos, e_phobos, i_phobos, Om_phobos, w_phobos, M0_phobos, n_phobos, ...
        R_phobos_mean, R_mars, rover_lat, rover_lon); %#ok<ASGLU>

    [sep_min, imin] = min(sep_f);
    is_transit = sep_min < ang_rad_f(imin);

    results(end+1).t_days = t_fine(imin);                                       %#ok<SAGROW>
    results(end).JD = JD0 + t_fine(imin);
    results(end).sep_min_deg = sep_min;
    results(end).phobos_radius_deg = ang_rad_f(imin);
    results(end).is_transit = is_transit;
    results(end).t_fine = t_fine; results(end).sep_f = sep_f; results(end).ang_rad_f = ang_rad_f;
    results(end).elev_ph_f = elev_ph_f; results(end).elev_ea_f = elev_ea_f;
    results(end).az_ph_f = az_ph_f; results(end).az_ea_f = az_ea_f;
    results(end).imin = imin;
end

%% ---------------- 10. Report results ----------------
fprintf('\n=== Refined Close-Approach Results ===\n');
best_i = 0; best_sep = inf;
for c = 1:length(results)
    dt_c = datetime(results(c).JD - 1721058.5, 'ConvertFrom','datenum');
    fprintf('%d) %s | min separation = %.3f deg (%.1f arcmin) | Phobos radius = %.3f deg | TRANSIT: %s\n', ...
        c, datestr(dt_c), results(c).sep_min_deg, results(c).sep_min_deg*60, ...
        results(c).phobos_radius_deg, string(results(c).is_transit));
    if results(c).sep_min_deg < best_sep
        best_sep = results(c).sep_min_deg; best_i = c;
    end
end

if isempty(results)
    fprintf('No close approaches found in this window - try extending search_days.\n');
else
    fprintf('\n>>> BEST EVENT: #%d, separation %.3f deg, transit = %s\n', ...
        best_i, results(best_i).sep_min_deg, string(results(best_i).is_transit));
end

%% ================= FIGURE 1: full-window separation timeline (+ zoom) =================
figure('Color','w','Position',[50 50 950 750]);

subplot(2,1,1);
plot(t_days, sep_deg, 'b-'); hold on;
plot(t_days, ang_radius_phobos_deg, 'r--');
xlabel(sprintf('days since %s', datestr(date_start)));
ylabel('angular separation [deg]');
legend('Phobos-Earth separation','Phobos angular radius (transit threshold)');
title('Full 90-Day Search (each spike = one Phobos pass, ~3/day)');
grid on;

subplot(2,1,2);
zoom_days = 3;                                   % show first few days so individual passes are readable
mask_zoom = t_days <= zoom_days;
plot(t_days(mask_zoom), sep_deg(mask_zoom), 'b-', 'LineWidth',1.2); hold on;
plot(t_days(mask_zoom), ang_radius_phobos_deg(mask_zoom), 'r--', 'LineWidth',1.2);
xlabel(sprintf('days since %s', datestr(date_start)));
ylabel('angular separation [deg]');
title(sprintf('Zoomed: First %d Days (individual Phobos passes visible)', zoom_days));
grid on;

sgtitle('Phobos-Earth Angular Separation, Jezero Crater');
saveas(gcf, fullfile(outdir,'fig1_separation_timeline.png'));

%% ================= FIGURE 1B: sky map (Phobos & Earth paths across the Martian sky) =====
% Polar "sky dome" projection: zenith at center, horizon at the outer ring,
% azimuth measured clockwise from North. Gaps in a track mean the object
% was below the horizon (not visible) at that time.
figure('Color','w','Position',[60 60 1100 550]);

subplot(1,2,1);
plot_sky_grid();
r_ph = 90 - elev_phobos(mask_zoom); r_ph(elev_phobos(mask_zoom)<0) = NaN;
r_ea = 90 - elev_earth(mask_zoom);  r_ea(elev_earth(mask_zoom)<0)  = NaN;
az_ph_z = az_phobos(mask_zoom); az_ea_z = az_earth(mask_zoom);
plot(r_ph.*sind(az_ph_z), r_ph.*cosd(az_ph_z), 'r-', 'LineWidth',1);
plot(r_ea.*sind(az_ea_z), r_ea.*cosd(az_ea_z), 'b-', 'LineWidth',2.5);
legend('Phobos path','Earth path','Location','southoutside','Orientation','horizontal');
title(sprintf('Multi-Pass Sky Track (first %d days)\nPhobos loops fast; Earth drifts like a slow star', zoom_days));

subplot(1,2,2);
plot_sky_grid();
if ~isempty(results)
    r = results(best_i);
    r_ph_c = 90 - r.elev_ph_f; r_ph_c(r.elev_ph_f<0) = NaN;
    r_ea_c = 90 - r.elev_ea_f; r_ea_c(r.elev_ea_f<0) = NaN;
    plot(r_ph_c.*sind(r.az_ph_f), r_ph_c.*cosd(r.az_ph_f), 'r-', 'LineWidth',2);
    plot(r_ea_c.*sind(r.az_ea_f), r_ea_c.*cosd(r.az_ea_f), 'b-', 'LineWidth',3);
    % mark the moment of closest approach on both tracks
    plot(r_ph_c(r.imin)*sind(r.az_ph_f(r.imin)), r_ph_c(r.imin)*cosd(r.az_ph_f(r.imin)), ...
        'kp', 'MarkerFaceColor','y', 'MarkerSize',14);
    legend('Phobos path','Earth path','closest approach','Location','southoutside','Orientation','horizontal');
    title(sprintf('Closest-Approach Detail (\\pm1 hr)\nmin separation = %.3f deg, transit = %s', ...
        r.sep_min_deg, string(r.is_transit)));
end
sgtitle('Sky Map: Phobos & Earth Tracks as Seen from Jezero Crater');
saveas(gcf, fullfile(outdir,'fig1b_sky_map.png'));

%% ================= FIGURE 2: zoom on best event =================
if ~isempty(results)
    r = results(best_i);
    figure('Color','w','Position',[100 100 800 500]);
    t_hr = (r.t_fine - r.t_fine(1))*24;
    plot(t_hr, r.sep_f, 'b-', 'LineWidth',1.5); hold on;
    plot(t_hr, r.ang_rad_f, 'r--', 'LineWidth',1.5);
    xlabel('hours from window start'); ylabel('angular separation [deg]');
    legend('Phobos-Earth separation','Phobos angular radius');
    title(sprintf('Closest Approach Detail (min sep = %.3f deg, transit = %s)', ...
        r.sep_min_deg, string(r.is_transit)));
    grid on;
    saveas(gcf, fullfile(outdir,'fig2_closest_approach_zoom.png'));
end

%% ================= FIGURE 3: eclipse view (to-scale disk plot) =================
if ~isempty(results)
    r = results(best_i);
    figure('Color','w','Position',[150 150 550 550]);
    theta = linspace(0,2*pi,100);
    plot(r.phobos_radius_deg*60*cos(theta), r.phobos_radius_deg*60*sin(theta), 'r-','LineWidth',2); hold on;
    % Earth offset placed along +x at the separation distance (arcmin)
    plot(r.sep_min_deg*60, 0, 'bo', 'MarkerFaceColor','b','MarkerSize',10);
    axis equal; grid on;
    xlabel('arcmin'); ylabel('arcmin');
    legend('Phobos disk (true angular scale)','Earth (point source)');
    title('Eclipse-View: Phobos Disk vs. Earth Position at Closest Approach');
    saveas(gcf, fullfile(outdir,'fig3_eclipse_view.png'));
end

%% ================= FIGURE 4: visibility / elevation over full window (+ zoom) =================
figure('Color','w','Position',[200 200 950 750]);

subplot(2,1,1);
plot(t_days, elev_phobos, 'r-'); hold on;
plot(t_days, elev_earth, 'b-');
yline(0,'k--');
xlabel(sprintf('days since %s', datestr(date_start))); ylabel('elevation above horizon [deg]');
legend('Phobos elevation','Earth elevation','horizon');
title('Full 90-Day Search');
grid on;

subplot(2,1,2);
plot(t_days(mask_zoom), elev_phobos(mask_zoom), 'r-', 'LineWidth',1.2); hold on;
plot(t_days(mask_zoom), elev_earth(mask_zoom), 'b-', 'LineWidth',1.2);
yline(0,'k--');
xlabel(sprintf('days since %s', datestr(date_start))); ylabel('elevation above horizon [deg]');
title(sprintf('Zoomed: First %d Days', zoom_days));
grid on;

sgtitle('Visibility of Phobos & Earth from Jezero Crater');
saveas(gcf, fullfile(outdir,'fig4_visibility.png'));

%% ================= FIGURE 5: interplanetary context =================
figure('Color','w','Position',[250 250 600 600]);
theta = linspace(0,2*pi,200);
plot(cos(theta), sin(theta), 'b:'); hold on;
plot(1.524*cos(theta), 1.524*sin(theta), 'r:');
plot(R_earth_helio(1,1), R_earth_helio(1,2), 'bo','MarkerFaceColor','b','MarkerSize',8);
plot(R_mars_helio(1,1), R_mars_helio(1,2), 'ro','MarkerFaceColor','r','MarkerSize',8);
plot(0,0,'y.','MarkerSize',30);
axis equal; grid on;
xlabel('AU'); ylabel('AU');
legend('Earth orbit','Mars orbit','Earth','Mars','Sun','Location','bestoutside');
title('Interplanetary Geometry at Search-Window Start');
saveas(gcf, fullfile(outdir,'fig5_interplanetary_context.png'));

fprintf('\nAll figures saved to %s\n', outdir);