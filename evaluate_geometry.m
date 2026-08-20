function [sep_deg, ang_rad_deg, elev_ph, elev_ea, d_ph, d_ea, az_ph, az_ea] = evaluate_geometry( ...
        t_days, JD0, JD_J2000, earth_el0, earth_rate, mars_el0, mars_rate, AU_km, ...
        a_phobos, e_phobos, i_phobos, Om_phobos, w_phobos, M0_phobos, n_phobos, ...
        R_phobos_mean, R_mars, rover_lat, rover_lon)
% EVALUATE_GEOMETRY  Re-usable full-geometry evaluation (heliocentric
% orbits, Mars orientation, Phobos orbit, rover position, angular
% separation, and local sky azimuth/elevation) at an arbitrary set of
% times. Used for the coarse full-window scan, fine-resolution
% refinement near candidate transit windows, and sky-map plotting.

    JD = JD0 + t_days;
    d_since_J2000 = JD - JD_J2000;
    T_cent = d_since_J2000/36525;

    R_earth_helio = heliocentric_position(earth_el0, earth_rate, T_cent);
    R_mars_helio  = heliocentric_position(mars_el0,  mars_rate,  T_cent);

    alpha0 = (317.269 - 0.108*T_cent) * pi/180;
    delta0 = (54.432  - 0.061*T_cent) * pi/180;
    W = mod(176.630 + 350.89198226*d_since_J2000, 360) * pi/180;

    p_hat0 = [cos(delta0(1))*cos(alpha0(1)), cos(delta0(1))*sin(alpha0(1)), sin(delta0(1))];
    N_hat0 = [-sin(alpha0(1)), cos(alpha0(1)), 0];
    PxN0   = cross(p_hat0, N_hat0);
    Xeq = N_hat0; Yeq = PxN0; Zeq = p_hat0;

    Xbf = cos(W).*N_hat0 + sin(W).*PxN0;
    Ybf = cross(repmat(p_hat0,length(W),1), Xbf, 2);
    Zbf = repmat(p_hat0, length(W), 1);

    lat_r = rover_lat*pi/180; lon_r = rover_lon*pi/180;
    r_bf = [cos(lat_r)*cos(lon_r), cos(lat_r)*sin(lon_r), sin(lat_r)];
    Rover_pos = R_mars*( r_bf(1).*Xbf + r_bf(2).*Ybf + r_bf(3).*Zbf );
    zenith_hat = Rover_pos / R_mars;

    % local topocentric East/North unit vectors (for azimuth), rotating with Mars
    East_bf  = [-sin(lon_r), cos(lon_r), 0];
    North_bf = [-sin(lat_r)*cos(lon_r), -sin(lat_r)*sin(lon_r), cos(lat_r)];
    East_hat  = East_bf(1).*Xbf  + East_bf(2).*Ybf  + East_bf(3).*Zbf;
    North_hat = North_bf(1).*Xbf + North_bf(2).*Ybf + North_bf(3).*Zbf;

    M_phobos = mod(M0_phobos + n_phobos*t_days, 2*pi);
    [xr,yr,zr] = kepler_orbit_to_frame(a_phobos, e_phobos, i_phobos, Om_phobos, w_phobos, M_phobos);
    Phobos_pos = xr.*Xeq + yr.*Yeq + zr.*Zeq;

    d_ph = Phobos_pos - Rover_pos;
    d_ea = (R_earth_helio - R_mars_helio)*AU_km - Rover_pos;

    dist_ph = vecnorm(d_ph,2,2); dist_ea = vecnorm(d_ea,2,2);
    cos_sep = sum(d_ph.*d_ea,2)./(dist_ph.*dist_ea);
    sep_deg = acosd(min(1,max(-1,cos_sep)));
    ang_rad_deg = asind(R_phobos_mean./dist_ph);
    elev_ph = asind(sum(zenith_hat.*d_ph,2)./dist_ph);
    elev_ea = asind(sum(zenith_hat.*d_ea,2)./dist_ea);
    az_ph = mod(atan2d(sum(East_hat.*d_ph,2), sum(North_hat.*d_ph,2)), 360);
    az_ea = mod(atan2d(sum(East_hat.*d_ea,2), sum(North_hat.*d_ea,2)), 360);
end