function R_eq = heliocentric_position(el0, rate, T_cent)
% HELIOCENTRIC_POSITION  Heliocentric equatorial (J2000) position [AU]
% from low-precision Keplerian elements.
%   el0, rate : 1x6 = [a,e,I,L,long_peri,long_node] (AU, deg, deg/century)
%   T_cent    : Nx1 Julian centuries since J2000
%   R_eq      : Nx3 position in AU, J2000 mean equatorial frame
el = el0 + rate.*T_cent;             % Nx6 broadcast
a = el(:,1); e = el(:,2); I = el(:,3)*pi/180;
L = el(:,4)*pi/180; lp = el(:,5)*pi/180; lo = el(:,6)*pi/180;
Om = lo; w = lp - lo; M = mod(L - lp + pi, 2*pi) - pi;

[xr,yr,zr] = kepler_orbit_to_frame_vec(a,e,I,Om,w,M);   % ecliptic, AU

eps0 = 23.43928*pi/180;
x_eq = xr;
y_eq = yr*cos(eps0) - zr*sin(eps0);
z_eq = yr*sin(eps0) + zr*cos(eps0);
R_eq = [x_eq, y_eq, z_eq];
end