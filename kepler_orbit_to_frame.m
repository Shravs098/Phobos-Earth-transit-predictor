function [xr,yr,zr] = kepler_orbit_to_frame(a,e,I,Om,w,M)
% KEPLER_ORBIT_TO_FRAME  Constant orbital elements (a,e,I,Om,w), M may be
% a vector (varies with time). Returns Cartesian position components in
% the reference frame whose fundamental plane/axes the elements are
% defined against (e.g. Mars-equatorial-inertial frame for Phobos).
E = kepler_solve(M, e);
x_orb = a*(cos(E)-e);
y_orb = a*sqrt(1-e^2)*sin(E);
cosOm=cos(Om); sinOm=sin(Om); cosw=cos(w); sinw=sin(w); cosI=cos(I); sinI=sin(I);
xr = (cosOm*cosw - sinOm*sinw*cosI)*x_orb + (-cosOm*sinw - sinOm*cosw*cosI)*y_orb;
yr = (sinOm*cosw + cosOm*sinw*cosI)*x_orb + (-sinOm*sinw + cosOm*cosw*cosI)*y_orb;
zr = (sinw*sinI)*x_orb + (cosw*sinI)*y_orb;
end