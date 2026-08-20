function [xr,yr,zr] = kepler_orbit_to_frame_vec(a,e,I,Om,w,M)
% KEPLER_ORBIT_TO_FRAME_VEC  Same as kepler_orbit_to_frame, but a,e,I,Om,w,M
% are all Nx1 arrays that vary in time (used for planetary heliocentric
% orbits, whose osculating elements drift slowly via secular rates).
E = kepler_solve(M, e);
x_orb = a.*(cos(E)-e);
y_orb = a.*sqrt(1-e.^2).*sin(E);
cosOm=cos(Om); sinOm=sin(Om); cosw=cos(w); sinw=sin(w); cosI=cos(I); sinI=sin(I);
xr = (cosOm.*cosw - sinOm.*sinw.*cosI).*x_orb + (-cosOm.*sinw - sinOm.*cosw.*cosI).*y_orb;
yr = (sinOm.*cosw + cosOm.*sinw.*cosI).*x_orb + (-sinOm.*sinw + cosOm.*cosw.*cosI).*y_orb;
zr = (sinw.*sinI).*x_orb + (cosw.*sinI).*y_orb;
end