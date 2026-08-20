function E = kepler_solve(M, e)
% KEPLER_SOLVE  Vectorized Newton-Raphson solution of Kepler's equation
%   E - e*sin(E) = M
M = mod(M + pi, 2*pi) - pi;
E = M;   % initial guess
for k = 1:10
    E = E - (E - e.*sin(E) - M)./(1 - e.*cos(E));
end
end