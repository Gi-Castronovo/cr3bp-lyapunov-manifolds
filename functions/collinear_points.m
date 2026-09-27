function [xL1, xL2] = collinear_points(mu)
%COLLINEAR_POINTS x-coordinates of the L1 and L2 libration points.
%   Computed numerically (grad U = 0 on the x-axis) instead of using
%   tabulated values: Lyapunov orbits are highly unstable, so the
%   equilibrium must be consistent with the value of mu used everywhere.

options = optimset('TolX', 1e-14);
x_earth = -mu;
x_moon  = 1 - mu;

dUdx = @(x) x - (1 - mu)*(x + mu)./abs(x + mu).^3 ...
              - mu*(x - x_moon)./abs(x - x_moon).^3;

xL1 = fzero(dUdx, [x_earth + 1e-6, x_moon - 1e-6], options);   % between Earth and Moon
xL2 = fzero(dUdx, [x_moon + 1e-6, 2.0], options);              % beyond the Moon
end
