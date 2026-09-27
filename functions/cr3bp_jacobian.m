function A = cr3bp_jacobian(x, y, mu)
%CR3BP_JACOBIAN Jacobian of the planar CR3BP vector field at (x, y).
%   A = CR3BP_JACOBIAN(x, y, mu) returns the 4x4 matrix used both for the
%   linear analysis around the libration points and for STM propagation.

r1 = sqrt((x + mu)^2     + y^2);
r2 = sqrt((x - 1 + mu)^2 + y^2);

% Second derivatives of the effective potential U(x, y)
Uxx = 1 - (1 - mu)/r1^3 - mu/r2^3 ...
        + 3*(1 - mu)*(x + mu)^2/r1^5 + 3*mu*(x - 1 + mu)^2/r2^5;
Uyy = 1 - (1 - mu)/r1^3 - mu/r2^3 ...
        + 3*(1 - mu)*y^2/r1^5        + 3*mu*y^2/r2^5;
Uxy = 3*(1 - mu)*(x + mu)*y/r1^5 + 3*mu*(x - 1 + mu)*y/r2^5;

A = [ 0,   0,   1, 0;
      0,   0,   0, 1;
      Uxx, Uxy, 0, 2;
      Uxy, Uyy, -2, 0];
end
