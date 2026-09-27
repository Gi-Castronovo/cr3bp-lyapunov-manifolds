function dX = cr3bp_dynamics(~, X, mu)
%CR3BP_DYNAMICS Planar CR3BP equations of motion augmented with the STM.
%   dX = CR3BP_DYNAMICS(t, X, mu) returns the time derivative of the
%   augmented state X = [x; y; u; v; Phi(:)] (20x1), expressed in the
%   Earth-Moon synodic frame and in nondimensional units.
%   The 4x4 State Transition Matrix Phi evolves as dPhi/dt = A(x,y) * Phi.

x = X(1);  y = X(2);  u = X(3);  v = X(4);
Phi = reshape(X(5:20), 4, 4);

r1 = sqrt((x + mu)^2     + y^2);   % distance from the Earth
r2 = sqrt((x - 1 + mu)^2 + y^2);   % distance from the Moon

% Gradient of the effective potential U(x, y)
Ux = x - (1 - mu)*(x + mu)/r1^3 - mu*(x - 1 + mu)/r2^3;
Uy = y - (1 - mu)*y/r1^3        - mu*y/r2^3;

dPhi = cr3bp_jacobian(x, y, mu) * Phi;

dX = [u;
      v;
      2*v + Ux;
     -2*u + Uy;
      dPhi(:)];
end
