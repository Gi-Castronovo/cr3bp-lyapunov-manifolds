function C = jacobi_constant(X, mu)
%JACOBI_CONSTANT Jacobi integral C = 2U - (u^2 + v^2) of a planar CR3BP state.
%   X = [x, y, u, v] (row or column). Lower C means higher energy.

x = X(1);  y = X(2);  u = X(3);  v = X(4);

r1 = sqrt((x + mu)^2     + y^2);
r2 = sqrt((x - 1 + mu)^2 + y^2);
U  = (x^2 + y^2)/2 + (1 - mu)/r1 + mu/r2;

C = 2*U - (u^2 + v^2);
end
