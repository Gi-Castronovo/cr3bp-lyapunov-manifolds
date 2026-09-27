function [Es, Eu, lambda_s, lambda_u] = monodromy_eigenvectors(X, T, mu)
%MONODROMY_EIGENVECTORS Stable/unstable directions of a periodic orbit at X.
%   Integrates the STM over one period starting from the orbit point X to
%   obtain the monodromy matrix M, then returns the eigenvectors associated
%   with the unstable (|lambda| > 1) and stable (|lambda| < 1) eigenvalues.
%   Both are normalised on their POSITION part, so that a perturbation
%   epsilon * E is a true displacement of epsilon in configuration space.

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[~, Y] = ode113(@(t, x) cr3bp_dynamics(t, x, mu), [0, T], [X(1:4); reshape(eye(4), 16, 1)], opts);
M = reshape(Y(end, 5:20), 4, 4);

[V, D] = eig(M);
lambda = real(diag(D));      % remove negligible imaginary round-off
V      = real(V);

[~, iu] = max(abs(lambda));  lambda_u = lambda(iu);  Eu = V(:, iu);
[~, is] = min(abs(lambda));  lambda_s = lambda(is);  Es = V(:, is);

Eu = Eu / norm(Eu(1:2));
Es = Es / norm(Es(1:2));
end
