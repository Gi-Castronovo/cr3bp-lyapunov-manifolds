function [X0_corr, T_corr] = single_shooting(X0, T, mu, mode)
%SINGLE_SHOOTING Differential corrector for symmetric planar Lyapunov orbits.
%   The orbit starts on the x-axis with a perpendicular velocity,
%   X0 = [x0, 0, 0, v0]. By the symmetry of the CR3BP it is periodic if it
%   crosses the x-axis perpendicularly at T/2, i.e. F = [y(T/2); u(T/2)] = 0.
%   Newton iterations use the STM to build the Jacobian DF of F.
%
%   mode = 'free_x'  : unknowns [x0; v0; T]  (used to refine the linear guess)
%   mode = 'fixed_x' : x0 is fixed, unknowns [v0; T]  (used in continuation)
%
%   Returns the corrected initial state (1x4 row) and the full period.

tol      = 1e-12;
max_iter = 50;
opts     = odeset('RelTol', 3e-14, 'AbsTol', 1e-14);
fixed_x  = strcmp(mode, 'fixed_x');

if fixed_x
    x0 = X0(1);
    z  = [X0(4); T];
else
    z  = [X0(1); X0(4); T];
end

F    = [1; 1];
iter = 0;
while norm(F) > tol
    iter = iter + 1;
    if iter > max_iter
        warning('single_shooting:noConvergence', ...
                'No convergence after %d iterations (|F| = %.2e).', max_iter, norm(F));
        break
    end

    if fixed_x
        Xi = [x0;   0; 0; z(1)];  Ti = z(2);
    else
        Xi = [z(1); 0; 0; z(2)];  Ti = z(3);
    end

    % Propagate state + STM for half a period
    [~, Y] = ode113(@(t, X) cr3bp_dynamics(t, X, mu), [0, Ti/2], [Xi; reshape(eye(4), 16, 1)], opts);
    Yf   = Y(end, :)';
    Phi  = reshape(Yf(5:20), 4, 4);
    dYf  = cr3bp_dynamics(0, Yf, mu);

    F     = [Yf(2); Yf(3)];             % [y(T/2); u(T/2)]
    dF_dT = 0.5*[Yf(4); dYf(3)];        % sensitivity to the period (evaluated at T/2)

    if fixed_x
        DF = [Phi(2,4), dF_dT(1);
              Phi(3,4), dF_dT(2)];
    else
        DF = [Phi(2,1), Phi(2,4), dF_dT(1);
              Phi(3,1), Phi(3,4), dF_dT(2)];
    end

    z = z - DF \ F;                     % Newton update
end

if fixed_x
    X0_corr = [x0,   0, 0, z(1)];  T_corr = z(2);
else
    X0_corr = [z(1), 0, 0, z(2)];  T_corr = z(3);
end
end
