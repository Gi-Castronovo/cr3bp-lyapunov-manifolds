%% LYAPUNOV ORBIT FAMILIES IN THE EARTH-MOON PCR3BP
% Part 1 of the project:
%   1. Linear analysis around L1 and L2 -> first-guess periodic orbits
%   2. Single-shooting differential correction (STM-based Newton method)
%   3. Natural-parameter continuation -> L1 and L2 Lyapunov families
%   4. Bisection -> pairs of L1/L2 orbits sharing the same Jacobi constant
%
% Output: data/lyapunov_pairs.mat, used by manifold_transfers.m
% Author: Girolamo Castronovo - Sapienza University of Rome, A.Y. 2025/26

clear; clc; close all;
addpath('functions');

SAVE_FIGURES = false;         % true -> export PNGs to /figures

%% System parameters (nondimensional units: DU = Earth-Moon distance)
GM_moon  = 4.902e3;            % km^3/s^2
GM_earth = 3.986e5;            % km^3/s^2
DU       = 384400;             % km
mu       = GM_moon/(GM_earth + GM_moon);

[xL1, xL2] = collinear_points(mu);
fprintf('mu = %.6f  |  xL1 = %.6f  |  xL2 = %.6f\n', mu, xL1, xL2);

amplitude = 50/DU;             % initial displacement from the equilibrium (50 km)
C_min     = 2.95;              % continuation stops below this Jacobi constant

%% 1. Linear guess around L1 and L2
[Xg_L1, Tg_L1] = linear_guess(xL1, mu, amplitude);
[Xg_L2, Tg_L2] = linear_guess(xL2, mu, amplitude);

Y = propagate(Xg_L1, Tg_L1, mu);
fig = figure('Name', 'Linear guess L1');
plot(Y(:,1), Y(:,2), 'b', 'LineWidth', 1.5); axis equal; grid on;
xlabel('x (nd)'); ylabel('y (nd)'); title('L1: propagation of the linear guess');
save_figure(fig, 'linear_guess_L1', SAVE_FIGURES);

Y = propagate(Xg_L2, Tg_L2, mu);
fig = figure('Name', 'Linear guess L2');
plot(Y(:,1), Y(:,2), 'r', 'LineWidth', 1.5); axis equal; grid on;
xlabel('x (nd)'); ylabel('y (nd)'); title('L2: propagation of the linear guess');
save_figure(fig, 'linear_guess_L2', SAVE_FIGURES);

%% 2. Differential correction (x0, v0 and T are all corrected)
[X_L1, T_L1] = single_shooting(Xg_L1, Tg_L1, mu, 'free_x');
[X_L2, T_L2] = single_shooting(Xg_L2, Tg_L2, mu, 'free_x');

plot_corrected_orbit(X_L1, T_L1, xL1, 'L1', 'b', mu, SAVE_FIGURES);
plot_corrected_orbit(X_L2, T_L2, xL2, 'L2', 'r', mu, SAVE_FIGURES);

%% 3. Natural-parameter continuation
% Each row of a family: [x0, v0, T, C]. Step sizes tuned by trial and error.
family_L1 = continue_family(X_L1, T_L1, -59e-5,  150, C_min, mu);
family_L2 = continue_family(X_L2, T_L2, 1.96e-4, 600, C_min, mu);

fprintf('L1 family: %d orbits, C in [%.4f, %.4f]\n', size(family_L1,1), min(family_L1(:,4)), max(family_L1(:,4)));
fprintf('L2 family: %d orbits, C in [%.4f, %.4f]\n', size(family_L2,1), min(family_L2(:,4)), max(family_L2(:,4)));

plot_family(family_L1, xL1, 'L1', mu, SAVE_FIGURES);
plot_family(family_L2, xL2, 'L2', mu, SAVE_FIGURES);

%% 4. Energy matching by bisection (master: L1 family, slave: L2 family)
% Every L1 orbit with C inside C_window gets an L2 "twin" with the same C.
C_window = [3.1370, 3.1493];
tol_C    = 1e-8;

C_masters = family_L1(family_L1(:,4) >= C_window(1) & family_L1(:,4) <= C_window(2), 4);

% pairs: [x0_L2, v0_L2, T_L2, C_L2, x0_L1, v0_L1, T_L1, C_L1]
pairs = zeros(numel(C_masters), 8);

for p = 1:numel(C_masters)
    C_target = C_masters(p);

    % Bracket C_target between two consecutive L2 orbits (C decreases along the family)
    k = find(family_L2(:,4) < C_target, 1, 'first');
    upper = family_L2(k,   :);     % C < C_target
    lower = family_L2(k-1, :);     % C > C_target

    C_mid = C_target + 1;
    while abs(C_mid - C_target) > tol_C
        x_mid = (lower(1) + upper(1))/2;
        v_mid = (lower(2) + upper(2))/2;
        T_mid = (lower(3) + upper(3))/2;

        % The interpolated state is not periodic: re-correct it
        [X_mid, T_mid] = single_shooting([x_mid, 0, 0, v_mid], T_mid, mu, 'fixed_x');
        C_mid = jacobi_constant(X_mid, mu);

        if C_mid > C_target
            lower = [X_mid(1), X_mid(4), T_mid, C_mid];
        elseif C_mid < C_target
            upper = [X_mid(1), X_mid(4), T_mid, C_mid];
        end
    end

    master   = family_L1(family_L1(:,4) == C_target, :);
    pairs(p, :) = [X_mid(1), X_mid(4), T_mid, C_mid, master];
end

fprintf('Matched %d L1/L2 orbit pairs.\n', size(pairs, 1));
if ~exist('data', 'dir'), mkdir('data'); end
save(fullfile('data', 'lyapunov_pairs.mat'), 'pairs', 'mu');

%% Plot of the matched pairs
fig = figure('Name', 'L1-L2 orbit pairs');
hold on; grid on; axis equal;
xlabel('x (nd)'); ylabel('y (nd)');
title('L1/L2 Lyapunov orbit pairs with equal Jacobi constant');

theta = linspace(0, 2*pi, 100);
R_moon = 1737/DU;
fill((1 - mu) + R_moon*cos(theta), R_moon*sin(theta), [0.6 0.6 0.6], 'EdgeColor', 'k', 'DisplayName', 'Moon');
plot(xL1, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'L1');
plot(xL2, 0, 'w+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'L2');

colors = lines(size(pairs, 1));
for p = 1:size(pairs, 1)
    Y2 = propagate([pairs(p,1), 0, 0, pairs(p,2)], pairs(p,3), mu);
    Y1 = propagate([pairs(p,5), 0, 0, pairs(p,6)], pairs(p,7), mu);
    plot(Y1(:,1), Y1(:,2), 'Color', colors(p,:), 'LineWidth', 2, 'DisplayName', sprintf('C = %.4f', pairs(p,4)));
    plot(Y2(:,1), Y2(:,2), 'Color', colors(p,:), 'LineWidth', 2, 'HandleVisibility', 'off');
end
legend('show', 'Location', 'best');
save_figure(fig, 'orbit_pairs', SAVE_FIGURES);


%% ===================== Local functions =====================

function [X_guess, T_guess] = linear_guess(xL, mu, amplitude)
% First-guess periodic orbit from the centre eigenvalue of the linearised
% dynamics around the collinear point (saddle x centre topology).
[V, D] = eig(cr3bp_jacobian(xL, 0, mu));
lambda = diag(D);
i = find(imag(lambda) > 1e-6, 1);       % positive imaginary eigenvalue
omega   = imag(lambda(i));
w       = V(:, i) / V(1, i);            % normalised on the x component
X_guess = [xL; 0; 0; 0] + amplitude*real(w);
T_guess = 2*pi/omega;
end

function Y = propagate(X0, T, mu)
% Propagates a state over one period (tight tolerances: the orbits are unstable)
opts = odeset('RelTol', 3e-14, 'AbsTol', 1e-14);
[~, Y] = ode113(@(t, X) cr3bp_dynamics(t, X, mu), [0, T], [X0(:); reshape(eye(4), 16, 1)], opts);
end

function plot_corrected_orbit(X0, T, xL, name, color, mu, save_flag)
Y   = propagate(X0, T, mu);
gap = norm(Y(end, 1:2) - Y(1, 1:2));   % closure error after one period
fprintf('%s: corrected period T = %.6f, closure error = %.3e\n', name, T, gap);
fig = figure('Name', ['Corrected orbit ' name]);
plot(Y(:,1), Y(:,2), color, 'LineWidth', 2); hold on;
plot(xL, 0, 'kx', 'MarkerSize', 8);
axis equal; grid on; xlabel('x (nd)'); ylabel('y (nd)');
title(sprintf('%s: corrected orbit (closure error %.2e)', name, gap));
save_figure(fig, ['corrected_orbit_' name], save_flag);
end

function family = continue_family(X, T, step, n_max, C_min, mu)
% Natural-parameter continuation in x0; v0 and T are re-corrected each step.
family = [X(1), X(4), T, jacobi_constant(X, mu)];
for k = 1:n_max
    [X, T] = single_shooting([X(1) + step, 0, 0, X(4)], T, mu, 'fixed_x');
    C = jacobi_constant(X, mu);
    if C < C_min
        fprintf('C = %.4f below C_min: continuation stopped.\n', C);
        break
    end
    family = [family; X(1), X(4), T, C]; %#ok<AGROW>
end
end

function plot_family(family, xL, name, mu, save_flag)
fig = figure('Name', ['Lyapunov family ' name]);
hold on; grid on; axis equal;
xlabel('x (nd)'); ylabel('y (nd)');
title(['Lyapunov family around ' name]);
plot(xL, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2);

n = size(family, 1);
colors = jet(n);
fprintf('Plotting %d %s orbits...\n', n, name);
for k = 1:n
    Y = propagate([family(k,1), 0, 0, family(k,2)], family(k,3), mu);
    plot(Y(:,1), Y(:,2), 'Color', colors(k,:), 'LineWidth', 1.5);
end
colormap(jet);
cb = colorbar;
cb.Label.String = 'Position along the family (normalised)';
save_figure(fig, ['lyapunov_family_' name], save_flag);
end

function save_figure(fig, name, save_flag)
if save_flag
    if ~exist('figures', 'dir'), mkdir('figures'); end
    exportgraphics(fig, fullfile('figures', [name '.png']), 'Resolution', 200);
end
end
