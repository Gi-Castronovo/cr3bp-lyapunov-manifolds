%% INVARIANT MANIFOLDS AND LOW-ENERGY TRANSFERS IN THE EARTH-MOON PCR3BP
% Part 2 of the project:
%   1. Stable/unstable directions of an L1 and an L2 Lyapunov orbit with the
%      same Jacobi constant (eigenvectors of the local monodromy matrix)
%   2. Globalisation of the manifold tubes and intersection with the
%      Poincare section x = 1 - mu (E -> M -> X configuration)
%   3. Propagation of ballistic trajectories for the four itineraries
%      X->M->E, E->M->X (transit) and X->M->X, E->M->E (non-transit)
%
% Input: data/lyapunov_pairs.mat, produced by lyapunov_families.m
% Author: Girolamo Castronovo - Sapienza University of Rome, A.Y. 2025/26

clear; clc; close all;
addpath('functions');

SAVE_FIGURES = true;          % true -> export PNGs to /figures
ITINERARY    = 'EMX';          % trajectory to propagate: 'XME', 'EMX', 'XMX', 'EME'

%% Parameters
load(fullfile('data', 'lyapunov_pairs.mat'), 'pairs');   % [x0_L2 v0_L2 T_L2 C_L2 x0_L1 v0_L1 T_L1 C_L1]

GM_moon  = 4.902e3;
GM_earth = 3.986e5;
DU       = 384400;                         % km
mu       = GM_moon/(GM_earth + GM_moon);
x_moon   = 1 - mu;

pair_idx = 3;                              % selected L1/L2 pair
C        = pairs(pair_idx, 4);             % common Jacobi constant
epsilon  = 50/DU;                          % manifold displacement (50 km)
N_points = 100;                            % points along each periodic orbit
fprintf('Selected pair %d: C = %.6f\n', pair_idx, C);

X0_L2 = [pairs(pair_idx,1); 0; 0; pairs(pair_idx,2)];   T_L2 = pairs(pair_idx,3);
X0_L1 = [pairs(pair_idx,5); 0; 0; pairs(pair_idx,6)];   T_L1 = pairs(pair_idx,7);

%% 1. Discretisation of the orbits and eigenvectors of the local monodromy matrix
opts_orb = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[~, Y_orb_L2] = ode113(@(t,x) cr3bp_dynamics(t,x,mu), linspace(0, T_L2, N_points+1), [X0_L2; reshape(eye(4),16,1)], opts_orb);
[~, Y_orb_L1] = ode113(@(t,x) cr3bp_dynamics(t,x,mu), linspace(0, T_L1, N_points+1), [X0_L1; reshape(eye(4),16,1)], opts_orb);

fprintf('Computing monodromy eigenvectors (%d points per orbit)...\n', N_points);
[Es_L2, Eu_L2] = eigenvectors_along_orbit(Y_orb_L2, T_L2, N_points, mu);
[Es_L1, Eu_L1] = eigenvectors_along_orbit(Y_orb_L1, T_L1, N_points, mu);

% Visual check of the perturbed initial conditions (the two branches must
% form two ordered rings around the periodic orbit)
plot_initial_conditions(Y_orb_L2, Eu_L2, epsilon, N_points, 'L2', SAVE_FIGURES);
plot_initial_conditions(Y_orb_L1, Eu_L1, epsilon, N_points, 'L1', SAVE_FIGURES);

%% 2. Manifold tubes and Poincare section x = 1 - mu
opts_tube    = odeset('RelTol', 1e-9,  'AbsTol', 1e-9,  'Events', @(t,y) event_box(t,y));
opts_section = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,y) event_vertical_section(t,y,x_moon,mu));

[X_mesh, Y_mesh, C_mesh] = zero_velocity_grid(mu);

fig_man = figure('Name', 'Invariant manifolds');
hold on; grid on; axis equal;
xlabel('x (nd)'); ylabel('y (nd)');
title('Invariant manifolds and Poincaré section x = 1 - \mu');
plot_forbidden_region(X_mesh, Y_mesh, C_mesh, C);

map_L2_S = [];     % L2 stable manifold (interior branch) on the section
map_L1_U = [];     % L1 unstable manifold (interior branch) on the section

% --- L2 orbit: exterior unstable branch (towards X), interior stable branch (from M)
T_tube = 15;
for k = 1:N_points
    Xk = Y_orb_L2(k, 1:4)';
    s  = branch_sign(Eu_L2(:,k), 'right');         % sign that moves the unstable start outwards
    X_U_ext = Xk + s*epsilon*Eu_L2(:,k);
    X_S_int = Xk + s*epsilon*Es_L2(:,k);

    [~, Y_U] = ode113(@(t,x) cr3bp_dynamics(t,x,mu), [0,  T_tube], [X_U_ext; reshape(eye(4),16,1)], opts_tube);
    [~, Y_S, ~, ye, ie] = ode113(@(t,x) cr3bp_dynamics(t,x,mu), [0, -T_tube], [X_S_int; reshape(eye(4),16,1)], opts_section);
    plot(Y_U(:,1), Y_U(:,2), 'c', 'LineWidth', 0.5);
    plot(Y_S(:,1), Y_S(:,2), 'b', 'LineWidth', 0.5);

    if ~isempty(ie) && ie(end) == 1 && ye(end,2) < 0      % lower half of the section
        map_L2_S = [map_L2_S; ye(end, 1:4)]; %#ok<AGROW>
    end
end

% --- L1 orbit: exterior stable branch (from E), interior unstable branch (towards M)
T_tube = 5;
for k = 1:N_points
    Xk = Y_orb_L1(k, 1:4)';
    s  = branch_sign(Es_L1(:,k), 'left');          % sign that moves the stable start towards the Earth
    X_S_ext = Xk + s*epsilon*Es_L1(:,k);
    X_U_int = Xk + s*epsilon*Eu_L1(:,k);

    [~, Y_S] = ode113(@(t,x) cr3bp_dynamics(t,x,mu), [0, -T_tube], [X_S_ext; reshape(eye(4),16,1)], opts_tube);
    [~, Y_U, ~, ye, ie] = ode113(@(t,x) cr3bp_dynamics(t,x,mu), [0,  T_tube], [X_U_int; reshape(eye(4),16,1)], opts_section);
    plot(Y_S(:,1), Y_S(:,2), 'r', 'LineWidth', 0.5);
    plot(Y_U(:,1), Y_U(:,2), 'y', 'LineWidth', 0.5);

    if ~isempty(ie) && ie(end) == 1 && ye(end,2) < 0
        map_L1_U = [map_L1_U; ye(end, 1:4)]; %#ok<AGROW>
    end
end

plot(Y_orb_L2(:,1), Y_orb_L2(:,2), 'w', 'LineWidth', 2, 'HandleVisibility', 'off');
plot(Y_orb_L1(:,1), Y_orb_L1(:,2), 'w', 'LineWidth', 2, 'HandleVisibility', 'off');
plot([x_moon, x_moon], [-0.12, 0.12], 'g', 'LineWidth', 4, 'HandleVisibility', 'off');

% Legend built from dummy handles (hundreds of lines were plotted above)
h = [plot(x_moon, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'DisplayName', 'Moon'), ...
     fill(nan(1,3), nan(1,3), [0.35 0.35 0.35], 'EdgeColor', 'k', 'DisplayName', 'Forbidden region (2U < C)'), ...
     plot(nan, nan, 'w', 'LineWidth', 2, 'DisplayName', 'Lyapunov orbits'), ...
     plot(nan, nan, 'r', 'LineWidth', 2, 'DisplayName', 'L1 stable W^S'), ...
     plot(nan, nan, 'y', 'LineWidth', 2, 'DisplayName', 'L1 unstable W^U'), ...
     plot(nan, nan, 'b', 'LineWidth', 2, 'DisplayName', 'L2 stable W^S'), ...
     plot(nan, nan, 'c', 'LineWidth', 2, 'DisplayName', 'L2 unstable W^U'), ...
     plot(nan, nan, 'g', 'LineWidth', 4, 'DisplayName', 'Poincaré section x = 1 - \mu')];
legend(h, 'Location', 'northeastoutside');
save_figure(fig_man, 'manifolds_EMX', SAVE_FIGURES);

%% 3. Poincare maps: every point inside the overlap of the two curves is a
%     state lying on both manifolds, i.e. a ballistic L1 -> L2 connection
fig = figure('Name', 'Poincare map (y, v)'); hold on; grid on; box on;
plot(map_L1_U(:,2), map_L1_U(:,4), 'y-', 'DisplayName', 'L1 unstable');
plot(map_L2_S(:,2), map_L2_S(:,4), 'b-', 'DisplayName', 'L2 stable');
xlabel('y (nd)'); ylabel('v (nd)'); title('Poincaré map at x = 1 - \mu: (y, v)');
legend('Location', 'best');
save_figure(fig, 'poincare_map_y_v', SAVE_FIGURES);

fig = figure('Name', 'Poincare map (y, u)'); hold on; grid on; box on;
plot(map_L1_U(:,2), map_L1_U(:,3), 'y-', 'DisplayName', 'L1 unstable');
plot(map_L2_S(:,2), map_L2_S(:,3), 'b-', 'DisplayName', 'L2 stable');
xlabel('y (nd)'); ylabel('u (nd)'); title('Poincaré map at x = 1 - \mu: (y, u)');
legend('Location', 'best');
save_figure(fig, 'poincare_map_y_u', SAVE_FIGURES);

%% 4. Ballistic trajectories for the four itineraries
% Initial conditions read from the Poincare maps. Position and velocity
% direction are prescribed; the speed follows from the Jacobi integral.
%              x        y      u-dir  v-dir   description
itineraries = {'XME', x_moon,  0.06,  -1,  0,  'X \rightarrow M \rightarrow E (transit)';
               'EMX', x_moon, -0.06,   1,  0,  'E \rightarrow M \rightarrow X (transit)';
               'XMX', 0.96,    0,      0, -1,  'X \rightarrow M \rightarrow X (non-transit)';
               'EME', 1.02,    0,      0,  1,  'E \rightarrow M \rightarrow E (non-transit)'};

row = find(strcmp(itineraries(:,1), ITINERARY));
if isempty(row), error('Unknown ITINERARY "%s".', ITINERARY); end
[x0, y0, du, dv, label] = itineraries{row, 2:6};

speed = sqrt(2*effective_potential(x0, y0, mu) - C);
X0    = [x0; y0; du*speed; dv*speed];
fprintf('%s initial state: [%.5f, %.5f, %.5f, %.5f]\n', ITINERARY, X0);

opts_final = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,y) event_box(t,y));
T_half = 2.5;
[~, Y_back] = ode113(@(t,x) cr3bp_dynamics(t,x,mu), [0, -T_half], [X0; reshape(eye(4),16,1)], opts_final);
[~, Y_fwd]  = ode113(@(t,x) cr3bp_dynamics(t,x,mu), [0,  T_half], [X0; reshape(eye(4),16,1)], opts_final);

fig = figure('Name', ['Trajectory ' ITINERARY]);
hold on; grid on; axis equal;
xlabel('x (nd)'); ylabel('y (nd)'); title(label);
plot_forbidden_region(X_mesh, Y_mesh, C_mesh, C);
plot(Y_orb_L2(:,1), Y_orb_L2(:,2), 'w', 'LineWidth', 2);
plot(Y_orb_L1(:,1), Y_orb_L1(:,2), 'w', 'LineWidth', 2);
plot(Y_back(:,1), Y_back(:,2), 'g', 'LineWidth', 2);
plot(Y_fwd(:,1),  Y_fwd(:,2),  'g', 'LineWidth', 2);
save_figure(fig, ['trajectory_' ITINERARY], SAVE_FIGURES);


%% ===================== Local functions =====================

function [Es, Eu] = eigenvectors_along_orbit(Y_orb, T, N, mu)
% Stable/unstable eigenvectors at N points of the orbit (4xN arrays)
Es = zeros(4, N);  Eu = zeros(4, N);
for k = 1:N
    [Es(:,k), Eu(:,k)] = monodromy_eigenvectors(Y_orb(k, 1:4)', T, mu);
end
end

function s = branch_sign(E, side)
% Sign of the perturbation that moves the start point to the requested side
% ('right' = larger x, 'left' = smaller x)
if strcmp(side, 'right'), s = 1 - 2*(E(1) <= 0);
else,                     s = 1 - 2*(E(1) >= 0);
end
end

function plot_initial_conditions(Y_orb, Eu, epsilon, N, name, save_flag)
starts_pos = zeros(N, 4);  starts_neg = zeros(N, 4);
for k = 1:N
    s = branch_sign(Eu(:,k), 'right');
    starts_pos(k,:) = Y_orb(k,1:4) + s*epsilon*Eu(:,k)';
    starts_neg(k,:) = Y_orb(k,1:4) - s*epsilon*Eu(:,k)';
end
fig = figure('Name', ['Manifold initial conditions ' name]);
hold on; axis equal; grid on;
plot(Y_orb(:,1), Y_orb(:,2), 'k', 'LineWidth', 2, 'DisplayName', ['Orbit ' name]);
plot(starts_pos(:,1), starts_pos(:,2), 'r.', 'MarkerSize', 8, 'DisplayName', 'Start +');
plot(starts_neg(:,1), starts_neg(:,2), 'b.', 'MarkerSize', 8, 'DisplayName', 'Start -');
xlabel('x (nd)'); ylabel('y (nd)'); legend;
title(['Manifold initial conditions around the ' name ' orbit']);
save_figure(fig, ['manifold_initial_conditions_' name], save_flag);
end

function U = effective_potential(x, y, mu)
r1 = sqrt((x + mu).^2 + y.^2);
r2 = sqrt((x - 1 + mu).^2 + y.^2);
U  = (x.^2 + y.^2)/2 + (1 - mu)./r1 + mu./r2;
end

function [X_mesh, Y_mesh, C_mesh] = zero_velocity_grid(mu)
[X_mesh, Y_mesh] = meshgrid(linspace(-2.5, 2.5, 800), linspace(-1.5, 1.5, 800));
C_mesh = 2*effective_potential(X_mesh, Y_mesh, mu);
end

function plot_forbidden_region(X_mesh, Y_mesh, C_mesh, C)
% Motion is only possible where v^2 = 2U - C >= 0: the region 2U < C is
% forbidden. contourf fills values ABOVE the level, hence the minus sign.
contourf(X_mesh, Y_mesh, -C_mesh, [-C, -C], 'LineColor', 'none', 'FaceColor', [0.35 0.35 0.35], 'HandleVisibility', 'off');
contour(X_mesh, Y_mesh, C_mesh, [C, C], 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
end

function save_figure(fig, name, save_flag)
if save_flag
    if ~exist('figures', 'dir'), mkdir('figures'); end
    exportgraphics(fig, fullfile('figures', [name '.png']), 'Resolution', 200);
end
end
