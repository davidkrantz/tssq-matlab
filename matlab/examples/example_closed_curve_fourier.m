% EXAMPLE_CLOSED_CURVE_FOURIER
% Compare standard singularity swap quadrature (SSQ) and translated SSQ 
% (TSSQ) on a closed curve using the Fourier basis for the Stokeslet 
% contribution of the slender-body kernel.
%
% AUTHOR: David Krantz (davkra@kth.se), September 2025

clear; close all; clc; rng(0);

% Options
opts = default_options('fourier');
opts.tol = 1e-4; % desired tolerance
opts.slender_eps = 0; % only Stokeslet: S(r) = I/|r| + (r r^T)/|r|^3

% Curve: closed circle gamma(t) = (cos(2*pi*t), sin(2*pi*t), 0), t in [0,1)
curve_ref.x  = @(t) cos(2*pi*t.');
curve_ref.y  = @(t) sin(2*pi*t.');
curve_ref.z  = @(t) 0*t.';
curve_ref.xp = @(t) -2*pi*sin(2*pi*t.');
curve_ref.yp = @(t) 2*pi*cos(2*pi*t.');
curve_ref.zp = @(t) 0*t.';
curve_ref.s  = @(t) sqrt(curve_ref.xp(t).^2 + curve_ref.yp(t).^2 + ...
    curve_ref.zp(t).^2);
% current Fourier implementation of SSQ/TSSQ requires t in [0,2*pi)
curve.x  = @(t) cos(t.');
curve.y  = @(t) sin(t.');
curve.z  = @(t) 0*t.';
curve.xp = @(t) -sin(t.');
curve.yp = @(t) cos(t.');
curve.zp = @(t) 0*t.';
curve.s  = @(t) sqrt(curve.xp(t).^2 + curve.yp(t).^2 + curve.zp(t).^2);

% Artificial layer density: [1+cos(2*pi*t),0,0], t in [0,1)
density_ref.f1 = @(t) 1 + cos(2*pi*t.');
density_ref.f2 = @(t) 0*t.';
density_ref.f3 = @(t) 0*t.';
density.f1 = @(t) 1 + cos(t.'); % t in [0,2*pi)
density.f2 = @(t) 0*t.';
density.f3 = @(t) 0*t.';

% Targets: a single line outward from (1,0,0)
Nt = 50; % number of targets
tline = 11*pi/3; % parameter location for anchor point
p = [curve.x(tline); curve.y(tline); curve.z(tline)]; % point on curve
nunit = [1;0;1];  % outward at (1,0,1)
distv = logspace(-8,-1,Nt); % distances along the line
targets = p + nunit .* distv; % 3 x Nt

% Reference: adaptive quadrature (panel-based, high accuracy)
disp('* Reference');
opts_ref = opts;
opts_ref.tol = 5e-14; opts_ref.nquad = 2*opts.nquad+2;
u1_ref = adaptive_quadrature_wrapper(curve_ref,density_ref,targets, ...
    opts_ref);

% SSQ
disp(' '); disp('* SSQ');
[u1_ssq,~,~,stats_ssq] = ssq_sbt(curve,density,targets,opts);

% TSSQ
disp(' '); disp('* TSSQ');
[u1_tssq,~,~,stats_tssq] = tssq_sbt(curve,density,targets,opts);

% Relative error (normalized by reference infinity norm)
err_ssq  = abs(u1_ref - u1_ssq)  ./ norm(u1_ref,inf);
err_tssq = abs(u1_ref - u1_tssq) ./ norm(u1_ref,inf);

% Plot curve and targets
figure;
tt = linspace(0,2*pi,400);
plot3(curve.x(tt), curve.y(tt), curve.z(tt), 'k-', 'LineWidth',1.2);
hold on;
plot3(targets(1,:), targets(2,:), targets(3,:), 'r.');
axis equal; grid on;
xlabel('x'); ylabel('y'); zlabel('z');
title('Closed curve (unit circle) and target line');

% Plot error vs distance to curve
figure;
loglog(distv, err_ssq, 'o-'); hold on;
loglog(distv, err_tssq, 'x-');
yline(opts.tol,'k','tolerance');
xlabel('Distance to curve'); ylabel('Relative error');
xlim([min(distv),max(distv)]);
grid on; legend('SSQ','TSSQ','tolerance');
title('Error vs distance (Fourier basis)');

alignfigs;
