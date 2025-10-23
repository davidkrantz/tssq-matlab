% EXAMPLE_OPEN_ARC_MONOMIAL
% Compare standard singularity swap quadrature (SSQ) translated SSQ (TSSQ)
% on a parabolic open arc using standard and translated monomial bases for 
% for the Stokeslet contribution of the slender-body kernel.
%
% AUTHOR: David Krantz (davkra@kth.se), September 2025

clear; close all; clc; rng(0);

% Options
opts = default_options('monomial');
opts.tol = 1e-4; % desired tolerance
opts.slender_eps = 0; % only Stokeslet: S(r) = I/|r| + (r r^T)/|r|^3

% Curve: parabolic arc: gamma(t) = (t, a*t^2, 0), t in [-1,1]
a = 0.9;
s_of_t = @(t) -1 + 2*t; % map [0,1] -> [-1,1]
curve.x  = @(t) s_of_t(t.');
curve.y  = @(t) a * (s_of_t(t.')).^2;
curve.z  = @(t) 0*t.';
curve.xp = @(t) 2; % derivatives wrt t
curve.yp = @(t) 2*a * s_of_t(t.') * 2;
curve.zp = @(t) 0*t.';
curve.s  = @(t) sqrt(curve.xp(t).^2 + curve.yp(t).^2 + curve.zp(t).^2);

% Artifical layer density: [1,0,0] for simplicity
density.f1 = @(t) 1 + 0*t.';
density.f2 = @(t) 0*t.';
density.f3 = @(t) 0*t.';

% Targets: a single line from curve
Nt = 50; % number of targets
tline = 0.3; % where on curve [0,1] the normal starts
p = [curve.x(tline); curve.y(tline); curve.z(tline)]; % point on curve
nvec = [1;1;1];
nunit = nvec/norm(nvec); % unit normal vector
distv = logspace(-8,0,Nt); % distances along the line
targets = p + nunit .* distv; % 3 x Nt

% Reference: adaptive quadrature
disp('* Reference');
opts_ref = opts;
opts_ref.tol = 5e-14; opts_ref.nquad = opts.nquad+2; opts_ref.Hlim = 1;
u1_ref = adaptive_quadrature_wrapper(curve,density,targets,opts_ref); 

% SSQ
disp(' '); disp('* SSQ');
u1_ssq = ssq_sbt(curve,density,targets,opts);

% TSSQ
disp(' '); disp('* TSSQ');
u1_tssq = tssq_sbt(curve,density,targets,opts);

% Relative error
err_ssq = abs(u1_ref-u1_ssq)./norm(u1_ref,inf);
err_tssq = abs(u1_ref-u1_tssq)./norm(u1_ref,inf);

% Plot curve and targets
figure;
tt = linspace(0,1,200);
plot3(curve.x(tt), curve.y(tt), curve.z(tt), 'k-', 'LineWidth',1.2);
hold on;
plot3(targets(1,:), targets(2,:), targets(3,:), 'r.');
axis equal; grid on;
xlabel('x'); ylabel('y'); zlabel('z');
title('Open arc and target line');

% Plot error vs distance to curve
figure;
loglog(distv,err_ssq,'o-'); hold on;
loglog(distv,err_tssq,'x-');
yline(opts.tol,'k','tolerance');
xlabel('Distance to curve'); ylabel('Relative error');
xlim([min(distv),max(distv)]);
grid on; legend('SSQ','TSSQ','Location','southwest');
title('Error vs distance (Monomial basis)');

alignfigs;
