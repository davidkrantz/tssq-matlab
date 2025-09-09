% LONG_FILAMENT  Example comparing SSQ and TSSQ on a tangled closed curve.
%
% This script reproduces results from the paper:
% - Evaluates Stokes slender-body potential
% - Compares SSQ and TSSQ in a monomial basis
% - Uses adaptive quadrature as reference
%
% AUTHOR: David Krantz (davkra@kth.se)

clear all;
close all;
format long;
rng(123);

opts = default_options('long_filament');
savefig = 0;

% Example parameters
tolv  = [1e-4;1e-6];
Neval = 1e3; % targets per distance
distv = fliplr(logspace(-8,-2,20)).'; % distances to test

% Setup curve and artificial layer density
curve = squiggle(); % x(t), y(t), z(t), xp(t), yp(t), zp(t), s(t)
density = struct('f1', @(t) curve.x(t), ...
                 'f2', @(t) curve.y(t), ...
                 'f3', @(t) curve.z(t));

% Target setup (random normals at fixed distances)
t = rand(1,Neval);
utang = [curve.xp(t); curve.yp(t); curve.zp(t)]./curve.s(t);
v = randn(3,Neval); v = v - utang.*sum(v.*utang,1);
v = v ./ sqrt(sum(v.^2,1));
X = []; Y = []; Z = [];
for d = distv.'
    X = [X curve.x(t) + d*v(1,:)];
    Y = [Y curve.y(t) + d*v(2,:)];
    Z = [Z curve.z(t) + d*v(3,:)];
end
targets = [X; Y; Z];

% Reference solution (adaptive quadrature)
disp('* Reference solution');
opts_ref = opts;
opts_ref.tol = 5e-14;
opts_ref.nquad = opts.nquad+2;
[uref1,uref2,uref3] = adaptive_quadrature_wrapper(curve,density,targets,opts_ref);

% Loop over tolerances
[adapquad_err_data,ssq_err_data,tssq_err_data] = deal(zeros(numel(distv),numel(tolv),3));
stats_tssq = cell(numel(tolv),1);
for ii = 1:numel(tolv)
    opts.tol = tolv(ii);

    % Adaptive quadrature (just for sanity check)
    disp(' ')
    disp('* Adaptive quadrature')
    [adquad1,adquad2,adquad3] = adaptive_quadrature_wrapper(curve,density,targets,opts);
    
    % Standard SSQ
    disp(' ')
    disp('* SSQ')
    [ussq1,ussq2,ussq3] = ssq_sbt(curve,density,targets,opts);

    % Translated SSQ (TSSQ)
    disp(' ')
    disp('* TSSQ')
    [utssq1,utssq2,utssq3,stats_tssq_tmp] = tssq_sbt(curve,density,targets,opts);
    stats_tssq{ii} = stats_tssq_tmp;

    % Compute errors
    [adquad_errmax,ssq_errmax,tssq_errmax] = deal(zeros(1,numel(distv)*Neval));
    for i = 1:numel(distv)
        idx = (i-1)*Neval+1:i*Neval;
        adquad_errmax(idx) = compute_error(uref1(idx), uref2(idx), uref3(idx), adquad1(idx), adquad2(idx), adquad3(idx));
        ssq_errmax(idx) = compute_error(uref1(idx), uref2(idx), uref3(idx), ussq1(idx), ussq2(idx), ussq3(idx));
        tssq_errmax(idx) = compute_error(uref1(idx), uref2(idx), uref3(idx), utssq1(idx), utssq2(idx), utssq3(idx));
        adapquad_err_data(i,ii,1) = min(adquad_errmax(idx));
        ssq_err_data(i,ii,1) = min(ssq_errmax(idx));
        tssq_err_data(i,ii,1) = min(tssq_errmax(idx));
        adapquad_err_data(i,ii,2) = mean(adquad_errmax(idx));
        ssq_err_data(i,ii,2) = mean(ssq_errmax(idx));
        tssq_err_data(i,ii,2) = mean(tssq_errmax(idx));
        adapquad_err_data(i,ii,3) = max(adquad_errmax(idx));
        ssq_err_data(i,ii,3) = max(ssq_errmax(idx));
        tssq_err_data(i,ii,3) = max(tssq_errmax(idx));
    end
end

% TSSQ timings
[kerevals_near,time_std_weights,time_ker_near,time_mod_weights,time_cancel_est,time_mod_vander_solve,time_mod_sigma_interp_a,time_mod_kernel_eval_a,time_mod_basis_int] = deal(zeros(numel(tolv),1));
for ii = 1:numel(tolv)
    kerevals_near(ii) = stats_tssq{ii}.kerevals_near;
    time_std_weights(ii) = stats_tssq{ii}.time_std_weights;
    time_ker_near(ii) = stats_tssq{ii}.time_ker_near;
    time_mod_weights(ii) = stats_tssq{ii}.time_mod_weights;
    time_cancel_est(ii) = stats_tssq{ii}.time_cancel_est;
    time_mod_basis_int(ii) = stats_tssq{ii}.time_mod_basis_int;
    time_mod_vander_solve(ii) = stats_tssq{ii}.time_mod_vander_solve;
    time_mod_sigma_interp_a(ii) = stats_tssq{ii}.time_mod_sigma_interp_a;
    time_mod_kernel_eval_a(ii) = stats_tssq{ii}.time_mod_kernel_eval_a;
end
table(tolv,kerevals_near,time_std_weights,time_mod_weights,time_cancel_est,time_ker_near,time_mod_basis_int,time_mod_vander_solve,time_mod_sigma_interp_a,time_mod_kernel_eval_a)

% Plots (hardcoded for paper example)
close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

FS = 15; % fontsize
LW = 1; % linewidth
MS = 8; % markersize
cols = colororder; % default Matlab color ordering
d2 = [distv; flipud(distv)]; % for filled region
falpha = 0.4; % opacity of filled region

figure; % to make figures same size
figure;

% tangled filament
figure('DefaultAxesFontSize',FS);
tj = adaptive_panelization(curve.s, opts.nquad, 1e-6);
xj = curve.x(tj); yj = curve.y(tj); zj = curve.z(tj);
plot3(xj, yj, zj, '.-k','MarkerSize',6);
hold on;
axis equal tight vis3d;
grid on;
box on;
xlabel('$x$','fontsize',FS,'interpreter','latex');
ylabel('$y$','fontsize',FS,'interpreter','latex');
zlabel('$z$','fontsize',FS,'interpreter','latex');

% tol = 1e-4
figure('DefaultAxesFontSize',FS);
loglog(distv,ssq_err_data(:,1,2),'Color',cols(1,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
hold on;
f = fill(d2,[ssq_err_data(:,1,1); flipud(ssq_err_data(:,1,3))],'g','FaceAlpha',falpha);
f.FaceColor = cols(1,:);
f.EdgeColor = cols(1,:);
loglog(distv,tssq_err_data(:,1,2),'Color',cols(2,:),'Marker','square','LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',cols(2,:));
f = fill(d2,[tssq_err_data(:,1,1); flipud(tssq_err_data(:,1,3))],'g','FaceAlpha',falpha);
f.FaceColor = cols(2,:);
f.EdgeColor = cols(2,:);
loglog(distv,5e-13./abs(distv).^2,'Color','k','Marker','none','LineStyle','--','LineWidth',LW,'MarkerSize',MS);
grid on;
xlabel('Distance to $\Gamma$, $d$','fontsize',FS,'interpreter','latex');
ylabel('Relative error','fontsize',FS,'interpreter','latex');
xlim([min(distv),max(distv)]);
legend('$\epsilon=10^{-4}$, SSQ','','$\epsilon=10^{-4}$, TSSQ','fontsize',FS,'interpreter','latex');
annotation('textarrow',[0.4 0.35],[0.78 0.67],'String','$\mathcal{O}(1/d^2)$','fontsize',FS,'interpreter','latex')
ylim([1e-15,1e8]);

% tol = 1e-6
figure('DefaultAxesFontSize',FS);
loglog(distv,ssq_err_data(:,2,2),'Color',cols(3,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
hold on;
f = fill(d2,[ssq_err_data(:,2,1); flipud(ssq_err_data(:,2,3))],'g','FaceAlpha',falpha);
f.FaceColor = cols(3,:);
f.EdgeColor = cols(3,:);
loglog(distv,tssq_err_data(:,2,2),'Color',cols(4,:),'Marker','square','LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',cols(4,:));
f = fill(d2,[tssq_err_data(:,2,1); flipud(tssq_err_data(:,2,3))],'g','FaceAlpha',falpha);
f.FaceColor = cols(4,:);
f.EdgeColor = cols(4,:);
loglog(distv,1e-15./abs(distv).^2,'Color','k','Marker','none','LineStyle','--','LineWidth',LW,'MarkerSize',MS);
grid on;
xlabel('Distance to $\Gamma$, $d$','fontsize',FS,'interpreter','latex');
ylabel('Relative error','fontsize',FS,'interpreter','latex');
xlim([min(distv),max(distv)]);
legend('$\epsilon=10^{-6}$, SSQ','','$\epsilon=10^{-6}$, TSSQ','fontsize',FS,'interpreter','latex');
annotation('textarrow',[0.41 0.35],[0.69 0.58],'String','$\mathcal{O}(1/d^2)$','fontsize',FS,'interpreter','latex')
ylim([1e-15,1e8]);

% sanity check
figure('DefaultAxesFontSize',FS);
markerstr = {'^','square','<'};
for ii = 1:numel(tolv)
    loglog(distv,ssq_err_data(:,ii,3),'Color',cols(ii,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
    hold on;
    loglog(distv,tssq_err_data(:,ii,3),'Color',cols(ii,:),'Marker',markerstr(ii),'LineStyle','-','LineWidth',LW,'MarkerSize',MS);
    loglog(distv,adapquad_err_data(:,ii,3),'Color','k','Marker','pentagram','LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor','k');
end
loglog(distv,1e-15./abs(distv).^2,'Color','k','Marker','none','LineStyle','--','LineWidth',LW,'MarkerSize',MS);
grid on;
xlabel('Distance to $\Gamma$, $d$','fontsize',FS,'interpreter','latex');
ylabel('Maximum relative error','fontsize',FS,'interpreter','latex');
xlim([distv(end),distv(1)]);
xticks([1e-8 1e-6 1e-4 1e-2]);
title('sanity check: compare with adaptive quad');

close(1);
close(2);

alignfigs;

if savefig
    disp('saving figures...');
    exportgraphics(figure(3),'figs/long_filament.pdf','Resolution',400);
    exportgraphics(figure(4),'figs/long_filament_err_vs_dist_tol4.pdf','Resolution',400);
    exportgraphics(figure(5),'figs/long_filament_err_vs_dist_tol6.pdf','Resolution',400);
    disp('sucessfully saved figures');
end
