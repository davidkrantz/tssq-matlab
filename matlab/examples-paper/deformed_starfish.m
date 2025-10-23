% DEFORMED_STARFISH  Example comparing SSQ and TSSQ on a deformed starfish.
%
% This script reproduces results from the paper:
% - Evaluates Stokes slender-body potential
% - Compares SSQ and TSSQ in a Fourier basis
% - Uses adaptive quadrature as reference
%
% AUTHOR: David Krantz (davkra@kth.se)

clear all;
close all;
format long;
rng(123);

opts = default_options('fourier');
savefig = 0;

% Example parameters
Neval = 5e3; % targets per distance
distv = fliplr(logspace(-8,-2,20)).'; % distances to test

% Setup curve and artificial layer density
% adaptive uses t in [0,1)
curve_adap = starfish3D(0.3,5,1); % x(t), y(t), z(t), xp(t), yp(t), zp(t), s(t)
density_adap.f1 = @(t) curve_adap.x(t);
density_adap.f2 = @(t) curve_adap.y(t);
density_adap.f3 = @(t) curve_adap.z(t);
% SSQ and TSSQ uses t in [0,2*pi)
curve = starfish3D_2pi(0.3,5,1);
density.f1 = @(t) curve.x(t)/(2*pi);
density.f2 = @(t) curve.y(t)/(2*pi);
density.f3 = @(t) curve.z(t)/(2*pi);

% Target setup (random normals at fixed distances)
t = rand(1,Neval);
utang = [curve_adap.xp(t); curve_adap.yp(t); curve_adap.zp(t)]./curve_adap.s(t);
v = randn(3,Neval); v = v - utang.*sum(v.*utang,1);
v = v ./ sqrt(sum(v.^2,1));
X = []; Y = []; Z = [];
for d = distv.'
    X = [X curve_adap.x(t) + d*v(1,:)];
    Y = [Y curve_adap.y(t) + d*v(2,:)];
    Z = [Z curve_adap.z(t) + d*v(3,:)];
end
targets = [X; Y; Z];

% Reference solution (adaptive quadrature)
disp('* Reference solution');
opts_ref = opts;
opts_ref.tol = 5e-14;
opts_ref.nquad = 2*opts.nquad+2;
[uref1,uref2,uref3] = adaptive_quadrature_wrapper(curve_adap,density_adap,targets,opts_ref);

% Adaptive quadrature (only for sanity check)
disp(' ')
disp('* Adaptive quadrature')
[adquad1,adquad2,adquad3] = adaptive_quadrature_wrapper(curve_adap,density_adap,targets,opts);

% Standard SSQ
disp(' ')
disp('* SSQ')
[ussq1,ussq2,ussq3,stats_ssq] = ssq_sbt(curve,density,targets,opts);

% TSSQ, correct both I_3 and I_5
disp(' ')
disp('* TSSQ')
[utssq1,utssq2,utssq3,stats_tssq] = tssq_sbt(curve,density,targets,opts);

% TSSQ, correct only I_5
disp(' ')
disp('* TSSQ: correct only I_5')
opts.corrR3 = false;
[utssq1_R5,utssq2_R5,utssq3_R5] = tssq_sbt(curve,density,targets,opts);

% Compute errors
[adapquad_err_data,ssq_err_data,tssq_err_data,tssq_err_R5_data] = deal(zeros(numel(distv),3));
for i = 1:numel(distv)
    idx = (i-1)*Neval+1:i*Neval;
    adquad_errmax = compute_error(uref1(idx), uref2(idx), uref3(idx), adquad1(idx), adquad2(idx), adquad3(idx));
    specquad_errmax = compute_error(uref1(idx), uref2(idx), uref3(idx), ussq1(idx), ussq2(idx), ussq3(idx));
    specquadsh_errmax = compute_error(uref1(idx), uref2(idx), uref3(idx), utssq1(idx), utssq2(idx), utssq3(idx));
    specquadsh_errmax_R5 = compute_error(uref1(idx), uref2(idx), uref3(idx), utssq1_R5(idx), utssq2_R5(idx), utssq3_R5(idx));
    adapquad_err_data(i,1) = min(adquad_errmax);
    ssq_err_data(i,1) = min(specquad_errmax);
    tssq_err_data(i,1) = min(specquadsh_errmax);
    tssq_err_R5_data(i,1) = min(specquadsh_errmax_R5);
    adapquad_err_data(i,2) = mean(adquad_errmax);
    ssq_err_data(i,2) = mean(specquad_errmax);
    tssq_err_data(i,2) = mean(specquadsh_errmax);
    tssq_err_R5_data(i,2) = mean(specquadsh_errmax_R5);
    adapquad_err_data(i,3) = max(adquad_errmax);
    ssq_err_data(i,3) = max(specquad_errmax);
    tssq_err_data(i,3) = max(specquadsh_errmax);
    tssq_err_R5_data(i,3) = max(specquadsh_errmax_R5);
end

% Print stats
stats_ssq
stats_tssq
nbr_pts_per_sec_ssq = stats_ssq.nbr_std_pts/stats_ssq.time_weights
nbr_pts_per_sec_tssq = stats_ssq.nbr_std_pts/stats_tssq.time_weights

% Plots
close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

FS = 15; % fontsize
LW = 1; % linewidth
MS = 8; % markersize
tmpcol = colororder; % default Matlab color ordering
falpha = 0.4; % transparency of area between min/max error
d2 = [distv; flipud(distv)];

figure;
figure;

figure('DefaultAxesFontSize',FS);
tj = adaptive_global_discretization(curve.s,5e-14);
xj = curve.x(tj); yj = curve.y(tj); zj = curve.z(tj);
plot3(xj, yj, zj, '.-k','MarkerSize',6);
hold on;
axis equal tight vis3d;
grid on;
box on;
xlabel('$x$','fontsize',FS,'interpreter','latex');
ylabel('$y$','fontsize',FS,'interpreter','latex');
zlabel('$z$','fontsize',FS,'interpreter','latex');

figure('DefaultAxesFontSize',FS);
markerstr = {'^','square','<'};
loglog(distv,ssq_err_data(:,2),'Color',tmpcol(4,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
hold on;
loglog(distv(1:12),tssq_err_data(1:12,2),'Color',tmpcol(5,:),'Marker',markerstr(2),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',tmpcol(5,:));
loglog(distv(13:end),tssq_err_data(13:end,2),'Color',tmpcol(5,:),'Marker',markerstr(2),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor','none');
loglog(distv(1:12),tssq_err_R5_data(1:12,2),'Color',tmpcol(6,:),'Marker',markerstr(1),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',tmpcol(6,:));
loglog(distv(13:end),tssq_err_R5_data(13:end,2),'Color',tmpcol(6,:),'Marker',markerstr(1),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor','none');
loglog(distv,1e-14./abs(distv).^2,'Color','k','Marker','none','LineStyle','--','LineWidth',LW,'MarkerSize',MS);
f1 = fill(d2,[ssq_err_data(:,1); flipud(ssq_err_data(:,3))],'g','FaceAlpha',falpha);
f1.FaceColor = tmpcol(4,:);
f2 = fill(d2,[tssq_err_data(:,1); flipud(tssq_err_data(:,3))],'g','FaceAlpha',falpha);
f2.FaceColor = tmpcol(5,:);
f3 = fill(d2,[tssq_err_R5_data(:,1); flipud(tssq_err_R5_data(:,3))],'g','FaceAlpha',falpha);
f3.FaceColor = tmpcol(6,:);
legend('mean, SSQ','mean, TSSQ ($I_3$ and $I_5$ corrected)','','mean, TSSQ (only $I_5$ corrected)','fontsize',FS-1,'interpreter','latex');
grid on;
xlabel('Distance to $\Gamma$, $d$','fontsize',FS,'interpreter','latex');
ylabel('Relative error','fontsize',FS,'interpreter','latex');
xlim([distv(end),distv(1)]);
xticks([1e-8 1e-6 1e-4 1e-2]);
annotation('textarrow',[0.32 0.38],[0.57 0.65],'String','$\mathcal{O}(1/d^2)$','fontsize',FS,'interpreter','latex')

figure('DefaultAxesFontSize',FS);
markerstr = {'^','square','<'};
loglog(distv,ssq_err_data(:,2),'Color',tmpcol(4,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
hold on;
loglog(distv,tssq_err_data(:,2),'Color',tmpcol(5,:),'Marker',markerstr(2),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',tmpcol(5,:));
loglog(distv,tssq_err_R5_data(:,2),'Color',tmpcol(6,:),'Marker',markerstr(1),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',tmpcol(6,:));
loglog(distv,adapquad_err_data(:,2),'Color','k','Marker','pentagram','LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor','k');
loglog(distv,1e-14./abs(distv).^2,'Color','k','Marker','none','LineStyle','--','LineWidth',LW,'MarkerSize',MS);
loglog(distv,1e-14./abs(distv).^1,'Color','k','Marker','none','LineStyle','--','LineWidth',LW,'MarkerSize',MS);
legend('std','mod both','mod only R5','fontsize',FS-3,'interpreter','latex');
grid on;
xlabel('Distance to $\Gamma$, $d$','fontsize',FS,'interpreter','latex');
ylabel('Mean relative error','fontsize',FS,'interpreter','latex');
xlim([distv(end),distv(1)]);
xticks([1e-8 1e-6 1e-4 1e-2]);
title('sanity check: compare against adaptive quad');

close(1);
close(2);

alignfigs;

if savefig
    if ~exist('figs', 'dir'); mkdir('figs'); end
    disp('saving figures...');
    exportgraphics(figure(3),'figs/deformed_starfish.pdf','Resolution',400);
    exportgraphics(figure(4),'figs/deformed_starfish_err_vs_dist_corrI3_corrI5.pdf','Resolution',400);
    disp('sucessfully saved figures');
end
