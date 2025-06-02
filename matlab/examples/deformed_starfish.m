% Evaluate Stokes slender body potential from closed-loop deformed thin
% starfish. Compared the standard and modified periodic SSQ methods, with 
% reference values obtained from an adaptive quadrature method.
%
% AUTHOR: David Krantz (davkra@kth.se), June 2025

clear all;
close all;
format long;
rng(123);

% Default setup
slender_eps = 1e-3;
narms = 5;
savefig = 0;

nquad_panel = 16;
Hlim = 4; % Adaptive quadrature distance rule
Hlim_ref = 4;
nquad_ref = 2*nquad_panel+2;

% Case setup
tol = 1e-10;
distv = fliplr(logspace(-8,-2,20)).';  % the dist of all pts from the curve

% Setup squiggle
[x, y, z, xp, yp, zp, ~, ~, ~] = starfish3D(0.3,narms,1);
s = @(t) sqrt(xp(t).^2 + yp(t).^2 + zp(t).^2);

% Artificial density (trivial)
f1 = @(t) x(t);
f2 = @(t) y(t);
f3 = @(t) z(t);

% all targs near curve, in random normal directions a fixed dist away
Ne = 5e3;  % # targs
t = rand(1,Ne);
v = randn(3,Ne); utang = [xp(t);yp(t);zp(t)]./s(t); % sloppy unit tangents
vdotutang = sum(v.*utang,1); v = v - utang.*vdotutang; % orthog v against the tangent
v = v./sqrt(sum(v.*v,1)); % normalize all the v vecs
X = []; Y = []; Z = [];
for i = 1:numel(distv)
    X = [X x(t) + distv(i)*v(1,:)]; % displace by v vecs from pts on curve
    Y = [Y y(t) + distv(i)*v(2,:)];
    Z = [Z z(t) + distv(i)*v(3,:)];
end

% Compute ref solution
disp('* Reference: Adaptive')
% Make sure that we get a different discretization for reference computations,
% otherwise errors get artifically small.
tol_ref = 5e-14;
[tj_ref, wj_ref, npan_ref] = adaptive_panelization(s, nquad_ref, tol_ref);
fprintf('Discretization: nquad=%d, npan=%d\n', nquad_ref, npan_ref);
[uref1, uref2, uref3] = adaptive_quadrature(...
    x(tj_ref), y(tj_ref), z(tj_ref), s(tj_ref), wj_ref, f1(tj_ref), ...
    f2(tj_ref), f3(tj_ref), X, Y, Z, nquad_ref, slender_eps, Hlim_ref);

% Panel-based discretization (only for testing)
fprintf('* Panel-based discretizing, tol=%.1e\n', tol)
[tjadap, wjadap, npan_adap, ~] = adaptive_panelization(s, nquad_panel, tol_ref);
fprintf('nquad=%d, npan=%d\n', nquad_panel, npan_adap);
xjadap = x(tjadap); yjadap = y(tjadap); zjadap = z(tjadap);
sjadap = s(tjadap);
f1jadap = f1(tjadap); f2jadap = f2(tjadap); f3jadap = f3(tjadap);

% Compute adaptive quadrature
disp(' ')
disp('* Adaptive quadrature')
[adquad1, adquad2, adquad3, adstats] = adaptive_quadrature(...
    xjadap, yjadap, zjadap, sjadap, wjadap, f1jadap, f2jadap, f3jadap, ...
    X, Y, Z, nquad_panel, slender_eps, Hlim);

% Global trapezoidal discretization
[d2, y2, z2, xp2, yp2, zp2, ~, ~, ~] = starfish3D_2pi(0.3,narms,1);
s2 = @(t) sqrt(xp2(t).^2 + yp2(t).^2 + zp2(t).^2);
f12 = @(t) d2(t)/(2*pi);
f22 = @(t) y2(t)/(2*pi);
f32 = @(t) z2(t)/(2*pi);
fprintf('* Global discretizing, tol=%.1e\n', 5e-14)
[tj,wj] = adaptive_global_discretization(s2,5e-14);
nquad_global = numel(tj);
fprintf('nquad_global=%d\n', nquad_global);
xj = d2(tj); yj = y2(tj); zj = z2(tj);
sj = s2(tj);
f1j = f12(tj); f2j = f22(tj); f3j = f32(tj);

% Compute special quadrature, correct both I_3 and I_5
disp(' ')
disp('* Interpolatory quadrature, correct I3 and I5')
[specquad1,specquad2,specquad3,specquadsh1,specquadsh2,specquadsh3,cancellation_errest,corr_needed,specstats] = interpolatory_quadrature(...
    xj, yj, zj, sj, tj, wj, f1j, f2j, f3j, X, Y, Z, ...
    d2, y2, z2, s2, nquad_global, slender_eps, tol, true, true);

% correct only I_5
disp(' ')
disp('* Interpolatory quadrature, correct I5')
[~,~,~,specquadsh1_R5,specquadsh2_R5,specquadsh3_R5,~,~,specstats_R5] = interpolatory_quadrature(...
    xj, yj, zj, sj, tj, wj, f1j, f2j, f3j, X, Y, Z, ...
    d2, y2, z2, s2, nquad_global, slender_eps, tol, false, true);

% no corrections (just to get timings)
disp(' ')
disp('* Interpolatory quadrature, no corrections')
[~,~,~,~,~,~,~,~,specstats_std] = interpolatory_quadrature(...
    xj, yj, zj, sj, tj, wj, f1j, f2j, f3j, X, Y, Z, ...
    d2, y2, z2, s2, nquad_global, slender_eps, tol, false, false);

% Compute errors
[adapquad_err_data,specquad_err_data,specquadsh_err_data,specquadsh_err_R5_data] = deal(zeros(numel(distv),3));
[adquad_errmax,specquad_errmax,specquadsh_errmax,specquadsh_errmax_R5] = deal(zeros(1,numel(distv)*Ne));
for i = 1:numel(distv)
    idx = (i-1)*Ne+1:i*Ne;
    adquad_errmax(idx) = compute_error(uref1(idx), uref2(idx), uref3(idx), adquad1(idx), adquad2(idx), adquad3(idx));
    specquad_errmax(idx) = compute_error(uref1(idx), uref2(idx), uref3(idx), specquad1(idx), specquad2(idx), specquad3(idx));
    specquadsh_errmax(idx) = compute_error(uref1(idx), uref2(idx), uref3(idx), specquadsh1(idx), specquadsh2(idx), specquadsh3(idx));
    specquadsh_errmax_R5(idx) = compute_error(uref1(idx), uref2(idx), uref3(idx), specquadsh1_R5(idx), specquadsh2_R5(idx), specquadsh3_R5(idx));
    adapquad_err_data(i,1) = min(adquad_errmax(idx));
    specquad_err_data(i,1) = min(specquad_errmax(idx));
    specquadsh_err_data(i,1) = min(specquadsh_errmax(idx));
    specquadsh_err_R5_data(i,1) = min(specquadsh_errmax_R5(idx));
    adapquad_err_data(i,2) = mean(adquad_errmax(idx));
    specquad_err_data(i,2) = mean(specquad_errmax(idx));
    specquadsh_err_data(i,2) = mean(specquadsh_errmax(idx));
    specquadsh_err_R5_data(i,2) = mean(specquadsh_errmax_R5(idx));
    adapquad_err_data(i,3) = max(adquad_errmax(idx));
    specquad_err_data(i,3) = max(specquad_errmax(idx));
    specquadsh_err_data(i,3) = max(specquadsh_errmax(idx));
    specquadsh_err_R5_data(i,3) = max(specquadsh_errmax_R5(idx));
end

% print errors
epsilon = tol*ones(size(distv));
table(distv,epsilon,adapquad_err_data,specquad_err_data,specquadsh_err_data,specquadsh_err_R5_data)

% plots
close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

FS = 15; % fontsize
LW = 1; % linewidth
MS = 8; % markersize
tmpcol = colororder; % default Matlab color ordering
falpha = 0.4; % transparency of area between min/max error

figure;
figure;

figure('DefaultAxesFontSize',FS);
plot3(xj, yj, zj, '.-k','MarkerSize',6);
hold on;
axis equal tight vis3d;
grid on;
box on;
xlabel('$x$','fontsize',FS,'interpreter','latex');
ylabel('$y$','fontsize',FS,'interpreter','latex');
zlabel('$z$','fontsize',FS,'interpreter','latex');

d2 = [distv; flipud(distv)];
figure('DefaultAxesFontSize',FS);
markerstr = {'^','square','<'};
loglog(distv,specquad_err_data(:,2),'Color',tmpcol(4,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
hold on;
f = fill(d2,[specquad_err_data(:,1); flipud(specquad_err_data(:,3))],'g','FaceAlpha',falpha);
f.FaceColor = tmpcol(4,:);
f.EdgeColor = tmpcol(4,:);
loglog(distv,specquadsh_err_data(:,2),'Color',tmpcol(5,:),'Marker',markerstr(2),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',tmpcol(5,:));
f = fill(d2,[specquadsh_err_data(:,1); flipud(specquadsh_err_data(:,3))],'g','FaceAlpha',falpha);
f.FaceColor = tmpcol(5,:);
f.EdgeColor = tmpcol(5,:);
loglog(distv,specquadsh_err_R5_data(:,2),'Color',tmpcol(6,:),'Marker',markerstr(1),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',tmpcol(6,:));
f = fill(d2,[specquadsh_err_R5_data(:,1); flipud(specquadsh_err_R5_data(:,3))],'g','FaceAlpha',falpha);
f.FaceColor = tmpcol(6,:);
f.EdgeColor = tmpcol(6,:);
loglog(distv,adapquad_err_data(:,2),'Color','k','Marker','pentagram','LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor','k');
loglog(distv,1e-14./abs(distv).^2,'Color','k','Marker','none','LineStyle','--','LineWidth',LW,'MarkerSize',MS);
loglog(distv,1e-14./abs(distv).^1,'Color','k','Marker','none','LineStyle','--','LineWidth',LW,'MarkerSize',MS);
legend('std','mod both','mod only R5','fontsize',FS-3,'interpreter','latex');
grid on;
xlabel('Distance to $\Gamma$, $d$','fontsize',FS,'interpreter','latex');
ylabel('Maximum relative error','fontsize',FS,'interpreter','latex');
xlim([distv(end),distv(1)]);
xticks([1e-8 1e-6 1e-4 1e-2]);

figure('DefaultAxesFontSize',FS);
markerstr = {'^','square','<'};
loglog(distv,specquad_err_data(:,2),'Color',tmpcol(4,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
hold on;
loglog(distv(1:12),specquadsh_err_data(1:12,2),'Color',tmpcol(5,:),'Marker',markerstr(2),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',tmpcol(5,:));
loglog(distv(13:end),specquadsh_err_data(13:end,2),'Color',tmpcol(5,:),'Marker',markerstr(2),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor','none');
loglog(distv(1:12),specquadsh_err_R5_data(1:12,2),'Color',tmpcol(6,:),'Marker',markerstr(1),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',tmpcol(6,:));
loglog(distv(13:end),specquadsh_err_R5_data(13:end,2),'Color',tmpcol(6,:),'Marker',markerstr(1),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor','none');
loglog(distv,1e-14./abs(distv).^2,'Color','k','Marker','none','LineStyle','--','LineWidth',LW,'MarkerSize',MS);
f1 = fill(d2,[specquad_err_data(:,1); flipud(specquad_err_data(:,3))],'g','FaceAlpha',falpha);
f1.FaceColor = tmpcol(4,:);
f2 = fill(d2,[specquadsh_err_data(:,1); flipud(specquadsh_err_data(:,3))],'g','FaceAlpha',falpha);
f2.FaceColor = tmpcol(5,:);
f3 = fill(d2,[specquadsh_err_R5_data(:,1); flipud(specquadsh_err_R5_data(:,3))],'g','FaceAlpha',falpha);
f3.FaceColor = tmpcol(6,:);
legend('mean, std','mean, mod ($I_3$ and $I_5$ corrected)','','mean, mod (only $I_5$ corrected)','fontsize',FS-1,'interpreter','latex');
grid on;
xlabel('Distance to $\Gamma$, $d$','fontsize',FS,'interpreter','latex');
ylabel('Relative error','fontsize',FS,'interpreter','latex');
xlim([distv(end),distv(1)]);
xticks([1e-8 1e-6 1e-4 1e-2]);
annotation('textarrow',[0.32 0.38],[0.57 0.65],'String','$\mathcal{O}(1/d^2)$','fontsize',FS,'interpreter','latex')

alignfigs;

if savefig
    disp('saving figures...');
    exportgraphics(figure(3),'figs/deformed_starfish.pdf','Resolution',400);
    exportgraphics(figure(5),'figs/deformed_starfish_err_vs_dist_corrI3_corrI5.pdf','Resolution',400);
    disp('sucessfully saved figures');
end

function [specquad1,specquad2,specquad3,specquadsh1,specquadsh2,specquadsh3,cancellation_errest,corrneeded,stats] = interpolatory_quadrature(...
    xj, yj, zj, sj, tj, wj, f1j, f2j, f3j, X, Y, Z, x, y, z, s, nquad, slender_eps, tol, correct_R3, correct_R5)

    % Setup
    [specquad1,specquad2,specquad3,specquadsh1,specquadsh2,specquadsh3] = deal(zeros(size(X)));
    cancellation_errest = nan*zeros(size(X));
    corrneeded = false(numel(X),1);
    time_weights = 0;
    kerevals_near = 0;
    maintic = tic();
    time_ker_near = 0;
    time_coeffs = 0;
    time_far = 0;

    % Compute quadrature weights
    atic = tic();
    [all_w1,all_w3,all_w5,all_mu1,all_mu3,all_mu5,specquad_needed,all_roots] = ...
        closed_curve_near_weights(tj, wj, xj, yj, zj, X, Y, Z, tol);

    % Evaluation count
    time_weights = time_weights + toc(atic);
    kerevals_near = kerevals_near + nquad*sum(specquad_needed(:));

    % Evaluate each panel-to-point pair
    for i=1:numel(X)
        Xi = X(i);
        Yi = Y(i);
        Zi = Z(i);
        if specquad_needed(i)
            atic = tic();
            [u1R1,u1R3,u1R5,u2R1,u2R3,u2R5,u3R1,u3R3,u3R5] = deal(zeros(1,nquad));
            for k=1:nquad
                r1k = xj(k)-Xi;
                r2k = yj(k)-Yi;
                r3k = zj(k)-Zi;
                [u1R1(k), u1R3(k), u1R5(k), u2R1(k), u2R3(k), u2R5(k), u3R1(k), u3R3(k), u3R5(k)] ...
                    = slender_body_kernel_split(r1k, r2k, r3k, f1j(k), f2j(k), f3j(k), ...
                                                slender_eps);
            end
            % integrands
            g1R1 = (sj.*u1R1).';
            g1R3 = (sj.*u1R3).';
            g1R5 = (sj.*u1R5).';
            g2R1 = (sj.*u2R1).';
            g2R3 = (sj.*u2R3).';
            g2R5 = (sj.*u2R5).';
            g3R1 = (sj.*u3R1).';
            g3R3 = (sj.*u3R3).';
            g3R5 = (sj.*u3R5).';
            % vectors to sum
            tmp1R1 = all_w1(:,i).*g1R1;
            tmp1R3 = all_w3(:,i).*g1R3;
            tmp1R5 = all_w5(:,i).*g1R5;
            tmp2R1 = all_w1(:,i).*g2R1;
            tmp2R3 = all_w3(:,i).*g2R3;
            tmp2R5 = all_w5(:,i).*g2R5;
            tmp3R1 = all_w1(:,i).*g3R1;
            tmp3R3 = all_w3(:,i).*g3R3;
            tmp3R5 = all_w5(:,i).*g3R5;
            I1R1 = sum(tmp1R1);
            I1R3 = sum(tmp1R3);
            I1R5 = sum(tmp1R5);
            I2R1 = sum(tmp2R1);
            I2R3 = sum(tmp2R3);
            I2R5 = sum(tmp2R5);
            I3R1 = sum(tmp3R1);
            I3R3 = sum(tmp3R3);
            I3R5 = sum(tmp3R5);
            I1R3sh = I1R3; I2R3sh = I2R3; I3R3sh = I3R3;
            I1R5sh = I1R5; I2R5sh = I2R5; I3R5sh = I3R5;
            if correct_R3 || correct_R5
                w3 = all_w3(:,i);
                w5 = all_w5(:,i);
                normw3 = norm(w3,inf);
                normw5 = norm(w5,inf);
                toln = tol;
                % estimate for 1/R^5
                corrR35 = zeros(3,2);
                corrR35(1,2) = cond_sum(normw5,g1R5,I1R5);
                corrR35(2,2) = cond_sum(normw5,g2R5,I2R5);
                corrR35(3,2) = cond_sum(normw5,g3R5,I3R5);
    %             corrR35(1,2) = sum(abs(tmp1R5))/abs(I1R5)*eps;
    %             corrR35(2,2) = sum(abs(tmp2R5))/abs(I2R5)*eps;
    %             corrR35(3,2) = sum(abs(tmp3R5))/abs(I3R5)*eps;
    %             corrR35(1,2) = sum(abs(tmp1R5))/abs(I1R5)*eps*nquad;
    %             corrR35(2,2) = sum(abs(tmp2R5))/abs(I2R5)*eps*nquad;
    %             corrR35(3,2) = sum(abs(tmp3R5))/abs(I3R5)*eps*nquad;
                corrR5 = sum(corrR35(:,2) > toln) > 0;
                if ~corrR5
                    % 1/R^5 ok, now check 1/R^3
                    corrR35(1,1) = cond_sum(normw3,g1R3,I1R3);
                    corrR35(2,1) = cond_sum(normw3,g2R3,I2R3);
                    corrR35(3,1) = cond_sum(normw3,g3R3,I3R3);
    %                 corrR35(1,1) = sum(abs(tmp1R3))/abs(I1R3)*eps;
    %                 corrR35(2,1) = sum(abs(tmp2R3))/abs(I2R3)*eps;
    %                 corrR35(3,1) = sum(abs(tmp3R3))/abs(I3R3)*eps;
    %                 corrR35(1,1) = sum(abs(tmp1R3))/abs(I1R3)*eps*nquad;
    %                 corrR35(2,1) = sum(abs(tmp2R3))/abs(I2R3)*eps*nquad;
    %                 corrR35(3,1) = sum(abs(tmp3R3))/abs(I3R3)*eps*nquad;
                    corrR3 = sum(corrR35(:,1) > toln) > 0;
                else
                    % 1/R^5 bad, assume same corr needed for 1/R^3
                    corrR3 = true;
                    corrR35(:,1) = corrR35(:,2);
                end
                cancellation_errest(i) = max(corrR35,[],'all');
    
                if corrR5 || corrR3
                    corrneeded(i) = true;
                    kstd = get_k_vec(nquad,2*pi).';
                    kmod = get_k_vec(nquad-2,2*pi).';
                    t0 = all_roots(i); % with real part in [0,2*pi)
                    a = real(t0);
                    b = imag(t0);
                    r = exp(-abs(b));
                    tdist = abs(exp(1i*tj)-exp(1i*t0));
                    estd = exp(1i*kstd.'*a);
                    emod = exp(1i*kmod*a);
                    mu1 = all_mu1(:,i);
                    mu3 = all_mu3(:,i);
                    mu5 = all_mu5(:,i);
                    if corrR5
                        p3mod = zeros(numel(kstd),1);
                        p5mod = zeros(numel(kstd),1);
                        p3modk = 0.5*(-(1-r)^2/(2*r*(1+r^2)).*mu3(2:end-1) + ((3/2+kmod-1)./(2*r).*mu1(2:end-1)-(3/2+kmod-2)./(1+r^2).*mu1(1:end-2))./(3/2-1));
                        p5modk = 0.5*(1-r)^(3-5) * (-(1-r)^2/(2*r*(1+r^2)).*mu5(2:end-1) + ((5/2+kmod-1)./(2*r).*mu3(2:end-1)-(5/2+kmod-2)./(1+r^2).*mu3(1:end-2))./(5/2-1));
                        if mod(nquad,2) == 0
                            p3mod(1) = 2*mu3(nquad/2+1)/(1-r)^(3-1);
                            p5mod(1) = 2*mu5(nquad/2+1)/(1-r)^(5-1);
                        else
                            p3mod(1) = 2*mu3((nquad-1)/2+1)/(1-r)^(3-1);
                            p5mod(1) = 2*mu5((nquad-1)/2+1)/(1-r)^(5-1);
                        end
                        p3mod(3:end) = emod.*p3modk;
                        p5mod(3:end) = emod.*p5modk;
                    else
                        p3mod = zeros(numel(kstd),1);
                        p3modk = 0.5*(-(1-r)^2/(2*r*(1+r^2)).*mu3(2:end-1) + ((3/2+kmod-1)./(2*r).*mu1(2:end-1)-(3/2+kmod-2)./(1+r^2).*mu1(1:end-2))./(3/2-1));
                        if mod(nquad,2) == 0
                            p3mod(1) = 2*mu3(nquad/2+1)/(1-r)^(3-1);
                        else
                            p3mod(1) = 2*mu3((nquad-1)/2+1)/(1-r)^(3-1);
                        end
                        p3mod(3:end) = emod.*p3modk;
                    end
    
                    % store all relevant values 
                    GR35 = zeros(3*nquad,2);
                    GR35(:,1) = [g1R3;g2R3;g3R3];
                    GR35(:,2) = [g1R5;g2R5;g3R5];
                    IR35 = zeros(3,2);
                    IR35(:,1) = [I1R3;I2R3;I3R3];
                    IR35(:,2) = [I1R5;I2R5;I3R5];
        
                    % eval layer dens at a=real(t0)
                    f1coeff = fftshift(fft(f1j)).'/nquad;
                    f2coeff = fftshift(fft(f2j)).'/nquad;
                    f3coeff = fftshift(fft(f3j)).'/nquad;
                    f1a = real(estd*f1coeff);
                    f2a = real(estd*f2coeff);
                    f3a = real(estd*f3coeff);
                    %f1a = x(a)/(2*pi); f2a = y(a)/(2*pi); f3a = z(a)/(2*pi); % analytic evaluation
                    % eval remaining integrand analytically
                    xa = x(a); ya = y(a); za = z(a);
                    r1a = xa-Xi; r2a = ya-Yi; r3a = za-Zi;
                    tdista = abs(exp(1i*a)-exp(1i*t0));
                    sa = s(a);
                    [~, u1R3a, u1R5a, ~, u2R3a, u2R5a, ~, u3R3a, u3R5a] ...
                        = slender_body_kernel_split(r1a, r2a, r3a, f1a, f2a, f3a, slender_eps);
                    uR35a = [u1R3a, u1R5a;
                             u2R3a, u2R5a;
                             u3R3a, u3R5a];
       
                    tdistMat = zeros(nquad,2);
                    tmp = tdist.^(2*1+1);
                    tdistMat(:,1) = tmp;
                    tdistMat(:,2) = tmp.*tdist.*tdist;
                    for ii = 1:3
                        for jj = 1:2
                            corr = corrR35(ii,jj) > toln;
                            if ~isnan(corr) && corr
                                g = GR35(((ii-1)*nquad+1):ii*nquad,jj);
                                h = g.*tdistMat(:,jj); % to expand in modified Fourier basis
                                ccoeff = fftshift(fft(h))/nquad;
                                btic = tic();
                                [a0,a1,bcoeff] = fourier2modcoeffs(ccoeff,a); % map c_k --> (a_0,a_1,b_k)
                                time_coeffs =  time_coeffs + toc(btic);
                                dcoeff = [a0;a1;bcoeff];
                                d1coeff = sa*uR35a(ii,jj)*tdista^(2*jj+1); % correction to d(1)
                                dcoeff(1) = d1coeff;
                                if jj == 1 && correct_R3
                                    IR35(ii,jj) = real(sum(dcoeff.*p3mod));
                                elseif jj == 2 && correct_R5
                                    IR35(ii,jj) = real(sum(dcoeff.*p5mod));
                                end
                            end
                        end
                    end
                    I1R3sh = IR35(1,1);
                    I2R3sh = IR35(2,1);
                    I3R3sh = IR35(3,1);
                    I1R5sh = IR35(1,2);
                    I2R5sh = IR35(2,2);
                    I3R5sh = IR35(3,2);
                end
            end

            q1 = I1R1 + I1R3 + I1R5;
            q2 = I2R1 + I2R3 + I2R5;
            q3 = I3R1 + I3R3 + I3R5;
            q1sh = I1R1 + I1R3sh + I1R5sh;
            q2sh = I2R1 + I2R3sh + I2R5sh;
            q3sh = I3R1 + I3R3sh + I3R5sh;
            time_ker_near =  time_ker_near + toc(atic);
        else            
            atic = tic();
            [q1, q2, q3] = quadsum(xj, yj, zj, sj, wj, f1j, f2j, f3j, ...
                                   Xi, Yi, Zi, nquad, slender_eps);
            q1sh = q1; q2sh = q2; q3sh = q3;
            time_far = time_far + toc(atic);
        end
        specquad1(i) = specquad1(i) + q1;
        specquad2(i) = specquad2(i) + q2;
        specquad3(i) = specquad3(i) + q3;
        specquadsh1(i) = specquadsh1(i) + q1sh;
        specquadsh2(i) = specquadsh2(i) + q2sh;        
        specquadsh3(i) = specquadsh3(i) + q3sh;
    end

    toc(maintic)
    fprintf('Near field kernel evals: %e\n', kerevals_near);
    fprintf('Near field kernel evals time: %f\n', time_ker_near)
    fprintf('Total time line3_near_weights: %f\n', time_weights)
    fprintf('Far field time: %f\n', time_far)
    fprintf('Number of target points: %d\n', numel(X))
    fprintf('Number special quadrature: %d\n', sum(specquad_needed))
    fprintf('Number corrections: %d\n', sum(corrneeded))
    stats.kerevals_near = kerevals_near;
    stats.time_weights = time_weights;
    stats.time_ker_near = time_ker_near;
    stats.time_coeffs = time_coeffs;
end
