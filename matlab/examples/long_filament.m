% Evaluate Stokes slender body potential from closed-loop squiggle curve.
% Compared the standard and modified SSQ (non-periodic) methods, with 
% reference values obtained from an adaptive quadrature method.
%
% AUTHOR: David Krantz (davkra@kth.se), May 2025
%
% NOTE: Based on demo_long_fiber.m from the GitHub repository
%   https://github.com/ludvigak/linequad

clear all;
close all;
format long;
rng(123);

% Default setup
slender_eps = 1e-3;
savefig = 0;
use_bjorck_pereyra = true;
use_mod = true;

UPSAMPLE = true; % Upsample before applying interpolatory quadrature
nquad = 16;
rho = 3; % Interpolatory quadrature limit rule, rho=1 --> disable specquad
Hlim = 1; % Adaptive quadrature distance rule 

% Case setup
tolv = [1e-6;1e-8;1e-10];
distv = fliplr(logspace(-8,-2,20)).';  % the dist of all pts from the curve

% Setup fiber
[x, y, z, xp, yp, zp] = squiggle();
s = @(t) sqrt(xp(t).^2 + yp(t).^2 + zp(t).^2);

% Artificial density (trivial)
f1 = @(t) x(t);
f2 = @(t) y(t);
f3 = @(t) z(t);

% all targs near curve, in random normal directions a fixed dist away
Ne = 5e1;  % # targs
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

% Compute ref solution using adaptive quadrature
% Make sure that we get a different discretization for reference 
% computations, otherwise errors get artifically small
disp('* Reference: Adaptive')
nquad_ref = nquad+2;
tol_ref = 5e-14;
Hlim_ref = 1;
[tj_ref, wj_ref, npan_ref] = adaptive_panelization(s, nquad_ref, tol_ref);
fprintf('Discretization: nquad=%d, npan=%d\n', nquad_ref, npan_ref);
[uref1, uref2, uref3] = adaptive_quadrature( ...
    x(tj_ref), y(tj_ref), z(tj_ref), s(tj_ref), wj_ref,f1(tj_ref), ...
    f2(tj_ref), f3(tj_ref), X, Y, Z, nquad_ref, slender_eps, Hlim_ref);

[adapquad_errmaxmax,specquad_errmaxmax,specquadsh_errmaxmax] = deal(zeros(numel(distv),numel(tolv)));
stats = cell(numel(tolv),1);
for ii = 1:numel(tolv)
    tol = tolv(ii);

    % Discretize
    fprintf('* Discretizing, tol=%.1e\n', tol)
    [tj, wj, npan, edges] = adaptive_panelization(s, nquad, tol);
    fprintf('nquad=%d, npan=%d\n', nquad, npan);
    xj = x(tj); yj = y(tj); zj = z(tj);
    sj = s(tj);
    f1j = f1(tj); f2j = f2(tj); f3j = f3(tj);

    % Compute adaptive quadrature (just for sanity check)
    disp(' ')
    disp('* Adaptive quadrature')
    [adquad1, adquad2, adquad3, adstats] = adaptive_quadrature(...
        xj, yj, zj, sj, wj, f1j, f2j, f3j, X, Y, Z, nquad, slender_eps, Hlim);
    
    % Compute special quadrature
    disp(' ')
    disp('* Interpolatory quadrature')
    [specquad1,specquad2,specquad3,specquadsh1,specquadsh2,specquadsh3,cancellation_errest,corr_needed,specstats] = interpolatory_quadrature(...
        xj, yj, zj, sj, tj, wj, f1j, f2j, f3j, X, Y, Z, ...
        x, y, z, s, nquad, edges, rho, UPSAMPLE, slender_eps, ...
        tol, use_bjorck_pereyra, use_mod);
    stats{ii} = specstats;
    
    % Compute errors
    [adquad_errmax,specquad_errmax,specquadsh_errmax] = deal(zeros(1,numel(distv)*Ne));
    for i = 1:numel(distv)
        idx = (i-1)*Ne+1:i*Ne;
        adquad_errmax(idx) = compute_error(uref1(idx), uref2(idx), uref3(idx), adquad1(idx), adquad2(idx), adquad3(idx));
        specquad_errmax(idx) = compute_error(uref1(idx), uref2(idx), uref3(idx), specquad1(idx), specquad2(idx), specquad3(idx));
        specquadsh_errmax(idx) = compute_error(uref1(idx), uref2(idx), uref3(idx), specquadsh1(idx), specquadsh2(idx), specquadsh3(idx));
        adapquad_errmaxmax(i,ii) = max(adquad_errmax(idx));
        specquad_errmaxmax(i,ii) = max(specquad_errmax(idx));
        specquadsh_errmaxmax(i,ii) = max(specquadsh_errmax(idx));
    end
end

% print errors
for ii = 1:numel(tolv)
    epsilon = tolv(ii)*ones(size(distv));
    table(distv,epsilon,adapquad_errmaxmax(:,ii),specquad_errmaxmax(:,ii),specquadsh_errmaxmax(:,ii))
end

% print timings
[kerevals_near,time_std_weights,time_ker_near,time_mod_weights,time_cancel_est,time_mod_vander_solve,time_mod_sigma_interp_a,time_mod_kernel_eval_a,time_mod_basis_int] = deal(zeros(numel(tolv),1));
for ii = 1:numel(tolv)
    kerevals_near(ii) = stats{ii}.kerevals_near;
    time_std_weights(ii) = stats{ii}.time_std_weights;
    time_ker_near(ii) = stats{ii}.time_ker_near;
    time_mod_weights(ii) = stats{ii}.time_mod_weights;
    time_cancel_est(ii) = stats{ii}.time_cancel_est;
    time_mod_basis_int(ii) = stats{ii}.time_mod_basis_int;
    time_mod_vander_solve(ii) = stats{ii}.time_mod_vander_solve;
    time_mod_sigma_interp_a(ii) = stats{ii}.time_mod_sigma_interp_a;
    time_mod_kernel_eval_a(ii) = stats{ii}.time_mod_kernel_eval_a;
end
table(tolv,kerevals_near,time_std_weights,time_mod_weights,time_cancel_est,time_ker_near,time_mod_basis_int,time_mod_vander_solve,time_mod_sigma_interp_a,time_mod_kernel_eval_a)
time_total = time_std_weights+time_ker_near

% plots
close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% plots
FS = 15; % fontsize
LW = 1; % linewidth
MS = 8; % markersize
tmpcol = colororder; % default Matlab color ordering

figure;
figure

[tj, wj, npan, edges] = adaptive_panelization(s, nquad, 1e-6);
fprintf('nquad=%d, npan=%d\n', nquad, npan);
xj = x(tj); yj = y(tj); zj = z(tj);

figure('DefaultAxesFontSize',FS);
plot3(xj, yj, zj, '.-k','MarkerSize',6);
hold on;
%plot3(X(1:Ne),Y(1:Ne),Z(1:Ne),'r.','MarkerSize',6);
axis equal tight vis3d;
grid on;
box on;
xlabel('$x$','fontsize',FS,'interpreter','latex');
ylabel('$y$','fontsize',FS,'interpreter','latex');
zlabel('$z$','fontsize',FS,'interpreter','latex');

figure('DefaultAxesFontSize',FS);
markerstr = {'^','square','<'};
for ii = 1:numel(tolv)
    loglog(distv,specquad_errmaxmax(:,ii),'Color',tmpcol(ii,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
    hold on;
    loglog(distv,specquadsh_errmaxmax(:,ii),'Color',tmpcol(ii,:),'Marker',markerstr(ii),'LineStyle','-','LineWidth',LW,'MarkerSize',MS);
end
loglog(distv,adapquad_errmaxmax(:,ii),'Color','k','Marker','pentagram','LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor','k');
loglog(distv,1e-15./abs(distv).^2,'Color','k','Marker','none','LineStyle','--','LineWidth',LW,'MarkerSize',MS);
legend('$\epsilon=10^{-6}$, std','$\epsilon=10^{-6}$, mod','$\epsilon=10^{-8}$, std','$\epsilon=10^{-8}$, mod','$\epsilon=10^{-10}$, std','$\epsilon=10^{-10}$, mod','','fontsize',FS-3,'interpreter','latex');
grid on;
xlabel('Distance to $\Gamma$, $d$','fontsize',FS,'interpreter','latex');
ylabel('Maximum relative error','fontsize',FS,'interpreter','latex');
xlim([distv(end),distv(1)]);
xticks([1e-8 1e-6 1e-4 1e-2]);
annotation('textarrow',[0.43 0.35],[0.78 0.65],'String','$\mathcal{O}(1/d^2)$','fontsize',FS,'interpreter','latex')

figure('DefaultAxesFontSize',FS);
markerstr = {'^','square','<'};
loglog(distv,specquad_errmaxmax(:,1),'Color',tmpcol(1,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
hold on;
loglog(distv,specquadsh_errmaxmax(:,1),'Color',tmpcol(1,:),'Marker',markerstr(1),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',tmpcol(1,:));
loglog(distv,specquad_errmaxmax(:,2),'Color',tmpcol(2,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
loglog(distv(1:14),specquadsh_errmaxmax(1:14,2),'Color',tmpcol(2,:),'Marker',markerstr(2),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',tmpcol(2,:));
loglog(distv(15:end),specquadsh_errmaxmax(15:end,2),'Color',tmpcol(2,:),'Marker',markerstr(2),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor','none');
loglog(distv,specquad_errmaxmax(:,3),'Color',tmpcol(3,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
loglog(distv(1:8),specquadsh_errmaxmax(1:8,3),'Color',tmpcol(3,:),'Marker',markerstr(3),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',tmpcol(3,:));
loglog(distv(9:end),specquadsh_errmaxmax(9:end,3),'Color',tmpcol(3,:),'Marker',markerstr(3),'LineStyle','-','LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor','none');
loglog(distv,1e-15./abs(distv).^2,'Color','k','Marker','none','LineStyle','--','LineWidth',LW,'MarkerSize',MS);
legend('$\epsilon=10^{-6}$, std','$\epsilon=10^{-6}$, mod','$\epsilon=10^{-8}$, std','$\epsilon=10^{-8}$, mod','','$\epsilon=10^{-10}$, std','$\epsilon=10^{-10}$, mod','','fontsize',FS-3,'interpreter','latex');
grid on;
xlabel('Distance to $\Gamma$, $d$','fontsize',FS,'interpreter','latex');
ylabel('Maximum relative error','fontsize',FS,'interpreter','latex');
xlim([distv(end),distv(1)]);
xticks([1e-8 1e-6 1e-4 1e-2]);
annotation('textarrow',[0.43 0.35],[0.78 0.65],'String','$\mathcal{O}(1/d^2)$','fontsize',FS,'interpreter','latex')

alignfigs;

if savefig
    disp('saving figures...');
    exportgraphics(figure(3),'figs/long_filament.pdf','Resolution',400);
    exportgraphics(figure(5),'figs/long_filament_err_vs_dist.pdf','Resolution',400);
    disp('sucessfully saved figures');
end

function [specquad1,specquad2,specquad3,specquadsh1,specquadsh2,specquadsh3,cancellation_errest,corrneeded,stats] = interpolatory_quadrature(...
    xj, yj, zj, sj, tj, wj, f1j, f2j, f3j, X, Y, Z, x, y, z, s, nquad, edges, rho, UPSAMPLE, slender_eps, tol, use_bjorck_pereyra, use_mod)
    [specquad1,specquad2,specquad3,specquadsh1,specquadsh2,specquadsh3] = deal(zeros(size(X)));
    cancellation_errest = nan*zeros(size(X));
    corrneeded = false(numel(X),1);
    npan = numel(xj)/nquad;    
    [tgl, wgl] = legendre.gauss(nquad);
    if UPSAMPLE
        nquad2 = 2*nquad;
        [tgl2, wgl2] = legendre.gauss(nquad2);
        B = bclag_interp_matrix(tgl, tgl2);
    else
        nquad2 = nquad;
        B = 1;
        tgl2 = tgl;
        wgl2 = wgl;
    end
    wbary = bclag_interp_weights(tgl2);
    time_std_weights = 0;
    time_mod_weights = 0;
    time_mod_basis_int = 0;
    time_mod_vander_solve = 0;
    time_mod_sigma_interp_a = 0;
    time_mod_kernel_eval_a = 0;
    time_cancel_est = 0;
    kerevals_near = 0;
    maintic = tic();
    time_ker_near = 0;
    time_far = 0;
    for j=1:npan
        % Load panel
        idx = (1:nquad) + nquad*(j-1);
        xjpan = xj(idx);
        yjpan = yj(idx);
        zjpan = zj(idx);
        sjpan = sj(idx);
        tjpan = tj(idx);
        wjpan = wj(idx);
        f1jpan = f1j(idx);
        f2jpan = f2j(idx);
        f3jpan = f3j(idx);
        % Upsample panel
        xjpan_up  = xjpan *B';
        yjpan_up  = yjpan *B';
        zjpan_up  = zjpan *B';
        sjpan_up  = sjpan *B';
        f1jpan_up = f1jpan*B';
        f2jpan_up = f2jpan*B';
        f3jpan_up = f3jpan*B';
        tjpan_up = tjpan.'*B';
        % current panel endpoints
        ta = edges(j);
        tb = edges(j+1);
        tsc = (tb-ta)/2;
        % Compute quadrature weights
        atic = tic();
        [all_w1, all_w3, all_w5, specquad_needed, all_roots] = line3_near_weights(tgl2, wgl2, xjpan_up, yjpan_up, zjpan_up, ...
                                                          X, Y, Z, rho);
        t0_all_roots = (ta+tb)/2+(tb-ta)/2*all_roots; % roots in [0,1]
        % Evaluation count
        time_std_weights = time_std_weights + toc(atic);
        kerevals_near = kerevals_near + nquad2*sum(specquad_needed(:));
        % Evaluate each panel-to-point pair
        for i=1:numel(X)    
            Xi = X(i);
            Yi = Y(i);
            Zi = Z(i);
            if specquad_needed(i)
                atic = tic();
                [u1R1,u1R3,u1R5,u2R1,u2R3,u2R5,u3R1,u3R3,u3R5] = deal(zeros(1,nquad2));
                for k=1:nquad2
                    r1k = xjpan_up(k)-Xi;
                    r2k = yjpan_up(k)-Yi;
                    r3k = zjpan_up(k)-Zi;
                    [u1R1(k), u1R3(k), u1R5(k), u2R1(k), u2R3(k), u2R5(k), u3R1(k), u3R3(k), u3R5(k)] ...
                        = slender_body_kernel_split(r1k, r2k, r3k, f1jpan_up(k), f2jpan_up(k), f3jpan_up(k), ...
                                                    slender_eps);
                end
                % integrands
                g1R1 = (sjpan_up.*u1R1).';
                g1R3 = (sjpan_up.*u1R3).';
                g1R5 = (sjpan_up.*u1R5).';
                g2R1 = (sjpan_up.*u2R1).';
                g2R3 = (sjpan_up.*u2R3).';
                g2R5 = (sjpan_up.*u2R5).';
                g3R1 = (sjpan_up.*u3R1).';
                g3R3 = (sjpan_up.*u3R3).';
                g3R5 = (sjpan_up.*u3R5).';
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
                scale_fac = sum(wjpan)/2;
                v0 = all_roots(i); % root having real part in [-1,1]
                if abs(real(v0)) <= 1.1 && use_mod
                    btic = tic();
                    w3 = all_w3(:,i);
                    w5 = all_w5(:,i);
                    normw3 = norm(w3,inf);
                    normw5 = norm(w5,inf);
                    % estimate for 1/R^5
                    errestR35 = zeros(3,2);
                    errestR35(1,2) = cond_sum(normw5,g1R5,I1R5);
                    errestR35(2,2) = cond_sum(normw5,g2R5,I2R5);
                    errestR35(3,2) = cond_sum(normw5,g3R5,I3R5);
%                     errestR35(1,2) = cond_sum(normw5,g1R5,I1R5)*scale_fac;
%                     errestR35(2,2) = cond_sum(normw5,g2R5,I2R5)*scale_fac;
%                     errestR35(3,2) = cond_sum(normw5,g3R5,I3R5)*scale_fac;
%                     errestR35(1,2) = sum(abs(tmp1R5))/abs(I1R5)*eps*nquad2;
%                     errestR35(2,2) = sum(abs(tmp2R5))/abs(I1R5)*eps*nquad2;
%                     errestR35(3,2) = sum(abs(tmp3R5))/abs(I1R5)*eps*nquad2;
                    corrR5 = sum(errestR35(:,2) > tol) > 0;
                    if ~corrR5
                        % 1/R^5 ok, now check 1/R^3
                        errestR35(1,1) = cond_sum(normw3,g1R3,I1R3);
                        errestR35(2,1) = cond_sum(normw3,g2R3,I2R3);
                        errestR35(3,1) = cond_sum(normw3,g3R3,I3R3);
%                         errestR35(1,1) = cond_sum(normw3,g1R3,I1R3)*scale_fac;
%                         errestR35(2,1) = cond_sum(normw3,g2R3,I2R3)*scale_fac;
%                         errestR35(3,1) = cond_sum(normw3,g3R3,I3R3)*scale_fac;
%                         errestR35(1,1) = sum(abs(tmp1R3))/abs(I1R3)*eps*nquad2;
%                         errestR35(2,1) = sum(abs(tmp2R3))/abs(I1R3)*eps*nquad2;
%                         errestR35(3,1) = sum(abs(tmp3R3))/abs(I1R3)*eps*nquad2;
                        corrR3 = sum(errestR35(:,1) > tol) > 0;
                    else
                        % 1/R^5 bad, assume same corr needed for 1/R^3
                        corrR3 = true;
                        errestR35(:,1) = errestR35(:,2);
                    end
                    time_cancel_est = time_cancel_est + toc(btic);
                    cancellation_errest(i) = max([max(errestR35,[],'all'),cancellation_errest(i)]);
    
                    if corrR5 || corrR3
                        ctic = tic();
                        corrneeded(i) = true;
                        sh = real(v0);
                        dtic = tic();
                        if corrR5
                            [~, p3sh, p5sh] = rsqrt_pow_integrals_shift(v0,nquad2);
                        else
                            [~, p3sh] = rsqrt_pow_integrals_shift(v0,nquad2);
                        end
                        time_mod_basis_int = time_mod_basis_int + toc(dtic);
    
                        % root w real part [0,1] (actual root of gamma(t), t in [0,1])
                        t0 = t0_all_roots(i);
                        a = real(t0);
                        tdist = abs(tjpan_up.'-t0);
        
                        % store all relevant values 
                        GR35 = zeros(3*nquad2,2);
                        GR35(:,1) = [g1R3;g2R3;g3R3];
                        GR35(:,2) = [g1R5;g2R5;g3R5];
                        IR35 = zeros(3,2);
                        IR35(:,1) = [I1R3;I2R3;I3R3];
                        IR35(:,2) = [I1R5;I2R5;I3R5];
        
                        % eval layer dens at real(v0) in [-1,1]
                        etic = tic();
                        [f1pana,f2pana,f3pana] = bclag_interp_sigma(f1jpan_up,f2jpan_up,f3jpan_up,tgl2,wbary,sh);
                        time_mod_sigma_interp_a = time_mod_sigma_interp_a + toc(etic);
                        %f1pana = x(a); f2pana = y(a); f3pana = z(a);
                        % eval remaining integrand analytically
                        ftic = tic();
                        xa = x(a); ya = y(a); za = z(a);
                        r1a = xa-Xi; r2a = ya-Yi; r3a = za-Zi;
                        tdista = abs(a-t0);
                        sa = s(a);
                        [~, u1R3a, u1R5a, ~, u2R3a, u2R5a, ~, u3R3a, u3R5a] ...
                            = slender_body_kernel_split(r1a, r2a, r3a, f1pana, f2pana, f3pana, slender_eps);
                        uR35a = [u1R3a, u1R5a;
                                 u2R3a, u2R5a;
                                 u3R3a, u3R5a];
                        time_mod_kernel_eval_a = time_mod_kernel_eval_a + toc(ftic);

                        % precompute factorization of Vandermonde matrix
                        alpha = tgl2-sh;
                        if ~use_bjorck_pereyra
                            gtic = tic();
                            V = ones(nquad2,nquad2); for k=2:nquad2, V(:,k) = V(:,k-1).*alpha; end
                            [L,U,P] = lu(V);
                            %[Q,R] = qr(V);
                            time_mod_vander_solve = time_mod_vander_solve + toc(gtic);
                        end

                        time_mod_weights = time_mod_weights + toc(ctic);

                        tdistMat = zeros(nquad2,2);
                        tmp = tdist.^(2*1+1);
                        tdistMat(:,1) = tmp;
                        tdistMat(:,2) = tmp.*tdist.*tdist;
                        for ii = 1:3
                            for jj = 1:2
                                corr = errestR35(ii,jj) > tol;
                                if corr
                                    dtic = tic();
                                    g = GR35(((ii-1)*nquad2+1):ii*nquad2,jj);
                                    h = g.*tdistMat(:,jj); % to expand in shifted monomial basis
                                    htic = tic();
                                    if use_bjorck_pereyra
                                        dcoeff = dvand(alpha,h); % Björck-Pereyra
                                    else
                                        dcoeff = U\(L\(P*h));
                                        %dcoeff = R\(Q\h);
                                    end
                                    %warning('off','MATLAB:nearlySingularMatrix');
                                    %dcoeff = V\h;
                                    %warning('on','MATLAB:nearlySingularMatrix')
                                    time_mod_vander_solve = time_mod_vander_solve + toc(htic);
                                    d1coeff = sa*uR35a(ii,jj)*tdista^(2*jj+1); % correction to d(1)
                                    dcoeff(1) = d1coeff;
                                    time_mod_weights = time_mod_weights + toc(dtic);
                                    if jj == 1
                                        IR35(ii,jj) = sum(dcoeff.*p3sh)/tsc^3;
                                    else
                                        IR35(ii,jj) = sum(dcoeff.*p5sh)/tsc^5;
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
                % Rescale (weights are for [-1,1])
                q1 = q1*scale_fac;
                q2 = q2*scale_fac;
                q3 = q3*scale_fac;
                q1sh = q1sh*scale_fac;
                q2sh = q2sh*scale_fac;
                q3sh = q3sh*scale_fac;
                time_ker_near =  time_ker_near + toc(atic);
            else            
                atic = tic();
                [q1, q2, q3] = quadsum(xjpan, yjpan, zjpan, sjpan, wjpan, f1jpan, f2jpan, f3jpan, ...
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
    end
    toc(maintic)
    fprintf('Near field kernel evals: %e\n', kerevals_near);
    fprintf('Near field kernel evals time: %f\n', time_ker_near)
    fprintf('Total time line3_near_weights: %f\n', time_std_weights)
    fprintf('Far field time: %f\n', time_far)
    stats.kerevals_near = kerevals_near;
    stats.time_std_weights = time_std_weights;
    stats.time_ker_near = time_ker_near;
    stats.time_mod_weights = time_mod_weights;
    stats.time_cancel_est = time_cancel_est;
    stats.time_mod_basis_int = time_mod_basis_int;
    stats.time_mod_sigma_interp_a = time_mod_sigma_interp_a;
    stats.time_mod_kernel_eval_a = time_mod_kernel_eval_a;
    stats.time_mod_vander_solve = time_mod_vander_solve;
end
