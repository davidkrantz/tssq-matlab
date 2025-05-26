function [errv,errshv,errestv,c,d,pmat,pshmat,irefv] = bases_comp_flat_panel(...
    m,sigma,r,delta,a,bv,n,adj_method,use_bjorck_pereyra,corr_coeff,errest_alt)
% BASES_COMP_FLAT_PANEL compares non-shifted and shifted monomial
%   basis functions to compute
%   
%   I^m(t_0) = \int_{-1}^{1} \frac{(t-a)^r+\delta}{|t-t_0|^m} dt, m=1,3,5
%
% [errv,errshv,c,d,pmat,pshmat,irefv] = bases_comp_flat_panel(...
%   m,sigma,r,delta,a,bv,n,adj_method,use_bjorck_pereyra,corr_coeff,errest_alt) 
%   returns the relative error of the two basis choices as well as the 
%   corresponding basis coefficients and basis integral values.
%
% Without arguments runs test defined below to generate figures in paper
%
% INPUTS:
%   m                  - power of singularity 1/r^m
%   sigma              - function handle for the layer density
%   r                  - vanishing order, (t-a)^r
%   delta              - small offset of integrand at t = a
%   a                  - real part of singularity location t0
%   bv                 - vector of imaginary parts of singularities t0 = a + ib
%   n                  - number of Gauss-Legendre nodes (basis function order)
%   adj_method         - if true, solves non-shifted using adjoint method
%   use_bjorck_pereyra - boolean, use Björck-Pereyra method to solve Vandermonde systems
%   errest_alt         - integer determines how to est cond number of sum
%
% OUTPUTS:
%   errv    - relative absolute error in non-shifted basis evaluation for each b
%   errshv  - relative absolute error in shifted basis evaluation for each b
%   errestv - cancellation error estimate
%   c       - monomial basis coefficients (standard, unshifted)
%   d       - monomial basis coefficients (modified, shifted)
%   pmat    - matrix of integrals P_k^m for standard basis
%   pshmat  - matrix of integrals \widetilde{P}_k^m for modified basis
%   irefv   - reference integral values computed via adaptive quadrature
%
% AUTHOR: David Krantz (davkra@kth.se), May 2025
%
% NOTE: See also matlab/test/test_flat_panel_basis_corr.m for more
%   extensive testing capabilities

if nargin == 0, run_comp_example; return; end

% numerator
h = @(t) t-a; % function small at t=a
f = @(t) sigma(t).*(h(t).^r+delta); % modified "density", RR^T style vanishing near t=a

% nodes
[tj,~] = gauss(n); % double precision

% samples
fj = f(tj); % standard double precision
sj = sigma(tj); % only layer density function, double precision

wbary = bclag_interp_weights(tj); % barycentric Lagrange interpolation wts

sh = a; % shift

% basis coefficients, non-shifted, i.e. c_k coeffs of f
if use_bjorck_pereyra
    c = dvand(tj,fj);
else
    V = ones(n,n); for k=2:n, V(:,k) = V(:,k-1).*(tj-0); end % Vandermonde not transp
    c = V\fj(:);
end

% basis coefficients, shifted, i.e. d_k coeffs of f
if use_bjorck_pereyra
    d = dvand(tj-sh,fj);
else
    Vsh = ones(n,n); for k=2:n, Vsh(:,k) = Vsh(:,k-1).*(tj-sh); end
    d = Vsh\fj(:);
end

% correct first coefficient in shifted basis
if corr_coeff
    sigma_interp = bclag_interp(sj,tj,wbary,sh); % barycentric Lagrange interp of layer dens at t=a
    d(1) = sigma_interp*(h(sh)^r+delta); % eval remaining part of kernel analytically
end

% loop over imag val of root
M = numel(bv);
irefv = zeros(M,1);
errv = zeros(M,1);
errestv = zeros(M,1);
errshv = zeros(M,1);
pmat = zeros(n,M);
pshmat = zeros(n,M);
for i = 1:M
    b = bv(i);
    t0 = a + 1i*b; % root in t

    R = @(t) abs(t-t0); % dist func
    integrand = @(t) f(t)./R(t).^m;
    if abs(a) < 1
        Ie1 = integral(integrand, -1, a, 'reltol', 0, 'abstol', 0);
        Ie2 = integral(integrand, a, 1, 'reltol', 0, 'abstol', 0);
        Ie = Ie1+Ie2;
    else
        Ie = integral(integrand, -1, 1, 'reltol', 0, 'abstol', 0);
    end
    irefv(i) = Ie;

    % basis integrals
    [p1, p3, p5] = rsqrt_pow_integrals(t0,n); % recurrences
    [p1sh, p3sh, p5sh] = rsqrt_pow_integrals_shift(t0,n); % w shift
    if m == 1
        p = p1;
        psh = p1sh;
    elseif m == 3
        p = p3;
        psh = p3sh;
    elseif m == 5
        p = p5;
        psh = p5sh;
    else
        error('invalid value m');
    end
    pmat(:,i) = p;
    pshmat(:,i) = psh;
    
    % non-shifted
    if adj_method
        % adj method
        if use_bjorck_pereyra
            lam = pvand(tj,p);
        else
            lam = V'\p;
        end
        I = sum(lam.*fj);
    else
        % plain non-adj method
        I = sum(c.*p);
    end
    errv(i) = abs((I-Ie)/Ie);

    % plain non-adj method, shifted
    Ish = sum(d.*psh);
    errshv(i) = abs((Ish-Ie)/Ie);
    
    % various estimation methods
    if adj_method
        if errest_alt == 1
            % standard
            kappa = sum(abs(lam.*fj))/abs(I);
            kappa = kappa*n;
        elseif errest_alt == 2
            % condition number of dot product using Euclidean norm
            kappa = norm(lam,inf)*norm(fj,inf)/abs(I);
        end
    else
        if errest_alt == 1
            % standard
            kappa = sum(abs(c.*p))/abs(I);
            kappa = kappa*n;
        elseif errest_alt == 2
            % condition number of dot product using Euclidean norm
            kappa = norm(c,inf)*norm(p,inf)/abs(I);
        end
    end
    % final cancellation error estimate
    errestv(i) = kappa*eps;
end

end

function run_comp_example

format long;
close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

savefig = 1; % saves figures to folder matlab/images

adj_method = 0; % solve non-shifted using adjoint method
corr_coeff = 1; % correct first poly coeff by exact value
use_bjorck_pereyra = 1; % or solve Vandermonde systems using "\"
errest_alt = 2; % which way to estimate cancellation error

a = 0.23; % real of root
bv = logspace(-5,0,20)'; % imag of root, "distance" to panel
delta = 1e-8; % offset
sigma = @(t) sin(t+1.53); % artifical layer density function
r = 2; % vanishing rate
n = 20; % nodes

% Example 1: m=1,3,5 (power of singularity 1/r^m), vary b, fix delta=1e-8
[errv1,errshv1,errestv1,c,d,pmat1,pshmat1,irefv1] = bases_comp_flat_panel(...
    1,sigma,r,delta,a,bv,n,adj_method,use_bjorck_pereyra,corr_coeff,errest_alt);
[errv3,errshv3,errestv3,~,~,pmat3,pshmat3,irefv3] = bases_comp_flat_panel(...
    3,sigma,r,delta,a,bv,n,adj_method,use_bjorck_pereyra,corr_coeff,errest_alt);
[errv5,errshv5,errestv5,~,~,pmat5,pshmat5,irefv5] = bases_comp_flat_panel(...
    5,sigma,r,delta,a,bv,n,adj_method,use_bjorck_pereyra,corr_coeff,errest_alt);

% Example 2: m=5, don't correct coefficient, vary b, fix delta=1e-8
[~,errshv5nocorr,~,~,~,~,~,~] = bases_comp_flat_panel(...
    5,sigma,r,delta,a,bv,n,adj_method,use_bjorck_pereyra,0,errest_alt);

% Example 3: m=5, fix b=1e-4, vary delta
deltav = logspace(-16,0,20);
ndelta = numel(deltav);
[errv5delta,errshv5delta,errestv5delta,irefv5delta] = deal(zeros(ndelta,1));
for i = 1:ndelta
    [errv5delta(i),errshv5delta(i),errestv5delta(i),~,~,~,~,irefv5delta(i)] = bases_comp_flat_panel(...
        5,sigma,r,deltav(i),a,1e-4,n,adj_method,use_bjorck_pereyra,corr_coeff,errest_alt);
end

% prints
fprintf('m=1: max relative error standard basis: %.14e\n',max(errv1));
fprintf('m=1: max relative error modified basis: %.14e\n',max(errshv1));
fprintf('m=3: max relative error standard basis: %.14e\n',max(errv3));
fprintf('m=3: max relative error modified basis: %.14e\n',max(errshv3));
fprintf('m=5: max relative error standard basis: %.14e\n',max(errv5));
fprintf('m=5: max relative error modified basis: %.14e\n',max(errshv5));

% prepare plots
M = numel(bv);
[tj,~] = gauss(n);
h = @(t) t-a;
f = @(t) sigma(t).*(h(t).^r+delta);
fj = f(tj);

% plots
FS = 15; % fontsize
LW = 1; % linewidth
MS = 7; % markersize
tmpcol = colororder; % default Matlab color ordering

% to make figure windows equal in size
figure; 
figure;
figure;

figure('DefaultAxesFontSize',FS);
loglog(abs(bv),errv1+eps,'Color',tmpcol(1,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
hold on;
loglog(abs(bv),errshv1+eps,'Color',tmpcol(1,:),'Marker','square','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
loglog(abs(bv),errv3+eps,'Color',tmpcol(2,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
loglog(abs(bv),errshv3+eps,'Color',tmpcol(2,:),'Marker','square','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
loglog(abs(bv),errv5+eps,'Color',tmpcol(3,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
loglog(abs(bv),errshv5+eps,'Color',tmpcol(3,:),'Marker','square','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
loglog(abs(bv),1e-16*1./abs(bv).^2,'Color','k','Marker','none','LineStyle','--','LineWidth',LW,'MarkerSize',MS);
loglog(abs(bv),errestv1+eps,'Color','k','Marker','.','LineStyle','-','LineWidth',LW,'MarkerSize',MS+1);
loglog(abs(bv),errestv3+eps,'Color','k','Marker','.','LineStyle','-','LineWidth',LW,'MarkerSize',MS+1);
loglog(abs(bv),errestv5+eps,'Color','k','Marker','.','LineStyle','-','LineWidth',LW,'MarkerSize',MS+1);
grid on;
legend('$m=1$, std','$m=1$, mod','$m=3$, std','$m=3$, mod','$m=5$, std','$m=5$, mod','fontsize',FS,'interpreter','latex');
xlabel('Distance to $\Gamma$, $b$','fontsize',FS,'interpreter','latex');
ylabel('Relative error','fontsize',FS,'interpreter','latex');
xticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e-0]);
annotation('textarrow',[0.35 0.5],[0.5 0.55],'String','$\mathcal{O}(1/b^2)~$','fontsize',FS,'interpreter','latex')
ylim([1e-16 1e-6]);
yticks([1e-16 1e-14 1e-12 1e-10 1e-8 1e-6]);

figure('DefaultAxesFontSize',FS);
loglog(nan,nan,'k*','LineStyle','none','MarkerSize',MS);
hold on;
loglog(nan,nan,'k','Marker','square','LineStyle','none','MarkerSize',MS);
loglog(nan,nan,'k','Marker','+','LineStyle','none','MarkerSize',MS);
loglog(abs(bv),vecnorm(c.*pmat1,inf,1),'Color',tmpcol(1,:),'Marker','*','LineStyle','none','LineWidth',LW,'MarkerSize',MS);
loglog(abs(bv),vecnorm(d.*pshmat1,inf,1),'Color',tmpcol(1,:),'Marker','square','LineStyle','none','LineWidth',LW,'MarkerSize',MS+2);
loglog(abs(bv),abs(irefv1),'Color',tmpcol(1,:),'Marker','+','LineStyle','none','LineWidth',LW,'MarkerSize',MS);
loglog(abs(bv),vecnorm(c.*pmat3,inf,1),'Color',tmpcol(2,:),'Marker','*','LineStyle','none','LineWidth',LW,'MarkerSize',MS);
loglog(abs(bv),vecnorm(d.*pshmat3,inf,1),'Color',tmpcol(2,:),'Marker','square','LineStyle','none','LineWidth',LW,'MarkerSize',MS+2);
loglog(abs(bv),abs(irefv3),'Color',tmpcol(2,:),'Marker','+','LineStyle','none','LineWidth',LW,'MarkerSize',MS);
loglog(abs(bv),vecnorm(c.*pmat5,inf,1),'Color',tmpcol(3,:),'Marker','*','LineStyle','none','LineWidth',LW,'MarkerSize',MS);
loglog(abs(bv),vecnorm(d.*pshmat5,inf,1),'Color',tmpcol(3,:),'Marker','square','LineStyle','none','LineWidth',LW,'MarkerSize',MS+2);
loglog(abs(bv),abs(irefv5),'Color',tmpcol(3,:),'Marker','+','LineStyle','none','LineWidth',LW,'MarkerSize',MS);
loglog(abs(bv),1./abs(bv).^2,'Color','k','Marker','none','LineStyle','--','LineWidth',LW,'MarkerSize',MS);
loglog(abs(bv),1./abs(bv).^4,'Color','k','Marker','none','LineStyle','--','LineWidth',LW,'MarkerSize',MS);
grid on;
legend('$\|\mathbf{c}\odot\mathbf{P}_m\|_{\infty}$','$\|\mathbf{d}\odot\widetilde{\mathbf{P}}_m\|_{\infty}$','$|I_m|$','fontsize',FS,'interpreter','latex');
xlabel('Distance to $\Gamma$, $b$','fontsize',FS,'interpreter','latex');
%ylabel('Value','fontsize',FS,'interpreter','latex');
xticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e-0]);
ylim([1e-2,1e20]);
yticks([1e0 1e5 1e10 1e15 1e20]);
annotation('textarrow',[0.2204 0.19],[0.62 0.55],'String','$\mathcal{O}(1/b^2)$','fontsize',FS,'interpreter','latex')
annotation('textarrow',[0.43 0.37],[0.8 0.72],'String','$\mathcal{O}(1/b^4)$','fontsize',FS,'interpreter','latex')

figure('DefaultAxesFontSize',FS);
loglog(abs(bv),errv5+eps,'Color',tmpcol(4,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
hold on;
loglog(abs(bv),errshv5+eps,'Color',tmpcol(5,:),'Marker','square','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
loglog(abs(bv),errshv5nocorr+eps,'Color',tmpcol(6,:),'Marker','o','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
grid on;
legend('$m=5$, std','$m=5$, mod ($d_1$ corrected)','$m=5$, mod ($d_1$ not corrected)','fontsize',FS-1,'interpreter','latex');
xlabel('Distance to $\Gamma$, $b$','fontsize',FS,'interpreter','latex');
ylabel('Relative error','fontsize',FS,'interpreter','latex');
xticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e-0]);
ylim([1e-16 1e-6]);
yticks([1e-16 1e-14 1e-12 1e-10 1e-8 1e-6]);

figure('DefaultAxesFontSize',FS);
loglog(deltav,errv5delta,'Color',tmpcol(4,:),'Marker','*','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
hold on;
loglog(deltav,errshv5delta,'Color',tmpcol(5,:),'Marker','square','LineStyle','-','LineWidth',LW,'MarkerSize',MS);
loglog(deltav,errestv5delta,'Color','k','Marker','.','LineStyle','-','LineWidth',LW,'MarkerSize',MS+1);
grid on;
legend('$m=5$, std','$m=5$, mod','fontsize',FS,'interpreter','latex');
xlabel('$\delta$','fontsize',FS,'interpreter','latex');
ylabel('Relative error','fontsize',FS,'interpreter','latex');
xlim([1e-16 1e0]);
xticks([1e-16 1e-12 1e-8 1e-4 1e0]);
ylim([1e-16 1e-6]);
yticks([1e-16 1e-14 1e-12 1e-10 1e-8 1e-6]);

% extra figures
if ~savefig
    figure;
    plot([-1,1],[0,0],'k');
    hold on;
    plot(tj,zeros(n,1),'k.');
    plot(a*ones(M,1),bv,'x');
    grid on;
    xlim([-1.5,1.5]);
    ylim([-1.5,1.5]);
    legend('E=[-1,1]','GL nodes','tar pts');
    xlabel('Real');
    ylabel('Imag');
    title(['a=' num2str(a)]);
    
    figure;
    tt = linspace(-1,1,1000);
    plot(tt,f(tt));
    hold on;
    plot(tj,fj,'.');
    xline(a,'k-','t=a');
    legend('numerator','samples');
    xlabel('t');
    ylabel('val');
    grid on;
    title(['delta=' num2str(delta)]);
    
    figure;
    semilogy(1:n,abs(c),'.--');
    hold on;
    semilogy(1:n,abs(d),'o--');
    grid on;
    xlim([1,n]);
    xlabel('k');
    ylabel('absolute value of coeff');
    legend('non-shifted: ck','shifted: dk');
end

alignfigs;

if savefig
    disp('saving figures...');
    exportgraphics(figure(4),'figs/ex1_err_vs_dist_std_mod.pdf','Resolution',400);
    exportgraphics(figure(5),'figs/ex1_inf_norm_quad_vec.pdf','Resolution',400);
    exportgraphics(figure(6),'figs/ex2_err_vs_dist_m5_test_corr_coeff.pdf','Resolution',400);
    exportgraphics(figure(7),'figs/ex3_err_vs_delta_m5.pdf','Resolution',400);
    disp('sucessfully saved figures');
end
end
