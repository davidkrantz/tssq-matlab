function [errv,errshv,c,d,p3mat,p3shmat,irefv] = test_flat_panel_basis_corr(...
    sigma,r,delta,a,bv,n,use_vpa,corr_coeff_exact,corr_coeff_interp_sig,...
    solve_shifted,use_bjorck_pereyra)
% TEST_FLAT_PANEL_BASIS_CORR compares non-shifted and shifted monomial
%   basis functions to compute
%   
%   I^3(t_0) = \int_{-1}^{1} \frac{f(t)}{|t-t_0|^3} dt
%
% [errv,errshv,c,d,p3mat,p3shmat,irefv] = test_flat_panel_basis_corr(...
%   sigma,r,delta,a,bv,n,use_vpa,corr_coeff_exact,corr_coeff_interp_sig,...
%   solve_shifted,use_bjorck_pereyra) returns the relative error of the two
%   basis choices as well as the corresponding basis coefficients and basis
%   integral values.
%
% Without arguments runs "test 1" with both corrections active
%
% INPUTS:
%   sigma                - function handle for the layer density
%   r                    - vanishing order of h(t) at t = a
%   delta                - small offset for location of vanishing, h(t) = t - a + delta
%   a                    - real part of singularity location t0
%   bv                   - vector of imaginary parts of singularities t0 = a + ib
%   n                    - number of Gauss-Legendre nodes (basis function order)
%   use_vpa              - whether to use vpa to compute basis integrals Pkm and Qkm
%   corr_coeff_exact     - if true, correct first basis coeff in shifted case with exact f(a)
%   corr_coeff_interp_sig- if true, correct first coeff using interpolated sigma(a)
%   solve_shifted        - string, 'dp'/'mp'/'qp' for double/mixed/quad precision solve
%   use_bjorck_pereyra   - boolean, use Bjorck-Pereyra method to solve Vandermonde
%
% OUTPUTS:
%   errv      - relative error in non-shifted basis evaluation for each b
%   errshv    - relative error in shifted basis evaluation for each b
%   c         - monomial basis coefficients (non-shifted)
%   d         - monomial basis coefficients (shifted)
%   p3mat     - matrix of integrals P_k^3 for non-shifted basis
%   p3shmat   - matrix of integrals Q_k^3 for shifted basis
%   irefv     - reference integral values computed via adaptive quadrature
%
% AUTHOR: David Krantz (davkra@kth.se) 2025-04-10

if nargin==0, test_corrections; return; end

% numerator
h = @(t) t-a+delta; % function small at t=a
f = @(t) sigma(t).*h(t).^r; % modified "density", RR^T style vanishing near t=a

% nodes
[tj,~] = gauss(n); % double precision
[tjqp,~] = gauss_qp(n); % quadruple precision

% samples
fj = f(tj); % standard double precision
fjqp = f(tjqp); % quadruple precision
fjmp = sigma(tj).*h(tjqp).^r; % "mixed" precision: sigma in dp and rest in qp
sj = sigma(tj); % only layer density function, double precision

wbary = bclag_interp_weights(tj); % barycentric Lagrange interpolation wts

sh = a; % shift

% basis coefficients, non-shifted, i.e. % c_k coeffs of f
if use_bjorck_pereyra
    c = dvand(tj,fj);
else
    V = ones(n,n); for k=2:n, V(:,k) = V(:,k-1).*(tj-0); end % Vandermonde not transp
    c = V\fj(:);
end

% basis coefficients, shifted, i.e. d_k coeffs of f
if strcmp(solve_shifted,'qp') % quadruple precision samples
    if use_bjorck_pereyra
        d = dvand(tjqp-sh,fjqp);
    else
        Vsh = ones(n,n); for k=2:n, Vsh(:,k) = Vsh(:,k-1).*(tjqp-sh); end % Vandermonde not transp, and shifted
        d = Vsh\fjqp(:);
    end
elseif strcmp(solve_shifted,'mp') % "mixed" precision samples
    if use_bjorck_pereyra
        d = dvand(tj-sh,fjmp);
    else
        Vsh = ones(n,n); for k=2:n, Vsh(:,k) = Vsh(:,k-1).*(tj-sh); end
        d = Vsh\fjmp(:);
    end
elseif strcmp(solve_shifted,'dp') % double precision samples
    if use_bjorck_pereyra
        d = dvand(tj-sh,fj);
    else
        Vsh = ones(n,n); for k=2:n, Vsh(:,k) = Vsh(:,k-1).*(tj-sh); end
        d = Vsh\fj(:);
    end
else
    error('invalid parameter solve_shifted');
end
% d = taylor_from_monomial(c,sh); % shifted coeffs from non-shifted coeffs c_k
% correct first coefficient in shifted basis
if corr_coeff_exact
    d(1) = f(sh); % try repair d(1) ! -> works (up to p3sh err), shows the culprit
end
if corr_coeff_interp_sig
    sigma_interp = bclag_interp(sj,tj,wbary,sh); % barycentric Lagrange interp of layer dens at t=a
    d(1) = sigma_interp*h(sh)^r; % eval remaining part of kernel analytically
end
if corr_coeff_exact && corr_coeff_interp_sig
    error('invalid parameters, choose one correction');
end

% loop over imag val of root
M = numel(bv);
irefv = zeros(M,1);
errv = zeros(M,1);
errshv = zeros(M,1);
p3mat = zeros(n,M);
p3shmat = zeros(n,M);
for i = 1:M
    b = bv(i);
    t0 = a + 1i*b; % root in t

    R = @(t) abs(t-t0); % dist func
    integrand3 = @(t) f(t)./R(t).^3;
    if abs(a) < 1
        I3e1 = integral(integrand3, -1, a, 'reltol', 0, 'abstol', 0);
        I3e2 = integral(integrand3, a, 1, 'reltol', 0, 'abstol', 0);
        I3e = I3e1+I3e2;
    else
        I3e = integral(integrand3, -1, 1, 'reltol', 0, 'abstol', 0);
    end
    irefv(i) = I3e;

    % basis integrals
    if use_vpa
        [p1, p3, p5] = rsqrt_pow_integrals(vpa(t0), n); % recurrences
        [p1sh, p3sh, p5sh] = rsqrt_pow_integrals_shift(vpa(t0), n); % w shift
    else
        [p1, p3, p5] = rsqrt_pow_integrals(t0, n);
        [p1sh, p3sh, p5sh] = rsqrt_pow_integrals_shift(t0, n);
    end
    p3mat(:,i) = p3;
    p3shmat(:,i) = p3sh;
    
    % plain non-adj method, non-shifted
    I3 = sum(c.*p3);
    errv(i) = abs((I3-I3e)/I3e);

    % plain non-adj method, shifted
    I3sh = sum(d.*p3sh);
    errshv(i) = abs((I3sh-I3e)/I3e);

    % print vectors
    if i == 1
        Pkm = p3;
        c_Pkm = c.*Pkm;
        sum_c_Pkm = nan*ones(numel(c_Pkm),1);
        sum_c_Pkm(1) = sum(c_Pkm);
        table(c,Pkm,c_Pkm,sum_c_Pkm)

        Qkm = p3sh;
        d_Qkm = d.*Qkm;
        sum_d_Qkm = nan*ones(numel(d_Qkm),1);
        sum_d_Qkm(1) = sum(d_Qkm);
        table(d,Qkm,d_Qkm,sum_d_Qkm)
    end
end

end

function test_corrections
savefig = 0; % saves figures to folder images

use_vpa = 1; % compute recurrences using vpa in digits(32)
corr_coeff_exact = 0; % correct first poly coeff by exact value
corr_coeff_interp_sig = 1; % correct using interp of layer dens
solve_shifted = 'dp'; % solve shifted Vandermonde system in dp/mp/qp
use_bjorck_pereyra = 1; % or solve Vandermonde systems using "\"

a = 0.23; % real of root
bv = logspace(-5,0,20)'; % imag of root, "distance" to panel
delta = 0; % offset 1e-3

sigma = @(t) sin(t+1); % artifical layer density function
r = 2; % vanishing rate

n = 16; % nodes

% run test
[errv,errshv,c,d,p3mat,p3shmat,irefv] = test_flat_panel_basis_corr(...
    sigma,r,delta,a,bv,n,use_vpa,corr_coeff_exact,corr_coeff_interp_sig,...
    solve_shifted,use_bjorck_pereyra);

% prepare plots
M = numel(bv);
[tj,~] = gauss(n);
h = @(t) t-a+delta;
f = @(t) sigma(t).*h(t).^r;
fj = f(tj);

% Plots
close all;

figure;
loglog(abs(bv),errv,'.');
hold on;
loglog(abs(bv),errshv,'o');
loglog(abs(bv),1e-16*1./abs(bv).^2,'k');
grid on;
legend('non-shifted','shifted','O(1/b^2)');
xlabel('b');
ylabel('rel err');
title(['use vpa=' num2str(use_vpa) ', corr coeff ex=' num2str(corr_coeff_exact), ', corr coeff interp=' num2str(corr_coeff_interp_sig) ', delta=' num2str(delta)]);

figure;
loglog(abs(bv),irefv,'.');
grid on;
xlabel('b');
ylabel('ref intval, I');
title(['delta=' num2str(delta)]);

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
semilogy(1:n,abs(p3mat),'.--');
grid on;
xlim([1,n]);
xlabel('k');
ylabel('|Pk3|');
bleg = cell(M,1);
for i = 1:M
    bleg{i} = ['b=' num2str(bv(i))];
end
legend(bleg,'location','northoutside','NumColumns',3);
title(['a=' num2str(a)]);

figure;
semilogy(1:n,abs(p3shmat),'o--');
grid on;
xlim([1,n]);
xlabel('k');
ylabel('|Qk3|');
legend(bleg,'location','northoutside','NumColumns',3);
title(['a=' num2str(a)]);

figure;
semilogy(1:n,abs(c),'.--');
hold on;
semilogy(1:n,abs(d),'o--');
grid on;
xlim([1,n]);
xlabel('k');
ylabel('absolute value of coeff');
legend('non-shifted: c_k','shifted: d_k');

figure;
semilogy(1:n,abs(d)./max(abs(d)),'.');
xlabel('Polynomial coefficient number');
ylabel('Relative magnitude of coefficient');
title(['Decay of shifted d_k coefficients, delta=' num2str(delta)]);
grid on;
xlim([1,n]);

alignfigs;

if savefig
    disp('saving figures...');
    exportgraphics(figure(1),['images/err_vpa' num2str(use_vpa) '_corrcoeff' num2str(corr_coeff_interp_sig) '_delta' num2str(delta) '_a' num2str(a) '.pdf'],'Resolution',400);
    exportgraphics(figure(2),['images/ivalref_delta' num2str(delta) '_a' num2str(a) '.pdf'],'Resolution',400);
    exportgraphics(figure(3),['images/interval_tarpts_a' num2str(a), '.pdf'],'Resolution',400);
    exportgraphics(figure(4),['images/numerator_delta' num2str(delta) '_a' num2str(a), '.pdf'],'Resolution',400);
    exportgraphics(figure(5),['images/Pk3_a=' num2str(a), '.pdf'],'Resolution',400);
    exportgraphics(figure(6),['images/Qk3_a=' num2str(a), '.pdf'],'Resolution',400);
    disp('sucessfully saved figures');
end
end
