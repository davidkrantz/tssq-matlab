function [errv,errshv,errestv,c,d,pmat,pshmat,irefv] = test_flat_panel_basis_corr(...
    m,sigma,r,delta,a,bv,n,use_vpa,no_hp_switch,adj_method,corr_coeff_exact,corr_coeff_interp_sig,...
    solve_nonshifted,solve_shifted,use_bjorck_pereyra,errest_alt)
% TEST_FLAT_PANEL_BASIS_CORR compares non-shifted and shifted monomial
%   basis functions to compute
%   
%   I^m(t_0) = \int_{-1}^{1} \frac{f(t)}{|t-t_0|^m} dt, m=1,3,5
%
% [errv,errshv,c,d,pmat,pshmat,irefv] = test_flat_panel_basis_corr(...
%   sigma,r,delta,a,bv,n,use_vpa,corr_coeff_exact,corr_coeff_interp_sig,...
%   solve_shifted,use_bjorck_pereyra) returns the relative error of the two
%   basis choices as well as the corresponding basis coefficients and basis
%   integral values.
%
% Without arguments runs test defined below
%
% INPUTS:
%   m                     - power of singularity 1/r^m
%   sigma                 - function handle for the layer density
%   r                     - vanishing order of h(t) at t = a
%   delta                 - small offset for location of vanishing, h(t) = t - a + delta
%   a                     - real part of singularity location t0
%   bv                    - vector of imaginary parts of singularities t0 = a + ib
%   n                     - number of Gauss-Legendre nodes (basis function order)
%   use_vpa               - string, 'none'/'all'/'init' indicating whether to use vpa to compute basis integrals Pkm and Qkm
%   no_hp_switch          - boolean, switch to disable half plane switch for I1(1)
%   adj_method            - if true, solves non-shifted using adjoint method
%   corr_coeff_exact      - if true, correct first basis coeff in shifted case with exact f(a)
%   corr_coeff_interp_sig - if true, correct first coeff using interpolated sigma(a)
%   solve_nonshifted      - string, 'dp'/'qp' for double/quad precision solve
%   solve_shifted         - string, 'dp'/'mp'/'qp' for double/mixed/quad precision solve
%   use_bjorck_pereyra    - boolean, use Björck-Pereyra method to solve Vandermonde systems
%   errest_alt            - integer determines how to est cond number of sum
%
% OUTPUTS:
%   errv    - relative absolute error in non-shifted basis evaluation for each b
%   errshv  - relative absolute error in shifted basis evaluation for each b
%   errestv - cancellation error estimate
%   c       - monomial basis coefficients (non-shifted)
%   d       - monomial basis coefficients (shifted)
%   pmat    - matrix of integrals P_k^m for non-shifted basis
%   pshmat  - matrix of integrals Q_k^m for shifted basis
%   irefv   - reference integral values computed via adaptive quadrature
%
% AUTHOR: David Krantz (davkra@kth.se), April 2025
%
% NOTE: Based on an idea by Alex Barnett. Some code has been taken from the 
%   GitHub repository https://github.com/ludvigak/linequad

if nargin == 0, test_corrections; return; end

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

% basis coefficients, non-shifted, i.e. c_k coeffs of f
if strcmp(solve_nonshifted,'qp')
    if use_bjorck_pereyra
        c = dvand(tjqp,fjqp);
    else
        V = ones(n,n); for k=2:n, V(:,k) = V(:,k-1).*(tjqp-0); end % Vandermonde not transp
        c = V\fjqp(:);
    end
elseif strcmp(solve_nonshifted,'dp')
    if use_bjorck_pereyra
        c = dvand(tj,fj);
    else
        V = ones(n,n); for k=2:n, V(:,k) = V(:,k-1).*(tj-0); end % Vandermonde not transp
        c = V\fj(:);
    end
else
    error('invalid parameter solve_nonshifted');
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
    d(1) = f(sh); % correct d(1) ! -> works (up to basis int err)
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
    if strcmp(use_vpa,'none')
        [p1, p3, p5] = rsqrt_pow_integrals(t0, n); % recurrences
        [p1sh, p3sh, p5sh] = rsqrt_pow_integrals_shift(t0, n, no_hp_switch); % w shift
    else
        [p1, p3, p5] = rsqrt_pow_integrals(vpa(t0), n); % use vpa for all terms
        [p1sh, p3sh, p5sh] = rsqrt_pow_integrals_shift(vpa(t0), n, no_hp_switch, use_vpa);
    end
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

    % estimate cancellation error estimate using adj method vectors
    if ~adj_method
        if use_bjorck_pereyra
            lam = pvand(tj,p);
        else
            lam = V'\p;
        end
    end
    adjv = lam.*fj;
    
    % various estimation methods
    if errest_alt == 1
        % standard
        kappa = sum(abs(adjv))/abs(I);
        kappa = kappa*n;
    elseif errest_alt == 2
        % upper bound on summation terms
        kappa = (max(abs(adjv))*n)/abs(I);
        kappa = kappa*n;
    elseif errest_alt == 3
        % average over k random samples
        k = 2;
        rind = randi(n,k,1);
        kappa = (n/k)*sum(abs(adjv(rind)))/abs(I);
        kappa = kappa*n;
    elseif errest_alt == 4
        % condition number of dot product using Euclidean norm
        kappa = norm(lam,2)*norm(fj,2)/abs(I);
    end
    % final cancellation error estimate
    errestv(i) = kappa*eps;

    % print vectors
    if i == 1
        Pkm = p;
        c_Pkm = c.*Pkm;
        sum_c_Pkm = nan*ones(numel(c_Pkm),1);
        sum_c_Pkm(1) = sum(c_Pkm);
        table(c,Pkm,c_Pkm,sum_c_Pkm)

        Qkm = psh;
        d_Qkm = d.*Qkm;
        sum_d_Qkm = nan*ones(numel(d_Qkm),1);
        sum_d_Qkm(1) = sum(d_Qkm);
        table(d,Qkm,d_Qkm,sum_d_Qkm)
    end
end


% condition number of Vandermonde matrices, just for testing
V = ones(n,n); for k=2:n, V(:,k) = V(:,k-1).*(tj-0); end
Vsh = ones(n,n); for k=2:n, Vsh(:,k) = Vsh(:,k-1).*(tj-sh); end
fprintf('cond nbr of unshifted Vandermonde matrix: %.14e\n',cond(V));
fprintf('cond nbr of shifted Vandermonde matrix: %.14e\n',cond(Vsh));

end

function test_corrections

test_no = 1; % which test to run

savefig = 0; % saves figures to folder matlab/images

m = 3; % power of singularity 1/r^m
use_vpa = 'none'; % compute recurrences using vpa in digits(32)
no_hp_switch = 0; % disables half plane switch for I1(1)
adj_method = 1; % solve non-shifted using adjoint method
corr_coeff_exact = 0; % correct first poly coeff by exact value
corr_coeff_interp_sig = 1; % correct using interp of layer dens
solve_nonshifted = 'dp'; % solve non-shifted Vandermonde system in dp/qp
solve_shifted = 'dp'; % solve shifted Vandermonde system in dp/mp/qp
use_bjorck_pereyra = 1; % or solve Vandermonde systems using "\"
errest_alt = 4; % which way to estimate cancellation error

switch test_no
    case 1
        a = 0.23; % real of root
        bv = logspace(-5,0,20)'; % imag of root, "distance" to panel
        delta = 0; % offset
        sigma = @(t) sin(t+1); % artifical layer density function
        r = 2; % vanishing rate
        n = 16; % nodes
    case 2
        a = -0.63;
        bv = -logspace(-5,0,20)';
        delta = 1e-3;
        sigma = @(t) sin(2*t+1)-3*cos(t);
        r = 3;
        n = 22;
    otherwise
        error('test does not exist');
end

% run test
[errv,errshv,errestv,c,d,pmat,pshmat,irefv] = test_flat_panel_basis_corr(...
    m,sigma,r,delta,a,bv,n,use_vpa,no_hp_switch,adj_method,...
    corr_coeff_exact,corr_coeff_interp_sig,solve_nonshifted,...
    solve_shifted,use_bjorck_pereyra,errest_alt);

% prepare plots
M = numel(bv);
[tj,~] = gauss(n);
h = @(t) t-a+delta;
f = @(t) sigma(t).*h(t).^r;
fj = f(tj);

% prints
fprintf('max relative error standard basis: %.14e\n', max(errv));
fprintf('max relative error modified basis: %.14e\n', max(errshv));

% plots
close all;

figure;
loglog(abs(bv),errv,'.');
hold on;
loglog(abs(bv),errshv,'o');
loglog(abs(bv),errestv,'mx-');
loglog(abs(bv),1e-16*1./abs(bv).^2);
grid on;
legend('non-shifted','shifted',['errest (type=' num2str(errest_alt) ')'],'O(1/b^2)');
xlabel('b');
ylabel('rel err');
title(['m=' num2str(m) ', ns adj=' num2str(adj_method) ', vpa=' use_vpa ', fix Q_1^1=' num2str(~no_hp_switch), ', corrcoeff=' num2str(corr_coeff_interp_sig) ', delta=' num2str(delta)]);

figure;
loglog(abs(bv),vecnorm(c.*pmat,inf,1),'.');
hold on;
loglog(abs(bv),vecnorm(d.*pshmat,inf,1),'o');
loglog(abs(bv),abs(irefv),'.');
loglog(abs(bv),1./abs(bv).^2);
if m == 5
    loglog(abs(bv),1./abs(bv).^4);
end
grid on;
if m == 5
    legend('inf-norm, quad vec nonshifted','inf-norm, quad vec shifted','|Iref|','O(1/b^2)','O(1/b^4)');
else
    legend('inf-norm, quad vec nonshifted','inf-norm, quad vec shifted','|Iref|','O(1/b^2)');
end
xlabel('b');
ylabel('value');
title(['m=' num2str(m) ', delta=' num2str(delta)]);

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
semilogy(1:n,abs(pmat),'.--');
grid on;
xlim([1,n]);
xlabel('k');
ylabel('|Pkm|');
bleg = cell(M,1);
for i = 1:M
    bleg{i} = ['b=' num2str(bv(i))];
end
legend(bleg,'location','northoutside','NumColumns',3);
title(['m=' num2str(m) ', a=' num2str(a)]);

figure;
semilogy(1:n,abs(pshmat),'o--');
grid on;
xlim([1,n]);
xlabel('k');
ylabel('|Qkm|');
legend(bleg,'location','northoutside','NumColumns',3);
title(['m=' num2str(m) ', a=' num2str(a)]);

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
    exportgraphics(figure(1),['matlab/images/err_m' num2str(m) '_vpa_' use_vpa '_fix_init' num2str(~no_hp_switch) '_corrcoeff' num2str(corr_coeff_interp_sig) '_delta' num2str(delta) '_a' num2str(a) '_erresttype' num2str(errest_alt) '.pdf'],'Resolution',400);
    exportgraphics(figure(2),['matlab/images/ivalref_m' num2str(m) '_delta' num2str(delta) '_a' num2str(a) '.pdf'],'Resolution',400);
    exportgraphics(figure(3),['matlab/images/interval_tarpts_a' num2str(a), '.pdf'],'Resolution',400);
    exportgraphics(figure(4),['matlab/images/numerator_delta' num2str(delta) '_a' num2str(a), '.pdf'],'Resolution',400);
    exportgraphics(figure(5),['matlab/images/Pk' num2str(m) '_a=' num2str(a), '.pdf'],'Resolution',400);
    exportgraphics(figure(6),['matlab/images/Qk' num2str(m) '_a=' num2str(a), '.pdf'],'Resolution',400);
    disp('sucessfully saved figures');
end
end
