function [errv,errmodv,errestv,c,d,pmat,pmodmat,irefv,kstd,kmod] = test_circle_basis_corr(...
    m,l,sigma,alpha,delta,a,bv,n,corr_coeff,add_shifted_sine,use_vpa,use_integral,use_fft)
% TEST_CIRCLE_BASIS_CORR compares the standard and modified Fourier basis
%   to compute
%   
%   I^m(t_0) = \int_{0}^{2\pi} \frac{f(t)}{|e^{it}-e^{it_0}|^m} dt, m=3,5
%
% [errv,errmodv,errestv,c,d,pmat,pmodmat,irefv,kstd,kmod] = ...
%   test_circle_basis_corr(m,l,sigma,alpha,delta,a,bv,n,...
%   corr_coeff,add_shifted_sine,use_vpa,use_integral,use_fft) returns the relative 
%   error of the two basis choices as well as the basis coefficients and 
%   basis integral values.
%
% Without arguments runs test defined below
%
% INPUTS:
%   m                - power of singularity 1/r^m
%   l                - number of equidistant nodes on [0,2*pi) to sample f at
%   sigma            - function handle for the layer density
%   alpha            - vanishing order of h(t) at t = a, i.e. h(t)^alpha
%   delta            - small offset for location of vanishing, h(t) = sin((t-a)/2) + delta
%   a                - real part of singularity location t0
%   bv               - vector of imaginary parts of singularities t0 = a + ib
%   corr_coeff       - string, 'exact'/'interp' for exact or Fourier interp
%   add_shifted_sine - boolean, if true, extends basis with sin(t-a)
%   use_vpa          - boolean, for certain calculations
%   use_integral     - boolean, if true, computes standard basis integrals 
%                      using "integral" instead of recurrence formulas
%   use_fft          - boolean, if true, computes modified basis
%                      coefficients from the standard ones computed via FFT
%
% OUTPUTS:
%   errv            - relative absolute error in standard basis evaluation for each b
%   errmodv         - relative absolute error in modified basis evaluation for each b
%   errestv         - cancellation error estimate
%   c               - standard Fourier coefficients
%   d               - modified Fourier coefficients
%   pmat            - matrix of integrals P_k^m for standard basis
%   pmodmat         - matrix of integrals Q_k^m for modified basis
%   irefv           - reference integral values computed via adaptive quadrature
%   kstd            - wavenumbers of standard Fourier basis
%   kmod            - wavenumbers of modified Fourier basis
%
% AUTHOR: David Krantz (davkra@kth.se), May 2025

if nargin == 0, test_corrections; return; end

% numerator
h = @(t) sin(t-a)+delta; % function small at t=a
f = @(t) sigma(t).*h(t).^alpha; % modified "density", RR^T style vanishing near t=a

tj = linspace(0,2*pi,n+1).'; tj(end) = []; % nodes in [0,2*pi)

kstd = get_k_vec(n,2*pi).'; % wavenumbers for standard Fourier basis
if add_shifted_sine % remove modes due to extended basis, keeps sys square
    kmod = get_k_vec(n-2,2*pi).'; % wavenumbers for modified Fourier basis
else
    kmod = get_k_vec(n-1,2*pi).'; 
end

fj = f(tj); % samples
sj = sigma(tj); % only layer density function

c = fftshift(fft(fj))/n; % standard Fourier basis coefficients
normc = norm(c,2); % for cancellation error estimate

if use_fft
    [a0,a1,b_coeffs] = fourier2modcoeffs(c,a); % map c_k --> (a_0,a_1,b_k)
    d = [a0;a1;b_coeffs];
else
    V = sin((tj-a)/2).^l.*exp(1i*kmod.'.*tj); % modified Fourier basis
    if add_shifted_sine
        V = [ones(n,1) sin(tj-a) V]; % add constant function 1 and shifted sine
    else
        V = [ones(n,1) V]; % add constant function 1
    end
    d = V\fj; % modified Fourier basis coefficients
end

% correct first coefficient in modified basis
if strcmp(corr_coeff,'exact')
    d(1) = f(a);
elseif strcmp(corr_coeff,'interp')
    sjk = fftshift(fft(sj))/n; % Fourier coefficients of layer dens
    sigma_interp = exp(1i*kstd.'*a)*sjk; % interpolate layer dens at t=a
    d(1) = sigma_interp*h(a)^alpha; % eval remaining part of kernel analytically
end

ell = @(r) -(4*r)./(1-r).^2; % function in basis integrals

% loop over imag val of root
M = numel(bv);
irefv = zeros(M,1);
errv = zeros(M,1);
errestv = zeros(M,1);
errmodv = zeros(M,1);
pmat = zeros(n,M);
pmodmat = zeros(n,M);
for i = 1:M
    if use_vpa
        b = vpa(bv(i));
    else
        b = bv(i);
    end
    r = exp(-abs(b));
    t0 = a + 1i*b; % root in t
    [K,E] = ellipke(r^2); % elliptic integrals of first and second kind

    % integrand and reference value
    R = @(t) abs(exp(1i*t)-exp(1i*double(t0))); % dist func
    integrand = @(t) f(t)./R(t).^m;
    warning off;
    Ie1 = integral(integrand,0,a,'AbsTol',0,'RelTol',0);
    Ie2 = integral(integrand,a,2*pi,'AbsTol',0,'RelTol',0);
    warning on;
    Ie = Ie1+Ie2;
    irefv(i) = Ie;

    % basis integrals, standard
    if use_integral
        gstd = @(t,k) cos(k*t)./(1-ell(double(r))*sin(t/2).^2).^(m/2);
        muk = zeros(numel(kstd),1);
        for j = 1:numel(kstd)
            muk(j) = integral(@(t) gstd(t,kstd(j)),0,pi,'AbsTol',0,'RelTol',0);
        end
        p = 2*exp(1i*kstd*a).*muk./(1-r).^m; % scale from rewriting og int
    else
        kmaxstd = max(abs(kstd));
        [mu1k,mu3k,mu5k] = periodic_basis_integrals(r,kmaxstd,n,K,E); % recurrences
        if m == 3
            muk = mu3k;
        else
            muk = mu5k;
        end
        p = (2*exp(1i*kstd*a).*muk)./(1-r).^(m-1);
    end
    pmat(:,i) = p;

    % basis integrals, modified
    pmod = zeros(numel(kstd),1);
    if use_integral || ~add_shifted_sine
        gmod = @(t,k) (sin(t/2).^l.*cos(k*t))./(1-ell(double(r))*sin(t/2).^2).^(m/2);
        pmodk = zeros(numel(kmod),1);
        for j = 1:numel(kmod)
            pmodk(j) = integral(@(t) gmod(t,kmod(j)),0,pi,'AbsTol',0,'RelTol',0);
        end
        pmodk = 2*pmodk./(1-r).^m;
    else
        % compute as combination of the old basis integral values
        if m == 3
            pmodk = 0.5*(1-r)^(3-m) * (-(1-r)^2/(2*r*(1+r^2)).*mu3k(2:end-1) + ((m/2+kmod-1)./(2*r).*mu1k(2:end-1)-(m/2+kmod-2)./(1+r^2).*mu1k(1:end-2))./(m/2-1));
        else
            pmodk = 0.5*(1-r)^(3-m) * (-(1-r)^2/(2*r*(1+r^2)).*mu5k(2:end-1) + ((m/2+kmod-1)./(2*r).*mu3k(2:end-1)-(m/2+kmod-2)./(1+r^2).*mu3k(1:end-2))./(m/2-1));
        end
    end
    % basis integral for constant term
    if m == 3
        if use_integral
            p0 = 2*(2/(1+r)*(2/(1+r)*E-(1-r)*K))/(1-r)^(m-1);
        else
            if mod(n,2) == 0
                p0 = 2*mu3k(n/2+1)/(1-r)^(m-1);
            else
                p0 = 2*mu3k((n-1)/2+1)/(1-r)^(m-1);
            end
        end
    else
        if use_integral
            p0 = 2*(2/(3*(1+r)^4)*(8*(1+r^2)*E-(1-r)*(1+r)*(5+3*r^2)*K))/(1-r)^(m-1);
        else
            if mod(n,2) == 0
                p0 = 2*mu5k(n/2+1)/(1-r)^(m-1);
            else
                p0 = 2*mu5k((n-1)/2+1)/(1-r)^(m-1);
            end
        end
    end
    pmod(1) = p0;
    if add_shifted_sine
        pmod(2) = 0; % basis integral for sin(t-a) term always equals zero
        pmod(3:end) = exp(1i*kmod*a).*pmodk;
    else
        pmod(2:end) = exp(1i*kmod*a).*pmodk;
    end
    pmodmat(:,i) = pmod;

    % plain non-adj method, standard
    I = real(sum(c.*p));
    errv(i) = abs((I-Ie)/Ie);

    % plain non-adj method, modified
    Imod = real(sum(d(2:end).*pmod(2:end))) + d(1)*p0; % separate cuz vpa
    errmodv(i) = abs((Imod-Ie)/Ie);

    % cancellation error estimate
    errestv(i) = cond_sum(normc,p,c.*p); % cond nbr of summation operation

    % print vectors
    if i == 1
        Pkm = p;
        c_Pkm = c.*Pkm;
        sum_c_Pkm = nan*ones(numel(c_Pkm),1);
        sum_c_Pkm(1) = sum(c_Pkm);
        table(c,Pkm,c_Pkm,sum_c_Pkm)

        Qkm = pmod;
        d_Qkm = d.*Qkm;
        sum_d_Qkm = nan*ones(numel(d_Qkm),1);
        sum_d_Qkm(1) = sum(d_Qkm);
        table(d,Qkm,d_Qkm,sum_d_Qkm)
    end
end
end

function test_corrections

test_no = 1; % which test to run

m = 3; % power of singularity 1/r^m
l = 2; % power of sine term in modified basis
add_shifted_sine = 1; % extends basis with sin(t-a)
corr_coeff = 'interp'; % determines how to correct the first modified coeff
use_vpa = 0; % use vpa for certain calculations
use_integral = 0; % compute basis integrals using "integral" instead of rec
use_fft = 1; % compute modified basis coefficients via FFT

switch test_no
    case 1
        n = 60; % nodes
        alpha = 2; % vanishing rate
        delta = 1e-6; % offset
        a = 4.23; % real of root
        bv = logspace(-5,0,20); % imag of root, "distance" to circle
        sigma = @(t) exp(cos(t)); % artifical layer density
    case 2
        n = 41;
        alpha = 2;
        delta = 0;
        a = 1.23;
        bv = logspace(-5,0,20);
        sigma = @(t) sin(2*t+1);
    otherwise
        error('invalid test number');
end

[errv,errmodv,errestv,c,d,pmat,pmodmat,irefv,kstd,kmod,] = test_circle_basis_corr(...
    m,l,sigma,alpha,delta,a,bv,n,corr_coeff,add_shifted_sine,use_vpa,use_integral,use_fft);

% prepare plots
M = numel(bv);
tj = linspace(0,2*pi,n+1).'; tj(end) = [];
h = @(t) sin(t-a)+delta;
f = @(t) sigma(t).*h(t).^alpha;
fj = f(tj);

% prints
fprintf('max relative error standard basis: %.14e\n', max(errv));
fprintf('max relative error modified basis: %.14e\n', max(errmodv));

% plots
close all;

figure;
loglog(abs(bv),errv,'.');
hold on;
loglog(abs(bv),errmodv,'o');
loglog(abs(bv),errestv,'mx-');
loglog(abs(bv),1e-16*1./abs(bv).^2);
grid on;
legend('standard','modified','errest','O(1/b^2)');
xlabel('b');
ylabel('rel err');
title(['m=' num2str(m) ', vpa=' num2str(use_vpa) ', corrcoeff=' corr_coeff ', delta=' num2str(delta)]);

figure;
loglog(abs(bv),abs(irefv),'.');
grid on;
xlabel('b');
ylabel('abs ivalref, |I|');
title(['m=' num2str(m) ', delta=' num2str(delta)]);

figure;
tt = linspace(0,2*pi,1000);
plot(cos(tt),sin(tt),'k');
hold on;
plot(cos(tj),sin(tj),'.k');
z0 = exp(1i*(a+1i*bv));
plot(real(z0),imag(z0),'x');
grid on;
axis equal;
legend('unit circle','trapz nodes','tar pts');
xlabel('Real');
ylabel('Imag');
title(['a=' num2str(a)]);

figure;
plot(tt,f(tt));
hold on;
plot(tj,fj,'.');
xlim([0,2*pi]);
xline(a,'k-','t=a');
legend('numerator','samples');
xlabel('t');
ylabel('val');
grid on;
title(['delta=' num2str(delta)]);

figure;
semilogy(kstd,abs(pmat),'.--');
grid on;
xlim([min(kstd),max(kstd)]);
xlabel('k');
ylabel('|Pkm|');
bleg = cell(M,1);
for i = 1:M
    bleg{i} = ['b=' num2str(bv(i))];
end
legend(bleg,'location','northoutside','NumColumns',3);
title(['m=' num2str(m) ', a=' num2str(a)]);

figure;
if add_shifted_sine
    semilogy(kmod,abs(pmodmat(3:end,:)),'o--');
else
    semilogy(kmod,abs(pmodmat(2:end,:)),'o--');
end
grid on;
hold on;
xlim([min(kmod),max(kmod)]);
xlabel('k');
ylabel('|Qkm|');
legend(bleg,'location','northoutside','NumColumns',3);
title(['m=' num2str(m) ', a=' num2str(a)]);

figure;
loglog(bv,abs(pmodmat(1,:)),'.');
grid on;
hold on;
loglog(bv,abs(d(1)*pmodmat(1,:)),'.');
legend('|p(1)|','|d(1)p(1)|');
xlabel('b');
ylabel('magnitude of coeff');
title('constant coeff in modified basis');

figure;
semilogy(kstd,abs(c),'.--');
hold on;
if add_shifted_sine
    semilogy(kmod,abs(d(3:end)),'o--');
else
    semilogy(kmod,abs(d(2:end)),'o--');
end
grid on;
xlim([min(kstd),max(kstd)]);
xlabel('k');
ylabel('absolute value of coeff');
legend('standard: c_k','modified: d_k');

figure;
if add_shifted_sine
    semilogy(kmod,abs(d(3:end))./max(abs(d(3:end))),'.');
else
    semilogy(kmod,abs(d(2:end))./max(abs(d(2:end))),'.');
end
xlabel('wavenumber, k');
ylabel('Relative magnitude of coefficient');
title(['Decay of modified d_k coefficients, delta=' num2str(delta)]);
grid on;
xlim([min(kmod),max(kmod)]);

alignfigs;

end

function est = cond_sum(normw,f,wf)
kappa = normw * norm(f,2) / abs(sum(wf));
est = eps*kappa;
end

% old function, keep here just in case
% function mu = new_rec(r,kmax,nph,m)
% kvec = 0:kmax;
% M = length(kvec);
% [K,E] = ellipke(r^2);
% 
% % p = 1/2
% P(1,1) = 2*K;
% P(2,1) = 2/r*(-E+K);
% % p = 3/2
% P(1,2) = 2/(1+r)*(2/(1+r)*E-(1-r)*K);
% % p = 5/2
% P(1,3) = 2/(3*(1+r)^4)*(8*(1+r^2)*E-(1-r)*(1+r)*(5+3*r^2)*K);
% 
% % Recurrence for p = 1/2
% for i = 3:M
%     k = kvec(i);
%     P(i,1) = (1+r^2)*2*(k-1)/(2*k-1)/r*P(i-1,1) - (2*k-3)/(2*k-1)*P(i-2,1);
% end
% 
% % Recurrence for p = 3/2, 5/2
% plist = [3/2;5/2];
% for j = 1:2
%     ptmp = plist(j);
%     for i = 2:M
%         k = kvec(i);
%         P(i,j+1) = (1+r^2)/2/r*P(i-1,j+1)-(1-r)^2*(ptmp+k-2)/(ptmp-1)/2/r*P(i-1,j);
%     end
% end
% 
% if m == 3
%     mu0 = 2.*r.^(-1).*(1+r).^(-2).*((-1).*E+K+K.*r);
%     muktmp = 0.5*(1-r)^(3-m) * (-(1-r)^2/(2*r*(1+r^2)).*P(2:end,2) + ((m/2+kvec(2:end).'-1)./(2*r).*P(2:end,1)-(m/2+kvec(2:end).'-2)./(1+r^2).*P(1:end-1,1))./(m/2-1));
% else
%     mu0 = (-2/3).*((-1)+r).^(-2).*r.^(-1).*(1+r).^(-4).*(E+E.*((-6)+r).*r+K.*((-1)+r.*(3+r+(-3).*r.^2)));
%     muktmp = 0.5*(1-r)^(3-m) * (-(1-r)^2/(2*r*(1+r^2)).*P(2:end,3) + ((m/2+kvec(2:end).'-1)./(2*r).*P(2:end,2)-(m/2+kvec(2:end).'-2)./(1+r^2).*P(1:end-1,2))./(m/2-1));
% end
% 
% muktmp = [mu0; muktmp];
% 
% if kmax ~= nph/2
%     mu = [conj(muktmp(M:-1:2)); muktmp];
% else
%     mu = [conj(muktmp(M:-1:2)); muktmp(1:end-1)];
% end
% end
