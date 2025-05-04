function [errv,errmodv,errestv,c,d,pmat,pmodmat,irefv,kstd,kmod,...
    maxerrbnd_Ie,maxerrbnd_p,maxerrbnd_pmod] = test_circle_basis_corr(...
    m,l,sigma,alpha,delta,a,bv,n,corr_coeff,add_shifted_sine,use_vpa,use_quadgk)
% TEST_CIRCLE_BASIS_CORR compares the standard and modified Fourier basis
%   to compute
%   
%   I^m(t_0) = \int_{0}^{2\pi} \frac{f(t)}{|e^{it}-e^{it_0}|^m} dt, m=1,3,5
%
% [errv,errmodv,errestv,c,d,pmat,pmodmat,irefv,maxerrbnd_Ie,...
%   maxerrbnd_p,maxerrbnd_pmod,kstd,kmod] = test_circle_basis_corr(...
%   m,l,sigma,alpha,delta,a,bv,n,corr_coeff,add_shifted_sine,use_vpa,use_quadgk)
%   returns the relative error of the two basis choices as well as the
%   basis coefficients and basis integral values.
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
%   use_quadgk       - boolean, if true, computes standard basis integrals 
%                      using quadgk instead of recurrence formulas
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
%   maxerrbnd_Ie    - max error estimate bound of quadgk for reference values
%   maxerrbnd_p     - max error estimate bound of quadgk for P_k^m
%   maxerrbnd_pmod  - same but for Q_k^m
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

V = sin((tj-a)/2).^l.*exp(1i*kmod.'.*tj); % modified Fourier basis
if add_shifted_sine
    V = [ones(n,1) sin(tj-a) V]; % add constant function 1 and shifted sine
else
    V = [ones(n,1) V]; % add constant function 1
end
d = V\fj; % modified Fourier basis coefficients

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
[maxerrbnd_Ie,maxerrbnd_p,maxerrbnd_pmod] = deal(1e-100);
warning off;
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
    [Ie,errbndtmp] = quadgk(integrand,0,2*pi,'MaxIntervalCount',1e7,'AbsTol',0,'RelTol',0);
    maxerrbnd_Ie = max(maxerrbnd_Ie,errbndtmp);
    irefv(i) = Ie;

    % basis integrals, standard
    if use_quadgk
        g = @(t,k) cos(k*t)./(1-ell(double(r))*sin(t/2).^2).^(m/2);
        muk = zeros(numel(kstd),1);
        for j = 1:numel(kstd)
            [muk(j),errbndtmp] = quadgk(@(t) g(t,kstd(j)),0,pi,'MaxIntervalCount',1e5,'AbsTol',0,'RelTol',0);
            maxerrbnd_p = max(maxerrbnd_p,errbndtmp);
        end
        p = 2*exp(1i*kstd*a).*muk./(1-r).^m; % scale from rewriting og int
    else
        kmaxstd = max(abs(kstd));
        muk = periodic_basis_integrals(r,kmaxstd,n,K,E,m); % recurrences
        p = (2*exp(1i*kstd*a).*muk)./(1-r).^(m-1);
    end
    pmat(:,i) = p;

    % basis integrals, modified
    gmod = @(t,k) (sin(t/2).^l.*cos(k*t))./(1-ell(double(r))*sin(t/2).^2).^(m/2);
    pmod = zeros(numel(kstd),1);
    pmodk = zeros(numel(kmod),1);
    for j = 1:numel(kmod)
        [pmodk(j),errbndtmp] = quadgk(@(t) gmod(t,kmod(j)),0,pi,'MaxIntervalCount',1e5,'AbsTol',0,'RelTol',0);
        maxerrbnd_pmod = max(maxerrbnd_pmod,errbndtmp);
    end
    p0 = 2*(2/(1+r)*(2/(1+r)*E-(1-r)*K))/(1-r)^(m-1); % hardcoded for m=3
    pmod(1) = p0; % basis integral for constant term
    if add_shifted_sine
        pmod(2) = 0; % basis integral for sin(t-a) term equals zero
        pmod(3:end) = 2*exp(1i*kmod*a).*pmodk./(1-r).^m;
    else
        pmod(2:end) = 2*exp(1i*kmod*a).*pmodk./(1-r).^m;
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
warning on;

end

function test_corrections

test_no = 1; % which test to run

m = 3; % power of singularity 1/r^m
l = 2; % power of sine term in modified basis
use_vpa = 0; % use vpa for certain calculations
use_quadgk = 0; % compute standard basis integrals using quadgk instead of recurrences
add_shifted_sine = 1; % extends basis with sin(t-a)
corr_coeff = 'interp'; % determines how to correct the first modified coeff

switch test_no
    case 1
        n = 40; % nodes
        alpha = 2; % vanishing rate
        delta = 1e-6; % offset
        a = 4.23; % real of root
        bv = logspace(-5,0,4); % imag of root, "distance" to circle
        sigma = @(t) exp(cos(t)); % artifical layer density
    case 2
        n = 40;
        alpha = 2;
        delta = 0;
        a = 1.23;
        bv = logspace(-5,0,4);
        sigma = @(t) sin(2*t+1);
    otherwise
        error('invalid test number');
end

[errv,errmodv,errestv,c,d,pmat,pmodmat,irefv,kstd,kmod,...
    maxerrbnd_Ie,maxerrbnd_p,maxerrbnd_pmod] = test_circle_basis_corr(...
    m,l,sigma,alpha,delta,a,bv,n,corr_coeff,add_shifted_sine,use_vpa,use_quadgk);

fprintf('maxerrbnd_Ie=%0.20e\n',maxerrbnd_Ie)
fprintf('maxerrbnd_p=%0.20e\n',maxerrbnd_p)
fprintf('maxerrbnd_pmod=%0.20e\n',maxerrbnd_pmod)

% prepare plots
M = numel(bv);
tj = linspace(0,2*pi,n+1).'; tj(end) = [];
h = @(t) sin(t-a)+delta;
f = @(t) sigma(t).*h(t).^alpha;
fj = f(tj);

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
yline(abs(d(1)),'k','|d(1)|');
if add_shifted_sine
    loglog(bv,abs(pmodmat(2,:)),'.');
    yline(abs(d(2)),'k','|d(2)|');
    legend('|p(1)|','','|p(2)|');
else
    legend('|p(1)|');
end
xlabel('b');
ylabel('magnitude of coeff');
title('magnitude of constant coeff in modified basis');

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
