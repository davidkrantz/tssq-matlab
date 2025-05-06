function [quadmat,recmat] = test_basis_int_recurrences(m,l,rvec,kvec)
% TEST_MOD_BASIS_INT_RECURRENCES tests the accuracy in the new and old
%   recurrence formulas for computing the standard and modified basis 
%   integrals, respectively. NOTE: This script does not compute the actual
%   basis integrals
%   
%   Q_k^m(t_0) = \int_{0}^{2\pi}
%   \frac{\sin((t-a)/2)^le^{ikt}}{|e^{it}-e^{it_0}|^m} dt, m=3,5,
%
%   but instead computes the related integral that only depends on 0<r<1,
%   
%   \tilde{Q}_k^m(r) = \frac{2}{(1-r)^m} \int_{0}^{\pi}
%   \frac{\sin(t/2)^2\cos(kt)}{(1-\ell(r)\sin(t/2)^2)^{m/2}}, m=3,5,
%
%   where \ell(r) = -(4r)/(1-r)^2.
%
% [quadmat,recmat] = test_mod_basis_int_recurrences(m,rvec,kvec)
%   returns the values of \tilde{Q}_k^m(r) for the r-values in rvec and
%   wavenumbers k in kvec.
%
% Without arguments runs test defined below
%
% INPUTS:
%   m    - power of singularity 1/r^m
%   l    - power of sine functions in basis (l=0 --> std, l=2 --> mod)
%   rvec - r-values, 0<r<1, to test
%   kvec - wavenumbers k to test, must be >0 in this script
%
% OUTPUTS:
%   quadmat   - integral values computed using "integral"
%   recmat    - integral values computed using recurrence formulas
%
% AUTHOR: David Krantz (davkra@kth.se), May 2025

if nargin == 0, test_recurrences; return; end

assert(sum(kvec<1)==0,'k must be >0 (not in general, just in this test)');

M = numel(kvec);
N = numel(rvec);
quadmat = zeros(M,N);
recmat = zeros(M,N);
parfor i = 1:M
    disp(['running k-value ' num2str(i) ' of ' num2str(M)]);
    k = kvec(i);
    for j = 1:N
        r = rvec(j);
        quadmat(i,j) = integral(@(t) 2*(sin(t./2).^l.*cos(k*t))./(1+(4*r)./(1-r).^2*sin(t./2).^2).^(m/2)./(1-r).^m,...
            0,pi,'AbsTol',0,'RelTol',0);
        recmat(i,j) = rec_test(r,k,m,l);
    end
end
end

function test_recurrences

test_no = 1; % which test to run

savefig = 0; % saves figures to folder matlab/images

switch test_no
    case 1
        m = 3; % power of singularity 1/r^m
        l = 0; % standard Fourier basis
        kmax = 80; % maximum wavenumber
        rvec = 1-logspace(-5,-1,100).'; % 0<r<1, log-spaced close to 1
    case 2
        m = 3;
        l = 2; % modified Fourier basis
        kmax = 80;
        rvec = 1-logspace(-5,-1,100).';
    case 3
        m = 5;
        l = 0;
        kmax = 80;
        rvec = 1-logspace(-5,-1,100).';
    case 4
        m = 5;
        l = 2;
        kmax = 80;
        rvec = 1-logspace(-5,-1,100).';
    otherwise
        error('test does not exist');
end

kvec = (1:kmax).'; % all wavenumbers

[quadmat,recmat] = test_basis_int_recurrences(m,l,rvec,kvec);

errmat = abs(quadmat-recmat)./abs(quadmat);

% plots
close all;

figure;
surf(rvec,kvec,quadmat)
xlabel('r');
ylabel('k');
shading interp;
title(['l=' num2str(l) ', m=', num2str(m) ': "integral"']);

figure;
surf(rvec,kvec,recmat)
xlabel('r');
ylabel('k');
shading interp;
title(['l=' num2str(l) ', m=', num2str(m) ': recurrence']);

figure;
surf(rvec,kvec,log10(errmat));
xlabel('r');
ylabel('k');
zlabel('rel abs err (log10)');
shading interp;
title(['l=' num2str(l) ', m=', num2str(m) ': rel abs err']);

leg = cell(numel(kvec),1);
for i = 1:numel(kvec)
    leg{i} = ['k=' num2str(kvec(i))];
end
figure;
loglog(1-rvec,errmat,'.');
grid on;
xlabel('1-r');
ylabel('rel abs err');
legend(leg,'Location','northoutside','NumColumns',6);
title(['l=' num2str(l) ', m=', num2str(m)]);

alignfigs;

if savefig
    disp('saving figures...');
    exportgraphics(figure(1),['matlab/images/basis_int_l' num2str(l) '_m' num2str(m) '_integral.pdf'],'Resolution',400);
    exportgraphics(figure(2),['matlab/images/basis_int_l' num2str(l) '_m' num2str(m) '_rec.pdf'],'Resolution',400);
    exportgraphics(figure(3),['matlab/images/basis_int_l' num2str(l) '_m' num2str(m) '_relabserr_surf.pdf'],'Resolution',400);
    exportgraphics(figure(4),['matlab/images/basis_int_l' num2str(l) '_m' num2str(m) '_relabserr_scatter.pdf'],'Resolution',400);
    disp('sucessfully saved figures');
end
end

function mu = rec_test(r,kk,m,l)
% naively implemented, will be rewritten later, new recurrences only valid
% for m>1 and k>0
kvec = 0:kk;
M = length(kvec);
[K,E] = ellipke(r^2);

% p = 1/2
P(1,1) = 2*K;
P(2,1) = 2/r*(-E+K);
% p = 3/2
P(1,2) = 2/(1+r)*(2/(1+r)*E-(1-r)*K);
% p = 5/2
P(1,3) = 2/(3*(1+r)^4)*(8*(1+r^2)*E-(1-r)*(1+r)*(5+3*r^2)*K);

% Recurrence for p = 1/2
for i = 3:M
    k = kvec(i);
    P(i,1) = (1+r^2)*2*(k-1)/(2*k-1)/r*P(i-1,1) - (2*k-3)/(2*k-1)*P(i-2,1);
end

% Recurrence for p = 3/2, 5/2
plist = [3/2;5/2];
for j = 1:2
    ptmp = plist(j);
    for i = 2:M
        k = kvec(i);
        P(i,j+1) = (1+r^2)/2/r*P(i-1,j+1)-(1-r)^2*(ptmp+k-2)/(ptmp-1)/2/r*P(i-1,j);
    end
end

% standard Fourier basis
if l == 0
    if m == 3
        mu = 2*P(end,2)/(1-r)^(m-1);
    else
        mu = 2*P(end,3)/(1-r)^(m-1);
    end
    return;
end

% modified Fourier basis, recurrences based on old ones
if m == 3
    mu = 0.5 * (-(1-r)^2/(2*r*(1+r^2))*P(end,2) + ((m/2+kk-1)/(2*r)*P(end,1)-(m/2+kk-2)/(1+r^2)*P(end-1,1))/(m/2-1));
else
    mu = 0.5*(1-r)^(3-m) * (-(1-r)^2/(2*r*(1+r^2))*P(end,3) + ((m/2+kk-1)/(2*r)*P(end,2)-(m/2+kk-2)/(1+r^2)*P(end-1,2))/(m/2-1));
end
end
