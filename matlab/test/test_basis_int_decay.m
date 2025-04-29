function [Pmat,Qmat,kvec] = test_basis_int_decay(m,a,bv,K,basis_alt)
% TEST_CIRCLE_BASIS_INT_DECAY computes the basis integrals
%
%   I^m(t_0) = \int_{0}^{2\pi} \frac{\Phi(t)}{|e^{it}-e^{it_0}|^m} dt, m=1,3,5
%
% using different basis functions \Phi(t), based on the variable basis_alt.
%
% [Pmat,Qmat,kvec] = test_basis_int_decay(m,a,bv,n,basis_alt) returns the
%   standard and modified basis integrals and well as the corresponding
%   wavenumbers.
%
% Without arguments runs test defined below
%
% INPUTS:
%   m         - power of denominator
%   a         - real part of root t_0
%   bv        - imag part of root t_0
%   K         - maximum wavenumber
%   basis_alt - which basis to test
%
% OUTPUTS:
%   Pmat - NxM matrix with standard basis integrals
%   Qmat - NxM matrix with modified basis integrals
%   kvec - wavenumbers of standard Fourier basis
%
% AUTHOR: David Krantz (davkra@kth.se), April 2025

if nargin == 0, test_basis; return; end

kvec = (-K:K).';
N = numel(kvec);

% loop over imag val of root
M = numel(bv);
Pmat = zeros(N,M);
Qmat = zeros(N,M);
for i = 1:M
    b = bv(i);
    t0 = a + 1i*b; % root in t

    R = @(t) abs(exp(1i*t)-exp(1i*t0)); % dist func
    int = @(t,k) exp(1i*k*t)./R(t).^m; % unshifted integrand

    switch basis_alt
        case 1
            intsh = @(t,k) (exp(1i*t)-exp(1i*a)).^k./R(t).^m; % shifted integrand
        case 2
            intsh = @(t,k) exp(1i*(t-a)).^k./R(t).^m;
        case 3
            intsh = @(t,k) exp(1i*(t-a-1i*b)).^k./R(t).^m;
        case 4
            intsh = @(t,k) (exp(1i*k*t)*exp(-k*b))./R(t).^m;
        case 5
            sigma = 1e-1;
            w = @(t) exp(-(t-a).^2./(2*sigma^2));
            intsh = @(t,k) (w(t).*exp(1i*k*t))./R(t).^m;
        case 6
            z0 = exp(1i*t0);
            intsh = @(t,k) ((exp(1i*t)-z0)./(1-conj(z0)*exp(1i*t))).^k./R(t).^m;
        case 7
            intsh = @(t,k) (exp(1i*t)-exp(1i*t0)).^k./R(t).^m;
        case 8
            intsh = @(t,k) sin((t-a)/2).^k./R(t).^m;
        case 9
            intsh = @(t,k) (exp(1i*k*t).*sin((t-a)/2).^2)./R(t).^m;
        otherwise
            error('basis alternative not supported')
    end

    warning off;
    for j = 1:N
        k = kvec(j);
        Pmat(j,i) = integral(@(t) int(t,k), 0, 2*pi, 'reltol', 0, 'abstol', 0);
        Qmat(j,i) = integral(@(t) intsh(t,k), 0, 2*pi, 'reltol', 0, 'abstol', 0);
    end
    warning on;
end

end

function test_basis

test_no = 1; % which test to run

m = 3; % power of singularity 1/r^m
basis_alt = 9; % which basis to use

switch test_no
    case 1
        a = 0.23; % real of root
        bv = logspace(-5,0,20)'; % imag of root, "distance" to panel
        K = 10; % max wavenumber to use
    case 2
        a = -0.63;
        bv = -logspace(-5,0,20)';
        K = 22;
    otherwise
        error('test does not exist');
end

% run test
[Pmat,Qmat,kvec] = test_basis_int_decay(m,a,bv,K,basis_alt);

% print
table(Pmat(:,1),Qmat(:,1))

% plot
close all;

leg = cell(numel(bv),1);
for i = 1:numel(bv)
    leg{i} = ['b=' num2str(bv(i))];
end

figure;
semilogy(kvec,abs(Pmat),'.--');
grid on;
xlabel('wavenumber, k');
ylabel('|Pkm|');
title(['standard Fourier basis: a=' num2str(a)]);
legend(leg,'location','northoutside','NumColumns',3);

figure;
semilogy(kvec,abs(Qmat),'.--');
grid on;
xlabel('wavenumber, k');
ylabel('|Qkm|');
title(['modified basis: a=' num2str(a)]);
legend(leg,'location','northoutside','NumColumns',3);

alignfigs;

end
