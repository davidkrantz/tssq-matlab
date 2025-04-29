function [fstd,fmod,ck,dk,kstd,kmod,Vmod,ttj,ffj] = test_periodic_modified_basis(f,N,l,a)
% TEST_PERIODIC_MODIFIED_BASIS compares the performances of a modified
%   Fourier basis with the standard one. The modified basis has to be
%   hardcoded below.
%
% [fstd,fmod,ck,dk,kstd,kmod,Vmod,ttj,ffj] = 
%   test_periodic_modified_basis(f,N,l,a) returns the two reconstructions
%   and the corresponding basis coefficients.
%
% Without arguments runs test defined below
%
% INPUTS:
%   f - function handle to function that typically vanishes at t=a
%   N - number of equidistant nodes on [0,2*pi) to sample f at
%   l - power of sine term in modified basis
%   a - point in (0,2*pi) where f(a) vanishes or almost vanishes
%
% OUTPUTS:
%   fstd - reconstructed function values at points ttj using standard Fourier basis
%   fmod - reconstructed function values at points ttj using modified Fourier basis
%   ck   - standard Fourier coefficients
%   dk   - modified Fourier coefficients
%   kstd - wavenumbers of standard Fourier basis
%   kmod - wavenumbers of modified Fourier basis
%   Vmod - Vandermonde matrix used to determine dk
%   ttj  - random points in [0,2*pi) where reconstructions are compared to true function values
%   ffj  - the true function values at ttj, f(ttj)
%
% AUTHOR: David Krantz (davkra@kth.se), April 2025

if nargin == 0, test_basis; return; end

rng(123); % fixes random seed
Neval = 301; % nbr pts to evaluate reconstruction

tj = linspace(0, 2*pi, N+1).'; tj(end) = []; % periodic grid in [0,2*pi)
fj = f(tj); % function values at grid points
kstd = get_k_vec(N,2*pi).'; % wavenumbers

ck = fftshift(fft(fj))/N; % standard Fourier basis coefficients

kmod = get_k_vec(N-1,2*pi).'; % N-1 to account for added const fcs, makes sys square
Vmod = sin((tj-a)/2).^l.*exp(1i*kmod.'.*tj); % modified Fourier basis
Vmod = [ones(N,1) Vmod]; % add constant function 1
dk = Vmod\fj; % modified Fourier basis coefficients

% check that we resolve the function
if min(abs(dk./max(abs(dk)))) > 1e-8
    warning('function might not be resolved by modified Fourier basis');
end

ttj = 2*pi*rand(Neval,1); % random test points in [0,2*pi)
ttj = sort(ttj);
ffj = f(ttj); % true function values

% reconstruct at new points using standard Fourier basis
Vstd = exp(1i*kstd.'.*ttj);
fstd = Vstd*ck;
if max(abs(imag(fstd)))>1e-8
    warning('recon imag part large for standard basis');
end
fstd = real(fstd);

% reconstruct at new points using modified Fourier basis
Vmod_eval = sin((ttj-a)/2).^l.*exp(1i*kmod.'.*ttj);
Vmod_eval = [ones(Neval,1) Vmod_eval];
fmod = Vmod_eval*dk;
if max(abs(imag(fmod)))>1e-8
    warning('recon imag part large for modified basis');
end
fmod = real(fmod);
end

function test_basis

close all;

test_no = 3; % which test to run

switch test_no
    case 1
        N = 40; % pts on grid
        l = 2; % vanishing rate of modified basis
        a = 1.5; % real part of root
        delta = 0; % constant offset at t=a
        f = @(t) (sin((t-a)/2)).^4.*exp(cos(t))+delta; % smooth func vanishing at t=a
    case 2
        N = 20;
        l = 2;
        a = 3.2512;
        delta = 1e-6;
        f = @(t) cos((t-a-pi)/2).^2+delta + (sin((t-a)/2)).^4.*cos(t);
    case 3
        N = 40;
        l = 2;
        a = 5.5;
        delta = 0.1;
        tmp = 1:5;
        f = @(t) sum(cos(tmp)./tmp.^2).*sin((t-a)/2).^4+delta;
    otherwise
        error('invalid test number');
end

[fstd,fmod,ck,dk,kstd,kmod,Vmod,ttj,ffj] = test_periodic_modified_basis(f,N,l,a);

err_std = abs(fstd-ffj);
err_mod = abs(fmod-ffj);

fprintf('max abs err standard Fourier recon: %.2e\n',max(err_std));
fprintf('max abs err modified Fourier recon: %.2e\n',max(err_mod));
fprintf('cond number of Vandermonde matrix:  %.2e\n',cond(Vmod));

% plot
figure;
semilogy(kstd,abs(ck),'.--');
hold on;
semilogy(kmod,abs(dk(2:end)),'.--');
grid on;
legend('Standard','Modified');
xlabel('wavenumber, k');
ylabel('Magnitude of coefficient');
title(['l=' num2str(l) ', cond(V)=' num2str(cond(Vmod)), ', |d(1)|=' num2str(abs(dk(1)))]);

figure;
semilogy(ttj,err_std,'.');
hold on;
semilogy(ttj,err_mod,'.');
grid on;
xlim([0,2*pi]);
legend('Standard','Modified');
xlabel('t');
ylabel('Absolute error');
title('reconstruction error');

figure;
plot(ttj,fstd,'.');
hold on;
plot(ttj,fmod,'o');
plot(ttj,ffj);
xlim([0,2*pi]);
legend('Standard','Modified','Exact');
grid on;
xlabel('t');
ylabel('function value');
title('true function and recons');

alignfigs;

end
