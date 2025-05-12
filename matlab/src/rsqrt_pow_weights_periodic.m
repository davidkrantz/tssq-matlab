function [w1,w3,w5,mu1,mu3,mu5] = rsqrt_pow_weights_periodic(tj,kmax,troot)
% RSQRT_POW_WEIGHTS_PERIODIC computes near evaluation quadrature weights
% for kernel 1/R^m, m=1,3,5.
%
% [w1,w3,w5,mu1,mu3,mu5] = rsqrt_pow_weights_periodic(tj,k,troot) returns 
%   the target-specific weights and Fourier basis integrals for wavenumbers
%   k and complex-valued root troot.
%
% INPUTS:
%   tj    - equdistant nodes in [0,2*pi)
%   kmax  - maximum wavenumber
%   troot - complex-valued root
%
% OUTPUTS:
%   w1, w3, w5    - matrix of size (numel(tj),numel(X)) with quadrature
%                   weights, replacing wj if needed
%   mu1, mu3, mu5 - Fourier basis integrals, zero-frequency shifted
%
% AUTHOR: David Krantz (davkra@kth.se), May 2025

n = numel(tj);

a = real(troot);
b = imag(troot);
r = exp(-abs(b)); % 0<r<1

[K,E] = ellipke(r^2);

% compute basis integrals for non-negative wavenumbers (0:kmax) through
% recurrences
[mu1pos,mu3pos,mu5pos] = periodic_basis_integrals(r,kmax,n,K,E,false);

% precompute flipped conjugate sections
mu1flip = conj(mu1pos(end:-1:2));
mu3flip = conj(mu3pos(end:-1:2));
mu5flip = conj(mu5pos(end:-1:2));

% decide index structure and handle symmetry
if mod(n,2) == 0
    mu1 = [mu1flip; mu1pos(1:end-1)];
    mu3 = [mu3flip; mu3pos(1:end-1)];
    mu5 = [mu5flip; mu5pos(1:end-1)];

    kstdord = [0:kmax-1, -kmax:-1]'; % standard Matlab orderings
    mu1stdord = [mu1pos(1:end-1); mu1flip];
    mu3stdord = [mu3pos(1:end-1); mu3flip];
    mu5stdord = [mu5pos(1:end-1); mu5flip];
else
    mu1 = [mu1flip; mu1pos];
    mu3 = [mu3flip; mu3pos];
    mu5 = [mu5flip; mu5pos];

    kstdord = [0:kmax, -kmax:-1]';
    mu1stdord = [mu1pos; mu1flip];
    mu3stdord = [mu3pos; mu3flip];
    mu5stdord = [mu5pos; mu5flip];
end

% compute near evaluation quadrature weights
ekstdord = exp(-1i*kstdord*a);

% correct scalings
p1 = 2*ekstdord.*mu1stdord;
p3 = 2*ekstdord.*mu3stdord./(1-r).^(3-1); % division is (1-r)^(m-1)
p5 = 2*ekstdord.*mu5stdord./(1-r).^(5-1);

W = ifft([p1 p3 p5],[],1,'symmetric'); % "adjoint" method

tdist = abs(exp(1i*tj)-exp(1i*troot)); % "regularizing factor"

w1 = W(:,1) .* tdist;
w3 = W(:,2) .* tdist.^3;
w5 = W(:,3) .* tdist.^5;

end