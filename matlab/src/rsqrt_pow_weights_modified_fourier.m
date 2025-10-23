function [w1,w3,w5,B03,B05,wbary] = rsqrt_pow_weights_modified_fourier(tj,kmax,troot)
% RSQRT_POW_WEIGHTS_MODIFIED_FOURIER computes near evaluation quadrature weights
%   for kernel 1/R^m, m=1,3,5. Uses a mix of standard and modified Fourier
%   basis.
%
% [w1,w3,w5,B03,B05,wbary] = rsqrt_pow_weights_periodic(tj,k,troot) returns 
%   the target-specific weights, constant coefficient Fourier basis
%   integrals and trigonometric barycentric interpolation weights.
%
% INPUTS:
%   tj    - equdistant nodes in [0,2*pi)
%   kmax  - maximum wavenumber
%   troot - complex-valued root
%
% OUTPUTS:
%   w1, w3, w5 - matrix of size (numel(tj),numel(X)) with quadrature
%                   weights, replacing wj if needed
%   B03, B05   - constant coeff Fourier basis integrals
%   wbary      - trig barycentric interp weights for interp at real(troot)
%
% AUTHOR: David Krantz (davkra@kth.se), May 2025

n = numel(tj);

a = real(troot);
b = imag(troot);
r = exp(-abs(b)); % 0<r<1

[K,E] = ellipke(r^2); % elliptic integrals of the 2nd kind

% compute basis integrals for non-negative wavenumbers (0:kmax) through
% recurrences
[mu1pos,mu3pos,mu5pos] = rsqrt_pow_integrals_fourier(r,kmax,n,K,E,false);

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
else
    mu1 = [mu1flip; mu1pos];
    mu3 = [mu3flip; mu3pos];
    mu5 = [mu5flip; mu5pos];

    kstdord = [0:kmax, -kmax:-1]';
    mu1stdord = [mu1pos; mu1flip];
end

% compute near evaluation quadrature weights
ekstdord = exp(-1i*kstdord*a);

% correct scalings
p1 = 2*ekstdord.*mu1stdord;

% modified basis integrals
kmod = get_k_vec(n-2,2*pi).';
emod = exp(1i*kmod*a);
p3osc = 0.5*(-(1-r)^2/(2*r*(1+r^2)).*mu3(2:end-1) + ((3/2+kmod-1)./(2*r).*mu1(2:end-1)-(3/2+kmod-2)./(1+r^2).*mu1(1:end-2))./(3/2-1));
p5osc = 0.5*(1-r)^(3-5) * (-(1-r)^2/(2*r*(1+r^2)).*mu5(2:end-1) + ((5/2+kmod-1)./(2*r).*mu3(2:end-1)-(5/2+kmod-2)./(1+r^2).*mu3(1:end-2))./(5/2-1));
if mod(n,2) == 0
    B03 = 2*mu3(n/2+1)/(1-r)^(3-1);
    B05 = 2*mu5(n/2+1)/(1-r)^(5-1);
else
    B03 = 2*mu3((n-1)/2+1)/(1-r)^(3-1);
    B05 = 2*mu5((n-1)/2+1)/(1-r)^(5-1);
end
Stilde3 = emod.*p3osc;
Stilde5 = emod.*p5osc;

% adjoint method for standard Fourier basis
W1 = ifft(p1,[],1,'symmetric');

% new adjoint method for modified Fourier basis
[W3,W5] = adjoint_modified_fourier(tj,a,Stilde3,Stilde5);

% barycentric interpolation weights
wbary = trig_barycentric_weights(tj,a);

tdist = abs(exp(1i*tj)-exp(1i*troot)); % "regularizing factor"
tdist3 = tdist.^3;
tdist5 = tdist3.*tdist.*tdist;

w1 = W1 .* tdist;
w3 = W3 .* tdist3;
w5 = W5 .* tdist5;

end