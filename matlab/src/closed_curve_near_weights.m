function [w1,w3,w5,mu1,mu3,mu5,specquad_needed,all_roots] = ...
    closed_curve_near_weights(tj, wj, xj, yj, zj, X, Y, Z, tol)
% CLOSED_CURVE_NEAR_WEIGHTS computes near evaluation quadrature weights for
%   kernel 1/R^m, m=1,3,5.
%
% [w1,w3,w5,mu1,mu3,mu5,specquad_needed,all_roots] = 
%   closed_curve_near_weights(tj, wj, xj, yj, zj, X, Y, Z, tol) returns the
%   target-specific weights, Fourier basis integrals and complex-valued
%   roots.
%
% INPUTS:
%   tj         - equdistant nodes in [0,2*pi)
%   wj         - quadrature weights accompanying tj
%   xj, yj, zj - coordinates on closed curve, corresponding to tj
%   X, Y, Z    - list of target points
%   tol        - desired error tolerance
%
% OUTPUTS:
%   w1, w3, w5      - matrix of size (numel(tj),numel(X)) with quadrature
%                     weights, replacing wj if needed
%   mu1, mu3, mu5   - Fourier basis integrals, zero-frequency shifted
%   specquad_needed - list of bool values indicating if special quadrature
%                     weights were computed
%   all_roots       - complex-valued root of R^2 for each target point
%
% AUTHOR: David Krantz (davkra@kth.se), May 2025

% reshape all inputs
tj = tj(:);
wj = wj(:);
xj = xj(:);
yj = yj(:);
zj = zj(:);

ntar = numel(X); % number of target points
n = length(xj); % number of discretization nodes

% Fourier expansion of discretization
xcoeff = fftshift(fft(xj))/n;
ycoeff = fftshift(fft(yj))/n;
zcoeff = fftshift(fft(zj))/n;
coeffs = [xcoeff,ycoeff,zcoeff]; % (n×3)
k = get_k_vec(n,2*pi).'; % wavenumber vector for FFT
kmax = max(abs(k)); % maximum wavenumber

% remove negligible Fourier modes, improves stability of root-finder for
% large n
max_vals = max(abs(coeffs), [], 1); % [max_x, max_y, max_z]
threshold = 1e-14;

% logical index for modes where all components are relatively small
remove_mask = ...
    abs(coeffs(:,1)) < threshold*max_vals(1) & ...
    abs(coeffs(:,2)) < threshold*max_vals(2) & ...
    abs(coeffs(:,3)) < threshold*max_vals(3);

% "truncate", i.e. remove modes
coeffstrunc = coeffs;
ktrunc = k;
coeffstrunc(remove_mask, :) = [];
ktrunc(remove_mask) = [];

% index and distance to closest grid node
kd = KDTree([xj yj zj]);
[tmpidx,dist] = kd.nn([X.' Y.' Z.']);

% first coarse filter: consider only target points with dist<cutdist
ds = (2*pi)/n; % grid spacing
cutdist = 8*ds; % cutoff distance
idxclose = dist<cutdist;
idxclosesest = tmpidx(idxclose);

% use closest point on curve as real part of initial guess
all_tinits = zeros(ntar,1);
all_tinits(idxclose) = tj(idxclosesest) + dist(idxclose)*1i; % add small imaginary part

% run root-finder
all_roots = complex(zeros(ntar,1));
rootfinder_converged = false(ntar,1);
target_idx = (1:ntar);
target_idx = target_idx(idxclose.');
for i = target_idx
    [troot,converged] = complex_fourier_newton(coeffstrunc,ktrunc,[X(i) Y(i) Z(i)],all_tinits(i));
    all_roots(i) = troot;
    rootfinder_converged(i) = converged;
end

% second filter: check which points need special quadrature based
% convergence of root-finder and imaginary part of root
all_rho = exp(-n*abs(imag(all_roots)));
specquad_needed = rootfinder_converged & (all_rho > tol);

% default is to return regular weights
[w1,w3,w5] = deal(repmat(wj,1,ntar));

% compute special weights where needed
[mu1,mu3,mu5] = deal(zeros(n,ntar));
for i = 1:ntar
    if specquad_needed(i)
        [w1(:,i),w3(:,i),w5(:,i),mu1(:,i),mu3(:,i),mu5(:,i)] = ...
            rsqrt_pow_weights_periodic(tj,kmax,all_roots(i));
    end
end
end
