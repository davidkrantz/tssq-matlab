function [w1,w3,w5,B03,B05,bary_weights,specquad_needed,std_needed,mod_needed,all_roots,timings] = ...
    near_weights_fourier(tj, wj, xj, yj, zj, X, Y, Z, tol, varargin)
% NEAR_WEIGHTS_FOURIER computes near evaluation quadrature weights for
%   kernel 1/R^m, m=1,3,5. Uses a mix of standard and modified Fourier
%   basis.
%
% [w1,w3,w5,B03,B05,bary_weights,specquad_needed,std_needed,mod_needed,all_roots,timings] = 
%   near_weights_fourier(tj, wj, xj, yj, zj, X, Y, Z, tol) returns the
%   target-specific weights, Fourier basis integrals and complex-valued
%   roots.
%
% INPUTS:
%   tj         - equdistant nodes in [0,2*pi)
%   wj         - quadrature weights accompanying tj
%   xj, yj, zj - coordinates on closed curve, corresponding to tj
%   X, Y, Z    - list of target points
%   tol        - desired error tolerance
%   bcrit      - switch from std to mod at |imag(root)| <= bcrit (default: 1e-2)
%   use_mod    - logical switch to use modified basis (default: true)
%   corrR3     - logical switch to correct 1/R^3 term (default: true)
%   corrR5     - logical switch to correct 1/R^3 term (default: true)
%
% OUTPUTS:
%   w1, w3, w5      - matrix of size (numel(tj),numel(X)) with quadrature
%                     weights
%   B03, B05        - constant coeff Fourier basis integrals
%   specquad_needed - list of bool values indicating if special quadrature
%                     weights were computed
%   std_needed      - list of bool values where standard basis is used
%   mod_needed      - list of bool values where modified basis is used
%   all_roots       - complex-valued root of R^2 for each target point
%   timings         - struct, timing of root-finding and quad weight gen
%
% AUTHOR: David Krantz (davkra@kth.se), May 2025

% create an input parser
p = inputParser;

% default values
default_bcrit = 1e-2;
default_use_mod = true;
default_corrR3 = true;
default_corrR5 = true;

% validation functions (optional but helpful)
validScalar = @(x) isnumeric(x) && isscalar(x);
validLogical = @(x) islogical(x) || islogical(logical(x));

% add named parameters with defaults and validation
addParameter(p,'bcrit',default_bcrit,validScalar);
addParameter(p,'use_mod',default_use_mod,validLogical);
addParameter(p,'corrR3',default_corrR3,validLogical);
addParameter(p,'corrR5',default_corrR5,validLogical);

% parse the inputs
parse(p, varargin{:});

% extract results
bcrit = p.Results.bcrit;
use_mod = p.Results.use_mod;
corrR3 = p.Results.corrR3;
corrR5 = p.Results.corrR5;

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
root_tic = tic();
all_roots = complex(zeros(ntar,1));
rootfinder_converged = false(ntar,1);
target_idx = (1:ntar);
target_idx = target_idx(idxclose.');
for i = target_idx
    [troot,converged] = complex_fourier_newton(coeffstrunc,ktrunc,[X(i) Y(i) Z(i)],all_tinits(i));
    all_roots(i) = troot;
    rootfinder_converged(i) = converged;
end
time_rootfinding = toc(root_tic);

% second filter: check which points need special quadrature based
% convergence of root-finder and imaginary part of root
all_rho = exp(-n*abs(imag(all_roots)));
specquad_needed = rootfinder_converged & (all_rho > tol);

% determine where to use standard and modified basis
if use_mod
    mod_needed = specquad_needed & (abs(imag(all_roots)) <= bcrit);
else
    mod_needed = false(size(specquad_needed));
end
std_needed = specquad_needed & ~mod_needed;

% init data
[B03,B05] = deal(zeros(ntar,1));
bary_weights = zeros(n,ntar);

% default is to return regular weights
weights_tic = tic();
[w1,w3,w5] = deal(repmat(wj,1,ntar));

% compute special weights where needed
for i = 1:ntar
    if mod_needed(i)
        [w1(:,i),w3(:,i),w5(:,i),B03(i),B05(i),bary_weights(:,i)] = ...
            rsqrt_pow_weights_modified_fourier(tj,kmax,all_roots(i));
    elseif std_needed(i)
        [w1(:,i),w3(:,i),w5(:,i)] = ...
            rsqrt_pow_weights_fourier(tj,kmax,all_roots(i));
    end
end
time_weights = toc(weights_tic);

% gather timings
timings.time_rootfinding = time_rootfinding;
timings.time_weights = time_weights;

% for testing purposes, allow to correct 1/R^3 or 1/R^5 only
if ~corrR3
    for i = 1:ntar
        if specquad_needed(i)
            [~,w3(:,i),~] = ...
                rsqrt_pow_weights_fourier(tj,kmax,all_roots(i));
        end
    end
end
if ~corrR5
    for i = 1:ntar
        if specquad_needed(i)
            [~,~,w5(:,i)] = ...
                rsqrt_pow_weights_fourier(tj,kmax,all_roots(i));
        end
    end
end
