function [w1,w3,w5,Ptilde03,Ptilde05,wbary,specquad_needed,std_needed,mod_needed,all_roots,timings] = ...
    near_weights_monomial(tj, wj, xj, yj, zj, X, Y, Z, varargin)
% NEAR_WEIGHTS_MONOMIAL computes near evaluation quadrature weights for
%   kernel 1/R^m, m=1,3,5. Uses a mix of standard and modified (translated)
%   polynomial basis.
%
% [w1,w3,w5,B03,B05,bary_weights,specquad_needed,std_needed,mod_needed,all_roots,timings] = 
%   near_weights_monomial(tj, wj, xj, yj, zj, X, Y, Z, tol) returns the
%   target-specific weights, Fourier basis integrals and complex-valued
%   roots.
%
% INPUTS:
%   tj           - nodes in parametrization space [-1, 1] (probably Gauss-Legendre points)
%   wj           - quadrature weights accompanying tj (e.g. Gauss-Legendre weights)
%   xj, yj, zj   - coordinates on curve, corresponding to tj
%   X, Y, Z      - list of target points
%   rho          - critical Bernstein radius, special quadrature weights computed for roots inside it (default: 4^(16/n))
%   acrit, bcrit - use mod basis for |real(root)|<=acrit and |imag(root)|<=bcrit (default: acrit=1.05, bcrit=1e-2)
%   use_mod      - logical switch to use modified basis (default: true)
%
% OUTPUTS:
%   w1, w3, w5         - matrix of size (numel(tj),numel(X)) with quadrature weights
%   Ptilde03, Ptilde05 - constant coeff monomial basis integrals
%   specquad_needed    - list of bool values indicating if special quadrature weights were computed
%   std_needed         - list of bool values where standard basis is used
%   mod_needed         - list of bool values where modified basis is used
%   all_roots          - complex-valued root of R^2 for each target point
%   timings            - struct, timing of root-finding and quad weight gen
%
% AUTHOR: David Krantz (davkra@kth.se)

USE_MEX_CODE = false; % Switch off to debug codes

n = length(xj); % nbr Gauss-Legendre points

% Create an input parser
p = inputParser;

% Default values
default_rho = 4^(16/n);
default_acrit = 1.05;
default_bcrit = 1e300;
default_use_mod = true;

% Validation functions (optional but helpful)
validScalar = @(x) isnumeric(x) && isscalar(x);
validLogical = @(x) islogical(x) || islogical(logical(x));

% Add named parameters with defaults and validation
addParameter(p,'rho',default_rho,validScalar);
addParameter(p,'acrit',default_acrit,validScalar);
addParameter(p,'bcrit',default_bcrit,validScalar);
addParameter(p,'use_mod',default_use_mod,validLogical);

% Parse the inputs
parse(p, varargin{:});

% Extract results
rho = p.Results.rho;
acrit = p.Results.acrit;
bcrit = p.Results.bcrit;
use_mod = p.Results.use_mod;

% Reshape all inputs
tj = tj(:);
wj = wj(:);
xj = xj(:);
yj = yj(:);
zj = zj(:);

% Legendre expansions of discretization
Legmat = legendre.matrix(n); % O(n^2) per panel
xhat = Legmat*xj;
yhat = Legmat*yj;
zhat = Legmat*zj;

% Cap expansion at 16 coefficients.
% This seems to be more stable near boundary
if n > 16
    xhat = xhat(1:16);
    yhat = yhat(1:16);
    zhat = zhat(1:16);
end

% Rootfinding: initial guesses
root_tic = tic();
if USE_MEX_CODE
    all_tinits = rootfinder_initial_guess_mex(tj, xj, yj, zj, X, Y, Z); % O(n) per point
else
    all_tinits = complex(zeros(size(X)));
    % Standard guess
    for i=1:numel(X)    
        tinit = rootfinder_initial_guess(tj, xj, yj, zj, X(i), Y(i), Z(i)); % O(n) per point
        all_tinits(i) = tinit;
    end
end

% First filter: Don't check points whose initial guess is far away
cp = (bernstein_radius(all_tinits) < 1.5*rho); % cp: Compute Points

% Rootfinding: Run
all_roots = deal(complex(zeros(size(X))));
rootfinder_converged = false(size(X));
if USE_MEX_CODE
    [all_roots(cp), rootfinder_converged(cp)] = rootfinder_mex(xhat, yhat, zhat, X(cp), Y(cp), Z(cp), all_tinits(cp));
else
    for i=find(cp(:)')
        tinit = all_tinits(i);
        [troot, converged] = rootfinder(xhat, yhat, zhat,  X(i), Y(i), Z(i), tinit); % O(n) per point        
        all_roots(i) = troot;
        rootfinder_converged(i) = converged;    
    end
end
time_rootfinding = toc(root_tic);

% Check which points need special attention
all_bernstein_radii = bernstein_radius(all_roots);
specquad_needed = rootfinder_converged & (all_bernstein_radii < rho);

% Determine where we use standard and modified basis
if use_mod
    mod_needed = specquad_needed & ...
                 (abs(real(all_roots)) <= acrit) & ...
                 (abs(imag(all_roots)) <= bcrit);
    bary_wts = bclag_interp_weights(tj); % Barycentric Lagrange interp wts
else
    mod_needed = false(size(specquad_needed));
end
std_needed = specquad_needed & ~mod_needed;

% Init data
[Ptilde03,Ptilde05] = deal(zeros(numel(X),1));
wbary = zeros(n,numel(X));

% Default is to return regular weights
weights_tic = tic();
[w1, w3, w5] = deal(repmat(wj(:), 1, numel(X)));

% Compute special weights where needed
if USE_MEX_CODE
    [tmp1, tmp3, tmp5] = rsqrt_pow_weights_mex(tj, all_roots(std_needed));
    w1(:,std_needed) = tmp1;
    w3(:,std_needed) = tmp3;
    w5(:,std_needed) = tmp5;
end
for i = 1:numel(X)
    if mod_needed(i)
        [w1(:,i), w3(:,i), w5(:,i), Ptilde03(i), Ptilde05(i), wbary(:,i)] = ...
            rsqrt_pow_weights_modified_monomial(tj, all_roots(i), bary_wts);
    elseif std_needed(i) && ~USE_MEX_CODE
            [w1(:,i), w3(:,i), w5(:,i)] = rsqrt_pow_weights(tj, all_roots(i));
    end
end
time_weights = toc(weights_tic);

% Timings
timings.time_rootfinding = time_rootfinding;
timings.time_weights = time_weights;
