function [mu1,mu3,mu5] = periodic_basis_integrals(r,kmax,n,K,E)
% PERIODIC_BASIS_INTEGRALS compute integrals mu of periodic basis functions
%
% [mu1 mu3,mu5] = PERIODIC_BASIS_INTEGRALS(r,kmax,n,K,E)
%   returns vectors of basis integrals of periodic basis functions with 
%   singularities of order p = 1/2, 3/2, 5/2, constructed via recurrence
%   relations
%
% INPUTS:
%   r    - scalar, 0 < r < 1, defined from root
%   kmax - maximum Fourier wavenumber (positive integer).
%   n    - number of equispaced discretization points in periodic direction
%   K    - complete elliptic integral of the 1st kind, eval at r^2, K(r^2)
%   E    - complete elliptic integral of the 2nd kind, eval at r^2, E(r^2)
%
% OUTPUTS:
%   mu1 - basis integrals for singularity p = 1/2
%   mu3 - basis integrals for singularity p = 3/2 (if requested)
%   mu5 - basis integrals for singularity p = 5/2 (if requested)
%
% NOTE: Recurrence relations from Lemma 2.3 in 
%   https://arxiv.org/pdf/2412.19575

klist = (0:kmax)';
M = kmax + 1;

% p=1/2 (col 1), 3/2 (col 2), 5/2 (col 3)
P = zeros(M, 3);

% constants
r2 = r^2;
rm1 = 1 - r;
rp1 = 1 + r;
rp1_4 = rp1^4;
inv_r = 1 / r;

% initial values
% p = 1/2
P(1,1) = 2 * K;
P(2,1) = 2 * inv_r * (-E + K);

% p = 3/2
P(1,2) = 2 / rp1 * (2 / rp1 * E - rm1 * K);

% p = 5/2
P(1,3) = 2 / (3 * rp1_4) * (8 * (1 + r2) * E - rm1 * rp1 * (5 + 3*r2) * K);

% recurrence for p=1/2
for i = 3:M
    k = klist(i);
    denom = 2*k - 1;
    P(i,1) = (1 + r2) * 2 * (k - 1) / denom * inv_r * P(i-1,1) ...
           - (2*k - 3) / denom * P(i-2,1);
end

if kmax ~= n / 2
    mu1 = [conj(P(M:-1:2,1)); P(:,1)];
else
    mu1 = [conj(P(M:-1:2,1)); P(1:end-1,1)];
end

if nargout == 1
    return;
end

% recurrence for p=3/2 and 5/2 (vectorized over p)
plist = [3/2, 5/2];
rm1_2 = rm1^2;
a = (1 + r2) / (2 * r);
b = rm1_2 / (2 * r);
for i = 2:M
    k = klist(i);
    c = (plist + (k - 2)) ./ (plist - 1);
    P(i,2:3) = a * P(i-1,2:3) - b * c .* P(i-1,1:2);
end

if kmax ~= n / 2
    mu3 = [conj(P(M-1:2,2)); P(:,2)];
    mu5 = [conj(P(M-1:2,3)); P(:,3)];
else
    mu3 = [conj(P(M:-1:2,2)); P(1:end-1,2)];
    mu5 = [conj(P(M:-1:2,3)); P(1:end-1,3)];
end
end
