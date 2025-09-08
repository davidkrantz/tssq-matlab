function [mu1,mu3,mu5] = periodic_basis_integrals(r,kmax,n,K,E,zero_order)
% PERIODIC_BASIS_INTEGRALS compute integrals mu of periodic basis functions
%
% [mu1 mu3,mu5] = PERIODIC_BASIS_INTEGRALS(r,kmax,n,K,E)
%   returns vectors of basis integrals of periodic basis functions with 
%   singularities of order p = 1/2, 3/2, 5/2, constructed via recurrence
%   relations
%
% INPUTS:
%   r          - scalar, 0 < r < 1, defined from root
%   kmax       - maximum Fourier wavenumber (positive integer).
%   n          - number of equispaced disc points in periodic direction
%   K          - complete elliptic integral of the 1st kind, eval at r^2, K(r^2)
%   E          - complete elliptic integral of the 2nd kind, eval at r^2, E(r^2)
%   zero_order - (optional, default: true) if true, output is zero-centered
%                [-kmax:-1, 0:kmax], else output is [0:kmax]
%
% OUTPUTS:
%   mu1 - basis integrals for singularity p = 1/2
%   mu3 - basis integrals for singularity p = 3/2
%   mu5 - basis integrals for singularity p = 5/2
%
% AUTHOR: David Krantz (davkra@kth.se)
%
% NOTE: Recurrence relations from Lemma 2.3 in 
%   https://arxiv.org/pdf/2412.19575

if nargin < 6
    zero_order = true;
end

% constants
M = kmax + 1;
r2 = r^2;
rm1 = 1 - r;
rp1 = 1 + r;
rp1_4 = rp1^4;
inv_r = 1 / r;

% allocate non-negative frequency vectors
mu1pos = zeros(M, 1);
mu3pos = zeros(M, 1);
mu5pos = zeros(M, 1);

% initial values
% p = 1/2
mu1pos(1) = 2*K;
mu1pos(2) = 2*inv_r*(-E+K);

% p = 3/2
mu3pos(1) = 2/rp1*(2/rp1*E-rm1*K);

% p = 5/2
mu5pos(1) = 2/(3*rp1_4)*(8*(1+r2)*E-rm1*rp1*(5+3*r2)*K);

% recurrence for p=1/2 (m=1)
for i = 3:M
    k = i - 1;
    denom = 2*k - 1;
    mu1pos(i) = (1 + r2) * 2 * (k - 1) / denom * inv_r * mu1pos(i-1) ...
           - (2*k - 3) / denom * mu1pos(i-2);
end

% recurrence for p=3/2 (m=3) and 5/2 (m=5)
rm1_2 = rm1^2;
a = (1 + r2) / (2 * r);
b = rm1_2 / (2 * r);
for i = 2:M
    k = i - 1;
    c3 = (3/2 + (k - 2)) ./ (3/2 - 1);
    c5 = (5/2 + (k - 2)) ./ (5/2 - 1);
    mu3pos(i) = a*mu3pos(i-1) - b*c3.*mu1pos(i-1);
    mu5pos(i) = a*mu5pos(i-1) - b*c5.*mu3pos(i-1);
end

if zero_order % zero-frequency order
    if kmax ~= n / 2
        mu1 = [conj(mu1pos(M:-1:2)); mu1pos];
        mu3 = [conj(mu3pos(M:-1:2)); mu3pos];
        mu5 = [conj(mu5pos(M:-1:2)); mu5pos];
    else
        mu1 = [conj(mu1pos(M:-1:2)); mu1pos(1:end-1)];
        mu3 = [conj(mu3pos(M:-1:2)); mu3pos(1:end-1)];
        mu5 = [conj(mu5pos(M:-1:2)); mu5pos(1:end-1)];
    end
else
    mu1 = mu1pos;
    mu3 = mu3pos;
    mu5 = mu5pos;
end
end
