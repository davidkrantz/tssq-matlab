function [w1, w3, w5, Ptilde03, Ptilde05, wbary] = rsqrt_pow_weights_modified_monomial(tj, troot, bary_wts)

n = numel(tj);
a = real(troot);

use_bjorck_pereyra = true;

p1 = rsqrt_pow_integrals(troot, n); % standard basis
[~, ptilde3, ptilde5] = rsqrt_pow_integrals_shift(troot, n); % translated basis

Ptilde03 = ptilde3(1);
Ptilde05 = ptilde5(1);
ptildek3 = [0; ptilde3(2:end)];
ptildek5 = [0; ptilde5(2:end)];

% Barycentric weights for interpolation at t=real(troot)
tmp = bary_wts./(a-tj);
wbary = tmp/sum(tmp);

% Standard and translated monomial basis
basis_std = tj;
basis_mod = tj-a;

% Compute "modified quadrature weights"
tdist = abs(tj-troot);
if n < 33 && use_bjorck_pereyra
    % Using Bjorck-Pereyra to solve Vandermonde system
    w1 = pvand(basis_std, p1) .* tdist; % O(n^2)
    w3 = pvand(basis_mod, ptildek3) .* tdist.^3;
    w5 = pvand(basis_mod, ptildek5) .* tdist.^5;
else
    % Direct solve with Vandermonde matrix is more stable for n>32
    % Still, n>40 is not a good idea    
    Astd = ones(n);
    Amod = ones(n);
    for j=2:n
        Astd(j,:) = Astd(j-1,:).*basis_std.'; % build transpose, O(n^2)
        Amod(j,:) = Amod(j-1,:).*basis_mod.';
    end
    warning('off', 'MATLAB:nearlySingularMatrix')
    W1 = Astd \ p1;
    W35 = Amod \ [ptildek3 ptildek5]; % O(n^3)
    W = [W1, W35];
    warning('on', 'MATLAB:nearlySingularMatrix')
    w1 = W(:,1) .* tdist;
    w3 = W(:,2) .* tdist.^3;
    w5 = W(:,3) .* tdist.^5;
end
