function [mu1,mu3,mu5] = periodic_basis_integrals(r,kmax,nph,K,E)
% quite naively implemented, should be rewritten
klist = (0:kmax)';
M = length(klist);
P = zeros(M,3);

% f = @(phi,r,k,m) cos(k*phi)./(1-2*r*cos(phi)+r^2).^(m/2);
% warning off;
% for i = 1:M
%     k = klist(i);
%     P(i,3) = integral(@(phi) f(phi,r,k,m),0,pi,'AbsTol',1e-14,'RelTol',1e-14);
% end
% warning on;
% muk = P(:,3);
% muk = (1-r).^(m-1)*muk;
% mu = [flipud(conj(muk(2:end))); muk];
% return;

% Check if we evaluate in root
if r ~= 1
    % p = 1/2
    P(1,1) = 2*K;
    P(2,1) = 2/r*(-E+K);
    % p = 3/2
    P(1,2) = 2/(1+r)*(2/(1+r)*E-(1-r)*K);
    % p = 5/2
    P(1,3) = 2/(3*(1+r)^4)*(8*(1+r^2)*E-(1-r)*(1+r)*(5+3*r^2)*K);
else
    % K(1) is not defined, compute instead as K(1-eps)
    K = log(4/sqrt(1-(1-eps)^2));
    % p = 1/2
    P(1,1) = 2*K;
    P(2,1) = 2/r*(-E+K);
    % p = 3/2
    P(1,2) = 1;
    % p = 5/2
    P(1,3) = 2/3;
end

% Recurrence for p = 1/2
for i = 3:M
    k = klist(i);
    P(i,1) = (1+r^2)*2*(k-1)/(2*k-1)/r*P(i-1,1) - (2*k-3)/(2*k-1)*P(i-2,1);
end

muk = P(:,1);
if kmax ~= nph/2
    mu1 = [flipud(conj(muk(2:end))); muk];
else
    mu1 = [flipud(conj(muk(2:end))); muk(1:end-1)];
end

if nargout == 1
    return;
end

% Recurrence for p = 3/2, 5/2
plist = [3/2;5/2];
for j = 1:2
    ptmp = plist(j);
    for i = 2:M
        k = klist(i);
        P(i,j+1) = (1+r^2)/2/r*P(i-1,j+1)-(1-r)^2*(ptmp+k-2)/(ptmp-1)/2/r*P(i-1,j);
    end
    if nargout > 1
        muk = P(:,2);
        if kmax ~= nph/2
            mu3 = [conj(muk(M:-1:2)); muk];
        else
            mu3 = [conj(muk(M:-1:2)); muk(1:end-1)];
        end
    end
    if nargout == 2
        return;
    end
end

muk = P(:,3);
if kmax ~= nph/2
    mu5 = [conj(muk(M:-1:2)); muk];
else
    mu5 = [conj(muk(M:-1:2)); muk(1:end-1)];
end
end
