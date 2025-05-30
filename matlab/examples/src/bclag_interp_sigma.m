function [ff1,ff2,ff3] = bclag_interp_sigma(f1,f2,f3,x,w,xx)
% [ff1,ff2,ff3] = bclag_interp(f1,f2,f3,x,w,xx)
% 
% Barycentric Lagrange Interpolation.
% Berrut, J.-P., & Trefethen, L. N. (2004).
% SIAM Review, 46(3), 501–517. doi:10.1137/S0036144502417715
%
% NOTE:
% Paper reads exact(xdiff==0) = 1, but it should be exact(xdiff==0) = j

numer1 = zeros(size(xx));
numer2 = zeros(size(xx));
numer3 = zeros(size(xx));
denom = zeros(size(xx));
exact = zeros(size(xx));
N = numel(x);
for j = 1:N
    xdiff = xx-x(j);
    temp = w(j)./xdiff;
    numer1 = numer1 + temp*f1(j);
    numer2 = numer2 + temp*f2(j);
    numer3 = numer3 + temp*f3(j);
    denom = denom + temp;
    exact(xdiff==0) = j;
end
ff1 = numer1./denom;
ff2 = numer2./denom;
ff3 = numer3./denom;
jj = find(exact);
ff1(jj) = f1(exact(jj));
ff2(jj) = f2(exact(jj));
ff3(jj) = f3(exact(jj));
end
