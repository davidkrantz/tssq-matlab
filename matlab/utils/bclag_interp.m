function ff = bclag_interp(f,x,w,xx)
% ff = bclag_interp(f,x,w,xx)
% 
% Barycentric Lagrange Interpolation.
% Berrut, J.-P., & Trefethen, L. N. (2004).
% SIAM Review, 46(3), 501–517. doi:10.1137/S0036144502417715
%
% NOTE:
% Paper reads exact(xdiff==0) = 1, but it should be exact(xdiff==0) = j

numer = zeros(size(xx));
denom = zeros(size(xx));
exact = zeros(size(xx));
N = numel(x);
for j = 1:N
    xdiff = xx-x(j);
    temp = w(j)./xdiff;
    numer = numer + temp*f(j);
    denom = denom + temp;
    exact(xdiff==0) = j;
end
ff = numer./denom;
jj = find(exact);
ff(jj) = f(exact(jj));
end
