function [x, w] = gauss_qp(N)
% Gauss-Legendre quadrature in quadruple precision
    digits(32);
    [x, w] = lgwt(N, vpa(-1), vpa(1));
end
