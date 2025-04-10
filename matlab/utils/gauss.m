function [x, w] = gauss(N)
% Gauss-Legendre quadrature 
    
% Golub-Welsch
    % k = 1:N-1;
    % b = sqrt(1./(4-1./k.^2));
    % A=zeros(N);
    % for n=2:N
    %     A(n,n-1)=b(n-1);
    %     A(n-1,n)=b(n-1);
    % end
    % [V, d] = eig(A, 'vector');
    % [x, idx] = sort(d);
    % w = 2*V(1, idx)'.^2;

% Newton method (seems to be better)
    [x, w] = lgwt(N, -1, 1);
end
