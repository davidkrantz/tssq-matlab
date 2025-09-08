function [tj,wj] = adaptive_global_discretization(s,tol)
n = 4;
res = tol+1;
itr = 0;
maxitr = 10;
while res > tol && itr < maxitr
    tj = linspace(0,2*pi,n+1).'; tj(end) = []; % periodic grid in [0,2*pi)
    sj = s(tj);
    c = fftshift(fft(sj))/n;
    res = max(abs(c(end-1:end))/max(abs(c)));
    n = 2*n;
    itr = itr+1;
end
wj = (tj(2)-tj(1))*ones(n/2,1);
end
