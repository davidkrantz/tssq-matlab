function a = taylor_from_monomial(c, t0)
% Computes Taylor expansion coefficients a of a polynomial f(t) = sum c_k * t^{k-1}
% centered at t0: f(t) ≈ sum a(m+1) * (t - t0)^m
n = length(c);
a = vpa(zeros(n, 1));

for m = 0:n-1
    sum_term = 0;
    for k = m+1:n
        coef = c(k) * factorial(k-1) / factorial(k-1 - m);
        sum_term = sum_term + coef * t0^(k - 1 - m);
    end
    a(m+1) = sum_term / factorial(m);
end
end
