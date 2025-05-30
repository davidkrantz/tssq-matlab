function [q1, q2, q3] = quadsum(xj, yj, zj, sj, wj, f1j, f2j, f3j, x, y, z, n, ...
                                slender_eps)
    q1 = 0; q2 = 0; q3 = 0;    
    for k=1:n
        r1 = xj(k)-x;
        r2 = yj(k)-y;
        r3 = zj(k)-z;                               
        [u1, u2, u3] = slender_body_kernel(r1, r2, r3, f1j(k), f2j(k), f3j(k), ...
                                           slender_eps);
        q1 = q1 + u1*sj(k)*wj(k);
        q2 = q2 + u2*sj(k)*wj(k);
        q3 = q3 + u3*sj(k)*wj(k);                    
    end
end
