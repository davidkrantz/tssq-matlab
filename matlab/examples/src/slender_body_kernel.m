function [u1, u2, u3] = slender_body_kernel(r1, r2, r3, f1, f2, f3, e)
    % Stokeslet + e^2/2*doublet
    R = sqrt(r1.^2 + r2.^2 + r3.^2);
    R3 = R.*R.*R;
    R5 = R3.*R.*R;
    rdotf = r1.*f1 + r2.*f2 + r3.*f3;
    u1 = f1./R + r1.*rdotf./R3 + e^2/2*(f1./R3 - 3*r1.*rdotf./R5);
    u2 = f2./R + r2.*rdotf./R3 + e^2/2*(f2./R3 - 3*r2.*rdotf./R5);
    u3 = f3./R + r3.*rdotf./R3 + e^2/2*(f3./R3 - 3*r3.*rdotf./R5);    
end
