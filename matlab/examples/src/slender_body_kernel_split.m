function [u1R1, u1R3, u1R5, u2R1, u2R3, u2R5, u3R1, u3R3, u3R5] = ...
        slender_body_kernel_split(r1, r2, r3, f1, f2, f3, e)
    R = sqrt(r1.^2 + r2.^2 + r3.^2);
    R3 = R.*R.*R;
    R5 = R3.*R.*R;    
    rdotf = r1.*f1 + r2.*f2 + r3.*f3;
    u1R1 = f1./R;
    u2R1 = f2./R;
    u3R1 = f3./R;
    u1R3 = (r1.*rdotf + e^2/2*f1)./R3;
    u2R3 = (r2.*rdotf + e^2/2*f2)./R3;
    u3R3 = (r3.*rdotf + e^2/2*f3)./R3;    
    u1R5 = -3*e^2/2*r1.*rdotf./R5;
    u2R5 = -3*e^2/2*r2.*rdotf./R5;
    u3R5 = -3*e^2/2*r3.*rdotf./R5;        
end
