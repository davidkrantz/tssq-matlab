function [tnew, wnew] = subdivide(xj, yj, zj, x0, y0, z0, Blo, Bhi, tgl, wgl, Rlim)
% Recursively subdivide panel until target point far away from all subpanels
    n = numel(xj);
    r1 = xj-x0;
    r2 = yj-y0;
    r3 = zj-z0;
    Rmin = sqrt( min(r1.^2 + r2.^2 + r3.^2) );
    if Rmin < Rlim
        dt = sum(wgl);
        tmid = (tgl(n)+tgl(1))/2;        
        tlo = (tgl+tmid-dt/2)/2;
        thi = (tgl+tmid+dt/2)/2;
        [tlo, wlo] = subdivide(Blo*xj, Blo*yj, Blo*zj, x0, y0, z0, ...
                               Blo, Bhi, tlo, wgl/2, Rlim/2);
        [thi, whi] = subdivide(Bhi*xj, Bhi*yj, Bhi*zj, x0, y0, z0, ...
                               Blo, Bhi, thi, wgl/2, Rlim/2);
        tnew = [tlo; thi];
        wnew = [wlo; whi];        
    else
        tnew = tgl;
        wnew = wgl;
    end    
end
