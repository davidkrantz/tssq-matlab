function [quad1, quad2, quad3, stats] = adaptive_quadrature(xj, yj, zj, sj, wj, f1j, f2j, f3j, ...
                                                     X, Y, Z, nquad, slender_eps, Hlim)
    npan = numel(xj)/nquad;
    [quad1, quad2, quad3] = deal(zeros(size(X)));
    [tgl,wgl] = legendre.gauss(nquad);
    tlo = (tgl-1)/2;
    thi = (tgl+1)/2;
    Blo = bclag_interp_matrix(tgl, tlo);
    Bhi = bclag_interp_matrix(tgl, thi);
    kerevals_near = 0;
    maintic = tic();
    time_sub = 0;
    time_interp = 0;
    time_ker_near = 0;
    time_far = 0;
    for j=1:npan
        % Load panel
        idx = (1:nquad) + nquad*(j-1);
        xjpan = xj(idx);
        yjpan = yj(idx);
        zjpan = zj(idx);
        sjpan = sj(idx);
        wjpan = wj(idx);        
        hpan = sjpan*wjpan;
        f1jpan = f1j(idx);
        f2jpan = f2j(idx);
        f3jpan = f3j(idx);    
        for i=1:numel(X)
            % Compute for each point
            Xi = X(i);
            Yi = Y(i);
            Zi = Z(i);
            r1 = xjpan-Xi;
            r2 = yjpan-Yi;
            r3 = zjpan-Zi;        
            Rmin = sqrt( min(r1.^2 + r2.^2 + r3.^2) );
            Rlim = Hlim*hpan;
            if Rmin < Rlim
                % Subdivide until target no closer than H*h from subpanel
                atic = tic();
                [tnew, wnew] = subdivide(xjpan(:), yjpan(:), zjpan(:), Xi, Yi, Zi, ...
                                         Blo, Bhi, tgl, wgl, Rlim);
                time_sub = time_sub + toc(atic);
                atic = tic();                
                B = bclag_interp_matrix_mex(tgl, tnew)';                
                nnew = numel(tnew);
                % This is faster than multiplying with 16x7 block and the unpacking
                xtmp = xjpan*B;
                ytmp = yjpan*B;
                ztmp = zjpan*B;
                f1tmp = f1jpan*B;
                f2tmp = f2jpan*B;
                f3tmp = f3jpan*B;
                stmp = sjpan*B;
                wnew = wnew*wjpan(1)/wgl(1);
                time_interp = time_interp + toc(atic);                                
                atic = tic();                
                [q1, q2, q3] = quadsum(xtmp, ytmp, ztmp, stmp, wnew, f1tmp, f2tmp, f3tmp, ...
                                       Xi, Yi, Zi, nnew, slender_eps);
                time_ker_near = time_ker_near + toc(atic);                                                
                kerevals_near = kerevals_near + nnew;
            else
                atic = tic();
                [q1, q2, q3] = quadsum(xjpan, yjpan, zjpan, sjpan, wjpan, f1jpan, f2jpan, f3jpan, ...
                                       Xi, Yi, Zi, nquad, slender_eps);
                time_far = time_far + toc(atic);                
            end                    
            quad1(i) =  quad1(i) + q1;
            quad2(i) =  quad2(i) + q2;
            quad3(i) =  quad3(i) + q3;            
        end
    end
    toc(maintic)
    fprintf('Near field kernel evals: %e\n', kerevals_near);
    fprintf('Near field kernel evals time: %f\n', time_ker_near)
    fprintf('Interpolation time: %f\n', time_interp)
    fprintf('Subdivision time: %f\n', time_sub)
    fprintf('Far field time: %f\n', time_far)
    stats.kerevals_near = kerevals_near;
    stats.time_interp = time_interp;    
    stats.time_ker_near = time_ker_near;    
end
