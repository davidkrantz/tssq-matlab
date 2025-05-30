function [tj, wj, npan, edges] = adaptive_panelization(s, nquad, tol)
    [tgl, wgl] = legendre.gauss(nquad); 
    tgl = (tgl+1)/2; wgl = wgl/2; % [0,1]
    L = legendre.matrix(nquad);
    function edges = recursion(ta, tb)
        assert(tb-ta > eps(), 'Recursion failed (tolerance too big?)')
        tj = ta + tgl*(tb-ta);
        sj = s(tj);
        c = abs(L*sj(:));
        res = max(c(end-1:end))/max(c);
        if res < tol
            % pass
            edges = [ta, tb];
        else 
            % subdivide
            tmid = (ta+tb)/2;
            edges = union( recursion(ta, tmid), recursion(tmid, tb) );
        end
    end    
    edges = recursion(0, 1);
    npan = numel(edges)-1;
    [tj, wj] = deal(zeros(npan*nquad, 1));
    for i=1:npan
        idx = (1:nquad) + nquad*(i-1);
        ta = edges(i);
        tb = edges(i+1);
        dt = tb-ta;
        tj(idx) = ta + tgl*dt;
        wj(idx) = wgl*dt;
    end
end
