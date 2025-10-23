function [specquad1,specquad2,specquad3,stats] = interpolatory_quadrature_monomial(...
    xj, yj, zj, sj, wj, f1j, f2j, f3j, X, Y, Z, nquad, edges, rho, acrit, bcrit, UPSAMPLE, slender_eps, use_mod)
% Based on https://github.com/ludvigak/linequad/blob/master/matlab/examples-paper/demo_long_fiber.m
    [specquad1,specquad2,specquad3] = deal(zeros(size(X)));
    npan = numel(xj)/nquad;    
    [tgl, wgl] = legendre.gauss(nquad);
    if UPSAMPLE
        nquad2 = 2*nquad;
        [tgl2, wgl2] = legendre.gauss(nquad2);
        B = bclag_interp_matrix(tgl, tgl2);
    else
        nquad2 = nquad;
        B = 1;
        tgl2 = tgl;
        wgl2 = wgl;
    end
    maintic = tic();
    time_weights = 0;
    time_rootfinding = 0;
    kerevals_near = 0;
    nbr_std_pts = 0;
    nbr_mod_pts = 0;
    time_near = 0;
    time_far = 0;
    % Loop over panels
    for j=1:npan
        % Load panel
        idx = (1:nquad) + nquad*(j-1);
        xjpan = xj(idx);
        yjpan = yj(idx);
        zjpan = zj(idx);
        sjpan = sj(idx);
        wjpan = wj(idx);
        f1jpan = f1j(idx);
        f2jpan = f2j(idx);
        f3jpan = f3j(idx);
        % Upsample panel
        xjpan_up  = xjpan *B';
        yjpan_up  = yjpan *B';
        zjpan_up  = zjpan *B';
        sjpan_up  = sjpan *B';
        f1jpan_up = f1jpan*B';
        f2jpan_up = f2jpan*B';
        f3jpan_up = f3jpan*B';
        % current panel endpoints
        ta = edges(j);
        tb = edges(j+1);
        tsc = (tb-ta)/2;
        tsc3 = tsc^3;
        tsc5 = tsc3*tsc*tsc;
        % Compute quadrature weights
        [all_w1, all_w3, all_w5, all_Ptilde03, all_Ptilde05, all_wbary_weights, specquad_needed, std_needed, mod_needed, all_roots, specquad_timings] = ...
            near_weights_monomial( ...
            tgl2, wgl2, xjpan_up, yjpan_up, zjpan_up, ...
            X, Y, Z, 'rho', rho, 'acrit', acrit, 'bcrit', bcrit, 'use_mod', use_mod);
        t0_all_roots = (ta+tb)/2+(tb-ta)/2*all_roots; % roots in [0,1]
        % Evaluation and timing count
        time_weights = time_weights + specquad_timings.time_weights;
        time_rootfinding = time_rootfinding + specquad_timings.time_rootfinding;
        kerevals_near = kerevals_near + nquad2*sum(specquad_needed(:));
        nbr_std_pts = nbr_std_pts + sum(std_needed);
        nbr_mod_pts = nbr_mod_pts + sum(mod_needed);
        % Evaluate each panel-to-point pair
        for i=1:numel(X)    
            Xi = X(i);
            Yi = Y(i);
            Zi = Z(i);
            if specquad_needed(i)
                atic = tic();
                [u1R1,u1R3,u1R5,u2R1,u2R3,u2R5,u3R1,u3R3,u3R5] = deal(zeros(1,nquad2));
                for k=1:nquad2
                    r1k = xjpan_up(k)-Xi;
                    r2k = yjpan_up(k)-Yi;
                    r3k = zjpan_up(k)-Zi;
                    [u1R1(k), u1R3(k), u1R5(k), u2R1(k), u2R3(k), u2R5(k), u3R1(k), u3R3(k), u3R5(k)] ...
                        = slender_body_kernel_split(r1k, r2k, r3k, f1jpan_up(k), f2jpan_up(k), f3jpan_up(k), ...
                                                    slender_eps);
                end
                % extract precomputed data
                w1 = all_w1(:,i);
                w3 = all_w3(:,i);
                w5 = all_w5(:,i);
                t0 = t0_all_roots(i); % root with real part in [0,1]
                a = real(t0);
                % integrands
                g1R1 = sjpan_up.*u1R1;
                g1R3 = sjpan_up.*u1R3;
                g1R5 = sjpan_up.*u1R5;
                g2R1 = sjpan_up.*u2R1;
                g2R3 = sjpan_up.*u2R3;
                g2R5 = sjpan_up.*u2R5;
                g3R1 = sjpan_up.*u3R1;
                g3R3 = sjpan_up.*u3R3;
                g3R5 = sjpan_up.*u3R5;
                % 1/R using adjoint standard SSQ
                I1R1 = g1R1*w1;
                I2R1 = g2R1*w1;
                I3R1 = g3R1*w1;
                if mod_needed(i)
                    btic = tic();
                    % adjoint TSSQ
                    Ptilde03 = all_Ptilde03(i);
                    Ptilde05 = all_Ptilde05(i);
                    wbary = all_wbary_weights(:,i);
                    % evaluate integrand in t=a
                    xa = xjpan_up*wbary; ya = yjpan_up*wbary; za = zjpan_up*wbary;
                    sa = sjpan_up*wbary;
                    f1pana = f1jpan_up*wbary; f2pana = f2jpan_up*wbary; f3pana = f3jpan_up*wbary;
                    r1a = xa-Xi; r2a = ya-Yi; r3a = za-Zi;
                    tdista = abs(a-t0);
                    tdista3 = tdista.^3;
                    tdista5 = tdista3.*tdista.*tdista;
                    [~, u1R3a, u1R5a, ~, u2R3a, u2R5a, ~, u3R3a, u3R5a] ...
                        = slender_body_kernel_split(r1a, r2a, r3a, f1pana, f2pana, f3pana, slender_eps);
                    % new adjoint method for 1/R^3 and 1/R^5
                    w03 = Ptilde03*tdista3*sa;
                    w05 = Ptilde05*tdista5*sa;
                    time_weights = time_weights + toc(btic);
                    I1R3 = g1R3*w3 + w03*u1R3a/tsc3;
                    I2R3 = g2R3*w3 + w03*u2R3a/tsc3;
                    I3R3 = g3R3*w3 + w03*u3R3a/tsc3;
                    I1R5 = g1R5*w5 + w05*u1R5a/tsc5;
                    I2R5 = g2R5*w5 + w05*u2R5a/tsc5;
                    I3R5 = g3R5*w5 + w05*u3R5a/tsc5;
                else
                    % adjoint standard SSQ for 1/R^3 and 1/R^5
                    I1R3 = g1R3*w3;
                    I2R3 = g2R3*w3;
                    I3R3 = g3R3*w3;
                    I1R5 = g1R5*w5;
                    I2R5 = g2R5*w5;
                    I3R5 = g3R5*w5;
                end
                scale_fac = sum(wjpan)/2;
                q1 = I1R1 + I1R3 + I1R5;
                q2 = I2R1 + I2R3 + I2R5;
                q3 = I3R1 + I3R3 + I3R5;
                % Rescale (weights are for [-1,1])
                q1 = q1*scale_fac;
                q2 = q2*scale_fac;
                q3 = q3*scale_fac;
                time_near =  time_near + toc(atic);
            else            
                atic = tic();
                [q1, q2, q3] = quadsum(xjpan, yjpan, zjpan, sjpan, wjpan, f1jpan, f2jpan, f3jpan, ...
                                       Xi, Yi, Zi, nquad, slender_eps);
                time_far = time_far + toc(atic);
            end
            specquad1(i) = specquad1(i) + q1;
            specquad2(i) = specquad2(i) + q2;
            specquad3(i) = specquad3(i) + q3;
        end
    end
    fprintf('Near field kernel evals: %e\n', kerevals_near);
    fprintf('Near field time: %f\n', time_near)
    fprintf('Time rootfinding: %f\n', time_rootfinding)
    fprintf('Time weights: %f\n', time_weights)
    fprintf('Far field time: %f\n', time_far)
    toc(maintic)
    stats.kerevals_near = kerevals_near;
    stats.time_weights = time_weights;
    stats.time_near = time_near;
    stats.time_far = time_far;
    stats.time_rootfinding = time_rootfinding;
    stats.nbr_std_pts = nbr_std_pts;
    stats.nbr_mod_pts = nbr_mod_pts;
end
