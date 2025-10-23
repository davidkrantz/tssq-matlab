function [specquad1,specquad2,specquad3,stats] = interpolatory_quadrature_fourier(...
    xj, yj, zj, sj, tj, wj, f1j, f2j, f3j, X, Y, Z, nquad, slender_eps, use_mod, bcrit, tol, corrR3, corrR5)
% Based on https://github.com/ludvigak/linequad/blob/master/matlab/examples-paper/demo_long_fiber.m

    % Setup
    [specquad1,specquad2,specquad3] = deal(zeros(size(X)));
    maintic = tic();
    time_near = 0;
    time_far = 0;

    % Compute quadrature weights
    [all_w1,all_w3,all_w5,all_B03,all_B05,all_bary_weights,specquad_needed,std_needed,mod_needed,all_roots,specquad_timings] = ...
        near_weights_fourier(tj, wj, xj, yj, zj, X, Y, Z, tol, ...
        'use_mod', use_mod, 'bcrit', bcrit, 'corrR3', corrR3, 'corrR5', corrR5);

    % Stats
    kerevals_near = nquad*sum(specquad_needed(:));
    time_weights = specquad_timings.time_weights;
    time_rootfinding = specquad_timings.time_rootfinding;
    nbr_std_pts = sum(std_needed);
    nbr_mod_pts = sum(mod_needed);

    % Evaluation point loop
    for i=1:numel(X)
        Xi = X(i);
        Yi = Y(i);
        Zi = Z(i);
        if specquad_needed(i)
            atic = tic();
            [u1R1,u1R3,u1R5,u2R1,u2R3,u2R5,u3R1,u3R3,u3R5] = deal(zeros(1,nquad));
            for k=1:nquad
                r1k = xj(k)-Xi;
                r2k = yj(k)-Yi;
                r3k = zj(k)-Zi;
                [u1R1(k), u1R3(k), u1R5(k), u2R1(k), u2R3(k), u2R5(k), u3R1(k), u3R3(k), u3R5(k)] ...
                    = slender_body_kernel_split(r1k, r2k, r3k, f1j(k), f2j(k), f3j(k), ...
                                                slender_eps);
            end
            % extract precomputed data
            w1 = all_w1(:,i);
            w3 = all_w3(:,i);
            w5 = all_w5(:,i);
            t0 = all_roots(i); % with real part in [0,2*pi)
            a = real(t0);
            % integrands
            g1R1 = sj.*u1R1;
            g1R3 = sj.*u1R3;
            g1R5 = sj.*u1R5;
            g2R1 = sj.*u2R1;
            g2R3 = sj.*u2R3;
            g2R5 = sj.*u2R5;
            g3R1 = sj.*u3R1;
            g3R3 = sj.*u3R3;
            g3R5 = sj.*u3R5;
            % 1/R using adjoint standard SSQ
            I1R1 = g1R1*w1;
            I2R1 = g2R1*w1;
            I3R1 = g3R1*w1;
            if mod_needed(i)
                btic = tic();
                % adjoint TSSQ
                wbary = all_bary_weights(:,i);
                B03 = all_B03(i); % const coeff basis integrals
                B05 = all_B05(i);
                % evaluate integrand in t=a
                xa = xj*wbary; ya = yj*wbary; za = zj*wbary;
                sa = sj*wbary;
                f1a = f1j*wbary; f2a = f2j*wbary; f3a = f3j*wbary;
                r1a = xa-Xi; r2a = ya-Yi; r3a = za-Zi;
                tdista = abs(exp(1i*a)-exp(1i*t0));
                tdista3 = tdista.^3;
                tdista5 = tdista3.*tdista.*tdista;
                [~, u1R3a, u1R5a, ~, u2R3a, u2R5a, ~, u3R3a, u3R5a] ...
                    = slender_body_kernel_split(r1a, r2a, r3a, f1a, f2a, f3a, slender_eps);
                % new adjoint method for 1/R^3 and 1/R^5
                w03 = B03*tdista3*sa;
                w05 = B05*tdista5*sa;
                time_weights = time_weights + toc(btic);
                I1R3 = g1R3*w3;
                I2R3 = g2R3*w3;
                I3R3 = g3R3*w3;
                I1R5 = g1R5*w5;
                I2R5 = g2R5*w5;
                I3R5 = g3R5*w5;
                if corrR3
                    I1R3 = I1R3 + w03*u1R3a;
                    I2R3 = I2R3 + w03*u2R3a;
                    I3R3 = I3R3 + w03*u3R3a;
                end
                if corrR5
                    I1R5 = I1R5 + w05*u1R5a;
                    I2R5 = I2R5 + w05*u2R5a;
                    I3R5 = I3R5 + w05*u3R5a;
                end
            else
                % adjoint standard SSQ for 1/R^3 and 1/R^5
                I1R3 = g1R3*w3;
                I2R3 = g2R3*w3;
                I3R3 = g3R3*w3;
                I1R5 = g1R5*w5;
                I2R5 = g2R5*w5;
                I3R5 = g3R5*w5;
            end
            % collect components
            q1 = I1R1 + I1R3 + I1R5;
            q2 = I2R1 + I2R3 + I2R5;
            q3 = I3R1 + I3R3 + I3R5;
            time_near =  time_near + toc(atic);
        else            
            atic = tic();
            [q1, q2, q3] = quadsum(xj, yj, zj, sj, wj, f1j, f2j, f3j, ...
                                   Xi, Yi, Zi, nquad, slender_eps);
            time_far = time_far + toc(atic);
        end
        specquad1(i) = specquad1(i) + q1;
        specquad2(i) = specquad2(i) + q2;
        specquad3(i) = specquad3(i) + q3;
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
