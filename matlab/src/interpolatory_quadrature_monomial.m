function [specquad1,specquad2,specquad3,cancellation_errest,corrneeded,stats] = interpolatory_quadrature_monomial(...
    xj, yj, zj, sj, tj, wj, f1j, f2j, f3j, X, Y, Z, x, y, z, s, nquad, edges, rho, UPSAMPLE, slender_eps, tol, use_bjorck_pereyra, use_mod, correct_R3, correct_R5)
% Based on https://github.com/ludvigak/linequad/blob/master/matlab/examples-paper/demo_long_fiber.m
    [specquad1,specquad2,specquad3] = deal(zeros(size(X)));
    cancellation_errest = nan*zeros(size(X));
    corrneeded = false(numel(X),1);
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
    wbary = bclag_interp_weights(tgl2);
    time_std_weights = 0;
    time_mod_weights = 0;
    time_mod_basis_int = 0;
    time_mod_vander_solve = 0;
    time_mod_sigma_interp_a = 0;
    time_mod_kernel_eval_a = 0;
    time_cancel_est = 0;
    kerevals_near = 0;
    maintic = tic();
    time_ker_near = 0;
    time_far = 0;
    for j=1:npan
        % Load panel
        idx = (1:nquad) + nquad*(j-1);
        xjpan = xj(idx);
        yjpan = yj(idx);
        zjpan = zj(idx);
        sjpan = sj(idx);
        tjpan = tj(idx);
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
        tjpan_up = tjpan.'*B';
        % current panel endpoints
        ta = edges(j);
        tb = edges(j+1);
        tsc = (tb-ta)/2;
        % Compute quadrature weights
        atic = tic();
        [all_w1, all_w3, all_w5, specquad_needed, all_roots] = line3_near_weights(tgl2, wgl2, xjpan_up, yjpan_up, zjpan_up, ...
                                                          X, Y, Z, rho);
        t0_all_roots = (ta+tb)/2+(tb-ta)/2*all_roots; % roots in [0,1]
        % Evaluation count
        time_std_weights = time_std_weights + toc(atic);
        kerevals_near = kerevals_near + nquad2*sum(specquad_needed(:));
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
                % integrands
                g1R1 = (sjpan_up.*u1R1).';
                g1R3 = (sjpan_up.*u1R3).';
                g1R5 = (sjpan_up.*u1R5).';
                g2R1 = (sjpan_up.*u2R1).';
                g2R3 = (sjpan_up.*u2R3).';
                g2R5 = (sjpan_up.*u2R5).';
                g3R1 = (sjpan_up.*u3R1).';
                g3R3 = (sjpan_up.*u3R3).';
                g3R5 = (sjpan_up.*u3R5).';
                % vectors to sum
                tmp1R1 = all_w1(:,i).*g1R1;
                tmp1R3 = all_w3(:,i).*g1R3;
                tmp1R5 = all_w5(:,i).*g1R5;
                tmp2R1 = all_w1(:,i).*g2R1;
                tmp2R3 = all_w3(:,i).*g2R3;
                tmp2R5 = all_w5(:,i).*g2R5;
                tmp3R1 = all_w1(:,i).*g3R1;
                tmp3R3 = all_w3(:,i).*g3R3;
                tmp3R5 = all_w5(:,i).*g3R5;
                I1R1 = sum(tmp1R1);
                I1R3 = sum(tmp1R3);
                I1R5 = sum(tmp1R5);
                I2R1 = sum(tmp2R1);
                I2R3 = sum(tmp2R3);
                I2R5 = sum(tmp2R5);
                I3R1 = sum(tmp3R1);
                I3R3 = sum(tmp3R3);
                I3R5 = sum(tmp3R5);
                I1R3sh = I1R3; I2R3sh = I2R3; I3R3sh = I3R3;
                I1R5sh = I1R5; I2R5sh = I2R5; I3R5sh = I3R5;
                scale_fac = sum(wjpan)/2;
                v0 = all_roots(i); % root having real part in [-1,1]
                toln = tol/npan;
                if abs(real(v0)) <= 1.1 && use_mod
                    btic = tic();
                    w3 = all_w3(:,i);
                    w5 = all_w5(:,i);
                    normw3 = norm(w3,inf);
                    normw5 = norm(w5,inf);
                    % estimate for 1/R^5
                    errestR35 = zeros(3,2);
                    errestR35(1,2) = cond_sum(normw5,g1R5,I1R5);
                    errestR35(2,2) = cond_sum(normw5,g2R5,I2R5);
                    errestR35(3,2) = cond_sum(normw5,g3R5,I3R5);
                    corrR5 = sum(errestR35(:,2) > toln) > 0;
                    if ~corrR5
                        % 1/R^5 ok, now check 1/R^3
                        errestR35(1,1) = cond_sum(normw3,g1R3,I1R3);
                        errestR35(2,1) = cond_sum(normw3,g2R3,I2R3);
                        errestR35(3,1) = cond_sum(normw3,g3R3,I3R3);
                        corrR3 = sum(errestR35(:,1) > toln) > 0;
                    else
                        % 1/R^5 bad, assume same corr needed for 1/R^3
                        corrR3 = true;
                        errestR35(:,1) = errestR35(:,2);
                    end
                    time_cancel_est = time_cancel_est + toc(btic);
                    cancellation_errest(i) = max([max(errestR35,[],'all'),cancellation_errest(i)]);
    
                    if corrR5 || corrR3
                        ctic = tic();
                        corrneeded(i) = true;
                        sh = real(v0);
                        dtic = tic();
                        if corrR5
                            [~, p3sh, p5sh] = rsqrt_pow_integrals_shift(v0,nquad2);
                        else
                            [~, p3sh] = rsqrt_pow_integrals_shift(v0,nquad2);
                        end
                        time_mod_basis_int = time_mod_basis_int + toc(dtic);
    
                        % root w real part [0,1] (actual root of gamma(t), t in [0,1])
                        t0 = t0_all_roots(i);
                        a = real(t0);
                        tdist = abs(tjpan_up.'-t0);
        
                        % store all relevant values 
                        GR35 = zeros(3*nquad2,2);
                        GR35(:,1) = [g1R3;g2R3;g3R3];
                        GR35(:,2) = [g1R5;g2R5;g3R5];
                        IR35 = zeros(3,2);
                        IR35(:,1) = [I1R3;I2R3;I3R3];
                        IR35(:,2) = [I1R5;I2R5;I3R5];
        
                        % eval layer dens at real(v0) in [-1,1]
                        if abs(sh) <= 1
                            etic = tic();
                            [f1pana,f2pana,f3pana] = bclag_interp_sigma(f1jpan_up,f2jpan_up,f3jpan_up,tgl2,wbary,sh);
                            time_mod_sigma_interp_a = time_mod_sigma_interp_a + toc(etic);
                            %f1pana = x(a); f2pana = y(a); f3pana = z(a);
                            % eval remaining integrand analytically
                            ftic = tic();
                            xa = x(a); ya = y(a); za = z(a);
                            r1a = xa-Xi; r2a = ya-Yi; r3a = za-Zi;
                            tdista = abs(a-t0);
                            sa = s(a);
                            [~, u1R3a, u1R5a, ~, u2R3a, u2R5a, ~, u3R3a, u3R5a] ...
                                = slender_body_kernel_split(r1a, r2a, r3a, f1pana, f2pana, f3pana, slender_eps);
                            uR35a = [u1R3a, u1R5a;
                                     u2R3a, u2R5a;
                                     u3R3a, u3R5a];
                            time_mod_kernel_eval_a = time_mod_kernel_eval_a + toc(ftic);
                        end

                        % precompute factorization of Vandermonde matrix
                        alpha = tgl2-sh;
                        if ~use_bjorck_pereyra
                            gtic = tic();
                            V = ones(nquad2,nquad2); for k=2:nquad2, V(:,k) = V(:,k-1).*alpha; end
                            [L,U,P] = lu(V);
                            %[Q,R] = qr(V);
                            time_mod_vander_solve = time_mod_vander_solve + toc(gtic);
                        end

                        time_mod_weights = time_mod_weights + toc(ctic);

                        tdistMat = zeros(nquad2,2);
                        tmp = tdist.^(2*1+1);
                        tdistMat(:,1) = tmp;
                        tdistMat(:,2) = tmp.*tdist.*tdist;
                        for ii = 1:3
                            for jj = 1:2
                                corr = errestR35(ii,jj) > toln;
                                if corr
                                    dtic = tic();
                                    g = GR35(((ii-1)*nquad2+1):ii*nquad2,jj);
                                    h = g.*tdistMat(:,jj); % to expand in shifted monomial basis
                                    htic = tic();
                                    if use_bjorck_pereyra
                                        dcoeff = dvand(alpha,h); % Björck-Pereyra
                                    else
                                        dcoeff = U\(L\(P*h));
                                        %dcoeff = R\(Q\h);
                                    end
                                    %warning('off','MATLAB:nearlySingularMatrix');
                                    %dcoeff = V\h;
                                    %warning('on','MATLAB:nearlySingularMatrix')
                                    time_mod_vander_solve = time_mod_vander_solve + toc(htic);
                                    if abs(sh) <= 1
                                        d1coeff = sa*uR35a(ii,jj)*tdista^(2*jj+1); % correction to d(1)
                                        dcoeff(1) = d1coeff;
                                    end
                                    time_mod_weights = time_mod_weights + toc(dtic);
                                    if jj == 1 && correct_R3
                                        IR35(ii,jj) = sum(dcoeff.*p3sh)/tsc^3;
                                    elseif jj == 2 && correct_R5
                                        IR35(ii,jj) = sum(dcoeff.*p5sh)/tsc^5;
                                    end
                                end
                            end
                        end
                        I1R3sh = IR35(1,1);
                        I2R3sh = IR35(2,1);
                        I3R3sh = IR35(3,1);
                        I1R5sh = IR35(1,2);
                        I2R5sh = IR35(2,2);
                        I3R5sh = IR35(3,2);
                    end
                end
                q1 = I1R1 + I1R3sh + I1R5sh;
                q2 = I2R1 + I2R3sh + I2R5sh;
                q3 = I3R1 + I3R3sh + I3R5sh;
                % Rescale (weights are for [-1,1])
                q1 = q1*scale_fac;
                q2 = q2*scale_fac;
                q3 = q3*scale_fac;
                time_ker_near =  time_ker_near + toc(atic);
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
    fprintf('Near field kernel evals time: %f\n', time_ker_near)
    fprintf('Time line3_near_weights: %f\n', time_std_weights)
    fprintf('Far field time: %f\n', time_far)
    toc(maintic)
    stats.kerevals_near = kerevals_near;
    stats.time_std_weights = time_std_weights;
    stats.time_ker_near = time_ker_near;
    stats.time_mod_weights = time_mod_weights;
    stats.time_cancel_est = time_cancel_est;
    stats.time_mod_basis_int = time_mod_basis_int;
    stats.time_mod_sigma_interp_a = time_mod_sigma_interp_a;
    stats.time_mod_kernel_eval_a = time_mod_kernel_eval_a;
    stats.time_mod_vander_solve = time_mod_vander_solve;
end
