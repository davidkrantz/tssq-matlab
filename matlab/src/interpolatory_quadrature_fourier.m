function [specquad1,specquad2,specquad3,cancellation_errest,corrneeded,stats] = interpolatory_quadrature_fourier(...
    xj, yj, zj, sj, tj, wj, f1j, f2j, f3j, X, Y, Z, x, y, z, s, nquad, slender_eps, tol, use_mod, correct_R3, correct_R5)
% Based on https://github.com/ludvigak/linequad/blob/master/matlab/examples-paper/demo_long_fiber.m

    % Setup
    [specquad1,specquad2,specquad3] = deal(zeros(size(X)));
    cancellation_errest = nan*zeros(size(X));
    corrneeded = false(numel(X),1);
    time_weights = 0;
    kerevals_near = 0;
    maintic = tic();
    time_ker_near = 0;
    time_coeffs = 0;
    time_far = 0;

    % Compute quadrature weights
    atic = tic();
    [all_w1,all_w3,all_w5,all_mu1,all_mu3,all_mu5,specquad_needed,all_roots] = ...
        closed_curve_near_weights(tj, wj, xj, yj, zj, X, Y, Z, tol);

    % Evaluation count
    time_weights = time_weights + toc(atic);
    kerevals_near = kerevals_near + nquad*sum(specquad_needed(:));

    % Evaluate each panel-to-point pair
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
            % integrands
            g1R1 = (sj.*u1R1).';
            g1R3 = (sj.*u1R3).';
            g1R5 = (sj.*u1R5).';
            g2R1 = (sj.*u2R1).';
            g2R3 = (sj.*u2R3).';
            g2R5 = (sj.*u2R5).';
            g3R1 = (sj.*u3R1).';
            g3R3 = (sj.*u3R3).';
            g3R5 = (sj.*u3R5).';
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
            if (correct_R3 || correct_R5) && use_mod
                w3 = all_w3(:,i);
                w5 = all_w5(:,i);
                normw3 = norm(w3,inf);
                normw5 = norm(w5,inf);
                toln = tol;
                % estimate for 1/R^5
                corrR35 = zeros(3,2);
                corrR35(1,2) = cond_sum(normw5,g1R5,I1R5);
                corrR35(2,2) = cond_sum(normw5,g2R5,I2R5);
                corrR35(3,2) = cond_sum(normw5,g3R5,I3R5);
                corrR5 = sum(corrR35(:,2) > toln) > 0;
                if ~corrR5
                    % 1/R^5 ok, now check 1/R^3
                    corrR35(1,1) = cond_sum(normw3,g1R3,I1R3);
                    corrR35(2,1) = cond_sum(normw3,g2R3,I2R3);
                    corrR35(3,1) = cond_sum(normw3,g3R3,I3R3);
                    corrR3 = sum(corrR35(:,1) > toln) > 0;
                else
                    % 1/R^5 bad, assume same corr needed for 1/R^3
                    corrR3 = true;
                    corrR35(:,1) = corrR35(:,2);
                end
                cancellation_errest(i) = max(corrR35,[],'all');
    
                if corrR5 || corrR3
                    corrneeded(i) = true;
                    kstd = get_k_vec(nquad,2*pi).';
                    kmod = get_k_vec(nquad-2,2*pi).';
                    t0 = all_roots(i); % with real part in [0,2*pi)
                    a = real(t0);
                    b = imag(t0);
                    r = exp(-abs(b));
                    tdist = abs(exp(1i*tj)-exp(1i*t0));
                    estd = exp(1i*kstd.'*a);
                    estdik = 1i*kstd.'.*estd;
                    emod = exp(1i*kmod*a);
                    mu1 = all_mu1(:,i);
                    mu3 = all_mu3(:,i);
                    mu5 = all_mu5(:,i);
                    if corrR5
                        p3mod = zeros(numel(kstd),1);
                        p5mod = zeros(numel(kstd),1);
                        p3modk = 0.5*(-(1-r)^2/(2*r*(1+r^2)).*mu3(2:end-1) + ((3/2+kmod-1)./(2*r).*mu1(2:end-1)-(3/2+kmod-2)./(1+r^2).*mu1(1:end-2))./(3/2-1));
                        p5modk = 0.5*(1-r)^(3-5) * (-(1-r)^2/(2*r*(1+r^2)).*mu5(2:end-1) + ((5/2+kmod-1)./(2*r).*mu3(2:end-1)-(5/2+kmod-2)./(1+r^2).*mu3(1:end-2))./(5/2-1));
                        if mod(nquad,2) == 0
                            p3mod(1) = 2*mu3(nquad/2+1)/(1-r)^(3-1);
                            p5mod(1) = 2*mu5(nquad/2+1)/(1-r)^(5-1);
                        else
                            p3mod(1) = 2*mu3((nquad-1)/2+1)/(1-r)^(3-1);
                            p5mod(1) = 2*mu5((nquad-1)/2+1)/(1-r)^(5-1);
                        end
                        p3mod(3:end) = emod.*p3modk;
                        p5mod(3:end) = emod.*p5modk;
                    else
                        p3mod = zeros(numel(kstd),1);
                        p3modk = 0.5*(-(1-r)^2/(2*r*(1+r^2)).*mu3(2:end-1) + ((3/2+kmod-1)./(2*r).*mu1(2:end-1)-(3/2+kmod-2)./(1+r^2).*mu1(1:end-2))./(3/2-1));
                        if mod(nquad,2) == 0
                            p3mod(1) = 2*mu3(nquad/2+1)/(1-r)^(3-1);
                        else
                            p3mod(1) = 2*mu3((nquad-1)/2+1)/(1-r)^(3-1);
                        end
                        p3mod(3:end) = emod.*p3modk;
                    end
    
                    % store all relevant values 
                    GR35 = zeros(3*nquad,2);
                    GR35(:,1) = [g1R3;g2R3;g3R3];
                    GR35(:,2) = [g1R5;g2R5;g3R5];
                    IR35 = zeros(3,2);
                    IR35(:,1) = [I1R3;I2R3;I3R3];
                    IR35(:,2) = [I1R5;I2R5;I3R5];
        
                    % eval layer dens at a=real(t0)
                    f1coeff = fftshift(fft(f1j)).'/nquad;
                    f2coeff = fftshift(fft(f2j)).'/nquad;
                    f3coeff = fftshift(fft(f3j)).'/nquad;
                    f1a = real(estd*f1coeff);
                    f2a = real(estd*f2coeff);
                    f3a = real(estd*f3coeff);
                    %f1a = x(a)/(2*pi); f2a = y(a)/(2*pi); f3a = z(a)/(2*pi); % analytic evaluation
                    % eval remaining integrand analytically
                    xa = x(a); ya = y(a); za = z(a);
                    r1a = xa-Xi; r2a = ya-Yi; r3a = za-Zi;
                    tdista = abs(exp(1i*a)-exp(1i*t0));
                    sa = s(a);
                    [~, u1R3a, u1R5a, ~, u2R3a, u2R5a, ~, u3R3a, u3R5a] ...
                        = slender_body_kernel_split(r1a, r2a, r3a, f1a, f2a, f3a, slender_eps);
                    uR35a = [u1R3a, u1R5a;
                             u2R3a, u2R5a;
                             u3R3a, u3R5a];
                    sin2fac = sin((tj-a)./2).^2;
       
                    tdistMat = zeros(nquad,2);
                    tmp = tdist.^(2*1+1);
                    tdistMat(:,1) = tmp;
                    tdistMat(:,2) = tmp.*tdist.*tdist;
                    for ii = 1:3
                        for jj = 1:2
                            corr = corrR35(ii,jj) > toln;
                            if ~isnan(corr) && corr
                                g = GR35(((ii-1)*nquad+1):ii*nquad,jj);
                                h = g.*tdistMat(:,jj); % to expand in modified Fourier basis
                                ccoeff = fftshift(fft(h))/nquad; % std Fourier coefficients
                                % compute modified Fourier coeffs
                                btic = tic();
                                % alt 1: transform coefficients (slow)
%                                 [a0,a1,bcoeff] = fourier2modcoeffs(ccoeff,a); % map c_k --> (a_0,a_1,b_k)
%                                 d1coeff = sa*uR35a(ii,jj)*tdista^(2*jj+1); % correction to d(1)
%                                 dcoeff = [d1coeff;a1;bcoeff];
                                % alt 2: interpolate and FFT (faster, ~2x)
                                a0 = sa*uR35a(ii,jj)*tdista^(2*jj+1);
                                a1 = real(estdik*ccoeff);
                                gj = (h-a0-a1*sin(tj-a))./sin2fac;
                                bk = fftshift(fft(gj)/nquad);
                                dcoeff = [a0;a1;bk(2:end-1)];
                                time_coeffs =  time_coeffs + toc(btic);
                                if jj == 1 && correct_R3
                                    IR35(ii,jj) = real(sum(dcoeff.*p3mod));
                                elseif jj == 2 && correct_R5
                                    IR35(ii,jj) = real(sum(dcoeff.*p5mod));
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
            time_ker_near =  time_ker_near + toc(atic);
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
    fprintf('Near field kernel evals time: %f\n', time_ker_near)
    fprintf('Total time line3_near_weights: %f\n', time_weights)
    fprintf('Far field time: %f\n', time_far)
    fprintf('Number of target points: %d\n', numel(X))
    fprintf('Number special quadrature: %d\n', sum(specquad_needed))
    fprintf('Number corrections: %d\n', sum(corrneeded))
    toc(maintic)
    stats.kerevals_near = kerevals_near;
    stats.time_weights = time_weights;
    stats.time_ker_near = time_ker_near;
    stats.time_coeffs = time_coeffs;
end