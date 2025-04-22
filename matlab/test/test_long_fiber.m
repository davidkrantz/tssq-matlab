% Evaluate Stokes slender body potential from closed-loop squiggle curve.
% Computes using interpolatory quadrature, and compares to semi-smart adaptive
% quadrature, which is reasonably fast and seems to compute the right thing.
% Slow computations of reference using integral() is also available.
%
% Matlab code attempts to be fast by avoiding inner loop vector ops, so that
% we can compute meaningful running times.

% Alex variant also includes all targs close to curve test - see targtype.

% Alex added (6) case for 3-digit slice test at nquad=4

function test_long_fiber(varargin) % Function makes Matlab JIT compile better

clear all;
format long;

SAVEPLOTS = false;

if nargin==0
    test_no = 7;
else
    test_no = varargin{1};
end

%% Default setup
slender_eps = 1e-3;
%WHICH_REFERENCE = 'integral';
WHICH_REFERENCE = 'adaptive';

UPSAMPLE = true; % Upsample before applying interpolatory quadrature
nquad = 16;
rho = 3; % Interpolatory quadrature limit rule
%Hlim = 1; % Adaptive quadrature distance rule
Hlim = 2;
%rho = 1 % Disable specquad

% When running near points
dist = 1e-3;            % the dist of all pts from the curve

%% Case setup
switch test_no
  case 0
    % Play
    targtype = 'n';   % 's' for original slice, or 'n' for all targs near curve at random locs.
    tol = 1e-6;
  case 1
    % Paper plot for field
    targtype = 's';
    tol = 1e-10;      
  case 2
    % Comparison for nearby points
    targtype = 'n';    
    tol = 1e-6;
    dist = 1e-2;
  case 3
    targtype = 'n';    
    tol = 1e-10;
    dist = 1e-2;
  case 4
    targtype = 'n';    
    tol = 1e-6;
    dist = 1e-4;
  case 5
    targtype = 'n';    
    tol = 1e-10;
    dist = 1e-4;
  case 6                    % try lower nquad for lower acc.
   nquad = 4;  % #digits + 1
   targtype = 's';    
   tol = 1e-3;
  case 7    % David
   targtype = 'n';
   tol = 1e-6;
   dist = 1e-8;
  otherwise
    error('Unknown test no')
end

%% Run

% Setup fiber
[x, y, z, xp, yp, zp] = squiggle();
s = @(t) sqrt(xp(t).^2 + yp(t).^2 + zp(t).^2);

% Artificial density (trivial)
f1 = @(t) x(t);
f2 = @(t) y(t);
f3 = @(t) z(t);

% Discretize
fprintf('* Discretizing, tol=%.1e\n', tol)
[tj, wj, npan, edges] = adaptive_panelization(s, nquad, tol);
fprintf('nquad=%d, npan=%d\n', nquad, npan);
xj = x(tj); yj = y(tj); zj = z(tj);
xpj = xp(tj); ypj = yp(tj); zpj = zp(tj);
sj = s(tj);
f1j = f1(tj); f2j = f2(tj); f3j = f3(tj);


if targtype=='s'     % Evaluation surface - slice
  Neval = 100;    % in each direction
  if SAVEPLOTS
      Neval = 200;
  end
  [X, Z, Y] = meshgrid( linspace(-1.4,1.4, Neval), ...                      
                        linspace(-1.4,1.4, Neval), ...
                        0.25);
elseif targtype=='n'  % all targs near curve, in random normal directions a fixed dist away
  Ne = 5e3;  % # targs
  t = rand(1,Ne); 
  v = randn(3,Ne); utang = [xp(t);yp(t);zp(t)]./s(t); % sloppy unit tangents
  vdotutang = sum(v.*utang,1); v = v - utang.*vdotutang;  % orthog v against the tangent
  v = v./sqrt(sum(v.*v,1));    % normalize all the v vecs
  X = x(t) + dist*v(1,:);   % displace by v vecs from pts on curve
  Y = y(t) + dist*v(2,:);
  Z = z(t) + dist*v(3,:);
end

% Compute ref solution
[uref1, uref2, uref3] = deal(zeros(size(X)));
switch WHICH_REFERENCE
  case 'integral'
    % Use integral() reference (slooooow, so only do one component)
    uref1 = zeros(size(X));
    disp('* Reference: integral()')
    tic
    parfor i=1:numel(X)
        uref1(i) = integral(@(t) s(t) .* ...
                            slender_body_kernel(x(t)-X(i), y(t)-Y(i), z(t)-Z(i), f1(t), f2(t), f3(t), slender_eps), ...
                            0, 1, 'abstol', 1e-15, 'reltol', 1e-15);
    end
    toc
  case 'adaptive'
    % Adaptive quadrature
    disp('* Reference: Adaptive')
    % Make sure that we get a different discretization for reference computations,
    % otherwise errors get artifically small.
    %nquad_ref = nquad+2;
    nquad_ref = 2*nquad;
    tol_ref = 5e-14;
    if nquad<16, tol_ref = 1e-2*tol;  end    % low-acc test cases
    %Hlim_ref = 1.0;
    Hlim_ref = 2.0;
    [tj_ref, wj_ref, npan_ref] = adaptive_panelization(s, nquad_ref, tol_ref);
    fprintf('Discretization: nquad=%d, npan=%d\n', nquad_ref, npan_ref);
    [uref1, uref2, uref3] = adaptive_quadrature(x(tj_ref), y(tj_ref), z(tj_ref), s(tj_ref), ...
                                                wj_ref, ...
                                                f1(tj_ref), f2(tj_ref), f3(tj_ref), ...
                                                X, Y, Z, nquad_ref, slender_eps, Hlim_ref);
end

% Compute adaptive quadrature
disp(' ')
disp('* Adaptive quadrature')
[adquad1, adquad2, adquad3, adstats] = adaptive_quadrature(xj, yj, zj, sj, wj, f1j, f2j, f3j, ...
                                             X, Y, Z, nquad, slender_eps, Hlim);

% Compute special quadrature
disp(' ')
disp('* Interpolatory quadrature')
[specquad1,specquad2,specquad3,specquadsh1,specquadsh2,specquadsh3,corr_needed,specstats] = interpolatory_quadrature(...
    xj, yj, zj, sj, tj, wj, f1j, f2j, f3j, X, Y, Z, ...
    x, y, z, xp, yp, zp, s, nquad, edges, rho, UPSAMPLE, slender_eps, tol);

% Compute errors
if strcmp(WHICH_REFERENCE, 'integral')
    % We only have reference for 1st component
    adquad2 = uref2;
    adquad3 = uref3;
    specquspec2 = uref2;
    specquspec3 = uref3;
end

adquad_errmax = compute_error(uref1, uref2, uref3, adquad1, adquad2, adquad3);
adquad_errmaxmax = max(adquad_errmax(:))

specquad_errmax = compute_error(uref1, uref2, uref3, specquad1, specquad2, specquad3);
specquad_errmaxmax = max(specquad_errmax(:))

specquadsh_errmax = compute_error(uref1, uref2, uref3, specquadsh1, specquadsh2, specquadsh3);
specquadsh_errmaxmax = max(specquadsh_errmax(:))

%%
close all;

% Plot
if SAVEPLOTS
    pubfig = @publication_fig;
else
    pubfig = @() false;
end

sfigure(1);
clf; pubfig();
plot3(xj, yj, zj, '.-k')
if targtype=='n', hold on; plot3(X,Y,Z,'r.'); end
axis equal tight vis3d
grid on
box on

sfigure(2);
clf; pubfig();
if targtype=='s'
  surf(X, Y, Z, log10(specquad_errmax))                  
  shading flat
  xlabel(colorbar(), 'log_{10} E_{rel}')  
  caxis auto
  hold on
  plot3(xj, yj, zj, '.-k')
  axis equal tight vis3d
  grid off
  box on
else
  semilogy(t(corr_needed),specquad_errmax(corr_needed),'o','markersize',3)
  hold on;
  semilogy(t,specquad_errmax,'.');
  xlabel('t'); ylabel('err')
  hold on
  semilogy(t,specquadsh_errmax,'.')
  semilogy(t,adquad_errmax,'.')
  yline(tol,'k',['tol=' num2str(tol)])
  legend('unshifted bad','unshifted all','shifted','adap')
  grid on
end
if ~SAVEPLOTS && targtype=='s'
    title('Interpolatory quadrature')
end

  
sfigure(3);
clf; pubfig();
if targtype=='s'
  surf(X, Y, Z, log10(adquad_errmax))                  
  shading flat
  colorbar
  caxis auto
  hold on
  plot3(xj, yj, zj, '.-k')
  axis equal tight vis3d
  grid off
  box on
else
  semilogy(t,adquad_errmax,'.'); xlabel('t'); ylabel('err')
  grid on
end
title('Adaptive quadrature')

% Test-specific plot saving
if SAVEPLOTS && test_no==1
    sfigure(2);
    axis off
    caxis([-16,-13])
    text(-1,0,-2, ['$\max E_{rel}=$', ...
                        sprintf('%.1e',specquad_errmaxmax)],'interpreter','latex')    
    saveplot(1, 'long_fiber_geo.png')
    saveplot(2, 'long_fiber_field.png')
end

% table line
adstats
fprintf(' %.1e &  %.1e & %.1e & %.2f & %.2f & %.1e & %.1e & %.2f & %.2f & %.1e \\\\ ', dist, tol, ...
        adstats.kerevals_near, adstats.time_ker_near, adstats.time_interp, adquad_errmaxmax,...
        specstats.kerevals_near, specstats.time_ker_near, specstats.time_weights, specquad_errmaxmax)
disp(['% case ' num2str(test_no)])

alignfigs;

keyboard;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveplot(num, filename)
    sfigure(num);
    publication_fig();
    path = ['../../linequadpaper/fig/', filename];
    print('-dpng', '-r200', path)
    system(['mogrify -trim ' path])    
    disp(['Wrote ' path])
end

function errmax = compute_error(uref1, uref2, uref3, q1, q2, q3)
    unorm = norm([uref1(:);uref2(:);uref3(:)], inf);
    err1 = abs(uref1-q1) ./ unorm;
    err2 = abs(uref2-q2) ./ unorm;
    err3 = abs(uref3-q3) ./ unorm;
    errmax = max(max(err1, err2), err3);
end

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

%%%%% INTERPOLATORY QUADRATURE

function [specquad1,specquad2,specquad3,specquadsh1,specquadsh2,specquadsh3,corrneeded,stats] = interpolatory_quadrature(...
    xj, yj, zj, sj, tj, wj, f1j, f2j, f3j, X, Y, Z, x, y, z, xp, yp, zp, s, nquad, edges, rho, UPSAMPLE, slender_eps, tol)
    [specquad1,specquad2,specquad3,specquadsh1,specquadsh2,specquadsh3] = deal(zeros(size(X)));
    corrneeded = false(numel(size(X)),1);
    npan = numel(xj)/nquad;    
    [tgl, wgl] = legendre.gauss(nquad);
    if UPSAMPLE
        nquad2 = 2*nquad;
        [tgl2, wgl2] = legendre.gauss(nquad2);
        B = bclag_interp_matrix(tgl, tgl2);
        %rho = sqrt(rho);
    else
        nquad2 = nquad;
        B = 1;
        tgl2 = tgl;
        wgl2 = wgl;
    end
    wbary = bclag_interp_weights(tgl2);
    time_weights = 0;
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
        % Compute quadrature weights
        atic = tic();
        [all_w1, all_w3, all_w5, specquad_needed, all_roots] = line3_near_weights(tgl2, wgl2, xjpan_up, yjpan_up, zjpan_up, ...
                                                          X, Y, Z, rho);
        % Evaluation count
        time_weights = time_weights + toc(atic);
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
                % estimate for 1/R^5
                corrR35 = false(3,2);
                corrR35(1,2) = cancellation_error_estimate(I1R5,tmp1R5) > tol;
                corrR35(2,2) = cancellation_error_estimate(I2R5,tmp2R5) > tol;
                corrR35(3,2) = cancellation_error_estimate(I3R5,tmp3R5) > tol;
                corrR5 = sum(corrR35(:,2)) > 0;
                if ~corrR5
                    % 1/R^5 ok, now check 1/R^3
                    corrR35(1,1) = cancellation_error_estimate(I1R3,tmp1R3) > tol;
                    corrR35(2,1) = cancellation_error_estimate(I2R3,tmp2R3) > tol;
                    corrR35(3,1) = cancellation_error_estimate(I3R3,tmp3R3) > tol;
                    corrR3 = logical(sum(corrR35(:,1)));
                else
                    % 1/R^5 bad, assume same corr needed for 1/R^3
                    corrR3 = true;
                    corrR35(:,1) = corrR35(:,2);
                end

                v0 = all_roots(i); % root having real part in [-1,1]
                %if false
                %if abs(imag(v0)) < 1e-3 && abs(real(v0)) < 1
                %    corrR35 = true(3,2);
                %if abs(real(v0)) < 1
                if (corrR5 || corrR3) && abs(real(v0)) < 1 %&& abs(imag(v0)) < 1e-3
                    corrneeded(i) = true;
                    sh = real(v0);
                    if corrR5
                        [~, p3sh, p5sh] = rsqrt_pow_integrals_shift(v0,nquad2);
                    else
                        [~, p3sh] = rsqrt_pow_integrals_shift(v0,nquad2);
                    end

                    % current panel endpoints
                    ta = edges(j);
                    tb = edges(j+1);
                    tsc = (tb-ta)/2;
                    %tsc = 1/2*sum(wjpan);

                    % root w real part [0,1] (actual root of gamma(t), t in [0,1])
                    t0 = (ta+tb)/2+(tb-ta)/2*v0;
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
                    f1pana = bclag_interp(f1jpan_up.',tgl2,wbary,sh);
                    f2pana = bclag_interp(f2jpan_up.',tgl2,wbary,sh);
                    f3pana = bclag_interp(f3jpan_up.',tgl2,wbary,sh);
                    %f1pana = x(a); f2pana = y(a); f3pana = z(a);
                    % eval remaining integrand analytically
                    xa = x(a); ya = y(a); za = z(a);
                    r1a = xa-Xi; r2a = ya-Yi; r3a = za-Zi;
                    tdista = abs(a-t0);
                    sa = s(a);
                    [~, u1R3a, u1R5a, ~, u2R3a, u2R5a, ~, u3R3a, u3R5a] ...
                        = slender_body_kernel_split(r1a, r2a, r3a, f1pana, f2pana, f3pana, slender_eps);
                    uR35a = [u1R3a, u1R5a;
                             u2R3a, u2R5a;
                             u3R3a, u3R5a];
    
                    for ii = 1:3
                        for jj = 1:2
                            corr = corrR35(ii,jj);
                            if corr
                                g = GR35(((ii-1)*nquad2+1):ii*nquad2,jj);
                                h = g.*tdist.^(2*jj+1); % to expand in shifted monomial basis
                                dcoeff = dvand(tgl2-sh,h);
                                d1coeff = sa*uR35a(ii,jj)*tdista^(2*jj+1); % correction to d(1)
                                %d1coeff = bclag_interp(h,tgl2,wbary,sh);
                                dcoeff(1) = d1coeff;
                                if jj == 1
                                    IR35(ii,jj) = sum(dcoeff.*p3sh)/tsc^3;
                                else
                                    IR35(ii,jj) = sum(dcoeff.*p5sh)/tsc^5;
                                end
                            end
                        end
                    end

%                     err1R3 = abs(I1R3-IR35(1,1))/abs(I1R3);
%                     err2R3 = abs(I2R3-IR35(2,1))/abs(I2R3);
%                     err3R3 = abs(I3R3-IR35(3,1))/abs(I3R3);
%                     err1R5 = abs(I1R5-IR35(1,2))/abs(I1R5);
%                     err2R5 = abs(I2R5-IR35(2,2))/abs(I2R5);
%                     err3R5 = abs(I3R5-IR35(3,2))/abs(I3R5);
%                     if err1R3 > 1e-5 || err2R3 > 1e-5 || err3R3 > 1e-5 || err1R5 > 1e-5 || err2R5 > 1e-5 || err3R5 > 1e-5
%                         keyboard;
%                     end

                    I1R3sh = IR35(1,1);
                    I2R3sh = IR35(2,1);
                    I3R3sh = IR35(3,1);
                    I1R5sh = IR35(1,2);
                    I2R5sh = IR35(2,2);
                    I3R5sh = IR35(3,2);
                    
                end

                q1 = I1R1 + I1R3 + I1R5;
                q2 = I2R1 + I2R3 + I2R5;
                q3 = I3R1 + I3R3 + I3R5;
                q1sh = I1R1 + I1R3sh + I1R5sh;
                q2sh = I2R1 + I2R3sh + I2R5sh;
                q3sh = I3R1 + I3R3sh + I3R5sh;
                % Rescale (weights are for [-1,1])
                q1 = q1/2*sum(wjpan);
                q2 = q2/2*sum(wjpan);
                q3 = q3/2*sum(wjpan);
                q1sh = q1sh/2*sum(wjpan);
                q2sh = q2sh/2*sum(wjpan);
                q3sh = q3sh/2*sum(wjpan);
                time_ker_near =  time_ker_near + toc(atic);
            else            
                atic = tic();
                [q1, q2, q3] = quadsum(xjpan, yjpan, zjpan, sjpan, wjpan, f1jpan, f2jpan, f3jpan, ...
                                       Xi, Yi, Zi, nquad, slender_eps);
                q1sh = q1; q2sh = q2; q3sh = q3;
                time_far = time_far + toc(atic);
            end
            specquad1(i) = specquad1(i) + q1;
            specquad2(i) = specquad2(i) + q2;
            specquad3(i) = specquad3(i) + q3;
            specquadsh1(i) = specquadsh1(i) + q1sh;
            specquadsh2(i) = specquadsh2(i) + q2sh;        
            specquadsh3(i) = specquadsh3(i) + q3sh;
        end
    end
    toc(maintic)
    fprintf('Near field kernel evals: %e\n', kerevals_near);
    fprintf('Near field kernel evals time: %f\n', time_ker_near)
    fprintf('Total time line3_near_weights: %f\n', time_weights)
    fprintf('Far field time: %f\n', time_far)
    stats.kerevals_near = kerevals_near;
    stats.time_weights = time_weights;
    stats.time_ker_near = time_ker_near;
end    

function est = cancellation_error_estimate(Stilde,x)
n = numel(x);
kappasum = sum(abs(x))/abs(Stilde);
est = kappasum*eps*n;
end

%%%%% ADAPTIVE QUADRATURE 2.5
    
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


%%%%% GENERAL FUNCTIONS

function [x, y, z, xp, yp, zp] = squiggle()
% random squiggle loop in t in [0,1]
    K = 20;   % max Fourier mode in squiggle
    K0 = 5;   % mode beyond which decay kicks in
    k = -K:K;
    nk = numel(k);
    rng(0);
    ampl = (K0+abs(k)).^(-1);   % Sobolev-type decay
    c = (randn(3,nk)+1i*randn(3,nk)) .* ampl;    % cmplx F coeffs in each coord
    x = @(t) real(c(1,:)*exp(2i*pi*k'*t(:)'));   % outer prod to eval the F series
    y = @(t) real(c(2,:)*exp(2i*pi*k'*t(:)'));
    z = @(t) real(c(3,:)*exp(2i*pi*k'*t(:)'));
    xp = @(t) real(c(1,:)*(2i*pi*k'.*exp(2i*pi*k'*t(:)')));
    yp = @(t) real(c(2,:)*(2i*pi*k'.*exp(2i*pi*k'*t(:)')));
    zp = @(t) real(c(3,:)*(2i*pi*k'.*exp(2i*pi*k'*t(:)')));
end

function [t, w] = panelization(npan, nquad)
    [X, W] = legendre.gauss(nquad); 
    X = (X+1)/2; W = W/2;   % on [0,1]
    t = bsxfun(@plus, X/npan, (0:npan-1)/npan);
    t = t(:);             % G-L panelization of param    
    w = repmat(W/npan, 1, npan);
    w = w(:);
end

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

%% SLENDER BODY KERNEL
% Stokeslet + e^2/2*doublet
function [u1, u2, u3] = slender_body_kernel(r1, r2, r3, f1, f2, f3, e)
    R = sqrt(r1.^2 + r2.^2 + r3.^2);
    R3 = R.*R.*R;
    R5 = R3.*R.*R;
    rdotf = r1.*f1 + r2.*f2 + r3.*f3;
    u1 = f1./R + r1.*rdotf./R3 + e^2/2*(f1./R3 - 3*r1.*rdotf./R5);
    u2 = f2./R + r2.*rdotf./R3 + e^2/2*(f2./R3 - 3*r2.*rdotf./R5);
    u3 = f3./R + r3.*rdotf./R3 + e^2/2*(f3./R3 - 3*r3.*rdotf./R5);    
end

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
