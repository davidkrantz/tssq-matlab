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
