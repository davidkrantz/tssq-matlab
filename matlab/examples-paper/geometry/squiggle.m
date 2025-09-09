function curve = squiggle_new()
%SQUIGGLE  Defines a test curve (closed loop squiggle)
%   Returns struct with x(t), y(t), z(t), xp(t), yp(t), zp(t), s(t)

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
    s = @(t) sqrt(xp(t).^2 + yp(t).^2 + zp(t).^2);

    curve.x = x;
    curve.y = y;
    curve.z = z;
    curve.xp = xp;
    curve.yp = yp;
    curve.zp = zp;
    curve.s = s;
end
