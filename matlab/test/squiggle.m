function [x,y,z,xp,yp,zp] = squiggle(L,K)
% SQUIGGLE generates a smooth closed squiggle curve in R^3 via Fourier series
%
% [x,y,z,xp,yp,zp] = squiggle(L,K) returns function handles defining the
%   closed loop in t in [0,(2*pi)/L] with max Fourier mode K
% 
% INPUTS:
%   L - parameter scale, t in [0,(2*pi)/L]
%   K - max Fourier mode in squiggle
%
% OUTPUTS:
%   x, y, z    - function handles for curve coordinates
%   xp, yp, zp - derivatives of the curve
%
% NOTE: Function modified from https://github.com/ludvigak/linequad

K0 = 5; % mode beyond which decay kicks in
k = -K:K;
nk = numel(k);
rng(0);
ampl = (K0+abs(k)).^(-1); % Sobolev-type decay
c = (randn(3,nk)+1i*randn(3,nk)).*ampl; % cmplx F coeffs in each coord
x = @(t) real(c(1,:)*exp(L*1i*k'*t(:).')); % outer prod to eval the F series
y = @(t) real(c(2,:)*exp(L*1i*k'*t(:).'));
z = @(t) real(c(3,:)*exp(L*1i*k'*t(:).'));
xp = @(t) real(c(1,:)*(L*1i*k'.*exp(L*1i*k'*t(:).')));
yp = @(t) real(c(2,:)*(L*1i*k'.*exp(L*1i*k'*t(:).')));
zp = @(t) real(c(3,:)*(L*1i*k'.*exp(L*1i*k'*t(:).')));
end
