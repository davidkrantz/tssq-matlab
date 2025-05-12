clear all;
close all;

[x, y, z, xp, yp, zp, xpp, ypp, zpp] = starfish3D(0.3,3,1);
s = @(t) sqrt(xp(t).^2 + yp(t).^2 + zp(t).^2);

tj = linspace(0,2*pi,1000).';

xj = x(tj); yj = y(tj); zj = z(tj);

tol = 1e-10;
[tj,wj] = adaptive_global_discretization(s,tol);
sj = s(tj);
numel(tj)

figure;
plot3(xj,yj,zj,'k.');
hold on;
plot3(xj(1),yj(1),zj(1),'ro');
plot3(xj(end),yj(end),zj(end),'gp');
axis equal;
grid on;

figure;
plot(tj,sj,'.');
grid on;

alignfigs;

function [tj,wj] = adaptive_global_discretization(s,tol)
n = 4;
res = tol+1;
itr = 0;
maxitr = 10;
while res > tol && itr < maxitr
    tj = linspace(0,2*pi,n+1).'; tj(end) = []; % periodic grid in [0,2*pi)
    sj = s(tj);
    c = fftshift(fft(sj))/n;
    res = max(abs(c(end-1:end))/max(abs(c)));
    n = 2*n;
    itr = itr+1;
end
wj = (tj(2)-tj(1))*ones(n,1);
end

function [x, y, z, xp, yp, zp, xpp, ypp, zpp] = starfish3D(amplitude, n_arms, radius)
if nargin < 3
    radius = 1.0;
end

% Smooth vertical wobble parameters
wobble_amp = 2; % vertical amplitude
wobble_freq = 1; % number of z-oscillations around the loop

% Position functions
x = @(t) radius * (1 + amplitude * cos(n_arms * t)) .* cos(t);
y = @(t) radius * (1 + amplitude * cos(n_arms * t)) .* sin(t);
z = @(t) wobble_amp * sin(wobble_freq * t);  % small 3D vertical variation

% First derivatives
common = @(t) radius * (1 + amplitude * cos(n_arms * t));
dcommon = @(t) -radius * amplitude * n_arms * sin(n_arms * t);

xp = @(t) dcommon(t) .* cos(t) - common(t) .* sin(t);
yp = @(t) dcommon(t) .* sin(t) + common(t) .* cos(t);
zp = @(t) wobble_amp * wobble_freq * cos(wobble_freq * t);

% Second derivatives
ddcommon = @(t) -radius * amplitude * n_arms^2 * cos(n_arms * t);

xpp = @(t) ddcommon(t) .* cos(t) - 2 * dcommon(t) .* sin(t) - common(t) .* cos(t);
ypp = @(t) ddcommon(t) .* sin(t) + 2 * dcommon(t) .* cos(t) - common(t) .* sin(t);
zpp = @(t) -wobble_amp * wobble_freq^2 * sin(wobble_freq * t);
end
