function curve = starfish3D(amplitude, n_arms, radius)
% STARFISH3D  Generate a smooth 3D starfish-shaped closed curve.
%
%   curve = starfish3D(amplitude, n_arms, radius)
%
%   INPUT:
%     amplitude : radial modulation amplitude
%     n_arms    : number of starfish arms
%     radius    : base radius (default 1)
%
%   OUTPUT:
%     curve struct with fields:
%       x(t), y(t), z(t)   : position
%       xp(t), yp(t), zp(t): derivatives wrt t
%       s(t)               : speed
%
%   NOTE:
%     t in [0, 1)

if nargin < 3
    radius = 1.0;
end

% Smooth vertical wobble parameters
wobble_amp = 2; % vertical amplitude
wobble_freq = 1; % number of z-oscillations around the loop

% Position functions
x = @(t) radius * (1 + amplitude * cos(n_arms * 2*pi*t(:)')) .* cos(2*pi*t(:)');
y = @(t) radius * (1 + amplitude * cos(n_arms * 2*pi*t(:)')) .* sin(2*pi*t(:)');
z = @(t) wobble_amp * sin(wobble_freq * 2*pi*t(:)');  % small 3D vertical variation

% First derivatives
common = @(t) radius * (1 + amplitude * cos(n_arms * 2*pi*t(:)'));
dcommon = @(t) -radius * amplitude * n_arms * sin(n_arms * 2*pi*t(:)');

xp = @(t) dcommon(t(:)') .* cos(2*pi*t(:)') - common(t(:)') .* sin(2*pi*t(:)');
yp = @(t) dcommon(t(:)') .* sin(2*pi*t(:)') + common(t(:)') .* cos(2*pi*t(:)');
zp = @(t) wobble_amp * wobble_freq * cos(wobble_freq * 2*pi*t(:)');

s = @(t) sqrt(xp(t).^2 + yp(t).^2 + zp(t).^2);

curve.x = x;
curve.y = y;
curve.z = z;
curve.xp = xp;
curve.yp = yp;
curve.zp = zp;
curve.s = s;
end
