function curve = starfish3D_2pi_new(amplitude, n_arms, radius)
% t in [0,2pi)
if nargin < 3
    radius = 1.0;
end

% Smooth vertical wobble parameters
wobble_amp = 2; % vertical amplitude
wobble_freq = 1; % number of z-oscillations around the loop

% Position functions
x = @(t) radius * (1 + amplitude * cos(n_arms * t(:)')) .* cos(t(:)');
y = @(t) radius * (1 + amplitude * cos(n_arms * t(:)')) .* sin(t(:)');
z = @(t) wobble_amp * sin(wobble_freq * t(:)');  % small 3D vertical variation

% First derivatives
common = @(t) radius * (1 + amplitude * cos(n_arms * t(:)'));
dcommon = @(t) -radius * amplitude * n_arms * sin(n_arms * t(:)');

xp = @(t) dcommon(t(:)') .* cos(t(:)') - common(t(:)') .* sin(t(:)');
yp = @(t) dcommon(t(:)') .* sin(t(:)') + common(t(:)') .* cos(t(:)');
zp = @(t) wobble_amp * wobble_freq * cos(wobble_freq * t(:)');

s = @(t) sqrt(xp(t).^2 + yp(t).^2 + zp(t).^2);

curve.x = x;
curve.y = y;
curve.z = z;
curve.xp = xp;
curve.yp = yp;
curve.zp = zp;
curve.s = s;
end