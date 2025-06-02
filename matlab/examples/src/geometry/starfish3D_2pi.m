function [x, y, z, xp, yp, zp, xpp, ypp, zpp] = starfish3D_2pi(amplitude, n_arms, radius)
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

% Second derivatives
ddcommon = @(t) -radius * amplitude * n_arms^2 * cos(n_arms * t(:)');

xpp = @(t) ddcommon(t(:)') .* cos(t(:)') - 2 * dcommon(t(:)') .* sin(t(:)') - common(t(:)') .* cos(t(:)');
ypp = @(t) ddcommon(t(:)') .* sin(t(:)') + 2 * dcommon(t(:)') .* cos(t(:)') - common(t(:)') .* sin(t(:)');
zpp = @(t) -wobble_amp * wobble_freq^2 * sin(wobble_freq * t(:)');
end