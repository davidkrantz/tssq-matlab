% plots f in basis integrand having numerator f*cos
%
% CONCLUSION: f is non-smooth, requires many points to compute integral in
% [0,2*pi] using FFT/trapz. Not feasible!
%
% ABORTED
clear all;
close all;

m = 3;
l = 2;
rv = linspace(0.9,1,500).';
tv = linspace(0,2*pi,500);
M = numel(rv);

ell = @(r) -(4*r)./(1-r).^2;
integrand = @(r,t) (sin(t./2).^l)./(1-ell(r).*sin(t./2).^2).^(m/2);

f = integrand(rv,tv);

%%
close all;

figure;
surf(tv,rv,f);
shading interp;
xlabel('t');
ylabel('r');
zlabel('integrand');
colorbar;