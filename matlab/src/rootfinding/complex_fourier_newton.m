function [troot,converged] = complex_fourier_newton(coeffs,k,X,tinit)
% COMPLEX_FOURIER_NEWTON solves for a complex troot such that |gamma(troot)-X|^2=0
%
% [troot,converged] = complex_fourier_newton(coeffs,k,X,tinit) uses the
% Fourier coefficients coeffs of the curve gamma with wavenumbers k and
% returns the complex-valued root of squared-distance function.
%
% Without arguments runs test defined below
%
% INPUTS:
%   coeffs - Nx3 matrix (x,y,z) of Fourier coefficients of gamma
%   k      - column vector of wavenumbers (Nx1)
%   X      - 1x3 vector, target point in R^3
%   tinit  - initial guess (complex scalar)
%
% OUTPUTS:
%   troot     - complex-valued root
%   converged - boolean flag for convergence
%
% AUTHOR: David Krantz (davkra@kth.se), May 2025

if nargin == 0, test_fourier_newton; return; end

% parameters
tol = 1e-14; % convergence tolerance
maxiter_newton = 40; % maximum iterations for Newton's method
VERBOSE = false;

t = tinit;
converged = false;
ik = 1i*k;

% Newton's method
for iter = 1:maxiter_newton
    coeffse = coeffs.*exp(ik*t);
    gamma = sum(coeffse);
    gamma_p = sum(ik.*coeffse);
    
    F = sum((gamma-X).^2);
    Fprime = 2*sum((gamma-X).*gamma_p);

    dt = -F/Fprime;
    t = t+dt;

     if abs(dt) < tol
        converged = true;
        if VERBOSE
            fprintf('Newton converged in %d iterations.\n', iter);
        end
        % ensure real part of root is in [0, 2pi)
        troot = mod(real(t),2*pi) + 1i*imag(t);
        return;
    end
end

% did not converge
if VERBOSE
    fprintf('Newton did not converge after %d iterations (abs(dt)=%g)\n', maxiter_newton, abs(dt));
end
troot = mod(real(t),2*pi) + 1i*imag(t);
end

function test_fourier_newton

clear; close all; clc;

% parameters for the squiggle curve
N = 1024; % number of equidistant nodes
t_nodes = linspace(0,2*pi,N+1); % parametric points in [0, 2pi)
t_nodes(end) = []; % remove duplicate endpoint

% construct test curve
[x,y,z,~,~,~] = squiggle(1,5); % parametric curve and its derivatives

% evaluate curve at nodes
x_vals = x(t_nodes);
y_vals = y(t_nodes);
z_vals = z(t_nodes);

% compute Fourier coefficients (column vectors)
x_coeffs = fftshift(fft(x_vals)).'/N;
y_coeffs = fftshift(fft(y_vals)).'/N;
z_coeffs = fftshift(fft(z_vals)).'/N;
coeffs = [x_coeffs, y_coeffs, z_coeffs]; % (num_nodes × 3)
k = get_k_vec(N,2*pi).'; % wavenumber vector for FFT

% remove negligible Fourier modes, makes it more stable if N is large
max_vals = max(abs(coeffs), [], 1); % [max_x, max_y, max_z]
threshold = 1e-14;

% logical index for modes where all components are relatively small
remove_mask = ...
    abs(coeffs(:,1)) < threshold * max_vals(1) & ...
    abs(coeffs(:,2)) < threshold * max_vals(2) & ...
    abs(coeffs(:,3)) < threshold * max_vals(3);
coeffs(remove_mask, :) = [];
k(remove_mask) = [];

X = [0.8,0.2,0.9]; % example 3D point

% use closest point on curve as real part of initial guess
[~, minidx] = min(sum(([x_vals.' y_vals.' z_vals.']-X).^2,2));
tinit = t_nodes(minidx) + 0.1i; % add small imaginary part

% run Newton solver
[troot, converged] = complex_fourier_newton(coeffs,k,X,tinit);

% evaluate squared-distance function
gamma_root = sum(coeffs.*exp(1i*k*troot),1); % 1×3 vector
residual = abs(sum((gamma_root-X).^2));

fprintf('converged: %d, residual: %.14e\n',converged,residual);

% visualization
figure;
plot3(x_vals, y_vals, z_vals, 'k.-'); hold on;
plot3(X(1), X(2), X(3), 'r.', 'MarkerSize', 20);
plot3(x(real(troot)), y(real(troot)), z(real(troot)), 'g.', 'MarkerSize', 20);
plot3(x(real(tinit)), y(real(tinit)), z(real(tinit)), 'mo', 'MarkerSize', 8);
grid on; axis image;
xlabel('x'); ylabel('y'); zlabel('z');
legend('Curve', 'Target', 'Real part of Newton Root', 'Initial Guess');

end
