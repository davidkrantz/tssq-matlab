function [a0,a1,b_coeffs] = fourier2modcoeffs(c_coeffs,a)
% FOURIER2MODCOEFFS converts Fourier basis coefficients to a modified
%   basis. Assuming N even, the regular expansion is
%
%   \sum_{k=-N/2}^{N/2-1} c_k e^{ikx}
%
%   while the modified expansion is
%
%   a_0 + a_1\sin(x-a) + \sum_{k=-N/2+1}^{N/2-2} b_k \sin((x-a)/2)^2 e^{ikx}.
%
% [a0,a1,b_coeffs] = fourier2modcoeffs(c_coeffs,a) takes the Fourier
%   coefficents c_coeffs (length-N vector) and maps them to the
%   coefficients of the modified basis defined by the shift a.
%
% Without arguments runs test defined below
%
% INPUTS:
%   c_coeffs - length-N vector of complex Fourier coefficients, assumed to be
%              ordered with symmetric wavenumbers using get_k_vec(N, 2*pi)
%   a        - real shift parameter defining the modified basis
%
% OUTPUTS:
%   a0       - scalar for constant term in the modified basis
%   a1       - scalar for sin(x - a) term in the modified basis
%   b_coeffs - length-(N-2) vector of complex coefficients for the modified
%              exponential basis modulated by sin((x-a)/2)^2
%
% AUTHOR: David Krantz (davkra@kth.se), May 2025
%
% NOTE: Function is based on the derivations and code 
%   test_coeffs_in_new_basis_AK.m by Anna-Karin Tornberg.

if nargin == 0, test_transform_coeffs; return; end

N = length(c_coeffs); % total number of standard Fourier modes

% constants
w = exp(1i*a); % exponential factors for shift a
winv = conj(w);
alpha1 = 2*w;
alpha2 = -w*w;
alpha3 = -4*w;
alphainv1 = 2*winv;
alphainv2 = -winv*winv;
alphainv3 = -4*winv;

b_coeffs = zeros(N-2,1); % modified basis coefficients output vector

% determine index of mode k=0 in c_coeffs and b_coeffs
if mod(N,2) == 0
    ind0_b = N/2; % wavenumers = -N/2+1 : N/2-2
    ind0_c = N/2+1; % wavenumbers = -N/2 : N/2 - 1
else
    ind0_b = (N-1)/2; % wavenumbers = -(N-1)/2+1 : (N-1)/2-1
    ind0_c = (N-1)/2+1; % wavenumbers = -(N-1)/2 : (N-1)/2
end

% positive modes
b_coeffs(N-2) = alpha3*c_coeffs(N);
b_coeffs(N-3) = alpha1*b_coeffs(N-2) + alpha3*c_coeffs(N-1);
for j = (N-4):-1:(ind0_b+1) % from high to low
    b_coeffs(j) = alpha1*b_coeffs(j+1) + alpha2*b_coeffs(j+2) + alpha3*c_coeffs(j+2);
end

% negative modes
b_coeffs(1) = alphainv3*c_coeffs(1);
b_coeffs(2) = alphainv1*b_coeffs(1) + alphainv3*c_coeffs(2);
for j = 3:(ind0_b-1) % from low to high
    b_coeffs(j) = alphainv1*b_coeffs(j-1) + alphainv2*b_coeffs(j-2) + alphainv3*c_coeffs(j);
end

% helper variables for a0, a1 and b0
d1 = -0.5*b_coeffs(ind0_b+1) + 0.25*w*b_coeffs(ind0_b+2) + c_coeffs(ind0_c + 1);
d2 = 0.25*winv*b_coeffs(ind0_b-1) + 0.25*w*b_coeffs(ind0_b+1) + c_coeffs(ind0_c);
d3 = -0.5*b_coeffs(ind0_b-1) + 0.25*winv*b_coeffs(ind0_b-2) + c_coeffs(ind0_c - 1);

b_coeffs(ind0_b) = -2*(w*d1 + winv*d3); % mode k=0, b_{0}
a1 = 1i*(w*d1 - winv*d3); % coeff of sin(t-a)
a0 = d2 - 0.5*b_coeffs(ind0_b); % coeff of constant function
end

function test_transform_coeffs

N = 128; % number of regular Fourier modes
a = pi/11; % shift in modified basis

kvec_b = get_k_vec(N-2,2*pi); % wavenumbers
kvec_c = get_k_vec(N,2*pi);

% random coefficients, but make sure they decay sufficiently
c_coeffs = (rand(size(kvec_c)) - 0.5)*2 + 1i*(rand(size(kvec_c)) - 0.5)*2;
c_coeffs = c_coeffs .* exp(-0.01 * kvec_c.^2);

[a0, a1, b_coeffs] = fourier2modcoeffs(c_coeffs, a); % transform coeffs

Np = 2^nextpow2(N*2); % evaluation grid
xv = (0:Np-1)/Np*2*pi;

Qsum = zeros(size(xv)); % evaluate original sum
for k = 1:N
    Qsum = Qsum + c_coeffs(k) * exp(1i*kvec_c(k)*xv);
end

Psum = a0+a1*sin(xv-a); % evaluate modified basis sum
for k = 1:(N-2)
    Psum = Psum + b_coeffs(k)*exp(1i*kvec_b(k)*xv).*sin((xv-a)/2).^2;
end

maxabserr = norm(Qsum-Psum,inf);
maxrelerr = norm((Qsum-Psum)./Qsum,inf);

fprintf('max absolute recon error: %.14e\n', maxabserr);
fprintf('max relative recon error: %.14e\n', maxrelerr);

end
