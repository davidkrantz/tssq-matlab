function [W3,W5] = adjoint_modified_fourier(tj,a,Stilde3,Stilde5)
% ADJOINT_MODIFIED_FOURIER returns the oscillatory part of the adjoint
%   quadrature weights

if nargin == 0, test_adjoint_weights; return; end

n = numel(tj);

% precompute exponentials
ea = exp(1i*a);
eam = conj(ea);
alpha1 = 2*eam;
alpha2 = -eam*eam;
alpha3 = -4*eam;

% determine index of mode k=0
if mod(n,2) == 0
    ind0 = n/2+1; % wavenumbers = -N/2 : N/2 - 1
    ind_flip_start = 2; % for symmetry fill later
    else
    ind0 = (n-1)/2+1; % wavenumbers = -(N-1)/2 : (N-1)/2
    ind_flip_start = 1; % for symmetry fill later
end

if nargout == 1
    % allocate
    What3 = zeros(n,1);

    % negative modes
    What3(ind0-1) = -2*eam*Stilde3(ind0-1);
    for j = (ind0-2):-1:1 % from low to high
        What3(j) = alpha1*What3(j+1)+alpha2*What3(j+2)+alpha3*Stilde3(j);
    end

    % positive modes through conjugate symmetry
    What3(ind0+1:end) = flipud(conj(What3(ind_flip_start:ind0-1)));

    % recover node weights (inverse of What_k = \sum_j W_j e^{+ik t_j})
    Whatshift3 = ifftshift(What3);
    W3 = ifft(conj(Whatshift3),'symmetric'); % slightly faster than below
    % W3 = real(fft(Whatshift3))/n;

    tdist3 = tdist.^3;
    W3 = W3*tdist3;
else
    % allocate
    What3 = zeros(n,1);
    What5 = zeros(n,1);

    % negative modes
    What3(ind0-1) = -2*eam*Stilde3(ind0-1);
    What5(ind0-1) = -2*eam*Stilde5(ind0-1);
    for j = (ind0-2):-1:1 % from low to high
        What3(j) = alpha1*What3(j+1)+alpha2*What3(j+2)+alpha3*Stilde3(j);
        What5(j) = alpha1*What5(j+1)+alpha2*What5(j+2)+alpha3*Stilde5(j);
    end

    % positive modes through conjugate symmetry
    What3(ind0+1:end) = flipud(conj(What3(ind_flip_start:ind0-1)));
    What5(ind0+1:end) = flipud(conj(What5(ind_flip_start:ind0-1)));
    
    % recover node weights (inverse of What_k = \sum_j W_j e^{+ik t_j})
    Whatshift = ifftshift([What3,What5],1);
    W = ifft(conj(Whatshift),[],1,'symmetric');
    W3 = W(:,1);
    W5 = W(:,2);
end

end

function test_adjoint_weights

rng(123);

nvec = [510, 511]; % try even and odd n

for i = 1:numel(nvec)
    n = nvec(i);
    % determine index of mode k=0
    if mod(n,2) == 0
        ind0 = n/2; 
        % wavenumbers = -N/2 : N/2 - 1
        ind_flip_start = 2;
        % for symmetry fill later
    else
        ind0 = (n-1)/2; 
        % wavenumbers = -(N-1)/2 : (N-1)/2
        ind_flip_start = 1;
        % for symmetry fill later
    end
    
    tj = linspace(0,2*pi,n+1).'; tj(end) = []; % nodes in [0,2*pi)
    
    % example root, try both off and on grid point
    if i == 1
        a = 2*pi*rand;
    else
        a = tj(11); % grid point
    end
    
    % example functions
    h = @(t) sin(t-a); % function small at t=a
    g = @(t) h(t).^2+1e-8;
    
    % setup variables
    ga = g(a);
    gj = g(tj);
    w0 = trig_barycentric_weights(tj,a);
    kmod = get_k_vec(n-2,2*pi).';
    Stilde3 = zeros(n-2,1);
    Stilde3(1:ind0-1) = rand(numel(1:ind0-1),1) + 1i*(rand(numel(1:ind0-1),1)-0.5);
    Stilde3(ind0) = rand;
    Stilde3(ind0+1:end) = flipud(conj(Stilde3(ind_flip_start:ind0-1)));
    Stilde5 = 2*Stilde3;
    pmod3 = [1e5*rand;0;Stilde3];
    pmod5 = [1e10*rand;0;Stilde5];
    
    % Reference weights (Vandermonde matrix approach)
    Vosc = sin((tj-a)/2).^2.*exp(1i*kmod.'.*tj); % oscilattory basis functions
    V = [ones(n,1) sin(tj-a) Vosc]; % add constant function 1 and shifted sine
    W3osc = V.'\[0; 0; pmod3(3:end)];
    W5osc = V.'\[0; 0; pmod5(3:end)];
    L3ref = real(W3osc.*gj+ga*pmod3(1)*w0);
    L5ref = real(W5osc.*gj+ga*pmod5(1)*w0);
    
    % FFT-accelerated weights
    [W3osc_fft,W5osc_fft] = adjoint_modified_fourier(tj,a,Stilde3,Stilde5);
    L3fft = real(W3osc_fft.*gj+ga*pmod3(1)*w0);
    L5fft = real(W5osc_fft.*gj+ga*pmod5(1)*w0);
    
    % errors
    maxabserr3 = norm(L3ref-L3fft,inf);
    maxrelerr3 = norm((L3ref-L3fft)./L3ref,inf);
    maxabserr5 = norm(L5ref-L5fft,inf);
    maxrelerr5 = norm((L5ref-L5fft)./L5ref,inf);
    
    fprintf('random values test with n=%.1i:\n', n);
    fprintf('max absolute error m=3: %.14e\n', maxabserr3);
    fprintf('max relative error m=3: %.14e\n', maxrelerr3);
    fprintf('max absolute error m=5: %.14e\n', maxabserr5);
    fprintf('max relative error m=5: %.14e\n', maxrelerr5);
end
end
