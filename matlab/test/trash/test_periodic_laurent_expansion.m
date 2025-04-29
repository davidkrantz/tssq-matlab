% THIS BASIS WAS BAD AT RESOLVING THE UNDERLYING FUNCTIONS
% 
% ABORTED
function test_periodic_laurent_expansion(f,n,t0)

if nargin == 0, run_test; return; end

% Nodes
tj = linspace(0, 2*pi, n+1);
tj = tj(1:end-1).'; % remove duplicate 2pi
zj = exp(1i*tj); % map to unit circle
fj = f(tj); % sampled values

ttj = 2*pi*rand(100,1); % new sample points
ttj = sort(ttj);
ffj = f(ttj);

N = numel(ttj);

z0 = exp(1i*real(t0)); % shift point on unit circle

% Fourier coefficients
kmax = floor((n-1)/2);
klist_f = (-kmax:kmax).';
ck = fftshift(fft(fj))/n

% Fourier reconstruction
f_fourier = zeros(N,1);
for j = 1:N
    f_fourier(j) = sum(ck.*exp(1i*klist_f*ttj(j)));
end

% Shifted Laurent expansion, all wavenumbers
klist_l_full = klist_f;
Vfull = zeros(n, length(klist_l_full));
for i = 1:length(klist_l_full)
    Vfull(:,i) = (zj - z0).^(klist_l_full(i));
end
cl_full = Vfull\fj % solve least-squares problem for coefficients

% Laurent reconstruction
f_laurent_full = zeros(N,1);
for j = 1:N
    f_laurent_full(j) = sum(cl_full.*(exp(1i*ttj(j))-z0).^klist_l_full);
end
%f_laurent_full = Vfull*cl_full;

% Shifted Laurent expansion, only positive wavenumbers
kmax = n-1;
klist_l_pos = (0:kmax)';
Vpos = zeros(n, length(klist_l_pos));
for i = 1:length(klist_l_pos)
    Vpos(:,i) = (zj - z0).^(klist_l_pos(i));
end
cl_pos = Vpos\fj % solve least-squares problem for coefficients
f_laurent_pos = zeros(N,1);
for j = 1:N
    f_laurent_pos(j) = sum(cl_pos.*(exp(1i*ttj(j))-z0).^klist_l_pos);
end
%f_laurent_pos = Vpos*cl_pos;

if isreal(fj)
    f_fourier = real(f_fourier);
    f_laurent_full = real(f_laurent_full);
    f_laurent_pos = real(f_laurent_pos);
end

err_fourier = norm(f_fourier-ffj,inf);
err_laurent_full = norm(f_laurent_full-ffj,inf);
err_laurent_pos = norm(f_laurent_pos-ffj,inf);

fprintf('Max error Fourier reconstruction: %.2e\n', err_fourier);
fprintf('Max error Laurent (full) reconstruction: %.2e\n', err_laurent_full);
fprintf('Max error Laurent (pos) reconstruction: %.2e\n', err_laurent_pos);

figure;
semilogy(ttj,ffj, 'k');
grid on;

% Compare reconstructions
figure;
plot(ttj,ffj, 'k');
hold on;
plot(ttj,f_fourier, '.')
plot(ttj,f_laurent_full, 'x')
plot(ttj,f_laurent_pos, 'd')
legend('True f', 'Fourier', 'Laurent (full)','Laurent (pos)')
title('Function and reconstructions')
xlabel('\phi'), grid on

figure;
semilogy(klist_f,abs(ck)./max(abs(ck)),'.--');
hold on;
semilogy(klist_l_full,abs(cl_full)./max(abs(cl_full)),'.--');
semilogy(klist_l_pos,abs(cl_pos)./max(abs(cl_pos)),'.--');
legend('Fourier','Laurent (full)','Laurent (pos)');
xlabel('wavenumber, k');
ylabel('relative magnitude of coeff');
grid on;

figure;
semilogy(ttj,abs(f_laurent_full-ffj),'.');
grid on;

alignfigs;

keyboard;
end

function run_test

clear all;
close all;
format long;

test_no = 1;
n = 21;

switch test_no
    case 1
        a = 2;
        b = 0.5;
        t0 = a + 1i*b;
        sigma = 3;
        m = 3;
        f = @(t) (1 - exp(-(1 - cos(t - a)) / (2 * sigma^2))).^m;
    case 2
        a = 1.23;
        b = 1e-4;
        t0 = a + 1i*b;
        m = 3;
        f = @(t) ((exp(1i*t)-exp(1i*t0)).^m);
        f = @(t) (1-cos(t-a)).^2;
        f = @(t) ((exp(1i*t)-exp(1i*t0)).^m)+2;
    case 3
        m = 3;
        x = [-1.32;1.1;0.1];
        R = @(t) sqrt((x(1) - cos(t)).^2 + (x(2) - sin(t)).^2 + x(3).^2);
        thbar = pi/2;
        t0 = find_phi_root(x(1),x(2),x(3),thbar,@(t) 1,@(t) 1);
        a = real(t0);
        sigma = @(t) sin(t-a).^3;
        f = @(t) (sigma(t)./R(t).^m).*abs(exp(1i*t)-exp(1i*t0)).^m;
    otherwise
        error('test does not exist');
end

test_periodic_expansion(f, n, t0);

end

function phi0 = find_phi_root(x,y,z,thbar,a,b)
atilde = a(thbar).*sin(thbar);
btilde = b(thbar).*cos(thbar);
lambda = (atilde.^2+x.^2+y.^2+(btilde-z).^2)./(2*atilde*sqrt(x.^2+y.^2));
phi0 = atan2(y,x) + 1i*log(lambda+sqrt(lambda.^2-1));
end
