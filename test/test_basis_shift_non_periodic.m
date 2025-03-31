% MATLAB script to demonstrate cancellation error in I_1(t0)
clear all;
close all;
clc;
format long;
rng(123);
addpath('utils');
addpath('chebfun-master');

% Define parameters
a = 0.5;   % Location where f_1 is small
b_values = logspace(-5, -1, 10);
delta = 1e-10;
m = 3;     % Order of singularity in g_1
N = 8;

% Shift basis
tshift = 0;

% Define the nearly singular function g_1
g1 = @(t, t0) 1 ./ abs(t - t0).^m;

% Define f_1(t) to be small near t = a
%f1 = @(t) (t - a).^4 + (t-a).^2; % Quadratic vanishing near a
f1 = @(t) (t-a).^4 + (t-a).^2 + delta;
%f1 = @(t) (t).^4 + (t).^2 + 1e-8;
frand = randnfun(2.8);
%f1 = @(t) (t - a).^4 + (t-a).^2 + frand(t);
%f1 = @(t) (t - a).^4.*frand(t) + (t-a).^2.*frand(t) + 1e-10;
%f1 = @(t) 1e-1*ones(size(t));

% Sample points for interpolation
xj = lgwt(N,-1,1);
%xj = lgwt(N,vpa(-1),vpa(1));
fj = f1(xj);

% Compute monomial expansion coefficients (solve Vandermonde system)
V = fliplr(vander(xj-tshift));
c = V \ fj

%c(idx) = 0;
%c = 0*c;
%c(1) = 1e-10;

% Plot error vs. b for different values
tt = linspace(-1,1,1000);
I1_refs = zeros(size(b_values));
I1_approxs = zeros(size(b_values));
Ik_values = zeros(N,length(b_values));
for i = 1:length(b_values)
    t0 = a + 1i * b_values(i);
    
    % Recompute Ik with updated troot
    Ik = zeros(N,1);
    for k = 1:N
        phi_k = @(t,tshift) (t-tshift).^(k-1);
        integrand = @(t) phi_k(t,tshift) .* g1(t, t0);
        %Ik(k) = integral(integrand, -1, 1,'AbsTol', 1e-12,'RelTol',1e-12);
        Ik1 = integral(integrand, -1, a,'AbsTol', 1e-12,'RelTol',1e-12);
        Ik2 = integral(integrand, a, 1,'AbsTol', 1e-12,'RelTol',1e-12);
        Ik(k) = Ik1 + Ik2;
    end
    
    % Compute I_1
    tmp = c .* Ik;
    Ik_values(:,i) = Ik;
    I1_approxs(i) = sum(c .* Ik);
    I1_refs(i) = integral(@(t) f1(t) .* g1(t, t0), -1, 1,'AbsTol', 1e-12,'RelTol',1e-12);

end

% Store error
abs_errors = abs(I1_approxs - I1_refs)+eps;
rel_errors = abs(I1_approxs - I1_refs) ./ abs(I1_refs)+eps;

%% Plots
close all;

% Plot error
figure;
loglog(b_values, abs_errors, '.--');
hold on;
loglog(b_values, rel_errors, 'o--');
loglog(b_values, 5*1e-16*b_values.^(-2),'k');
loglog(b_values, 1e-16*b_values.^(-4),'k');
xlabel('b (imaginary part of t0)');
ylabel('Error');
title('Cancellation error vs. b');
legend('Abs','Rel','O(1/b^2)','O(1/b^4)');
grid on;

% Plot integral value
figure;
loglog(b_values,I1_refs,'.--');
hold on;
loglog(b_values,I1_approxs,'o--');
legend('Ref','Approx');
title('Integral value');
xlabel('b (imaginary part of t0)');
ylabel('Value');
grid on;

% Plot function f1
figure;
plot(tt,f1(tt));
hold on;
xline(a,'-','t=a');
grid on;
title('Function f1');

% Plot interpolation error
kk = (1:N).'-1;
f1exact = f1(tt);
f1interp = zeros(1,numel(tt));
for i = 1:numel(tt)
    f1interp(i) = sum(c.*(tt(i)-tshift).^kk);
end
figure;
semilogy(tt,abs(f1exact-f1interp)+eps);
hold on;
semilogy(tt,abs(f1exact-f1interp)./abs(f1exact)+eps);
grid on;
legend('Abs','Rel');
xlabel('t');
ylabel('Error');
title('Error in interpolation');

% Plot Ik values
figure;
for i = 1:N
    semilogy(1:N,abs(Ik_values(:,i)),'x--');
    hold on;
end
grid on;
leg = {};
for i = 1:N
    leg{i} = num2str(b_values(i));
end
legend(leg);
title('|Ik| for different values of b');
xlabel('N');
semilogy(1:N,abs(c),'ro--');
ylabel('val');
legend('|Ik|','|ck|');

% Plot tmp vec
figure;
for i = 1:N
    plot(1:N,c.*Ik_values(:,i),'x--');
    hold on;
end
grid on;
legend(leg);
ylabel('val');
xlabel('N');
title('vector to sum');

alignfigs;
