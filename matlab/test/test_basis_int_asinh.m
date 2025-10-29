clear all;
close all;
clc;

M = 100;
N = 100;

a = linspace(0.9,1.1,M).';
b = logspace(-8,0,N).';

[AA,BB] = meshgrid(a,b);
z = AA+1i*BB;

I1 = zeros(M,N);
I2 = zeros(M,N);
I3 = zeros(M,N);
Iref = vpa(zeros(M,N));
parfor i = 1:M
    disp([i M]);
    for j = 1:N
        [I1(i,j),I2(i,j),I3(i,j)] = rsqrt_pow_integrals_test(z(i));
        [~,Iref(i,j),] = rsqrt_pow_integrals_test(vpa(z(i)));
    end
end

err1 = double(norm(I1-Iref,'inf')/norm(Iref,'inf'))
err2 = double(norm(I2-Iref,'inf')/norm(Iref,'inf'))
err3 = double(norm(I3-Iref,'inf')/norm(Iref,'inf'))

E1 = abs(I1-Iref)./abs(Iref);
E2 = abs(I2-Iref)./abs(Iref);
E3 = abs(I3-Iref)./abs(Iref);

close all;

figure;
mesh(AA,log10(BB),log10(E1+eps));
colorbar;
xlabel('a'); ylabel('b'); zlabel('rel err');
title('I1');

figure;
mesh(AA,log10(BB),log10(E2+eps));
colorbar;
xlabel('a'); ylabel('b'); zlabel('rel err');
title('I2');

figure;
mesh(AA,log10(BB),log10(E3+eps));
colorbar;
xlabel('a'); ylabel('b'); zlabel('rel err');
title('I3');

alignfigs

function [I1,I2,I3] = rsqrt_pow_integrals_test(z)
NO_POWER_SERIES = false;

zr = real(z);
zi = imag(z);

u1 = sqrt((1+zr)^2 + zi^2);
u2 = sqrt((1-zr)^2 + zi^2);

% Vanilla expression
I1 = log(1-zr+u2)-log(-1-zr+u1);

if nargout == 1
    return;
end

% Evaluate after substitution zr -> -|zr|
arg2 = 1+abs(zr) + sqrt((1+abs(zr))^2 + zi^2);          
in_rhomb = 4*abs(zi) < 1-abs(zr);    
if ~in_rhomb || NO_POWER_SERIES
    arg1 = -1+abs(zr) + sqrt((-1+abs(zr))^2 + zi^2);        
else
    % Series evaluation needed inside 
    % rhombus [-1, i/4, 1, -i/4, -1].
    % Here arg1 has cancellation due to structure
    % -x + sqrt(x^2+b^2)
    Ns = 11;
    coeffs = coeffs_I1(Ns);
    arg1 = (1-abs(zr))*eval_series(coeffs, 1-abs(zr), zi, Ns);    
end   
I2 = log(arg2)-log(arg1);

% New expression
I3 = asinh((1-zr)/abs(zi)) + asinh((1+zr)/abs(zi));
end
