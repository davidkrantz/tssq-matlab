clear all;
close all;
clc;
format long;

n = 16;

use_vpa = 0;

Na = 100;
Nb = 100;

av = linspace(-1.5,1.5,Na);
bv = linspace(-1.5,1.5,Nb);

av = linspace(1e-5,1.5,Na);
bv = linspace(1e-5,1.5,Nb);

[A,B] = meshgrid(av,bv);
T0 = A+1i*B;

[P1,P3,P5,P1sh,P3sh,P5sh] = deal(nan*zeros(Na,Nb,n));
for i = 1:Na
    for j = 1:Nb
        t0 = T0(i,j);
        if use_vpa
            [P1(i,j,:),P3(i,j,:),P5(i,j,:)] = rsqrt_pow_integrals(vpa(t0), n);
            [P1sh(i,j,:),P3sh(i,j,:),P5sh(i,j,:)] = rsqrt_pow_integrals_shift(vpa(t0), n);
        else
            [P1(i,j,:),P3(i,j,:),P5(i,j,:)] = rsqrt_pow_integrals(t0, n);
            [P1sh(i,j,:),P3sh(i,j,:),P5sh(i,j,:)] = rsqrt_pow_integrals_shift(t0, n);
        end
    end
end

%%
m = 3;
k = 4;

F = chebfun2(@(a,b) myfunc(a,b,m,k),[1e-5,1.5,1e-5,1.5])

close all;
figure;
plot(F);
figure;
plotcoeffs(F);

%% Plot
close all;
for i = 1:n
    figure;
    surf(A,B,log10(abs(P3sh(:,:,i))));
    colorbar;
    shading interp;
    hold on;
    plot3([-1 1],[0 0],1e100*[1 1],'k','linewidth',1.5);
    title(['k=' num2str(i)]);
    xlabel('a');
    ylabel('b');
    xlim([min(av) max(av)]);
    ylim([min(bv) max(bv)]);
    view(0,90);
end

alignfigs;

%%
close all;

figure;
for i = 2:n
    loglog(bv,(abs(P3sh(:,1,i))))
    hold on;
end
grid on;

function Q = myfunc(a,b,m,k)
t0 = a+1i*b;
[Q1,Q3,Q5] = rsqrt_pow_integrals_shift(t0,k);

if m == 1
    Q = Q1(end);
elseif m == 3
    Q = Q3(end);
else
    Q = Q5(end);
end

Q = log10(abs(Q));

end
