% tests to check the feasibility of directly interpolating the modified
% periodic basis function

% aborted, probably use recurrences instead
clear all;
close all;
clc;

m = 3;
l = 2;
kmax = 30;
kvec = (0:kmax).';
M = numel(kvec);

%integrand = @(r,t,k) (sin(t./2).^l.*cos(k.*t))./(1-2*r.*cos(t)+r.^2).^(m/2);
ell = @(r) -(4*r)./(1-r).^2;
integrand = @(r,t,k) (sin(t./2).^l.*cos(k*t))./(1-ell(r)*sin(t./2).^2).^(m/2);

chebv = cell(M,1);
warning off;
for i = 1:M
    i
    k = kvec(i);
    f = @(r) quadgk(@(t) integrand(r,t,k),0,2*pi,'MaxIntervalCount',1e6,'AbsTol',1e-14,'RelTol',1e-14);
    %f = @(r) rec_test(r,k);
    fcheb = chebfun(@(r) f(r),[0.9,1-1e-8],'splitting','off');
    chebv{i} = fcheb;
end
warning on;

%%
close all;

figure;
hold on;
kleg = cell(M,1);
for i = 1:M
    plot(chebv{i});
    grid on;
    xlabel('r');
    ylabel('val');
    kleg{i} = ['k=' num2str(kvec(i))];
end
title(['m=' num2str(m) ', l=' num2str(l)]);
legend(kleg);
xlim([0.5,1]);

for i = 1:M
    figure;
    subplot(1,2,1);
    plot(chebv{i});
    grid on;
    xlabel('r');
    ylabel('val');
    title(['k=' num2str(kvec(i))]);
    subplot(1,2,2);
    plotcoeffs(chebv{i});
    grid on;
    title(['k=' num2str(kvec(i))]);
end

alignfigs;

function muk3 = rec_test(r,kk)
kvec = 0:(kk+1);
M = length(kvec);
[K,E] = ellipke(r^2);

% p = 1/2
P(1,1) = 2*K;
P(2,1) = 2/r*(-E+K);
% p = 3/2
P(1,2) = 2/(1+r)*(2/(1+r)*E-(1-r)*K);
% p = 5/2
P(1,3) = 2/(3*(1+r)^4)*(8*(1+r^2)*E-(1-r)*(1+r)*(5+3*r^2)*K);

% Recurrence for p = 1/2
for i = 3:M
    k = kvec(i);
    P(i,1) = (1+r^2)*2*(k-1)/(2*k-1)/r*P(i-1,1) - (2*k-3)/(2*k-1)*P(i-2,1);
end

% Recurrence for p = 3/2, 5/2
plist = [3/2;5/2];
for j = 1:2
    ptmp = plist(j);
    for i = 2:M
        k = kvec(i);
        P(i,j+1) = (1+r^2)/2/r*P(i-1,j+1)-(1-r)^2*(ptmp+k-2)/(ptmp-1)/2/r*P(i-1,j);
    end
end

Ikm1 = 2*P(end-2,2)*(1-r);
Ik = 2*P(end-1,2)*(1-r);
Ikp1 = 2*P(end,2)*(1-r);
muk3 = 0.5*Ik-0.25*(Ikp1+Ikm1);
end
