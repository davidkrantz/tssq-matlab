%Given c-coeffs for regular Fourier basis, want to determine b-coeffs, plus a0 and a1. 
%Reg sum is 
%sum from -(N+1) to (N+1) c_k exp(ikx)
%Modified basis is
%a0+a1*sin(x-a)+sum from -N to N b_k exp(ikx) sin((x-a)/2).^2

N=128; 
a=pi/11; 
%%Shift in basis. 
w=exp(1i*a);
winv=exp(-1i*a);

kvec_c=-(N+1):(N+1);
%k incdecies for c-coeffs
kvec_b=-N:N;
%k incdecies for b-coeffs



Np=256; 
xv=(0:Np-1)/Np*2*pi; 
c_coeffs=(rand(size(kvec_c))-0.5)*2+1i*(rand(size(kvec_c))-0.5)*2;
c_coeffs=c_coeffs.*exp(-0.01*kvec_c.^2);
%c_coeffs=zeros(size(kvec_c)); 
%c_coeffs(N+2)=1; 
%c_coeffs(N)=1; 

Qsum=zeros(size(xv)); 
for l=1:length(kvec_c)
    Qsum=Qsum+c_coeffs(l)*exp(1i*kvec_c(l)*xv); 
end;

b_coeffs=zeros(size(kvec_b)); 
indB=2*N+1; 
indC=2*N+3; 
%b_N:
b_coeffs(indB)=-4*w*c_coeffs(indC); 
%b_{N-1}
b_coeffs(indB-1)=4*w*(0.5*b_coeffs(indB)-c_coeffs(indC-1));

cind=indC-2; 
%Set b_{N-2} down to b_1
for kind=indB-2:-1:N+2
    b_coeffs(kind)=4*w*(0.5*b_coeffs(kind+1)-0.25*w*b_coeffs(kind+2)-c_coeffs(cind));
    cind=cind-1; 
end;

%b_{-N}:
b_coeffs(1)=-4*winv*c_coeffs(1); 
%b_{-N+1}
b_coeffs(2)=4*winv*(0.5*b_coeffs(1)-c_coeffs(2));
cind=3; 
%Set b_{-N+2} up to b_{-1}
for kind=3:N
    b_coeffs(kind)=4*winv*(0.5*b_coeffs(kind-1)-0.25*winv*b_coeffs(kind-2)-c_coeffs(cind));
    cind=cind+1; 
end;

%d1=-0.5*b_1+0.25*w*b_2+c_1
d1=-0.5*b_coeffs(N+2)+0.25*w*b_coeffs(N+3)+c_coeffs(N+3); 
%d2=0.25*winv*b_{-1}+0.25*w*b_1+c_0
d2=0.25*winv*b_coeffs(N)+0.25*w*b_coeffs(N+2)+c_coeffs(N+2); 
%d3=-0.25*b_{-1}+0.25*winv*b_{-2}+c_{-1}
d3=-0.5*b_coeffs(N)+0.25*winv*b_coeffs(N-1)+c_coeffs(N+1); 

%b_0
b_coeffs(N+1)=-2*(w*d1+winv*d3);
%a_1
a1=1i*(w*d1-winv*d3); 
%a_0=d_2-0.5*b_0
a0=d2-0.5*b_coeffs(N+1);


%%Now, add up to Psum
Psum=a0*ones(size(xv))+a1*sin(xv-a);
for l=1:length(kvec_b)
    Psum=Psum+b_coeffs(l)*exp(1i*kvec_b(l)*xv).*sin((xv-a)/2).^2; 
end;


clf;
subplot(211)
plot(xv,real(Qsum)); 
hold on; 
plot(xv,imag(Qsum)); 
plot(xv,real(Psum),'--'); 
plot(xv,imag(Psum),'--'); 
legend('Qsum - real','Qsum - imag','Psum - real','Psum - imag'); 

subplot(212)
semilogy(xv,abs(real(Qsum-Psum))); 
hold on; 
semilogy(xv,abs(imag(Qsum-Psum))); 
legend('Diff real part','Diff imag part'); 
