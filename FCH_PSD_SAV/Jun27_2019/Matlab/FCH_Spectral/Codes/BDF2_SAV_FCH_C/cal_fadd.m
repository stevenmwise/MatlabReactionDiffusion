function [phie,fadd]=cal_fadd(M,epsilon,eta,A,C,xx,yy,time,lap)

phie=A*sin(C*xx).*cos(C*yy)*sin(time);
phie_t=A*sin(C*xx).*cos(C*yy)*cos(time);
phie_xx=-A*C^2*sin(C*xx).*cos(C*yy)*sin(time);
phie_yy=-A*C^2*sin(C*xx).*cos(C*yy)*sin(time);

fpe=phie.^3-phie;
fppe=3*phie.^2-1;
mu=epsilon^2*(phie_xx+phie_yy)-fpe;
mue=epsilon^2*real(ifft2(lap.*fft2(mu)))+(eta-fppe).*mu;

fadd=phie_t-M*real(ifft2(lap.*fft2(mue)));
