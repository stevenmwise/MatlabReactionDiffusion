function  [u,Z]=cal_UZ(epsilon,eta2,gamma,h,B,psi,phi,phi_x,phi_y,phi_lap,phi3_lap)

F=0.5*(psi+1).^2.*(0.5*(psi-1).^2+gamma/3.0*(psi-2.0));
Fp=(psi+1).*(psi-1).*(psi+0.5*gamma);
Fpp=(psi+0.5*gamma).*(2*psi)+(psi+1).*(psi-1);


U=3*epsilon^2*phi.^2.*(phi_x.^2+phi_y.^2)+0.5*Fp.^2-eta2*F;

u=sum(sum(U))*h*h;
u=sqrt(u+B);

Z=-epsilon^2*(phi3_lap+3*phi.^2.*phi_lap)+Fpp.*Fp-eta2*Fp;
Z=Z/u;









