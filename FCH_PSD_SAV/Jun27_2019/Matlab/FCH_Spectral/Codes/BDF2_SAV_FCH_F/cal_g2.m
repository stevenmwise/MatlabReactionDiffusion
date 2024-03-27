function [g]=cal_g2(M,dt,epsilon,eta1,gamma,h,S1,S2,Ucross,Zstar,Zstar_lap,phicross,phistar_lap,phistar_lap2)

sau=sum(sum(Zstar.*phicross/3))*h*h;

g=phicross/(2*M*dt)+S1*phistar_lap2-S2*phistar_lap+Ucross/3*Zstar_lap...
    -0.5*Zstar_lap*sau+(2+eta1+2.0/12.0*gamma^2)*epsilon^2*phistar_lap2 ;

