function [g]=cal_g1(epsilon,eta1,gamma,M,dt,h,S1,S2,U,Z,Z_lap,phi,phi_lap,phi_lap2)

sau=sum(sum(Z.*phi))*h*h;

g=phi/(M*dt)+S1*phi_lap2-S2*phi_lap+U*Z_lap-0.5*Z_lap*sau+(2+eta1+2.0/12.0*gamma^2)*epsilon^2*phi_lap2;
