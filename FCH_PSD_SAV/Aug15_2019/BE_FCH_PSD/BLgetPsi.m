function [BLPsi,BLPsip]=BLevans(dlam)
global eps eta bm bp a Xbl Phibl lamb
lamb=dlam;
tf=Xbl(end);
options= odeset('Reltol',1e-7,'AbsTol', 1e-7);
%Solve for Yost functions
%Note that Psi-(t)=Psi+(-t) so Evans function is
% E = Psi+(0)Psi+'(0)
mu=sqrt(Wpp(bm)+lamb);
init=[1 mu]*exp(-tf*mu);
tspan=[-tf 0.2];
[T,Psi]=ode45('BLEvansfield',tspan,init,options);
%     BLPsi=interp1(T,Psi(:,1).*(1-tanh(5*(T-9)))/2,-Xbl);
%     BLPsip=interp1(T,Psi(:,2).*(1-tanh(5*(T-9)))/2,-Xbl);
BLPsi=interp1(T,Psi(:,1),-Xbl);
BLPsip=-interp1(T,Psi(:,2),-Xbl);

end

