function [Ub,z,data] = BilayerNumerics
%
global Xbl Phibl Phiblp tag tagwell bllam0 delta Psi0Bl Psi0Blp gam bm bp
%
% Initial parameters
Xmax = 20;  % fix max x value
delta = 0.01;
Xbl = [1e-10:delta:Xmax];
Xbl = [Xbl Xmax];
%
% Obtain Bilayer Profile (Phibl is Ub)%
MakeBilayer; % return solution in Phibl & Phiblp
Intcheck = trapz(Xbl,1/2*(Phiblp).^2-W(Phibl));
if (abs(Intcheck)>1e-6)
   sprintf('Intcheck>1e-6 in Bilayer Profile')
   Intcheck
end

Mbl = 2*trapz(Xbl,Phibl-bm);
data(1) = Mbl;
B1 = 2*trapz(Xbl,(Phibl-bm).^2);
data(2) = B1;
sigmab = 2*trapz(Xbl,(Phiblp).^2);
data(3) = sigmab;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Bilayer ground state eigenvalue and eigenfunction %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% fix range for fzero to look for lambda
%lambdaRange=[abs(gam)/1.2 abs(gam)*1.2]; %good for Keith's W
lambdaRange=[0.005 15]; %good for Arjen's W
options=optimset('TolX',1e-4);
bllam0=fzero('BLEvans', lambdaRange ,options);
data(4)=bllam0;
%
% Find the ground eigenfunction
[Psi0Bl,Psi0Blp] = BLgetPsi(bllam0);
Psi0Bl2p = CalcDerivative(-Psi0Blp, delta);
%
% Verify solution
Psi0Resid = trapz(Xbl,Psi0Bl2p-Wpp(Phibl).*Psi0Bl-bllam0.*Psi0Bl);
if (abs(Psi0Resid)>1e-4)
   sprintf('Psi0 Residual>1e-4 in Bilayer Psi0, %0.5g',Psi0Resid)
end
%
% Normalize eigenfunc
Psi0Bl = Psi0Bl/sqrt((2*trapz(Xbl,Psi0Bl.^2)));
Psi0Blp = Psi0Blp/sqrt((2*trapz(Xbl,Psi0Bl.^2)));
Psi0Bl2p = Psi0Bl2p/sqrt((2*trapz(Xbl,Psi0Bl.^2)));


%These are GLOBAL variables do not delete
%Ubl2=2*trapz(Xbl,Phibl.^2);
%Ubpl2=2*trapz(Xbl, Phiblp.^2);
%Psibll2=2*trapz(Xbl,Psi0Bl.^2);
%Psiblpl2=2*pi*trapz(Xbl,Psi0Blp.^2);
%rhobl=Ubpl2/Mbl;
%End define global variables

%
% Calc S:=\int \PSI1W'''(U_b)\psi_0\,dz
z = [fliplr(-Xbl) Xbl];
Ub = [fliplr(Phibl), Phibl];
[Phi1,BlPhi1p,BlPhi12p] = BlPhi1;
%
%   Verify Phi_1
% VerifyPhi1 = norm(BlPhi12p-(Wpp(Phibl).*Phi1)-1);
Intcheck = trapz(Xbl,BlPhi12p-(Wpp(Phibl).*Phi1)-1);
if (abs(Intcheck)>5e-5)
   sprintf('Intcheck>5e-5 in Bilayer Phi1')
   Intcheck
end
%
Sbl = 2*trapz(Xbl, Psi0Bl.^2.*Wppp(Phibl).*Phi1);
data(5) = Sbl;
end
