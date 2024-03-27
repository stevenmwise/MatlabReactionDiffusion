function [PSI1,PSI1Blp,PSI1Bl2p] = BlPhi1

% SOLVE  d2x/dt2+W''(U_b)x = 1
% initial conditions: x(0) = 0, x'(0)=0

global bm bp Xbl Phibl Phiblp Psi0Bl Psi0Blp delta

%close all

epsilon=1;%e-10;
clear solinit;
%Phimax=fzero('W',[0.9*bm -b]);
%options= odeset('Reltol',epsilon,'AbsTol',epsilon);
tspan=Xbl;
%tf=Xbl(end);
%mu=sqrt(Wpp(b));
%init=[1 mu]*exp(-tf*mu);%[Phimax 0];
% chi=transpose(1-(1-tanh((tspan-12)*7))/2);% cut-off function to make guess smooth

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up the initial guess %
%
BVPInitGuess();
%solinit = bvpinit(Xbl,@BVPInitGuess);
%figure
%plot(solinit.x, solinit.y(1,:))
%hold on
%plot(solinit.x, solinit.y(2,:), 'r')
%txt=sprintf('Initial Guess');
%title(txt)
%xlabel('z')
%ylabel('Y1,Y2')

% solinit.x=tspan;
% solinit.y(1,:)=Y(:,1).*chi+b.*(1-chi);
% solinit.y(2,:)=Y(:,2).*chi;
options=bvpset('RelTol',epsilon);
sol=bvp4c('PhiBlODE','PhiBlBC',solinit,options);
PSI1= interp1(sol.x,sol.y(1,:),Xbl)-1/Wpp(bm);%sol.y(1,:);
PSI1Blp=interp1(sol.x,sol.y(2,:),Xbl);
temp=(diff(sol.y(2,:))/0.01);
PSI1Bl2p=CalcDerivative(PSI1Blp,delta);
PSI1Bl2p = interp1(sol.x,PSI1Bl2p,Xbl);

% figure
% plot(Xbl, PSI1,'Linewidth',3)
% hold on
% % plot(sol.x, PSI1Blp,'r','Linewidth',3)
% % plot(sol.x, PSI1Bl2p,'k','Linewidth',3)
% 
% txt=sprintf( 'Bilayer \\Phi_1^b function');
% title(txt,'FontSize', 20)
% xlabel('z','FontSize', 16)
% ylabel('\Phi_1','FontSize', 16)
% print -depsc BlPhi1.eps
  

    function InitGuess=BVPInitGuess()
        G =fliplr( InitProfile(tspan));
        Gp_res=fliplr(Gp(tspan));
        InitGuess = [G Gp_res];
        solinit.x = Xbl;
        solinit.y(1,:) = Psi0Bl;%Phibl+1;%G;
        solinit.y(2,:)=-Psi0Blp;%Phiblp; %Gp_res;

    end

    function G = InitProfile(u)
       MaxValPhi = 1000;
G = exp(Wpp(bm).*(u/2));
MaxInd = length(u);
for Ind=1:1:MaxInd
    if (G(Ind) >= 688)
        G(Ind) = (1-(1-tanh((u(Ind)-9.5)))/2).*MaxValPhi;
    end
end
% figure
%  plot(Xbl,G,'b');
    end


    function Gp_res=Gp(u)
        MaxValPhi = 1000;
        del=1e-9;
%         G = @(x) exp(Wpp(b).*(x/2));
%         chi = @(x)(1-(1-tanh((x-9.5)))/2).*MaxValPhi;
%         Lp=(G(u+del)-G(u-del))/2/del;
%         MaxInd = length(u);
%         for Ind=1:1:MaxInd
%             if (Lp(Ind) >= 688)
%                 Lp(Ind) = (chi(u(Ind)+del)-chi(u(Ind)-del))/2/del;
%             end
%         end
%         
%         chip = (chi(u+del)-chi(u-del))/2/del;
        Gp_res=(InitProfile(u+del)-InitProfile(u-del))/2/del;
%         for i=1:1:length(u)
%         if (InitProfile(u(i)+del) < 670)
%             Gp_res(i)=(InitProfile(u(i)+del)-InitProfile(u(i)-del))/2/del;
%         else
%             Gp_res(i)=0;
%         end
%         end
    end

end

