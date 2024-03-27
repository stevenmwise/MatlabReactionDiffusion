% Solves  Y'' = W'(Y)
%   Y(\pm infty)=b
%  By integrating the ODE from x=0
% with initial data Y(0)=Ym, Y'(0)=0
% where W(Ym)=0
function MakeBilayer;
 global bm bp gam Xbl Phibl Phiblp tag tagbl
 Phimax=fzero('W',[0.9*bm bp]); % For proper choice of W these bound zero
 epsilon=1e-10;
 clear solinit;
 options= odeset('Reltol',epsilon,'AbsTol',epsilon);
 tspan=Xbl;
 init=[Phimax 0];
 [T,Y]=ode45('Wblfield', tspan, init, options);
 Y1=Y(:,1);
 Y2=Y(:,2);
 
 temp= (diff(Y(:,2))/0.01);
 Y3=transpose([transpose(temp), temp(end)]);
 chi=(1-tanh((T-12)*7))/2;
 Phibl=interp1(T,bm+(Y(:,1)-bm).*chi,Xbl, 'linear');
 Phiblp=interp1(T,Y(:,2).*chi,Xbl, 'linear');
 %%%%%%%%%%%%%%%%%%%%%%%%%
 %  Verify
 %%%%%%%%%%%%%%%%%%%%%%%%%
 Phiblpp = interp1(T,Y3.*chi,Xbl);
%norm(Phiblpp-Wp(Phibl))
 
end %compute Phibl
