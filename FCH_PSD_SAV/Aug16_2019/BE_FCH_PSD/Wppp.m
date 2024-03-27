function Wppp=Wppp(u)
  global eps eta a bp bm gam delta
%   epsilon=1e-4;
%   Wppp=(Wpp(u+epsilon)-Wpp(u-epsilon))/(2*epsilon);
% %   p =polyfit(u,Wppp,1)
% Wppp = CalcDerivative(Wpp(u),delta);
%    Wppp=6*u;
     Wppp=6*u-3*(bp+bm)+gam;
end
