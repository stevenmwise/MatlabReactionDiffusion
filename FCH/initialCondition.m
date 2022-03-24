function [phi] = initialCondition(param)
%
N = param.N;
L = param.L;
epsilon = param.epsilon;
phiAve = param.phiAve;
initType = param.initType;
%
h = L/N;
xx = zeros(N,N);
yy = zeros(N,N);
for j = 1:N
  for i = 1:N
    xx(i,j) = h*i;
    yy(i,j) = h*j;
  end
end
%
phi = zeros(N,N);
%
switch initType
  case 1
%
    phi = 2.0*exp(sin(xx)+sin(yy)-2) ...
      + 2.2*exp(-sin(xx)-sin(yy)-2.0)-1.0;
%
  case 2
%
    phi = phiAve+0.001*(2*rand(N,N)-1);
%
  case 3
%
    ro = 0.34;
    for i = 1:N
      for j = 1:N
        if     (xx(i,j) > sin(4.0*pi*yy(i,j)/L)+L/2+ro)
          phi(i,j) = -1.0;
        elseif (xx(i,j) < sin(4.0*pi*yy(i,j)/L)+L/2-ro)
          phi(i,j) = -1.0;
        else
          phi(i,j) = 1.0;
            end
        end
    end
%
  case 4
%
    ro = 0.34;
    [phi] = smoothedSCurve(N,L,ro);
%
  case 5
%
    flag = 1;
    [phi] = smoothedAnnulus1(N,L,epsilon,flag);
%
  case 6
%
    stretch = 2.0;
    [phi] = smoothedAnnulus2(N,L,stretch);
%
  otherwise
%
    Error("Error in initialCondition: No such initial conditon!")
%
end