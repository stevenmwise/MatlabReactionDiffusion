function [phi] = initialCondition(N,xx,yy,L,ratio)
%
phi = zeros(N,N);
%
for j = 1:N
  for i = 1:N
    if (xx(i,j) <= L/2.0);
      phi(i,j) = -0.38;
    else
      phi(i,j) = -0.28+(0.050*rand-0.025);
    end
  end
end
%
end