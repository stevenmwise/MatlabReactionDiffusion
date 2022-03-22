function [phi] = smoothedAnnulus1(N,L,epsilon,flag)
%
M = 2^10;
%
if N > M
  error('Error in smoothedAnnulus1: N is too large for smoothing.')
end
%
phiPrime = zeros(M,M);
hPrime = L/M;
xPrime = zeros(M,M);
yPrime = zeros(M,M);
for j = 1:M
  for i = 1:M
    xPrime(i,j) = hPrime*i;
    yPrime(i,j) = hPrime*j;
  end
end
%
alphaFilt = 50*log(10);
%
alphaAlias = [0:M/2-1,M/2,-M/2+1:-1];
%    
filtX = zeros(M,M);
filtY = zeros(M,M);
for i=1:M
  filtX(i,:) = alphaAlias;
  filtY(:,i) = alphaAlias';
end
%
ita1 = 2*filtX/M;
ita2 = 2*filtY/M;
%
filtX = exp(-alphaFilt*ita1.^2);
filtY = exp(-alphaFilt*ita2.^2);
%
phi1 = zeros(M,M);           
phi2 = zeros(M,M);
%  
switch flag
  case 1 % circular ring
    a = 1.0;
    b = 1.0;
  case  2 % elliptical ring
    a = 1.0;
    b = 2.0;
  otherwise
    error("Error in smoothedAnnulus1: flag error!")
end
%
r1 = L/3;
r2 = r1*0.95;
%
phiPrime = zeros(M,M);
for i = 1:M
  for j = 1:M
    dist1 = sqrt((xPrime(i,j)-L/2)^2/a^2+(yPrime(i,j)-L/2)^2/b^2);
    phi1(i,j) = tanh((abs(dist1)-r1)/(sqrt(2.0)*epsilon));
    phi2(i,j) = tanh((abs(dist1)-r2)/(sqrt(2.0)*epsilon));
    phiPrime(i,j) = -phi1(i,j)*phi2(i,j);
  end
end
%
phiPrime = 2.0*(phiPrime-min(min(phiPrime)))/(max(max(phiPrime)) ...
  -min(min(phiPrime)))-1.0;
%
phiPrimeHat = fft2(phiPrime);
phiPrime = real(ifft2(filtX.*filtY.*phiPrimeHat));
%
phi(1:N,1:N) = phiPrime(M/N:M/N:M,M/N:M/N:M);
%
end