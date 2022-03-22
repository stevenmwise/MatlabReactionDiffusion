function [phi] = smoothedAnnulus2(N,L,stretch)
%
M = 2^10;
%
if N > M
  error('Error in smoothedAnnulus2: N is too large for smoothing.')
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
phiPrime = zeros(M,M);
for i=1:M
  for j=1:M
    r = sqrt((xPrime(i,j)-L/2)^2+((yPrime(i,j)-L/2)^2)/stretch);
    if     (r > L/4+0.2)
      phiPrime(i,j) = -1.0;
    elseif (r < L/4-0.2)
      phiPrime(i,j) = -1.0;
    else
      phiPrime(i,j) =  1.0;
    end
  end
end
%
phiPrimeHat = fft2(phiPrime);
phiPrime = real(ifft2(filtX.*filtY.*phiPrimeHat));
%
phi(1:N,1:N) = phiPrime(M/N:M/N:M,M/N:M/N:M);
%
end