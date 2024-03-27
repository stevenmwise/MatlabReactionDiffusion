function [U] = FCH_Init_SixFCirc(gamma,epsilon,eta1,eta2,N,len)
global gam bm bp
%
bm = -1;
bp = 1;
gam = gamma;
eps = epsilon;
%
% Calculate Internal Parameters
[Up,ZZ] = BilayerNumerics;
%
L = len/2.0;
%
dx = len/N;
XX = zeros(N,N);
YY = zeros(N,N);
for i = 1:N
    for j = 1:N
        XX(i,j) = -L+dx*(i-1);
        YY(i,j) = -L+dx*(j-1);
    end
end
%
U = getinit(ZZ,Up,XX,YY,eps,bm,bp,L);
%
end

function [U]=getinit(ZZ,Up,XX,YY,eps,bm,bp,L)
%
% Circle of radius 3 with 6-fold wiggle, and moderate boost
alpham = Wpp(bm);
R1 = sqrt(XX.^2+YY.^2);
theta = atan(YY./XX);
R1 = R1+eps*cos(6*theta); % wiggle the circle
R1s = min((R1-3.0)/eps,20);
R1s = max(R1s,-20);
U = interp1(ZZ,Up,R1s);
%
% Add small boost to background IC
U = U+0.5/alpham^2*eps;
%
end

