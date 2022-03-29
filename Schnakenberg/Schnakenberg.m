% This program calculates approximate solutions to the Schnakenberg
% system in 2D assuming periodic boundary conditions. We use 
% a simple backward Euler IMEX scheme. Space is discretized via the
% pseudo-spectral method. This system experiences Turing instability
% whenever
%
%
%
% See the reaction terms Schnak1 and Schnak2 below.
%
% In 1D the system reads
%
%     phi_t = D1*(phi_xx+phi_yy)+Schnak1(phi,psi),
%
%     psi_t = D2*(psi_xx+psi_yy)+Schnak2(phi,psi),
%
% where
%
%     Schnak1(phi,psi) =  phi.*phi.*psi-phi+(b-a)/2.0,
%
%     Schnak2(phi,psi) = -phi.*phi.*psi+(b+a)/2.0.
%
% See the paper
%
% H.K. Khudhair, Y. Zhang, N. Fukawa, "Pattern Selection in the
%   Schnakenberg equations: From normal to anomalous diffusion,"
%   Numer. Method. Partial Diff. Eq. (2021).
%   https://doi.org/10.1002/num.22842
%
% 
% for more details. The IMEX scheme employed here is
%
%     phi^{k+1} - phi^k = dt*(D1*phi^{k+1}_xx
%                        + Schnak1(phi^k,psi^k)),
%
%     psi^{k+1} - psi^k = dt*(D2*psi^{k+1}_xx
%                        + Schnak2(phi^{k+1},psi^k)).
%
%
clear;
clc;

dt = 1.0e-02;
stepsPerPlot = 500;
stepsPerScreenPlot = 50;
stepsPerReport = 100;
plotFrames = 50;
maxSteps = stepsPerPlot*plotFrames;

N  = 256;
D1 = 1.000e-02;
D2 = 1.000e+00;
a  = 3.200e+00;
b  = 3.500e+00;
L  = 8.000e+00;

% Grid and laplacian matrices:
kx = 2.0*pi*[0:N/2-1 N/2 -N/2+1:-1]/L;
ky = 2.0*pi*[0:N/2-1 N/2 -N/2+1:-1]/L;
kx2 = kx.*kx;
ky2 = ky.*ky;

h = L/N;
xx  = zeros(N,N);
yy  = zeros(N,N);
lap = zeros(N,N);
for j = 1:N
  for i = 1:N
    xx(i,j) = h*i;
    yy(i,j) = h*j;
    lap(i,j) = -(kx2(i)+ky2(j));
  end
end
%
% Parameters:

param.D1 = D1;
param.D2 = D2;
param.a = a;
param.b = b;
param.L  = L;
param.N  = N;
%
% Initialization:
phi = zeros(N,N);
psi = zeros(N,N);
for i = 1:N
  for j = 1:N
%
    phi(i,j) = b;
    psi(i,j) = (a+b)/(2.0*b*b)+0.01*(rand-0.5);
%
  end
end
%
% Print out initial frame:
%
printField('phi',phi,0,dt,0,param)
%
figure(1);
pcolor(xx,yy,phi);
axis equal;
shading interp;
colormap(jet);
title('Initial Conditions');
colorbar;
getframe;

for k = 1:maxSteps
%
  time = k*dt;
%
  if (mod(k,stepsPerReport) == 0)
    fprintf('step : %6d       time : %8.3f \n', k, time)
  end
%
  q = phi+dt*Schnak1(phi,psi,param);
%
  coef = 1.0-dt*D1*lap;
%
  phi = real(ifft2(fft2(q)./coef));
%
  q = psi+dt*Schnak2(phi,psi,param);
%
  coef = 1.0-dt*D2*lap;
%
  psi = real(ifft2(fft2(q)./coef));
%
  if (mod(k,stepsPerScreenPlot) == 0)
    figure(1);
    pcolor(xx,yy,phi);
    axis equal
    shading interp; 
    colormap(jet);
    colorbar;
    title(['Schnakenberg-IMEX: dt = ',num2str(dt), ...
      '  time = ',num2str(time)]);
    getframe;
  end
%
  if (mod(k,stepsPerPlot) == 0)
%
    frame = round(k/stepsPerPlot);
    printField('phi',phi,frame,dt,time,param)
%
  end
end
%
disp('program done')
%
% Embedded function(s) below:
%
function [fphi] = Schnak1(phi,psi,param)
%
a = param.a;
b = param.b;
%
fphi =  phi.*phi.*psi-phi+(b-a)/2.0;
%
end % function Schnak1
%
function [fphi] = Schnak2(phi,psi,param)
%
a = param.a;
b = param.b;
%
fphi = -phi.*phi.*psi    +(b+a)/2.0;
%
end % function Schnak2
%
function [ ] = printField(fieldName,fieldArray,frame,dt,time,param)
%
N = param.N;
L = param.L;
%
s1 = ['0000000' num2str(frame)];
s2 = s1((length(s1)-4):length(s1));
fid = fopen(['./OUT/',fieldName,s2,'.dat'],'w');
fprintf(fid,'%25.15e %25.15e %10i %25.15e\n',time,dt,N,L);
for j = 1:N
  for i = 1:N
    fprintf(fid,'%25.15e\n',fieldArray(i,j));
  end
end
fclose(fid);
%
end % function printField

