% This program calculates approximate solutions to the 
% Functionalized Cahn-Hilliard (FCH) equation using a first-order
% IMEX scheme, with stabilization. The method was designed by Chen
% and Shen in the paper
% 
% Feng Chen, Jie Shen, "Efficient spectral-Galerkin methods for 
%   systems of coupled second-order equations and their 
%   applications," J. Comput. Phys. 231 (15) (2012) 5016-5028.
% 
% A descprition of the model can be found in the paper
%
% C. Zhang, J. Ouyang, C. Wang, and S.M. Wise, "Numerical 
%   Comparison of Modified-Energy Stable SAV-Type Schemes 
%   and Classical BDF Methods on Benchmark Problems for the 
%   Functionalized Cahn-Hilliard Equation," J. Comput. Phys. 
%   423 (2020) 109772 (35 pages).
%
% Most of the test problems below come from that 2020 paper.
%
% initType = 1  ::  High accuracy benchmark
% initType = 2  ::  phiAve + random noise
% initType = 3  ::  S shape
%	initType = 4  ::  S shape with smoothing filter
%	initType = 5  ::  annulus 1 with smoothing filter
%	initType = 6  ::  annulus 2 with smoothing filter
%
clear;
clc;
clf;
initType = 2;
%
% Chen-Shen stabilization parameters:
sig1 = 1.0;
sig2 = 1.0;
%
% Time stepping and printout parameters:
dt = 1e-4;
stepsPerPlot = 10000;
stepsPerScreenPlot = 1000;
stepsPerReport = 100;
plotFrames = 100;
maxSteps = stepsPerPlot*plotFrames;
phiAve = 0.0;
%
switch initType
  case 1
    M = 1.0;
    L = 2*pi; 
    N = 2^7;
    epsilon = 0.18;
    eta1 = epsilon^2;
    eta2 = epsilon^2;
    gamma = 0.00;
  case 2
    M = 1.0;
    L = 2*pi; 
    N = 2^8;
    epsilon = 0.04;
    eta1 = 5.0*epsilon;
    eta2 = 3.0*epsilon;
    gamma = 0.0;
    phiAve = 0.5;
  case {3, 4}
    M = 1.0;
    L = 12.8;
    N = 2^8;
    epsilon = 0.1;
    eta1 = 0.2;
    eta2 = 0.2;
    gamma = 0.00;
  case 5
    M = 1.0;
    L = 2*pi;
    N = 2^8;
    epsilon = 0.04;
    eta1 = 2.0*epsilon;
    eta2 = 2.0*epsilon;
    gamma = 0.25;
  case 6
    M = 1.0;
    L = 4*pi;
    N = 2^8;
    epsilon = 0.1;
    eta1 = 1.45*epsilon;
    eta2 = 2.0*epsilon;
    gamma = 0.25;
  otherwise
    error("Error in initType: No such initial conditon!")
end
%
param.L        = L;
param.N        = N;
param.epsilon  = epsilon;
param.eta1     = eta1;
param.eta2     = eta2;
param.gamma    = gamma;
param.phiAve   = phiAve;
param.initType = initType;
%
% Grid and laplacian matrices:
kx = 2.0*pi*[0:N/2-1 N/2 -N/2+1:-1]/L;
ky = 2.0*pi*[0:N/2-1 N/2 -N/2+1:-1]/L;
kx2 = kx.*kx;
ky2 = ky.*ky;
%
h = L/N;
xx  = zeros(N,N);
yy  = zeros(N,N);
lap = zeros(N,N);
%
for j = 1:N
  for i = 1:N
    xx(i,j) = h*i;
    yy(i,j) = h*j;
    lap(i,j) = -(kx2(i)+ky2(j));
  end
end
%
% Get initial condition
[phi] = initialCondition(param);
%
% Print out initial frame:
%
printField('phi',phi,0,dt,0,param)
%
figure(1);
pcolor(xx,yy,phi);
axis equal
shading interp;
colormap('jet');
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
  fp  = (phi-1.0).*(phi+1.0).*(phi+0.5*gamma);
  fpp = (phi+1.0).*(phi+0.5*gamma)+(phi-1.0).*(phi+0.5*gamma) ...
    +(phi-1.0).*(phi+1.0);
%
  wt = epsilon*epsilon*real(ifft2(fft2(phi).*lap))-fp;
  q1 = (sig1-fpp).*wt+(eta1-eta2)*fp;
  q2 = sig2*phi-fp;
  q3 = epsilon*epsilon*real(ifft2(fft2(q2).*lap))-(sig1-eta1)*q2+q1;
  q4 = phi+dt*M*real(ifft2(fft2(q3).*lap));
%
  coef = 1.0-dt*M*epsilon^4*lap.^3 ...
    + dt*M*epsilon^2*(sig2+sig1-eta1)*lap.^2 ...
    - dt*M*(sig2*(sig1-eta1))*lap;
%
  phi = real(ifft2(fft2(q4)./coef));
%
  if (mod(k,stepsPerScreenPlot) == 0)
    figure(1);
    pcolor(xx,yy,phi);
    axis equal
    shading interp; colormap(jet);
    colorbar;
    title(['FCH-IMEX: dt = ',num2str(dt),'  time = ', ...
      num2str(time)]);
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


