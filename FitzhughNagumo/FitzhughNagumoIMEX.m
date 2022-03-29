% This program calculates approximate solutions to the Fitzhugh-
% Nagumo equation in 2D assuming periodic boundary conditions. We use 
% a simple backward Euler IMEX scheme. Space is discretized via the
% pseudo-spectral method. 
%
% In 1D the equation reads
%
%     phi_t = eps2*phi_xx-Allee(phi),
%
% The 2D version is analogous. The IMEX scheme is
%
%     phi^{k+1} - phi^k = dt*(eps2*phi^{k+1}_xx-Allee(phi^k)),
%
% or, in other words,
%
%     phi^{k+1}-dt*eps2*phi^{k+1}_xx = phi^k-dt*Allee(phi^k)
%
%
clear;
clc;
clf;
%
dt = 1.0e-01;
stepsPerPlot = 500;
stepsPerScreenPlot = 100;
stepsPerReport = 50;
plotFrames = 50;
maxSteps = stepsPerPlot*plotFrames;

N = 256;
eps2 = 2.000e-05;
mid  = 5.000e-01;
L = 1.0;

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
param.eps2 = eps2;
param.mid  = mid;
param.L    = L;
param.N    = N;
%
% Initialization:
phi = zeros(N,N);
for i = 1:N
  for j = 1:N
%
%    phi(i,j) = mid+0.05*(rand-0.5);
%
    if (xx(i,j)>0.15 && xx(i,j)<0.25) && ...
        (yy(i,j)>0.15 && yy(i,j)<0.85)
      phi(i,j) = 1.0;
    end
    if (xx(i,j)>0.75 && xx(i,j)<0.85) && ...
        (yy(i,j)>0.15 && yy(i,j)<0.85)
      phi(i,j) = 1.0;
    end
    if (xx(i,j)>0.15 && xx(i,j)<0.85) && ...
        (yy(i,j)>0.75 && yy(i,j)<0.85)
      phi(i,j) = 1.0;
    end
    if (xx(i,j)>0.35 && xx(i,j)<0.85) && ...
        (yy(i,j)>0.15 && yy(i,j)<0.25)
      phi(i,j) = 1.0;
    end
    if (xx(i,j)>0.35 && xx(i,j)<0.45) && ...
        (yy(i,j)>0.15 && yy(i,j)<0.65)
      phi(i,j) = 1.0;
    end
    if (xx(i,j)>0.35 && xx(i,j)<0.65) && ...
        (yy(i,j)>0.55 && yy(i,j)<0.65)
      phi(i,j) = 1.0;
    end
    if (xx(i,j)>0.55 && xx(i,j)<0.65) && ...
        (yy(i,j)>0.35 && yy(i,j)<0.65)
      phi(i,j) = 1.0;
    end
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
%
for k = 1:maxSteps
%
  time = k*dt;
%
  if (mod(k,stepsPerReport) == 0)
    fprintf('step : %6d       time : %8.3f \n', k, time)
  end
%
  q = phi-dt*Allee(phi,param);
%
  coef = 1.0-dt*eps2*lap;
%
  phi = real(ifft2(fft2(q)./coef));
%
  if (mod(k,stepsPerScreenPlot) == 0)
    figure(1);
    pcolor(xx,yy,phi);
    axis equal
    shading interp; 
    colormap(jet);
    colorbar;
    title(['FN-IMEX: dt = ',num2str(dt),'  time = ',num2str(time)]);
    getframe;
  end
%
  if (mod(k,stepsPerPlot) == 0)
%
    frame = round(k/stepsPerPlot);
    printField('phi',phi,frame,dt,time,param)
%
  end
%
end
%
disp('program done')
%
% Embedded function(s) below:
%
function [fphi] = Allee(phi,param)
%
% The Allee function reaction term.
%
mid = param.mid;
%
fphi = phi.*(phi-mid).*(phi-1.0);
%
end % function Allee
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

