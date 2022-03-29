% This program calculates approximate solutions to the Cahn-Hilliard 
% equation in 2D assuming periodic boundary conditions. We use a 
% simple backward Euler IMEX scheme. Space is discretized via the
% pseudo-spectral method. 
%
% In 1D, the equation reads
%
%     phi_t = (dF(phi)-eps2*phi_xx)_xx.
%
% The 2D version is analogous. The IMEX scheme is
%
%     phi^{k+1} - phi^k = dt*(dF(phi^k)-eps2*phi^{k+1}_xx)_xx,
%
% or, in other words,
%
%     phi^{k+1}+dt*eps2*phi^{k+1}_xxxx = phi^k+dt(dF(phi^k))_xx
%
%
clear;
clc;
clf;
%
dt = 1.0e-04;
stepsPerPlot = 100;
stepsPerScreenPlot = 10;
stepsPerReport = 100;
plotFrames = 100;
maxSteps = stepsPerPlot*plotFrames;
%
N = 512;
eps2 = 1.000e-05;
ave = 3.000e-02;
%
L = 1.0;
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
param.eps2 = eps2;
param.ave  = ave;
param.L    = L;
param.N    = N;
%
% Initialization:
phi = zeros(N,N);
for i = 1:N
  for j = 1:N
%
%    phi(i,j) = ave+0.05*(rand-0.5);
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
s1 = ['0000000' num2str(0)];
s2 = s1((length(s1)-4):length(s1));
fid = fopen(['./OUT/phi',s2,'.dat'],'w');
fprintf(fid,'%25.15e %25.15e %10i %25.15e\n',0,dt,N,L);
for j = 1:N
  for i = 1:N
    fprintf(fid,'%25.15e\n',phi(i,j));
  end
end
fclose(fid);
%
figure(1);
pcolor(xx,yy,phi);
shading interp; colormap(jet);
title('Initial Conditions');
colorbar;
getframe;
axis equal
%
for k = 1:maxSteps
%
  time = k*dt;
%
  if (mod(k,stepsPerReport) == 0)
    fprintf('step : %6d       time : %8.3f \n', k, time)
  end
%
  fphi = dF(phi,param);
  q = phi+dt*real(ifft2(fft2(fphi).*lap));
%
  coef = 1.0+dt*eps2*lap.*lap;
%
  phi = real(ifft2(fft2(q)./coef));
%
  if (mod(k,stepsPerScreenPlot) == 0)
    figure(1);
    pcolor(xx,yy,phi);
    axis equal
    shading interp; colormap(jet);
    colorbar;
    title(['CH-IMEX: dt = ',num2str(dt),'  time = ',num2str(time)]);
    getframe;
  end
%
  if (mod(k,stepsPerPlot) == 0)
%
    s1 = ['0000000' num2str(round(k/stepsPerPlot))];
    s2 = s1((length(s1)-4):length(s1));
    fid = fopen(['./OUT/phi',s2,'.dat'],'w');
    fprintf(fid,'%25.15e %25.15e %10i %25.15e\n',time,dt,N,L);
    for j = 1:N
      for i = 1:N
        fprintf(fid,'%25.15e\n',phi(i,j));
      end
    end
    fclose(fid);
  end
end
%
disp('program done')
%
% Embedded function(s) below:
%
function [fphi] = dF(phi,param)
%
% Derivative of the homogeneous free energy density.
%
fphi = phi.*(phi-0.5).*(phi-1.0);
%
end % function dF


