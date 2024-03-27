%   This program calculates the solution to the PFC equation using the 
%     Backward Euler scheme with the PSD solver.
%
%   index == 1        ::  lap = [0:N/2-1 N/2 -N/2+1:-1];
%   index /= 1        ::  lap = [0:N/2-1   0 -N/2+1:-1];
%--------------------------------------------------------------------------
clear;
clc;
index = 1;
%--------------------------------------------------------------------------
%
% Set this to true to suppress screen output:
suppressOutput = false;
%--------------------------------------------------------------------------
M = 1.0; 
dt = 1e+04;
phi_0 = 0.0;
%
% Approximate line search parameter:
stepsize = 1.5;
%
plotFrames = 10;
stepsPerPlot = 10;
maxSteps = stepsPerPlot*plotFrames;
%
stepsPerScreenPlot = 10;
stepsPerReport = 1;
%
h = 0.785;
N = 512;
L = N*h; 
r = -0.4;
%
ratio = L/(2*pi);
xx = zeros(N,N);
yy = zeros(N,N);
%
% Avoid Matlab's meshgrid function:
for j = 1:N
  for i = 1:N
    xx(i,j) = h*i;
    yy(i,j) = h*j;
  end
end
%
[derivx,derivy,lap,lap0] = initOperators(N,ratio,index);
%
% Get initial condition
[phi] = initialCondition(N,xx,yy,L,ratio);
%
% Print initial data:
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
if ~suppressOutput
  figure(1);
  pcolor(xx,yy,phi);
  shading interp;
  title('Initial Conditions');
  colorbar;
  getframe;
end
phio = phi;
%
for k=1:maxSteps
%
  time=k*dt;
%
  phioo = phio;
  phio = phi;
  phi = 2*phio-phioo;
%
  nrmSD = 1.0;
%    
  if (mod(k,stepsPerReport) == 0)
    disp(['step :',num2str(k),'         time :',num2str(time)]);
  end
%
% PSD iteration:
  ell = 0;
  ellMax = 199;
  while nrmSD > 1.0e-09 && ell <= ellMax
    ell = ell+1;
%
    [residual] = GetResidual(phi,phio,lap,lap0,h,dt,r,M,L);
%
    [searchdir] = GetSearchDirection(residual,lap,lap0,h,dt,r,M,L);
%
% Update with a stepsize:
%
    phi = phi+stepsize*searchdir;
%
    nrmSD = max(max(abs(searchdir)));
%       
    if (mod(k,stepsPerReport) == 0 && ~suppressOutput)
      fprintf('    PSD iteration %d, Residual %e\n',ell,nrmSD);
    end
  end
%
  if ell >= ellMax
    disp(['step :',num2str(k),'         time :',num2str(time)]);
    disp('Too many PSD iterations!')
    return
  end
%
  if (mod(k,stepsPerScreenPlot) == 0 && ~suppressOutput)
    figure(1);
    pcolor(xx,yy,phi);
    shading interp; colormap(jet);
    colorbar;
    title(['FCH-BE-PSD: dt = ',num2str(dt),'  time = ',num2str(time)]);
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


