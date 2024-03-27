%   This program calculate the FCH equation using the Backward Euler scheme with
%       the PSD solver.
%
%   init_type = 1     ::  High accuracy benchmark
%   init_type = 2     ::  phi_0 random number
%   init_type = 3     ::  S shape
%	init_type = 4     ::  S shape  with filter by 2^10
%	init_type = 5     ::  circle ring or elliptical ring for pearling
%	init_type = 6     ::  circle
%
%   index = 1         ::  lap = [0:N/2-1 N/2 -N/2+1:-1];
%   index /= 1        ::  lap = [0:N/2-1   0 -N/2+1:-1];
%
%   step_plot         ::  plot the picture every 'step_plot' steps.
%   Energy_max        ::  Energy=zeros(1,Energy_max)
%   ie                ::  the index of  'Energy'.
%--------------------------------------------------------------------------
clear;
clc;
index = 1;
init_type = 4;
%--------------------------------------------------------------------------
M = 1;
phi_0 = 0.0;
dt = 1e-3;
steps_per_plot = 10000;
steps_per_screen_plot = 1000;
steps_per_report = 100;
plot_frames = 100;
max_steps = steps_per_plot*plot_frames;
%
if     init_type==1
    L=2*pi; N=2^7;  epsilon=0.18;
    eta1=epsilon^2;      eta2=epsilon^2;    gamma = 0.00;
elseif init_type==2
    L=2*pi; N=2^7;  epsilon=0.04;
    eta1=5.0*epsilon;    eta2=3.0*epsilon;  gamma = 0.00;
elseif init_type==3 || init_type==4
    L=12.8; N=2^8;  epsilon=0.1;
    eta1=0.2;            eta2=0.2;          gamma = 0.00;
elseif init_type==5
    L=2*pi; N=2^8;  epsilon=0.04;
    eta1=2.0*epsilon;    eta2=2.0*epsilon;  gamma = 0.25;
elseif init_type==6
    L=4*pi; N=2^8;  epsilon=0.1;
    eta1=1.45*epsilon;   eta2=2.0*epsilon;  gamma = 0.25;
else
    disp("No such initial conditon!!!")
    return
end
%
% Grid and  derivative matrix
ratio = L/(2*pi);
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
[derivx,derivy,lap,lap0] = init_operators(N,ratio,index);
%
% Get initial condition
[phi] = initial_condition(N,xx,yy,L,ratio,epsilon,phi_0,init_type);
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
figure(1);
pcolor(xx,yy,phi);
shading interp;
title('Initial Conditions');
colorbar;
getframe;
phio = phi;
%
stepsize = 0.72;
%
for k=1:max_steps
%
    time=k*dt;
%
    phioo = phio;
    phio = phi;
    phi = 2*phio-phioo;
%
    nrm_sd = 1.0;
%    
    if (mod(k,steps_per_report) == 0)
        disp(['step :',num2str(k),'         time :',num2str(time)]);
    end
%
% PSD iteration:
    ell = 0;
    ell_max = 199;
    while nrm_sd > 1.0e-09 && ell <= ell_max
        ell = ell+1;
%
        [residual] = GetResidual(phi,phio,phioo,lap,lap0,h,dt,epsilon,gamma, ...
                eta1,eta2,M,L,k);
%
        [searchdir] = GetSearchDirection(residual,lap,lap0,h,dt,epsilon, ...
                gamma,eta1,eta2,M,L,k);
%
% Update with a stepsize of stepp = 1/Lip:
%
        phi = phi+stepsize*searchdir;
%
%        nrm_rd = max(max(abs(residual)));
        nrm_sd = max(max(abs(searchdir)));
%       
        if (mod(k,steps_per_report) == 0)
            fprintf('    PSD iteration %d, Residual %e\n',ell,nrm_sd);
        end
    end
%
    if ell >= ell_max
        disp(['step :',num2str(k),'         time :',num2str(time)]);
        display('Too many PSD iterations!')
        return
    end
%
    if (mod(k,steps_per_screen_plot) == 0)
        figure(1);
        pcolor(xx,yy,phi);
        shading interp; colormap(jet);
        colorbar;
        title(['FCH-BDF2: dt = ',num2str(dt),'  time = ',num2str(time)]);
        getframe;
    end
%
    if (mod(k,steps_per_plot) == 0)
%
        s1 = ['0000000' num2str(round(k/steps_per_plot))];
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


