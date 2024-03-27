%   This program calculate the FCH equation using the BDF2 scheme with
%       the PSD solver.
%
%   init_type = 1     ::  High accuracy benchmark
%   init_type = 2     ::  random field with average phi_0
%   init_type = 3     ::  S shape
%	init_type = 4     ::  S shape with smoothing filter: N_fine = 2^10
%	init_type = 5     ::  S shape with smoothing filter: N_fine = 2^10
%                           and different parameters.
%	init_type = 6     ::  annulus with smoothing filter: N_fine = 2^10
%	init_type = 7     ::  Keith's six-fold circle: splitting
%	init_type = 8     ::  Keith's six-fold circle
%
%   index == 1        ::  lap = [0:N/2-1 N/2 -N/2+1:-1];
%   index /= 1        ::  lap = [0:N/2-1   0 -N/2+1:-1];
%--------------------------------------------------------------------------
clear;
clc;
index = 1;
init_type = 5;
%--------------------------------------------------------------------------
%
% Set this to true to suppress screen output:
suppressOutput = false;
%--------------------------------------------------------------------------
M = 1; dt = 1e-03;
phi_0 = 0.0;
%
plot_frames = 100;
steps_per_plot = 100;
max_steps = steps_per_plot*plot_frames;
%
steps_per_screen_plot = 10;
steps_per_report = 1;
%
if     init_type==1
    L = 2*pi; N = 2^7;  epsilon = 0.18;
    eta1 = epsilon^2;      eta2 = epsilon^2;    gamma = 0.00;
elseif init_type==2
    L = 2*pi; N = 2^7;  epsilon = 0.04;
    eta1 = 5.0*epsilon;    eta2 = 3.0*epsilon;  gamma = 0.00;
elseif init_type==3 || init_type==4
    L = 12.8; N = 2^8;  epsilon = 0.10;
    eta1 = 0.2;            eta2 = 0.2;          gamma = 0.00;
elseif init_type==5
    L = 12.8; N = 2^8;  epsilon = 0.10;
    eta1 = 1.45*epsilon;   eta2 = 3.0*epsilon;  gamma = 0.30;
elseif init_type==6
    L = 4*pi; N = 2^8;  epsilon = 0.10;
    eta1 = 1.45*epsilon;   eta2 = 2.0*epsilon;  gamma = 0.25;
elseif init_type==7
    L = 4*pi; N = 2^8;  epsilon = 0.10;
    eta1 = 1.45*epsilon;   eta2 = 3.0*epsilon;  gamma = 0.30;
elseif init_type==8
    L = 4*pi; N = 2^8;  epsilon = 0.10;
    eta1 = 1.45*epsilon;   eta2 = 2.0*epsilon;  gamma = 0.30;
else
    disp("No such initial conditon!!!")
    return
end
%
ratio = L/(2*pi);
h = L/N;
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
[derivx,derivy,lap,lap0] = init_operators(N,ratio,index);
%
% Get initial condition
[phi] = initial_condition(N,xx,yy,L,ratio,epsilon,gamma,eta1,eta2,phi_0, ...
    init_type);
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
% Approximate line search parameter:
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
% Update with a stepsize:
%
        phi = phi+stepsize*searchdir;
%
%        nrm_rd = max(max(abs(residual)));
        nrm_sd = max(max(abs(searchdir)));
%       
        if (mod(k,steps_per_report) == 0 && ~suppressOutput)
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
    if (mod(k,steps_per_screen_plot) == 0 && ~suppressOutput)
        figure(1);
        pcolor(xx,yy,phi);
        shading interp; colormap(jet);
        colorbar;
        title(['FCH-BDF2-PSD: dt = ',num2str(dt),'  time = ',num2str(time)]);
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


