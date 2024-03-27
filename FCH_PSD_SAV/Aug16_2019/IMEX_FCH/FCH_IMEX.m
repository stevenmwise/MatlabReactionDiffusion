%   This program calculate the FCH equation using a linear IMEX scheme.
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
init_type = 7;
%--------------------------------------------------------------------------
%
% Set this to true to suppress screen output:
suppressOutput = false;
%--------------------------------------------------------------------------
M = 1.0; dt = 1e-03;
phi_0 = 0.0;
%
plot_frames = 100;
steps_per_plot = 100;
max_steps = steps_per_plot*plot_frames;
%
steps_per_screen_plot = 10;
steps_per_report = 1;
%
sig1 = 1.0; sig2 = 1.0;
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
%
for k=1:max_steps
%
    time = k*dt;
%
    if (mod(k,steps_per_report) == 0 && ~suppressOutput)
        disp(['step :',num2str(k),'         time :',num2str(time)]);
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
    coef = 1.0-dt*M*epsilon^4*lap.^3+dt*M*epsilon^2*(sig2+sig1-eta1)*lap.^2 ...
            -dt*M*(sig2*(sig1-eta1))*lap;
%
    phi = real(ifft2(fft2(q4)./coef));
%
    if (mod(k,steps_per_screen_plot) == 0 && ~suppressOutput)
        figure(1);
        pcolor(xx,yy,phi);
        shading interp; colormap(jet);
        colorbar;
        title(['FCH-IMEX: dt = ',num2str(dt),'  time = ',num2str(time)]);
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


