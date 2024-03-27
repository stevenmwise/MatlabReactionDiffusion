%   This program calculate the FCH equation using the BDF2 scheme by SAV
%    method.
%
%   init_type=1     ::  High accuracy benchmark
%   init_type=2     ::  phi_0 random number
%   init_type=3     ::  S shape
%	init_type=4     ::  S shape  with filter by 2^10
%	init_type=5     ::  circle ring or elliptical ring for pearling
%	init_type=6     ::  circle
%-------------------------------------------------------------------------------
clear ;clc;
index     = 1;
init_type = 6;
%-------------------------------------------------------------------------------
M=1;
gamma = 0.25;
phi_0 = 0.6;
dt = 1e-4;
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
S1=5*epsilon^1;    S2=0;      B=L^2;
%
%-------------------------------------------------------------------------------
%
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
[psi] = initial_condition(N,xx,yy,L,ratio,epsilon,phi_0,init_type);
%
%-------------------------------------------------------------------------------
% Get initial condition

s1 = ['0000000' num2str(0)];
s2 = s1((length(s1)-4):length(s1));
fid = fopen(['./OUT/phi',s2,'.dat'],'w');
fprintf(fid,'%25.15e %25.15e %10i %25.15e\n',0,dt,N,L);
for j = 1:N
    for i = 1:N
        fprintf(fid,'%25.15e\n',psi(i,j));
    end
end
fclose(fid);
%
figure(1);
pcolor(xx,yy,psi);
shading interp;
title('Initial Conditions');
colorbar;
getframe;
%
%-------------------------------------------------------------------------------
%
phi = psi+gamma/6.0;
%
for k=1:max_steps
%
    time = k*dt;
%
    if k==1
%
        phin = phi;
%
        [phin_x,phin_y,phin_lap,phin_lap2,phin_lap3,phin3_lap] ...
            = cal_phi(phin,derivx,derivy,lap);
%
        [U0,Zn] = cal_UZ(epsilon,eta2,gamma,h,B,psi,phin,phin_x,phin_y, ...
            phin_lap,phin3_lap);
%
        Un = U0;
%
        Zn_hat = fft2(Zn);
        Zn_lap = real(ifft2(lap.*Zn_hat));
%
        [g] = cal_g1(epsilon,eta1,gamma,M,dt,h,S1,S2,Un,Zn,Zn_lap,phin,phin_lap,phin_lap2);
%
        coef = 1.0/(M*dt)-epsilon^4*lap.^3+S1*lap.^2-S2*lap;
% 
        X1 = real(ifft2((lap.*Zn_hat)./coef));
        X2 = real(ifft2(fft2(g)./coef));
%
        Y1 = sum(sum(Zn.*X1))*h*h;
        Y2 = sum(sum(Zn.*X2))*h*h;
        Y3 = Y2/(1-0.5*Y1);
%
        phi = X2+0.5*Y3*X1;
%
        Un1 = Un;
        Un = Un+0.5*sum(sum(Zn.*(phi-phin)))*h*h;
%
    else
%
        Zstar = 2*Zn-Zn1;
        Zstar_hat = fft2(Zstar);
        Zstar_lap = real(ifft2(lap.*Zstar_hat));
%
        phicross = 4*phin-phin1;
        Ucross = 4*Un-Un1;
%
        [g] = cal_g2(M,dt,epsilon,eta1,gamma,h,S1,S2,Ucross,Zstar, ...
            Zstar_lap,phicross, ...
            phistar_lap,phistar_lap2);
%
        coef = 3.0/(2.0*M*dt)-epsilon^4*lap.^3+S1*lap.^2-S2*lap;
%
        X1 = real(ifft2((lap.*Zstar_hat)./coef));
        X2 = real(ifft2(fft2(g)./coef));
%
        Y1 = sum(sum(Zstar.*X1))*h*h;
        Y2 = sum(sum(Zstar.*X2))*h*h;
        Y3 = Y2/(1-0.5*Y1);
%
        phi = X2+0.5*Y3*X1;
%
        Un1 = Un;
        Un = Ucross/3+0.5*sum(sum(Zstar.*(3*phi-4*phin+phin1)/3))*h*h;
%     
    end
%  
    psi = phi-gamma/6.0;
%  
    phin1 = phin;
    phin   = phi;
    phin1_x = phin_x;
    phin1_y = phin_y;
    phistar = 2.0*phin-phin1;
%
    [phin_x,phin_y,phin_lap,phin_lap2,phin_lap3,phin3_lap] ...
        = cal_phi(phin,derivx,derivy,lap);
%
    Zn1 = Zn;
%
    [U0,Zn] = cal_UZ(epsilon,eta2,gamma,h,B,psi,phin,phin_x,phin_y,phin_lap,phin3_lap);
%
    [phistar_x,phistar_y,phistar_lap,phistar_lap2,phistar_lap3,phistar3_lap] ...
        = cal_phi(phistar,derivx,derivy,lap);
%
    if (mod(k,steps_per_report) == 0)
        disp(['step :',num2str(k),'         time :',num2str(time)]);
    end
        
    if (mod(k,steps_per_screen_plot) == 0)
        figure(1);
        pcolor(xx,yy,psi);
        shading interp; colormap(jet);
        colorbar;
        title(['FCH-SAV-BDF2-',num2str(dt),'-time :',num2str(time)]);
        getframe;
    end
%
    if (mod(k,steps_per_plot) == 0)
        s1 = ['0000000' num2str(round(k/steps_per_plot))];
        s2 = s1((length(s1)-4):length(s1));
        fid = fopen(['./OUT/phi',s2,'.dat'],'w');
        fprintf(fid,'%25.15e %25.15e %10i %25.15e\n',time,dt,N,L);
        for j = 1:N
            for i = 1:N
                fprintf(fid,'%25.15e\n',psi(i,j));
            end
        end
        fclose(fid);
    end
    
end
%
disp('program done')


