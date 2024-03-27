function [phi]=initial_condition(N,xx,yy,L,ratio,epsilon,phi_0,init_type)
%
phi = zeros(N,N);
flag = 1;
%
if  init_type==1
%
    phi=2*exp(sin(xx)+sin(yy)-2)+2.2*exp(-sin(xx)-sin(yy)-2)-1;
%
elseif init_type==2
%
    %rand('state',0);
    phi = phi_0 + 0.001*(2*rand(N,N)-1);
%
elseif init_type==3
    r = 0.34;
%
    for i=1:N
        for j=1:N
            if (xx(i,j)>sin(4*pi*yy(i,j)/L)+L/2+r)
                phi(i,j)=-1.0;
            elseif (xx(i,j)<sin(4*pi*yy(i,j)/L)+L/2-r)
                phi(i,j)=-1.0;
            else
                phi(i,j)=1.0;
            end
        end
    end
%
elseif init_type==4
%
    [phi] = smoothed_s_curve(N,L);
%
elseif init_type==5
%
    phi1= zeros(N,N);           
    phi2= zeros(N,N);
%  
    if flag == 1      %circle ring
        a = 1.0;
        b = 1.0;
    elseif flag == 2  %Elliptical ring
        a=1.0;
        b=2.0;
    else
        disp("No such circle for pearling!!!")
    end
%
    r1=L/3;
    r2=r1*0.95;
%
    for i=1:N
        for j=1:N
            dist1=sqrt((xx(i,j)-L/2)^2/a^2+(yy(i,j)-L/2)^2/b^2);
            phi1(i,j)=tanh((abs(dist1)-r1)/(sqrt(2.0)*epsilon));
            phi2(i,j)=tanh((abs(dist1)-r2)/(sqrt(2.0)*epsilon));
            phi(i,j)=-phi1(i,j)*phi2(i,j);
        end
    end
%
    phi=2*( phi-min(min(phi)))/( max(max(phi))-min(min(phi)))-1;
%
 elseif init_type==6
%
    [phi] = smoothed_annulus(N,L);
%
else
%
    disp(["No such initial conditon"])
%
end