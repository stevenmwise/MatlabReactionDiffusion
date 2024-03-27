function [phi]=initial_condition(N,xx,yy,L,ratio,epsilon,gamma,eta1,eta2, ...
    phi_0,init_type)
%
phi = zeros(N,N);
%
if  init_type == 1
%
    phi=2*exp(sin(xx)+sin(yy)-2)+2.2*exp(-sin(xx)-sin(yy)-2)-1;
%
elseif init_type == 2
%
    phi = phi_0 + 0.001*(2*rand(N,N)-1);
%
elseif init_type == 3
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
    [phi] = smoothed_s_curve(N,L);
%
elseif init_type == 6
%
    [phi] = smoothed_annulus(N,L);
%
elseif init_type == 7 || init_type == 8
%
    [phi] = FCH_Init_SixFCirc(gamma,epsilon,eta1,eta2,N,L);
%
else
%
    disp(["Error: No such initial conditon"])
%
end