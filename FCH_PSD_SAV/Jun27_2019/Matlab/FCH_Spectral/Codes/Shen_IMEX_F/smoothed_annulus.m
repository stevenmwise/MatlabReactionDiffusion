function [phi] = smoothed_s_curve(N,L)
%
M = 2^10;
phi_prime = zeros(M,M);
h_prime = L/M;
%
% The use of meshgrid should be avoided:
%[x_prime,y_prime] = meshgrid((1:M)*h_prime);
%
% Instead, use
x_prime = zeros(M,M);
y_prime = zeros(M,M);
for j = 1:M
    for i = 1:M
        x_prime(i,j) = h_prime*i;
        y_prime(i,j) = h_prime*j;
    end
end
%
alpha_filt = 50*log(10);
%
alpha_alias = [0:M/2-1,M/2,-M/2+1:-1];
%    
filt_x = zeros(M,M);
filt_y = zeros(M,M);
for i=1:M
    filt_x(i,:) = alpha_alias;
    filt_y(:,i) = alpha_alias';
end
%
ita1 = 2*filt_x/M;
ita2 = 2*filt_y/M;
%
filt_x = exp(-alpha_filt*ita1.^2);
filt_y = exp(-alpha_filt*ita2.^2);
%
stretch = 2.0;
%
for i=1:M
    for j=1:M
        r = sqrt((x_prime(i,j)-L/2)^2+((y_prime(i,j)-L/2)^2)/stretch);
        if     (r>L/4+0.2)
            phi_prime(i,j) = -1.0;
        elseif (r<L/4-0.2)
            phi_prime(i,j) = -1.0;
        else
            phi_prime(i,j) =  1.0;
        end
    end
end
%
%figure(1);
%pcolor(x_prime,y_prime,phi_prime);
%shading interp;
%title('Unfiltered Function');
%colorbar;
%
phi_prime_hat = fft2(phi_prime);
%phi_prime_o = phi_prime;
phi_prime = real(ifft2(filt_x.*filt_y.*phi_prime_hat));
%
phi(1:N,1:N) = phi_prime(M/N:M/N:M,M/N:M/N:M);
%
%N1 = 2^9;
%phi1(1:N1,1:N1) = phi_prime(M/N1:M/N1:M,M/N1:M/N1:M);
%N2 = 2^7;
%phi2(1:N2,1:N2) = phi_prime(M/N2:M/N2:M,M/N2:M/N2:M);
%N3 = 2^6;
%phi3(1:N3,1:N3) = phi_prime(M/N3:M/N3:M,M/N3:M/N3:M);
%
%figure(2);
%plot((1:M)*h_prime,phi_prime_o(:,M/2),'-s')
%hold 'on'
%plot([L/N1:L/N1:L],phi1(:,N1/2),'o')
%hold on;
%plot([L/N2:L/N2:L],phi2(:,N2/2),'+')
%hold on;
%plot([L/N3:L/N3:L],phi3(:,N3/2),'x')
%hold on;
%title('Filtered Functions');
%hold off;
%
%figure(3);
%pcolor(x_prime,y_prime,phi_prime);
%shading interp;
%title('Unfiltered Function');
%colorbar;
%
%display('Done')
end