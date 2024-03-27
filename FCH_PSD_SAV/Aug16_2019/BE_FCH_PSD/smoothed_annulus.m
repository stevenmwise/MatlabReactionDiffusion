function [phi] = smoothed_s_curve(N,L)
%
M = 2^10;
phi_prime = zeros(M,M);
h_prime = L/M;
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
phi_prime_hat = fft2(phi_prime);
%
phi_prime = real(ifft2(filt_x.*filt_y.*phi_prime_hat));
%
phi(1:N,1:N) = phi_prime(M/N:M/N:M,M/N:M/N:M);
%
end