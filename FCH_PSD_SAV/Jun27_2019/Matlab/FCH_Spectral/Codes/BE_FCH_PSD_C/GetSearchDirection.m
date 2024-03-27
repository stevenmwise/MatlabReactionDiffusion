function [searchdir] = GetSearchDirection(residual,lap,lap0,h,dt,epsilon, ...
                gamma,eta1,eta2,M,L)
%
% Both the residual and search direction should be of mean zero.
%
coef = -1.0./(M*dt*lap0)+1.0+epsilon^4*lap.^2;
%
searchdir = real(ifft2(fft2(residual)./coef));
%
searchdir = searchdir-(sum(sum(searchdir))*h*h)/(L*L);
end

