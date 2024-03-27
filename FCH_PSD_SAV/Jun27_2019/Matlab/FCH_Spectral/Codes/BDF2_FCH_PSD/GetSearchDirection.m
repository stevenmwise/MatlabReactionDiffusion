function [searchdir] = GetSearchDirection(residual,lap,lap0,h,dt,epsilon, ...
                gamma,eta1,eta2,M,L,k)
%
% Both the residual and search direction should be of mean zero.
%
% The residual must be made mean zero in the function GetResidual.
%
if k == 1
    coef = -1.0./(    M*dt*lap0)+1.0+epsilon^4*lap.^2;
else
    coef = -3.0./(2.0*M*dt*lap0)+1.0+epsilon^4*lap.^2;
end
%
searchdir = real(ifft2(fft2(residual)./coef));
%
% The next step is not really necessary.
searchdir = searchdir-(sum(sum(searchdir))*h*h)/(L*L);
end

