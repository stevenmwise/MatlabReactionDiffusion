function [searchdir] = GetSearchDirection(residual,lap,lap0,h,dt,r,M,L)
%
% Both the residual and search direction should be of mean zero.
%
coef = -1.0./(M*dt*lap0)+1.0+(1.0+lap).*(1.0+lap);
%
searchdir = real(ifft2(fft2(residual)./coef));
%
searchdir = searchdir-(sum(sum(searchdir))*h*h)/(L*L);
end

