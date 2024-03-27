function [residual] = GetResidual(phi,phio,lap,lap0,h,dt,r,M,L)
%
omega = real(ifft2(fft2(phi).*(1.0+lap).*(1.0+lap)))+phi.*phi.*phi+r*phio;
%
phiDiff = phi-phio;
massPhiDiff = (sum(sum(phiDiff))*h*h)/(L*L);
phiDiff = phiDiff-massPhiDiff;
%
residual = -omega-real(ifft2(fft2(phiDiff)./(-dt*M*lap0)));
%
% Make the residual mean zero:
residual = residual-(sum(sum(residual))*h*h)/(L*L);
end

