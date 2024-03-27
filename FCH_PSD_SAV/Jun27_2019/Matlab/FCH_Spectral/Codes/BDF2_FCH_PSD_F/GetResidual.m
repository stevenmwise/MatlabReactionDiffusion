function [residual] = GetResidual(phi,phio,phioo,lap,lap0,h,dt,epsilon,gamma,eta1, ...
                eta2,M,L,k)
%
fp  = (phi-1.0).*(phi+1.0).*(phi+0.5*gamma);
fpp = (phi+1.0).*(phi+0.5*gamma)+(phi-1.0).*(phi+0.5*gamma) ...
            +(phi-1.0).*(phi+1.0);
%
omega = epsilon^2*real(ifft2(fft2(phi).*lap))-fp;
%
if k == 1
    phi_diff = phi-phio;
else
    phi_diff = (3.0*phi-4.0*phio+phioo)/2;
end
%
mass_phi_diff = (sum(sum(phi_diff))*h*h)/(L*L);
phi_diff = phi_diff-mass_phi_diff;
%
residual = -(epsilon^2*real(ifft2(fft2(omega).*lap))-omega.*fpp ...
            +eta1*omega+(eta1-eta2)*fp ...
            +real(ifft2(fft2(phi_diff)./(-dt*M*lap0))));
%
% Make the residual mean zero. This is absolutely necessary:
residual = residual-(sum(sum(residual))*h*h)/(L*L);
end

