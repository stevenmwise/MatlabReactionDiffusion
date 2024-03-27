function [phi_x,phi_y,phi_lap,phi_lap2,phi_lap3,phi3_lap]=cal_phi(phi,derivx,derivy,lap)

phi_hat=fft2(phi);
phi_x=real(ifft2(derivx.*phi_hat));
phi_y=real(ifft2(derivy.*phi_hat));
phi_lap=real(ifft2(lap.*phi_hat));
phi_lap2=real(ifft2(lap.^2.*phi_hat));
phi_lap3=real(ifft2(lap.^3.*phi_hat));
phi3_lap=real(ifft2(lap.*fft2(phi.^3)));
