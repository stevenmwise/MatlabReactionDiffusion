function W=W(u)
global bm bp gam 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%555
%   Keith's Well
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W=(u-bm).^2/2.*((u-bp).^2/2+gam.*(u-(3*bp-bm)/2)/3);
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%555
% %   Arjen's Well
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p = 3;
% Tildemp = (p/2)^(1/(p-2));
% TildeWp =@(x) 1/(p-2)*(p.*x.^2-2.*x.^p);
% W = TildeWp(u+1)+20.*(u-Tildemp+1).^(p+1).*heaviside(u-Tildemp+1);
end
