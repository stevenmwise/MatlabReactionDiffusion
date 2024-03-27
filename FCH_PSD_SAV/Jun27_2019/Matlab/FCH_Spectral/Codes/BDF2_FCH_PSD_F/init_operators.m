function [derivx,derivy,lap,lap0] = init_operators(N,ratio,index)
%
if (mod(N,2)==0)
    alpha_alias  = [0:N/2-1 N/2 -N/2+1:-1];
    alpha_alias2 = [0:N/2-1   0 -N/2+1:-1];
else
    alpha_alias  = [0:(N-1)/2 -(N-1)/2:-1];
    alpha_alias2 = [0:(N-1)/2 -(N-1)/2:-1];
end
%
alpha  = zeros(N,N); beta  = zeros(N,N);
alpha2 = zeros(N,N); beta2 = zeros(N,N);
%
for i=1:N
    alpha( i,:) = alpha_alias  ;
    beta(  :,i) = alpha_alias' ;
    alpha2(i,:) = alpha_alias2 ;
    beta2( :,i) = alpha_alias2';
end
%
derivx = sqrt(-1)*alpha2/ratio;
derivy = sqrt(-1)*beta2/ratio;
%
if index == 1
    lap = - (alpha .^2+beta .^2)/(ratio^2);
else
    lap = - (alpha2.^2+beta2.^2)/(ratio^2);
end
%
lap0 = lap; lap0(1,1) = 1.0;