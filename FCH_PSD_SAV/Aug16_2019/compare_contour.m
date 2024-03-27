function [l_inf_diff,l_2_diff] = compare_contour(k,plot_name)
%
s1 = ['000000' num2str(k)];
s2 = s1((length(s1)-4):length(s1));
fid = fopen(['BDF2_FCH_PSD/OUT/phi',s2,'.dat'],'r');
%
[time] =  fscanf(fid, '%f', 1);
[dt  ] =  fscanf(fid, '%f', 1);
[N   ] =  fscanf(fid, '%i', 1);
[L   ] =  fscanf(fid, '%f', 1);
%
A1 = zeros(1,N*N,'double');
[A1]=fscanf(fid,'%f',[1,N*N]);
%
fclose(fid);
%
s1 = ['000000' num2str(k)];
s2 = s1((length(s1)-4):length(s1));
fid = fopen([plot_name,'/OUT/phi',s2,'.dat'],'r');
%
[time] =  fscanf(fid, '%f', 1);
[dt  ] =  fscanf(fid, '%f', 1);
[N   ] =  fscanf(fid, '%i', 1);
[L   ] =  fscanf(fid, '%f', 1);
%
A2 = zeros(1,N*N,'double');
[A2]=fscanf(fid,'%f',[1,N*N]);
%
fclose(fid);
%
h = L/N;
xx  = zeros(N,N,'double');
yy  = zeros(N,N,'double');
phi1 = zeros(N,N,'double');
phi2 = zeros(N,N,'double');
%
for j = 1:N
  for i = 1:N
    xx(i,j) = h*i;
    yy(i,j) = h*j;
    phi1(i,j) = A1(1,(j-1)*N+i);
    phi2(i,j) = A2(1,(j-1)*N+i);
  end
end
%
l_inf_diff = max(max(abs(phi1-phi2)));
l_2_diff = sqrt((h/L)*(h/L)*sum(sum((phi1-phi2).*(phi1-phi2))));
%
clf
contour(xx,yy,phi1,[-0.2 0.2],'Color','black','LineWidth',1.5,'LineStyle','-.')
hold on
[C,h] = contour(xx,yy,phi2,[-0.2 0.2],'Color','red','LineWidth',0.5);
%clabel(C,h,'FontWeight','bold','Color','green');
axis([L/2,0.8*L,L/2,0.8*L])
axis equal
axis([L/2,0.8*L,L/2,0.8*L])
plot_title = strrep(plot_name, '_', ' ');
title([plot_title, '  Comparison @ time = ', num2str(time)]);
legend('BDF2 FCH PSD F',plot_title)
err = sprintf('%2.3e',l_inf_diff);
text(7,7.25,['L_{\infty} diff =  ',err]);
err = sprintf('%6.3e',l_2_diff);
text(7,7.00,['L_2 diff =  ',err]);
%
s3 = ['ComparisonPlots/',plot_name,s2,'.jpg'];
print('-djpeg','-r400',s3)
end