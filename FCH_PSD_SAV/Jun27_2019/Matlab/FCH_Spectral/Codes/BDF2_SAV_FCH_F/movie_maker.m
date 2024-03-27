
function [ ] = movie_maker(video_name,frame_rate,plot_frames)
%
% This program makes image files and a movie from the saved data files.
%
s0 = [video_name, '.mp4'];
vidObj = VideoWriter(s0,'MPEG-4');
vidObj.FrameRate = frame_rate;
open(vidObj);
%
for k = 0:plot_frames
    display(k)
    plot_frame(k,video_name)
    set(gca,'nextplot','replacechildren');
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end
%
close(vidObj);
%
end

function [ ] = plot_frame(k,plot_name)

s1 = ['000000' num2str(k)];
s2 = s1((length(s1)-4):length(s1));
fid = fopen(['./OUT/phi',s2,'.dat'],'r');
%
[time] =  fscanf(fid, '%f', 1);
[dt  ] =  fscanf(fid, '%f', 1);
[N   ] =  fscanf(fid, '%i', 1);
[L   ] =  fscanf(fid, '%f', 1);
%
h = L/N;
%
xx  = zeros(N,N,'double');
yy  = zeros(N,N,'double');
phi = zeros(N,N,'double');
%
A = zeros(1,N*N,'double');
[A]=fscanf(fid,'%f',[1,N*N]);
%
fclose(fid);
%
for j = 1:N
  for i = 1:N
    xx(i,j) = h*i;
    yy(i,j) = h*j;
    phi(i,j) = A(1,(j-1)*N+i);
  end
end
%
pcolor(xx,yy,phi);
axis([0,L,0,L])
axis equal
axis([0,L,0,L])
shading interp;
caxis([-1 1])
colormap(jet);
caxis([-1 1])
colorbar;
title([plot_name, '   dt = ', num2str(dt), '   time = ', num2str(time)]);
%
s3 = ['./OUT/phi',s2,'.jpg']
print('-djpeg','-r400',s3)
end