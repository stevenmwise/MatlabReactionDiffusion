
function [ ] = movieMaker(videoName,frameRate,plotFrames,cmin,cmax)
%
% This program makes image files and a movie from the saved data 
% files in the ./OUT directory. The program assumes that the files
% in that directory are named phi00ddd.dat, where ddd is a non-
% negative integer, and have a certain format.
%
% The typical calling structure is
%
%     movieMaker('CahnHilliard',10,20,0,1)
%
% videoName is a string or character vector
%
% frameRate is the number of frames per second
%
% plotFrames is the total number of frames (not counting the initial
%   condition) to be plotted
%
% cmin is the desired minimum value of the colorbar
%
% cmax is the desired maximum value of the colorbar
%
s0 = [videoName, '.mp4'];
vidObj = VideoWriter(s0,'MPEG-4');
vidObj.FrameRate = frameRate;
open(vidObj);
%
for k = 0:plotFrames
  display(k)
  plotFrame(k,videoName,cmin,cmax)
  set(gca,'nextplot','replacechildren');
  currFrame = getframe(gcf);
  writeVideo(vidObj,currFrame);
end
%
close(vidObj);
%
end % function movieMaker
%
% Embedded function(s) below:
%
function [ ] = plotFrame(k,plotName,cmin,cmax)
%
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
caxis([cmin,cmax])
colormap(jet);
caxis([cmin,cmax])
colorbar;
title([plotName, '   dt = ', num2str(dt), '   time = ', ...
  num2str(time)]);
%
s3 = ['./OUT/phi',s2,'.jpg']
print('-djpeg','-r400',s3)
end % function PlotFrame