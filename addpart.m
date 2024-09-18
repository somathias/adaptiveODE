function [nc1,xc1,yc1,zc1]=addpart(nc0,xc0,yc0,zc0,r0)
% initiate proliferation of cell ipart at random location
xc1=xc0;
yc1=yc0;
zc1=zc0;
% a random cell
ipart=floor(1+(nc0-0.1)*rand);
dirpart(1)=2*rand-1;
dirpart(2)=2*rand-1;
dirpart(3)=2*rand-1;
dirnorm=norm(dirpart);
dirpart(:)=0.5*r0*dirpart(:)/dirnorm;
xc1(ipart)=xc0(ipart)+dirpart(1);
yc1(ipart)=xc0(ipart)+dirpart(2);
zc1(ipart)=xc0(ipart)+dirpart(3);
nc1=nc0+1;
xc1(nc1)=xc0(ipart)-dirpart(1);
yc1(nc1)=xc0(ipart)-dirpart(2);
zc1(nc1)=xc0(ipart)-dirpart(3);

