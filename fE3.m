function [xd,yd,zd,k1,fev]=fE3(i,xc,yc,zc,dt,N,fev0,s,a,ra,mu)
% forward Euler fE, one step in 3D
% input
% i: time point to advance to
% (xc(j,i-1),xc(j,i-1),xc(j,i-1)) j=1:N, initial coordinates
% dt: time step
% N: number of cells
% fev0: accumulated F evaluations
% s,a,ra,mu: force parameters
% output
% (xd(j),yd(j),zd(j)) j=1:N, end coordinates
% k1: force vector at t(i-1)
% fev: updated F evaluations
x(1,:)=xc(:,i-1);
x(2,:)=yc(:,i-1);
x(3,:)=zc(:,i-1);
% compute the force
k1=fcmp3(N,x,s,a,ra,mu);  
for j=1:N
   xd(j)=xc(j,i-1)+dt*k1(1,j);
   yd(j)=yc(j,i-1)+dt*k1(2,j);
   zd(j)=zc(j,i-1)+dt*k1(3,j);
end
fev=fev0+1;
