function [xd,yd,zd,k1,k2,fev]=RK23(i,xc,yc,zc,k1cmp,k1,dt,N,fev0,s,a,ra,mu)
% Runge's method RK2, order 2, one step in 3D
% input
% i: time point to advance to
% (xc(j,i-1),xc(j,i-1),xc(j,i-1)) j=1:N, initial coordinates
% if k1cmp then compute new k1
% k1: force vector at t(i-1)
% dt: time step
% N: number of cells
% fev0: accumulated F evaluations
% s,a,ra,mu: force parameters
% output
% (xd(j),yd(j),zd(j)) j=1:N, end coordinates
% k1: force vector at t(i-1)
% fev: updated F evaluations
fev1=fev0;
if k1cmp
  % compute forces at t(i-1)
  x(1,:)=xc(:,i-1);
  x(2,:)=yc(:,i-1);
  x(3,:)=zc(:,i-1);
  k1=fcmp3(N,x,s,a,ra,mu);
  fev1=fev1+1;
end
% intermediate coordinates, see (24)
for j=1:N
   x(1,j)=xc(j,i-1)+0.5*dt*k1(1,j);
   x(2,j)=yc(j,i-1)+0.5*dt*k1(2,j);
   x(3,j)=zc(j,i-1)+0.5*dt*k1(3,j);
end
% corresponding force vector 
k2=fcmp3(N,x,s,a,ra,mu);
fev1=fev1+1;
% solution at t(i)
for j=1:N
   xd(j)=xc(j,i-1)+dt*k2(1,j);
   yd(j)=yc(j,i-1)+dt*k2(2,j);
   zd(j)=zc(j,i-1)+dt*k2(3,j);
end
fev=fev1;
