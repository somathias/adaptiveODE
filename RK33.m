function [xd,yd,zd,fev]=RK33(i,xc,yc,zc,k1,k2,dt,N,fev0,s,a,ra,mu)
% explicit Runge-Kutta method, order 3, one step in 3D, see (30)
% input
% i: time point to advance to
% (xc(j,i-1),xc(j,i-1),xc(j,i-1)) j=1:N, initial coordinates
% k1: force vector at t(i-1)
% k2: force vector at intermediate time
% dt: time step
% N: number of cells
% fev0: accumulated F evaluations
% s,a,ra,mu: force parameters
% output
% (xd(j),yd(j),zd(j)) j=1:N, end coordinates
% fev: updated F evaluations
for j=1:N
   x(1,j)=xc(j,i-1)+0.75*dt*k2(1,j);
   x(2,j)=yc(j,i-1)+0.75*dt*k2(2,j);
   x(3,j)=zc(j,i-1)+0.75*dt*k2(3,j);
end
% compute the force at intermediate time point
k3=fcmp3(N,x,s,a,ra,mu);
fev=fev0+1;
for j=1:N
   xd(j)=xc(j,i-1)+dt*((2/9)*k1(1,j)+(1/3)*k2(1,j)+(4/9)*k3(1,j));
   yd(j)=yc(j,i-1)+dt*((2/9)*k1(2,j)+(1/3)*k2(2,j)+(4/9)*k3(2,j));
   zd(j)=zc(j,i-1)+dt*((2/9)*k1(3,j)+(1/3)*k2(3,j)+(4/9)*k3(3,j));
end
