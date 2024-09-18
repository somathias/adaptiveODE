function [xd,yd,zd,k1,iter,fev]=bE3(i,xc,yc,zc,dt,N,fev0,epsN,s,a,ra,mu)
% backward Euler, one step
% input
% i: time point to advance to
% (xc(j,i-1),xc(j,i-1),xc(j,i-1)) j=1:N, initial coordinates
% dt: time step
% N: number of cells
% fev0: accumulated F evaluations
% epsN: tolerance in nonlinear iterations
% s,a,ra,mu: force parameters
% output
% (xd(j),yd(j),zd(j)) j=1:N, end coordinates
% k1: force vector at t(i)
% iter: number of iterations
% fev: updated F evaluations
x(1,:)=xc(:,i-1);
x(2,:)=yc(:,i-1);
x(3,:)=zc(:,i-1);
% solve k1=F(x+dt*k1), x1=x+dt*k1 by nonlinear conjugate gradient method
[x1,k1,iter,fev]=nonlincg3(N,x,x,dt,1,fev0,epsN,s,a,ra,mu);
xd(:)=x1(1,:);
yd(:)=x1(2,:);
zd(:)=x1(3,:);
end
