function [xd,yd,zd,k1,k2,iter,fev]=TM3(i,xc,yc,zc,dt,N,fev0,epsN,s,a,ra,mu)
% trapezoidal method, order 2, one step in 3D
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
% k2: force vector at intermediate time point
% iter: number of iterations
% fev: updated F evaluations
x(1,:)=xc(:,i-1);
x(2,:)=yc(:,i-1);
x(3,:)=zc(:,i-1);
% force vector at t(i-1)
k1=fcmp3(N,x,s,a,ra,mu);
fev1=fev0+1;
% solve k2=F(x+0.5*dt*(k1+k2)), x1=x+0.5*dt*(k1+k2) by nonlinear conjugate gradient method
f0=x+0.5*dt*k1;
[x1,k2,iter,fev]=nonlincg3(N,x,f0,dt,0.5,fev1,epsN,s,a,ra,mu);
xd(:)=x1(1,:);
yd(:)=x1(2,:);
zd(:)=x1(3,:);
end
