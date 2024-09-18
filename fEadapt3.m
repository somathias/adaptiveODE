function [xd,yd,zd,fev]=fEadapt3(i,xc,yc,zc,k0,dt1,mstep,fev0,N,K0,s,a,ra,mu)
% forward Euler, one step with multirate method in 3D
% input
% i: time point to advance to
% (xc(j,i-1),xc(j,i-1),xc(j,i-1)) j=1:N, initial coordinates
% k0: force vector evaluated at t(i-1)
% dt1: time step
% mstep: number of substeps
% fev0: accumulated F evaluations
% K0: set of fast particles
% N: number of cells
% s,a,ra,mu: force parameters
% output
% (xd(j),yd(j),zd(j)) j=1:N, end coordinates
% fev: updated F evaluations
k1=zeros(3,N);
dt0=dt1/mstep;
x(1,:)=xc(:,i-1);
x(2,:)=yc(:,i-1);
x(3,:)=zc(:,i-1);
% compute forces for slow particles in K1 at t(i-1)
K1=ones(N,1)-K0;
fev1=0;
sumK0=sum(K0);
if sumK0>0
   % there are fast particles
   % step one long step forward with the fast particles
   for k=1:mstep
      % take mstep short steps
      % compute forces for fast particles in K0
      if k>1
         % compute forces selected by K0
         k1=fcmpadapt3(N,K0,x,s,a,ra,mu);
      else
	 % use forces in k0 in the first substep
         for j=1:N
	    if K0(j)
	       k1(:,j)=k0(:,j);
            end
         end
      end
      x1=x+dt0*k1;
      x=x1;
   end
   % update F evaluations, only part of the force vector has been changed
   fev1=fev1+(mstep-1)*sumK0/N;
end
% step one long step with the slow particles 
for j=1:N
  if K1(j)
      % integrate slow particles
      xd(j)=xc(j,i-1)+dt1*k0(1,j);
      yd(j)=yc(j,i-1)+dt1*k0(2,j);
      zd(j)=zc(j,i-1)+dt1*k0(3,j);
   else
      % fast particles have already been integrated
      xd(j)=x(1,j);
      yd(j)=x(2,j);
      zd(j)=x(3,j);
   end
end
fev=fev0+fev1;
