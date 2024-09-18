function [xc,yc,zc,t,nt,fev,dt1]=exp3Dord23(ntim,T,N,x0,y0,z0,dt0,errest,intst,eps,s,a,ra,mu)
% CBM in 3D using RK2 and error estimate by RK3
% input
% T: final time
% N: number of cells
% (x0(i),y0(i),z0(i), i=1:N) initial cell coordinates
% if errest then local errors are estimated and the time step varied, otherwise the time step is constant dt0
% if intst then the integration is started by finding the first time step
% eps: local error parameter
% s, a, ra, mu: force parameters
% output
% (xc(i,j),yc(i,j),zc(i,j), i=1:N, j=1:nt) cell coordinates in nt time steps
% (t(j), j=1:nt) time points
% fev is number of function evaluations
ntim1=ntim+1;
xc=zeros(N,ntim1);
yc=zeros(N,ntim1);
zc=zeros(N,ntim1);
t=zeros(ntim1,1);
xc(:,1)=x0;
yc(:,1)=y0;
zc(:,1)=z0;
i=2;
ti=0;
fev0=0;
if errest && intst
   % find first step, dt0 initial guess, integrate by RK2
   [xc(:,i),yc(:,i),zc(:,i),k1,k2,fev1]=RK23(i,xc,yc,zc,1,x0,dt0,N,fev0,s,a,ra,mu);
   fev0=fev1;
   % integrate by RK3
   [xd,yd,zd,fev1]=RK33(i,xc,yc,zc,k1,k2,dt0,N,fev0,s,a,ra,mu);
   % determine the first time step
   dt=newdt3(xc(:,i),yc(:,i),zc(:,i),xd',yd',zd',dt0,eps,2);
else
   dt=dt0;
   k1=zeros(3,ntim1);
   fev1=0;
end
while ti<T
   % integrate until final time T
   ti=t(i-1)+dt;
   if ti>T
      dt=T-t(i-1);
      ti=T;
   else
      dtsav=dt;
   end
   t(i)=ti;
   fev0=fev1;
   % take one step forward with RK2 from (xc(:,i-1),yc(:,i-1),zc(:,i-1)) to (xc(:,i),yc(:,i),zc(:,i)) with time step dt
   [xc(:,i),yc(:,i),zc(:,i),k1,k2,fev1]=RK23(i,xc,yc,zc,1,k1,dt,N,fev0,s,a,ra,mu);
   if errest
       fev0=fev1;
      % take the step with RK3
      [xd,yd,zd,fev1]=RK33(i,xc,yc,zc,k1,k2,dt,N,fev0,s,a,ra,mu);
      % determine a new time step
      dt1=newdt3(xc(:,i),yc(:,i),zc(:,i),xd',yd',zd',dt,eps,2);
      dt=dt1;
   end
   i=i+1;
end
dt1=dtsav;
nt=i-1;
fev=fev1;

