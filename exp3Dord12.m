function [xc,yc,zc,t,nt,fev,dt1]=exp3Dord12(ntim,T,N,x0,y0,z0,dt0,errest,intst,eps,s,a,ra,mu)
% CBM in 3D using forward Euler (fE) and error estimate by RK2
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
fev=0;
i=2;
ti=0;
if errest && intst
   fev0=fev;
   % find first step, dt0 initial guess, integrate by fE
   [xc(:,i),yc(:,i),zc(:,i),k1,fev]=fE3(i,xc,yc,zc,dt0,N,fev0,s,a,ra,mu);
   fev0=fev;
   % integrate by RK2
   [xd,yd,zd,k1,k2,fev]=RK23(i,xc,yc,zc,0,k1,dt0,N,fev0,s,a,ra,mu);
   % determine the first time step
   dt=newdt3(xc(:,i),yc(:,i),zc(:,i),xd',yd',zd',dt0,eps,1);
else
   dt=dt0;
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
   fev0=fev;
   % take one step forward with fE from (xc(:,i-1),yc(:,i-1),zc(:,i-1)) to (xc(:,i),yc(:,i),zc(:,i)) with time step dt
   [xc(:,i),yc(:,i),zc(:,i),k1,fev]=fE3(i,xc,yc,zc,dt,N,fev0,s,a,ra,mu);
   % k1 evaluated at t(i-1)
   if errest && (ti<T)
      fev0=fev;
      % take the step with RK2
      [xd,yd,zd,k1,k2,fev]=RK23(i,xc,yc,zc,0,k1,dt,N,fev0,s,a,ra,mu);
      % determine the new time step
      dt=newdt3(xc(:,i),yc(:,i),zc(:,i),xd',yd',zd',dt,eps,1);
   end
   i=i+1;
end
dt1=dtsav;
nt=i-1;
 
