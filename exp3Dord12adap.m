function [xc,yc,zc,t,nt,fev,dt1]=exp3Dord12adap(ntim,T,mstep,N,x0,y0,z0,dt0,eps,s,a,ra,mu)
% CBM in 3D using forward Euler, multirate time stepping and error estimate by Jacobian
% input
% T: final time
% mstep: number of substeps
% N: number of cells
% (x0(i),y0(i),z0(i), i=1:N) initial cell coordinates
% if errest then local errors are estimated and the time step varied, otherwise the time step is constant dt0
% eps: local error parameter
% s, a, ra, mu: force parameters
% output
% (xc(i,j),yc(i,j),zc(i,j), i=1:N, j=1:nt) cell coordinates in nt time steps
% (t(j), j=1:nt) time points
% fev is the number of function evaluations
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
fev=0;
while ti<T
   fev0=fev;
   % determine the new time step dt1
   [dt1,K0,k0,fev]=newdtadapt3(N,xc(:,i-1),yc(:,i-1),zc(:,i-1),mstep,fev0,eps,s,a,ra,mu);
   % k0 is the force at t(i-1)
   dtsav=dt1;
   fev0=fev;
   ti=t(i-1)+dt1;
   if ti>T
      dt1=T-t(i-1);
      ti=T;
   end
   t(i)=ti;
   % take one step forward from (xc(:,i-1),yc(:,i-1),zc(:,i-1)) to (xc(:,i),yc(:,i),zc(:,i)) with time step dt1
   [xc(:,i),yc(:,i),zc(:,i),fev]=fEadapt3(i,xc,yc,zc,k0,dt1,mstep,fev0,N,K0,s,a,ra,mu);
   i=i+1;
end
dt1=dtsav;
nt=i-1;
