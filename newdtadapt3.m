function [dt1,K0,k1,fev]=newdtadapt3(N,x0,y0,z0,mstep,fev0,eps,s,a,ra,mu)
% compute new time step dt1 for the slow particles in the multirate forward Euler method in 3D
% fev is number of force function evaluations
% input
% N: number of cells
% (x0(j),y0(j),z0(j)) j=1:N, coordinates at t(i-1)
% mstep: number of substeps
% fev0: accumulated F evaluations
% eps: error tolerance
% s,r,ra,mu: force parameters
% output
% dt1: new time step
% K0: the set of fast particles
% k1 is evaluated by (x0,y0,z0) at t(i-1)
% fev: F evaluations at t(i)
epsdt=1e-4;
x=zeros(3,N);
% sets of slow (K1) and fast (K0) cells 
K1=zeros(N,1);
K0=ones(N,1);
for j=1:N
   x(1,j)=x0(j);
   x(2,j)=y0(j);
   x(3,j)=z0(j);
end
% compute AF, the second derivative estimate, see (25)
k1=fcmp3(N,x,s,a,ra,mu);
x=x+epsdt*k1;
F=fcmp3(N,x,s,a,ra,mu);
% difference approximation of AF
AF=(F-k1)/epsdt;
fev=fev0+2;
% compute chi for error estimate
for j=1:N
   maxAF(j)=max(abs(AF(1,j)),max(abs(AF(2,j)),abs(AF(3,j))));
end
k=1;
chia0=max(maxAF);
chia1=chia0/mstep;
% determine sets K0 and K1
for j=1:N
   if maxAF(j) < chia1
      % the set with long time steps
      K1(j)=1;
   end
end
% the set with short time steps, see (43)
K0=K0-K1;
% the long time step
dt1=sqrt(2*eps/chia1);
