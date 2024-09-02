% Generation of Fig 4 in the paper 'Numerical integration of mechanical forces in
% center-based models for biological cell populations' by P Lotstedt and S Mathias
spheroid=1;
% initalize force parameters
[s,a,ra,mu]=initiasph(1);
% load center coordinates (xc00(1:nc00),yc00(1:nc00),zc00(1:nc00))
% spheroid coordinates
load coords_spheroid_217cells.txt
% number of cells
nc00=217;
xc00=coords_spheroid_217cells(:,1);
yc00=coords_spheroid_217cells(:,2);
zc00=coords_spheroid_217cells(:,3);
% dt0 is initial provisional time step
ntim=100; dt0=0.0001; kmax=10;
dte12=zeros(ntim,kmax);
% number of function evaluations for the methods
fevfEconst=zeros(kmax,1); % fEc
fev12=zeros(kmax,1); % fE
fev12m=fev12; % mfE
fev23=fev12; % RK2
fevi12=fev12; % bE
fevi23=fev12; % TM
% number of proliferations for each \Delta t_{div}
ndtmax=10;
% number of \Delta t_{div}
ndelmax=14;
% total cpu time for each \Delta t_{div} for 6 methods
cputot=zeros(ndelmax,6);
fevtot=zeros(ndelmax,6);
cput=zeros(ndtmax,6);
fev=zeros(ndtmax,6);
tedim=100*(ndtmax+1);
% local error bound
eps=0.005;
% time between proliferations, abscissa in Fig 4
tdel=zeros(ndelmax,1);
tdeltot=zeros(ndelmax,1);
tdel(1)=0.02; tdel(2)=0.05; tdel(3)=0.1; tdel(4)=0.15; tdel(5)=0.2;
tdel(6)=0.25; tdel(7)=0.3; tdel(8)=0.4; tdel(9)=0.5; tdel(10)=0.7;
tdel(11)=1.0; tdel(12)=1.3; tdel(13)=1.6; tdel(14)=2.0;
%
for kdelta=1:ndelmax
%
tdelta=zeros(ndtmax,1);
for j=1:ndtmax
   tdelta(j)=tdel(kdelta);
end
% measure cpu time for different proliferation rates 
% (xc(n,j),yc(n,j)) coordinates of n:th cell at time t_j, n=1:nc, j=1:ntim, 
% time points in t(1:ntim)
% solution with forward Euler (fE)
nte12tot=0;
% initial cell coordinates
xc0=xc00;
yc0=yc00;
zc0=zc00;
nc=nc00;
te12gl=zeros(tedim,1);
rseed=rng;
for ndt=1:ndtmax
   dtim=tdelta(ndt);
   intst=(ndt==1);
%  integrate between 0 and dtim from (xc0,yc0,zc0) to (xce12,yce12,zce12) with fE
   tic
   [xce12,yce12,zce12,te12,nte12,fev12(ndt),dt1]=exp3Dord12(ntim,dtim,nc,xc0,yc0,zc0,dt0,1,intst,eps,s,a,ra,mu);
%  cpu time for this \Delta t_{div}
   cput(ndt,1)=toc;
%  F evaluations for this \Delta t_{div}
   fev(ndt,1)=fev12(ndt);
%  a new cell is added in proliferation to final solution above, separation between cells is 0.3*s
   [nc1,xc0,yc0,zc0]=addpart(nc,xce12(:,nte12),yce12(:,nte12),zce12(:,nte12),0.3*s);
   if ndt>1
       te12end=te12gl(nte12tot);
   else
       te12end=0;
       te12first=te12;
       nte12first=nte12;
   end
%  time points for the solution
   for j=1:nte12
      te12gl(nte12tot+j)=te12(j)+te12end;
   end
   nte12tot=nte12tot+nte12;
   nc=nc1;
   dt0=dt1;
end
% determine the time steps in the fE solution
dte120=dtextr(nte12first,te12first);
nte120=nte12first-2;
% forward Euler with constant time step (fEc) and at least the same
% accuracy as fE, constant time step is chosen to be the smallest time step for adaptive fE
dtconst=min(dte120(1:nte120));
nte1consttot=0;
xc0=xc00;
yc0=yc00;
zc0=zc00;
nc=nc00;
te1constgl=zeros(tedim,1);
rng(rseed);
for ndt=1:ndtmax
   dtim=tdelta(ndt);
   intst=0;
%  integrate between 0 and dtim from (xc0,yc0,zc0) to (xce1const,yce1const,zce1const) with fEc
   tic
   [xce1const,yce1const,zce1const,te1const,nte1const,fevfEconst(ndt),dt1]=exp3Dord12(ntim,dtim,nc,xc0,yc0,zc0,dtconst,0,intst,eps,s,a,ra,mu);
%  cpu time for this \Delta t_{div}
   cput(ndt,2)=toc;
%  F evaluations for this \Delta t_{div}
   fevfEconst(ndt)=ceil(dtim/dtconst);
   fev(ndt,2)=fevfEconst(ndt);
%  a new cell is added in proliferation to final solution above
   [nc1,xc0,yc0,zc0]=addpart(nc,xce1const(:,nte1const),yce1const(:,nte1const),zce1const(:,nte1const),0.3*s);
   if ndt>1
       te1constend=te1constgl(nte1consttot);
   else
       te1constend=0;
   end
%  time points for the solution
   for j=1:nte1const
      te1constgl(nte1consttot+j)=te1const(j)+te1constend;
   end
   nte1consttot=nte1consttot+nte1const;
   nc=nc1;
end
% solution with multirate forward Euler (mfE)
nte12mtot=0;
xc0=xc00;
yc0=yc00;
zc0=zc00;
nc=nc00;
te12mgl=zeros(tedim,1);
rng(rseed);
for ndt=1:ndtmax
   dtim=tdelta(ndt);
%  integrate between 0 and dtim from (xc0,yc0,zc0) to (xce12m,yce12m,zce12m) with mfE
   tic
   [xce12m,yce12m,zce12m,te12m,nte12m,fev12m(ndt),dt1]=exp3Dord12adap(ntim,dtim,6,nc,xc0,yc0,zc0,dt0,eps,s,a,ra,mu);
%  cpu time for this \Delta t_{div}
   cput(ndt,3)=toc;
%  F evaluations for this \Delta t_{div}
   fev(ndt,3)=fev12m(ndt);
%  a new cell is added in proliferation to final solution above
   [nc1,xc0,yc0,zc0]=addpart(nc,xce12m(:,nte12m),yce12m(:,nte12m),zce12m(:,nte12m),0.3*s);
   if ndt>1
       te12mend=te12mgl(nte12mtot);
   else
       te12mend=0;
   end
%  time points for the solution
   for j=1:nte12m
      te12mgl(nte12mtot+j)=te12m(j)+te12mend;
   end
   nte12mtot=nte12mtot+nte12m;
   nc=nc1;
end
% solution with Runge-Kutta, 2nd order
nte23tot=0;
xc0=xc00;
yc0=yc00;
zc0=zc00;
nc=nc00;
te23gl=zeros(tedim,1);
rng(rseed);
for ndt=1:ndtmax
   dtim=tdelta(ndt);
   intst=(ndt==1);
%  integrate between 0 and dtim from (xc0,yc0,zc0) to (xce23,yce23,zce23) with RK2
   tic
   [xce23,yce23,zce23,te23,nte23,fev23(ndt),dt1]=exp3Dord23(ntim,dtim,nc,xc0,yc0,zc0,dt0,1,intst,eps,s,a,ra,mu);
%  cpu time for this \Delta t_{div}
   cput(ndt,4)=toc;
%  F evaluations for this \Delta t_{div}
   fev(ndt,4)=fev23(ndt);
%  a new cell is added in proliferation to final solution above
   [nc1,xc0,yc0,zc0]=addpart(nc,xce23(:,nte23),yce23(:,nte23),zce23(:,nte23),0.3*s);
   if ndt>1
       te23end=te23gl(nte23tot);
   else
       te23end=0;
   end
%  time points for the solution
   for j=1:nte23
      te23gl(nte23tot+j)=te23(j)+te23end;
   end
   nte23tot=nte23tot+nte23;
   nc=nc1;
   dt0=dt1;
end
% solution with backward Euler
nti12tot=0;
xc0=xc00;
yc0=yc00;
zc0=zc00;
nc=nc00;
ti12gl=zeros(tedim,1);
rng(rseed);
for ndt=1:ndtmax
   dtim=tdelta(ndt);
   intst=(ndt==1);
%  integrate between 0 and dtim from (xc0,yc0,zc0) to (xci12,yci12,zci12) with bE
   tic
   [xci12,yci12,zci12,ti12,nti12,iteri12,fevi12(ndt),dt1]=imp3Dord12(ntim,dtim,nc,xc0,yc0,zc0,dt0,1,intst,eps,s,a,ra,mu);
%  cpu time for this \Delta t_{div}
   cput(ndt,5)=toc;
%  F evaluations for this \Delta t_{div}
   fev(ndt,5)=fevi12(ndt);
%  a new cell is added in proliferation to final solution above
   [nc1,xc0,yc0,zc0]=addpart(nc,xci12(:,nti12),yci12(:,nti12),zci12(:,nti12),0.3*s);
   if ndt>1
       ti12end=ti12gl(nti12tot);
   else
       ti12end=0;
   end
%  time points for the solution
   for j=1:nti12
      ti12gl(nti12tot+j)=ti12(j)+ti12end;
   end
   nti12tot=nti12tot+nti12;
   nc=nc1;
   dt0=dt1;
end
% solution with trapezoidal method
nti23tot=0;
xc0=xc00;
yc0=yc00;
zc0=zc00;
nc=nc00;
ti23gl=zeros(tedim,1);
rng(rseed);
for ndt=1:ndtmax
   dtim=tdelta(ndt);
   intst=(ndt==1);
%  integrate between 0 and dtim from (xc0,yc0,zc0) to (xci23,yci23,zci23) with TM
   tic
   [xci23,yci23,zci23,ti23,nti23,iteri23,fevi23(ndt),dt1]=imp3Dord23(ntim,dtim,nc,xc0,yc0,zc0,dt0,1,intst,eps,s,a,ra,mu);
%  cpu time for this \Delta t_{div}
   cput(ndt,6)=toc;
%  F evaluations for this \Delta t_{div}
   fev(ndt,6)=fevi23(ndt);
%  a new cell is added in proliferation to final solution above
   [nc1,xc0,yc0,zc0]=addpart(nc,xci23(:,nti23),yci23(:,nti23),zci23(:,nti23),0.3*s);
   if ndt>1
       ti23end=ti23gl(nti23tot);
   else
       ti23end=0;
   end
%  time points for the solution
   for j=1:nti23
      ti23gl(nti23tot+j)=ti23(j)+ti23end;
   end
   nti23tot=nti23tot+nti23;
   nc=nc1;
   dt0=dt1;
end
% sum contribution from ndtmax intervals
[cputot1,fevtot1,tdeltot1]=cpufevsum(kdelta,ndtmax,cputot,fevtot,tdeltot,cput,fev,tdel(kdelta));
cputot=cputot1;
fevtot=fevtot1;
tdeltot=tdeltot1;
kdelta
% end of delta loop
end
%
% plot summed cpu time and force evals for each delta and ODE method
res=plotcpufev(ndelmax,cputot,fevtot,tdeltot);

function [nc1,xc1,yc1,zc1]=addpart(nc0,xc0,yc0,zc0,r0)
% initiate proliferation of cell ipart at random location
xc1=xc0;
yc1=yc0;
zc1=zc0;
% a random cell
ipart=floor(1+(nc0-0.1)*rand);
dirpart(1)=2*rand-1;
dirpart(2)=2*rand-1;
dirpart(3)=2*rand-1;
dirnorm=norm(dirpart);
dirpart(:)=0.5*r0*dirpart(:)/dirnorm;
xc1(ipart)=xc0(ipart)+dirpart(1);
yc1(ipart)=xc0(ipart)+dirpart(2);
zc1(ipart)=xc0(ipart)+dirpart(3);
nc1=nc0+1;
xc1(nc1)=xc0(ipart)-dirpart(1);
yc1(nc1)=xc0(ipart)-dirpart(2);
zc1(nc1)=xc0(ipart)-dirpart(3);

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

function [cputot1,fevtot1,tdeltot1]=cpufevsum(k,ndtmax,cputot,fevtot,tdeltot,cput,fev,tdelta)
% sum cpu time and F evaluations for each method and divide by the number of intervals between proliferations
cputot1=cputot;
fevtot1=fevtot;
tdeltot1=tdeltot;
tdeltot1(k)=tdelta;
for j=1:6
   cputot1(k,j)=sum(cput(:,j))/ndtmax;
end
for j=1:6
   fevtot1(k,j)=sum(fev(:,j))/ndtmax;
end

function dt=dtextr(nt,t)
% determine the time steps from nt time points
nt0=nt-2;
dt=zeros(nt0,1);
for j=1:nt0
   dt(j)=t(j+1)-t(j);
end

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

function F=fcmp3(N,x,s,a,ra,mu)
% force calculation for CBM in 3D
% input
% N: number of cells
% x(1:3,j), j=1:N, cell coordinates
% s,a,ra,mu: force parameters
% output
% F(1:3,j), j=1:N, force vector
    F=zeros(3,N);
    for i=1:N
        fsum=zeros(3,1);
        for j=1:N
          if i~=j
	    % diff is direction of the force between cell i and j
            diff=x(:,j)-x(:,i);
            ndiff=norm(diff);
	    % magnitude of force between two cells i and j
            fabs=force(ndiff,s,a,ra,mu);
            f=fabs*diff/ndiff;
            fsum=fsum+f;
          end
        end
	% sum of forces on cell i
        F(:,i)=fsum;
    end
     
function F=fcmpadapt3(N,K,x,s,a,ra,mu)
% force calculation for selected CBM in 3D
% input
% N: number of cells
% K: set of selected cells 
% x(1:3,j), j=1:N, cell coordinates
% s,a,ra,mu: force parameters
% output
% F(1:3,j), j=1:N, force vector
   F=zeros(3,N);
   for i=1:N
      if K(i)
        % if K(i) then compute force, otherwise 0
        fsum=zeros(3,1);
        for j=1:N
           if i~=j 
	      % diff is direction of the force between cell i and j
              diff=x(:,j)-x(:,i);
              ndiff=norm(diff);
	      % magnitude of force between two cells i and j
              fabs=force(ndiff,s,a,ra,mu);
              f=fabs*diff/ndiff;
              fsum=fsum+f;
           end
        end
        F(:,i)=fsum;
      end
   end

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

function f=force(r,s,a,ra,mu)
% the force between two cells
% input
% r: distance between cell centers
% s,a,ra,mu: force parameters in (4)
% output
% f: force magnitude
    if r<ra 
        f=mu*(r-ra)^2*(r-s);
    else
        f=0;
    end

function [xc,yc,zc,t,nt,iter,fev,dt1]=imp3Dord12(ntim,T,N,x0,y0,z0,dt0,errest,intst,eps,s,a,ra,mu)
% CBM in 3D using backward Euler (bE) and error estimate by RK2
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
% (iter(j), j=1:nt) is the number of iterations in each time step
% fev is number of function evaluations
% dt1: final dt
ntim1=ntim+1;
xc=zeros(N,ntim1);
yc=zeros(N,ntim1);
zc=zeros(N,ntim1);
t=zeros(ntim1,1);
iter=zeros(ntim1,1);
% tolerance in nonlinear iterations
epsN=1e-2*eps;
xc(:,1)=x0;
yc(:,1)=y0;
zc(:,1)=z0;
i=2;
ti=0;
fev0=0;
if errest
   % find first step, dt0 initial guess, integrate with bE
   [xc(:,i),yc(:,i),zc(:,i),k1,iter(i-1),fev1]=bE3(i,xc,yc,zc,dt0,N,fev0,epsN,s,a,ra,mu);
   fev0=fev1;
   % integrate with RK2
   [xd,yd,zd,k1,k2,fev1]=RK23(i,xc,yc,zc,1,k1,dt0,N,fev0,s,a,ra,mu);
   % determine the first time step
   dt=newdt3(xc(:,i),yc(:,i),zc(:,i),xd',yd',zd',dt0,eps,1);
   fev0=fev1;
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
   % take one step forward with bE from (xc(:,i-1),yc(:,i-1),zc(:,i-1)) to (xc(:,i),yc(:,i),zc(:,i)) with time step dt
   [xc(:,i),yc(:,i),zc(:,i),k11,iter(i),fev1]=bE3(i,xc,yc,zc,dt,N,fev0,epsN,s,a,ra,mu);
   % k11 at t(i)
   if errest && (ti<T)
      fev0=fev1;
      % take the step with RK2
      % k1 at t(i-1)
      [xd,yd,zd,k1,k2,fev1]=RK23(i,xc,yc,zc,1,k1,dt,N,fev0,s,a,ra,mu);
      k1=k11;
      % determine the new time step
      dt1=newdt3(xc(:,i),yc(:,i),zc(:,i),xd',yd',zd',dt,eps,1);
      dt=dt1;
   end
   fev0=fev1;
   i=i+1;
end
dt1=dtsav;
nt=i-1;
fev=fev1;

function [xc,yc,zc,t,nt,iter,fev,dt1]=imp3Dord23(ntim,T,N,x0,y0,z0,dt0,errest,intst,eps,s,a,ra,mu)
% CBM in 3D using trapezoidal method (TM) and error estimate by RK*
% input
% ntim: upper bound on the number of time steps
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
% (iter(j), j=1:nt) is the number of iterations in each time step
% fev is number of function evaluations
% dt1: final dt
ntim1=ntim+1;
xc=zeros(N,ntim1);
yc=zeros(N,ntim1);
zc=zeros(N,ntim1);
t=zeros(ntim1,1);
iter=zeros(ntim1,1);
% tolerance in nonlinear iterations
epsN=1e-2*eps;
xc(:,1)=x0;
yc(:,1)=y0;
zc(:,1)=z0;
i=2;
ti=0;
fev0=0;
if errest
   % find first step, dt0 initial guess, integrate with TM
   [xc(:,i),yc(:,i),zc(:,i),k1,k2,iter(i-1),fev1]=TM3(i,xc,yc,zc,dt0,N,fev0,epsN,s,a,ra,mu);
   fev0=fev1;
   % integrate with RK*
   [xd,yd,zd,fev1]=RKS3(i,xc,yc,zc,k1,k2,dt0,N,fev0,s,a,ra,mu);
   % determine the first time step
   dt=newdt3(xc(:,i),yc(:,i),zc(:,i),xd',yd',zd',dt0,eps,2);
   fev0=fev1;
else
   dt=dt0;
end
dtsav=dt;
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
   % take one step forward with TM from (xc(:,i-1),yc(:,i-1),zc(:,i-1)) to (xc(:,i),yc(:,i),zc(:,i)) with time step dt
   [xc(:,i),yc(:,i),zc(:,i),k1,k2,iter(i),fev1]=TM3(i,xc,yc,zc,dt,N,fev0,epsN,s,a,ra,mu);
   if errest && (ti<T)
      fev0=fev1;
      % take the step with RK*
      [xd,yd,zd,fev1]=RKS3(i,xc,yc,zc,k1,k2,dt,N,fev0,s,a,ra,mu);
      % determine the new time step
      dt1=newdt3(xc(:,i),yc(:,i),zc(:,i),xd',yd',zd',dt,eps,2);
      dt=dt1;
   end
   fev0=fev1;
   i=i+1;
end
dt1=dtsav;
nt=i-1;
fev=fev1;

function [s,a,ra,mu]=initiasph(N)
% initialize constants for CBM spheroid
    s=1;
    a=10;
    ra=1.5;
    mu=5.7;

function dt=newdt3(x0,y0,z0,x1,y1,z1,dt0,eps,iord)
% compute new time step by comparing (x0,y0,z0) with (x1,y1,z1) in max norm
% input
% (x0(j),y0(j),z0(j)) j=1:N, method of low order  
% (x1(j),y1(j),z1(j)) j=1:N, method of high order  
% dt0: previous time step
% eps: error tolerance
% iord: order of low order method
% output
% dt: new time step
xerr=max(abs(x0-x1));
yerr=max(abs(y0-y1));
zerr=max(abs(z0-z1));
err=max(xerr,max(yerr,zerr));
if iord==1
   dt=dt0*sqrt(eps/err);
else
   dt=dt0*(eps/err)^(1/3);
end

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

function [x1,f1,iter,fev]=nonlincg3(N,x0,f0,dt,c,fev0,eps,s,a,ra,mu)
% solve nonlinf3 by nonlinear conjugate gradient iterations
% to minimize F(x) with gradient f(x), theory according to Hager and Zhang
% input
% N: number of cells
% x0(1:3,j), j=1:N, cell coordinates
% f0(1:3,j), j=1:N, constant right hand side, b in (72), (73)
% c=1: backward Euler, c=1/2: trapezoidal method
% fev0: accumulated number of force function evaluations
% eps: iteration tolerance
% s,a,ra,mu: force parameters
% output
% x1(1:3,j), j=1:N, solution
% f1: the force at x1
% iter: number of iterations
% fev: updated F evaluations
epscg=1e-3;
p=zeros(3,N);
N2=2*N;
x=x0;
r0=nonlinf3(N,x,f0,dt,c,s,a,ra,mu);
p0=r0;
iter=0;
while (norm(r0)>eps)&&(iter<20)
% iterate until the l2 norm is sufficiently small or the number of iterations is 20
   iter=iter+1;
   for j=1:N
      p(1,j)=p0(j);
      p(2,j)=p0(N+j);
      p(3,j)=p0(N2+j);
   end
   % approximate Jacobian A times p (A is Hessian of F)
   x1=x+epscg*p;
   r=nonlinf3(N,x1,f0,dt,c,s,a,ra,mu);
   Ap=(r-r0)/epscg;
   % minimize F(x) with alpha, F(x1)=F(x-alpha*p)<F(x), assuming
   % quadratic approximation of F around x
   alpha=r0'*r0/(p0'*Ap);
   x1=x-alpha*p;
   r1=nonlinf3(N,x1,f0,dt,c,s,a,ra,mu);
   % beta according to Fletcher-Reeves [16]
   beta=r1'*r1/(r0'*r0);
   p0=r1+beta*p0;
   r0=r1;
   x=x1;
end
x1=x;
% the force at the final x
f1=fcmp3(N,x1,s,a,ra,mu);
fev=fev0+1+2*iter;

function [nlf]=nonlinf3(N,x0,f0,dt,c,s,a,ra,mu)
% compute the nonlinear function to be solved for x^{n+1} in bE and TM
% input
% N: number of cells
% x0(1:3,j), j=1:N, cell coordinates
% f0(1:3,j), j=1:N, constant right hand side in function
% dt: time step
% c: bE: c=1, TM: c=0.5
% s,a,ra,mu: force parameters
% output
% nlf(k), k=1:(3N) value of nonlinear function
nlf=zeros(3*N,1);
% compute the force at x0
frc=fcmp3(N,x0,s,a,ra,mu);
% value of nonlinear function
nlf0=x0-c*dt*frc-f0;
k=N;
l=2*N;
% reformat the output
for j=1:N
   k=k+1;
   l=l+1;
   nlf(j)=nlf0(1,j);
   nlf(k)=nlf0(2,j);
   nlf(l)=nlf0(3,j);
end

function res=plotcpufev(ndelmax,cputot,fevtot,tdeltot)
% plot log10 of cpu time for 5 methods divided by the cpu time for fEc 
cputot6=cputot(:,2);
fevtot6=fevtot(:,2);
% scaling by data for constant time step
for j=1:ndelmax
   for k=1:6
      cputot2(j,k)=cputot(j,k)/cputot6(j);
      fevtot2(j,k)=fevtot(j,k)/fevtot6(j);
   end
end
%subplot(2,1,1)
% Fig 4
plot(tdeltot(1:ndelmax),log10(cputot2(1:ndelmax,1)),'g*', ...
    tdeltot(1:ndelmax),log10(cputot2(1:ndelmax,2)),'k-.',tdeltot(1:ndelmax),log10(cputot2(1:ndelmax,3)),'bo', ...
     tdeltot(1:ndelmax),log10(cputot2(1:ndelmax,4)),'m^',tdeltot(1:ndelmax),log10(cputot2(1:ndelmax,5)),'c+', ...
     tdeltot(1:ndelmax),log10(cputot2(1:ndelmax,6)),'rx')
% plot(tdeltot(1:ndelmax),cputot2(1:ndelmax,1),'g*', ...
%     tdeltot(1:ndelmax),cputot2(1:ndelmax,2),'k-.',tdeltot(1:ndelmax),cputot2(1:ndelmax,3),'bo', ...
%      tdeltot(1:ndelmax),cputot2(1:ndelmax,4),'m^',tdeltot(1:ndelmax),cputot2(1:ndelmax,5),'c+', ...
%      tdeltot(1:ndelmax),cputot2(1:ndelmax,6),'rx')
%subplot(2,1,2)
%plot(tdeltot(1:ndelmax),log10(fevtot2(1:ndelmax,1)),'g*',tdeltot(1:ndelmax),log10(fevtot2(1:ndelmax,2)),'k-.', ...
%    tdeltot(1:ndelmax),log10(fevtot2(1:ndelmax,3)),'bo', ...
%	tdeltot(1:ndelmax),log10(fevtot2(1:ndelmax,4)),'m^',tdeltot(1:ndelmax),log10(fevtot2(1:ndelmax,5)),'c+', ...
%    tdeltot(1:ndelmax),log10(fevtot2(1:ndelmax,6)),'rx')
res=0;

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

% initial coordinates for 217 cells in a spheroid tissue
0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00
5.000000000000000000e-01 2.886751345948128655e-01 8.164965809277259234e-01
0.000000000000000000e+00 0.000000000000000000e+00 1.632993161855451847e+00
5.000000000000000000e-01 2.886751345948128655e-01 2.449489742783177881e+00
0.000000000000000000e+00 0.000000000000000000e+00 3.265986323710903694e+00
5.000000000000000000e-01 2.886751345948128655e-01 4.082482904638630394e+00
5.000000000000000000e-01 8.660254037844385966e-01 0.000000000000000000e+00
0.000000000000000000e+00 1.154700538379251462e+00 8.164965809277259234e-01
5.000000000000000000e-01 8.660254037844385966e-01 1.632993161855451847e+00
0.000000000000000000e+00 1.154700538379251462e+00 2.449489742783177881e+00
5.000000000000000000e-01 8.660254037844385966e-01 3.265986323710903694e+00
0.000000000000000000e+00 1.154700538379251462e+00 4.082482904638630394e+00
0.000000000000000000e+00 1.732050807568877193e+00 0.000000000000000000e+00
5.000000000000000000e-01 2.020725942163690281e+00 8.164965809277259234e-01
0.000000000000000000e+00 1.732050807568877193e+00 1.632993161855451847e+00
5.000000000000000000e-01 2.020725942163690281e+00 2.449489742783177881e+00
0.000000000000000000e+00 1.732050807568877193e+00 3.265986323710903694e+00
5.000000000000000000e-01 2.020725942163690281e+00 4.082482904638630394e+00
5.000000000000000000e-01 2.598076211353316012e+00 0.000000000000000000e+00
0.000000000000000000e+00 2.886751345948128655e+00 8.164965809277259234e-01
5.000000000000000000e-01 2.598076211353316012e+00 1.632993161855451847e+00
0.000000000000000000e+00 2.886751345948128655e+00 2.449489742783177881e+00
5.000000000000000000e-01 2.598076211353316012e+00 3.265986323710903694e+00
0.000000000000000000e+00 2.886751345948128655e+00 4.082482904638630394e+00
0.000000000000000000e+00 3.464101615137754386e+00 0.000000000000000000e+00
5.000000000000000000e-01 3.752776749732567030e+00 8.164965809277259234e-01
0.000000000000000000e+00 3.464101615137754386e+00 1.632993161855451847e+00
5.000000000000000000e-01 3.752776749732567030e+00 2.449489742783177881e+00
0.000000000000000000e+00 3.464101615137754386e+00 3.265986323710903694e+00
5.000000000000000000e-01 3.752776749732567030e+00 4.082482904638630394e+00
5.000000000000000000e-01 4.330127018922192761e+00 0.000000000000000000e+00
0.000000000000000000e+00 4.618802153517005848e+00 8.164965809277259234e-01
5.000000000000000000e-01 4.330127018922192761e+00 1.632993161855451847e+00
0.000000000000000000e+00 4.618802153517005848e+00 2.449489742783177881e+00
5.000000000000000000e-01 4.330127018922192761e+00 3.265986323710903694e+00
0.000000000000000000e+00 4.618802153517005848e+00 4.082482904638630394e+00
1.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00
1.500000000000000000e+00 2.886751345948128655e-01 8.164965809277259234e-01
1.000000000000000000e+00 0.000000000000000000e+00 1.632993161855451847e+00
1.500000000000000000e+00 2.886751345948128655e-01 2.449489742783177881e+00
1.000000000000000000e+00 0.000000000000000000e+00 3.265986323710903694e+00
1.500000000000000000e+00 2.886751345948128655e-01 4.082482904638630394e+00
1.500000000000000000e+00 8.660254037844385966e-01 0.000000000000000000e+00
1.000000000000000000e+00 1.154700538379251462e+00 8.164965809277259234e-01
1.500000000000000000e+00 8.660254037844385966e-01 1.632993161855451847e+00
1.000000000000000000e+00 1.154700538379251462e+00 2.449489742783177881e+00
1.500000000000000000e+00 8.660254037844385966e-01 3.265986323710903694e+00
1.000000000000000000e+00 1.154700538379251462e+00 4.082482904638630394e+00
1.000000000000000000e+00 1.732050807568877193e+00 0.000000000000000000e+00
1.500000000000000000e+00 2.020725942163690281e+00 8.164965809277259234e-01
1.000000000000000000e+00 1.732050807568877193e+00 1.632993161855451847e+00
1.500000000000000000e+00 2.020725942163690281e+00 2.449489742783177881e+00
1.000000000000000000e+00 1.732050807568877193e+00 3.265986323710903694e+00
1.500000000000000000e+00 2.020725942163690281e+00 4.082482904638630394e+00
1.500000000000000000e+00 2.598076211353316012e+00 0.000000000000000000e+00
1.000000000000000000e+00 2.886751345948128655e+00 8.164965809277259234e-01
1.500000000000000000e+00 2.598076211353316012e+00 1.632993161855451847e+00
1.000000000000000000e+00 2.886751345948128655e+00 2.449489742783177881e+00
1.500000000000000000e+00 2.598076211353316012e+00 3.265986323710903694e+00
1.000000000000000000e+00 2.886751345948128655e+00 4.082482904638630394e+00
1.000000000000000000e+00 3.464101615137754386e+00 0.000000000000000000e+00
1.500000000000000000e+00 3.752776749732567030e+00 8.164965809277259234e-01
1.000000000000000000e+00 3.464101615137754386e+00 1.632993161855451847e+00
1.500000000000000000e+00 3.752776749732567030e+00 2.449489742783177881e+00
1.000000000000000000e+00 3.464101615137754386e+00 3.265986323710903694e+00
1.500000000000000000e+00 3.752776749732567030e+00 4.082482904638630394e+00
1.500000000000000000e+00 4.330127018922192761e+00 0.000000000000000000e+00
1.000000000000000000e+00 4.618802153517005848e+00 8.164965809277259234e-01
1.500000000000000000e+00 4.330127018922192761e+00 1.632993161855451847e+00
1.000000000000000000e+00 4.618802153517005848e+00 2.449489742783177881e+00
1.500000000000000000e+00 4.330127018922192761e+00 3.265986323710903694e+00
1.000000000000000000e+00 4.618802153517005848e+00 4.082482904638630394e+00
2.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00
2.500000000000000000e+00 2.886751345948128655e-01 8.164965809277259234e-01
2.000000000000000000e+00 0.000000000000000000e+00 1.632993161855451847e+00
2.500000000000000000e+00 2.886751345948128655e-01 2.449489742783177881e+00
2.000000000000000000e+00 0.000000000000000000e+00 3.265986323710903694e+00
2.500000000000000000e+00 2.886751345948128655e-01 4.082482904638630394e+00
2.500000000000000000e+00 8.660254037844385966e-01 0.000000000000000000e+00
2.000000000000000000e+00 1.154700538379251462e+00 8.164965809277259234e-01
2.500000000000000000e+00 8.660254037844385966e-01 1.632993161855451847e+00
2.000000000000000000e+00 1.154700538379251462e+00 2.449489742783177881e+00
2.500000000000000000e+00 8.660254037844385966e-01 3.265986323710903694e+00
2.000000000000000000e+00 1.154700538379251462e+00 4.082482904638630394e+00
2.000000000000000000e+00 1.732050807568877193e+00 0.000000000000000000e+00
2.500000000000000000e+00 2.020725942163690281e+00 8.164965809277259234e-01
2.000000000000000000e+00 1.732050807568877193e+00 1.632993161855451847e+00
2.500000000000000000e+00 2.020725942163690281e+00 2.449489742783177881e+00
2.000000000000000000e+00 1.732050807568877193e+00 3.265986323710903694e+00
2.500000000000000000e+00 2.020725942163690281e+00 4.082482904638630394e+00
2.500000000000000000e+00 2.598076211353316012e+00 0.000000000000000000e+00
2.000000000000000000e+00 2.886751345948128655e+00 8.164965809277259234e-01
2.500000000000000000e+00 2.598076211353316012e+00 1.632993161855451847e+00
2.000000000000000000e+00 2.886751345948128655e+00 2.449489742783177881e+00
2.500000000000000000e+00 2.598076211353316012e+00 3.265986323710903694e+00
2.000000000000000000e+00 2.886751345948128655e+00 4.082482904638630394e+00
2.000000000000000000e+00 3.464101615137754386e+00 0.000000000000000000e+00
2.500000000000000000e+00 3.752776749732567030e+00 8.164965809277259234e-01
2.000000000000000000e+00 3.464101615137754386e+00 1.632993161855451847e+00
2.500000000000000000e+00 3.752776749732567030e+00 2.449489742783177881e+00
2.000000000000000000e+00 3.464101615137754386e+00 3.265986323710903694e+00
2.500000000000000000e+00 3.752776749732567030e+00 4.082482904638630394e+00
2.500000000000000000e+00 4.330127018922192761e+00 0.000000000000000000e+00
2.000000000000000000e+00 4.618802153517005848e+00 8.164965809277259234e-01
2.500000000000000000e+00 4.330127018922192761e+00 1.632993161855451847e+00
2.000000000000000000e+00 4.618802153517005848e+00 2.449489742783177881e+00
2.500000000000000000e+00 4.330127018922192761e+00 3.265986323710903694e+00
2.000000000000000000e+00 4.618802153517005848e+00 4.082482904638630394e+00
3.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00
3.500000000000000000e+00 2.886751345948128655e-01 8.164965809277259234e-01
3.000000000000000000e+00 0.000000000000000000e+00 1.632993161855451847e+00
3.500000000000000000e+00 2.886751345948128655e-01 2.449489742783177881e+00
3.000000000000000000e+00 0.000000000000000000e+00 3.265986323710903694e+00
3.500000000000000000e+00 2.886751345948128655e-01 4.082482904638630394e+00
3.500000000000000000e+00 8.660254037844385966e-01 0.000000000000000000e+00
3.000000000000000000e+00 1.154700538379251462e+00 8.164965809277259234e-01
3.500000000000000000e+00 8.660254037844385966e-01 1.632993161855451847e+00
3.000000000000000000e+00 1.154700538379251462e+00 2.449489742783177881e+00
3.500000000000000000e+00 8.660254037844385966e-01 3.265986323710903694e+00
3.000000000000000000e+00 1.154700538379251462e+00 4.082482904638630394e+00
3.000000000000000000e+00 1.732050807568877193e+00 0.000000000000000000e+00
3.500000000000000000e+00 2.020725942163690281e+00 8.164965809277259234e-01
3.000000000000000000e+00 1.732050807568877193e+00 1.632993161855451847e+00
3.500000000000000000e+00 2.020725942163690281e+00 2.449489742783177881e+00
3.000000000000000000e+00 1.732050807568877193e+00 3.265986323710903694e+00
3.500000000000000000e+00 2.020725942163690281e+00 4.082482904638630394e+00
3.500000000000000000e+00 2.598076211353316012e+00 0.000000000000000000e+00
3.000000000000000000e+00 2.886751345948128655e+00 8.164965809277259234e-01
3.500000000000000000e+00 2.598076211353316012e+00 1.632993161855451847e+00
2.868896213024025421e+00 2.851572818896319195e+00 2.385659747678811726e+00
3.500000000000000000e+00 2.598076211353316012e+00 3.265986323710903694e+00
3.000000000000000000e+00 2.886751345948128655e+00 4.082482904638630394e+00
3.000000000000000000e+00 3.464101615137754386e+00 0.000000000000000000e+00
3.500000000000000000e+00 3.752776749732567030e+00 8.164965809277259234e-01
3.000000000000000000e+00 3.464101615137754386e+00 1.632993161855451847e+00
3.500000000000000000e+00 3.752776749732567030e+00 2.449489742783177881e+00
3.000000000000000000e+00 3.464101615137754386e+00 3.265986323710903694e+00
3.500000000000000000e+00 3.752776749732567030e+00 4.082482904638630394e+00
3.500000000000000000e+00 4.330127018922192761e+00 0.000000000000000000e+00
3.000000000000000000e+00 4.618802153517005848e+00 8.164965809277259234e-01
3.500000000000000000e+00 4.330127018922192761e+00 1.632993161855451847e+00
3.000000000000000000e+00 4.618802153517005848e+00 2.449489742783177881e+00
3.500000000000000000e+00 4.330127018922192761e+00 3.265986323710903694e+00
3.000000000000000000e+00 4.618802153517005848e+00 4.082482904638630394e+00
4.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00
4.500000000000000000e+00 2.886751345948128655e-01 8.164965809277259234e-01
4.000000000000000000e+00 0.000000000000000000e+00 1.632993161855451847e+00
4.500000000000000000e+00 2.886751345948128655e-01 2.449489742783177881e+00
4.000000000000000000e+00 0.000000000000000000e+00 3.265986323710903694e+00
4.500000000000000000e+00 2.886751345948128655e-01 4.082482904638630394e+00
4.500000000000000000e+00 8.660254037844385966e-01 0.000000000000000000e+00
4.000000000000000000e+00 1.154700538379251462e+00 8.164965809277259234e-01
4.500000000000000000e+00 8.660254037844385966e-01 1.632993161855451847e+00
4.000000000000000000e+00 1.154700538379251462e+00 2.449489742783177881e+00
4.500000000000000000e+00 8.660254037844385966e-01 3.265986323710903694e+00
4.000000000000000000e+00 1.154700538379251462e+00 4.082482904638630394e+00
4.000000000000000000e+00 1.732050807568877193e+00 0.000000000000000000e+00
4.500000000000000000e+00 2.020725942163690281e+00 8.164965809277259234e-01
4.000000000000000000e+00 1.732050807568877193e+00 1.632993161855451847e+00
4.500000000000000000e+00 2.020725942163690281e+00 2.449489742783177881e+00
4.000000000000000000e+00 1.732050807568877193e+00 3.265986323710903694e+00
4.500000000000000000e+00 2.020725942163690281e+00 4.082482904638630394e+00
4.500000000000000000e+00 2.598076211353316012e+00 0.000000000000000000e+00
4.000000000000000000e+00 2.886751345948128655e+00 8.164965809277259234e-01
4.500000000000000000e+00 2.598076211353316012e+00 1.632993161855451847e+00
4.000000000000000000e+00 2.886751345948128655e+00 2.449489742783177881e+00
4.500000000000000000e+00 2.598076211353316012e+00 3.265986323710903694e+00
4.000000000000000000e+00 2.886751345948128655e+00 4.082482904638630394e+00
4.000000000000000000e+00 3.464101615137754386e+00 0.000000000000000000e+00
4.500000000000000000e+00 3.752776749732567030e+00 8.164965809277259234e-01
4.000000000000000000e+00 3.464101615137754386e+00 1.632993161855451847e+00
4.500000000000000000e+00 3.752776749732567030e+00 2.449489742783177881e+00
4.000000000000000000e+00 3.464101615137754386e+00 3.265986323710903694e+00
4.500000000000000000e+00 3.752776749732567030e+00 4.082482904638630394e+00
4.500000000000000000e+00 4.330127018922192761e+00 0.000000000000000000e+00
4.000000000000000000e+00 4.618802153517005848e+00 8.164965809277259234e-01
4.500000000000000000e+00 4.330127018922192761e+00 1.632993161855451847e+00
4.000000000000000000e+00 4.618802153517005848e+00 2.449489742783177881e+00
4.500000000000000000e+00 4.330127018922192761e+00 3.265986323710903694e+00
4.000000000000000000e+00 4.618802153517005848e+00 4.082482904638630394e+00
5.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00
5.500000000000000000e+00 2.886751345948128655e-01 8.164965809277259234e-01
5.000000000000000000e+00 0.000000000000000000e+00 1.632993161855451847e+00
5.500000000000000000e+00 2.886751345948128655e-01 2.449489742783177881e+00
5.000000000000000000e+00 0.000000000000000000e+00 3.265986323710903694e+00
5.500000000000000000e+00 2.886751345948128655e-01 4.082482904638630394e+00
5.500000000000000000e+00 8.660254037844385966e-01 0.000000000000000000e+00
5.000000000000000000e+00 1.154700538379251462e+00 8.164965809277259234e-01
5.500000000000000000e+00 8.660254037844385966e-01 1.632993161855451847e+00
5.000000000000000000e+00 1.154700538379251462e+00 2.449489742783177881e+00
5.500000000000000000e+00 8.660254037844385966e-01 3.265986323710903694e+00
5.000000000000000000e+00 1.154700538379251462e+00 4.082482904638630394e+00
5.000000000000000000e+00 1.732050807568877193e+00 0.000000000000000000e+00
5.500000000000000000e+00 2.020725942163690281e+00 8.164965809277259234e-01
5.000000000000000000e+00 1.732050807568877193e+00 1.632993161855451847e+00
5.500000000000000000e+00 2.020725942163690281e+00 2.449489742783177881e+00
5.000000000000000000e+00 1.732050807568877193e+00 3.265986323710903694e+00
5.500000000000000000e+00 2.020725942163690281e+00 4.082482904638630394e+00
5.500000000000000000e+00 2.598076211353316012e+00 0.000000000000000000e+00
5.000000000000000000e+00 2.886751345948128655e+00 8.164965809277259234e-01
5.500000000000000000e+00 2.598076211353316012e+00 1.632993161855451847e+00
5.000000000000000000e+00 2.886751345948128655e+00 2.449489742783177881e+00
5.500000000000000000e+00 2.598076211353316012e+00 3.265986323710903694e+00
5.000000000000000000e+00 2.886751345948128655e+00 4.082482904638630394e+00
5.000000000000000000e+00 3.464101615137754386e+00 0.000000000000000000e+00
5.500000000000000000e+00 3.752776749732567030e+00 8.164965809277259234e-01
5.000000000000000000e+00 3.464101615137754386e+00 1.632993161855451847e+00
5.500000000000000000e+00 3.752776749732567030e+00 2.449489742783177881e+00
5.000000000000000000e+00 3.464101615137754386e+00 3.265986323710903694e+00
5.500000000000000000e+00 3.752776749732567030e+00 4.082482904638630394e+00
5.500000000000000000e+00 4.330127018922192761e+00 0.000000000000000000e+00
5.000000000000000000e+00 4.618802153517005848e+00 8.164965809277259234e-01
5.500000000000000000e+00 4.330127018922192761e+00 1.632993161855451847e+00
5.000000000000000000e+00 4.618802153517005848e+00 2.449489742783177881e+00
5.500000000000000000e+00 4.330127018922192761e+00 3.265986323710903694e+00
5.000000000000000000e+00 4.618802153517005848e+00 4.082482904638630394e+00
3.131103786975974579e+00 2.921929872999938116e+00 2.513319737887544036e+00

  
