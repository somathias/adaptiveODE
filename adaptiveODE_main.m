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
