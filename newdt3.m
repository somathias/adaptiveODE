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
