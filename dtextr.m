function dt=dtextr(nt,t)
% determine the time steps from nt time points
nt0=nt-2;
dt=zeros(nt0,1);
for j=1:nt0
   dt(j)=t(j+1)-t(j);
end

