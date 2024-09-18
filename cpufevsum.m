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
