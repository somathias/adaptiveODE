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
