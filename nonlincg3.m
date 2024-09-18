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

  
