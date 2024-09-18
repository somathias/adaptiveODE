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
    

    

    
