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
    
