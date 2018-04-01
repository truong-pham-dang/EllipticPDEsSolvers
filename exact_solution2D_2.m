function u_ex=exact_solution2D_2(x,y,k)

% Homogenous Dirichlet bc

if (k==1)  
    u_ex=x*y*(x - 1)*(y - 1);
end

if (k == 10)
    u_ex=x*y*(x - 1)*(y - 1);
end

% Nonhomogenous Dirichlet bc

%u_ex=(x-1/3)^2*(y-2/3)^2;
