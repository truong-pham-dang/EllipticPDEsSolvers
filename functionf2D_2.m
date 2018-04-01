function f=functionf2D_2(x,y,k)

% Homogenous Dirichlet bc

if (k == 1)
    f=x*y*(x - 1)*(y - 1) - 2*y*(y - 1) - 2*x*(x - 1);
end

if (k == 10)
    f=10*x*y*(x - 1)*(y - 1) - 2*y*(y - 1) - 2*x*(x - 1);
end

% Nonhomogenous Dirichlet bc

%f=(-2)*(y-2/3)^2-2*(x-1/3)^2+(x-1/3)^2*(y-2/3)^2;