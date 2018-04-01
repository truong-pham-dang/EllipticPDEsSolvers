syms x y;
u = x*(x-1)*y*(y-1)                 % Homogenous Dirichlet bc
%u  = (x-1/3)^2*(y-2/3)^2            % Nonhomogenous Dirichlet bc
diffx = diff(u,x)
diffxx = diff(diffx,x)
diffy = diff(u,y) 
diffyy = diff(diffy,y)
-(diffxx+diffyy)+1*u
 
