% Author:              DANG Truong.
% Date (dd/mm/yyyy):   01/04/2018.
% Solve the following PDE by finite volume method:
% - u_xx(x,y) - u_yy(x,y) + alpha * u(x,y) = f(x,y) with homogenous Dirichlet boundary conditions.
% i.e u = 0 (on 4 boundaries).
% Features:
% 1. Evalute mean value of f over control volume by MATLAB intrinsic function integral2.
% 2. Uniformy regular Cartesian mesh.
clear all
close all
N     = 4;
k1    = 10;
alpha = 0.0001;
if (alpha == 1)
    fun   = @(x,y)x.*y.*(x - 1).*(y - 1) - 2.*y.*(y - 1) - 2.*x.*(x - 1);
end

if (alpha == 10)
    fun   = @(x,y)10.*x.*y.*(x - 1).*(y - 1) - 2.*y.*(y - 1) - 2.*x.*(x - 1);
end

if (alpha == 0.001)
    fun   = @(x,y)0.001.*x.*y.*(x - 1).*(y - 1) - 2.*y.*(y - 1) - 2.*x.*(x - 1);
end

if (alpha == 0.0001)
    fun   = @(x,y)0.0001.*x.*y.*(x - 1).*(y - 1) - 2.*y.*(y - 1) - 2.*x.*(x - 1);
end
for j=1:4
    a    = 0;
    b    = 1;
    c    = 0;
    d    = 1;
    N    = 2*N;
    M(j) = N;
    
    %%Make mesh points x(i+1/2) 
    h    = (b-a)/(N-1);
    x(1) = a;
    for i=2:N-1
        x(i) = (i-1)*h + x(1);
    end
    x(N) = b;
    
    %%Make control points xcp(i) 
    xcp=zeros(N+1,1);
    xcp(1)=x(1);
    for i=2:N+1
        if (i==N+1)
            xcp(i)=x(i-1);
        else
            xcp(i)=(x(i)+x(i-1))/2;
        end
    end
    
    %%Make mesh points y(i+1/2) 
    k=(d-c)/(N-1);
    y(1)=c;
    for i=2:N-1
        y(i)=(i-1)*k + y(1);
    end
    y(N)=d;
    
    %%Make control points ycp(i) 
    ycp=zeros(N+1,1);
    ycp(1)=y(1);
    for i=2:N+1
        if (i==N+1)
            ycp(i)=y(i-1);
        else
            ycp(i)=(y(i)+y(i-1))/2;
        end
    end
    %%
    B=zeros((N-1)*(N-1),(N-1)*(N-1));
    
    %%
    for i = 1:N-1
        struct_matrix_A(i).A = zeros(N-1,N-1);
        for jj = 1:N-1
            ai = 1/((x(jj+1)-x(jj))*(xcp(jj+1)-xcp(jj)));
            bi = 1/((x(jj+1)-x(jj))*(xcp(jj+2)-xcp(jj+1)));
            cj = 1/((y(jj+1)-y(jj))*(ycp(jj+1)-ycp(jj)));
            dj = 1/((y(jj+1)-y(jj))*(ycp(jj+2)-ycp(jj+1)));
            sij=(ai+bi+cj+dj);
            if (jj==1)
                struct_matrix_A(i).A(jj,jj)=sij+alpha/((x(i+1)-x(i))*(y(jj+1)-y(jj)));
                struct_matrix_A(i).A(jj,jj+1)=-bi;
            else
                if(jj==N-1)
                struct_matrix_A(i).A(jj,jj)=sij+alpha/((x(i+1)-x(i))*(y(jj+1)-y(jj)));
                struct_matrix_A(i).A(jj,jj-1)=-ai;
            else
                struct_matrix_A(i).A(jj,jj)=sij+alpha/((x(i+1)-x(i))*(y(jj+1)-y(jj)));
                struct_matrix_A(i).A(jj,jj-1)=-ai;
                struct_matrix_A(i).A(jj,jj+1)=-bi;
            end
        end
        end
    end
            
    %%
    B = zeros(N-1,N-1);      
    for i=1:N-1
        if(i==1)
            B((i-1)*(N-1)+1:i*(N-1),(i-1)*(N-1)+1:i*(N-1))=struct_matrix_A(i).A;
            B((i-1)*(N-1)+1:i*(N-1),i*(N-1)+1:(i+1)*(N-1))=-1/((y(i+1)-y(i))*(ycp(i+2)-ycp(i+1))) * eye(N-1); 
        else
            if(i==N-1)
                B((i-1)*(N-1)+1:i*(N-1),(i-1)*(N-1)+1:i*(N-1))=struct_matrix_A(i).A;
                B((i-1)*(N-1)+1:i*(N-1),(i-2)*(N-1)+1:(i-1)*(N-1))=-1/((y(i+1)-y(i))*(ycp(i+1)-ycp(i))) * eye(N-1);
            else
                B((i-1)*(N-1)+1:i*(N-1),(i-1)*(N-1)+1:i*(N-1))=struct_matrix_A(i).A;
                B((i-1)*(N-1)+1:i*(N-1),i*(N-1)+1:(i+1)*(N-1))=-1/((y(i+1)-y(i))*(ycp(i+2)-ycp(i+1))) * eye(N-1);
                B((i-1)*(N-1)+1:i*(N-1),(i-2)*(N-1)+1:(i-1)*(N-1))=-1/((y(i+1)-y(i))*(ycp(i+1)-ycp(i))) * eye(N-1);
            end
        end
    end
    %%
   
    F=zeros((N-1)*(N-1),1);
    for i1 = 1:N-1
        for i2 = 1:N-1
            F((i1-1)*(N-1)+i2)=integral2(fun,x(i1),x(i1+1),y(i2),y(i2+1))/((x(i1+1)-x(i1))*(y(i2+1)-y(i2)));
        end
    end
    
    u=B\F;
    
    u_dis=zeros((N+1),(N+1));
    
     for i1=1:N-1
         for i2=1:N-1
             u_dis(i1+1,i2+1)=u((i1-1)*(N-1)+i2);
         end
     end
    
    u_exact=zeros((N+1),(N+1));
    for i1=1:N+1
        for i2=1:N+1
            u_exact(i1,i2)=exact_solution2D_2(xcp(i2),ycp(i1),k1);
        end
    end
    
    figure
    surf(xcp,ycp,u_dis)
    title('u discrete', 'FontWeight','bold')


    figure
    surf(xcp,ycp,u_exact)
    title('u exact', 'FontWeight','bold')


    %% Odd indices figures: u discrete
    %% Even indices figures: u exact
    
    %%
    norm_l2(j)=0;
    for i=1:N-1
        for t=1:N-1
            norm_l2(j)=norm_l2(j)+(u_dis(i,t)-u_exact(i,t))^2*(x(i+1)-x(i))*(y(t+1)-y(t));
        end
    end
    norm_l2(j)=(norm_l2(j))^(1/2);

    %%
    %%
    norm_h1(j)=0;
     for i=1:N
        for t=2:N-1
            norm_h1(j)=norm_h1(j)+((abs(u_dis(i+1,t)-u_exact(i+1,t))-abs(u_dis(i,t)-u_exact(i,t)))^2*(y(t+1)-y(t)))/(xcp(i+1)-xcp(i));
        end
     end
     for i=2:N-1
        for t=1:N
            norm_h1(j)=norm_h1(j)+((abs(u_dis(i,t+1)-u_exact(i,t+1))-abs(u_dis(i,t)-u_exact(i,t)))^2*(x(i+1)-x(i)))/(ycp(t+1)-ycp(t));
        end
     end
     norm_h1(j)=(norm_h1(j))^(1/2);

    %%
    
end
figure
plot( log(M), -log(norm_l2), 'red', log(M), 1*log(M),'black', log(M), -log(norm_h1), 'blue', log(M), (1/2)*log(M), 'cyan')
legend('|| ||_{L^2}','Slope = 1','|| ||_{H^1_0}','Slope = 1/2')
legend('Location','northwest')
xlabel('Log(N)')
ylabel('-Log(err)')