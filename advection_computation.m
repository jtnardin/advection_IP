%advection_computation  written 1-29-18 by JTN to compute numerical solution
%to advection equation

function [umodel,err] = advection_computation(q,g,dx,xn,x_int,xbd_0,xbd_1,dt,...
    tn,IC,A,x,xc,sampNo,num_meth,t,tc,soln)


%     q = [alpha,beta];
    
    %evaluate Velocity
    v = g(x);

    if strcmp(num_meth,'upwind')
        Acomp = A(v(x_int),v(x_int-1),x_int,1,dt/dx);
    elseif strcmp(num_meth,'laxfried')
        Acomp = A(v(x_int+1),v(x_int-1),x_int,1,dt/dx);
    elseif strcmp(num_meth,'laxwend')
        Acomp = A(v(x_int+1),v(x_int-1),v(x_int),x_int,1,dt/dx);
    else
        error('Incorrect numerical method specified')
    end
   


    %initialize u
    u = zeros(xn,tn);
    u(:,1) = IC;
   
    %computation over time
    for i = 2:tn

        
        u(:,i) = Acomp*u(:,i-1);

        
    end
    
    
    
    %model simulation
    
%     umodel = u(1:sampNo:end,1:sampNo*2*400/10:end);
    

    [X,T] = meshgrid(x,t);
    [XC,TC] = meshgrid(xc,tc);

    umodel = interp2(X,T,u',XC,TC,'cubic');

    umodel = umodel';
    
    % to compute grid norm
%     [T,X] = meshgrid(t,x);
%     soln_grid = soln(T,X);
%     soln_grid(isnan(soln_grid))=0;
%     
%     err = dx*dt*sum(sum(abs(u-soln_grid)));       

    err = 1;

end
