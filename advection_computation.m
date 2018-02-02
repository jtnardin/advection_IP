%advection_computation  written 1-29-18 by JTN to compute numerical solution
%to advection equation

function umodel = advection_computation(q,g,dx,xn,x_int,xbd_0,xbd_1,dt,...
    tn,IC,A,Abd,x,xc,num_meth,t,tc)


%     q = [alpha,beta];
    
    %evaluate Velocity
    v = g(x);

    if strcmp(num_meth,'upwind')
        Acomp = A(v(x_int),v(x_int-1),x_int,1,dt/dx);
    elseif strcmp(num_meth,'laxfried')
        Acomp = A(v(x_int+1),v(x_int-1),x_int,1,dt/dx);
    elseif strcmp(num_meth,'laxwend')
        Acomp = A(v(x_int+1),v(x_int-1),v(x_int),x_int,1,dt/dx);
    elseif strcmp(num_meth,'beamwarm')
        x_int_m = [x_int(2:end) xn];
        Acomp = A(v(x_int_m),v(x_int_m-2),v(x_int_m-1),x_int_m,1,dt/dx) + ...
            Abd(v(2),v(1),2,1,dt/dx);
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
    
    [X,T] = meshgrid(x,t);
    [XC,TC] = meshgrid(xc,tc);

    umodel = interp2(X,T,u',XC,TC,'cubic');

    umodel = umodel';


end
