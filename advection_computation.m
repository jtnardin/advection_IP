%advection_computation  written 1-29-18 by JTN to compute numerical solution
%to advection equation

function umodel = advection_computation(q,g,dx,xn,x_int,xbd_0,xbd_1,dt,...
    tn,IC,A,Abd,x,xc,num_meth,t,tc)


%     q = [alpha,beta];
    
    %evaluate Velocity
    v = g(x);

    if strcmp(num_meth,'upwind')
        Acomp = A(v(x_int),v(x_int-1),x_int,1,dt/dx) + Abd;
    elseif strcmp(num_meth,'laxfried')
        Acomp = A(v(x_int+1),v(x_int-1),x_int,1,dt/dx) + Abd;
    elseif strcmp(num_meth,'laxwend')
        Acomp = A(v(x_int+1),v(x_int-1),v(x_int),x_int,1,dt/dx) + Abd;
    elseif strcmp(num_meth,'beamwarm')
        x_int_m = [x_int(2:end)];
        Acomp = A(v(x_int_m),v(x_int_m-2),v(x_int_m-1),x_int_m,1,dt/dx) + ...
            Abd(v(2),v(1),2,1,dt/dx);
    else
        error('Incorrect numerical method specified')
    end
   
    
    %which computational values should we save?
    write_index = zeros(length(tc),1);
    
    for i = 1:length(tc)
        write_index(i) = find(t==tc(i));
    end

    %initialize u
    u = zeros(xn,length(tc));
    u(:,1) = IC;
    write_count = 1;
    
    utmp = IC';
    
    %computation over time
    for i = 2:tn

        
%         u(:,i) = Acomp*u(:,i-1);

        utmp = Acomp*utmp;
        
        %write for correct indices
        if any(i==write_index)
           write_count = write_count + 1; 
           u(:,write_count) = utmp; 
        end

        
    end
    
    
    
    %model simulation
    
    [X,TC1] = meshgrid(x,tc);
    [XC,TC2] = meshgrid(xc,tc);

    umodel = interp2(X,TC1,u',XC,TC2,'cubic');

    umodel = umodel';

    

end
