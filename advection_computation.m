%advection_computation  written 1-29-18 by JTN to compute numerical solution
%to advection equation

function [umodel u] = advection_computation(q,g,dx,xn,x_int,xbd_0,xbd_1,dt,...
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
    elseif strcmp(num_meth,'upwindfl')
        Acomp = @(se,sw,se1,sw1) A(se,sw,v(x_int)',v(x_int-1)',x_int,1,dt/dx) + ...
            Abd(se1,sw1,v(2),v(end-1),dt/dx);
        
    else
        error('Incorrect numerical method specified')
    end
   
    
    %which computational values should we save?
    write_index = zeros(length(tc),1);
    
    for i = 1:length(tc)
        write_index(i) = find(abs(t-tc(i))==min(abs(t-tc(i))));
    end

    %initialize u
    u = zeros(length(x),length(tc));
    u(:,1) = IC;
    write_count = 1;
    
    utmp = IC';
    
    %computation over time
    
    if strcmp(num_meth,'upwindfl')
        
        sigma = @(r) (r+abs(r))./(1+abs(r));
        
        for i = 2:tn

            %create sensors
            [u_e,u_w] = positive_sensor(utmp,x_int,x_int==2,1);
            
            utmp = Acomp(sigma(u_e),sigma(u_w),sigma(u_e(2)),sigma(u_w(end-1)))*utmp;

            %write for correct indices
            if any(i==write_index)
               write_count = write_count + 1; 
               u(:,write_count) = utmp; 
            end


        end
        
    else
    
        for i = 2:tn

            utmp = Acomp*utmp;

            %write for correct indices
            if any(i==write_index)
               write_count = write_count + 1; 
               u(:,write_count) = utmp; 
            end


        end
    end
    
    
    
    %model simulation
    
    [X,TC1] = meshgrid(x,tc);
    [XC,TC2] = meshgrid(xc,tc);

    umodel = interp2(X,TC1,u',XC,TC2,'cubic');

    umodel = umodel';

    

end
