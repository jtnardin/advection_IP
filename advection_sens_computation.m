%advection_computation  written 1-29-18 by JTN to compute numerical solution
%to advection equation

function [umodel,s1model,s2model] = advection_sens_computation(g,galpha,gbeta,dx,xn,x_int,xbd_0,xbd_1,dt,...
    tn,IC,A,Abd,x,xc,num_meth,t,tc)


%     q = [alpha,beta];
    
    %evaluate Velocity
    v = g(x);
    vs1 = galpha(x);
    vs2 = gbeta(x);
    
    
    if strcmp(num_meth,'upwindsens')
       Acomp = A(v(x_int),v(x_int-1),vs1(x_int),vs1(x_int-1),vs2(x_int),vs2(x_int-1),x_int,1,dt/dx) + ...
            Abd;
    elseif strcmp(num_meth,'laxfriedsens')
       Acomp = A(v(x_int+1),v(x_int-1),vs1(x_int+1),vs1(x_int-1),vs2(x_int+1),vs2(x_int-1),x_int,1,dt/dx) + ...
            Abd;
    elseif strcmp(num_meth,'laxwendsens')
       Acomp = A(v(x_int+1),v(x_int-1),v(x_int),vs1(x_int+1),vs1(x_int-1),...
           vs1(x_int),vs2(x_int+1),vs2(x_int-1),vs2(x_int),x_int,1,dt/dx) + ...
            Abd;
    elseif strcmp(num_meth,'beamwarmsens')
       x_int_m = x_int(2:end);
       Acomp = A(v(x_int_m+1),v(x_int_m-1),v(x_int_m),vs1(x_int_m+1),vs1(x_int_m-1),...
           vs1(x_int_m),vs2(x_int_m+1),vs2(x_int_m-1),vs2(x_int_m),x_int_m,1,dt/dx) + ...
            Abd(v(2),v(1),vs1(2),vs1(1),vs2(2),vs2(1),2,1,dt/dx);
    else
        error('Incorrect numerical method specified')
    end
   
    
    %which computational values should we save?
    write_index = zeros(length(tc),1);
    
    for i = 1:length(tc)
        write_index(i) = find(abs(t-tc(i))==min(abs(t-tc(i))));
    end

    %initialize u
    u = zeros(3*xn,length(tc));
    u(:,1) = IC;
    write_count = 1;
    
    utmp = IC';
    
%     figure
    
    %computation over time
    for i = 2:tn

        
%         u(:,i) = Acomp*u(:,i-1);

        utmp = Acomp*utmp;
        
        %write for correct indices
        if any(i==write_index)
           write_count = write_count + 1; 
           u(:,write_count) = utmp; 
        end
        
%         if mod(i,250) == 0
%            disp('hey') 
%         end

%         if mod(i,10)==0
%            subplot(3,1,1)
%            plot(x,utmp(1:xn))
%            subplot(3,1,2)
%            plot(x,utmp(xn+1:2*xn))
%            subplot(3,1,3)
%            plot(x,utmp(2*xn+1:3*xn))
%            
%            title(num2str(t(i)))
%            pause(.125)
%            
%         end
        
    end
    
    v = u(1:xn,:);
    s1 = u(xn+1:2*xn,:);
    s2 = u(2*xn+1:3*xn,:);
    
    %model simulation
    
    [X,TC1] = meshgrid(x,tc);
    [XC,TC2] = meshgrid(xc,tc);

    umodel = interp2(X,TC1,v',XC,TC2,'cubic');
    
    s1model = interp2(X,TC1,s1',XC,TC2,'cubic');
    
    s2model = interp2(X,TC1,s2',XC,TC2,'cubic');
    
    

end
