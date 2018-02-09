%function advection_rate.m written 2-1-18 by JTN. Takes the input string to
%determine and return the rate of advection and corresponding sigma, sigma_inv
%functions.

function phi = IC_spec(input)

    if strcmp(input,'step')
        
        phi = @(x) 5*(x>.0999).*(x<=.3);
        
        
    elseif strcmp(input,'gauss')

        phi = @(x) (exp(-(x-.2).^2/.005));
  
     elseif strcmp(input,'front')

        phi = @(x) 5*(x<=.2);
  
        
    elseif strcmp(input,'abs')
        
        phi = @(x) (1-10*abs(x-.2)).*(x<=.3).*(x>.0999);

    end
    
end