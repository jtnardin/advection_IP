%function advection_rate.m written 2-1-18 by JTN. Takes the input string to
%determine and return the rate of advection and corresponding sigma, sigma_inv
%functions.

function [g,sigma,sigma_inv,galpha,gbeta] = advection_rate(input,alpha,beta)

    if strcmp(input,'root')
        
        g = @(x) alpha*x.^(1/beta);
        galpha = @(x) x.^(1/beta);
        gbeta = @(x) alpha/beta*x.^(1/beta-1);

        %sigma, sigma_inv terms
        sigma = @(x,x0) 1/(alpha*(1-1/beta))*(x.^(1-1/beta)-x0.^(1-1/beta));
        sigma_inv = @(t,x0) (alpha*(1-1/beta)*t+x0.^(1-1/beta)).^(beta/(beta-1));
        
    elseif strcmp(input,'const')



        %rate of advection
        g = @(x) alpha*ones(size(x));

        %sigma, sigma_inv terms
        sigma = @(x,x0) (x-x0)/alpha;
        sigma_inv = @(t,x0) x0 + alpha*t;

    end
    
end