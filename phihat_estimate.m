function phi = phihat_estimate(res)

    phi = sum(res(1:end-1).*res(2:end))/sum(res.^2);
    
    phi(isnan(phi)) = 0;

end