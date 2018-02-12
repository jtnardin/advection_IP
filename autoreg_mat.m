function [B1,B2] = autoreg_mat(phi1,phi2,ind1,ind2,xn)

    B1 = sparse([ind1 ind1(2:end)],[ind1 ind1(1:end-1)],[sqrt(1-phi1^2) ones(1,length(ind1)-1) -phi1*ones(1,length(ind1)-1)]...
            ,xn,xn)/sqrt(1-phi1^2);
        
    B2 = sparse([ind2 ind2(1:end-1)],[ind2 ind2(2:end)],[sqrt(1-phi2^2) ones(1,length(ind2)-1) -phi2*ones(1,length(ind2)-1)]...
            ,xn,xn)/sqrt(1-phi2^2);
end