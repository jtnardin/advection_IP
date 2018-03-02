function [B1,B2] = autoreg_mat(phi1,phi2,ind1,ind2,xn)


    %autocorrelation points past shock
    if length(ind1)==1
        B1 = sparse([ind1],[ind1],1,xn,xn);
    else
        B1 = sparse([ind1 ind1(2:end) ind1(1)],[ind1 ind1(1:end-1) ind1(2)],...
            [-phi1 ones(1,length(ind1)-1) -phi1*ones(1,length(ind1)-1) 1]...
                ,xn,xn)/sqrt(1-phi1^2);
    end


    %autocorrelation points before shock
    if length(ind2) == 1
        B2 = sparse(ind2,ind2,1,xn,xn);
    else
        B2 = sparse([ind2 ind2(1:end-1) ind2(end)],[ind2 ind2(2:end) ind2(end-1)]...
            ,[ones(1,length(ind2)-1) -phi2 -phi2*ones(1,length(ind2)-1) 1]...
            ,xn,xn)/sqrt(1-phi2^2);
    end
    

end