function [A,Abd] = aMatrixupwind(xn,method)


    if strcmp(method,'upwind')
        %create upwind matrix
        A = @(ve,vw,ind,dn,l) sparse([ind ind],[ind ind-dn],[1-l*ve l*vw],xn,xn);
        
        Abd = sparse([1 xn],[1 xn],[1 1],xn,xn);
    elseif strcmp(method,'laxfried')
        %create lax-friedrich method
        A = @(ve,vw,ind,dn,l) sparse([ind ind],[ind-dn ind+dn],...
            [1+l*vw 1-l*ve]/2,xn,xn);
        
        Abd = sparse([1 xn],[1 xn],[1 1],xn,xn);
    elseif strcmp(method,'laxwend')
        %create lax-wendroff method
        A = @(ve,vw,vp,ind,dn,l) sparse([ind ind ind],[ind-dn ind ind+dn],...
            [.5*(l*vw+l^2*vp.*vw) 1-(l^2)/2*(vp.*ve + vp.^2) ...
            0.5*(-l*ve+l^2*ve.^2)],xn,xn);
        
        Abd = sparse([1 xn],[1 xn],[1 1],xn,xn);
    elseif strcmp(method,'beamwarm')
        %create beam warming method
        A = @(ve,vw,vp,ind,dn,l) sparse([ind ind ind],[ind-2*dn ind-dn ind],...
            [-l/2*vw+.5*l^2*vw.*vp 2*l*vp-1/2*l^2*vp.*ve-1/2*l^2*vp.^2 ...
            1-3/2*l*ve+1/2*l^2*ve.^2],xn,xn);
        
        Abd = @(ve,vp,ind,dn,l) sparse([ind ind 1 xn],[ind-dn ind 1 xn],...
            [2*l*vp-1/2*l^2*vp.*ve-1/2*l^2*vp.^2 1-3/2*l*ve+1/2*l^2*ve.^2 1 1],xn,xn);
        
        
    else
        error('incorrect method specified')
    end
    
end