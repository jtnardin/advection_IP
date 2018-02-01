function A = aMatrixupwind(xn,method)


    if strcmp(method,'upwind')
        %create upwind matrix
        A = @(ve,vw,ind,dn,l) sparse([ind ind],[ind ind-dn],[1-l*ve l*vw],xn,xn);
    elseif strcmp(method,'laxfried')
        %create lax-friedrich method
        A = @(ve,vw,ind,dn,l) sparse([ind ind],[ind-dn ind+dn],...
            [1+l*vw 1-l*ve]/2,xn,xn);
    elseif strcmp(method,'laxwend')
        %create lax-wendroff method
        A = @(ve,vw,vp,ind,dn,l) sparse([ind ind ind],[ind-dn ind ind+dn],...
            [.5*(l*vw+l^2*vp.*vw) 1-(l^2)/2*(vp.*ve + vp.^2) ...
            0.5*(-l*ve+l^2*ve.^2)],xn,xn);
    %             [0.5*l*(1+l)*vw 1-vp*(l^2) 0.5*l*(-1+l)*ve],xn,xn);
    else
        error('incorrect method specified')
    end
    
end