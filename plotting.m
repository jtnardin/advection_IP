figure

j = 3;

[g,sigma,sigma_inv] = advection_rate('root',q{j,m,num_meth}(1),...
    q{j,m,num_meth}(2));

for i = 1:6
    subplot(3,2,i)
    
    for j = 4
        plot(xdata,res{j}(i,:),'*')
        hold on
        plot(xdata,res_mod{j}(i,:),'s')
    end
%     plot(xdata,model_sims{j}(i,:))
%     plot(repmat(sigma_inv(tdata(i),.1),2,1),[-2 6],'k--')
%     plot(repmat(sigma_inv(tdata(i),.3),2,1),[-2 6],'k--')
%     axis([0 1 -2 6])
    plot([0 1],[0 0],'k--')
    axis([0 1 -4 4])
end