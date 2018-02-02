

% j = 4;
% 
% x = linspace(0,1,10*2^(j-1)+1);
% 
% figure('units','normalize','outerposition',[0 0 1 1])
% 
% for i = 1:6
%     subplot(3,2,i)
%     
%        
%     plot(x,model_sims{j}(:,i),'-')
%     
%     hold on
%     
%     plot(x,soln(tdata(i),x),'--')
%     
%     if i == 1
%         legend('Numerical','True','location','northeast')
%     end
%     
% end

% exportfig(gcf,['advection_compare_' num2str(j) '.eps'],'color','rgb','fontsize',1.5);
% saveas(gcf,['advection_compare_' num2str(j) '.fig']);



% j = 5;

x = linspace(0,1,size(model_sims,1));

figure('units','normalize','outerposition',[0 0 1 1])

for i = 1:6
    subplot(3,2,i)
    hold on
    
    for j = 2:2:6
       
        plot(x,model_sims(:,i,j),'-')
    
    end
    
    plot(x,soln(tdata(i),x),'--')
    
    if i == 1
        h=legend(['Numerical, dx = ' num2str(1/(xnsize(2)-1))],['Numerical, dx = ' num2str(1/(xnsize(4)-1))],...
            ['Numerical, dx = ' num2str(1/(xnsize(6)-1))],'True Soln','location','northeast');
        set(h,'interpreter','latex')
    end
    
end

% exportfig(gcf,['advection_compare_' num2str(j) '_' num_method '.eps'],'color','rgb','fontsize',1.5);
% saveas(gcf,['advection_compare_' num2str(j) '_' num_method '.fig']);