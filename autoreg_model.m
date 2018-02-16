
m = 3;
num_meth = 1;

phihat1 = zeros(7,6);
phihat2 = zeros(7,6);

res_mod = cell(6,1);


for xni = 1:7

    res_mod{xni} = zeros(size(res{xni}));
    
    for tj = 2:6
        
        %get relevant residual values    
        a = res{xni}(tj,:);

%         %find analytical location of discontinuity
%         [g,sigma,sigma_inv] = advection_rate('root',q{xni+1,m,num_meth}(1),...
%             q{xni+1,m,num_meth}(2));
% 
%         %find indices before/after discontinuity
%         res_past_shock = a(xdata>sigma_inv(tdata(tj),.2));
%         res_before_shock = a(xdata<sigma_inv(tdata(tj),.2));
%         
%         ind_past_shock = find(xdata>sigma_inv(tdata(tj),.2));
%         ind_before_shock = find(xdata<=sigma_inv(tdata(tj),.2));
        
        a_max = max(abs(a));
        a_max_loc = xdata(abs(a)==a_max);
%         
        if a(abs(a)==a_max) >= 0 %max res value positive

            res_past_shock = a(xdata>=a_max_loc);
            res_before_shock = a(xdata<a_max_loc);

            ind_past_shock = find(xdata>=a_max_loc);
            ind_before_shock = find(xdata<a_max_loc);
            
        elseif a(abs(a)==a_max) < 0 %max res value negative
           
            
            res_past_shock = a(xdata>a_max_loc);
            res_before_shock = a(xdata<=a_max_loc);

            ind_past_shock = find(xdata>a_max_loc);
            ind_before_shock = find(xdata<=a_max_loc);
           
            
        end

        
%         phihat(xni,tj) = sum(x_past_shock(1:end-1).*x_past_shock(2:end))/sum(a(past_ind:end).^2);
        phihat1(xni,tj) = abs(phihat_estimate(res_past_shock));
        phihat2(xni,tj) = abs(phihat_estimate(res_before_shock));
        [B1,B2] = autoreg_mat(phihat1(xni,tj),phihat2(xni,tj),ind_past_shock,ind_before_shock,length(xdata));
        

        res_mod{xni}(tj,:) = (B1+B2)*a';

    %     figure
    %     subplot(2,1,1)
    %     plot(a,'.')
    %     axis([ind-1 50 -.1 2])
    % 
    %     subplot(2,1,2)
    %     plot(B*a','.')
    %     axis([ind-1 50 -.1 2])

    end

end

for i = 1:7

    figure('units','normalized','outerposition',[0 0 1 1])
    for j = 2:6
        subplot(3,2,j)

%         plot(res{i}(j,:),['b.'])
        hold on
        plot(xdata,res_mod{i}(j,:),[colors(j) '*'])
        plot([xdata(1) xdata(end)],[0 0],'k')
        axis([0 1 -1 1])
        
        xlabel('x')
        ylabel('A*(model - data)')
        title(['Autoregessive Residual, t = ' num2str(tdata(j))])
    end
    
        
        
    
end


exportfig(gcf,['mod_residual_' IC_str '_' num2str(i) '_noise_' num2str(eta_str(m)) '.eps'],'fontsize',2,'color','rgb')
saveas(gcf,['mod_residual_' IC_str '_' num2str(i) '_noise_' num2str(eta_str(m)) '.fig'])