
IC_str = '_front';
stat_meth = '_autor';

% load(['CI' IC_str '_OLS.mat'])
load(['advection_art_data' IC_str '.mat'])

xnsize = [21,41,81,161,321,641,2*640+1];

num_meth_cell = cell(4,1);
num_meth_cell{1} = 'Upwind';
num_meth_cell{2} = 'Lax-Friedrich';
num_meth_cell{3} = 'Lax-Wendroff';
num_meth_cell{4} = 'Beam warming';


num_meth_short_cell = cell(4,1);
num_meth_short_cell{1} = 'UW';
num_meth_short_cell{2} = 'LF';
num_meth_short_cell{3} = 'LW';
num_meth_short_cell{4} = 'BW';

for i = 1:length(xd)
    xndata(i) = length(xd{i});
end
for i = 1:length(eta)
    eta_vec(i) = eta(i);
end


xnstr = zeros(1,8);
eta_str = zeros(1,8);


for m = 1:8
    

    xdi = ceil(m/length(eta_vec));
    sigmaj = mod(m,length(eta_vec));
    
    if sigmaj == 0
        sigmaj = length(eta_vec);
    end
    
    xnstr(m) = xndata(xdi);
    eta_str(m) = eta_vec(sigmaj);
    
end





for j = [1 4]

    figure('units','normalized','outerposition',[0 0 1 1])
    
    for i = 1:8

        subplot(2,4,i)
        hold on
        plot(log2(q0(1)),log2(q0(2)),'k*')



        for l = 1:7
            C = CI{l,i,j};
            if (~any(isnan(C(:))))&&(~isempty(C))
                rectangle('position',[log2(C(1,1)) log2(C(2,1)) log2(C(1,2))-log2(C(1,1)) log2(C(2,2))-log2(C(2,1))],...
                    'edgecolor',[(l-0)/7 0 1-(l-0)/7])%repmat(1-(l+1)/8,1,3))
            end
        end

        xlabel('$\log_2(\alpha)$','interpreter','latex')
        ylabel('$\log_2(\beta)$','interpreter','latex')

        title(strcat('95\% CI,  ',num_meth_short_cell{j},' method, $N$ = ',num2str(xnstr(i)),', $\eta^2 = $',...
                num2str(eta_str(i)^2)),'interpreter','latex')

        if strcmp(IC_str,'_front')
            if j == 1
                axis([-3 -1.5 -1.1 -.7])
            elseif j == 3
                axis([-2.5 .5 -1.7 -.7])
            elseif j == 4
                axis([-2.7 -1.5 -1.3 -.8])
            end
        elseif strcmp(IC_str,'_gauss')
            if j == 1
                axis([-3.5 -1.5 -1.4 -.75])
            elseif j == 2
                axis([-2 -1.5 -1.35 -1.26])
            elseif j ==3
%                 axis([-1.5 -2 -1.35 -1.26])
            elseif j == 4
                axis([-2.4 -1.65 -1.4 -1.15])
            end
        end
        
        
    end

% %     exportfig(gcf,['CI_h_plot' IC_str '_' num2str(j) stat_meth '.eps'],'fontsize',1.5,'color','rgb')
% %     saveas(gcf,['CI_h_plot' IC_str '_' num2str(j) stat_meth '.fig'])

    
end