clear all; clc

load('advection_rates_2_3.mat')
load('advection_art_data.mat')


num_meth_cell = cell(4,1);
num_meth_cell{1} = 'upwind';
num_meth_cell{2} = 'Lax-Friedrich';
num_meth_cell{3} = 'Lax-Wendroff';
num_meth_cell{4} = 'Beam warming';


xndata = [10, 50];
eta = [0 0.1];

xnstr = zeros(1,4);
eta_str = zeros(1,4);

for m = 1:4
    
    xdi = ceil(m/2);
    sigmaj = mod(m,2);
    
    if sigmaj == 0
        sigmaj = 2;
    end
    
    xnstr(m) = xndata(xdi);
    eta_str(m) = eta(sigmaj);
    
end


figure
for i = 1:4
    
    subplot(2,2,i)
    plot(J(:,:,i))
    
    if i == 1
        
        h=legend(strcat('$x_n$ = ',num2str(xnstr'),', $\sigma^2 = $',...
            num2str(eta_str')),'location','northeast');
        
        set(h,'interpreter','latex');
        
    end
    
    title(['J, ' num_meth_cell{i}])
    
end

qnorm = zeros(7,4,4);

for i = 1:7
    for j = 1:4
        for k = 1:4
            
            if J(i,j,k)>0
                q_norm(i,j,k) = norm(q{i,j,k}-q0);
            end
            
        end
    end
end



figure
for i = 1:4
    
    subplot(2,2,i)
    plot(q_norm(:,:,i))
    
    if i == 1
        
        h=legend(strcat('$x_n$ = ',num2str(xnstr'),', $\sigma^2 = $',...
            num2str(eta_str')),'location','northeast');
        
        set(h,'interpreter','latex');
        
    end
    
    axis([0 7 0 2])
    
    title(['$\| \hat{q} - q_0 \|_2$, ' num_meth_cell{i}],'interpreter','latex')
    
end
