%optim_running written 2-3-18 by JTN to run optimization routine
%several times over

q_ols = cell(7,4,4);
J_ols = zeros(7,4,4);

q_autoreg = cell(7,4,4);
J_autoreg = zeros(7,4,4);

phi1 = cell(7,4,4);
phi2 = cell(7,4,4);

IC_str = '_front';

tic
for i = 7%1:7
    for j = 4%1:4
        for k = 1%1:4
            [i,j,k]
            
            if k == 1
                num_meth = 'upwind';
            elseif k == 2
                num_meth = 'laxfried';
            elseif k == 3
                num_meth = 'laxwend';
            elseif k == 4
                num_meth = 'beamwarm';
            end

%             [q_ols{i,j,k},J_ols(i,j,k)]= art_advec_fitting_f(i,j,num_meth,IC_str);


            [q_ols{i,j,k},J_ols(i,j,k),q_autoreg{i,j,k},J_autoreg(i,j,k)...
                ,phi1{i,j,k},phi2{i,j,k}]...
                = autoreg_art_advec_fitting_f(i,j,num_meth,IC_str);

            
        end
    end
end
toc

save(['advection_rates_autoreg_IC' IC_str '.mat'],'q_ols','J_ols',...
    'q_autoreg','J_autoreg','phi1','phi2')