%optim_running written 2-3-18 by JTN to run optimization routine
%several times over

q_ols = cell(7,12,4);
J_ols = zeros(7,12,4);

IC_str = '_gauss';

tic
for i = 1:7
    for j = 1:12
        for k = 1:4
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



             [q_ols{i,j,k},J_ols(i,j,k)]= art_advec_fitting_f(i,j,num_meth,IC_str);
            
        end
    end
end
% end
toc

save(['advection_rates' IC_str '_IC_3_6.mat'],'q_ols','J_ols');%,...

