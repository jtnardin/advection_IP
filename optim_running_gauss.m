%optim_running written 2-3-18 by JTN to run optimization routine
%several times over

IC_str = '_gauss';
data_str = ['advection_art_data' IC_str '_all_3_26.mat'];
model_str = 'root';

load(data_str)

load(['advection_rates' IC_str '_IC_all_3_26.mat']);%,...

tic
for i = 1:7
    for j = 1:numel(data)
        for k = 4
            [i,j,k]
            
            if isempty(q_ols{i,j,k})


                if k == 1
                    num_meth = 'upwind';
                elseif k == 2
                    num_meth = 'laxwend';
                elseif k == 3
                    num_meth = 'beamwarm';
                elseif k == 4
                    num_meth = 'upwindfl';
                end



                 [q_ols{i,j,k},J_ols(i,j,k)]= art_advec_fitting_f(i,j,num_meth,IC_str,model_str,data_str);
                 save(['advection_rates' IC_str '_IC_all_3_26.mat'],'q_ols','J_ols')
            end
        end
    end
end
% end
toc

save(['advection_rates' IC_str '_IC_all_3_26.mat'],'q_ols','J_ols');%,...

