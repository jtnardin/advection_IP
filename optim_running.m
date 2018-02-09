%optim_running written 2-3-18 by JTN to run optimization routine
%several times over

q = cell(7,4,4);
J = zeros(7,4,4);

IC_str = '_front';

tic
for i = 1:4
    for j = 1:4
        for k = 1
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


            [q{i,j,k},J(i,j,k)] = art_advec_fitting_f(i,j,num_meth,IC_str);

            
        end
    end
end
toc

save(['advection_rates_IC' IC_str '.mat'],'J','q')