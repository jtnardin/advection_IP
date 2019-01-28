IC_str = '_gauss';

load(['CI' IC_str '_OLS_2.mat'])

CI_new = CI;

load(['CI' IC_str '_OLS.mat'])

for i = 1:7
    for j = 1:8
        for k = 1:4
            
            if isempty(CI{i,j,k})
                CI{i,j,k} = CI_new{i,j,k};
            end
        end
    end
end

save(['CI' IC_str '_OLS.mat'],'CI')