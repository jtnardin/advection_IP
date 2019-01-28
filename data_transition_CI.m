clear all; clc

IC_str = '_gauss';

load(['CI' IC_str '_OLS_3_1_newton.mat'])

CI1 = CI;

load(['CI' IC_str '_OLS_3_1_newton2.mat'])

CI2 = CI;

CI = cell(7,8,4);

for i = 1:7
    for j = 1:8
        for k = [1 4]
            
            CI{i,j,k} = CI1{i,j,k};
            
        end
    end
end



for i = 1:7
    for j = 1:8
        for k = 3
            
            CI{i,j,k} = CI2{i,j,k};
            
        end
    end
end

save(['CI' IC_str '_OLS.mat'],'CI')