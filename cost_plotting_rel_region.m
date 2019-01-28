clear all; clc

rel_range_gauss = cell(8,4);
rel_range_front = cell(8,4);

% default
for i = 1:8
    for j = 1:4
        rel_range_gauss{i,j} = 1:7;
        rel_range_front{i,j} = 1:7;
    end
end
        
%exceptions
% rel_range_front{1} = 1:7;
% rel_range_front{2} = 1:6;
% rel_range_front{3} = 1:7;
% rel_range_front{4} = 1:4;
% 


rel_range_gauss{1,1} = 1:5;
rel_range_gauss{1,3} = 1:4;
rel_range_gauss{1,4} = 1:5;

rel_range_gauss{2,1} = 1:5;
rel_range_gauss{2,3} = 1:3;
rel_range_gauss{2,4} = 1:4;

rel_range_gauss{3,1} = 1:4;
rel_range_gauss{3,3} = 1:3;
rel_range_gauss{3,4} = 1:3;


rel_range_gauss{4,1} = 2:7;
rel_range_gauss{4,3} = 2:7;
rel_range_gauss{4,4} = 2:7;


rel_range_gauss{5,1} = 1:6;
rel_range_gauss{5,3} = 1:4;
rel_range_gauss{5,4} = 1:6;


rel_range_gauss{6,1} = 1:6;
rel_range_gauss{6,3} = 1:4;
rel_range_gauss{6,4} = 1:4;

rel_range_gauss{7,1} = 1:4;
rel_range_gauss{7,3} = 1:3;
rel_range_gauss{7,4} = 1:3;


rel_range_gauss{8,1} = 2:7;
rel_range_gauss{8,3} = 2:7;
rel_range_gauss{8,4} = 2:7;


rel_range_front{1,3} = 2:7;

rel_range_front{2,1} = 1:6;
rel_range_front{2,3} = 2:7;
rel_range_front{2,4} = 1:3;

rel_range_front{3,1} = 1:5;
rel_range_front{3,3} = 2:4;
rel_range_front{3,4} = 1:3;



rel_range_front{5,1} = 2:7;
rel_range_front{5,3} = 1:6;
rel_range_front{5,4} = [1:5 7];


rel_range_front{6,1} = 2:7;
rel_range_front{6,3} = 1:7;
rel_range_front{6,4} = 1:5;

rel_range_front{7,1} = 2:6;
rel_range_front{7,3} = 1:4;
rel_range_front{7,4} = 2:5;


save('rel_range_cost_sims.mat');

