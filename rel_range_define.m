clear all; clc

rel_range_gauss = cell(4,1);
rel_range_front = cell(4,1);

rel_range_front{1} = 1:7;
rel_range_front{2} = 1:6;
rel_range_front{3} = 1:7;
rel_range_front{4} = 1:4;


rel_range_gauss{1} = 1:7;
rel_range_gauss{2} = 4:7;
rel_range_gauss{3} = 1:4;
rel_range_gauss{4} = 2:5;

save('rel_range_sims.mat');

