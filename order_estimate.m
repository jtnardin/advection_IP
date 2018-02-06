%%%%%% order_estimate.m written 2-5-18 by JTN
%%% mean to find some order estimates for J(h), q(h)



load('advection_rates_gauss_IC.mat')
load('advection_art_data_gauss.mat')

xnsize = [21,41,81,161,321,641,2*640+1];

num_meth_cell = cell(4,1);
num_meth_cell{1} = 'upwind';
num_meth_cell{2} = 'Lax-Friedrich';
num_meth_cell{3} = 'Lax-Wendroff';
num_meth_cell{4} = 'Beam warming';

%%%%%% J computation

%estimate order using best-fit line
order_table = zeros(4,4);
poly_table = cell(4,4);

for i = 1:4
    for j = 1:4
        
        poly_table{i,j} = polyfit(log(1./xnsize)',log(squeeze(J(:,i,j))),1);
        order_table(i,j) = poly_table{i,j}(1);
             
    end
end


%estimate order using 2-point methods
R=abs(J(1:end-1,:,:)./J(2:end,:,:));

phat = log2(R);


%%%%%%now do the same for q

%get norm of ||q-q0||_2
q_norm = zeros(7,4,4);

for i = 1:7
    for j = 1:4
        for k = 1:4
            
            if J(i,j,k)>0
                q_norm(i,j,k) = norm(q{i,j,k}-q0);
            end
            
        end
    end
end

%%%%% write to latex table.

input.data = order_table;
input.dataFormat = {'%.3f'};
input.makeCompleteLatexDocument = 1;
text = latexTable(input);
q_text = ' ';
for i = 1:length(text)
    q_text = [q_text text{i}];
end
fileID = fopen('J_1.tex','w');
fprintf(fileID,'%s',q_text);
fclose(fileID)


input2.data = squeeze(mean(phat));
input2.dataFormat = {'%.3f'};
input2.makeCompleteLatexDocument = 1;
text2 = latexTable(input2);
q_text2 = ' ';
for i = 1:length(text2)
    q_text2 = [q_text2 text2{i}];
end
fileID2 = fopen('J_2.tex','w');
fprintf(fileID2,'%s',q_text2);
fclose(fileID2)

%estimate order using best-fit line
q_order_table = zeros(4,4);
q_poly_table = cell(4,4);

for i = 1:4
    for j = 1:4
        
        q_poly_table{i,j} = polyfit(log(1./xnsize)',log(squeeze(q_norm(:,i,j))),1);
        q_order_table(i,j) = q_poly_table{i,j}(1);
             
    end
end


%estimate order using 2-point methods
qR=abs(q_norm(1:end-1,:,:)./q_norm(2:end,:,:));

qphat = log2(qR);

