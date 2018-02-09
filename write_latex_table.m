function write_latex_table(name,table)

    
    input.data = table;
    input.dataFormat = {'%.3f'};
    input.makeCompleteLatexDocument = 1;
    text = latexTable(input);
    q_text = ' ';
    for i = 1:length(text)
        q_text = [q_text text{i}];
    end
    fileID = fopen(name,'w');
    fprintf(fileID,'%s',q_text);
    fclose(fileID)




end