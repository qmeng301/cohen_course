clear
close all
clc

row_n = 4;
col_n = 8;

rand_out_matrix  = rand_num_check_output(row_n,col_n);

fid = fopen('rand_matrix_output.txt','w');

% variable labels
variable_labels = {'row number';'column number';'result of test'};

for vari=1:length(variable_labels)
    fprintf(fid,'%s\t',variable_labels{vari});
    % the %s is for string; %g is for number.
end

% insert a new-line character
fprintf(fid,'\n');

for rowi=1:size(rand_out_matrix,1)% loop through rows (variables)
    
    % loop through columns (variables)
    %for columni=1:size(rand_out_matrix,2)
    %    fprintf(fid,'%g\t',rand_out_matrix(rowi,columni));
    %end
    %fprintf(fid,'\n'); % end-of-line 
    
    % You could also do this in one line:
     fprintf(fid,'%g\t%g\t%g\n',rand_out_matrix(rowi,1),rand_out_matrix(rowi,2),rand_out_matrix(rowi,3));
    
    fprintf('Finished writing line %g of %g\n',rowi,size(rand_out_matrix,1));
end

fclose(fid);
