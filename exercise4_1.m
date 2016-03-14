close all
clear
clc

rand_num = rand (4,8);

for i = 1: length(rand_num(:,1))    
    for j = 1:length(rand_num(1,:))
        
        if rand_num(i,j) <= 0.5
            
            switch i
                case 1
                    suff_row = 'st';                    
                case 2 
                    suff_row = 'nd';                    
                case 3
                    suff_row = 'rd';
                otherwise
                    suff_row = 'th';
            end
                 
            switch j
                case 1
                    suff_col = 'st';                    
                case 2
                    suff_col = 'nd';                    
                case 3
                    suff_col = 'rd';
                otherwise
                    suff_col = 'th'; 
            end              
            disp (['The ', num2str(i), suff_row, ' row and ', num2str(j), suff_col, ' column has a value of ',...
                num2str(rand_num(i,j)), ' and is not bigger than 0.5.'])
        end
    end
end

