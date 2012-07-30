function f_plot_matrix_std(matrix)
%% plot the array #vec, where y(k) is X and y(k+1)

[M N]=size(matrix);


for i=1:M
    str=sprintf('%s \t & %s \t & %s \t & %s \t & %s \t \\\\',num2str(matrix(i,1),5), num2str(matrix(i,2),5), num2str(matrix(i,3),5), num2str(matrix(i,4),5),num2str(matrix(i,5),5));
    disp(str);    
end

end

