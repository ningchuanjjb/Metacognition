profile on
t0_cat2 = tic;
for tempi=1:1
    [filename_output,length_list_output,error_flag] = fun_concatenate(file_list,M1,'hdf5');
end
fprintf('\nConcatenate     ending, time = %.1f seconds.\n',toc(t0_cat2));
profile viewer
