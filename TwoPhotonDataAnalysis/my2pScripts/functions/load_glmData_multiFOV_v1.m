function [data, MAT_file_name] = load_glmData_multiFOV_v1(searchName, targetPATH)
% cd(targetPATH)

%MAT_file=dir([pwd,targetPATH,'\',searchName,'*']);
MAT_file=dir([targetPATH,'\','*',searchName,'*']);
MAT_file_name = cell(1, max(size(MAT_file)));
MAT_file_load_cell = cell(1, max(size(MAT_file)));
rowHeadings = cell(1, max(size(MAT_file)));
for i=1:size(MAT_file)
    %MAT_file_name(i) = {MAT_file(i).name};
    MAT_file_name(i) = {[MAT_file(i).folder,'\',MAT_file(i).name]};
    MAT_file_load_cell(i) = {load(cell2mat(MAT_file_name(i)))};
    rowHeadings(i) = {sprintf('file%d', i)};
end
MAT_file_load = cell2struct(MAT_file_load_cell,rowHeadings,2);
data = MAT_file_load;
