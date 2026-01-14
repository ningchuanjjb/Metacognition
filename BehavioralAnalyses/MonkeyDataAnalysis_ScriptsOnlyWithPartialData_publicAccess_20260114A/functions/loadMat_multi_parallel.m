function [data, MAT_file_name] = loadMat_multi_parallel(searchName, targetPATH, temp_if_parallel)
cd(targetPATH)

MAT_file=dir([searchName, '*']);
MAT_file_name = cell(1, max(size(MAT_file)));
MAT_file_load_cell = cell(1, max(size(MAT_file)));
rowHeadings = cell(1, max(size(MAT_file)));

for i=1:size(MAT_file,1)
    MAT_file_name(i) = {[MAT_file(i).folder,'\',MAT_file(i).name]};
    rowHeadings(i) = {sprintf('file%d', i)};
end
size_MAT_file = size(MAT_file,1);  
batchSize = 2;
batchNum = ceil(size_MAT_file/batchSize);


if temp_if_parallel == 0
    for i=1:size_MAT_file
        MAT_file_name(i) = {[MAT_file(i).folder,'\',MAT_file(i).name]};
        %MAT_file_load_cell(i) = {load(cell2mat(MAT_file_name(i)))};
        MAT_file_load_cell{i} = load(MAT_file_name{i});
        rowHeadings(i) = {sprintf('file%d', i)};
    end
elseif temp_if_parallel == 1
    for tempi=1:batchNum
        tempj = (tempi-1)*batchSize + 1;
        parfor i=tempj:min((tempj+batchSize),size_MAT_file)            
            MAT_file_load_cell{i} = load(MAT_file_name{i});
        end
    end
end
MAT_file_load = cell2struct(MAT_file_load_cell,rowHeadings,2);
data = MAT_file_load;
