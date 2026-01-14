function data = loadMat_single(searchName, targetPATH)
cd(targetPATH)

MAT_file=dir([searchName, '*']);
% MAT_file_name = {MAT_file.name};
MAT_file_name = {[MAT_file.folder,'\',MAT_file.name]};
MAT_file_load_cell = {load(cell2mat(MAT_file_name))};
MAT_file_load = MAT_file_load_cell{:};

data = MAT_file_load;
