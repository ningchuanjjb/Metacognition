function auto_file_name = autoGetMarkerFileName(searchName, targetPATH)
%cd(targetPATH)

MAT_file=dir([targetPATH, '\',searchName, '*']);
file_num = length(MAT_file);
temp_id = zeros(1, file_num);
for tempi=1:file_num
    temp_name = MAT_file(tempi).name;
    temp_id(tempi) = str2double(temp_name(22:24));
    
end
[M,I]=max(temp_id); %#ok<*ASGLU>

target_id = I;

auto_file_name = [MAT_file(target_id).folder,'\',MAT_file(target_id).name];
