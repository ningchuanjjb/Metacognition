function auto_file_name = autoGetFileName_general(searchName, targetPATH, if_max0_min1)
%cd(targetPATH)

if ~exist('if_max0_min1','var')
    if_max0_min1 = 0;
end

MAT_file=dir([targetPATH, '\',searchName, '*']);
file_num = length(MAT_file);

temp_name = cell(1, file_num);
for tempi=1:file_num
    temp_name{tempi} = MAT_file(tempi).name;  
end
temp_name_str = string(temp_name);
[B,sortIndex] = sort(temp_name_str,2,'descend'); %#ok<*ASGLU>
if isempty(sortIndex) == false
    if if_max0_min1 == 0
        I = sortIndex(1);
    elseif if_max0_min1 == 1
        I = sortIndex(end);
    end
    target_id = I;
    
    auto_file_name = [MAT_file(target_id).folder,'\',MAT_file(target_id).name];
else 
    auto_file_name = '';
end
