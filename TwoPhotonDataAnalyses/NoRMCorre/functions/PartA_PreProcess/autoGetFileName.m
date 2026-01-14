function auto_file_name = autoGetFileName(searchName, targetPATH, if_max0_min1)
%cd(targetPATH)

if ~exist('if_max0_min1','var')
    if_max0_min1 = 0;
end

MAT_file=dir([targetPATH, '\',searchName, '*']);
file_num = length(MAT_file);
temp_date = zeros(1, file_num);
temp_time = zeros(1, file_num);
temp_date_and_time = zeros(1, file_num);
for tempi=1:file_num
    temp_name = MAT_file(tempi).name;
    %temp_date(tempi) = str2double(temp_name(15:end-10));
    %temp_time(tempi) = str2double(temp_name(end-8:end-3));
    %temp_date_and_time(tempi) = str2double([temp_name(15:end-10) temp_name(end-8:end-3)]);
    temp_date(tempi) = str2double(temp_name(end-17:end-10));
    temp_time(tempi) = str2double(temp_name(end-8:end-3));
    temp_date_and_time(tempi) = str2double([temp_name(end-17:end-10) temp_name(end-8:end-3)]);    
end
if if_max0_min1 == 0
    [M,I]=max(temp_date_and_time); %#ok<*ASGLU>
elseif if_max0_min1 == 1
    [M,I]=min(temp_date_and_time); %#ok<*ASGLU>
end
target_id = I;

auto_file_name = [MAT_file(target_id).folder,'\',MAT_file(target_id).name];
