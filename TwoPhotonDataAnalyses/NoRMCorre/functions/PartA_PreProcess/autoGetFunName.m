function searchName_v = autoGetFunName(searchName, targetPATH)
%cd(targetPATH)

%MAT_file=dir([pwd,'\','functions','\',searchName,'*','_v','*','.m']);
MAT_file=dir([targetPATH,'\','\',searchName,'*','_v','*','.m']);
MAT_file_name = {MAT_file.name};
temp_fileNum = zeros(1, length(MAT_file_name));
for tempi=1:length(MAT_file_name)
    temp_name = MAT_file_name{tempi};
    temp_v = strfind(temp_name,'v');
    temp_fileNum(tempi) = str2num(temp_name(temp_v+1:end-2));%#ok<ST2NM>
end
[tempM,I] = max(temp_fileNum); %#ok<*ASGLU>
searchName_v = MAT_file_name{I}(1:end-2);