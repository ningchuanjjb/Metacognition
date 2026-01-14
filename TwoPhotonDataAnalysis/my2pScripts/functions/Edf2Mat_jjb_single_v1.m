function file_name = Edf2Mat_jjb_single_v1(searchName,edfPATH_data,matPATH_save,originPATH)

% edfPATH_data = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\EDF\edfdata\edf';
% matPATH_save = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\EDF\edfdata\mat';

cd(edfPATH_data)


% searchName = '113TEST1';

MAT_file=dir([edfPATH_data, '\', searchName, '*', 'edf']);


file_size = max(size(MAT_file));
for file_index=1:file_size
    file_name = MAT_file(file_index).name(1:end-4);
    
    if exist([file_name,'.mat'],'file') ~= 0
        fprintf('Target eye movement mat file has already existed.\n');
        continue
    end    
    
    fprintf("Now is processing file %d!\n", file_index);
    fprintf("File name is '%s.edf'!\n", file_name);
    
    edf0 = Edf2Mat([file_name '.edf']);
    % plot(edf0);
    
    edf0_fieldnames = fieldnames(edf0);
    edf0_Samples_fieldnames = fieldnames(edf0.Samples);
    
    edf0_struct = struct;
    for tempi=1:size(edf0_fieldnames,1)
        if strcmp(edf0_fieldnames{tempi},'fails') == 0
            edf0_struct = setfield(edf0_struct, edf0_fieldnames{tempi}, getfield(edf0, edf0_fieldnames{tempi}));
        end
    end
    clear edf0
    
    edf0_struct;
    edf0_saved = struct;
    
    edf0_saved.gx = edf0_struct.Samples.gx;
    edf0_saved.gy = edf0_struct.Samples.gy;
    edf0_saved.pupilSize = edf0_struct.Samples.pupilSize;
    edf0_saved.time = edf0_struct.Samples.time;
    edf0_saved.Messages = edf0_struct.Events.Messages;
    edf0_saved.fixationEvent = edf0_struct.Events.Efix;
    clear edf0_struct
    
    
    save([matPATH_save, '\', file_name '.mat'], 'edf0_saved');
    
end

% cd 'C:\ASDROOT\STUDY\MonkeyDataAnalysis'
cd(originPATH)