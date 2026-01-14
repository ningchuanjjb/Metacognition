% Chuan's 32th script (20260114)
%This script: To preprocess edf files from eyelink and save them as mat files.
close all
clear all %#ok<CLALL>

% edfPATH = 'C:\ASDROOT\STUDY\edfdata\edf';
% matPATH = 'C:\ASDROOT\STUDY\edfdata\mat';

edfPATH = 'D:\edfdata\edf';
matPATH = 'D:\edfdata\mat';

cd(edfPATH)


% searchName = 'EYETEST1_';
% searchName = 'edf2022-07-12trainZelku_Step13_forRecord_niMarker_202005027C_27-2';
% searchName = 'Z_';
% searchName = 'ZTrain1_';
% searchName = 'D_';
% searchName = 'D2_';
% searchName = ' DTrain2_';
% searchName = 'test_';
% searchName = '113TEST1';

% searchName = 'D2p_';
searchName = 'Z2p_';


MAT_file=dir([edfPATH, '\', searchName, '*', 'edf']);


file_size = max(size(MAT_file));
for file_index=1:file_size
    file_name = MAT_file(file_index).name(1:end-4);
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
    
    
    save([matPATH, '\', file_name '.mat'], 'edf0_saved');
    a = 1;
end

cd 'C:\ASDROOT\STUDY\MonkeyDataAnalysis'