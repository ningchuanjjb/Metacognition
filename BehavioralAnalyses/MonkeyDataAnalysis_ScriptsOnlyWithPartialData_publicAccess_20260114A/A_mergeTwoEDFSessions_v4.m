% Please mannual load session a and session b ( eye link data with .mat)
% This script is to merge a and b

% a = edf0_saved;
% b = edf0_saved;
a;
b;
%a_trial_count = trial_para.trial_count;
%b_trial_count = trial_para.trial_count;

c = struct;
c.gx = [a.gx; b.gx];
c.gy = [a.gy; b.gy];
c.pupilSize = [a.pupilSize; b.pupilSize];
c.time = [a.time; b.time];

c.Messages = a.Messages;
c.Messages.time = [c.Messages.time b.Messages.time];
c.Messages.info = [c.Messages.info b.Messages.info];

c.fixationEvent = a.fixationEvent;
c.fixationEvent.eye = [c.fixationEvent.eye b.fixationEvent.eye];
c.fixationEvent.start = [c.fixationEvent.start b.fixationEvent.start];
c.fixationEvent.end = [c.fixationEvent.end b.fixationEvent.end];
c.fixationEvent.duration = [c.fixationEvent.duration b.fixationEvent.duration];
c.fixationEvent.posX = [c.fixationEvent.posX b.fixationEvent.posX];
c.fixationEvent.posY = [c.fixationEvent.posY b.fixationEvent.posY];
c.fixationEvent.pupilSize = [c.fixationEvent.pupilSize b.fixationEvent.pupilSize];

edf0_saved = c;

edf0_saved.Messages.info;
last_marker_count = 0;
last_marker_count_index = 0;
flag_isSecondFile = 0;
delete_index = struct;
for tempi=1:length(edf0_saved.Messages.info)
    edf0_saved.Messages.info{tempi};
    k = strfind(edf0_saved.Messages.info{tempi}, 'TRIALID ');% if find, k should be 1
    if k == 1
        marker_count = str2double(edf0_saved.Messages.info{tempi}(length('TRIALID ')+1:end));
        if flag_isSecondFile == 1
            marker_count = marker_count + a_trial_count;
            edf0_saved.Messages.info{tempi}(length('TRIALID ')+1:end) = [];
            edf0_saved.Messages.info{tempi} = [edf0_saved.Messages.info{tempi} num2str(marker_count)];            
        end
        
        marker_count_index = tempi;
        if marker_count < last_marker_count
            marker_count_index;
            last_marker_count_index;
            delete_index.marker_count_index = marker_count_index;
            delete_index.last_marker_count_index = last_marker_count_index;
            
            flag_isSecondFile = 1;
            marker_count = marker_count + a_trial_count;
            edf0_saved.Messages.info{tempi}(length('TRIALID ')+1:end) = [];
            edf0_saved.Messages.info{tempi} = [edf0_saved.Messages.info{tempi} num2str(marker_count)];                        
        end
        last_marker_count = marker_count;
        last_marker_count_index = tempi;
    end
end
edf0_saved.Messages.time(...
    delete_index.last_marker_count_index:(delete_index.marker_count_index-1)...
    ) = [];
edf0_saved.Messages.info(...
    delete_index.last_marker_count_index:(delete_index.marker_count_index-1)...
    ) = [];

if false
    %file_name = 'D_edf2022-06-29trainDing_Step12_fixation_202005027C_forBehavior_6-merge'
    file_name = 'D_edf2022-06-30trainDing_Step12_fixation_202005027C_forBehavior_6-merge'
    
    matPATH = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\EDF\edfdata\mat'; %#ok<*UNRCH>               
    save([matPATH, '\', file_name '.mat'], 'edf0_saved');
end