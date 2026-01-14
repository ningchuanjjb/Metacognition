%% Initialization
% clear
close all

if_load = 0;

currentSessionIndex_AB = 2;



currentABSession_multi = string;

currentABSession_multi = [currentABSession_multi; '20230512A_and_20230513A'];%1
currentABSession_multi = [currentABSession_multi; '20230527A_and_20230528A'];%2
currentABSession_multi = [currentABSession_multi; '20230530A_and_20230531A'];%3

currentABSession_multi(1) = [];
num_FOV_AB = length(currentABSession_multi);


currentABSession = currentABSession_multi{currentSessionIndex_AB};

fprintf('currentABSession = %s.\n',currentABSession);



if if_load == 1
    %temp_load = load(['D:\twoPhotonData_motionCorrected\twoSessions\' '20230527A_and_20230528A']);
    temp_load = load(['D:\twoPhotonData_motionCorrected\twoSessions\' currentABSession]);    
    decodingDataSimplified_AB = temp_load.decodingDataSimplified_AB;
    clear temp_load
end

a = 1;
decodingDataSimplified_AB;









%% End