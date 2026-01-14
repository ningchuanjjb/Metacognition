% Chuan's 31th script (20251214)
% This script: To prepare for pseudo-population analyses where data from all sessions were combined
%% Initialization
% clear
close all

if_compute = 1;
if_save = 0;

% if_monkey_D0_Z1;


if_resample_seqBased0_locBased1 = 1;


resampleTrialCount_min = 7;
resampleTrialCount_length = [14 7 7];%Only for seqBased mode, [12 6 6]
resampleIterCount = 32;

if if_resample_seqBased0_locBased1 == 1
    resampleTrialCount_min = 1;
end

targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)


currentABSession_multi = string;

if if_monkey_D0_Z1 == 0
    % currentABSession_multi = [currentABSession_multi; '20230426A_and_20230427A'];%1 few
    currentABSession_multi = [currentABSession_multi; '20230502A_and_20230503A'];%2
    currentABSession_multi = [currentABSession_multi; '20230504A_and_20230508A_and_20230509A'];%3
    currentABSession_multi = [currentABSession_multi; '20230510A_and_20230510B_and_20230511A'];%4
    % currentABSession_multi = [currentABSession_multi; '20230512A_and_20230513A'];%5 few
    currentABSession_multi = [currentABSession_multi; '20230515A_and_20230516A'];%6
    currentABSession_multi = [currentABSession_multi; '20230525A_and_20230526A'];%7
    currentABSession_multi = [currentABSession_multi; '20230527A_and_20230528A'];%8
    currentABSession_multi = [currentABSession_multi; '20230530A_and_20230531A'];%9
    currentABSession_multi = [currentABSession_multi; '20230601A_and_20230602A'];%10
    currentABSession_multi = [currentABSession_multi; '20230604A_and_20230605A'];%11
    currentABSession_multi = [currentABSession_multi; '20230614A_and_20230615A'];%12

elseif if_monkey_D0_Z1 == 1
        currentABSession_multi = [currentABSession_multi; '20240122A_and_20240123A'];%1
%         currentABSession_multi = [currentABSession_multi; '20240124A_and_20240126A'];%2, few        
        currentABSession_multi = [currentABSession_multi; '20240129A_and_20240131A'];%3        
%         currentABSession_multi = [currentABSession_multi; '20240202A_and_20240203A'];%4, few           
%         currentABSession_multi = [currentABSession_multi; '20240207A_and_20240208A'];%5, few           
%         currentABSession_multi = [currentABSession_multi; '20240210A_and_20240211A'];%6, few           
        currentABSession_multi = [currentABSession_multi; '20240216A_and_20240218A'];%7        
        currentABSession_multi = [currentABSession_multi; '20240220A_and_20240221A'];%8        
        currentABSession_multi = [currentABSession_multi; '20240226A_and_20240227A'];%9        
        currentABSession_multi = [currentABSession_multi; '20240229A_and_20240301A'];%10        
        currentABSession_multi = [currentABSession_multi; '20240304A_and_20240305A'];%11        
        currentABSession_multi = [currentABSession_multi; '20240307A_and_20240308A'];%12        
%         currentABSession_multi = [currentABSession_multi; '20240312A_and_20240315A'];%13, few           
        currentABSession_multi = [currentABSession_multi; '20240319A_and_20240320A'];%14        
        currentABSession_multi = [currentABSession_multi; '20240322A_and_20240323A'];%15        
%         currentABSession_multi = [currentABSession_multi; '20240329A_and_20240330A'];%16, few           
        %currentABSession_multi = [currentABSession_multi; '20240402A_and_20240403A'];%17, bad      
        %currentABSession_multi = [currentABSession_multi; '20240410A_and_20240411A'];%18, bad
        
end


currentABSession_multi(1) = [];
num_FOV_AB = length(currentABSession_multi);


output_path = 'D:\twoPhotonData_motionCorrected\twoSessions';





numFrames = 6;
pointKindsNum = 4;
target_seqSet = get_target_seqSet_v1(numFrames,pointKindsNum);
numSeq = zeros(1,pointKindsNum);
for tempi=1:pointKindsNum
    numSeq(tempi) = length(target_seqSet{tempi});
end


%% Get resampleTrialCount_seq
if if_resample_seqBased0_locBased1 == 0
    resampleTrialCount_seq = zeros(1,sum(numSeq(1:3)));
    for tempi=1:3
        temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
        resampleTrialCount_seq(temp_range) = resampleTrialCount_length(tempi);
    end
    
elseif if_resample_seqBased0_locBased1 == 1
    if if_compute == 1
        seqCount_valid_allSessions = cell(1,num_FOV_AB);
        
        for currentSessionIndex=1:num_FOV_AB
            % for currentSessionIndex=1:2
            currentABSession = currentABSession_multi{currentSessionIndex};
            fprintf('currentABSession = %s.\n',currentABSession);
            
            temp_load = load([output_path '\' currentABSession]);
            temp_decodingDataSimplified = temp_load.decodingDataSimplified_AB;
            
            seqIndex = temp_decodingDataSimplified.seqIndex;
            trialIndex_bool_memoryCorrect = temp_decodingDataSimplified.trialIndex_bool_memoryCorrect;
            target_seqSet_inOne = temp_decodingDataSimplified.target_seqSet_inOne;
            
            temp_trialIndex_bool = trialIndex_bool_memoryCorrect; % use memory correct trial
            seqIndex_valid = seqIndex(temp_trialIndex_bool);
            
            temp_seqCount = zeros(1,sum(numSeq(1:3)));
            for tempi=1:sum(numSeq(1:3))
                temp_seqCount(tempi) = sum(seqIndex_valid==tempi,'all');
            end
            seqCount_valid_allSessions{currentSessionIndex} = temp_seqCount;
        end
        
    end
    
    
    seqCount_valid_allSessions_collapsed = zeros(length(seqCount_valid_allSessions{1}),num_FOV_AB);
    for tempi=1:num_FOV_AB
        seqCount_valid_allSessions_collapsed(:,tempi) = seqCount_valid_allSessions{tempi};
    end
    
    seqCount_valid_allSessions_collapsed_min = min(seqCount_valid_allSessions_collapsed,[],2);
    
    a = 1;
    
    seqIndex_valid_fromAllSessionsMin = zeros(sum(seqCount_valid_allSessions_collapsed_min),1);
    for tempi=1:length(seqCount_valid_allSessions_collapsed_min)
        temp_range = (sum(seqCount_valid_allSessions_collapsed_min(1:tempi-1))+1):sum(seqCount_valid_allSessions_collapsed_min(1:tempi));
        
        seqIndex_valid_fromAllSessionsMin(temp_range) = tempi;
    end
    
    boolIndex_location_seq = false(numFrames,length(target_seqSet_inOne));
    for tempi=1:length(target_seqSet_inOne)
        currentSequence = target_seqSet_inOne{tempi};
        boolIndex_location_seq(currentSequence,tempi) = true;
    end
    boolIndex_location_seq_T = boolIndex_location_seq';
    
    
    threshold_locRatio = 0.055;%0.06-->0.05-->0.045
    threshold_locStd = 0.026;%0.0278-->0.025-->0.026
    
    fun_locBasedResample_seqCount_Name_v = autoGetFunName_myScripts('fun_locBasedResample_seqCount', [targetPATH '\functions']);
    fun_locBasedResample_seqCount = str2func(fun_locBasedResample_seqCount_Name_v);
    
    options_locBased.threshold_locRatio = threshold_locRatio;
    options_locBased.threshold_locStd = threshold_locStd;
    options_locBased.numSeq = numSeq;
    options_locBased.numFrames = numFrames;
    options_locBased.seqIndex_valid = seqIndex_valid_fromAllSessionsMin;
    options_locBased.boolIndex_location_seq_T = boolIndex_location_seq_T;
    
    resampleTrialCount_seq = fun_locBasedResample_seqCount(options_locBased);    
end
resampleTrialCount_seq;

a1 = sum(resampleTrialCount_seq>0);

a = 1;

if if_compute == 1
    temp_F_dff_decisionBin_resample_allSessions = cell(1,num_FOV_AB);
    seqIndex_valid_resample_allSessions = cell(1,num_FOV_AB);
    roiNum_allSessions = zeros(1,num_FOV_AB);    
    seqCount_valid_allSessions = cell(1,num_FOV_AB);
    
    for currentSessionIndex=1:num_FOV_AB
        % for currentSessionIndex=1:2
        currentABSession = currentABSession_multi{currentSessionIndex};
        fprintf('currentABSession = %s.\n',currentABSession);
        
        temp_load = load([output_path '\' currentABSession]);
        temp_decodingDataSimplified = temp_load.decodingDataSimplified_AB;
        
        seqIndex = temp_decodingDataSimplified.seqIndex;
        trialIndex_bool_memoryCorrect = temp_decodingDataSimplified.trialIndex_bool_memoryCorrect;
        decisionPeriodA_interval = temp_decodingDataSimplified.decisionPeriodA_interval;
        F_dff_decisionPeriodA = temp_decodingDataSimplified.F_dff_decisionPeriodA;
        roiNum_allSessions(currentSessionIndex) = size(temp_decodingDataSimplified.F_dff_decisionPeriodA,1);
        
        
        temp_trialIndex_bool = trialIndex_bool_memoryCorrect; % use memory correct trial
        F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
        F_dff_decisionBin1 = double(F_dff_decisionBin1);
        F_dff_decisionBin1 = F_dff_decisionBin1 + eps;
        
        
        temp_F_dff_decisionBin = F_dff_decisionBin1(:,temp_trialIndex_bool);
        seqIndex_valid = seqIndex(temp_trialIndex_bool);
        
        
        temp_seqCount = zeros(1,sum(numSeq(1:3)));
        for tempi=1:sum(numSeq(1:3))
            temp_seqCount(tempi) = sum(seqIndex_valid==tempi,'all');
        end
        
        %temp_F_dff_decisionBin_resample = zeros(resampleIterCount,size(temp_F_dff_decisionBin,1),sum(numSeq(1:3))*resampleTrialCount);
        %seqIndex_valid_resample = zeros(resampleIterCount,sum(numSeq(1:3))*resampleTrialCount);
        for tempIter=1:resampleIterCount
            
            %tempIndex_resample = zeros(resampleTrialCount,sum(numSeq(1:3)));
            tempIndex_resample_cell = cell(1,sum(numSeq(1:3)));            
            for tempi=1:sum(numSeq(1:3))
                %temp_resampleTrialCount = resampleTrialCount;
                temp_resampleTrialCount = resampleTrialCount_seq(tempi);
                
                tempIndex_resample_cell{tempi} = zeros(1,temp_resampleTrialCount);
                tempIndex = find(seqIndex_valid==tempi);
                temp_length = length(tempIndex);
                if temp_length == 0
                    continue
                end
                if temp_length >= temp_resampleTrialCount
                    tempI = sort(randperm(temp_length,temp_resampleTrialCount));
                elseif temp_length < temp_resampleTrialCount
                    temp_resampleTrialCount2 = temp_resampleTrialCount;
                    temp_intLoopCount = floor(temp_resampleTrialCount2/temp_length);
                    temp_resampleTrialCount2 = temp_resampleTrialCount2-temp_length*temp_intLoopCount;
                    tempI1 = [];
                    for temptempi=1:temp_intLoopCount
                        tempI1 = [tempI1 (1:temp_length)]; %#ok<*AGROW>
                    end
                    tempI2 = sort(randperm(temp_length,temp_resampleTrialCount2));
                    
                    tempI = [tempI1 tempI2]; %#ok<*NASGU>
                end
                tempI = tempI(randperm(length(tempI))); % shuffle the resampled trial index
                
                %tempIndex_resample(:,tempi) = tempIndex(tempI);
                tempIndex_resample_cell{tempi} = tempIndex(tempI);
            end
            
            tempIndex_resample_cell;            
            %tempIndex_resample;
            %tempIndex_resample_1d = reshape(tempIndex_resample,1,[]);
            tempIndex_resample_1d = [];
            for tempi=1:length(tempIndex_resample_cell)
                tempIndex_resample_1d = [tempIndex_resample_1d,tempIndex_resample_cell{tempi}];                
            end
            a = 1;
            
            temptemp_isvalid = tempIndex_resample_1d > 0;
            
            % -1 means no trial in this sequence
            temptemp_F_dff_decisionBin_resample = -1*ones(size(temp_F_dff_decisionBin,1),length(tempIndex_resample_1d));
            tempseqIndex_valid_resample = -1*ones(1,length(tempIndex_resample_1d));
            
            temptemp_F_dff_decisionBin_resample(:,temptemp_isvalid) = temp_F_dff_decisionBin(:,tempIndex_resample_1d(temptemp_isvalid));
            tempseqIndex_valid_resample(temptemp_isvalid) = seqIndex_valid(tempIndex_resample_1d(temptemp_isvalid));
            
            if tempIter == 1
                %temp_F_dff_decisionBin_resample = zeros(resampleIterCount,size(temp_F_dff_decisionBin,1),sum(numSeq(1:3))*temp_resampleTrialCount);
                %seqIndex_valid_resample = zeros(resampleIterCount,sum(numSeq(1:3))*temp_resampleTrialCount);
                temp_F_dff_decisionBin_resample = zeros(resampleIterCount,size(temp_F_dff_decisionBin,1),length(tempIndex_resample_1d));
                seqIndex_valid_resample = zeros(resampleIterCount,length(tempIndex_resample_1d));                
            end
            temp_F_dff_decisionBin_resample(tempIter,:,:) = temptemp_F_dff_decisionBin_resample;
            seqIndex_valid_resample(tempIter,:) = tempseqIndex_valid_resample;
        end
        temp_F_dff_decisionBin_resample;
        seqIndex_valid_resample;
        
        temp_F_dff_decisionBin_resample_allSessions{currentSessionIndex} = temp_F_dff_decisionBin_resample;
        seqIndex_valid_resample_allSessions{currentSessionIndex} = seqIndex_valid_resample(1,:);
        seqCount_valid_allSessions{currentSessionIndex} = temp_seqCount;
    end
end
a = 1;


%% Collapse all sessions
seqIndex_valid_resample_allSessions_collapsed = zeros(length(seqIndex_valid_resample_allSessions{1}),num_FOV_AB);
seqCount_valid_allSessions_collapsed = zeros(length(seqCount_valid_allSessions{1}),num_FOV_AB);
for tempi=1:num_FOV_AB
    seqIndex_valid_resample_allSessions_collapsed(:,tempi) = seqIndex_valid_resample_allSessions{tempi};
    seqCount_valid_allSessions_collapsed(:,tempi) = seqCount_valid_allSessions{tempi};    
end
a = 1;
temp_F_dff_decisionBin_resample_allSessions_collaped = zeros(resampleIterCount,sum(roiNum_allSessions),length(seqIndex_valid_resample_allSessions{1}));
for tempi=1:num_FOV_AB
    temp_range_roi = (sum(roiNum_allSessions(1:tempi-1))+1):sum(roiNum_allSessions(1:tempi));
    [temp_1,temp_2] = bounds(temp_range_roi);
    
    temp_F_dff_decisionBin_resample_allSessions_collaped(:,temp_range_roi,:) = ...
        temp_F_dff_decisionBin_resample_allSessions{tempi};    
end
a = 1;
seqIndex_valid_resample_allSessions_collapsed_raw = seqIndex_valid_resample_allSessions_collapsed;
temp_F_dff_decisionBin_resample_allSessions_collaped_raw = temp_F_dff_decisionBin_resample_allSessions_collaped;

seqCount_valid_allSessions_collapsed;


seqCount_valid_allSessions_collapsed_min = min(seqCount_valid_allSessions_collapsed,[],2);

trialNum_resapmled_raw = size(seqIndex_valid_resample_allSessions_collapsed,1);
trialValidBoolIndex = false(trialNum_resapmled_raw,1);

% seq_valid_allSessions = seqCount_valid_allSessions_collapsed_min > 3;%0
seq_valid_allSessions = seqCount_valid_allSessions_collapsed_min >= resampleTrialCount_min;
for tempi=1:sum(numSeq(1:3))
    %temp_resampleTrialCount = resampleTrialCount;
    %temp_resampleTrialCount = resampleTrialCount_seq(tempi);

    temptempBool = seq_valid_allSessions(tempi);
    %temp_range_resampledTrialIndex = (temp_resampleTrialCount*(tempi-1)+1):temp_resampleTrialCount*tempi;
    temp_range_resampledTrialIndex = (sum(resampleTrialCount_seq(1:tempi-1))+1):sum(resampleTrialCount_seq(1:tempi));
    
    if temptempBool == true
        trialValidBoolIndex(temp_range_resampledTrialIndex) = true;
    end
end
a = 1;
a1 = sum(seq_valid_allSessions)

seqIndex_valid_resample_allSessions_collapsed = seqIndex_valid_resample_allSessions_collapsed_raw(trialValidBoolIndex,1);
temp_F_dff_decisionBin_resample_allSessions_collaped = temp_F_dff_decisionBin_resample_allSessions_collaped_raw(:,:,trialValidBoolIndex);

% temp1 = temp_F_dff_decisionBin_resample_allSessions_collaped_raw(:,:,~trialValidBoolIndex);
% temp2 = sum(temp1,'all');
% a1 = squeeze(temp1(1,:,:));


seqIndex_valid_resample_allSessions_collapsed;
temp_F_dff_decisionBin_resample_allSessions_collaped;

roiNum_allSessions;
roiNum_allSessions_sum = sum(roiNum_allSessions);
num_FOV_AB;
resampleIterCount;
% resampleTrialCount;
resampleTrialCount_seq;
seq_range = 1:sum(numSeq(1:3));

decodingDataSimplified_allSessions_resampled = struct;
decodingDataSimplified_allSessions_resampled.seqIndex_valid_resample_allSessions_collapsed = seqIndex_valid_resample_allSessions_collapsed;
decodingDataSimplified_allSessions_resampled.temp_F_dff_decisionBin_resample_allSessions_collaped = temp_F_dff_decisionBin_resample_allSessions_collaped;
decodingDataSimplified_allSessions_resampled.seqCount_valid_allSessions_collapsed = seqCount_valid_allSessions_collapsed;
decodingDataSimplified_allSessions_resampled.seqCount_valid_allSessions_collapsed_min = seqCount_valid_allSessions_collapsed_min;

decodingDataSimplified_allSessions_resampled.roiNum_allSessions = roiNum_allSessions;
decodingDataSimplified_allSessions_resampled.roiNum_allSessions_sum = roiNum_allSessions_sum;
decodingDataSimplified_allSessions_resampled.num_FOV = num_FOV_AB;
decodingDataSimplified_allSessions_resampled.resampleIterCount = resampleIterCount;
% decodingDataSimplified_allSessions_resampled.resampleTrialCount = resampleTrialCount;
decodingDataSimplified_allSessions_resampled.resampleTrialCount_seq = resampleTrialCount_seq;
decodingDataSimplified_allSessions_resampled.seq_range = seq_range;

decodingDataSimplified_allSessions_resampled.seqIndex_valid_resample_allSessions_collapsed_raw = seqIndex_valid_resample_allSessions_collapsed_raw;
decodingDataSimplified_allSessions_resampled.temp_F_dff_decisionBin_resample_allSessions_collaped_raw = temp_F_dff_decisionBin_resample_allSessions_collaped_raw;
decodingDataSimplified_allSessions_resampled.trialValidBoolIndex = trialValidBoolIndex;

decodingDataSimplified_allSessions_resampled.target_seqSet_inOne = temp_decodingDataSimplified.target_seqSet_inOne;
decodingDataSimplified_allSessions_resampled.target_seqSet = temp_decodingDataSimplified.target_seqSet;


seqCount_min_length = cell(3,1);
for tempi=1:3
    temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
    
    seqCount_min_length{tempi} = seqCount_valid_allSessions_collapsed_min(temp_range)';
end
seqCount_min_length

%% Save
if if_save == 1    
    if if_resample_seqBased0_locBased1 == 0
        temp_fileName = 'allSessions_resampled_seqBased';        
    elseif if_resample_seqBased0_locBased1 == 1
        temp_fileName = 'allSessions_resampled_locBased';        
    end
    
    if if_monkey_D0_Z1 == 0
        temp_fileName = [temp_fileName '_monkeyD'];
    elseif if_monkey_D0_Z1 == 1
        temp_fileName = [temp_fileName '_monkeyZ'];        
    end    
    
    save([output_path,'\',temp_fileName],'decodingDataSimplified_allSessions_resampled');
end

%% End
