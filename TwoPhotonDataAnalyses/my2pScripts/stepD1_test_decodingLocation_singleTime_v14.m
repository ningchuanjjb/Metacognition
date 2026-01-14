%% Initialization
% clear
% close all


% currentSession = '113Recording_20230426A_Ding_Site16';
% currentSession = '113Recording_20230427A_Ding_Site16_sameFOV0426';
% currentSession = '113Recording_20230502A_Ding_Site13';
% currentSession = '113Recording_20230503A_Ding_Site13_sameFOV0502';
% currentSession = '113Recording_20230504A_Ding_Site02';
% currentSession = '113Recording_20230509A_Ding_Site02';
% currentSession = '113Recording_20230508A_Ding_Site02_sameFOV0509';
% currentSession = '113Recording_20230511A_Ding_Site05';
% currentSession = '113Recording_20230510A_Ding_Site05_sameFOV0511';
% currentSession = '113Recording_20230510B_Ding_Site05_sameFOV0511';
% currentSession = '113Recording_20230512A_Ding_Site09';
% currentSession = '113Recording_20230513A_Ding_Site09_sameFOV0512';
% currentSession = '113Recording_20230515A_Ding_Site24_sameFOV0516';
% currentSession = '113Recording_20230516A_Ding_Site24';
% currentSession = '113Recording_20230517A_Ding_Site16B';
% currentSession = '113Recording_20230522A_Ding_Site05B';
% currentSession = '113Recording_20230525A_Ding_Site09B';
% currentSession = '113Recording_20230526A_Ding_Site09B_sameFOV0525';
% currentSession = '113Recording_20230527A_Ding_Site02B';
currentSession = '113Recording_20230528A_Ding_Site02B_sameFOV0527';
% currentSession = '113Recording_20230530A_Ding_Site05C';
% currentSession = '113Recording_20230531A_Ding_Site05C_sameFOV0530';
% currentSession = '113Recording_20230601A_Ding_Site13B';
% currentSession = '113Recording_20230602A_Ding_Site13B_sameFOV0601';
% currentSession = '113Recording_20230604A_Ding_Site07';
% currentSession = '113Recording_20230605A_Ding_Site07_sameFOV0604';
% currentSession = '113Recording_20230612A_Ding_Site14';
% currentSession = '113Recording_20230614A_Ding_Site05D';
% currentSession = '113Recording_20230615A_Ding_Site05D_sameFOV0614';
% currentSession = '113Recording_20230619A_Ding_Site02C';


targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)


% output_shortPath = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected';
output_shortPath = 'D:\twoPhotonData_motionCorrected';
temp_currentSession_path = [output_shortPath '\' currentSession];
temp_if_max0_min1 = 0;
output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);


%% Params setting
if_compute = 1;%1
if_plot = 1;

if_profile = 0;

if_seqLevel0_trialLevel1 = 1;
if_decoder_6location0_56seq1 = 0;%0

behavior_rProb0_gAcc1_trialCount2 = 1;
if_plot_all0_cluster1 = 0;
clusterNum = 10;

numFrames = 6;
rng(1); % For reproducibility
% t_decoder = templateSVM('Standardize',true); % To standardise input data
% t_decoder = templateSVM('Standardize',true,'KernelFunction','linear'); % To standardise input data
t_decoder = templateSVM('Standardize',true,'KernelFunction','linear'); % To standardise input data
% t_decoder = templateSVM('Standardize',true,'KernelFunction','polynomial','PolynomialOrder',2); % To standardise input data
% t_decoder = templateKNN('NumNeighbors',30,'Standardize',1); % To standardise input data


if if_profile == 1
    profile on
end

%% Load data
trial_para = []; %#ok<*NASGU>

% load decodingData
% load trial_para
% load markerParse_trialLevel
load([output_path,'\','decodingData.mat'],'decodingData');
load([output_path,'\','trial_para.mat'],'trial_para');
load([output_path,'\','markerParse_trialLevel.mat'],'markerParse_trialLevel');


temp_fields = fieldnames(decodingData);
for tempi=1:length(temp_fields)
    eval([temp_fields{tempi},'=decodingData.(temp_fields{tempi});']);
end


%% Others
a = 1;


numFrames = 6;
pointKindsNum = 4;
target_seqSet = get_target_seqSet_v1(numFrames,pointKindsNum);
numSeq = zeros(1,pointKindsNum);
for tempi=1:pointKindsNum
    numSeq(tempi) = length(target_seqSet{tempi});
end

F_dff_decisionBin1 = mean(F_dff_decisionPeriod(:,:,decisionPeriod_interval(1):decisionPeriod_interval(2)),3);
% F_dff_decisionBin1 = mean(F_dff_decisionPeriod(:,:,1:27),3);
% F_dff_decisionBin1 = mean(F_dff_decisionPeriod(:,:,decisionPeriod_interval(2):decisionPeriod_interval(3)),3);


% temp_dff = F_dff_decisionPeriod(:,:,decisionPeriod_interval(1):decisionPeriod_interval(2)-1);
% temp_dff2 = reshape(temp_dff,size(temp_dff,1),size(temp_dff,2),[],5);
temp_dff = F_dff_decisionPeriod(:,:,decisionPeriod_interval(1):decisionPeriod_interval(2));
temp_dff2 = reshape(temp_dff,size(temp_dff,1),size(temp_dff,2),[],2);
temp_dff2_mean = squeeze(mean(temp_dff2,3));
F_dff_decisionBin1_littleBin = temp_dff2_mean;


% topX = 50;
% F_dff_decisionBin1 = F_dff_decisionBin1(I_sorted(1:topX),:);


temp_F_dff_decisionBin = F_dff_decisionBin1;

roiNum = size(temp_F_dff_decisionBin,1);

a = 1;
if if_seqLevel0_trialLevel1 == 0
    temp_F_dff_decisionBin_seqMerged = zeros(roiNum,sum(numSeq));
    temp_F_dff_decisionBin;
    validSeqBoolIndex = true(1,sum(numSeq));
    for tempi=1:sum(numSeq)
        
        temp_seqBoolIndex = seqIndex==tempi;
        temp_boolIndex = temp_seqBoolIndex & trialIndex_bool_memoryCorrect;
        %temp_boolIndex = temp_seqBoolIndex;
        
        %temp_dff = temp_F_dff_decisionBin(:,seqIndex==tempi);
        temp_dff = temp_F_dff_decisionBin(:,temp_boolIndex);
        
        if isempty(temp_dff) == true
            %temp_dff = temp_F_dff_decisionBin(:,seqIndex==tempi);
            validSeqBoolIndex(tempi) = false;
        end
        
        temp_F_dff_decisionBin_seqMerged(:,tempi) = mean(temp_dff,2);
    end
    temp_F_dff_decisionBin_raw = temp_F_dff_decisionBin;
    
    temp_F_dff_decisionBin = temp_F_dff_decisionBin_seqMerged;
    
    
    
    boolIndex_location = false(numFrames,length(target_seqSet_inOne));
    for tempi=1:length(target_seqSet_inOne)
        currentSequence = target_seqSet_inOne{tempi};
        boolIndex_location(currentSequence,tempi) = true;
    end
    if size(temp_F_dff_decisionBin,2) ~= sum(validSeqBoolIndex)
        temp_F_dff_decisionBin_raw = temp_F_dff_decisionBin;
        temp_F_dff_decisionBin = temp_F_dff_decisionBin_raw(:,validSeqBoolIndex);
        boolIndex_location_raw = boolIndex_location;
        boolIndex_location = boolIndex_location_raw(:,validSeqBoolIndex);
    end
    
elseif if_seqLevel0_trialLevel1 == 1
    
%     temp_trialIndex_bool = trialIndex_lengthx_bool(2,:);
    %temp_trialIndex_bool = trialIndex_lengthx_bool_memoryCorrect(1,:);
%     temp_trialIndex_bool = trialIndex_lengthx_bool_memoryCorrect(2,:);
    %temp_trialIndex_bool = trialIndex_lengthx_bool_memoryCorrect(1,:) | trialIndex_lengthx_bool_memoryCorrect(2,:) | trialIndex_lengthx_bool_memoryCorrect(3,:);    
    temp_trialIndex_bool = trialIndex_bool_memoryCorrect; % use memory correct trial
%     temp_trialIndex_bool = true(1,size(F_dff_decisionBin1,2)); % use all trial
    %temp_trialIndex_bool = trialIndex_bool_choiceMemory; % use choice memory trial
    
    temp_F_dff_decisionBin = F_dff_decisionBin1(:,temp_trialIndex_bool);
    
    temp_F_dff_decisionBin_littleBin = F_dff_decisionBin1_littleBin(:,temp_trialIndex_bool,:);
    
    seqIndex_valid = seqIndex(temp_trialIndex_bool);
    
    boolIndex_location_trial = false(numFrames,length(seqIndex_valid));
    for tempi=1:length(seqIndex_valid)
        currentSequence = target_seqSet_inOne{seqIndex_valid(tempi)};
        boolIndex_location_trial(currentSequence,tempi) = true;
    end
    
    boolIndex_location_seq = false(numFrames,length(target_seqSet_inOne));
    for tempi=1:length(target_seqSet_inOne)
        currentSequence = target_seqSet_inOne{tempi};
        boolIndex_location_seq(currentSequence,tempi) = true;
    end
    
end


% validSeqBoolIndex = ~isnan(sum(temp_F_dff_decisionBin,1));

a = 1;
% if if_profile == 1
%    profile on
% end
if if_compute == 1
    t_fun_seqProbSVM_6location = tic;
    if if_seqLevel0_trialLevel1 == 0
        seqProbSVM_6location_seqLevel_Name_v = autoGetFunName_myScripts('seqProbSVM_6location_seqLevel', [targetPATH '\functions']);
        fun_seqProbSVM_6location_seqLevel = str2func(seqProbSVM_6location_seqLevel_Name_v);
        svm_seqPosterior = fun_seqProbSVM_6location_seqLevel(temp_F_dff_decisionBin,boolIndex_location,t_decoder,numFrames);
    elseif if_seqLevel0_trialLevel1 == 1
        if if_decoder_6location0_56seq1 == 0
            seqProbSVM_6location_trialLevel_Name_v = autoGetFunName_myScripts('seqProbSVM_6location_trialLevel', [targetPATH '\functions']);
            fun_seqProbSVM_6location_trialLevel = str2func(seqProbSVM_6location_trialLevel_Name_v);
            
            if_resample = 0;
            resampleTrialCount = 8;%4-->2-->8
            resampleIterCount = 64;%10-->64

            if if_resample == 0
                [svm_seqPosterior_raw,validSeqBoolIndex] = fun_seqProbSVM_6location_trialLevel(temp_F_dff_decisionBin,boolIndex_location_trial,boolIndex_location_seq,t_decoder,numFrames);
                svm_seqPosterior = svm_seqPosterior_raw(validSeqBoolIndex);
            elseif if_resample == 1                
                
                
                
                uniqueSeq = unique(seqIndex_valid);
                temp1 = 1:sum(numSeq);
                validSeqBoolIndex = ismember(temp1,uniqueSeq);
                
                seqCount_all = zeros(1,sum(numSeq));
                tempj = 0;
                for tempi=1:sum(numSeq)
                    %if validSeqBoolIndex(tempi) == false
                    %if ismember(tempi,uniqueSeq) == 0
                    %    continue
                    %end
                    %tempj = tempj + 1;
                    %seqCount_all(tempi) = sum(seqIndex_valid==uniqueSeq(tempj),'all');
                    seqCount_all(tempi) = sum(seqIndex_valid==tempi,'all');
                end
                
                validSeqBoolIndex2 = seqCount_all>=resampleTrialCount;
                temp1 = 1:sum(numSeq);
                uniqueSeq2 = temp1(validSeqBoolIndex2);
                
                
                %svm_seqPosterior_resampleIter = zeros(resampleIterCount,sum(validSeqBoolIndex2));   
                svm_seqPosterior_resampleIter = zeros(resampleIterCount,sum(numSeq));   
                %for tempIter=1:resampleIterCount
                parfor tempIter=1:resampleIterCount
                    
                    tempIndex_resample = zeros(resampleTrialCount,sum(validSeqBoolIndex2));
                    tempj = 0;
                    for tempi=1:sum(numSeq)
                        if validSeqBoolIndex2(tempi) == false
                            continue
                        end
                        tempj = tempj + 1;
                        tempIndex = find(seqIndex_valid==uniqueSeq2(tempj));
                        tempI = sort(randperm(length(tempIndex),resampleTrialCount));
                        tempIndex_resample(:,tempj) = tempIndex(tempI);
                    end
                    tempIndex_resample_1d = reshape(tempIndex_resample,1,[]);
                    
                    temp_F_dff_decisionBin_resample = temp_F_dff_decisionBin(:,tempIndex_resample_1d);
                    boolIndex_location_trial_resample = boolIndex_location_trial(:,tempIndex_resample_1d);
                    
                    [svm_seqPosterior_raw_resample,validSeqBoolIndex_resample] = ...
                        fun_seqProbSVM_6location_trialLevel(temp_F_dff_decisionBin_resample,boolIndex_location_trial_resample,boolIndex_location_seq,t_decoder,numFrames);
                    
                    %svm_seqPosterior_resampleIter(tempIter,:) = svm_seqPosterior_raw_resample(validSeqBoolIndex2);
                    svm_seqPosterior_resampleIter(tempIter,:) = svm_seqPosterior_raw_resample;
                    a = 1;
                end
                svm_seqPosterior = mean(svm_seqPosterior_resampleIter,1);
                validSeqBoolIndex_raw = validSeqBoolIndex;
                validSeqBoolIndex = validSeqBoolIndex2;
                svm_seqPosterior_raw = svm_seqPosterior;
                svm_seqPosterior = svm_seqPosterior_raw(validSeqBoolIndex);
            end
        elseif if_decoder_6location0_56seq1 == 1
            seqProbSVM_56seq_trialLevel_Name_v = autoGetFunName_myScripts('seqProbSVM_56seq_trialLevel', [targetPATH '\functions']);
            fun_seqProbSVM_56seq_trialLevel = str2func(seqProbSVM_56seq_trialLevel_Name_v);
            %[svm_seqPosterior_raw,validSeqBoolIndex] = fun_seqProbSVM_56seq_trialLevel(temp_F_dff_decisionBin,boolIndex_location_trial,boolIndex_location_seq,t_decoder,numFrames,temp_F_dff_decisionBin_littleBin);            
            %svm_seqPosterior = svm_seqPosterior_raw(validSeqBoolIndex); 
            
                        
            threshold_trialCount = 10;%10, 7, 
            
            temp_svm_Y_raw = boolIndex_location_trial;
            temp_svm_Y_raw2 = zeros(1,size(temp_svm_Y_raw,2));
            for tempi=1:size(boolIndex_location_seq,2)
                tempSeq_bool = boolIndex_location_seq(:,tempi)';
                temp1_boolIndex = sum(temp_svm_Y_raw' == tempSeq_bool,2) == numFrames;
                temp_svm_Y_raw2(temp1_boolIndex) = tempi;
            end
            uniqueSeq = unique(temp_svm_Y_raw2);
            seqCount = zeros(1,length(uniqueSeq));
            for tempi=1:length(uniqueSeq)
                seqCount(tempi) = sum(temp_svm_Y_raw2==uniqueSeq(tempi),'all');
            end
            
            lowTrialCount_seqIndex = uniqueSeq(seqCount<threshold_trialCount);
            
            validSeqBoolIndex = false(1,sum(numSeq));
            validSeqBoolIndex(uniqueSeq) = true;
            validSeqBoolIndex(lowTrialCount_seqIndex) = false;
            
            validTrialBoolIndex = ismember(temp_svm_Y_raw2,find(validSeqBoolIndex==true));
            validTrialCount = sum(validTrialBoolIndex);
            
            [svm_seqPosterior_raw,validSeqBoolIndex] = fun_seqProbSVM_56seq_trialLevel(temp_F_dff_decisionBin(:,validTrialBoolIndex),boolIndex_location_trial(:,validTrialBoolIndex),boolIndex_location_seq,t_decoder,numFrames,temp_F_dff_decisionBin_littleBin(:,validTrialBoolIndex,:));
            svm_seqPosterior = svm_seqPosterior_raw(validSeqBoolIndex);                        
                        
        end
    end
    % toc(t_fun_seqProbSVM_6location);
    fprintf('t_fun_seqProbSVM_6location = %.1f secs.\n',toc(t_fun_seqProbSVM_6location));
elseif if_compute == 0
    %validSeqBoolIndex = validSeqBoolIndex2;
    %validSeqBoolIndex(7:end) = false;
    %svm_seqPosterior = svm_seqPosterior_raw(validSeqBoolIndex);
end
a = 1;
% if if_profile == 1
%    profile viewer
% end

if behavior_rProb0_gAcc1_trialCount2 == 0
    temp_behaviorValue = offloadingProb_inOne;
elseif behavior_rProb0_gAcc1_trialCount2 == 1
    temp_behaviorValue = 1-gAcc_noChoice_collapsed_inOne;
elseif behavior_rProb0_gAcc1_trialCount2 == 2
    
    seqIndex_valid;
    uniqueSeq = unique(seqIndex_valid);
    seqCount_all = zeros(1,sum(numSeq));
    tempj = 0;
    for tempi=1:sum(numSeq)
        %if validSeqBoolIndex(tempi) == false
        %    continue
        %end
        %tempj = tempj + 1;
        %seqCount_all(tempi) = sum(seqIndex_valid==uniqueSeq(tempj),'all');
        seqCount_all(tempi) = sum(seqIndex_valid==tempi,'all');
    end
    temp_behaviorValue = seqCount_all;
end

a = 1;

if_shuffle = 1;
shuffleNum = 30;

if if_shuffle == 1
    if if_compute == 1
        t_shuffle = tic;
        svm_seqPosterior_shuffle = zeros(shuffleNum,sum(numSeq));
        if if_seqLevel0_trialLevel1 == 0
            parfor tempi=1:shuffleNum
                shuffleIndex = randperm(size(boolIndex_location,2));
                svm_seqPosterior_shuffle(tempi,:) = fun_seqProbSVM_6location_seqLevel (...
                    temp_F_dff_decisionBin(:,shuffleIndex),...
                    boolIndex_location(:,shuffleIndex),t_decoder,numFrames); %#ok<*PFBNS>
            end
        elseif if_seqLevel0_trialLevel1 == 1
            parfor tempi=1:shuffleNum
                
                % shuffle in seq level
                tempBoolIndex = find(validSeqBoolIndex==true);
                shuffleIndex_1 = randperm(sum(validSeqBoolIndex));
                shuffleIndex_seq = tempBoolIndex(shuffleIndex_1);
                
                shuffleIndex = zeros(1,size(boolIndex_location_trial,2));
                for tempj=1:length(seqIndex_valid)
                    currentSeqIndex = seqIndex_valid(tempj);
                    currentSeqIndex_shuffled = shuffleIndex_seq(find(tempBoolIndex==currentSeqIndex));                 %#ok<*FNDSB>
                    shuffleIndex(tempj) = currentSeqIndex_shuffled;
                end
                
                % shuffle in trial level
                %shuffleIndex = randperm(size(boolIndex_location_trial,2));
                
                [temp_svm_seqPosterior_raw,temp_validSeqBoolIndex] = fun_seqProbSVM_6location_trialLevel (...
                    temp_F_dff_decisionBin,...
                    boolIndex_location_trial(:,shuffleIndex),...
                    boolIndex_location_seq,t_decoder,numFrames); %#ok<*PFBNS>
                svm_seqPosterior_shuffle(tempi,:) = temp_svm_seqPosterior_raw;
            end
        end
        svm_seqPosterior_shuffle_median = median(svm_seqPosterior_shuffle,1);
        fprintf('t_shuffle = %.1f secs.\n',toc(t_shuffle));
        
        
        % svm_seqPosterior_n11n = svm_seqPosterior./svm_seqPosterior_shuffle_median;
        % svm_seqPosterior_n11n = (svm_seqPosterior-svm_seqPosterior_shuffle_median)...
        %     ./svm_seqPosterior_shuffle_median;
        svm_seqPosterior_n11n = svm_seqPosterior_raw-svm_seqPosterior_shuffle_median;
        
    elseif if_compute == 0
        svm_seqPosterior_n11n = svm_seqPosterior_raw-svm_seqPosterior_shuffle_median;
    end
elseif if_shuffle == 0
    svm_seqPosterior_n11n = svm_seqPosterior_raw;
end

% svm_seqPosterior_n11n = svm_seqPosterior_raw(validSeqBoolIndex)-svm_seqPosterior_shuffle_median;


% x = svm_seqPosterior(validSeqBoolIndex);
x = svm_seqPosterior_n11n(validSeqBoolIndex);
y = temp_behaviorValue(validSeqBoolIndex);
mdl_all = fitglm(x,y);
% mdl_all = fitglm(x,y,'Link','log');
r2_all = mdl_all.Rsquared.Adjusted;
fprintf('r2_all = %.4f.\n',r2_all);

[r_corr,p_corr] = corr(x',y');

numSeq_cumu = zeros(1,length(numSeq));
for tempi=1:length(numSeq)
    numSeq_cumu(tempi) = sum(numSeq(1:tempi));    
end
lengthxSeqBoolIndex = false(length(numSeq),sum(numSeq));
lengthxSeqBoolIndex_valid = false(length(numSeq),sum(numSeq));
for tempi=1:length(numSeq)
    if tempi == 1
        lengthxSeqBoolIndex(tempi,1:numSeq_cumu(tempi)) = true;
    else
        lengthxSeqBoolIndex(tempi,numSeq_cumu(tempi-1)+1:numSeq_cumu(tempi)) = true;
    end
    lengthxSeqBoolIndex_valid(tempi,:) = lengthxSeqBoolIndex(tempi,:) & validSeqBoolIndex;
end
lengthxSeqBoolIndex;
lengthxSeqBoolIndex_valid;

x_length = zeros(1,length(numSeq));
y_length = zeros(1,length(numSeq));
for tempi=1:length(numSeq)
    x_length(tempi) = mean(svm_seqPosterior_raw(lengthxSeqBoolIndex_valid(tempi,:)));
    y_length(tempi) = mean(temp_behaviorValue(lengthxSeqBoolIndex_valid(tempi,:)));
end
% mdl_length = fitglm(x_length,y_length);
% r2_length = mdl_length.Rsquared.Adjusted;


[temp_behaviorValue_sort,temp_I] = sort(temp_behaviorValue);
temp1 = 1:sum(numSeq);
validSeqIndex = temp1(validSeqBoolIndex);
temp2BoolIndex = ismember(temp_I,validSeqIndex);
temp_I_valid = temp_I(temp2BoolIndex);
temp_behaviorValue_sort_valid = temp_behaviorValue_sort(temp2BoolIndex);


% clusterNum = 10;

clusterxSeqBoolIndex = false(clusterNum,sum(numSeq));
clusterSize = floor(length(temp_I_valid)/clusterNum)*ones(1,clusterNum);
for tempi=1:(length(temp_I_valid)-sum(clusterSize))
    clusterSize(clusterNum+1-tempi) = clusterSize(clusterNum+1-tempi) + 1;
end


x_cluster = zeros(1,clusterNum);
y_cluster = zeros(1,clusterNum);
for tempi=1:clusterNum
    temp_range = sum(clusterSize(1:tempi-1))+1:sum(clusterSize(1:tempi));
    temp_range2 = temp_I_valid(temp_range);
    x_cluster(tempi) = mean(svm_seqPosterior_raw(temp_range2));    
    y_cluster(tempi) = mean(temp_behaviorValue(temp_range2));
    a = 1;
end
% mdl_cluster = fitglm(x_cluster,y_cluster);
% r2_cluster = mdl_cluster.Rsquared.Adjusted;

% if_plot_all0_cluster1 = 0;



%% Plot
if if_plot == 1
    close all
    
    fig1 = figure('Name','Fig1','NumberTitle','off');
    set(gcf,'Position',[35+0 35+0 600 450]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    
    nexttile
    if if_plot_all0_cluster1 == 0
        x = svm_seqPosterior_n11n(validSeqBoolIndex);
        %x = svm_seqPosterior_raw(validSeqBoolIndex);
        y = temp_behaviorValue(validSeqBoolIndex);
    elseif if_plot_all0_cluster1 == 1
        x = x_cluster;
        y = y_cluster;
    end
    
    %x = svm_seqPosterior_raw(validSeqBoolIndex);
    %y = seqCount_all(validSeqBoolIndex)/max(seqCount_all);
    
    
    mdl = fitglm(x,y);
    r2 = mdl.Rsquared.Adjusted;

    scatter(x, y, 15, 'filled',...
        'MarkerFaceColor', [0 0 0], 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.8);%0.4
    hold on
    
    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y);

    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    
    xlim([temp_xmin temp_xmax]);
    ylim([temp_ymin temp_ymax]);
    
    set(gca, 'FontSize', 11)
    set(gca,'box','off');% 取消右、上边框
    %xlabel(sprintf('svm posterior probability, r2=%.3f',r2_all), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(sprintf('decoding accuracy, r2=%.3f',r2), 'FontSize', 14, 'FontWeight', 'bold');
    if if_plot_all0_cluster1 == 0
        title(sprintf('%d seqs',sum(validSeqBoolIndex)), 'FontSize', 14, 'FontWeight', 'bold');
    elseif if_plot_all0_cluster1 == 1
        title(sprintf('%d clusters',clusterNum), 'FontSize', 14, 'FontWeight', 'bold');
    end
    if behavior_rProb0_gAcc1_trialCount2 == 0        
        ylabel(sprintf('offloading rate'), 'FontSize', 14, 'FontWeight', 'bold');
    elseif behavior_rProb0_gAcc1_trialCount2 == 1
        ylabel(sprintf('Error rate'), 'FontSize', 14, 'FontWeight', 'bold');        
    elseif behavior_rProb0_gAcc1_trialCount2 == 2
        ylabel(sprintf('Trial count'), 'FontSize', 14, 'FontWeight', 'bold');                
    end
end

if if_profile == 1
    profile viewer
end


%% End