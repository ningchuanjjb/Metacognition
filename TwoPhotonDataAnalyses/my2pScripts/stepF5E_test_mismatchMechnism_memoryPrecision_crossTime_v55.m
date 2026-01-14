% Chuan's 14th script (20251214)
% This script: Time courses of memory strength, related to figure 4.
%% Initialization
close all

if_plot = 1;

if_plot_multiPanel = 0;
if_plot_singlePanel = 0;
if_plot_multiPanelB = 0;%1
if_plot_multiPanelC = 0;

if_dynamics_precision0_locCorr1 = 1;%0,1

if_compute = 1;%1

if_smooth = 0;%0, 1 is only temopory for within-seq plot
if_StatQuant = 1;%1, 0 is only temopory for within-seq plot

if_fastCompute = 0;%0, for within-seq plot

if if_fastCompute == 1
    if_plot = 0;
    
    if_smooth = 1;    
    if_StatQuant = 0;
end



if_delay1_forward0_backward1 = 0; % new

% if if_monkey_D0_Z1 == 1
%     if_delay1_forward0_backward1 = 1;
% end

if exist('if_compute_summary','var') == 1
    if_delay1_forward0_backward1 = 1;
end

if_ttest_crossTime = 1;

if_allTrial0_matchMismatch1_allTrialWithError2_correct3 = 2;%3

% if_trainPeriod_eachBin0_delay1Bin1 = 0;

if_plot_fineTuningMismatch = 0;

p_threshold = 0.01;


color_choiceMemoryHigh = [128,205,193]/255; %[146,197,222]/255
color_choiceMemoryLow = [1,133,113]/255; %[5,113,176]/255
color_choiceOffloadHigh = [166,97,26]/255; %[202,0,32]/255
color_choiceOffloadLow = [223,194,125]/255; %[244,165,130]/255
color_choiceMemory = [1,133+20,113]/255;
color_choiceOffload = [166+40,97,26]/255;
color_choiceMemoryError = [0.3 0.3 0.3];

if if_memeoryPrecision_stimuli0_response1 == 0
    temp_seqIndex = seqIndex;
elseif if_memeoryPrecision_stimuli0_response1 == 1
    temp_seqIndex = seqIndex_response;
end


if_4typesFilter_all0_cMC1_cF2_cME3 = 0;%0

if if_4typesFilter_all0_cMC1_cF2_cME3 == 0
    tempBoolIndex_4typesFilter = true(length(choiceMemoryCorrectBoolIndex),1);
    
elseif if_4typesFilter_all0_cMC1_cF2_cME3 == 1
    tempBoolIndex_4typesFilter = choiceMemoryCorrectBoolIndex';
    
elseif if_4typesFilter_all0_cMC1_cF2_cME3 == 2
    tempBoolIndex_4typesFilter = choiceOffloadBoolIndex';
    
elseif if_4typesFilter_all0_cMC1_cF2_cME3 == 3
    tempBoolIndex_4typesFilter = choiceMemoryErrorBoolIndex';
end


memoryPrecision_trialLevel;
memoryPrecision_trialLevel;

if_profileOn = 0;
if if_profileOn == 1
    profile on
end
%% Test memoryPrecision
t_test = tic;
if if_compute == 1
    fun_testMemoryPrecision_currentTime_Name_v = autoGetFunName_myScripts('fun_testMemoryPrecision_currentTime', [targetPATH '\functions']);
    fun_testMemoryPrecision_currentTime = str2func(fun_testMemoryPrecision_currentTime_Name_v);
    
    options_testMemoryPrecision = struct;
    
    options_testMemoryPrecision.trialIndex_bool_memoryCorrect = trialIndex_bool_memoryCorrect;
    options_testMemoryPrecision.allMemoryCorrectBoolIndex = allMemoryCorrectBoolIndex;
    options_testMemoryPrecision.temp_trialIndex_valid_resample = temp_trialIndex_valid_resample;
    options_testMemoryPrecision.svm_train_length1_outputs = svm_train_length1_outputs;
    options_testMemoryPrecision.svm_train_length2_outputs = svm_train_length2_outputs;
    options_testMemoryPrecision.svm_train_length3_outputs = svm_train_length3_outputs;
    options_testMemoryPrecision.numSeq = numSeq;
    options_testMemoryPrecision.resampleIterCount = resampleIterCount;
    %options_testMemoryPrecision.resampleIterCount = 1;%resampleIterCount-->2-->1
    options_testMemoryPrecision.seqIndex = seqIndex;
    options_testMemoryPrecision.boolIndex_location_seq = boolIndex_location_seq;
    options_testMemoryPrecision.seqIndex_valid = seqIndex_valid;
    
    options_testMemoryPrecision.allMemoryErrorBoolIndex = allMemoryErrorBoolIndex;
    options_testMemoryPrecision.seqIndex_response = seqIndex_response;
    options_testMemoryPrecision.KFold_num = KFold_num;
    options_testMemoryPrecision.numFrames = numFrames;
    options_testMemoryPrecision.if_memeoryPrecision_stimuli0_response1 = if_memeoryPrecision_stimuli0_response1;
    %options_testMemoryPrecision.if_memeoryPrecision_stimuli0_response1 = 1;
    options_testMemoryPrecision.if_inferOffloadResponse = if_inferOffloadResponse;
    
    options_testMemoryPrecision.choiceOffloadBoolIndex = choiceOffloadBoolIndex;
    
    options_testMemoryPrecision.if_compute = if_compute;
    options_testMemoryPrecision.if_fastCompute = if_fastCompute;
    options_testMemoryPrecision.if_memoryPrecision_accuracy0_sigma1 = if_memoryPrecision_accuracy0_sigma1;
    options_testMemoryPrecision.score_stimuli_to_response = score_stimuli_to_response;
    options_testMemoryPrecision.if_precision_meanProb0_sumProb1 = if_precision_meanProb0_sumProb1;
    options_testMemoryPrecision.if_entropy = if_entropy;
    options_testMemoryPrecision.fun_sigmaBased_singleTrialPrecision = fun_sigmaBased_singleTrialPrecision;    
    
    F_dff_baselinePeriod;
    F_dff_length1_sample;
    F_dff_length2_sample;
    F_dff_length3_sample;
    F_dff_decisionPeriodA;
    
    for loopCount=0:4
    %for loopCount=0:0
    
    %memoryPrecision_trialLevel_crossTime_baseline = [];
    %for loopCount=1:1
        %if true
        %    loopCount = 4;
        
        fprintf('loopCount=%d.\n',loopCount);
        
        if loopCount==0
            F_dff_current = F_dff_baselinePeriod;
        elseif loopCount==1
            F_dff_current = F_dff_length1_sample;
        elseif loopCount==2
            F_dff_current = F_dff_length2_sample;
        elseif loopCount==3
            F_dff_current = F_dff_length3_sample;
        elseif loopCount==4
            if if_delay1_forward0_backward1 == 1
                F_dff_current = F_dff_decisionPeriodA;
            elseif if_delay1_forward0_backward1 == 0
                F_dff_current = F_dff_decisionPeriodA;
                F_dff_current(:,:,1:decisionPeriodB_interval(2)) = F_dff_decisionPeriodB;
            end            
        end
        
        if if_memoryPrecision_accuracy0_sigma1 == 1
            if if_smooth == 1
                %Fast trial: 10(698,677),20(683,653),30(680,644),40(680,639),50(680,634),60(680,629)
                %Full trial: 10(707,677),20(697,653),30(695,641)
                F_dff_current = smoothdata(F_dff_current,3,'gaussian',20);
            end
        end        
        
        
        temptemp_range = 1:size(F_dff_current,3);
        
        temp_memoryPrecision_trialLevel_crossTime = zeros(length(seqIndex),length(temptemp_range));
        temp_locCorr_trialLevel_crossTime = zeros(length(seqIndex),length(temptemp_range));
        
        temp_memoryPrecision_trialLevel_resampleIter_crossTime = zeros(length(seqIndex),resampleIterCount,length(temptemp_range));
        temp_locCorr_trialLevel_resampleIter_crossTime = zeros(length(seqIndex),resampleIterCount,length(temptemp_range));
        
        for tempi=temptemp_range(1):temptemp_range(end)        
        %for tempi=temptemp_range(1):temptemp_range(1)
        %parfor tempi=temptemp_range(1):temptemp_range(end)
        %    rng(5);%5
            F_dff_currentTime = F_dff_current(:,:,temptemp_range(tempi)); %#ok<*PFBNS>
            testMemoryPrecision_output = fun_testMemoryPrecision_currentTime(F_dff_currentTime,options_testMemoryPrecision);
            temp_memoryPrecision_trialLevel_crossTime(:,tempi) = testMemoryPrecision_output.memoryPrecision_trialLevel_currentTime;
            temp_locCorr_trialLevel_crossTime(:,tempi) = testMemoryPrecision_output.locCorr_trialLevel_currentTime;            
            
            temp_memoryPrecision_trialLevel_resampleIter_crossTime(:,:,tempi) = testMemoryPrecision_output.memoryPrecision_trialLevel_resampleIter_matrix_currentTime;
            temp_locCorr_trialLevel_resampleIter_crossTime(:,:,tempi) = testMemoryPrecision_output.locCorr_trialLevel_resampleIter_matrix_currentTime;            
        end
        
        if loopCount==0
            memoryPrecision_trialLevel_crossTime_baseline_raw2 = temp_memoryPrecision_trialLevel_crossTime;
            locCorr_trialLevel_crossTime_baseline = temp_locCorr_trialLevel_crossTime;
            
            memoryPrecision_trialLevel_resampleIter_crossTime_baseline_raw2 = temp_memoryPrecision_trialLevel_resampleIter_crossTime;
            locCorr_trialLevel_resampleIter_crossTime_baseline = temp_locCorr_trialLevel_resampleIter_crossTime;
            
        elseif loopCount==1
            memoryPrecision_trialLevel_crossTime_length1_raw2 = temp_memoryPrecision_trialLevel_crossTime;
            locCorr_trialLevel_crossTime_length1 = temp_locCorr_trialLevel_crossTime;            
            
            memoryPrecision_trialLevel_resampleIter_crossTime_length1_raw2 = temp_memoryPrecision_trialLevel_resampleIter_crossTime;
            locCorr_trialLevel_resampleIter_crossTime_length1 = temp_locCorr_trialLevel_resampleIter_crossTime;
            
        elseif loopCount==2
            memoryPrecision_trialLevel_crossTime_length2_raw2 = temp_memoryPrecision_trialLevel_crossTime;
            locCorr_trialLevel_crossTime_length2 = temp_locCorr_trialLevel_crossTime;            
            
            memoryPrecision_trialLevel_resampleIter_crossTime_length2_raw2 = temp_memoryPrecision_trialLevel_resampleIter_crossTime;
            locCorr_trialLevel_resampleIter_crossTime_length2 = temp_locCorr_trialLevel_resampleIter_crossTime;
            
        elseif loopCount==3
            memoryPrecision_trialLevel_crossTime_length3_raw2 = temp_memoryPrecision_trialLevel_crossTime;
            locCorr_trialLevel_crossTime_length3 = temp_locCorr_trialLevel_crossTime;            
            
            memoryPrecision_trialLevel_resampleIter_crossTime_length3_raw2 = temp_memoryPrecision_trialLevel_resampleIter_crossTime;
            locCorr_trialLevel_resampleIter_crossTime_length3 = temp_locCorr_trialLevel_resampleIter_crossTime;
            
        elseif loopCount==4
            memoryPrecision_trialLevel_crossTime_delay1_raw2 = temp_memoryPrecision_trialLevel_crossTime;
            locCorr_trialLevel_crossTime_delay1 = temp_locCorr_trialLevel_crossTime;            
            
            memoryPrecision_trialLevel_resampleIter_crossTime_delay1_raw2 = temp_memoryPrecision_trialLevel_resampleIter_crossTime;
            locCorr_trialLevel_resampleIter_crossTime_delay1 = temp_locCorr_trialLevel_resampleIter_crossTime;
            
        end
        
        fprintf('t = %.1f secs.\n',toc(t_test));
    end
    
    
    
%     memoryPrecision_trialLevel_crossTime_baseline_raw;
%     
%     temp_baseline_mean = mean(memoryPrecision_trialLevel_crossTime_baseline_raw,'all','omitnan');
%     temp_baseline_std = std(memoryPrecision_trialLevel_crossTime_baseline_raw,0,'all','omitnan');
%     
%     memoryPrecision_trialLevel_crossTime_baseline_raw = (memoryPrecision_trialLevel_crossTime_baseline_raw-temp_baseline_mean)./temp_baseline_std;    
%     memoryPrecision_trialLevel_crossTime_length1_raw = (memoryPrecision_trialLevel_crossTime_length1_raw-temp_baseline_mean)./temp_baseline_std;
%     memoryPrecision_trialLevel_crossTime_length2_raw = (memoryPrecision_trialLevel_crossTime_length2_raw-temp_baseline_mean)./temp_baseline_std;
%     memoryPrecision_trialLevel_crossTime_length3_raw = (memoryPrecision_trialLevel_crossTime_length3_raw-temp_baseline_mean)./temp_baseline_std;
%     memoryPrecision_trialLevel_crossTime_delay1_raw = (memoryPrecision_trialLevel_crossTime_delay1_raw-temp_baseline_mean)./temp_baseline_std;    
    

    memoryPrecision_trialLevel_crossTime_baseline_raw = memoryPrecision_trialLevel_crossTime_baseline_raw2;
    memoryPrecision_trialLevel_crossTime_length1_raw = memoryPrecision_trialLevel_crossTime_length1_raw2;
    memoryPrecision_trialLevel_crossTime_length2_raw = memoryPrecision_trialLevel_crossTime_length2_raw2;
    memoryPrecision_trialLevel_crossTime_length3_raw = memoryPrecision_trialLevel_crossTime_length3_raw2;
    memoryPrecision_trialLevel_crossTime_delay1_raw = memoryPrecision_trialLevel_crossTime_delay1_raw2;
    
    
    %% Z-score
    temp_seqIndex;    
    for tempi=1:3
        temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
        
        temptempBoolIndex = ismember(temp_seqIndex,temp_range);
        
        temp_baseline_mean = mean(memoryPrecision_trialLevel_crossTime_baseline_raw2(temptempBoolIndex,:),'all','omitnan');
        temp_baseline_std = std(memoryPrecision_trialLevel_crossTime_baseline_raw2(temptempBoolIndex,:),0,'all','omitnan');
        
        memoryPrecision_trialLevel_crossTime_baseline_raw(temptempBoolIndex,:) = (memoryPrecision_trialLevel_crossTime_baseline_raw2(temptempBoolIndex,:)-temp_baseline_mean)./temp_baseline_std;
        memoryPrecision_trialLevel_crossTime_length1_raw(temptempBoolIndex,:) = (memoryPrecision_trialLevel_crossTime_length1_raw2(temptempBoolIndex,:)-temp_baseline_mean)./temp_baseline_std;
        memoryPrecision_trialLevel_crossTime_length2_raw(temptempBoolIndex,:) = (memoryPrecision_trialLevel_crossTime_length2_raw2(temptempBoolIndex,:)-temp_baseline_mean)./temp_baseline_std;
        memoryPrecision_trialLevel_crossTime_length3_raw(temptempBoolIndex,:) = (memoryPrecision_trialLevel_crossTime_length3_raw2(temptempBoolIndex,:)-temp_baseline_mean)./temp_baseline_std;
        memoryPrecision_trialLevel_crossTime_delay1_raw(temptempBoolIndex,:) = (memoryPrecision_trialLevel_crossTime_delay1_raw2(temptempBoolIndex,:)-temp_baseline_mean)./temp_baseline_std;        
    end

    
end

if false
    memoryPrecision_trialLevel_crossTime_baseline_raw = memoryPrecision_trialLevel_crossTime_baseline;
    memoryPrecision_trialLevel_crossTime_length1_raw = memoryPrecision_trialLevel_crossTime_length1;
    memoryPrecision_trialLevel_crossTime_length2_raw = memoryPrecision_trialLevel_crossTime_length2;
    memoryPrecision_trialLevel_crossTime_length3_raw = memoryPrecision_trialLevel_crossTime_length3;
    memoryPrecision_trialLevel_crossTime_delay1_raw = memoryPrecision_trialLevel_crossTime_delay1;
end

if if_fastCompute == 0
    if if_dynamics_precision0_locCorr1 == 0
        memoryPrecision_trialLevel_crossTime_baseline = memoryPrecision_trialLevel_crossTime_baseline_raw;
        memoryPrecision_trialLevel_crossTime_length1 = memoryPrecision_trialLevel_crossTime_length1_raw;
        memoryPrecision_trialLevel_crossTime_length2 = memoryPrecision_trialLevel_crossTime_length2_raw;
        memoryPrecision_trialLevel_crossTime_length3 = memoryPrecision_trialLevel_crossTime_length3_raw;
        memoryPrecision_trialLevel_crossTime_delay1 = memoryPrecision_trialLevel_crossTime_delay1_raw;
        
        memoryPrecision_trialLevel_resampleIter_crossTime_baseline = memoryPrecision_trialLevel_resampleIter_crossTime_baseline_raw2;
        memoryPrecision_trialLevel_resampleIter_crossTime_length1 = memoryPrecision_trialLevel_resampleIter_crossTime_length1_raw2;
        memoryPrecision_trialLevel_resampleIter_crossTime_length2 = memoryPrecision_trialLevel_resampleIter_crossTime_length2_raw2;
        memoryPrecision_trialLevel_resampleIter_crossTime_length3 = memoryPrecision_trialLevel_resampleIter_crossTime_length3_raw2;
        memoryPrecision_trialLevel_resampleIter_crossTime_delay1 = memoryPrecision_trialLevel_resampleIter_crossTime_delay1_raw2;
        
    elseif if_dynamics_precision0_locCorr1 == 1
        memoryPrecision_trialLevel_crossTime_baseline = locCorr_trialLevel_crossTime_baseline;
        memoryPrecision_trialLevel_crossTime_length1 = locCorr_trialLevel_crossTime_length1;
        memoryPrecision_trialLevel_crossTime_length2 = locCorr_trialLevel_crossTime_length2;
        memoryPrecision_trialLevel_crossTime_length3 = locCorr_trialLevel_crossTime_length3;
        memoryPrecision_trialLevel_crossTime_delay1 = locCorr_trialLevel_crossTime_delay1;
        
        memoryPrecision_trialLevel_resampleIter_crossTime_baseline = locCorr_trialLevel_resampleIter_crossTime_baseline;
        memoryPrecision_trialLevel_resampleIter_crossTime_length1 = locCorr_trialLevel_resampleIter_crossTime_length1;
        memoryPrecision_trialLevel_resampleIter_crossTime_length2 = locCorr_trialLevel_resampleIter_crossTime_length2;
        memoryPrecision_trialLevel_resampleIter_crossTime_length3 = locCorr_trialLevel_resampleIter_crossTime_length3;
        memoryPrecision_trialLevel_resampleIter_crossTime_delay1 = locCorr_trialLevel_resampleIter_crossTime_delay1;
        
    end
end

if if_fastCompute == 1
    fun_testMemoryPrecision_currentTime_Name_v = autoGetFunName_myScripts('fun_testMemoryPrecision_currentTime', [targetPATH '\functions']);
    fun_testMemoryPrecision_currentTime = str2func(fun_testMemoryPrecision_currentTime_Name_v);
    
    options_testMemoryPrecision = struct;
    
    options_testMemoryPrecision.trialIndex_bool_memoryCorrect = trialIndex_bool_memoryCorrect;
    options_testMemoryPrecision.allMemoryCorrectBoolIndex = allMemoryCorrectBoolIndex;
    options_testMemoryPrecision.temp_trialIndex_valid_resample = temp_trialIndex_valid_resample;
    options_testMemoryPrecision.svm_train_length1_outputs = svm_train_length1_outputs;
    options_testMemoryPrecision.svm_train_length2_outputs = svm_train_length2_outputs;
    options_testMemoryPrecision.svm_train_length3_outputs = svm_train_length3_outputs;
    options_testMemoryPrecision.numSeq = numSeq;
    options_testMemoryPrecision.resampleIterCount = resampleIterCount;
    options_testMemoryPrecision.seqIndex = seqIndex;
    options_testMemoryPrecision.boolIndex_location_seq = boolIndex_location_seq;
    options_testMemoryPrecision.seqIndex_valid = seqIndex_valid;
    
    options_testMemoryPrecision.allMemoryErrorBoolIndex = allMemoryErrorBoolIndex;
    options_testMemoryPrecision.seqIndex_response = seqIndex_response;
    options_testMemoryPrecision.KFold_num = KFold_num;
    options_testMemoryPrecision.numFrames = numFrames;
    options_testMemoryPrecision.if_memeoryPrecision_stimuli0_response1 = if_memeoryPrecision_stimuli0_response1;
    %options_testMemoryPrecision.if_memeoryPrecision_stimuli0_response1 = 0;%?
    %options_testMemoryPrecision.if_memeoryPrecision_stimuli0_response1 = 1;
    options_testMemoryPrecision.if_inferOffloadResponse = if_inferOffloadResponse;
    
    options_testMemoryPrecision.choiceOffloadBoolIndex = choiceOffloadBoolIndex;
    
    options_testMemoryPrecision.if_compute = if_compute;
    options_testMemoryPrecision.if_fastCompute = if_fastCompute;
    options_testMemoryPrecision.if_memoryPrecision_accuracy0_sigma1 = if_memoryPrecision_accuracy0_sigma1;
    options_testMemoryPrecision.score_stimuli_to_response = score_stimuli_to_response;
    options_testMemoryPrecision.if_precision_meanProb0_sumProb1 = if_precision_meanProb0_sumProb1;
    options_testMemoryPrecision.if_entropy = if_entropy;
    options_testMemoryPrecision.fun_sigmaBased_singleTrialPrecision = fun_sigmaBased_singleTrialPrecision;    
    
    
    F_dff_baselinePeriod;
    F_dff_length1_sample;
    F_dff_length2_sample;
    F_dff_length3_sample;
    F_dff_decisionPeriodA;
    
    %for loopCount=0:4
    %if true
    loopCount = 4;
    
    %fprintf('loopCount=%d.\n',loopCount);
    
    if loopCount==0
        F_dff_current = F_dff_baselinePeriod;
    elseif loopCount==1
        F_dff_current = F_dff_length1_sample;
    elseif loopCount==2
        F_dff_current = F_dff_length2_sample;
    elseif loopCount==3
        F_dff_current = F_dff_length3_sample;
    elseif loopCount==4
        if if_delay1_forward0_backward1 == 1
            F_dff_current = F_dff_decisionPeriodA;
        elseif if_delay1_forward0_backward1 == 0
            F_dff_current = F_dff_decisionPeriodA;
            F_dff_current(:,:,1:decisionPeriodB_interval(2)) = F_dff_decisionPeriodB;
        end
        
        F_dff_current;
        if if_memoryPrecision_accuracy0_sigma1 == 1
            %Fast trial: 10(698,677),20(683,653),30(680,644),40(680,639),50(680,634),60(680,629)
            %Full trial: 10(707,677),20(697,653),30(695,641)
            F_dff_current = smoothdata(F_dff_current,3,'gaussian',20);
        end

    end
    
    temptemp_range = 1:size(F_dff_current,3);
    
    temp_memoryPrecision_trialLevel_crossTime = zeros(length(seqIndex),length(temptemp_range));
    for tempi=temptemp_range(1):temptemp_range(end)
        %parfor tempi=temptemp_range(1):temptemp_range(end)
        %for tempi=temptemp_range(1):temptemp_range(1)
        
        F_dff_currentTime = F_dff_current(:,:,temptemp_range(tempi)); %#ok<*PFBNS>
        %F_dff_currentTime = F_dff_decisionBin1; % only for test
        %F_dff_currentTime = mean(F_dff_current,3); % only for test
        
        testMemoryPrecision_output = fun_testMemoryPrecision_currentTime(F_dff_currentTime,options_testMemoryPrecision);
        if if_memoryPrecision_accuracy0_sigma1 == 0
            temp_memoryPrecision_trialLevel_crossTime(:,tempi) = testMemoryPrecision_output.memoryPrecision_trialLevel_currentTime;
        elseif if_memoryPrecision_accuracy0_sigma1 == 1
            temp_memoryPrecision_trialLevel_crossTime(:,tempi) = testMemoryPrecision_output.memoryPrecision_trialLevel_currentTime;
            %temp_memoryPrecision_trialLevel_crossTime(:,tempi) = testMemoryPrecision_output.locCorr_trialLevel_currentTime;
        end
    end
    
    if loopCount==0
        memoryPrecision_trialLevel_crossTime_baseline = temp_memoryPrecision_trialLevel_crossTime;
    elseif loopCount==1
        memoryPrecision_trialLevel_crossTime_length1 = temp_memoryPrecision_trialLevel_crossTime;
    elseif loopCount==2
        memoryPrecision_trialLevel_crossTime_length2 = temp_memoryPrecision_trialLevel_crossTime;
    elseif loopCount==3
        memoryPrecision_trialLevel_crossTime_length3 = temp_memoryPrecision_trialLevel_crossTime;
    elseif loopCount==4
        memoryPrecision_trialLevel_crossTime_delay1 = temp_memoryPrecision_trialLevel_crossTime;
    end
    
    fprintf('t = %.1f secs.\n',toc(t_test));
    %end
    
    memoryPrecision_trialLevel_crossTime_baseline = [];
    memoryPrecision_trialLevel_crossTime_length1 = [];
    memoryPrecision_trialLevel_crossTime_length2 = [];
    memoryPrecision_trialLevel_crossTime_length3 = [];
    
end

memoryPrecision_trialLevel_crossTime_baseline;
memoryPrecision_trialLevel_crossTime_length1;
memoryPrecision_trialLevel_crossTime_length2;
memoryPrecision_trialLevel_crossTime_length3;
memoryPrecision_trialLevel_crossTime_delay1;




%% Statistical quantification of temporal dynamics
% if if_compute == 1
if true
    if if_StatQuant == 1
        memoryPrecision_trialLevel_resampleIter_crossTime_baseline;
        memoryPrecision_trialLevel_resampleIter_crossTime_length1;
        memoryPrecision_trialLevel_resampleIter_crossTime_length2;
        memoryPrecision_trialLevel_resampleIter_crossTime_length3;
        memoryPrecision_trialLevel_resampleIter_crossTime_delay1;
        
        a = 1;
        trialBoolIndex_length1;
        
        memoryPrecision_trialLevel_resampleIter_crossTime_lastT = nan(size(memoryPrecision_trialLevel_resampleIter_crossTime_length1));
        
        memoryPrecision_trialLevel_resampleIter_crossTime_lastT(trialBoolIndex_length1,:,:) = ...
            memoryPrecision_trialLevel_resampleIter_crossTime_length1(trialBoolIndex_length1,:,end-17:end);
        
        memoryPrecision_trialLevel_resampleIter_crossTime_lastT(trialBoolIndex_length2,:,:) = ...
            memoryPrecision_trialLevel_resampleIter_crossTime_length2(trialBoolIndex_length2,:,end-17:end);
        
        memoryPrecision_trialLevel_resampleIter_crossTime_lastT(trialBoolIndex_length3,:,:) = ...
            memoryPrecision_trialLevel_resampleIter_crossTime_length3(trialBoolIndex_length3,:,end-17:end);        
        
        
        temp_size_baseline = size(memoryPrecision_trialLevel_resampleIter_crossTime_baseline,3);        
        temp_size = size(memoryPrecision_trialLevel_resampleIter_crossTime_lastT,3) + ...
            size(memoryPrecision_trialLevel_resampleIter_crossTime_delay1(:,:,1:33),3);
        
        temp_r_baseline = nan(resampleIterCount,temp_size_baseline);
        temp_p = nan(resampleIterCount,temp_size);
        temp_r = nan(resampleIterCount,temp_size);
        for tempi=1:resampleIterCount
            
            temp1_baseline = squeeze(memoryPrecision_trialLevel_resampleIter_crossTime_baseline(:,tempi,:));
            temp1_baseline_choiceMemoryCorrect = temp1_baseline(choiceMemoryCorrectBoolIndex,:);
            
            temp1_baseline_choiceMemoryCorrect_median = median(temp1_baseline_choiceMemoryCorrect,'all','omitnan');
            if isnan(temp1_baseline_choiceMemoryCorrect_median)
                continue
            end
            
            temp1_baseline_choiceMemoryCorrect_median2 = mean(temp1_baseline_choiceMemoryCorrect,1,'omitnan');
            
            temp1_lastT = squeeze(memoryPrecision_trialLevel_resampleIter_crossTime_lastT(:,tempi,:));
            temp1_lastT_choiceMemoryCorrect = temp1_lastT(choiceMemoryCorrectBoolIndex,:);

            temp1_delay1 = squeeze(memoryPrecision_trialLevel_resampleIter_crossTime_delay1(:,tempi,:));
            temp1_delay1_choiceMemoryCorrect = temp1_delay1(choiceMemoryCorrectBoolIndex,:);
            

            temp2 = [temp1_lastT_choiceMemoryCorrect,temp1_delay1_choiceMemoryCorrect(:,1:33)];
            
            %% Baseline
            temp_r_baseline(tempi,:) = mean(temp1_baseline,1,'omitnan');
            
            
            %% LastT ~ Delay1
            temptemp_p = nan(1,size(temp2,2));
            for temptempi=1:size(temp2,2)
                %[~,temptemp_p(temptempi)] = ttest(temp2(:,temptempi),temp1_baseline_choiceMemoryCorrect_median);
                [~,temptemp_p(temptempi)] = ttest(temp1_baseline_choiceMemoryCorrect_median2,mean(temp2(:,temptempi),'omitnan'),'tail','left');                
            end
            temp_p(tempi,:) = temptemp_p;
            temp_r(tempi,:) = mean(temp2,1,'omitnan');
        end
        temp_p_valid = temp_p(~isnan(temp_p(:,1)),:);        
                
        temp_signifTimeStamp = temp_p_valid < p_threshold;
        
        temp_signifTimeStamp_memory = temp_signifTimeStamp;
        
        temp_r_dynamics_memory_baseline = temp_r_baseline;
        temp_r_dynamics_memory = temp_r;  
        
        
        temp_signifTimeStamp_memory;
        
        temp_persistentTimeStamp = 5;%5
        
        temp_signifTimeStampCount_memory = nan(size(temp_signifTimeStamp_memory,1),1);
        for tempi=1:length(temp_signifTimeStampCount_memory)
            temp1 = temp_signifTimeStamp_memory(tempi,:);
            
            for tempj=1:(length(temp1)-temp_persistentTimeStamp+1)
               if sum(temp1(tempj:tempj+temp_persistentTimeStamp-1)) == temp_persistentTimeStamp
                   temp_signifTimeStampCount_memory(tempi) = tempj;
                   break
               end                
            end
            
        end
        temp_signifTimeStampCount_memory;
        
        if isnan(temp_signifTimeStampCount_memory)
            temp_signifTimeStampCount_memory = length(temp_signifTimeStamp_memory);
        end
        
        temp_signifTimeStampCount_memory_mean = mean(temp_signifTimeStampCount_memory,'omitnan');
        temp_signifTimeStampCount_memory_sem = std(temp_signifTimeStampCount_memory,0,'omitnan')./sqrt(length(temp_signifTimeStampCount_memory));

        temp_r_dynamics_memory_baseline_mean = mean(temp_r_dynamics_memory_baseline,1,'omitnan');        
        temp_r_dynamics_memory_mean = mean(temp_r_dynamics_memory,1,'omitnan');
        
        
        
        %% Plot
        temp_signifTimeStampCount_memory_mean;
        
        temp_r_dynamics_memory_whole = [temp_r_dynamics_memory_baseline,temp_r_dynamics_memory];
        
        temp_r_dynamics_memory_whole_mean = mean(temp_r_dynamics_memory_whole,1,'omitnan');
        temp_r_dynamics_memory_whole_sem = std(temp_r_dynamics_memory_whole,0,'omitnan')./sqrt(size(temp_r_dynamics_memory_whole,1));
        temp_r_dynamics_memory_whole_std = std(temp_r_dynamics_memory_whole,0,'omitnan');
        
        temp1 = temp_r_dynamics_memory_whole;
        temp1_mean = temp_r_dynamics_memory_whole_mean;
        temp1_sem = temp_r_dynamics_memory_whole_sem;
        temp1_std = temp_r_dynamics_memory_whole_std;
        
        temp2_timeStamp = temp_signifTimeStampCount_memory_mean + size(temp_r_dynamics_memory_baseline,2);        
        temp2 = zeros(1,size(temp1,2));
        temp2(ceil(temp2_timeStamp):end) = 1;
        
        if if_plot == 1
            backgounrdColor = [1 1 1];
            
            temp_blankSize = 10;            

            fig = figure('Name','asd','NumberTitle','off');
            temptemp1 = 360*0.975*0.74*1.37*0.97*1.02*1.015;
            temptemp2 = 200*0.8*0.8*0.87;
            
            %set(gcf,'Position',[600 50+350 temptemp1 temptemp2]);
            set(gcf,'Position',[600 50+350 temptemp1*0.84 temptemp2]);
            
            t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
            
            nexttile
            

            
            x = 1:baselinePeriod_interval(end);
            y = temp1_mean(1:baselinePeriod_interval(end));
            plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1.5);
            hold on
            
%             y_std = temp1_std(1:baselinePeriod_interval(end));
%             patch([x(:); flipud(x(:))]', [y(:)+y_std(:); flipud(y(:)-y_std(:))]',[0.3 0.3 0.3],'FaceAlpha',0.05,'EdgeColor',[1 1 1]*0.725);%[1 1 1]*0.725           
%             hold on
            
            x = ((baselinePeriod_interval(end)+1):size(temp1,2))+temp_blankSize;            
            y = temp1_mean((baselinePeriod_interval(end)+1):size(temp1,2));
            plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1.5);
            hold on        
            
%             y_std = temp1_std((baselinePeriod_interval(end)+1):size(temp1,2));
%             patch([x(:); flipud(x(:))]', [y(:)+y_std(:); flipud(y(:)-y_std(:))]',[0.3 0.3 0.3],'FaceAlpha',0.05,'EdgeColor',[1 1 1]*0.725);            
%             hold on
            
            
            x = (1:size(temp1,2))+temp_blankSize;
            y = temp1_mean;
            
            [y_min,y_max] = bounds(y);
            
            if if_ttest_crossTime == 1
                scatter(x(temp2==1),y_max+(y_max-y_min)*0.15,8,[0 0 0],'*');
            end
            
            
            temp_interval_all = [baselinePeriod_interval(1),...
                (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)),...
                (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)) + temp_blankSize,...
                (baselinePeriod_interval(end) + temp_blankSize + length1_sample_interval(end) + decisionPeriodA_interval)];
                        
            text(mean(temp_interval_all(2:3)),y_min-(y_max-y_min)*0.1,...
                sprintf('%s','//'),'HorizontalAlignment','center','fontsize',8);
                        
            xticks(temp_interval_all([1 2 3 4]));
            xticklabels({'Fixation','T1','LastT','Delay-on'});
            
            yticks([0 0.5]);
            if if_monkey_D0_Z1 == 1
                yticks([0 0.2]);
            end
            
            xlim([temp_interval_all(1) temp_interval_all(end-1)]);
            
            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 7.5)%12
            set(gca,'color',backgounrdColor);
            set(gca,'box','off');% 取消右、上边框
            
            
            ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.2]);
            ylabel('Correlation', 'FontSize', 8, 'FontWeight', 'normal');%10
            
            %title('Memory quality time course','FontSize',9);
            title('Memory time course','FontSize',9);
            subtitle('Correct trials, location correlation','FontSize',7.5);
            
            xtickangle(0);
            
            
        end
        
        %% Compare temp_signifTimeStampCount_meta & temp_signifTimeStampCount_memory
        if if_plot == 1
            if true
                
                
                temp_1 = temp_signifTimeStampCount_memory;
                temp_2 = temp_signifTimeStampCount_meta;

                [~,temp_p_signifTimeStamp_memoryMeta] = ttest2(temp_1,temp_2);
                
                
                fig = figure('Name','example seq distribution','NumberTitle','off');
                set(gcf,'Position',[680 450 201.4*1.08*0.95*0.98*0.97*0.65*0.8*1.5 179.5*1.08*0.97*1.5*1.3*1.05*0.43*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                
                t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
                
                nexttile
                
                temp_p12 = temp_p_signifTimeStamp_memoryMeta;
                
                temp_1 = (temp_signifTimeStampCount_memory-1).*33.3+rand(16,1)./1000;
                %temp_1 = temp_signifTimeStampCount_meta.*33;
                
                temp_2 = (temp_signifTimeStampCount_meta-1).*33.3;
                
                fprintf('signifTimeStampCount_memory: mean = %.3f, std = %.3f.\n',mean(temp_1),std(temp_1));
                fprintf('signifTimeStampCount_meta: mean = %.3f, std = %.3f.\n',mean(temp_2),std(temp_2));
                
                temp_y_min = min([temp_1;temp_2]);
                temp_y_max = max([temp_1;temp_2]);
                
                
                temp_data = [temp_1;temp_2];
                
                g1 = repmat({'A'},length(temp_1),1);
                g2 = repmat({'B'},length(temp_2),1);
                
                temp_label = [g1;g2];
                
                temptemp_color1 = [1 1 1]*0.5;
                temptemp_color2 = repmat(temptemp_color1, 2, 1);
                
                h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                    'GroupOrder',[{'A'};{'B'}]);
                h(1).ViolinPlot.FaceAlpha = 0.1;
                h(2).ViolinPlot.FaceAlpha = 0.1;
                
                
                tempTxt = sprintf('');
                if temp_p12 < 0.001
                    tempTxt = sprintf('***');
                elseif temp_p12 < 0.01
                    tempTxt = sprintf('**');
                elseif temp_p12 < 0.05
                    tempTxt = sprintf('*');
                end
                text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.20,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
                    'HorizontalAlignment','center');
                plot([1.1 1.9],temp_y_max+(temp_y_max-temp_y_min)*0.16*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
                hold on
                
                
                set(gca,'linewidth',1.5)
                xlim([0.45 2.55])
                ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.30]);
                %ylim([-0.125 1.30])
                set(gca, 'FontSize', 8);%10
                
                xticks([1 2]);
                
                xtl = ["Memory";"Meta"];
                xt=get(gca,'XTick');
                yt=get(gca,'YTick');
                xtext_xp=xt;                
                xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.28;
                text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',8);%9
                set(gca,'xticklabel','');
                
                
                %yticks([0 1]);
                
                set(gca,'box','off');% 取消右、上边框
                ylabel('Latency time (ms)', 'FontSize', 8, 'FontWeight', 'normal');%10
                %     if if_memeoryPrecision_stimuli0_response1 == 0
                %         title(sprintf('Stimuli-labeled'),'fontsize',9);
                %     elseif if_memeoryPrecision_stimuli0_response1 == 1
                %         title(sprintf('Response-labeled'),'fontsize',9);
                %     end
                
            end
            
        end
        
        
    end
    
end


% if_plot_multiPanel=1;
%% Plot
if if_plot == 1
    %close all
    
    if if_plot_multiPanel == 1
        backgounrdColor = [1 1 1]*1;%0.875-->0.915-->1
        
        fig = figure('Name','asd','NumberTitle','off'); %#ok<*NASGU>
        %set(gcf,'Position',[10 50 1650 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 50 720 200]);
        t = tiledlayout(1,5,'TileSpacing','Compact','Padding','Compact');
        
        t.Title.String = sprintf('MemoryPrecision dynamics, %s, p_threshold=%.3f',FOVName_currentFOV2,p_threshold);
        t.Title.FontSize = 10;%12
        t.Title.Interpreter = 'none';
        
        [y_min,y_max] = bounds([memoryPrecision_trialLevel_crossTime_baseline(:);...
            memoryPrecision_trialLevel_crossTime_length1(:);memoryPrecision_trialLevel_crossTime_length2(:);memoryPrecision_trialLevel_crossTime_length3(:);...
            memoryPrecision_trialLevel_crossTime_delay1(:)]);
        
        for loopCount=0:4
            %for loopCount=1:4
            
            if loopCount==0
                temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_baseline;
                temp_str = 'Baseline';
                %nexttile(1);%4
                nexttile
                temp_interval = baselinePeriod_interval;
            elseif loopCount==1
                temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_length1;
                temp_str = 'Sample, length1';
                %nexttile(4);%2
                nexttile
                temp_interval = length1_sample_interval;
            elseif loopCount==2
                temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_length2;
                temp_str = 'Sample, length2';
                %nexttile(5);%5
                nexttile
                temp_interval = length2_sample_interval;
            elseif loopCount==3
                temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_length3;
                temp_str = 'Sample, length3';
                %nexttile(6);%8
                nexttile
                temp_interval = length3_sample_interval;
            elseif loopCount==4
                temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_delay1;
                temp_str = 'Delay1';
                %nexttile(3);%6
                nexttile
                temp_interval = decisionPeriodA_interval;
            end
            
            
            % mismatch trials
            if if_plot_fineTuningMismatch == 1
                temp_memoryPrecision_trialLevel_crossTime_overMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh,:);
                temp_memoryPrecision_trialLevel_crossTime_underMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow,:);
            elseif if_plot_fineTuningMismatch == 0
                %temp_memoryPrecision_trialLevel_crossTime_overMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLowError_choiceMemory,:);
                %temp_memoryPrecision_trialLevel_crossTime_underMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_choiceOffload,:);
                temp_memoryPrecision_trialLevel_crossTime_overMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaHigh_choice,:);
                temp_memoryPrecision_trialLevel_crossTime_underMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaLow_choice,:);
            end
            
            % match trials
            if if_plot_fineTuningMismatch == 1
                temp_memoryPrecision_trialLevel_crossTime_highMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh,:);
                temp_memoryPrecision_trialLevel_crossTime_lowMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_choiceOffloadLow,:);
            elseif if_plot_fineTuningMismatch == 0
                %temp_memoryPrecision_trialLevel_crossTime_highMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory,:);
                %temp_memoryPrecision_trialLevel_crossTime_lowMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_choiceOffload,:);
                temp_memoryPrecision_trialLevel_crossTime_highMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaHigh_choice,:);
                temp_memoryPrecision_trialLevel_crossTime_lowMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaLow_choice,:);
            end
            
            % all trials
            temp_memoryPrecision_trialLevel_crossTime_choiceMemory = temp_memoryPrecision_trialLevel_crossTime(choiceMemoryBoolIndex,:);
            temp_memoryPrecision_trialLevel_crossTime_choiceOffload = temp_memoryPrecision_trialLevel_crossTime(choiceOffloadBoolIndex,:);
            
            % all trials with error
            temp_memoryPrecision_trialLevel_crossTime_choiceMemoryCorrect = temp_memoryPrecision_trialLevel_crossTime(choiceMemoryCorrectBoolIndex,:);
            temp_memoryPrecision_trialLevel_crossTime_choiceMemoryError = temp_memoryPrecision_trialLevel_crossTime(choiceMemoryErrorBoolIndex,:);
            
            a1 = temp_memoryPrecision_trialLevel_crossTime_highMatch;
            a2 = temp_memoryPrecision_trialLevel_crossTime_lowMatch;
            a3 = temp_memoryPrecision_trialLevel_crossTime_overMismatch;
            a4 = temp_memoryPrecision_trialLevel_crossTime_underMismatch;
            a5 = temp_memoryPrecision_trialLevel_crossTime_choiceMemory;
            a6 = temp_memoryPrecision_trialLevel_crossTime_choiceOffload;
            a7 = temp_memoryPrecision_trialLevel_crossTime_choiceMemoryCorrect;
            a8 = temp_memoryPrecision_trialLevel_crossTime_choiceMemoryError;
            
            a1_mean = mean(a1,1,'omitnan');
            a2_mean = mean(a2,1,'omitnan');
            a3_mean = mean(a3,1,'omitnan');
            a4_mean = mean(a4,1,'omitnan');
            a5_mean = mean(a5,1,'omitnan');
            a6_mean = mean(a6,1,'omitnan');
            a7_mean = mean(a7,1,'omitnan');
            a8_mean = mean(a8,1,'omitnan');            
            
            a1_sem = std(a1,1,1,'omitnan')./sqrt(size(a1,1));
            a2_sem = std(a2,1,1,'omitnan')./sqrt(size(a2,1));
            a3_sem = std(a3,1,1,'omitnan')./sqrt(size(a3,1));
            a4_sem = std(a4,1,1,'omitnan')./sqrt(size(a4,1));
            a5_sem = std(a5,1,1,'omitnan')./sqrt(size(a5,1));
            a6_sem = std(a6,1,1,'omitnan')./sqrt(size(a6,1));
            a7_sem = std(a7,1,1,'omitnan')./sqrt(size(a7,1));
            a8_sem = std(a8,1,1,'omitnan')./sqrt(size(a8,1));
            
            [~,p_a5_a6] = ttest2(a5,a6);
            
            range_crossTime = 1:size(temp_memoryPrecision_trialLevel_crossTime,2);
            x = range_crossTime;
            
            h_line = [];
            
            if if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 1
                y = a1_mean;
                h_line = [h_line plot(x,y,'color',color_choiceMemoryHigh,'linewidth',1)]; %#ok<*AGROW>
                hold on
                y_sem = a1_sem;
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryHigh,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                hold on
                
                y = a2_mean;
                h_line = [h_line plot(x,y,'color',color_choiceOffloadLow,'linewidth',1)];
                hold on
                y_sem = a2_sem;
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                hold on
                
                y = a3_mean;
                h_line = [h_line plot(x,y,'color',color_choiceMemoryLow,'linewidth',1)];
                hold on
                y_sem = a3_sem;
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                hold on
                
                y = a4_mean;
                h_line = [h_line plot(x,y,'color',color_choiceOffloadHigh,'linewidth',1)];
                hold on
                y_sem = a4_sem;
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadHigh,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                hold on
                
            elseif if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 0
                y = a5_mean;
                h_line = [h_line plot(x,y,'color',color_choiceMemory,'linewidth',1)]; %#ok<*AGROW>
                hold on
                y_sem = a5_sem;
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                hold on
                
                y = a6_mean;
                h_line = [h_line plot(x,y,'color',color_choiceOffload,'linewidth',1)];
                hold on
                y_sem = a6_sem;
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                hold on
                
                p_a5_a6;
                
                %scatter(x(p_a5_a6<0.001),y_max-(y_max-y_min)*0.61,8,[0 0 0],'*');
                %scatter(x(p_a5_a6<p_threshold),y_max-(y_max-y_min)*0.61,8,[0 0 0],'*');
                scatter(x(p_a5_a6<p_threshold),y_max-(y_max-y_min)*0.56,8,[0 0 0],'*');
            elseif if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 2
                y = a7_mean;
                h_line = [h_line plot(x,y,'color',color_choiceMemory,'linewidth',1)]; %#ok<*AGROW>
                hold on
                y_sem = a7_sem;
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                hold on
                
                
                y = a8_mean;
                h_line = [h_line plot(x,y,'color',color_choiceMemoryError,'linewidth',1)]; %#ok<*AGROW>
                hold on
                y_sem = a8_sem;
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryError,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                hold on     
                
                
                y = a6_mean;
                h_line = [h_line plot(x,y,'color',color_choiceOffload,'linewidth',1)];
                hold on
                y_sem = a6_sem;
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                hold on                                
                
            end
            
            temptempTrialNum_a1 = sum(~isnan(a1(:,1)));
            temptempTrialNum_a2 = sum(~isnan(a2(:,1)));
            temptempTrialNum_a3 = sum(~isnan(a3(:,1)));
            temptempTrialNum_a4 = sum(~isnan(a4(:,1)));
            temptempTrialNum_a5 = sum(~isnan(a5(:,1)));
            temptempTrialNum_a6 = sum(~isnan(a6(:,1)));
            
            
            if loopCount==1 || loopCount==2 || loopCount==3 || loopCount== 4
                for tempi=1:(length(temp_interval)-1)
                    plot(temp_interval(tempi)*[1 1],[y_min y_max],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                if loopCount==1 || loopCount==2 || loopCount==3
                    for tempi=1:(length(temp_interval)-1)
                        x = [temp_interval(1+tempi-1) y_min;...
                            temp_interval(1+tempi-1) y_max;...
                            temp_interval(1+tempi-1)+6 y_max;...
                            temp_interval(1+tempi-1)+6 y_min];
                        y = [1 2 3 4];
                        patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.1,'EdgeColor','none');
                        hold on
                    end
                end
            end
            
            if loopCount==0
                xticks(temp_interval(1:end-1));
                %xticklabels({'PreFixation','Fixation'});
                xticklabels({'Fixation',''});
            elseif loopCount==1
                xticks(temp_interval(1:end-1));
                xticklabels({'T1'});
            elseif loopCount==2
                xticks(temp_interval(1:end-1));
                xticklabels({'T1','T2'});
            elseif loopCount==3
                xticks(temp_interval(1:end-1));
                xticklabels({'T1','T2','T3'});
            elseif loopCount==4
                xticks(temp_interval(1:end-1));
                xticklabels({'Delay1','ChoiceCue'});
            end
            
            
            
            %if loopCount == 0
            %if loopCount == 0
            if false
                if if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 1 %#ok<*UNRCH>
                    le = legend(h_line,sprintf('highMatch'),...
                        sprintf('lowMatch'),...
                        sprintf('overMismatch'),...
                        sprintf('underMismatch'),...
                        'Location','northeast','fontsize',10,'NumColumns',4);
                    
                    le.ItemTokenSize = ones(1,4)*14;%20
                    le.Color = backgounrdColor;
                elseif if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 0
                    le = legend(h_line,sprintf('choiceMemory'),...
                        sprintf('choiceOffload'),...
                        'Location','northeast','fontsize',10,'NumColumns',4);
                end
            end
            
            xlim([range_crossTime(1) range_crossTime(end)]);
            %ylim([y_min-(y_max-y_min)*0.0 y_max-(y_max-y_min)*0.4]);%0.1
            %if if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 1
                %ylim([y_min-(y_max-y_min)*0.0 y_max-(y_max-y_min)*0.45]);%0.1
            %elseif if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 0
                %ylim([y_min-(y_max-y_min)*0.0 y_max-(y_max-y_min)*0.58]);%0.1
            %end
            ylim([y_min-(y_max-y_min)*0.0 y_max-(y_max-y_min)*0.5]);
            %ylim([0 1]);
            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 10)
            set(gca,'color',backgounrdColor);
            set(gca,'box','off');% 取消右、上边框
            %xlabel('Frame', 'FontSize', 14, 'FontWeight', 'bold');
            %if loopCount == 1
            if loopCount == 0
                ylabel('Memory Precision', 'FontSize', 10, 'FontWeight', 'normal');
            end
            
            xtickangle(0);
            
            %if if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 1
                %temp_title = title(sprintf('%s, [%d,%d,%d,%d]',...
                %    temp_str,temptempTrialNum_a1,temptempTrialNum_a2,temptempTrialNum_a3,temptempTrialNum_a4),...
                %    'FontSize',14,'FontWeight','normal');%'bold'
                %temp_title = title(sprintf('%s',temp_str),'FontSize',10,'FontWeight','normal');%'bold'
            %elseif if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 0
                %temp_title = title(sprintf('%s, [%d,%d]',...
                %    temp_str,temptempTrialNum_a5,temptempTrialNum_a6),...
                %    'FontSize',14,'FontWeight','normal');%'bold
                %temp_title = title(sprintf('%s',temp_str),'FontSize',10,'FontWeight','normal');%'bold
            %end
            temp_title = title(sprintf('%s',temp_str),'FontSize',10,'FontWeight','normal');
            temp_title.Interpreter = 'none';
            
            
        end
    end
    
    
    %% single panel
    if if_plot_singlePanel == 1
        %close all
        backgounrdColor = [1 1 1]*1;%0.875
        
        fig = figure('Name','asda','NumberTitle','off');
        %set(gcf,'Position',[10 50 350 230]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 340 130*1.4]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 340 130*1.15]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[35+0 42+0 340 130*1.15*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        %t = tiledlayout(1,2,'TileSpacing','tight','Padding','tight'); %#ok<*NASGU>
        t = tiledlayout(1,2,'TileSpacing','compact','Padding','tight'); %#ok<*NASGU>
        
        %t.Title.String = sprintf('Population');
        %t.Title.String = sprintf('Population, n=%d',size(F_dff_decisionPeriodA,1));
        %t.Title.FontSize = 11;
        %t.Title.Interpreter = 'none';
        
        
        nexttile
        
        [y_min,y_max] = bounds([memoryPrecision_trialLevel_crossTime_delay1(:)]); %#ok<*NBRAK>
        
        
        temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_delay1;
        temp_str = 'Delay1';
        temp_interval = decisionPeriodA_interval;
        
        
        % mismatch trials
        if if_plot_fineTuningMismatch == 1
            temp_memoryPrecision_trialLevel_crossTime_overMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh,:);
            temp_memoryPrecision_trialLevel_crossTime_underMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow,:);
        elseif if_plot_fineTuningMismatch == 0
            %temp_memoryPrecision_trialLevel_crossTime_overMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLowError_choiceMemory,:);
            %temp_memoryPrecision_trialLevel_crossTime_underMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_choiceOffload,:);
            temp_memoryPrecision_trialLevel_crossTime_overMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaHigh_choice,:);
            temp_memoryPrecision_trialLevel_crossTime_underMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaLow_choice,:);
        end
        
        % match trials
        if if_plot_fineTuningMismatch == 1
            temp_memoryPrecision_trialLevel_crossTime_highMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh,:);
            temp_memoryPrecision_trialLevel_crossTime_lowMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_choiceOffloadLow,:);
        elseif if_plot_fineTuningMismatch == 0
            %temp_memoryPrecision_trialLevel_crossTime_highMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory,:);
            %temp_memoryPrecision_trialLevel_crossTime_lowMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_choiceOffload,:);
            temp_memoryPrecision_trialLevel_crossTime_highMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaHigh_choice,:);
            temp_memoryPrecision_trialLevel_crossTime_lowMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaLow_choice,:);
        end
        
        % all trials
        temp_memoryPrecision_trialLevel_crossTime_choiceMemory = temp_memoryPrecision_trialLevel_crossTime(choiceMemoryBoolIndex,:);
        temp_memoryPrecision_trialLevel_crossTime_choiceOffload = temp_memoryPrecision_trialLevel_crossTime(choiceOffloadBoolIndex,:);
        
        % all trials with error
        temp_memoryPrecision_trialLevel_crossTime_choiceMemoryCorrect = temp_memoryPrecision_trialLevel_crossTime(choiceMemoryCorrectBoolIndex,:);
        temp_memoryPrecision_trialLevel_crossTime_choiceMemoryError = temp_memoryPrecision_trialLevel_crossTime(choiceMemoryErrorBoolIndex,:);
        
        a1 = temp_memoryPrecision_trialLevel_crossTime_highMatch;
        a2 = temp_memoryPrecision_trialLevel_crossTime_lowMatch;
        a3 = temp_memoryPrecision_trialLevel_crossTime_overMismatch;
        a4 = temp_memoryPrecision_trialLevel_crossTime_underMismatch;
        a5 = temp_memoryPrecision_trialLevel_crossTime_choiceMemory;
        a6 = temp_memoryPrecision_trialLevel_crossTime_choiceOffload;
        a7 = temp_memoryPrecision_trialLevel_crossTime_choiceMemoryCorrect;
        a8 = temp_memoryPrecision_trialLevel_crossTime_choiceMemoryError;
        
        a1_mean = mean(a1,1,'omitnan');
        a2_mean = mean(a2,1,'omitnan');
        a3_mean = mean(a3,1,'omitnan');
        a4_mean = mean(a4,1,'omitnan');
        a5_mean = mean(a5,1,'omitnan');
        a6_mean = mean(a6,1,'omitnan');
        a7_mean = mean(a7,1,'omitnan');
        a8_mean = mean(a8,1,'omitnan');
        
        a1_sem = std(a1,1,1,'omitnan')./sqrt(size(a1,1));
        a2_sem = std(a2,1,1,'omitnan')./sqrt(size(a2,1));
        a3_sem = std(a3,1,1,'omitnan')./sqrt(size(a3,1));
        a4_sem = std(a4,1,1,'omitnan')./sqrt(size(a4,1));
        a5_sem = std(a5,1,1,'omitnan')./sqrt(size(a5,1));
        a6_sem = std(a6,1,1,'omitnan')./sqrt(size(a6,1));
        a7_sem = std(a7,1,1,'omitnan')./sqrt(size(a7,1));
        a8_sem = std(a8,1,1,'omitnan')./sqrt(size(a8,1));
        
        range_crossTime = 1:size(temp_memoryPrecision_trialLevel_crossTime,2);
        x = range_crossTime;
        
        h_line = [];
        
        if if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 1
            y = a1_mean;
            h_line = [h_line plot(x,y,'color',color_choiceMemoryHigh,'linewidth',1)]; %#ok<*AGROW>
            hold on
            y_sem = a1_sem;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryHigh,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
            hold on
            
            y = a4_mean;
            h_line = [h_line plot(x,y,'color',color_choiceOffloadHigh,'linewidth',1)];
            hold on
            y_sem = a4_sem;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadHigh,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
            hold on
            
            
            y = a3_mean;
            h_line = [h_line plot(x,y,'color',color_choiceMemoryLow,'linewidth',1)];
            hold on
            y_sem = a3_sem;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
            hold on
            
            
            y = a2_mean;
            h_line = [h_line plot(x,y,'color',color_choiceOffloadLow,'linewidth',1)];
            hold on
            y_sem = a2_sem;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
            hold on
            
        elseif if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 0
            y = a5_mean;
            h_line = [h_line plot(x,y,'color',color_choiceMemory,'linewidth',1)]; %#ok<*AGROW>
            hold on
            y_sem = a5_sem;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
            hold on
            
            y = a6_mean;
            h_line = [h_line plot(x,y,'color',color_choiceOffload,'linewidth',1)];
            hold on
            y_sem = a6_sem;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
            hold on
        elseif if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 2
            y = a7_mean;
            h_line = [h_line plot(x,y,'color',color_choiceMemory,'linewidth',1)]; %#ok<*AGROW>
            hold on
            y_sem = a7_sem;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
            hold on
            
            
            y = a8_mean;
            h_line = [h_line plot(x,y,'color',color_choiceMemoryError,'linewidth',1)]; %#ok<*AGROW>
            hold on
            y_sem = a8_sem;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryError,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
            hold on
            
            
            y = a6_mean;
            h_line = [h_line plot(x,y,'color',color_choiceOffload,'linewidth',1)];
            hold on
            y_sem = a6_sem;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
            hold on
            
        end
        
        temptempTrialNum_a1 = sum(~isnan(a1(:,1)));
        temptempTrialNum_a2 = sum(~isnan(a2(:,1)));
        temptempTrialNum_a3 = sum(~isnan(a3(:,1)));
        temptempTrialNum_a4 = sum(~isnan(a4(:,1)));
        temptempTrialNum_a5 = sum(~isnan(a5(:,1)));
        temptempTrialNum_a6 = sum(~isnan(a6(:,1)));
        
        
        for tempi=1:(length(temp_interval)-1)
            plot(temp_interval(tempi)*[1 1],[y_min y_max],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        xticks(temp_interval(1:end-1));
        xticklabels({'Delay1','ChoiceCue'});
        xtickangle(0);
        
        yticks([0 0.5]);
        
        set(gca,'linewidth',1.5)
        xlim([range_crossTime(1) range_crossTime(end)]);
        %ylim([y_min-(y_max-y_min)*0.0 y_max+(y_max-y_min)*0]);%0.1
        %ylim([y_min-(y_max-y_min)*0.01 y_max-(y_max-y_min)*0.45]);
        ylim([0 y_max-(y_max-y_min)*0.45]);
        %ylim([0 1]);
        set(gca, 'FontSize', 10)
        set(gca,'color',backgounrdColor);
        set(gca,'box','off');% 取消右、上边框
        %xlabel('Frame', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('Quality', 'FontSize', 10, 'FontWeight', 'normal');
        
        %temp_title = title(sprintf('Population'),'FontSize',12,'FontWeight','normal');%'bold'
        %temp_title.Interpreter = 'none';
        
        
        
        %% bar plot
        nexttile
        
        memoryPrecision_trialLevel;
        %memoryMetaMismatch;
        trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh;
        trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh;
        trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow;
        trialBoolIndex_memoryPrecisionLow_choiceOffloadLow;
        
        temp_memoryPrecision_highMatch = memoryPrecision_trialLevel(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh);
        temp_memoryPrecision_OverMismatch = memoryPrecision_trialLevel(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh);
        temp_memoryPrecision_underMismatch = memoryPrecision_trialLevel(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow);
        temp_memoryPrecision_lowMatch = memoryPrecision_trialLevel(trialBoolIndex_memoryPrecisionLow_choiceOffloadLow);
        
        
        %         temp_1 = temp_memoryPrecision_OverMismatch;
        %         temp_2 = temp_memoryPrecision_highMatch;
        %         temp_3 = temp_memoryPrecision_lowMatch;
        %         temp_4 = temp_memoryPrecision_underMismatch;
        temp_1 = temp_memoryPrecision_highMatch;
        temp_2 = temp_memoryPrecision_underMismatch;
        temp_3 = temp_memoryPrecision_lowMatch;
        temp_4 = temp_memoryPrecision_OverMismatch;
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        temp3_SEM = std(temp_3)/sqrt(length(temp_3));
        temp4_SEM = std(temp_4)/sqrt(length(temp_4));
        
        temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM,mean(temp_3)-temp3_SEM,mean(temp_4)-temp4_SEM]);
        temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM,mean(temp_4)+temp4_SEM]);
        
        temp_y_max12 = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
        temp_y_max13 = max([mean(temp_1)+temp1_SEM,mean(temp_3)+temp3_SEM]);
        temp_y_max14 = max([mean(temp_1)+temp1_SEM,mean(temp_4)+temp4_SEM]);
        temp_y_max23 = max([mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM]);
        temp_y_max24 = max([mean(temp_2)+temp2_SEM,mean(temp_4)+temp4_SEM]);
        temp_y_max34 = max([mean(temp_3)+temp3_SEM,mean(temp_4)+temp4_SEM]);
        
        
        temp_bar = bar([0 1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3) mean(temp_4)], ...
            'FaceColor','flat');
        
        %         temp_bar.CData(1,:) = color_choiceMemoryLow;
        %         temp_bar.CData(2,:) = color_choiceMemoryHigh;
        %         temp_bar.CData(3,:) = color_choiceOffloadLow;
        %         temp_bar.CData(4,:) = color_choiceOffloadHigh;
        temp_bar.CData(1,:) = color_choiceMemoryHigh;
        temp_bar.CData(2,:) = color_choiceOffloadHigh;
        temp_bar.CData(3,:) = color_choiceOffloadLow;
        temp_bar.CData(4,:) = color_choiceMemoryLow;
        
        hold on
        
        %         errorbar([0], [mean(temp_1)],[temp1_SEM],'.','Color',color_choiceMemoryLow*0.65,...
        %             'LineWidth',2,'CapSize',6); %#ok<*NBRAK>
        %         hold on
        %         errorbar([1], [mean(temp_2)],[temp2_SEM],'.','Color',color_choiceMemoryHigh*0.65,...
        %             'LineWidth',2,'CapSize',6);
        %         hold on
        %         errorbar([2], [mean(temp_3)],[temp3_SEM],'.','Color',color_choiceOffloadLow*0.65,...
        %             'LineWidth',2,'CapSize',6);
        %         hold on
        %         errorbar([3], [mean(temp_4)],[temp4_SEM],'.','Color',color_choiceOffloadHigh*0.65,...
        %             'LineWidth',2,'CapSize',6);
        %         hold on
        errorbar([0], [mean(temp_1)],[temp1_SEM],'.','Color',color_choiceMemoryHigh*0.65,...
            'LineWidth',2,'CapSize',6); %#ok<*NBRAK>
        hold on
        errorbar([1], [mean(temp_2)],[temp2_SEM],'.','Color',color_choiceOffloadHigh*0.65,...
            'LineWidth',2,'CapSize',6);
        hold on
        errorbar([2], [mean(temp_3)],[temp3_SEM],'.','Color',color_choiceOffloadLow*0.65,...
            'LineWidth',2,'CapSize',6);
        hold on
        errorbar([3], [mean(temp_4)],[temp4_SEM],'.','Color',color_choiceMemoryLow*0.65,...
            'LineWidth',2,'CapSize',6);
        hold on
        
        
        yticks([0 0.5]);
        
        set(gca,'linewidth',1.5)
        ylim([0 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
        set(gca, 'FontSize', 10) %14
        set(gca,'XTickLabel','');
        set(gca,'box','off');% 取消右、上边框
        %ylabel('Memory precision', 'FontSize', 10, 'FontWeight', 'normal');
        
        
    end
    
end



%% Plot multiPanelB
if if_plot == 1
    if if_plot_multiPanelB == 1
        backgounrdColor = [1 1 1];
        
        for plot_lengthFlag=1:3
            temp_range = sum(numSeq(1:plot_lengthFlag-1))+1:sum(numSeq(1:plot_lengthFlag));
            temptempBoolIndex = ismember(seqIndex,temp_range)';
            
            if plot_lengthFlag==1
                trialBoolIndex_length1 = temptempBoolIndex;
            elseif plot_lengthFlag==2
                trialBoolIndex_length2 = temptempBoolIndex;
            elseif plot_lengthFlag==3
                trialBoolIndex_length3 = temptempBoolIndex;
            end
        end
        trialBoolIndex_length1;
        trialBoolIndex_length2;
        trialBoolIndex_length3;
        trialBoolIndex_length123 = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;

        %for plot_lengthFlag=1:3
        for plot_lengthFlag=1:7
            
            
            if plot_lengthFlag == 1
                memoryPrecision_trialLevel_crossTime_lengthx = memoryPrecision_trialLevel_crossTime_length1;
                lengthx_sample_interval = length1_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length1;
            elseif plot_lengthFlag == 2
                memoryPrecision_trialLevel_crossTime_lengthx = memoryPrecision_trialLevel_crossTime_length2;
                lengthx_sample_interval = length2_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length2;
            elseif plot_lengthFlag == 3
                memoryPrecision_trialLevel_crossTime_lengthx = memoryPrecision_trialLevel_crossTime_length3;
                lengthx_sample_interval = length3_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length3;
            elseif plot_lengthFlag == 4
                lengthx_sample_interval = length3_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
                
        
                temp_memoryPrecision_trialLevel_crossTime = nan(size(memoryPrecision_trialLevel_crossTime_length3));
                
                for temptempi=1:size(temp_memoryPrecision_trialLevel_crossTime,1)
                    if trialBoolIndex_length1(temptempi) == true
                        temp1 = memoryPrecision_trialLevel_crossTime_length1(temptempi,:);
                        temp_interp_factor = 3;
                        
                        if isnan(temp1(1)) == true
                            continue
                        end
                        
                        x = 1:length(temp1);
                        x_interp = (1:temp_interp_factor:length(temp1)*temp_interp_factor)+(temp_interp_factor-1);
                        
                        temp1_interp = interp1(x_interp,temp1,1:length(temp1)*temp_interp_factor,'makima');
                        
                        for temptempj=1:ceil((temp_interp_factor-1))
                            temp1_interp(temptempj) = temp1_interp(ceil(temp_interp_factor));
                        end
                        
                        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = temp1_interp;
                    end
                    if trialBoolIndex_length2(temptempi) == true
                        temp1 = memoryPrecision_trialLevel_crossTime_length2(temptempi,:);
                        temp_interp_factor = 1.5;
                        
                        if isnan(temp1(1)) == true
                            continue
                        end
                        
                        x = 1:length(temp1);
                        x_interp = (1:temp_interp_factor:length(temp1)*temp_interp_factor)+(temp_interp_factor-1);
                        
                        temp1_interp = interp1(x_interp,temp1,1:length(temp1)*temp_interp_factor,'makima');
                        
                        for temptempj=1:ceil((temp_interp_factor-1))
                            temp1_interp(temptempj) = temp1_interp(ceil(temp_interp_factor));
                        end
                        
                        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = temp1_interp;
                    end
                    if trialBoolIndex_length3(temptempi) == true
                        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = memoryPrecision_trialLevel_crossTime_length3(temptempi,:);
                    end
                    
                end
                temp_memoryPrecision_trialLevel_crossTime;
                
                memoryPrecision_trialLevel_crossTime_lengthx = temp_memoryPrecision_trialLevel_crossTime;
                
            elseif plot_lengthFlag == 5
                lengthx_sample_interval = length3_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
                
                memoryPrecision_trialLevel_crossTime_lengthx = memoryPrecision_trialLevel_crossTime_lengthx;%#ok<*ASGSL>                
                
            elseif plot_lengthFlag == 6
                lengthx_sample_interval = length1_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
                
        
                temp_memoryPrecision_trialLevel_crossTime = nan(size(memoryPrecision_trialLevel_crossTime_length1));
                
                for temptempi=1:size(temp_memoryPrecision_trialLevel_crossTime,1)
                    if trialBoolIndex_length1(temptempi) == true
                        temp1 = memoryPrecision_trialLevel_crossTime_length1(temptempi,:);
                        if isnan(temp1(1)) == true
                            continue
                        end                                                
                        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = temp1;
                    end
                    if trialBoolIndex_length2(temptempi) == true
                        temp1 = memoryPrecision_trialLevel_crossTime_length2(temptempi,:);
                        if isnan(temp1(1)) == true
                            continue
                        end
                        temp1 = temp1(end-lengthx_sample_interval(end)+1:end);
                        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = temp1;
                    end
                    if trialBoolIndex_length3(temptempi) == true
                        temp1 = memoryPrecision_trialLevel_crossTime_length3(temptempi,:);
                        if isnan(temp1(1)) == true
                            continue
                        end
                        temp1 = temp1(end-lengthx_sample_interval(end)+1:end);
                        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = temp1;
                    end
                    
                end
                temp_memoryPrecision_trialLevel_crossTime;
                
                memoryPrecision_trialLevel_crossTime_lengthx = temp_memoryPrecision_trialLevel_crossTime;
                             
            elseif plot_lengthFlag == 7
                lengthx_sample_interval = length1_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
                
                memoryPrecision_trialLevel_crossTime_lengthx = memoryPrecision_trialLevel_crossTime_lengthx;%#ok<*ASGSL>
                               
            end
            
            temp_trialBoolIndex_lengthx = temp_trialBoolIndex_lengthx & tempBoolIndex_4typesFilter;
            
            
            fig = figure('Name','asd','NumberTitle','off');
            %set(gcf,'Position',[10 50+200*(plot_lengthFlag-1) 360*0.975 200*0.8*0.8]);
            %set(gcf,'Position',[10 50+200*(plot_lengthFlag-1) 360*0.975*0.74 200*0.8*0.8]);  
            
            %temptemp1 = 360*0.975*0.74;
            temptemp1 = 360*0.975*0.74*1.37*0.97*1.02;
            %temptemp2 = 200*0.8*0.8;
            temptemp2 = 200*0.8*0.8*0.87;
            
            if if_dynamics_precision0_locCorr1 == 1
                temptemp1 = temptemp1 * 1.02;
            end
            
            if plot_lengthFlag <=4
                set(gcf,'Position',[10 50+200*(plot_lengthFlag-1) temptemp1 temptemp2]);
            elseif plot_lengthFlag >= 5
                %set(gcf,'Position',[10+400 50+200*(plot_lengthFlag-2) temptemp1 temptemp2]);
                set(gcf,'Position',[10+400*(plot_lengthFlag-4) 50+600 temptemp1 temptemp2]);
            end

            t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
            
            nexttile
            
            [y_min,y_max] = bounds([memoryPrecision_trialLevel_crossTime_baseline(:);...
                memoryPrecision_trialLevel_crossTime_length1(:);memoryPrecision_trialLevel_crossTime_length2(:);memoryPrecision_trialLevel_crossTime_length3(:);...
                memoryPrecision_trialLevel_crossTime_delay1(:)]);
            
            if if_dynamics_precision0_locCorr1 == 0
                if if_monkey_D0_Z1 == 0
                    if currentSessionIndex_AB == 4
                        y_max = 0.4;
                    end
                    if currentSessionIndex_AB == 7
                        y_max = 0.225;
                    end
                    if currentSessionIndex_AB == 3
                        y_max = 0.18;
                    end
                    if currentSessionIndex_AB == 9
                        y_max = 0.23;
                    end
                    if currentSessionIndex_AB == 11
                        y_max = 0.17;
                    end
                    if currentSessionIndex_AB == 12
                        y_max = 0.25;
                    end
                    
                elseif if_monkey_D0_Z1 == 1
                    y_max = 0.13;
                    
                    if currentSessionIndex_AB == 2
                        y_max = 0.13;
                    end
                    if currentSessionIndex_AB == 4
                        y_max = 0.15;
                    end
                    if currentSessionIndex_AB == 10
                        y_max = 0.18;
                    end
                    if currentSessionIndex_AB == 14
                        y_max = 0.20;
                    end
                    if currentSessionIndex_AB == 16
                        y_max = 0.19;
                    end
                    
                end
            end
            
            if if_dynamics_precision0_locCorr1 == 1
                if if_monkey_D0_Z1 == 0
                    
                    y_max = 1;
                    y_min =-0.5;
                    
                    if currentSessionIndex_AB == 8
                        y_max = 1.4;%2
                        y_min =-0.2;%-0.5
                    end                          
                    
                elseif if_monkey_D0_Z1 == 1
                    y_max = 1;
                    y_min =-0.5;
                end
                
            end
            
            %             temp1 = memoryPrecision_trialLevel_crossTime_delay1;
            %             temp1 = temp1(:,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
            %             [y_min,y_max] = bounds([memoryPrecision_trialLevel_crossTime_baseline(:);...
            %                 memoryPrecision_trialLevel_crossTime_lengthx(:);...
            %                 temp1(:)]);
            
            
            for loopCount=0:2
                
                if loopCount==0
                    temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_baseline;
                    temp_interval = baselinePeriod_interval;
                elseif loopCount==1
                    temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_lengthx;
                    temp_interval = baselinePeriod_interval(end) + lengthx_sample_interval;
                elseif loopCount==2
                    temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_delay1;
                    temp_interval = baselinePeriod_interval(end) + lengthx_sample_interval(end) + decisionPeriodA_interval;
                end
                
                
                % mismatch trials
                if if_plot_fineTuningMismatch == 1
                    temp_memoryPrecision_trialLevel_crossTime_overMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh & temp_trialBoolIndex_lengthx,:);
                    temp_memoryPrecision_trialLevel_crossTime_underMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow & temp_trialBoolIndex_lengthx,:);
                elseif if_plot_fineTuningMismatch == 0
                    temp_memoryPrecision_trialLevel_crossTime_overMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaHigh_choice & temp_trialBoolIndex_lengthx,:);
                    temp_memoryPrecision_trialLevel_crossTime_underMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaLow_choice & temp_trialBoolIndex_lengthx,:);
                end
                
                % match trials
                if if_plot_fineTuningMismatch == 1
                    temp_memoryPrecision_trialLevel_crossTime_highMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh & temp_trialBoolIndex_lengthx,:);
                    temp_memoryPrecision_trialLevel_crossTime_lowMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_choiceOffloadLow & temp_trialBoolIndex_lengthx,:);
                elseif if_plot_fineTuningMismatch == 0
                    temp_memoryPrecision_trialLevel_crossTime_highMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaHigh_choice & temp_trialBoolIndex_lengthx,:);
                    temp_memoryPrecision_trialLevel_crossTime_lowMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaLow_choice & temp_trialBoolIndex_lengthx,:);
                end
                
                % all trials
                temp_memoryPrecision_trialLevel_crossTime_choiceMemory = temp_memoryPrecision_trialLevel_crossTime(choiceMemoryBoolIndex' & temp_trialBoolIndex_lengthx,:);
                temp_memoryPrecision_trialLevel_crossTime_choiceOffload = temp_memoryPrecision_trialLevel_crossTime(choiceOffloadBoolIndex' & temp_trialBoolIndex_lengthx,:);
                
                % all trials with error
                temp_memoryPrecision_trialLevel_crossTime_choiceMemoryCorrect = temp_memoryPrecision_trialLevel_crossTime(choiceMemoryCorrectBoolIndex' & temp_trialBoolIndex_lengthx,:);
                temp_memoryPrecision_trialLevel_crossTime_choiceMemoryError = temp_memoryPrecision_trialLevel_crossTime(choiceMemoryErrorBoolIndex' & temp_trialBoolIndex_lengthx,:);
                
                a1 = temp_memoryPrecision_trialLevel_crossTime_highMatch;
                a2 = temp_memoryPrecision_trialLevel_crossTime_lowMatch;
                a3 = temp_memoryPrecision_trialLevel_crossTime_overMismatch;
                a4 = temp_memoryPrecision_trialLevel_crossTime_underMismatch;
                a5 = temp_memoryPrecision_trialLevel_crossTime_choiceMemory;
                a6 = temp_memoryPrecision_trialLevel_crossTime_choiceOffload;
                a7 = temp_memoryPrecision_trialLevel_crossTime_choiceMemoryCorrect;
                a8 = temp_memoryPrecision_trialLevel_crossTime_choiceMemoryError;
                
                a1_mean = mean(a1,1,'omitnan');
                a2_mean = mean(a2,1,'omitnan');
                a3_mean = mean(a3,1,'omitnan');
                a4_mean = mean(a4,1,'omitnan');
                a5_mean = mean(a5,1,'omitnan');
                a6_mean = mean(a6,1,'omitnan');
                a7_mean = mean(a7,1,'omitnan');
                a8_mean = mean(a8,1,'omitnan');
                
                a1_sem = std(a1,1,1,'omitnan')./sqrt(size(a1,1));
                a2_sem = std(a2,1,1,'omitnan')./sqrt(size(a2,1));
                a3_sem = std(a3,1,1,'omitnan')./sqrt(size(a3,1));
                a4_sem = std(a4,1,1,'omitnan')./sqrt(size(a4,1));
                a5_sem = std(a5,1,1,'omitnan')./sqrt(size(a5,1));
                a6_sem = std(a6,1,1,'omitnan')./sqrt(size(a6,1));
                a7_sem = std(a7,1,1,'omitnan')./sqrt(size(a7,1));
                a8_sem = std(a8,1,1,'omitnan')./sqrt(size(a8,1));
                
                [~,p_a5_a6] = ttest2(a5,a6);
                
                %if true
                if false
                    if loopCount == 0
                        diff_baseline = a7_mean-a6_mean;
                        diff_baseline_median = median(diff_baseline);
                    else
                        [~,p_a7_a6] = ttest2(a7,a6+diff_baseline_median);
                    end
                end
                
                if true
                %if false
                    if loopCount == 0
                        diff_baseline = a7_mean;
                        diff_baseline_median = median(diff_baseline);
                        %diff_baseline_median = mean(diff_baseline);
                        
                        diff_baseline_trial = a7;
                    else
                        p_a7_a6 = nan(1,size(a7,2));
                        for temptempi=1:size(a7,2)
                            %[~,p_a7_a6(temptempi)] = ttest2(a7(:,temptempi),diff_baseline');
                            %[~,p_a7_a6(temptempi)] = ttest(a7(:,temptempi),diff_baseline_median);
                            %[~,p_a7_a6(temptempi)] = ttest(a7(:,temptempi),diff_baseline_median,'tail','right');
                            %[~,p_a7_a6(temptempi)] = ttest(a7(:,temptempi),0,'tail','right');
                            [~,p_a7_a6(temptempi)] = ttest2(diff_baseline,mean(a7(:,temptempi),'omitnan'),'tail','left');                            
                            %[~,p_a7_a6(temptempi)] = ttest2(diff_baseline_trial(:),mean(a7(:,temptempi),'omitnan'),'tail','left');
                                                        
                        end
                    end                    
                end
                
                temp_memoryPrecision_seqLevel_x_crossTime = zeros(sum(numSeq(1:3)),size(temp_memoryPrecision_trialLevel_crossTime,2));
                for temptempi=1:sum(numSeq(1:3))
                    temptempBoolIndex = temp_seqIndex==temptempi;
                    %temptempBoolIndex_x = temptempBoolIndex;
                    temptempBoolIndex_x = temptempBoolIndex & trialIndex_bool_memoryCorrect;
                    %temptempBoolIndex_x = temptempBoolIndex & choiceMemoryCorrectBoolIndex;   
                    %temptempBoolIndex_x = temptempBoolIndex & (trialIndex_bool_memoryCorrect&~choiceBoolIndex);  
                    temp_memoryPrecision_seqLevel_x_crossTime(temptempi,:) = mean(temp_memoryPrecision_trialLevel_crossTime(temptempBoolIndex_x,:),'omitnan');
                end
                
                temp_x = temp_memoryPrecision_seqLevel_x_crossTime;
                temp_y = seqPrecision_behavior(1:sum(numSeq(1:3)));
                
                temp_r = nan(1,size(temp_memoryPrecision_seqLevel_x_crossTime,2));
                temp_p = nan(1,size(temp_memoryPrecision_seqLevel_x_crossTime,2));
                for temptempi=1:size(temp_memoryPrecision_seqLevel_x_crossTime,2)
                    temptemp_x = temp_x(:,temptempi);
                    [temp_r(temptempi),temp_p(temptempi)] = corr(temptemp_x(~isnan(temptemp_x)),temp_y(~isnan(temptemp_x)));
                end
                
                p_a9 = temp_p;  
                
                a9 = temp_r;
                if loopCount == 0
                    a9_baseline = a9;
                    a9_baseline_median = median(a9);
                end
                    
                range_crossTime = (1:size(temp_memoryPrecision_trialLevel_crossTime,2))+temp_interval(1)-1;
                x = range_crossTime;
                
                h_line = [];
                
                if plot_lengthFlag <= 4 || plot_lengthFlag == 6
                    if if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 1
                        y = a3_mean;
                        h_line = [h_line plot(x,y,'color',color_choiceMemoryLow,'linewidth',1)];
                        hold on
                        y_sem = a3_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        y = a1_mean;
                        h_line = [h_line plot(x,y,'color',color_choiceMemoryHigh,'linewidth',1)]; %#ok<*AGROW>
                        hold on
                        y_sem = a1_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryHigh,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        y = a2_mean;
                        h_line = [h_line plot(x,y,'color',color_choiceOffloadLow,'linewidth',1)];
                        hold on
                        y_sem = a2_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        y = a4_mean;
                        h_line = [h_line plot(x,y,'color',color_choiceOffloadHigh,'linewidth',1)];
                        hold on
                        y_sem = a4_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadHigh,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                    elseif if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 0
                        y = a5_mean;
                        h_line = [h_line plot(x,y,'color',color_choiceMemory,'linewidth',1)]; %#ok<*AGROW>
                        hold on
                        y_sem = a5_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        y = a6_mean;
                        h_line = [h_line plot(x,y,'color',color_choiceOffload,'linewidth',1)];
                        hold on
                        y_sem = a6_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        p_a5_a6;
                        
                        if if_ttest_crossTime == 1
                            scatter(x(p_a5_a6<p_threshold),y_max-(y_max-y_min)*0.15,8,[0 0 0],'*');
                        end
                        
                    elseif if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 2
                        y = a7_mean;
                        h_line = [h_line plot(x,y,'color',color_choiceMemory,'linewidth',1)]; %#ok<*AGROW>
                        hold on
                        y_sem = a7_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        y = a6_mean;
                        h_line = [h_line plot(x,y,'color',color_choiceOffload,'linewidth',1)];
                        hold on
                        y_sem = a6_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        y = a8_mean;
                        h_line = [h_line plot(x,y,'color',color_choiceMemoryError,'linewidth',1)]; %#ok<*AGROW>
                        hold on
                        y_sem = a8_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryError,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        
                        if if_ttest_crossTime == 1
                            
                            if plot_lengthFlag <= 3
                                if loopCount >= 1
                                    scatter(x(p_a7_a6<p_threshold),y_max-(y_max-y_min)*0.42,8,[0 0 0],'*');
                                end
                            end
                            if plot_lengthFlag >= 4
                                if loopCount >= 1
                                    scatter(x(p_a7_a6<p_threshold),y_max-(y_max-y_min)*0.59,8,[0 0 0],'*');
                                end
                            end
                            
                        end
                    
                        elseif if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 3
                            y = a7_mean;
                            h_line = [h_line plot(x,y,'color',[1 1 1]*0.3,'linewidth',1)]; %#ok<*AGROW>
                            hold on
                            
                            if if_ttest_crossTime == 1
                                
                                if plot_lengthFlag <= 3
                                    if loopCount >= 1
                                        scatter(x(p_a7_a6<p_threshold),y_max-(y_max-y_min)*0.42,8,[0 0 0],'*');
                                    end
                                end
                                if plot_lengthFlag >= 4
                                    if loopCount >= 1
                                        %scatter(x(p_a7_a6<p_threshold),y_max-(y_max-y_min)*0.59,8,[0 0 0],'*');
                                        if plot_lengthFlag ~= 6
                                            scatter(x(p_a7_a6<p_threshold),y_max-(y_max-y_min)*0.575,8,[0 0 0],'*');
                                        elseif plot_lengthFlag == 6
                                            if if_dynamics_precision0_locCorr1 == 0
                                                scatter(x(p_a7_a6<p_threshold),y_max-(y_max-y_min)*0.555,8,[0 0 0],'*');
                                            elseif if_dynamics_precision0_locCorr1 == 1
                                                scatter(x(p_a7_a6<p_threshold),y_max-(y_max-y_min)*0.355,8,[0 0 0],'*');
                                            end
                                        end
                                    end
                                end
                                
                            end
                        
                    end                    
                    
                elseif plot_lengthFlag == 5 || plot_lengthFlag == 7
                    y_max = 1;
                    y_min =-0.5;
                    
                    y = a9;
                    h_line = [h_line plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1)]; %#ok<*AGROW>
                    hold on
                    
                    if loopCount == 2
                        plot([1 x(end)],[1 1]*a9_baseline_median,':','color',[0.3 0.3 0.3],'linewidth',1);
                        hold on
                    end
                    
                    if if_ttest_crossTime == 1
                        if loopCount >= 1
                            scatter(x(p_a9<p_threshold),y_max+(y_max-y_min)*0.0,8,[0 0 0],'*');
                        end
                    end
                    
                end
                
                if loopCount==1 && plot_lengthFlag<=3
                    for tempi=1:(length(temp_interval)-1)
                        x = [temp_interval(1+tempi-1) y_min;...
                            temp_interval(1+tempi-1) y_max;...
                            temp_interval(1+tempi-1)+6 y_max;...
                            temp_interval(1+tempi-1)+6 y_min];
                        y = [1 2 3 4];
                        patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.1,'EdgeColor','none');
                        hold on
                    end
                end
                
            end
            
            %if false
                if if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 2
                    if plot_lengthFlag == 1 || plot_lengthFlag == 4
                        le = legend(h_line,'Choice-memory correct','Choice-offload','Choice-memory error',...
                            'Location','northwest','fontsize',6.5);%10
                        le.ItemTokenSize = ones(1,3)*10;
                        le.Color = backgounrdColor;
                    end
                end
            %end
            
            
            temp_interval_all = [baselinePeriod_interval(1) ...
                (baselinePeriod_interval(end) + lengthx_sample_interval(1:end-1)) ...
                (baselinePeriod_interval(end) + lengthx_sample_interval(end) + decisionPeriodA_interval)];
            
            xticks(temp_interval_all(1:end-2));
            if plot_lengthFlag == 1
                xticklabels({'Fixation','T1','Delay-on'});
            elseif plot_lengthFlag == 2
                xticklabels({'Fixation','T1','T2','Delay-on'});
            elseif plot_lengthFlag == 3
                xticklabels({'Fixation','T1','T2','T3','Delay-on'});
            elseif plot_lengthFlag == 4
                xticks(temp_interval_all([1 2 5]));
                xticklabels({'Fixation','Sample','Delay-on'});
            elseif plot_lengthFlag == 5
                xticks(temp_interval_all([1 2 5]));
                xticklabels({'Fixation','Sample','Delay-on'});
            elseif plot_lengthFlag == 6
                xticks(temp_interval_all([1 2 3]));
                xticklabels({'Fixation','Last T','Delay-on'});
            elseif plot_lengthFlag == 7
                xticks(temp_interval_all([1 2 3]));
                xticklabels({'Fixation','Last T','Delay-on'});
            end
            
            xlim([temp_interval_all(1) temp_interval_all(end-1)]);
            
            
            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 8)%12
            set(gca,'color',backgounrdColor);
            set(gca,'box','off');% 取消右、上边框
            xtickangle(0);
            
            if plot_lengthFlag <= 4 || plot_lengthFlag == 6
                ylim([y_min+(y_max-y_min)*0.0 y_max-(y_max-y_min)*0.38]);
                if plot_lengthFlag == 4 || plot_lengthFlag == 6
                    ylim([y_min+(y_max-y_min)*0.0 y_max-(y_max-y_min)*0.565]);
                    %ylim([y_min+(y_max-y_min)*0.0 y_max-(y_max-y_min)*0.545]);
                end
                
                ylabel('Quality', 'FontSize', 9, 'FontWeight', 'normal');%10
                
                if plot_lengthFlag == 6
                    if if_dynamics_precision0_locCorr1 == 0                        
                        ylim([y_min+(y_max-y_min)*0.0 y_max-(y_max-y_min)*0.525]);
                    elseif if_dynamics_precision0_locCorr1 == 1
                        ylim([y_min+(y_max-y_min)*0.0 y_max-(y_max-y_min)*0.325]);
                        ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');%10
                    end
                end
                
                
                
                if if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 3
                    title('Memory quality time course','FontSize',9);
                    if if_dynamics_precision0_locCorr1 == 0
                        subtitle('Correct trials','FontSize',7.5);
                    elseif if_dynamics_precision0_locCorr1 == 1
                       subtitle('Correct trials, location correlation','FontSize',7.5); 
                    end
                end
            elseif plot_lengthFlag == 5 || plot_lengthFlag == 7
                ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);
                
                %set(gca,'ydir','reverse');
                ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');%10
                title('Memory quality VS. accuracy','FontSize',9);
            end
            
            
        end
        
    end
    
end


%% Plot multiPanelC
if if_plot == 1
    if if_plot_multiPanelC == 1
        backgounrdColor = [1 1 1];
        
        for plot_lengthFlag=1:3
            temp_range = sum(numSeq(1:plot_lengthFlag-1))+1:sum(numSeq(1:plot_lengthFlag));
            temptempBoolIndex = ismember(seqIndex,temp_range)';
            
            if plot_lengthFlag==1
                trialBoolIndex_length1 = temptempBoolIndex;
            elseif plot_lengthFlag==2
                trialBoolIndex_length2 = temptempBoolIndex;
            elseif plot_lengthFlag==3
                trialBoolIndex_length3 = temptempBoolIndex;
            end
        end
        trialBoolIndex_length1;
        trialBoolIndex_length2;
        trialBoolIndex_length3;
        trialBoolIndex_length123 = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;

        for plot_lengthFlag=1:3
                        
            if plot_lengthFlag == 1
                memoryPrecision_trialLevel_crossTime_lengthx = memoryPrecision_trialLevel_crossTime_length1;
                lengthx_sample_interval = length1_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length1;
            elseif plot_lengthFlag == 2
                memoryPrecision_trialLevel_crossTime_lengthx = memoryPrecision_trialLevel_crossTime_length2;
                lengthx_sample_interval = length2_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length2;
            elseif plot_lengthFlag == 3
                memoryPrecision_trialLevel_crossTime_lengthx = memoryPrecision_trialLevel_crossTime_length3;
                lengthx_sample_interval = length3_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length3;
            end
            
            temp_trialBoolIndex_lengthx = temp_trialBoolIndex_lengthx & tempBoolIndex_4typesFilter;
            
            
            fig = figure('Name','asd','NumberTitle','off');
            
            temptemp1 = 360*0.975*0.74*1.37*0.97*1.02;
            temptemp2 = 200*0.8*0.8;
            
            set(gcf,'Position',[410 50+200*(plot_lengthFlag-1) temptemp1 temptemp2]);

            t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
            
            nexttile
                        
            y_min = -0.5;%-1
            y_max = 1;
            
            for loopCount=0:2
                
                if loopCount==0
                    temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_baseline;
                    temp_interval = baselinePeriod_interval;
                elseif loopCount==1
                    temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_lengthx;
                    temp_interval = baselinePeriod_interval(end) + lengthx_sample_interval;
                elseif loopCount==2
                    temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_delay1;
                    temp_interval = baselinePeriod_interval(end) + lengthx_sample_interval(end) + decisionPeriodA_interval;
                end
                

                temp_memoryPrecision_seqLevel_x_crossTime = zeros(sum(numSeq(1:3)),size(temp_memoryPrecision_trialLevel_crossTime,2));
                for temptempi=1:sum(numSeq(1:3))
                    temptempBoolIndex = temp_seqIndex==temptempi;
                    %temptempBoolIndex_x = temptempBoolIndex;
                    temptempBoolIndex_x = temptempBoolIndex & trialIndex_bool_memoryCorrect;
                    %temptempBoolIndex_x = temptempBoolIndex & choiceMemoryCorrectBoolIndex;   
                    %temptempBoolIndex_x = temptempBoolIndex & (trialIndex_bool_memoryCorrect&~choiceBoolIndex);  
                    temp_memoryPrecision_seqLevel_x_crossTime(temptempi,:) = mean(temp_memoryPrecision_trialLevel_crossTime(temptempBoolIndex_x,:),'omitnan');
                end
                
                temp_range = sum(numSeq(1:plot_lengthFlag-1))+1:sum(numSeq(1:plot_lengthFlag));
                
                temptempBoolIndex = true(size(temp_memoryPrecision_seqLevel_x_crossTime,1),1);
                temptempBoolIndex(temp_range) = false;
                
                temp_memoryPrecision_seqLevel_x_crossTime(temptempBoolIndex,:) = nan;
                
                temp_x = temp_memoryPrecision_seqLevel_x_crossTime;
                temp_y = seqPrecision_behavior(1:sum(numSeq(1:3)));
                
                temp_r = nan(1,size(temp_memoryPrecision_seqLevel_x_crossTime,2));
                temp_p = nan(1,size(temp_memoryPrecision_seqLevel_x_crossTime,2));
                for temptempi=1:size(temp_memoryPrecision_seqLevel_x_crossTime,2)
                    temptemp_x = temp_x(:,temptempi);
                    [temp_r(temptempi),temp_p(temptempi)] = corr(temptemp_x(~isnan(temptemp_x)),temp_y(~isnan(temptemp_x)));
                end
                
                p_a9 = temp_p;
                
                a9 = temp_r;
                if loopCount == 0
                    a9_baseline = a9;
                    a9_baseline_median = median(a9);
                end
                    
                range_crossTime = (1:size(temp_memoryPrecision_trialLevel_crossTime,2))+temp_interval(1)-1;
                x = range_crossTime;
                
                h_line = [];
                
                y = a9;
                h_line = [h_line plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1)]; %#ok<*AGROW>
                hold on
                
                if loopCount == 2
                    plot([1 x(end)],[1 1]*a9_baseline_median,':','color',[0.3 0.3 0.3],'linewidth',1);
                    hold on
                end
                
                if loopCount >= 1
                    scatter(x(p_a9<p_threshold),y_max+(y_max-y_min)*0.0,8,[0 0 0],'*');
                end                
                
                    
                if loopCount==1 && plot_lengthFlag<=3
                    for tempi=1:(length(temp_interval)-1)
                        x = [temp_interval(1+tempi-1) y_min;...
                            temp_interval(1+tempi-1) y_max;...
                            temp_interval(1+tempi-1)+6 y_max;...
                            temp_interval(1+tempi-1)+6 y_min];
                        y = [1 2 3 4];
                        patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.1,'EdgeColor','none');
                        hold on
                    end
                end
                
            end            
            
            temp_interval_all = [baselinePeriod_interval(1) ...
                (baselinePeriod_interval(end) + lengthx_sample_interval(1:end-1)) ...
                (baselinePeriod_interval(end) + lengthx_sample_interval(end) + decisionPeriodA_interval)];
            
            xticks(temp_interval_all(1:end-2));
            if plot_lengthFlag == 1
                xticklabels({'Fixation','T1','Delay1'});
            elseif plot_lengthFlag == 2
                xticklabels({'Fixation','T1','T2','Delay1'});
            elseif plot_lengthFlag == 3
                xticklabels({'Fixation','T1','T2','T3','Delay1'});              
            end
            
            xlim([temp_interval_all(1) temp_interval_all(end-1)]);
            
            %ylim([-0.7 0.7]);
            ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);
            
            
            ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');%10
            title('Memory quality VS. accuracy','FontSize',10);
            
            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 8)%12
            set(gca,'color',backgounrdColor);
            set(gca,'box','off');% 取消右、上边框
            xtickangle(0);
            
            
        end
        
    end
    
end


if if_profileOn == 1
    profile viewer
end

%% End