% Chuan's 13th script (20251214)
% This script: Time courses of meta-memory, related to figure 4.
%% Initialization
close all

if_plot = 1;%0

if_plot_multiPanel = 0;
if_plot_singlePanel = 0;

if_plot_multiPanelB = 0; % 1, choice level
if_plot_multiPanelC = 0; % sequence level correlation
if_plot_multiPanelD = 0; % decoding accuracy

if_compute = 1;%1

if_smooth = 0;

if_StatQuant = 1;


if_delay1_forward0_backward1 = 0; % new

% if if_monkey_D0_Z1 == 1
%     if_delay1_forward0_backward1 = 1;
% end

if exist('if_compute_summary','var') == 1
    if_delay1_forward0_backward1 = 1;
end

if_resample_meta = if_resample_meta;
resampleIterCount_meta = resampleIterCount_meta;

if_ttest_crossTime = 1;%1


if_allTrial0_matchMismatch1_allTrialWithError2_correct3 = 2;%3

if_trainPeriod_eachBin0_delay1Bin1 = 1;%0

if if_resample_meta == 1
    if_trainPeriod_eachBin0_delay1Bin1 = 1;
end

if_plot_fineTuningMismatch = 0;


p_threshold = 0.01;%0.05,0.001,0.01


meta_trialLevel;
meta_trialLevel;


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

%% To categorize match and mismatch trials
% only neuron label
if true
    trialBoolIndex_memoryPrecisionLow_metaLow_choice;
    trialBoolIndex_memoryPrecisionHigh_metaHigh_choice;
    trialBoolIndex_memoryPrecisionLow_metaHigh_choice;
    trialBoolIndex_memoryPrecisionHigh_metaLow_choice;
    
    b1 = sum(trialBoolIndex_memoryPrecisionLow_metaLow_choice);
    b2 = sum(trialBoolIndex_memoryPrecisionHigh_metaHigh_choice);
    b3 = sum(trialBoolIndex_memoryPrecisionLow_metaHigh_choice);
    b4 = sum(trialBoolIndex_memoryPrecisionHigh_metaLow_choice);
end


%% Test meta
if if_compute == 1
    fun_testMeta_currentTime_Name_v = autoGetFunName_myScripts('fun_testMeta_currentTime', [targetPATH '\functions']);
    fun_testMeta_currentTime = str2func(fun_testMeta_currentTime_Name_v);
    
    options_testMeta = struct;
    options_testMeta.KFold_num = KFold_num_meta;
    options_testMeta.choiceBoolIndex = choiceBoolIndex_validLength;
    options_testMeta.noChoiceBoolIndex = noChoiceBoolIndex_validLength;
    
    F_dff_baselinePeriod;
    F_dff_length1_sample;
    F_dff_length2_sample;
    F_dff_length3_sample;
    F_dff_decisionPeriodA;
    
    fun_metaDecoder_Name_v = autoGetFunName_myScripts('metaDecoder', [targetPATH '\functions']);
    fun_metaDecoder = str2func(fun_metaDecoder_Name_v);
    t_decoder = templateSVM('Standardize',true,'KernelFunction','linear'); % To standardise input data
    
    if if_trainPeriod_eachBin0_delay1Bin1 == 1
        options_testMeta.svm_Meta = svm_Meta;
    end
    
    for loopCount=0:4
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
                F_dff_current = smoothdata(F_dff_current,3,'gaussian',20);
            end
        end         
        
        if if_trainPeriod_eachBin0_delay1Bin1 == 0
            if loopCount<4
                F_dff_current_bin = mean(F_dff_current,3);
            else
                F_dff_current_bin = mean(F_dff_current(:,:,1:33),3);
            end
            
            svm_options = struct;
            svm_options.numSeq = numSeq;
            svm_options.seqIndex_choice = seqIndex_choice;
            svm_options.KFold_num = KFold_num;
            svm_options.t_decoder = t_decoder;
            
            x_raw = F_dff_current_bin';
            x = x_raw(choiceBoolIndex,:);
            y_raw = choiceMemoryBoolIndex';
            y = y_raw(choiceBoolIndex);
            
            svm_Meta = fun_metaDecoder(x,y,svm_options);
            options_testMeta.svm_Meta = svm_Meta;
        end
        
        
        temptemp_range = 1:size(F_dff_current,3);
        
        if if_resample_meta == 0
            temp_meta_trialLevel_crossTime = zeros(length(seqIndex),length(temptemp_range));
            parfor tempi=temptemp_range(1):temptemp_range(end)
                F_dff_currentTime = F_dff_current(:,:,temptemp_range(tempi)); %#ok<*PFBNS>
                testMeta_output = fun_testMeta_currentTime(F_dff_currentTime,options_testMeta);
                temp_meta_trialLevel_crossTime(:,tempi) = testMeta_output.meta_trialLevel_currentTime;
            end
            
        elseif if_resample_meta == 1
            
            a1 = validChoiceMemoryTrialIndex_resampled;
            a2 = validChoiceOffloadTrialIndex_resampled;
            validTrialIndex_resampled = [a1,a2];
            
            
            x_raw_time = F_dff_current;
            
            timeBinSize = 1;
            timeBinNum = floor(size(x_raw_time,3)/timeBinSize);
            
            Posterior_2d_time_raw = nan(size(validTrialIndex_resampled,2),timeBinNum,resampleIterCount_meta);
            
            parfor tempIter=1:resampleIterCount_meta
                
                a1 = [];
                
                temp_svm_Meta = svm_Meta.temp_svm_Meta_resampled{tempIter};
                
                temp_validTrialIndex = validTrialIndex_resampled(tempIter,:);
                
                for tempi=1:timeBinNum
                    
                    temp_range = ((tempi-1)*timeBinSize+1):tempi*timeBinSize;
                    
                    if timeBinSize > 1
                        a1 = mean(x_raw_time(:,:,temp_range),3);
                    else
                        a1 = x_raw_time(:,:,temp_range);
                    end
                    
                    temp_x_raw = a1';
                    
                    %%  Test meta-memory on training set trials (choice trials)
                    multi_Posterior_cell = cell(1, KFold_num_meta);
                    temp_Mdl_CV_binary = temp_svm_Meta.temp_Mdl_CV_binary;
                    temp_x = temp_x_raw(temp_validTrialIndex,:);
                    
                    F_dff_T = temp_x;
                    
                    for tempk=1:KFold_num_meta
                        tempTrialBoolIndex_fold = temp_Mdl_CV_binary.ModelParameters.Generator.UseObsForIter(:,tempk);
                        temp_F_dff_T = F_dff_T(~tempTrialBoolIndex_fold,:); %#ok<*PFBNS>
                        temp_F_dff_T_2d = temp_F_dff_T;
                        [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary.Trained{tempk},temp_F_dff_T_2d);
                        if size(tempPosterior,2) == 1
                            tempPosterior_2 = tempPosterior(:,1);
                        else
                            tempPosterior_2 = tempPosterior(:,2);
                        end
                        multi_Posterior_cell{tempk} = tempPosterior_2;
                    end
                    
                    temp_Posterior_2d = zeros(size(F_dff_T,1),1);
                    for tempk=1:KFold_num_meta
                        temp_Posterior = multi_Posterior_cell{tempk};
                        tempTrialBoolIndex_fold = temp_Mdl_CV_binary.ModelParameters.Generator.UseObsForIter(:,tempk);
                        temp_Posterior_2d(~tempTrialBoolIndex_fold) = temp_Posterior;
                    end
                    temp_Posterior_2d;
                    
                    Posterior_2d_time_raw(:,tempi,tempIter) = temp_Posterior_2d;
                    
                end
                
            end
            
            
            Posterior_2d_time_raw;
            Posterior_2d_time = nan(length(choiceBoolIndex),timeBinNum,resampleIterCount_meta);
            
            
            a1 = validChoiceMemoryTrialIndex_resampled;
            a2 = validChoiceOffloadTrialIndex_resampled;
            a12 = [a1,a2];
            
            for tempi=1:resampleIterCount_meta
                temp_posterior = Posterior_2d_time_raw(:,:,tempi);
                
                temp_trialIndex = a12(tempi,:);
                
                Posterior_2d_time(temp_trialIndex,:,tempi) = temp_posterior;
            end
            
            Posterior_2d_time_mean = mean(Posterior_2d_time,3,'omitnan');
            
            temp_meta_trialLevel_crossTime = Posterior_2d_time_mean;
            
            Posterior_2d_time_T = permute(Posterior_2d_time,[1 3 2]);
            
        end
        
        if loopCount==0
            meta_trialLevel_crossTime_baseline = temp_meta_trialLevel_crossTime;            
            meta_trialLevel_resampleIter_crossTime_baseline = Posterior_2d_time_T;
            
        elseif loopCount==1
            meta_trialLevel_crossTime_length1 = temp_meta_trialLevel_crossTime;
            meta_trialLevel_resampleIter_crossTime_length1 = Posterior_2d_time_T;
            
        elseif loopCount==2
            meta_trialLevel_crossTime_length2 = temp_meta_trialLevel_crossTime;
            meta_trialLevel_resampleIter_crossTime_length2 = Posterior_2d_time_T;
            
        elseif loopCount==3
            meta_trialLevel_crossTime_length3 = temp_meta_trialLevel_crossTime;
            meta_trialLevel_resampleIter_crossTime_length3 = Posterior_2d_time_T;
            
        elseif loopCount==4
            meta_trialLevel_crossTime_delay1 = temp_meta_trialLevel_crossTime;
            meta_trialLevel_resampleIter_crossTime_delay1 = Posterior_2d_time_T;
            
        end
        
    end
end
meta_trialLevel_crossTime_baseline;
meta_trialLevel_crossTime_length1;
meta_trialLevel_crossTime_length2;
meta_trialLevel_crossTime_length3;
meta_trialLevel_crossTime_delay1;

meta_trialLevel_resampleIter_crossTime_baseline;
meta_trialLevel_resampleIter_crossTime_length1;
meta_trialLevel_resampleIter_crossTime_length2;
meta_trialLevel_resampleIter_crossTime_length3;
meta_trialLevel_resampleIter_crossTime_delay1;

%% Meta-memory decoding accuracy
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

y_raw = nan((length(choiceMemoryBoolIndex)),1);
y_raw(choiceMemoryBoolIndex_validLength) = true;
y_raw(choiceOffloadBoolIndex_validLength) = false;

% for plot_lengthFlag=1:4
%     if plot_lengthFlag == 1
%         temp_trialBoolIndex_lengthx = trialBoolIndex_length1;
%     elseif plot_lengthFlag == 2
%         temp_trialBoolIndex_lengthx = trialBoolIndex_length2;
%     elseif plot_lengthFlag == 3
%         temp_trialBoolIndex_lengthx =  trialBoolIndex_length3;
%     elseif plot_lengthFlag == 4
%         temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
%     end
%
%     for tempi=1:5
%         if tempi == 1
%             x_raw = meta_trialLevel_crossTime_baseline;
%         elseif tempi == 2
%             x_raw = meta_trialLevel_crossTime_length1;
%         elseif tempi == 3
%             x_raw = meta_trialLevel_crossTime_length2;
%         elseif tempi == 4
%             x_raw = meta_trialLevel_crossTime_length3;
%         elseif tempi == 5
%             x_raw = meta_trialLevel_crossTime_delay1;
%         end
%
%         temptempBoolIndex = (~isnan(x_raw(:,1))) & temp_trialBoolIndex_lengthx;
%
%         x = x_raw(temptempBoolIndex,:);
%
%         y = y_raw(temptempBoolIndex);
%         y = y==true;
%
%         y_predict = x>=metaDecoderThreshold_delay1;
%
%         temp_meta_crossTime_decodingAcc = nan(1,size(x,2));
%         for tempj=1:size(x,2)
%             temp_meta_crossTime_decodingAcc(tempj) = sum(y==y_predict(:,tempj))/length(y);
%         end
%
%         if plot_lengthFlag == 1
%             if tempi == 1
%                 meta_crossTime_decodingAcc_baseline = temp_meta_crossTime_decodingAcc;
%             elseif tempi == 2
%                 meta_crossTime_decodingAcc_length1 = temp_meta_crossTime_decodingAcc;
%             elseif tempi == 3
%                 meta_crossTime_decodingAcc_length2 = temp_meta_crossTime_decodingAcc;
%             elseif tempi == 4
%                 meta_crossTime_decodingAcc_length3 = temp_meta_crossTime_decodingAcc;
%             elseif tempi == 5
%                 meta_crossTime_decodingAcc_delay1 = temp_meta_crossTime_decodingAcc;
%             end
%         elseif plot_lengthFlag == 2
%
%         elseif plot_lengthFlag == 3
%
%         elseif plot_lengthFlag == 4
%
%         end
%     end
% end


for plot_lengthFlag=1:4
    if plot_lengthFlag == 1
        temp_trialBoolIndex_lengthx = trialBoolIndex_length1;
    elseif plot_lengthFlag == 2
        temp_trialBoolIndex_lengthx = trialBoolIndex_length2;
    elseif plot_lengthFlag == 3
        temp_trialBoolIndex_lengthx =  trialBoolIndex_length3;
    elseif plot_lengthFlag == 4
        temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
    end
    
    for tempi=1:2
        if tempi == 1
            x_raw = meta_trialLevel_crossTime_baseline;
        elseif tempi == 2
            x_raw = meta_trialLevel_crossTime_delay1;
        end
        
        temptempBoolIndex = (~isnan(x_raw(:,1))) & temp_trialBoolIndex_lengthx;
        
        x = x_raw(temptempBoolIndex,:);
        
        y = y_raw(temptempBoolIndex);
        y = y==true;
        
        y_predict = x>=metaDecoderThreshold_delay1;
        
        temp_meta_crossTime_decodingAcc = nan(1,size(x,2));
        for tempj=1:size(x,2)
            temp_meta_crossTime_decodingAcc(tempj) = sum(y==y_predict(:,tempj))/length(y);
        end
        
        if plot_lengthFlag == 1
            if tempi == 1
                meta_crossTime_decodingAcc_baseline_length1 = temp_meta_crossTime_decodingAcc;
            elseif tempi == 2
                meta_crossTime_decodingAcc_delay1_length1 = temp_meta_crossTime_decodingAcc;
            end
        elseif plot_lengthFlag == 2
            if tempi == 1
                meta_crossTime_decodingAcc_baseline_length2 = temp_meta_crossTime_decodingAcc;
            elseif tempi == 2
                meta_crossTime_decodingAcc_delay1_length2 = temp_meta_crossTime_decodingAcc;
            end
        elseif plot_lengthFlag == 3
            if tempi == 1
                meta_crossTime_decodingAcc_baseline_length3 = temp_meta_crossTime_decodingAcc;
            elseif tempi == 2
                meta_crossTime_decodingAcc_delay1_length3 = temp_meta_crossTime_decodingAcc;
            end
        elseif plot_lengthFlag == 4
            if tempi == 1
                meta_crossTime_decodingAcc_baseline = temp_meta_crossTime_decodingAcc;
            elseif tempi == 2
                meta_crossTime_decodingAcc_delay1 = temp_meta_crossTime_decodingAcc;
            end
        end
    end
end
meta_crossTime_decodingAcc_baseline_length1;
meta_crossTime_decodingAcc_delay1_length1;
meta_crossTime_decodingAcc_baseline_length2;
meta_crossTime_decodingAcc_delay1_length2;
meta_crossTime_decodingAcc_baseline_length3;
meta_crossTime_decodingAcc_delay1_length3;
meta_crossTime_decodingAcc_baseline;
meta_crossTime_decodingAcc_delay1;

meta_crossTime_decodingAcc_baseline_cell = cell(1,4);
meta_crossTime_decodingAcc_baseline_cell{1} = meta_crossTime_decodingAcc_baseline_length1;
meta_crossTime_decodingAcc_baseline_cell{2} = meta_crossTime_decodingAcc_baseline_length2;
meta_crossTime_decodingAcc_baseline_cell{3} = meta_crossTime_decodingAcc_baseline_length3;
meta_crossTime_decodingAcc_baseline_cell{4} = meta_crossTime_decodingAcc_baseline;

meta_crossTime_decodingAcc_delay1_cell = cell(1,4);
meta_crossTime_decodingAcc_delay1_cell{1} = meta_crossTime_decodingAcc_delay1_length1;
meta_crossTime_decodingAcc_delay1_cell{2} = meta_crossTime_decodingAcc_delay1_length2;
meta_crossTime_decodingAcc_delay1_cell{3} = meta_crossTime_decodingAcc_delay1_length3;
meta_crossTime_decodingAcc_delay1_cell{4} = meta_crossTime_decodingAcc_delay1;

for plot_lengthFlag=1:4
    if plot_lengthFlag == 1
        temp_trialBoolIndex_lengthx = trialBoolIndex_length1;
        x_raw = meta_trialLevel_crossTime_length1;
    elseif plot_lengthFlag == 2
        temp_trialBoolIndex_lengthx = trialBoolIndex_length2;
        x_raw = meta_trialLevel_crossTime_length2;
    elseif plot_lengthFlag == 3
        temp_trialBoolIndex_lengthx =  trialBoolIndex_length3;
        x_raw = meta_trialLevel_crossTime_length3;
    elseif plot_lengthFlag == 4
        temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
        %x_raw = meta_trialLevel_crossTime_length3;
        
        temp_meta_trialLevel_crossTime = nan(size(meta_trialLevel_crossTime_length1));
        
        for temptempi=1:size(temp_meta_trialLevel_crossTime,1)
            if trialBoolIndex_length1(temptempi) == true
                temp1 = meta_trialLevel_crossTime_length1(temptempi,:);
                if isnan(temp1(1)) == true
                    continue
                end
                temp_meta_trialLevel_crossTime(temptempi,:) = temp1;
            end
            if trialBoolIndex_length2(temptempi) == true
                temp1 = meta_trialLevel_crossTime_length2(temptempi,:);
                if isnan(temp1(1)) == true
                    continue
                end
                temp1 = temp1(end-length1_sample_interval(end)+1:end);%lengthx_sample_interval
                temp_meta_trialLevel_crossTime(temptempi,:) = temp1;
            end
            if trialBoolIndex_length3(temptempi) == true
                temp1 = meta_trialLevel_crossTime_length3(temptempi,:);
                if isnan(temp1(1)) == true
                    continue
                end
                temp1 = temp1(end-length1_sample_interval(end)+1:end);%lengthx_sample_interval
                temp_meta_trialLevel_crossTime(temptempi,:) = temp1;
            end
            
        end
        temp_meta_trialLevel_crossTime;
        
        x_raw = temp_meta_trialLevel_crossTime;
        
    end
    
    temptempBoolIndex = (~isnan(x_raw(:,1))) & temp_trialBoolIndex_lengthx;
    
    x = x_raw(temptempBoolIndex,:);
    
    y = y_raw(temptempBoolIndex);
    y = y==true;
    
    y_predict = x>=metaDecoderThreshold_delay1;
    
    temp_meta_crossTime_decodingAcc = nan(1,size(x,2));
    for tempj=1:size(x,2)
        temp_meta_crossTime_decodingAcc(tempj) = sum(y==y_predict(:,tempj))/length(y);
    end
    
    if plot_lengthFlag == 1
        meta_crossTime_decodingAcc_length1 = temp_meta_crossTime_decodingAcc;
    elseif plot_lengthFlag == 2
        meta_crossTime_decodingAcc_length2 = temp_meta_crossTime_decodingAcc;
    elseif plot_lengthFlag == 3
        meta_crossTime_decodingAcc_length3 = temp_meta_crossTime_decodingAcc;
    elseif plot_lengthFlag == 4
        meta_crossTime_decodingAcc_lastT = temp_meta_crossTime_decodingAcc;
    end
    
end
meta_crossTime_decodingAcc_length1;
meta_crossTime_decodingAcc_length2;
meta_crossTime_decodingAcc_length3;
meta_crossTime_decodingAcc_lastT;

meta_crossTime_decodingAcc_sample_cell = cell(1,4);
meta_crossTime_decodingAcc_sample_cell{1} = meta_crossTime_decodingAcc_length1;
meta_crossTime_decodingAcc_sample_cell{2} = meta_crossTime_decodingAcc_length2;
meta_crossTime_decodingAcc_sample_cell{3} = meta_crossTime_decodingAcc_length3;
meta_crossTime_decodingAcc_sample_cell{4} = meta_crossTime_decodingAcc_lastT;


%% Split meta_trialLevel_crossTime
% mismatch trials
trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh;
trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow;
trialBoolIndex_memoryPrecisionLowError_choiceMemory;
trialBoolIndex_memoryPrecisionHigh_choiceOffload;

a1 = sum(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh);
a2 = sum(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow);


% match trials
trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh;
trialBoolIndex_memoryPrecisionLow_choiceOffloadLow;

a3 = sum(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh);
a4 = sum(trialBoolIndex_memoryPrecisionLow_choiceOffloadLow);


%% Statistical quantification of temporal dynamics
% if if_compute == 1
if true
    if if_StatQuant == 1
        meta_trialLevel_resampleIter_crossTime_baseline;
        meta_trialLevel_resampleIter_crossTime_length1;
        meta_trialLevel_resampleIter_crossTime_length2;
        meta_trialLevel_resampleIter_crossTime_length3;
        meta_trialLevel_resampleIter_crossTime_delay1;
        
        a = 1;
        trialBoolIndex_length1;
        
        meta_trialLevel_resampleIter_crossTime_lastT = nan(size(meta_trialLevel_resampleIter_crossTime_length1));
        
        meta_trialLevel_resampleIter_crossTime_lastT(trialBoolIndex_length1,:,:) = ...
            meta_trialLevel_resampleIter_crossTime_length1(trialBoolIndex_length1,:,end-17:end);
        
        meta_trialLevel_resampleIter_crossTime_lastT(trialBoolIndex_length2,:,:) = ...
            meta_trialLevel_resampleIter_crossTime_length2(trialBoolIndex_length2,:,end-17:end);
        
        meta_trialLevel_resampleIter_crossTime_lastT(trialBoolIndex_length3,:,:) = ...
            meta_trialLevel_resampleIter_crossTime_length3(trialBoolIndex_length3,:,end-17:end);        
        
        
        temp_size_baseline = size(meta_trialLevel_resampleIter_crossTime_baseline,3);
        temp_size = size(meta_trialLevel_resampleIter_crossTime_lastT,3) + ...
            size(meta_trialLevel_resampleIter_crossTime_delay1(:,:,1:33),3);
        
        temp_r_baseline = nan(resampleIterCount_meta,temp_size_baseline);        
        temp_p = nan(resampleIterCount_meta,temp_size);
        temp_r = nan(resampleIterCount_meta,temp_size);        
        for tempi=1:resampleIterCount_meta
            
            temp1_baseline = squeeze(meta_trialLevel_resampleIter_crossTime_baseline(:,tempi,:));
            temp1_baseline_choiceMemory = temp1_baseline(choiceMemoryBoolIndex,:);
            
            temp1_baseline_choiceMemory_median = median(temp1_baseline_choiceMemory,'all','omitnan');
            if isnan(temp1_baseline_choiceMemory_median)
                continue
            end
                        
            temp1_lastT = squeeze(meta_trialLevel_resampleIter_crossTime_lastT(:,tempi,:));

            temp1_delay1 = squeeze(meta_trialLevel_resampleIter_crossTime_delay1(:,tempi,:));
            
            temp2 = [temp1_lastT,temp1_delay1(:,1:33)];

            %% Baseline
            temp_meta_seqLevel_choice_crossTime = zeros(sum(numSeq(1:3)),size(temp1_baseline,2));
            for temptempi=1:sum(numSeq(1:3))
                temptempBoolIndex = temp_seqIndex==temptempi;
                temptempBoolIndex_choice = temptempBoolIndex & choiceBoolIndex;
                temp_meta_seqLevel_choice_crossTime(temptempi,:) = mean(temp1_baseline(temptempBoolIndex_choice,:),'omitnan');
            end
            
            temp_x = temp_meta_seqLevel_choice_crossTime;
            temp_y = offloadingProb_inOne(1:sum(numSeq(1:3)))';
            
            temptemp_r = nan(1,size(temp_meta_seqLevel_choice_crossTime,2));
            temptemp_p = nan(1,size(temp_meta_seqLevel_choice_crossTime,2));
            for temptempi=1:size(temp_meta_seqLevel_choice_crossTime,2)
                temptemp_x = temp_x(:,temptempi);
                [temptemp_r(temptempi),temptemp_p(temptempi)] = corr(temptemp_x(~isnan(temptemp_x)),temp_y(~isnan(temptemp_x)));
            end
            
            temp_r_baseline(tempi,:) = temptemp_r;
            temptemp_r_baseline = temptemp_r;
            
            %% LastT ~ Delay1
            temp_meta_seqLevel_choice_crossTime = zeros(sum(numSeq(1:3)),size(temp2,2));
            for temptempi=1:sum(numSeq(1:3))
                temptempBoolIndex = temp_seqIndex==temptempi;
                temptempBoolIndex_choice = temptempBoolIndex & choiceBoolIndex;
                temp_meta_seqLevel_choice_crossTime(temptempi,:) = mean(temp2(temptempBoolIndex_choice,:),'omitnan');
            end
            
            temp_x = temp_meta_seqLevel_choice_crossTime;
            temp_y = offloadingProb_inOne(1:sum(numSeq(1:3)))';
            
            temptemp_r = nan(1,size(temp_meta_seqLevel_choice_crossTime,2));
            temptemp_p = nan(1,size(temp_meta_seqLevel_choice_crossTime,2));
            for temptempi=1:size(temp_meta_seqLevel_choice_crossTime,2)
                temptemp_x = temp_x(:,temptempi);
                %[temptemp_r(temptempi),temptemp_p(temptempi)] = corr(temptemp_x(~isnan(temptemp_x)),temp_y(~isnan(temptemp_x)));
                
                
                [temptemp_r(temptempi),~] = corr(temptemp_x(~isnan(temptemp_x)),temp_y(~isnan(temptemp_x)));
                [~,temptemp_p(temptempi)] = ttest(temptemp_r_baseline,temptemp_r(temptempi),'tail','right');

            end
            
            a = 1;
            
            
            temp_p(tempi,:) = temptemp_p;
            temp_r(tempi,:) = temptemp_r;
        end
        
        temp_p_valid = temp_p(~isnan(temp_p(:,1)),:);        
                
        temp_signifTimeStamp = temp_p_valid < p_threshold;
        
        temp_signifTimeStamp_meta = temp_signifTimeStamp;
        
        temp_r_dynamics_meta_baseline = temp_r_baseline;
        temp_r_dynamics_meta = temp_r;
        
        
        temp_signifTimeStamp_meta;
        
        temp_persistentTimeStamp = 5;
        
        temp_signifTimeStampCount_meta = nan(size(temp_signifTimeStamp_meta,1),1);
        for tempi=1:length(temp_signifTimeStampCount_meta)
            temp1 = temp_signifTimeStamp_meta(tempi,:);
            
            for tempj=1:(length(temp1)-temp_persistentTimeStamp+1)
               if sum(temp1(tempj:tempj+temp_persistentTimeStamp-1)) == temp_persistentTimeStamp
                   temp_signifTimeStampCount_meta(tempi) = tempj;
                   break
               end                
            end
            
        end
        temp_signifTimeStampCount_meta;
        
        temp_signifTimeStampCount_meta_mean = mean(temp_signifTimeStampCount_meta,'omitnan');
        temp_signifTimeStampCount_meta_sem = std(temp_signifTimeStampCount_meta,0,'omitnan')./sqrt(length(temp_signifTimeStampCount_meta));
        
        temp_r_dynamics_meta_baseline_mean = mean(temp_r_dynamics_meta_baseline,1,'omitnan');
        temp_r_dynamics_meta_mean = mean(temp_r_dynamics_meta,1,'omitnan');
        
        
        
        %% Plot
        temp_signifTimeStampCount_meta_mean;
        
        temp_r_dynamics_meta_whole = [temp_r_dynamics_meta_baseline,temp_r_dynamics_meta];
        
        temp_r_dynamics_meta_whole_mean = mean(temp_r_dynamics_meta_whole,1,'omitnan');
        temp_r_dynamics_meta_whole_sem = std(temp_r_dynamics_meta_whole,0,'omitnan')./sqrt(size(temp_r_dynamics_meta_whole,1));
        temp_r_dynamics_meta_whole_std = std(temp_r_dynamics_meta_whole,0,'omitnan');
        
        temp1 = temp_r_dynamics_meta_whole;
        temp1_mean = temp_r_dynamics_meta_whole_mean;
        temp1_sem = temp_r_dynamics_meta_whole_sem;
        temp1_std = temp_r_dynamics_meta_whole_std;
        
        temp2_timeStamp = temp_signifTimeStampCount_meta_mean + size(temp_r_dynamics_meta_baseline,2);        
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
%             patch([x(:); flipud(x(:))]', [y(:)+y_std(:); flipud(y(:)-y_std(:))]',[0.3 0.3 0.3],'FaceAlpha',0.05,'EdgeColor',[1 1 1]*0.725);            
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
            
            %y_sem = temp1_sem;
            %patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.3 0.3 0.3],'FaceAlpha',0.05,'EdgeColor',[1 1 1]*0.725);
            %y_std = temp1_std;
            %patch([x(:); flipud(x(:))]', [y(:)+y_std(:); flipud(y(:)-y_std(:))]',[0.3 0.3 0.3],'FaceAlpha',0.05,'EdgeColor',[1 1 1]*0.725);            
            %hold on
                
            [y_min,y_max] = bounds(y);
            
            if if_ttest_crossTime == 1
                scatter(x(temp2==1),y_min-(y_max-y_min)*0.24,8,[0 0 0],'*');
            end
            
            
            temp_interval_all = [baselinePeriod_interval(1),...
                (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)),...
                (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)) + temp_blankSize,...
                (baselinePeriod_interval(end) + temp_blankSize + length1_sample_interval(end) + decisionPeriodA_interval)];
                        
            text(mean(temp_interval_all(2:3)),y_max+(y_max-y_min)*0.17,...
                sprintf('%s','//'),'HorizontalAlignment','center','fontsize',8);
            
            xticks(temp_interval_all([1 2 3 4]));
            xticklabels({'Fixation','T1','LastT','Delay-on'});
            
            yticks([-1 0 0.5]);
            %if if_monkey_D0_Z1 == 1
                %yticks([-1 -0.5 0 0.5]);
            %end
            
            xlim([temp_interval_all(1) temp_interval_all(end-1)]);
            
            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 7.5)%12,8
            set(gca,'color',backgounrdColor);
            set(gca,'box','off');% 取消右、上边框
                        
            ylim([y_min-(y_max-y_min)*0.29 y_max+(y_max-y_min)*0.17]);
            if if_monkey_D0_Z1 == 1
                ylim([y_min-(y_max-y_min)*0.29 y_max+(y_max-y_min)*0.235]);
            end
            
            set(gca,'ydir','reverse');
            ylabel('Correlation', 'FontSize', 8, 'FontWeight', 'normal');%
            title('Meta-memory time course','FontSize',9);
            subtitle('Correlation between meta-memory and offloading rate','FontSize',7.5);
            
            xtickangle(0);
            
            
        end
        
    end
    
end


%% Plot
% p_threshold = 0.001;%0.05
if if_plot == 1
    %close all
    if if_plot_multiPanel == 1
        %backgounrdColor = [1 1 1]*0.915;%0.875
        backgounrdColor = [1 1 1];
        
        fig = figure('Name','asd','NumberTitle','off');
        %set(gcf,'Position',[10 50 1650 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 50 720 200]);
        t = tiledlayout(1,5,'TileSpacing','Compact','Padding','Compact');
        
        if if_trainPeriod_eachBin0_delay1Bin1 == 0
            t.Title.String = sprintf('Meta-memory dynamics, (trainEachBin), %s, p_threshold=%.3f',FOVName_currentFOV2,p_threshold);
        elseif if_trainPeriod_eachBin0_delay1Bin1 == 1
            t.Title.String = sprintf('Meta-memory dynamics, (trainDelay1Bin), %s, p_threshold=%.3f',FOVName_currentFOV2,p_threshold);
        end
        t.Title.FontSize = 10;%12
        t.Title.Interpreter = 'none';
        
        [y_min,y_max] = bounds([meta_trialLevel_crossTime_baseline(:);...
            meta_trialLevel_crossTime_length1(:);meta_trialLevel_crossTime_length2(:);meta_trialLevel_crossTime_length3(:);...
            meta_trialLevel_crossTime_delay1(:)]);
        
        for loopCount=0:4
            %for loopCount=1:4
            
            if loopCount==0
                temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_baseline;
                temp_str = 'Baseline';
                %nexttile(1);%4
                nexttile
                temp_interval = baselinePeriod_interval;
            elseif loopCount==1
                temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_length1;
                temp_str = 'Sample, length1';
                %nexttile(4);%2
                nexttile
                temp_interval = length1_sample_interval;
            elseif loopCount==2
                temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_length2;
                temp_str = 'Sample, length2';
                %nexttile(5);%5
                nexttile
                temp_interval = length2_sample_interval;
            elseif loopCount==3
                temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_length3;
                temp_str = 'Sample, length3';
                %nexttile(6);%8
                nexttile
                temp_interval = length3_sample_interval;
            elseif loopCount==4
                temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_delay1;
                temp_str = 'Delay1';
                %nexttile(3);%6
                nexttile
                temp_interval = decisionPeriodA_interval;
            end
            
            
            % mismatch trials
            if if_plot_fineTuningMismatch == 1
                temp_meta_trialLevel_crossTime_overMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh,:);
                temp_meta_trialLevel_crossTime_underMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow,:);
            elseif if_plot_fineTuningMismatch == 0
                %temp_meta_trialLevel_crossTime_overMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLowError_choiceMemory,:);
                %temp_meta_trialLevel_crossTime_underMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_choiceOffload,:);
                temp_meta_trialLevel_crossTime_overMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaHigh_choice,:);
                temp_meta_trialLevel_crossTime_underMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaLow_choice,:);
            end
            
            % match trials
            if if_plot_fineTuningMismatch == 1
                temp_meta_trialLevel_crossTime_highMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh,:);
                temp_meta_trialLevel_crossTime_lowMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_choiceOffloadLow,:);
            elseif if_plot_fineTuningMismatch == 0
                %temp_meta_trialLevel_crossTime_highMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory,:);
                %temp_meta_trialLevel_crossTime_lowMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_choiceOffload,:);
                temp_meta_trialLevel_crossTime_highMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaHigh_choice,:);
                temp_meta_trialLevel_crossTime_lowMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaLow_choice,:);
            end
            
            % all trials
            temp_meta_trialLevel_crossTime_choiceMemory = temp_meta_trialLevel_crossTime(choiceMemoryBoolIndex,:);
            temp_meta_trialLevel_crossTime_choiceOffload = temp_meta_trialLevel_crossTime(choiceOffloadBoolIndex,:);
            
            % all trials with error
            temp_meta_trialLevel_crossTime_choiceMemoryCorrect = temp_meta_trialLevel_crossTime(choiceMemoryCorrectBoolIndex,:);
            temp_meta_trialLevel_crossTime_choiceMemoryError = temp_meta_trialLevel_crossTime(choiceMemoryErrorBoolIndex,:);
            
            a1 = temp_meta_trialLevel_crossTime_highMatch;
            a2 = temp_meta_trialLevel_crossTime_lowMatch;
            a3 = temp_meta_trialLevel_crossTime_overMismatch;
            a4 = temp_meta_trialLevel_crossTime_underMismatch;
            a5 = temp_meta_trialLevel_crossTime_choiceMemory;
            a6 = temp_meta_trialLevel_crossTime_choiceOffload;
            a7 = temp_meta_trialLevel_crossTime_choiceMemoryCorrect;
            a8 = temp_meta_trialLevel_crossTime_choiceMemoryError;
            
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
            
            
            range_crossTime = 1:size(temp_meta_trialLevel_crossTime,2);
            x = range_crossTime;
            
            h_line = [];
            
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
                
                %scatter(x(p_a5_a6<0.001),y_max-(y_max-y_min)*0.15,8,[0 0 0],'*');
                scatter(x(p_a5_a6<p_threshold),y_max-(y_max-y_min)*0.15,8,[0 0 0],'*');
                
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
            %if loopCount == 1
            if false
                if if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 1 %#ok<*UNRCH>
                    %                 le = legend(h_line,sprintf('highMatch'),...
                    %                     sprintf('lowMatch'),...
                    %                     sprintf('overMismatch'),...
                    %                     sprintf('underMismatch'),...
                    %                     'Location','northeast','fontsize',10,'NumColumns',4);
                    le = legend(h_line,...
                        sprintf('Over-mismatch'),...
                        sprintf('High-match'),...
                        sprintf('Low-match'),...
                        sprintf('Under-mismatch'),...
                        'Location','northeast','fontsize',8,'NumColumns',4);
                    
                elseif if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 0
                    le = legend(h_line,sprintf('Choice-memory'),...
                        sprintf('Choice-offload'),...
                        'Location','northeast','fontsize',10,'NumColumns',4);
                end
                
                le.ItemTokenSize = ones(1,4)*14;%20
                le.Color = backgounrdColor;
            end
            
            xlim([range_crossTime(1) range_crossTime(end)]);
            %if if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 1
            %ylim([y_min-(y_max-y_min)*0.0 y_max+(y_max-y_min)*0]);%0.1
            %ylim([y_min+(y_max-y_min)*0.25 y_max-(y_max-y_min)*0.08]);
            %elseif if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 0
            %ylim([y_min+(y_max-y_min)*0.33 y_max-(y_max-y_min)*0.13]);
            %end
            ylim([y_min+(y_max-y_min)*0.25 y_max-(y_max-y_min)*0.08]);
            
            %ylim([0 1]);
            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 10)%12
            set(gca,'color',backgounrdColor);
            set(gca,'box','off');% 取消右、上边框
            %xlabel('Frame', 'FontSize', 14, 'FontWeight', 'bold');
            %if loopCount == 1
            if loopCount == 0
                %ylabel('Meta-memory', 'FontSize', 12, 'FontWeight', 'bold');
                ylabel('Meta-memory', 'FontSize', 10, 'FontWeight', 'normal');
            end
            xtickangle(0);
            
            %if if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 1
            %temp_title = title(sprintf('%s, [%d,%d,%d,%d]',...
            %    temp_str,temptempTrialNum_a1,temptempTrialNum_a2,temptempTrialNum_a3,temptempTrialNum_a4),...
            %    'FontSize',12,'FontWeight','normal');%'bold'
            %    temp_title = title(sprintf('%s',temp_str),'FontSize',10,'FontWeight','normal');%'bold'
            %elseif if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 0
            %temp_title = title(sprintf('%s, [%d,%d]',...
            %    temp_str,temptempTrialNum_a5,temptempTrialNum_a6),...
            %    'FontSize',12,'FontWeight','normal');%'bold
            %    temp_title = title(sprintf('%s',temp_str),'FontSize',10,'FontWeight','normal');%'bold
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
        %t = tiledlayout(1,2,'TileSpacing','loose','Padding','loose'); %#ok<*NASGU>
        
        %t.Title.String = sprintf('Population, n=%d',size(F_dff_decisionPeriodA,1));
        %t.Title.FontSize = 11;
        %t.Title.Interpreter = 'none';
        
        
        nexttile
        
        [y_min,y_max] = bounds([meta_trialLevel_crossTime_delay1(:)]);
        
        
        temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_delay1;
        temp_str = 'Delay1';
        temp_interval = decisionPeriodA_interval;
        
        
        % mismatch trials
        if if_plot_fineTuningMismatch == 1
            temp_meta_trialLevel_crossTime_overMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh,:);
            temp_meta_trialLevel_crossTime_underMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow,:);
        elseif if_plot_fineTuningMismatch == 0
            %temp_meta_trialLevel_crossTime_overMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLowError_choiceMemory,:);
            %temp_meta_trialLevel_crossTime_underMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_choiceOffload,:);
            temp_meta_trialLevel_crossTime_overMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaHigh_choice,:);
            temp_meta_trialLevel_crossTime_underMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaLow_choice,:);
        end
        
        % match trials
        if if_plot_fineTuningMismatch == 1
            temp_meta_trialLevel_crossTime_highMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh,:);
            temp_meta_trialLevel_crossTime_lowMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_choiceOffloadLow,:);
        elseif if_plot_fineTuningMismatch == 0
            %temp_meta_trialLevel_crossTime_highMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory,:);
            %temp_meta_trialLevel_crossTime_lowMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_choiceOffload,:);
            temp_meta_trialLevel_crossTime_highMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaHigh_choice,:);
            temp_meta_trialLevel_crossTime_lowMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaLow_choice,:);
        end
        
        % all trials
        temp_meta_trialLevel_crossTime_choiceMemory = temp_meta_trialLevel_crossTime(choiceMemoryBoolIndex,:);
        temp_meta_trialLevel_crossTime_choiceOffload = temp_meta_trialLevel_crossTime(choiceOffloadBoolIndex,:);
        
        % all trials with error
        temp_meta_trialLevel_crossTime_choiceMemoryCorrect = temp_meta_trialLevel_crossTime(choiceMemoryCorrectBoolIndex,:);
        temp_meta_trialLevel_crossTime_choiceMemoryError = temp_meta_trialLevel_crossTime(choiceMemoryErrorBoolIndex,:);
        
        a1 = temp_meta_trialLevel_crossTime_highMatch;
        a2 = temp_meta_trialLevel_crossTime_lowMatch;
        a3 = temp_meta_trialLevel_crossTime_overMismatch;
        a4 = temp_meta_trialLevel_crossTime_underMismatch;
        a5 = temp_meta_trialLevel_crossTime_choiceMemory;
        a6 = temp_meta_trialLevel_crossTime_choiceOffload;
        a7 = temp_meta_trialLevel_crossTime_choiceMemoryCorrect;
        a8 = temp_meta_trialLevel_crossTime_choiceMemoryError;
        
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
        
        
        range_crossTime = 1:size(temp_meta_trialLevel_crossTime,2);
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
            h_line = [h_line plot(x,y,'color',color_choiceOffloadLow,'linewidth',1)];
            hold on
            y_sem = a3_sem;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
            hold on
            
            
            y = a2_mean;
            h_line = [h_line plot(x,y,'color',color_choiceMemoryLow,'linewidth',1)];
            hold on
            y_sem = a2_sem;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
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
        
        
        %         le = legend(h_line,...
        %             sprintf('overMismatch'),...
        %             sprintf('highMatch'),...
        %             sprintf('lowMatch'),...
        %             sprintf('underMismatch'),...
        %             'Location','northeast','fontsize',10,'NumColumns',4);
        %
        %         le.ItemTokenSize = ones(1,4)*13;%20
        %         le.Color = backgounrdColor;
        
        %         le = legend(h_line,...
        %             sprintf('High-match'),...
        %             sprintf('Under-mismatch'),...
        %             sprintf('Low-match'),...
        %             sprintf('Over-mismatch'),...
        %             'Location','northeast','fontsize',7,'NumColumns',4);
        %         le.ItemTokenSize = ones(1,4)*13;%20
        %         le.Color = backgounrdColor;
        
        xticks(temp_interval(1:end-1));
        xticklabels({'Delay1','ChoiceCue'});
        xtickangle(0);
        
        yticks([0.2 0.8]);
        
        set(gca,'linewidth',1.5)
        xlim([range_crossTime(1) range_crossTime(end)]);
        %ylim([y_min-(y_max-y_min)*0.0 y_max+(y_max-y_min)*0]);%0.1
        ylim([y_min+(y_max-y_min)*0.15 y_max-(y_max-y_min)*0.05]);
        set(gca, 'FontSize', 10)
        set(gca,'color',backgounrdColor);
        set(gca,'box','off');% 取消右、上边框
        %xlabel('Frame', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('Meta-memory', 'FontSize', 10, 'FontWeight', 'normal');
        
        %temp_title = title(sprintf('Population'),'FontSize',12,'FontWeight','normal');%'bold'
        %temp_title.Interpreter = 'none';
        
        
        %% bar plot
        nexttile
        
        meta_trialLevel;
        %memoryMetaMismatch;
        trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh;
        trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh;
        trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow;
        trialBoolIndex_memoryPrecisionLow_choiceOffloadLow;
        
        temp_meta_highMatch = meta_trialLevel(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh);
        temp_meta_OverMismatch = meta_trialLevel(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh);
        temp_meta_underMismatch = meta_trialLevel(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow);
        temp_meta_lowMatch = meta_trialLevel(trialBoolIndex_memoryPrecisionLow_choiceOffloadLow);
        
        
        %         temp_1 = temp_meta_OverMismatch;
        %         temp_2 = temp_meta_highMatch;
        %         temp_3 = temp_meta_lowMatch;
        %         temp_4 = temp_meta_underMismatch;
        
        temp_1 = temp_meta_highMatch;
        temp_2 = temp_meta_underMismatch;
        temp_3 = temp_meta_lowMatch;
        temp_4 = temp_meta_OverMismatch;
        
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
        %ylabel('Meta-memory', 'FontSize', 10, 'FontWeight', 'normal');
        
        
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
        
        %for plot_lengthFlag=1:3
        for plot_lengthFlag=1:7
            
            
            if plot_lengthFlag == 1
                meta_trialLevel_crossTime_lengthx = meta_trialLevel_crossTime_length1;
                lengthx_sample_interval = length1_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length1;
            elseif plot_lengthFlag == 2
                meta_trialLevel_crossTime_lengthx = meta_trialLevel_crossTime_length2;
                lengthx_sample_interval = length2_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length2;
            elseif plot_lengthFlag == 3
                meta_trialLevel_crossTime_lengthx = meta_trialLevel_crossTime_length3;
                lengthx_sample_interval = length3_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length3;
            elseif plot_lengthFlag == 4
                lengthx_sample_interval = length3_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
                
                
                temp_meta_trialLevel_crossTime = nan(size(meta_trialLevel_crossTime_length3));
                
                for temptempi=1:size(temp_meta_trialLevel_crossTime,1)
                    if trialBoolIndex_length1(temptempi) == true
                        temp1 = meta_trialLevel_crossTime_length1(temptempi,:);
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
                        
                        temp_meta_trialLevel_crossTime(temptempi,:) = temp1_interp;
                    end
                    if trialBoolIndex_length2(temptempi) == true
                        temp1 = meta_trialLevel_crossTime_length2(temptempi,:);
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
                        
                        temp_meta_trialLevel_crossTime(temptempi,:) = temp1_interp;
                    end
                    if trialBoolIndex_length3(temptempi) == true
                        temp_meta_trialLevel_crossTime(temptempi,:) = meta_trialLevel_crossTime_length3(temptempi,:);
                    end
                    
                end
                temp_meta_trialLevel_crossTime;
                
                meta_trialLevel_crossTime_lengthx = temp_meta_trialLevel_crossTime;
                
            elseif plot_lengthFlag == 5
                lengthx_sample_interval = length3_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
                
                meta_trialLevel_crossTime_lengthx = meta_trialLevel_crossTime_lengthx;%#ok<*ASGSL>
                
            elseif plot_lengthFlag == 6
                lengthx_sample_interval = length1_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
                
                
                temp_meta_trialLevel_crossTime = nan(size(meta_trialLevel_crossTime_length1));
                
                for temptempi=1:size(temp_meta_trialLevel_crossTime,1)
                    if trialBoolIndex_length1(temptempi) == true
                        temp1 = meta_trialLevel_crossTime_length1(temptempi,:);
                        if isnan(temp1(1)) == true
                            continue
                        end
                        temp_meta_trialLevel_crossTime(temptempi,:) = temp1;
                    end
                    if trialBoolIndex_length2(temptempi) == true
                        temp1 = meta_trialLevel_crossTime_length2(temptempi,:);
                        if isnan(temp1(1)) == true
                            continue
                        end
                        temp1 = temp1(end-lengthx_sample_interval(end)+1:end);
                        temp_meta_trialLevel_crossTime(temptempi,:) = temp1;
                    end
                    if trialBoolIndex_length3(temptempi) == true
                        temp1 = meta_trialLevel_crossTime_length3(temptempi,:);
                        if isnan(temp1(1)) == true
                            continue
                        end
                        temp1 = temp1(end-lengthx_sample_interval(end)+1:end);
                        temp_meta_trialLevel_crossTime(temptempi,:) = temp1;
                    end
                    
                end
                temp_meta_trialLevel_crossTime;
                
                meta_trialLevel_crossTime_lengthx = temp_meta_trialLevel_crossTime;
                
            elseif plot_lengthFlag == 7
                lengthx_sample_interval = length1_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
                
                meta_trialLevel_crossTime_lengthx = meta_trialLevel_crossTime_lengthx;%#ok<*ASGSL>
                
            end
            
            %% Plot each length
            if true
                fig = figure('Name','asd','NumberTitle','off');
                %set(gcf,'Position',[10 50+200*(plot_lengthFlag-1) 360*0.975 200*0.8*0.8]);
                
                %temptemp1 = 360*0.975*0.74;
                %temptemp1 = 360*0.975*0.74*1.37*0.97*1.02;
                temptemp1 = 360*0.975*0.74*1.37*0.97*1.02*1.015;                
                %temptemp2 = 200*0.8*0.8;
                temptemp2 = 200*0.8*0.8*0.87;
                
                if plot_lengthFlag <=4
                    set(gcf,'Position',[10 50+200*(plot_lengthFlag-1) temptemp1 temptemp2]);
                elseif plot_lengthFlag >= 5
                    %set(gcf,'Position',[10+400 50+200*(plot_lengthFlag-2) temptemp1 temptemp2]);
                    set(gcf,'Position',[10+400*(plot_lengthFlag-4) 50+600 temptemp1 temptemp2]);
                end
                t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
                
                nexttile
                
                [y_min,y_max] = bounds([meta_trialLevel_crossTime_baseline(:);...
                    meta_trialLevel_crossTime_length1(:);meta_trialLevel_crossTime_length2(:);meta_trialLevel_crossTime_length3(:);...
                    meta_trialLevel_crossTime_delay1(:)]);
                
                %                 temp1 = meta_trialLevel_crossTime_delay1;
                %                 temp1 = temp1(:,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
                %                 [y_min,y_max] = bounds([meta_trialLevel_crossTime_baseline(:);...
                %                     meta_trialLevel_crossTime_lengthx(:);...
                %                     temp1(:)]);
                
                for loopCount=0:2
                    
                    if loopCount==0
                        temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_baseline;
                        temp_interval = baselinePeriod_interval;
                    elseif loopCount==1
                        temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_lengthx;
                        temp_interval = baselinePeriod_interval(end) + lengthx_sample_interval;
                    elseif loopCount==2
                        temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_delay1;
                        temp_interval = baselinePeriod_interval(end) + lengthx_sample_interval(end) + decisionPeriodA_interval;
                    end
                    
                    
                    % mismatch trials
                    if if_plot_fineTuningMismatch == 1
                        temp_meta_trialLevel_crossTime_overMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh & temp_trialBoolIndex_lengthx,:);
                        temp_meta_trialLevel_crossTime_underMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow & temp_trialBoolIndex_lengthx,:);
                    elseif if_plot_fineTuningMismatch == 0
                        temp_meta_trialLevel_crossTime_overMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaHigh_choice & temp_trialBoolIndex_lengthx,:);
                        temp_meta_trialLevel_crossTime_underMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaLow_choice & temp_trialBoolIndex_lengthx,:);
                    end
                    
                    % match trials
                    if if_plot_fineTuningMismatch == 1
                        temp_meta_trialLevel_crossTime_highMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh & temp_trialBoolIndex_lengthx,:);
                        temp_meta_trialLevel_crossTime_lowMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_choiceOffloadLow & temp_trialBoolIndex_lengthx,:);
                    elseif if_plot_fineTuningMismatch == 0
                        temp_meta_trialLevel_crossTime_highMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaHigh_choice & temp_trialBoolIndex_lengthx,:);
                        temp_meta_trialLevel_crossTime_lowMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaLow_choice & temp_trialBoolIndex_lengthx,:);
                    end
                    
                    % all trials
                    temp_meta_trialLevel_crossTime_choiceMemory = temp_meta_trialLevel_crossTime(choiceMemoryBoolIndex' & temp_trialBoolIndex_lengthx,:);
                    temp_meta_trialLevel_crossTime_choiceOffload = temp_meta_trialLevel_crossTime(choiceOffloadBoolIndex' & temp_trialBoolIndex_lengthx,:);
                    
                    % all trials with error
                    temp_meta_trialLevel_crossTime_choiceMemoryCorrect = temp_meta_trialLevel_crossTime(choiceMemoryCorrectBoolIndex' & temp_trialBoolIndex_lengthx,:);
                    temp_meta_trialLevel_crossTime_choiceMemoryError = temp_meta_trialLevel_crossTime(choiceMemoryErrorBoolIndex' & temp_trialBoolIndex_lengthx,:);
                    
                    a1 = temp_meta_trialLevel_crossTime_highMatch;
                    a2 = temp_meta_trialLevel_crossTime_lowMatch;
                    a3 = temp_meta_trialLevel_crossTime_overMismatch;
                    a4 = temp_meta_trialLevel_crossTime_underMismatch;
                    a5 = temp_meta_trialLevel_crossTime_choiceMemory;
                    a6 = temp_meta_trialLevel_crossTime_choiceOffload;
                    a7 = temp_meta_trialLevel_crossTime_choiceMemoryCorrect;
                    a8 = temp_meta_trialLevel_crossTime_choiceMemoryError;
                    
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
                    
                    if if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 3
                        if loopCount == 0
                            diff_baseline_56 = a5_mean-a6_mean;
                            diff_baseline_median_56 = median(diff_baseline_56);
                        else
                            [~,p_a5_a6] = ttest2(a5,a6+diff_baseline_median_56);
                        end
                    end
                    
                    if loopCount == 0
                        diff_baseline = a7_mean-a6_mean;
                        diff_baseline_median = median(diff_baseline);
                    else
                        [~,p_a7_a6] = ttest2(a7,a6+diff_baseline_median);
                    end
                    
                    temp_meta_seqLevel_choice_crossTime = zeros(sum(numSeq(1:3)),size(temp_meta_trialLevel_crossTime,2));
                    for temptempi=1:sum(numSeq(1:3))
                        temptempBoolIndex = temp_seqIndex==temptempi;
                        temptempBoolIndex_choice = temptempBoolIndex & choiceBoolIndex;
                        temp_meta_seqLevel_choice_crossTime(temptempi,:) = mean(temp_meta_trialLevel_crossTime(temptempBoolIndex_choice,:),'omitnan');
                    end
                    
                    temp_x = temp_meta_seqLevel_choice_crossTime;
                    temp_y = offloadingProb_inOne(1:sum(numSeq(1:3)))';
                    
                    temp_r = nan(1,size(temp_meta_seqLevel_choice_crossTime,2));
                    temp_p = nan(1,size(temp_meta_seqLevel_choice_crossTime,2));
                    for temptempi=1:size(temp_meta_seqLevel_choice_crossTime,2)
                        temptemp_x = temp_x(:,temptempi);
                        [temp_r(temptempi),temp_p(temptempi)] = corr(temptemp_x(~isnan(temptemp_x)),temp_y(~isnan(temptemp_x)));
                    end
                    
                    a9 = temp_r;
                    
                    
                    
                    if loopCount == 0
                        a9_baseline = a9;                  
                        a9_baseline_median = median(a9);
                    else
                        p_a9 = nan(size(a9));
                        for temptempi=1:length(a9)
                            [~,p_a9(temptempi)] = ttest(a9_baseline,a9(temptempi),'tail','right');
                            
                            %[~,p_a9(temptempi)] = ttest(a9_baseline,a9(temptempi));                            
                        end
                    end
                    
                    p_a9 = temp_p;
                    
                    
                    range_crossTime = (1:size(temp_meta_trialLevel_crossTime,2))+temp_interval(1)-1;
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
                            
                            scatter(x(p_a5_a6<p_threshold),y_max-(y_max-y_min)*0.15,8,[0 0 0],'*');
                            
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
                            
                            

                            
                            if loopCount >= 1
                                if if_ttest_crossTime == 1
                                    scatter(x(p_a7_a6<p_threshold),y_max-(y_max-y_min)*0.15,8,[0 0 0],'*');
                                end
                            end
                            
                            
                        elseif if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 3
                            %y = a7_mean-a6_mean;
                            y = a5_mean-a6_mean;
                            h_line = [h_line plot(x,y,'color',[1 1 1]*0.3,'linewidth',1)]; %#ok<*AGROW>
                            hold on
                            
                            y_min = -0.25;
                            y_max = 0.4;
                            
                            if loopCount >= 1
                                if if_ttest_crossTime == 1
                                    %scatter(x(p_a7_a6<p_threshold),y_max-(y_max-y_min)*0.15,8,[0 0 0],'*');
                                    scatter(x(p_a5_a6<p_threshold),y_max-(y_max-y_min)*0.15,8,[0 0 0],'*');
                                end
                            end
                            
                        end
                        
                    elseif plot_lengthFlag == 5 || plot_lengthFlag == 7
                        y = a9;
                        h_line = [h_line plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1)]; %#ok<*AGROW>
                        hold on
                        
                        y_min = -1.1;%-1
                        y_max = 0.5;%0.75
                        
                        %if if_monkey_D0_Z1 == 0
                        %    if currentSessionIndex_AB == 8
                        %        y_min = -1;
                        %        y_max = 0.5;
                        %    end
                        %end
                            
                        %if loopCount == 2
                        %    plot([1 x(end)],[1 1]*a9_baseline_median,':','color',[0.3 0.3 0.3],'linewidth',1);
                        %    hold on
                        %end
                        
                        if loopCount >= 1
                            if if_ttest_crossTime == 1
                                %scatter(x(p_a9<p_threshold),y_max-(y_max-y_min)*0.15,8,[0 0 0],'*');
                                %scatter(x(p_a9<p_threshold),0.5,8,[0 0 0],'*');
                                scatter(x(p_a9<p_threshold),-1.05,8,[0 0 0],'*');
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
                
                if if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 2
                    if plot_lengthFlag == 4
                        %le = legend(h_line,'Choice-memory correct','Choice-offload','Choice-memory error',...
                        %    'Location','northwest','fontsize',6.5);%10
                        le = legend(h_line,'Memory-correct','Offload','Memory-error',...
                           'Location','northwest','fontsize',6.5);%10                        
                        le.ItemTokenSize = ones(1,3)*10;
                        le.Color = backgounrdColor;
                        legend('boxoff');
                    end
                end
                
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
                    %yticks([-1 0]);
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
                
                if plot_lengthFlag <= 4 || plot_lengthFlag == 6
                    ylim([y_min+(y_max-y_min)*0.3 y_max-(y_max-y_min)*0.08]);
                    ylabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');%
                    
                    %if if_allTrial0_matchMismatch1_allTrialWithError2_correct3 == 3
                    %    ylim([-0.05 0.3]);
                    %end
                    if plot_lengthFlag == 6
                        title('Meta-memory time course','FontSize',9);
                        subtitle('Choice-memory VS. Choice-offload','FontSize',7.5);
                    end
                    
                elseif plot_lengthFlag == 5 || plot_lengthFlag == 7
                    ylim([y_min y_max]);
                    
                    set(gca,'ydir','reverse');
                    ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');%
                    %title('Meta-memory VS. offloading rate','FontSize',9);%10
                    title('Meta-memory time course','FontSize',9);
                    subtitle('Correlation between meta-memory and offloading rate','FontSize',7.5);
                end
                
                xtickangle(0);
                
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
                meta_trialLevel_crossTime_lengthx = meta_trialLevel_crossTime_length1;
                lengthx_sample_interval = length1_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length1;
            elseif plot_lengthFlag == 2
                meta_trialLevel_crossTime_lengthx = meta_trialLevel_crossTime_length2;
                lengthx_sample_interval = length2_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length2;
            elseif plot_lengthFlag == 3
                meta_trialLevel_crossTime_lengthx = meta_trialLevel_crossTime_length3;
                lengthx_sample_interval = length3_sample_interval;
                temp_trialBoolIndex_lengthx = trialBoolIndex_length3;
            end
            
            %temp_trialBoolIndex_lengthx = temp_trialBoolIndex_lengthx & tempBoolIndex_4typesFilter;
            temp_trialBoolIndex_lengthx = temp_trialBoolIndex_lengthx;
            
            
            fig = figure('Name','asd','NumberTitle','off');
            
            temptemp1 = 360*0.975*0.74*1.37*0.97*1.02;
            temptemp2 = 200*0.8*0.8;
            
            set(gcf,'Position',[410 50+200*(plot_lengthFlag-1) temptemp1 temptemp2]);
            
            t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
            
            nexttile
            
            y_min = -1;
            y_max = 1;
            
            for loopCount=0:2
                
                if loopCount==0
                    temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_baseline;
                    temp_interval = baselinePeriod_interval;
                elseif loopCount==1
                    temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_lengthx;
                    temp_interval = baselinePeriod_interval(end) + lengthx_sample_interval;
                elseif loopCount==2
                    temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_delay1;
                    temp_interval = baselinePeriod_interval(end) + lengthx_sample_interval(end) + decisionPeriodA_interval;
                end
                
                
                temp_meta_seqLevel_choice_crossTime = zeros(sum(numSeq(1:3)),size(temp_meta_trialLevel_crossTime,2));
                for temptempi=1:sum(numSeq(1:3))
                    temptempBoolIndex = temp_seqIndex==temptempi;
                    temptempBoolIndex_choice = temptempBoolIndex & choiceBoolIndex;
                    temp_meta_seqLevel_choice_crossTime(temptempi,:) = mean(temp_meta_trialLevel_crossTime(temptempBoolIndex_choice,:),'omitnan');
                end
                
                temp_x = temp_meta_seqLevel_choice_crossTime;
                temp_y = offloadingProb_inOne(1:sum(numSeq(1:3)))';
                
                temp_r = nan(1,size(temp_meta_seqLevel_choice_crossTime,2));
                temp_p = nan(1,size(temp_meta_seqLevel_choice_crossTime,2));
                for temptempi=1:size(temp_meta_seqLevel_choice_crossTime,2)
                    temptemp_x = temp_x(:,temptempi);
                    [temp_r(temptempi),temp_p(temptempi)] = corr(temptemp_x(~isnan(temptemp_x)),temp_y(~isnan(temptemp_x)));
                end
                
                p_a9 = temp_p;
                
                a9 = temp_r;
                if loopCount == 0
                    a9_baseline = a9;
                    a9_baseline_median = median(a9);
                end
                
                %                 p_a9 = nan(size(a9));
                %                 if loopCount >= 1
                %                     for temptempi=1:length(p_a9)
                %                         %[~,p_a9(temptempi)] = ttest(a9_baseline,a9(temptempi));
                %                         [~,p_a9(temptempi)] = ttest(a9_baseline,a9(temptempi),'tail','right');
                %                     end
                %                 end
                
                
                
                range_crossTime = (1:size(temp_meta_trialLevel_crossTime,2))+temp_interval(1)-1;
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
                    scatter(x(p_a9<p_threshold),y_max-(y_max-y_min)*0.15,8,[0 0 0],'*');
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
            
            ylim([-1 1]);
            set(gca,'ydir','reverse');
            
            ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
            title('Meta-memory VS. offloading rate','FontSize',9);%10
            
            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 8)%12
            set(gca,'color',backgounrdColor);
            set(gca,'box','off');% 取消右、上边框
            xtickangle(0);
            
            
        end
        
    end
    
end


%% Plot multiPanelD
if if_plot == 1
    if if_plot_multiPanelD == 1
        backgounrdColor = [1 1 1];
        
        y_min = 0.2;
        y_max = 0.9;
        
        for plot_lengthFlag=1:4
            if plot_lengthFlag == 1
                lengthx_sample_interval = length1_sample_interval;
            elseif plot_lengthFlag == 2
                lengthx_sample_interval = length2_sample_interval;
            elseif plot_lengthFlag == 3
                lengthx_sample_interval = length3_sample_interval;
            elseif plot_lengthFlag == 4
                lengthx_sample_interval = length1_sample_interval;
            end
            
            %% Plot each length
            if true
                fig = figure('Name','asd','NumberTitle','off');
                
                temptemp1 = 360*0.975*0.74*1.37*0.97*1.02;
                temptemp2 = 200*0.8*0.8;
                
                set(gcf,'Position',[10 50+200*(plot_lengthFlag-1) temptemp1 temptemp2]);
                t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
                
                nexttile
                
                
                for loopCount=0:2
                    
                    if loopCount==0
                        temp_meta_trialLevel_crossTime = meta_crossTime_decodingAcc_baseline_cell{plot_lengthFlag};
                        temp_interval = baselinePeriod_interval;
                    elseif loopCount==1
                        temp_meta_trialLevel_crossTime = meta_crossTime_decodingAcc_sample_cell{plot_lengthFlag};
                        temp_interval = baselinePeriod_interval(end) + lengthx_sample_interval;
                    elseif loopCount==2
                        temp_meta_trialLevel_crossTime = meta_crossTime_decodingAcc_delay1_cell{plot_lengthFlag};
                        temp_interval = baselinePeriod_interval(end) + lengthx_sample_interval(end) + decisionPeriodA_interval;
                    end
                    
                    
                    a1 = temp_meta_trialLevel_crossTime;
                    
                    
                    if loopCount == 0
                        diff_baseline = a1;
                        diff_baseline_median = median(diff_baseline);
                        diff_baseline_prctile = prctile(diff_baseline,99.9);
                    else
                        p_a1 = nan(1,length(a1));
                        for temptempi=1:length(a1)
                            [~,p_a1(temptempi)] = ttest(diff_baseline,a1(temptempi));
                        end
                    end
                    
                    
                    range_crossTime = (1:size(temp_meta_trialLevel_crossTime,2))+temp_interval(1)-1;
                    x = range_crossTime;
                    
                    y = a1;
                    h_line = [h_line plot(x,y,'color',[1 1 1]*0.3,'linewidth',1)]; %#ok<*AGROW>
                    hold on
                    
                    
                    if loopCount >= 1
                        if if_ttest_crossTime == 1
                            %scatter(x(p_a1<p_threshold),y_max-(y_max-y_min)*0.15,8,[0 0 0],'*');
                            %scatter(x(p_a1<0.001),y_max-(y_max-y_min)*0.15,8,[0 0 0],'*');
                            
                            scatter(x(a1>diff_baseline_prctile),y_max-(y_max-y_min)*0.15,8,[0 0 0],'*');                            
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
                elseif plot_lengthFlag == 4
                    xticks(temp_interval_all([1 2 3]));
                    xticklabels({'Fixation','Last T','Delay1'});
                end
                
                xlim([temp_interval_all(1) temp_interval_all(end-1)]);
                
                set(gca,'linewidth',1.5)
                set(gca, 'FontSize', 8)%12
                set(gca,'color',backgounrdColor);
                set(gca,'box','off');% 取消右、上边框
                
                ylim([y_min+(y_max-y_min)*0.3 y_max-(y_max-y_min)*0.08]);
                ylabel(sprintf('Meta-memory\nDecoding accuracy'), 'FontSize', 9, 'FontWeight', 'normal');%
                
                xtickangle(0);
                
            end
            
        end
        
    end
    
end



%% New



%% End