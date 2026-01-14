function output = fun_wrapped_locSVM_v4(options)

if_normal0_resample1_leaveOneSeqOut2 = options.if_normal0_resample1_leaveOneSeqOut2;
seqCount = options.seqCount;
if_resample_seqBased0_locBased1 = options.if_resample_seqBased0_locBased1;
resampleTrialCount = options.resampleTrialCount;
numSeq = options.numSeq;
targetPATH = options.targetPATH;

numFrames = options.numFrames;
seqIndex_valid = options.seqIndex_valid;
boolIndex_location_seq_T = options.boolIndex_location_seq_T;
resampleIterCount = options.resampleIterCount;
temp_F_dff_decisionBin = options.temp_F_dff_decisionBin;

gAcc_target_collapsed_inOne = options.gAcc_target_collapsed_inOne;
boolIndex_location_seq = options.boolIndex_location_seq;
KFold_num = options.KFold_num;
if_locDecoder_meanPosterior0_meanAcc1 = options.if_locDecoder_meanPosterior0_meanAcc1;
t_decoder = options.t_decoder;

if_decodingAcc_threshold0_sort1 = options.if_decodingAcc_threshold0_sort1;
boolIndex_location_trial = options.boolIndex_location_trial;


if if_normal0_resample1_leaveOneSeqOut2 == 0 || if_normal0_resample1_leaveOneSeqOut2 == 1
    %% Get resampleTrialCount_seq
    if if_normal0_resample1_leaveOneSeqOut2 == 0
        resampleTrialCount_seq = seqCount;
    elseif if_normal0_resample1_leaveOneSeqOut2 == 1
        
        if if_resample_seqBased0_locBased1 == 0
            resampleTrialCount_seq = resampleTrialCount*ones(1,sum(numSeq(1:3)));
            for target_length=1:3
                temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
                temp_seqCount = seqCount(temp_range);
                temp_resampleTrialCount = resampleTrialCount;
                
                %if target_length == 1
                %    if min(temp_seqCount) > temp_resampleTrialCount
                %        temp_resampleTrialCount = min(temp_seqCount);
                %    end
                %end
                if target_length > 1
                    a = 1; %#ok<*NASGU>
                    temptempBoolIndex = temp_seqCount >= temp_resampleTrialCount;
                    while sum(temptempBoolIndex) < numSeq(target_length)*0.4
                        temp_resampleTrialCount = temp_resampleTrialCount - 1;
                        temptempBoolIndex = temp_seqCount >= temp_resampleTrialCount;
                        if temp_resampleTrialCount == 4
                            break
                        end
                    end
                end
                fprintf('temp_resampleTrialCount=%d, for length %d.\n',temp_resampleTrialCount,target_length);
                
                resampleTrialCount_seq(temp_range) = temp_resampleTrialCount;
            end
            resampleTrialCount_seq(seqCount<resampleTrialCount) = 0;
            
        elseif if_resample_seqBased0_locBased1 == 1
            threshold_locRatio = 0.055;%0.06-->0.05-->0.045
            threshold_locStd = 0.026;%0.0278-->0.025-->0.026
            
            %threshold_locRatio = 0.1;%0.055
            %threshold_locStd = 0.05;%0.026
            
            %threshold_locRatio = 1;%0.055
            %threshold_locStd = 1;%0.026

            
            fun_locBasedResample_seqCount_Name_v = autoGetFunName_myScripts('fun_locBasedResample_seqCount', [targetPATH '\functions']);
            fun_locBasedResample_seqCount = str2func(fun_locBasedResample_seqCount_Name_v);
            
            options_locBased.threshold_locRatio = threshold_locRatio;
            options_locBased.threshold_locStd = threshold_locStd;
            options_locBased.numSeq = numSeq;
            options_locBased.numFrames = numFrames;
            options_locBased.seqIndex_valid = seqIndex_valid;
            options_locBased.boolIndex_location_seq_T = boolIndex_location_seq_T;
            
            resampleTrialCount_seq = fun_locBasedResample_seqCount(options_locBased);
        end
    end
    
    
    %% Get resampled XY
    fun_resampleXY_seqCountBased_Name_v = autoGetFunName_myScripts('fun_resampleXY_seqCountBased', [targetPATH '\functions']);
    fun_resampleXY_seqCountBased = str2func(fun_resampleXY_seqCountBased_Name_v);
    
    options_resampleXY.numSeq = numSeq;
    options_resampleXY.resampleIterCount = resampleIterCount;
    options_resampleXY.resampleTrialCount_seq = resampleTrialCount_seq;
    options_resampleXY.seqIndex_valid = seqIndex_valid;
    options_resampleXY.temp_F_dff_decisionBin = temp_F_dff_decisionBin;
    
    output_resampleXY = fun_resampleXY_seqCountBased(options_resampleXY);
    
    temp_F_dff_decisionBin_resampled = output_resampleXY.temp_F_dff_decisionBin_resample;
    seqIndex_valid_resampled = output_resampleXY.seqIndex_valid_resample;
    % temp_trialIndex_valid_resample here is for single sessions or two sessions, to trace trial-level memory precision
    temp_trialIndex_valid_resample = output_resampleXY.temp_trialIndex_valid_resample;
    
    %% Get resampled SVM
    fun_resampled_svm_train_lengthx_outputs_Name_v = autoGetFunName_myScripts('fun_resampled_svm_train_lengthx_outputs', [targetPATH '\functions']);
    fun_resampled_svm_train_lengthx_outputs = str2func(fun_resampled_svm_train_lengthx_outputs_Name_v);
    
    options_resampledSVM = struct;
    options_resampledSVM.numSeq = numSeq;
    options_resampledSVM.gAcc_target_collapsed_inOne = gAcc_target_collapsed_inOne;
    options_resampledSVM.seqIndex_valid_resampled = seqIndex_valid_resampled;
    options_resampledSVM.temp_F_dff_decisionBin_resampled = temp_F_dff_decisionBin_resampled;
    options_resampledSVM.boolIndex_location_seq = boolIndex_location_seq;
    options_resampledSVM.KFold_num = KFold_num;
    options_resampledSVM.resampleIterCount = resampleIterCount;
    options_resampledSVM.numFrames = numFrames;
    options_resampledSVM.if_locDecoder_meanPosterior0_meanAcc1 = if_locDecoder_meanPosterior0_meanAcc1;
    options_resampledSVM.targetPATH = targetPATH;
    options_resampledSVM.t_decoder = t_decoder;
    options_resampledSVM.if_normal0_resample1_leaveOneSeqOut2 = if_normal0_resample1_leaveOneSeqOut2;
    options_resampledSVM.if_decodingAcc_threshold0_sort1 = if_decodingAcc_threshold0_sort1;
    
    output_resampledSVM = fun_resampled_svm_train_lengthx_outputs(options_resampledSVM);
    
    svm_train_length1_outputs = output_resampledSVM.svm_train_length1_outputs;
    svm_train_length2_outputs = output_resampledSVM.svm_train_length2_outputs;
    svm_train_length3_outputs = output_resampledSVM.svm_train_length3_outputs;
    
    a = 1;
    
elseif if_normal0_resample1_leaveOneSeqOut2 == 2
    seqProbSVM_6location_trialLevel_train_Name_v = autoGetFunName_myScripts('seqProbSVM_6location_trialLevel_train', [targetPATH '\functions']);
    fun_seqProbSVM_6location_trialLevel_train = str2func(seqProbSVM_6location_trialLevel_train_Name_v);
    
    svm_options = struct;
    svm_options.boolIndex_location_seq = boolIndex_location_seq;
    svm_options.t_decoder = t_decoder;
    svm_options.numFrames = numFrames;
    svm_options.targetPATH = targetPATH;
    svm_options.if_normal0_resample1_leaveOneSeqOut2 = if_normal0_resample1_leaveOneSeqOut2;
    svm_options.if_decodingAcc_threshold0_sort1 = if_decodingAcc_threshold0_sort1;
    
    
    for target_length=1:3
        %for target_length=3:3
        temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
        gAcc = gAcc_target_collapsed_inOne(temp_range);
        
        svm_options.gAcc = gAcc;
        svm_options.seq_range = temp_range;
        svm_options.target_length = target_length;
        
        tempBoolIndex = ismember(seqIndex_valid,temp_range);
        x = temp_F_dff_decisionBin(:,tempBoolIndex);
        y = boolIndex_location_trial(:,tempBoolIndex);
        
        %% leave-one-seq-out mode
        a = 1;
        
        seqProbSVM_leaveOneSeqOut_Name_v = autoGetFunName_myScripts('seqProbSVM_6location_trialLevel_train_leaveOneSeqOut', [targetPATH '\functions']);
        fun_seqProbSVM_leaveOneSeqOut = str2func(seqProbSVM_leaveOneSeqOut_Name_v);
        
        svm_train_outputs = fun_seqProbSVM_leaveOneSeqOut(x,y,svm_options);
        
        
        temp_svm = svm_train_outputs;
        
        %% Compute seq-level mean location distribution
        temp_trialNum = size(temp_svm.temp_svm_Y_valid_T,1);
        temp_seqIndex_trial = zeros(temp_trialNum,1);
        for tempi=1:temp_trialNum
            temp_seqBoolIndex_trial = temp_svm.temp_svm_Y_valid_T(tempi,:);
            for tempj=1:numSeq(target_length)
                temp_seqIndex = temp_svm.seq_range(tempj);
                temp_seqBoolIndex = temp_svm.boolIndex_location_seq(:,temp_seqIndex)';
                
                if sum(temp_seqBoolIndex==temp_seqBoolIndex_trial) == numFrames
                    temp_seqIndex_trial(tempi) = temp_svm.seq_range(tempj);
                    break
                end
            end
        end
        temp_svm.seqIndex_trial = temp_seqIndex_trial;
        
        
        Posterior_2d_n11n_lengthx_mean = zeros(numSeq(target_length),numFrames);
        Posterior_2d_w_lengthx_mean = zeros(numSeq(target_length),numFrames);
        for tempi=1:numSeq(target_length)
            temp_seqIndex = temp_svm.seq_range(tempi);
            tempSeqBoolIndex = temp_svm.boolIndex_location_seq(:,temp_seqIndex)';
            
            tempBoolIndex = temp_seqIndex_trial == temp_seqIndex;
            temp_Posterior_2d_n11n = temp_svm.Posterior_2d_n11n(tempBoolIndex,:);
            temp_Posterior_2d_w = temp_svm.Posterior_2d_w(tempBoolIndex,:);
            
            if if_locDecoder_meanPosterior0_meanAcc1 == 0
                Posterior_2d_n11n_lengthx_mean(tempi,:) = mean(temp_Posterior_2d_n11n,1);
                Posterior_2d_w_lengthx_mean(tempi,:) = mean(temp_Posterior_2d_w,1);
            elseif if_locDecoder_meanPosterior0_meanAcc1 == 1
                tempBoolIndex2 = temp_Posterior_2d_n11n>0.5;
                for tempj=1:numFrames
                    tempBoolIndex3 = tempBoolIndex2(:,tempj) == tempSeqBoolIndex(tempj);
                    Posterior_2d_n11n_lengthx_mean(tempi,tempj) = sum(tempBoolIndex3)/length(tempBoolIndex3);
                end
                tempBoolIndex2 = temp_Posterior_2d_w>0.5;
                for tempj=1:numFrames
                    tempBoolIndex3 = tempBoolIndex2(:,tempj) == tempSeqBoolIndex(tempj);
                    Posterior_2d_w_lengthx_mean(tempi,tempj) = sum(tempBoolIndex3)/length(tempBoolIndex3);
                end
            end
        end
        temp_svm.Posterior_2d_n11n_lengthx_mean = Posterior_2d_n11n_lengthx_mean;
        temp_svm.Posterior_2d_w_lengthx_mean = Posterior_2d_w_lengthx_mean;
        
        svm_train_outputs = temp_svm;
        
        if target_length == 1
            svm_train_length1_outputs = svm_train_outputs;
        elseif target_length == 2
            svm_train_length2_outputs = svm_train_outputs;
        elseif target_length == 3
            svm_train_length3_outputs = svm_train_outputs;
        end
        
    end
    temp_trialIndex_valid_resample = [];
end

output = struct;
output.svm_train_length1_outputs = svm_train_length1_outputs;
output.svm_train_length2_outputs = svm_train_length2_outputs;
output.svm_train_length3_outputs = svm_train_length3_outputs;
output.temp_trialIndex_valid_resample = temp_trialIndex_valid_resample;

