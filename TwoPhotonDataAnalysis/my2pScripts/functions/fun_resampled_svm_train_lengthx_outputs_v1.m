function output = fun_resampled_svm_train_lengthx_outputs_v1(options)

numSeq = options.numSeq;

gAcc_target_collapsed_inOne = options.gAcc_target_collapsed_inOne;
seqIndex_valid_resampled = options.seqIndex_valid_resampled;
temp_F_dff_decisionBin_resampled = options.temp_F_dff_decisionBin_resampled;
boolIndex_location_seq = options.boolIndex_location_seq;
KFold_num = options.KFold_num;

resampleIterCount = options.resampleIterCount;
numFrames = options.numFrames;
if_locDecoder_meanPosterior0_meanAcc1 = options.if_locDecoder_meanPosterior0_meanAcc1;

targetPATH = options.targetPATH;
t_decoder = options.t_decoder;
if_normal0_resample1_leaveOneSeqOut2 = options.if_normal0_resample1_leaveOneSeqOut2;
if_decodingAcc_threshold0_sort1 = options.if_decodingAcc_threshold0_sort1;

% KFold_num = options.KFold_num;
% KFold_num = options.KFold_num;





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
    svm_options.if_weighted = 0;%1
    
    tempBoolIndex = ismember(seqIndex_valid_resampled,temp_range);
    x = temp_F_dff_decisionBin_resampled(:,:,tempBoolIndex);
    y = boolIndex_location_seq(:,seqIndex_valid_resampled(tempBoolIndex));
    
    %temp_resampleTrialCount = resampleTrialCount;
    
    c = cvpartition(size(y,2),'KFold',KFold_num);
    
    
    temp_svm_resample = cell(1,resampleIterCount);
    % It seems that parfor would inevitably cause some minor variance of results, even the rng seed is same
    %for tempIter=1:resampleIterCount
    %for tempIter=1:1
    parfor tempIter=1:resampleIterCount
        %warning('off','all');
        x_resampled = squeeze(x(tempIter,:,:));
        y_resampled = y; %#ok<*PFTUSE>
        
        svm_options2 = struct;
        svm_options2.c = c;
        svm_options2.KFold_num = KFold_num;
        
        temptemp_svm_resample = fun_seqProbSVM_6location_trialLevel_train(x_resampled,y_resampled,svm_options,svm_options2); %#ok<*PFBNS>
        temp_svm_resample{tempIter} = temptemp_svm_resample;
        %warning('on','all');
    end
    temp_svm = temp_svm_resample{1};
    
    
    for tempIter=1:resampleIterCount
        temp_svm = temp_svm_resample{tempIter};
        
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
        
        temp_svm_resample{tempIter} = temp_svm;
        a = 1; %#ok<*NASGU>
    end
    
    
    svm_train_outputs = struct;
    
    svm_train_outputs.coeff_w_power = temp_svm.coeff_w_power;
    svm_train_outputs.boolIndex_location_seq = temp_svm.boolIndex_location_seq;
    svm_train_outputs.seq_range = temp_svm.seq_range;
    
    svm_posterior_lengthx = 0;
    Posterior_2d_n11n_lengthx_mean = 0;
    Posterior_2d_w_lengthx_mean = 0;
    
    for tempIter=1:resampleIterCount
        temp_svm = temp_svm_resample{tempIter};
        svm_posterior_lengthx = svm_posterior_lengthx + temp_svm.svm_posterior_lengthx;
        Posterior_2d_n11n_lengthx_mean = Posterior_2d_n11n_lengthx_mean + temp_svm.Posterior_2d_n11n_lengthx_mean;
        Posterior_2d_w_lengthx_mean = Posterior_2d_w_lengthx_mean + temp_svm.Posterior_2d_w_lengthx_mean;
    end
    svm_posterior_lengthx = svm_posterior_lengthx ./ resampleIterCount;
    Posterior_2d_n11n_lengthx_mean = Posterior_2d_n11n_lengthx_mean ./ resampleIterCount;
    Posterior_2d_w_lengthx_mean = Posterior_2d_w_lengthx_mean ./ resampleIterCount;
    
    svm_train_outputs.svm_posterior_lengthx = svm_posterior_lengthx;
    svm_train_outputs.Posterior_2d_n11n_lengthx_mean = Posterior_2d_n11n_lengthx_mean;
    svm_train_outputs.Posterior_2d_w_lengthx_mean = Posterior_2d_w_lengthx_mean;
    
    svm_train_outputs.temp_svm_resample = temp_svm_resample;
    
    
    if target_length == 1
        svm_train_length1_outputs = svm_train_outputs;
    elseif target_length == 2
        svm_train_length2_outputs = svm_train_outputs;
    elseif target_length == 3
        svm_train_length3_outputs = svm_train_outputs;
    end
end

output = struct;
output.svm_train_length1_outputs = svm_train_length1_outputs;
output.svm_train_length2_outputs = svm_train_length2_outputs;
output.svm_train_length3_outputs = svm_train_length3_outputs;
