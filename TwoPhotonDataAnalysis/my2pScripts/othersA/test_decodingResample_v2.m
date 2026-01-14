if_resample = 1;

resampleTrialCount = 10;%6-->10
resampleTrialCount_min = 4;
resampleIterCount = 5;%2-->10



svm_options = struct;
svm_options.boolIndex_location_seq = boolIndex_location_seq;
svm_options.t_decoder = t_decoder;
svm_options.numFrames = numFrames;
svm_options.targetPATH = targetPATH;


for target_length=1:3
% for target_length=1:1
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    gAcc = gAcc_target_collapsed_inOne(temp_range);
    
    svm_options.gAcc = gAcc;
    svm_options.seq_range = temp_range;
    svm_options.target_length = target_length;
    
    tempBoolIndex = ismember(seqIndex_valid,temp_range);
    x = temp_F_dff_decisionBin(:,tempBoolIndex);
    y = boolIndex_location_trial(:,tempBoolIndex);
    
    
    if if_resample == 1
        temp_uniqueSeq = unique(seqIndex_valid(tempBoolIndex));
        validSeqBoolIndex = ismember(temp_range,temp_uniqueSeq);
        
        temp_seqCount = zeros(1,numSeq(target_length));
        for tempi=1:numSeq(target_length)
            tempj = temp_range(tempi);
            temp_seqCount(tempi) = sum(seqIndex_valid==tempj,'all');
        end
        
        validSeqBoolIndex2 = temp_seqCount>=resampleTrialCount_min;
        temp_uniqueSeq2 = temp_range(validSeqBoolIndex2);
        
        temp_svm_resample = cell(1,tempIter);
        for tempIter=1:resampleIterCount
            %parfor tempIter=1:resampleIterCount
            
            tempIndex_resample = zeros(resampleTrialCount,sum(validSeqBoolIndex2));
            tempj = 0;
            for tempi=1:numSeq(target_length)
                if validSeqBoolIndex2(tempi) == false
                    continue
                end
                tempj = tempj + 1;
                tempIndex = find(seqIndex_valid==temp_uniqueSeq2(tempj));
                
                temp_length = length(tempIndex);
                if temp_length >= resampleTrialCount
                    tempI = sort(randperm(temp_length,resampleTrialCount));
                elseif temp_length < resampleTrialCount
                    temp_intLoopCount = floor(resampleTrialCount/temp_length);
                    temp_resampleTrialCount = resampleTrialCount-temp_length*temp_intLoopCount;
                    tempI1 = [];
                    for temptempi=1:temp_intLoopCount
                        tempI1 = [tempI1 (1:temp_length)]; %#ok<*AGROW>
                    end
                    tempI2 = sort(randperm(temp_length,temp_resampleTrialCount));
                    
                    %tempI = sort([tempI1 tempI2]);
                    tempI = [tempI1 tempI2]; %#ok<*NASGU>
                end
                
                tempIndex_resample(:,tempj) = tempIndex(tempI);
            end
            tempIndex_resample_1d = reshape(tempIndex_resample,1,[]);
            
            x_resample = temp_F_dff_decisionBin(:,tempIndex_resample_1d);
            y_resample = boolIndex_location_trial(:,tempIndex_resample_1d);
            
            x = x_resample;
            y = y_resample;
            
            
            temptemp_svm_resample = fun_seqProbSVM_6location_trialLevel_train(x,y,svm_options);
            temp_svm_resample{tempIter} = temptemp_svm_resample;
            
        end
        a = 1;
        
        
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
            a = 1;
        end
        a = 1;
        
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
        
        a = 1;
        
    elseif if_resample == 0
        svm_train_outputs = fun_seqProbSVM_6location_trialLevel_train(x,y,svm_options);
        
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
        
    end
    
    
    
    if target_length == 1
        svm_train_length1_outputs = svm_train_outputs;
    elseif target_length == 2
        svm_train_length2_outputs = svm_train_outputs;
    elseif target_length == 3
        svm_train_length3_outputs = svm_train_outputs;
    end
    
end
a = 1;