%% Initialization


if_compute = 0;

if_fastCompute = 1;


%% Test memoryPrecision
t_test = tic;

if if_fastCompute == 1   
    
    decodingDataSimplified;
    
    seqIndex_valid_resample_allSessions_collapsed;
    temp_F_dff_decisionBin_resample_allSessions_collaped;
    trialValidBoolIndex;
    
    
    trialIndex_bool_memoryCorrect = trialValidBoolIndex;
    allMemoryCorrectBoolIndex = trialValidBoolIndex;
    temp_trialIndex_valid_resample = repmat(seqIndex_valid_resample_allSessions_collapsed,[1,32])';
    seqIndex = seqIndex_valid_resample_allSessions_collapsed;
    seqIndex_valid = seqIndex;

    
    %F_dff_baselinePeriod;
    %F_dff_length1_sample;
    %F_dff_length2_sample;
    %F_dff_length3_sample;
    %F_dff_decisionPeriodA;
    
    temp_dff = squeeze(temp_F_dff_decisionBin_resample_allSessions_collaped(1,:,:));
    F_dff_decisionPeriodA = repmat(temp_dff,[1,1,5]);
    
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
        F_dff_current = F_dff_decisionPeriodA;
    end
    
    temptemp_range = 1:size(F_dff_current,3);
    
    temp_memoryPrecision_trialLevel_crossTime = zeros(length(seqIndex),length(temptemp_range));
    for tempTimeIndex=temptemp_range(1):temptemp_range(end)
        %parfor tempi=temptemp_range(1):temptemp_range(end)
        %for tempi=temptemp_range(1):temptemp_range(1)
        F_dff_currentTime = F_dff_current(:,:,temptemp_range(tempTimeIndex)); %#ok<*PFBNS>
        
        
        if_accelerate = 1;
        
        if if_accelerate == 0
            step_iter = 1;
        elseif if_accelerate == 1
            %step_iter = 4;%time of 'par for': 2,3,4(170),8(170)
            %step_iter = 40;%time of 'for': 4(293),30(75),40(38)
            step_iter = 40;%4,40
        end
                
        
        boolIndex_location_seq_T = boolIndex_location_seq';
        
        memoryPrecision_trialLevel_resampleIter = cell(length(trialIndex_bool_memoryCorrect),1);
        
        %% Sub-step A: Compute memory precision on training set trials (noChoiceCorrect, choiceMemoryCorrect trials)
        temp_trialIndex_bool = allMemoryCorrectBoolIndex;
        
        temp_trialIndex = find(temp_trialIndex_bool==1);
        temp_trialIndex_resample = temp_trialIndex(temp_trialIndex_valid_resample);
        
        for target_length=1:3
            if target_length == 1
                svm_train_lengthx_outputs = svm_train_length1_outputs;
            elseif target_length == 2
                svm_train_lengthx_outputs = svm_train_length2_outputs;
            elseif target_length == 3
                svm_train_lengthx_outputs = svm_train_length3_outputs;
            end
            temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
            
            temp_svm_resample = svm_train_lengthx_outputs.temp_svm_resample;
            for tempIter=1:step_iter:resampleIterCount
                temptemp_svm = temp_svm_resample{tempIter};
                temptemp_trialIndex = temp_trialIndex_resample(tempIter,:);
                
                temptemp_seqIndex = seqIndex(temptemp_trialIndex);
                temptempBoolIndex = ismember(temptemp_seqIndex,temp_range);
                
                temptemp_seqIndex_target_length = temptemp_seqIndex(temptempBoolIndex);
                temptemp_trialIndex_target_length = temptemp_trialIndex(temptempBoolIndex);
                
                temptemp_trialIndex_target_length; %#ok<*VUNUS>
                
                temptemp_Posterior_2d_n11n = []; %#ok<*PFTUSW>
                                
                
                x = F_dff_currentTime(:,temptemp_trialIndex_target_length);
                y = boolIndex_location_seq(:,temptemp_seqIndex_target_length);
                x_T = x';
                y_T = y';
                
                % multi_Posterior_cell
                multi_Posterior_cell = cell(1, KFold_num);
                temp_Mdl_CV_binary = temptemp_svm.temp_Mdl_CV_binary;
                F_dff_T = x_T;
                
                for temploc=1:numFrames
                    for tempk=1:KFold_num
                        tempTrialBoolIndex_fold = temp_Mdl_CV_binary{temploc}.ModelParameters.Generator.UseObsForIter(:,tempk);
                        temp_F_dff_T = F_dff_T(~tempTrialBoolIndex_fold,:); %#ok<*PFBNS>
                        temp_F_dff_T_2d = temp_F_dff_T;
                        [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary{temploc}.Trained{tempk},temp_F_dff_T_2d);
                        %tempPosterior_2 = tempPosterior(:,2);
                        if size(tempPosterior,2) == 1
                            tempPosterior_2 = tempPosterior(:,1);
                        else
                            tempPosterior_2 = tempPosterior(:,2);
                        end
                        multi_Posterior_cell{tempk}(:,temploc) = tempPosterior_2;
                    end
                end
                a = 1;
                
                % Posterior_2d
                Posterior_2d = zeros(size(x,2),numFrames);
                for tempk=1:KFold_num
                    temp_Posterior = multi_Posterior_cell{tempk};
                    for temploc=1:numFrames
                        tempTrialBoolIndex_fold = temp_Mdl_CV_binary{temploc}.ModelParameters.Generator.UseObsForIter(:,tempk);
                        Posterior_2d(~tempTrialBoolIndex_fold,temploc) = temp_Posterior(:,temploc);
                    end
                end
                temptemp_Posterior_2d_n11n = Posterior_2d;
                
                
                temp_boolIndex_location_seq = boolIndex_location_seq(:,temptemp_seqIndex_target_length)';
                temp_p = temptemp_Posterior_2d_n11n;
                temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
                temp_p_production = prod(temp_p,2);
                
                for tempi=1:length(temptemp_trialIndex_target_length)
                    memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                        [memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_p_production(tempi)];
                end
                
            end
        end
        temptempTrialResampleNum = zeros(length(memoryPrecision_trialLevel_resampleIter),1);
        for tempi=1:length(memoryPrecision_trialLevel_resampleIter)
            temptempBoolIndex(tempi) = ~isempty(memoryPrecision_trialLevel_resampleIter{tempi});
            temptempTrialResampleNum(tempi) = length(memoryPrecision_trialLevel_resampleIter{tempi});
        end
        aa1 = sum(seqIndex_valid<=41); %#ok<*NASGU>
        aa2 = sum(temptempTrialResampleNum>0);
        
        a1 = temptempTrialResampleNum(temptempTrialResampleNum>0);
        [M,I] = sort(a1); %#ok<*ASGLU>
        
        memoryPrecision_trialLevel_resampleMean = nan(length(memoryPrecision_trialLevel_resampleIter),1);
        for tempi=1:length(memoryPrecision_trialLevel_resampleIter)
            if temptempTrialResampleNum(tempi) == 0
                continue
            end
            memoryPrecision_trialLevel_resampleMean(tempi) = mean(memoryPrecision_trialLevel_resampleIter{tempi});
        end
        a1 = sum(~isnan(memoryPrecision_trialLevel_resampleMean));
        memoryPrecision_trialLevel_resampleMean;
        
        
        temp_memoryPrecision_trialLevel_crossTime(:,tempi) = memoryPrecision_trialLevel_resampleMean;

    end
    
    if loopCount==0
        memoryPrecision_trialLevel_crossTime_baseline = temp_memoryPrecision_trialLevel_crossTime; %#ok<*NASGU>
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


%% End