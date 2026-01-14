function testMemoryPrecision_output = fun_testMemoryPrecision_currentTime_v10(F_dff_currentTime,options_testMemoryPrecision)
%% Initialization


trialIndex_bool_memoryCorrect = options_testMemoryPrecision.trialIndex_bool_memoryCorrect;
allMemoryCorrectBoolIndex = options_testMemoryPrecision.allMemoryCorrectBoolIndex;
temp_trialIndex_valid_resample = options_testMemoryPrecision.temp_trialIndex_valid_resample;
svm_train_length1_outputs = options_testMemoryPrecision.svm_train_length1_outputs;
svm_train_length2_outputs = options_testMemoryPrecision.svm_train_length2_outputs;
svm_train_length3_outputs = options_testMemoryPrecision.svm_train_length3_outputs;
numSeq = options_testMemoryPrecision.numSeq;
resampleIterCount = options_testMemoryPrecision.resampleIterCount;
seqIndex = options_testMemoryPrecision.seqIndex;
boolIndex_location_seq = options_testMemoryPrecision.boolIndex_location_seq;
seqIndex_valid = options_testMemoryPrecision.seqIndex_valid;

allMemoryErrorBoolIndex = options_testMemoryPrecision.allMemoryErrorBoolIndex;
seqIndex_response = options_testMemoryPrecision.seqIndex_response;
KFold_num = options_testMemoryPrecision.KFold_num;
numFrames = options_testMemoryPrecision.numFrames;
if_memeoryPrecision_stimuli0_response1 = options_testMemoryPrecision.if_memeoryPrecision_stimuli0_response1;
if_inferOffloadResponse = options_testMemoryPrecision.if_inferOffloadResponse;

choiceOffloadBoolIndex = options_testMemoryPrecision.choiceOffloadBoolIndex;

if_compute = options_testMemoryPrecision.if_compute;
if_fastCompute = options_testMemoryPrecision.if_fastCompute;
if_memoryPrecision_accuracy0_sigma1 = options_testMemoryPrecision.if_memoryPrecision_accuracy0_sigma1;
score_stimuli_to_response = options_testMemoryPrecision.score_stimuli_to_response;
if_precision_meanProb0_sumProb1 = options_testMemoryPrecision.if_precision_meanProb0_sumProb1;
if_entropy = options_testMemoryPrecision.if_entropy;
fun_sigmaBased_singleTrialPrecision = options_testMemoryPrecision.fun_sigmaBased_singleTrialPrecision;


if_accelerate = 1;%1

if if_accelerate == 0
    step_iter = 1;
    
    if_testKfold = 1;

elseif if_accelerate == 1
    %step_iter = 4;%time of 'par for': 2,3,4(170),8(170)
    %step_iter = 40;%time of 'for': 4(293),30(75),40(38)
    
    
    if if_fastCompute == 1
        step_iter = 4;
        %step_iter = 40;        
    end
    if if_compute == 1
        step_iter = 40;
        %step_iter = 2;
    end    
    
    %step_iter = 4;%4,40    
    
    if_testKfold = 0;
end



boolIndex_location_seq_T = boolIndex_location_seq';

memoryPrecision_trialLevel_resampleIter = cell(length(trialIndex_bool_memoryCorrect),1);
locCorr_trialLevel_resampleIter = cell(length(trialIndex_bool_memoryCorrect),1);

memoryPrecision_trialLevel_resampleIter_matrix = nan(length(trialIndex_bool_memoryCorrect),resampleIterCount);
locCorr_trialLevel_resampleIter_matrix = nan(length(trialIndex_bool_memoryCorrect),resampleIterCount);

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
    %for tempIter=1:resampleIterCount
    for tempIter=1:step_iter:resampleIterCount
    %tempIterA=1:step_iter:resampleIterCount;
    %for tempIterB=1:length(tempIterA)
    %    tempIter = tempIterA(tempIterB);
        temptemp_svm = temp_svm_resample{tempIter};
        temptemp_trialIndex = temp_trialIndex_resample(tempIter,:);
        
        temptemp_seqIndex = seqIndex(temptemp_trialIndex);
        temptempBoolIndex = ismember(temptemp_seqIndex,temp_range);
        
        temptemp_seqIndex_target_length = temptemp_seqIndex(temptempBoolIndex);
        temptemp_trialIndex_target_length = temptemp_trialIndex(temptempBoolIndex);
        
        temptemp_trialIndex_target_length; %#ok<*VUNUS>
        %temptemp_Posterior_2d_n11n = temptemp_svm.Posterior_2d_n11n;
        
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
        
        
            temp_precision = nan(length(temptemp_trialIndex_target_length),1);
            parfor tempi=1:length(temptemp_trialIndex_target_length)
                
                temp_locDistri = temptemp_Posterior_2d_n11n(tempi,:);
                tempSeqIndex = temptemp_seqIndex_target_length(tempi);
                valid_length = 3;
                
                temp_score = score_stimuli_to_response(tempSeqIndex,1:sum(numSeq(1:valid_length)));
                
                temp_options = struct;
                temp_options.boolIndex_location_seq_T = boolIndex_location_seq_T;
                temp_options.numSeq = numSeq;
                temp_options.valid_length = valid_length;
                temp_options.if_precision_meanProb0_sumProb1 = if_precision_meanProb0_sumProb1;
                
                temp_options.temp_locDistri = temp_locDistri;
                temp_options.temp_score = temp_score;
                temp_options.if_entropy = if_entropy;
                temp_options.tempSeqIndex = tempSeqIndex;
                temp_options.tempSeqLength = target_length;
                
                temp_precision(tempi) = fun_sigmaBased_singleTrialPrecision(temp_options); %#ok<*PFOUS,*PFBNS>                
            end        
        
        
        temp_r_loc = nan(size(temptemp_Posterior_2d_n11n,1),1);
        for tempi=1:size(temptemp_Posterior_2d_n11n,1)
            [temp_r_loc(tempi),~] = corr(temptemp_Posterior_2d_n11n(tempi,:)',temp_boolIndex_location_seq(tempi,:)');
        end
        for tempi=1:length(temptemp_trialIndex_target_length)
            locCorr_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                [locCorr_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_r_loc(tempi)];      
            
            locCorr_trialLevel_resampleIter_matrix(temptemp_trialIndex_target_length(tempi),tempIter) = temp_r_loc(tempi);
        end        
        
        for tempi=1:length(temptemp_trialIndex_target_length)
            if if_memoryPrecision_accuracy0_sigma1 == 0
                memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_p_production(tempi)];
                
                memoryPrecision_trialLevel_resampleIter_matrix(temptemp_trialIndex_target_length(tempi),tempIter) = temp_p_production(tempi);
            elseif if_memoryPrecision_accuracy0_sigma1 == 1
                memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_precision(tempi)];
                
                memoryPrecision_trialLevel_resampleIter_matrix(temptemp_trialIndex_target_length(tempi),tempIter) = temp_precision(tempi);                
            end
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

%memoryPrecision_trialLevel = nan(length(memoryPrecision_trialLevel_resampleIter),1);
memoryPrecision_trialLevel_resampleMean = nan(length(memoryPrecision_trialLevel_resampleIter),1);
for tempi=1:length(memoryPrecision_trialLevel_resampleIter)
    if temptempTrialResampleNum(tempi) == 0
        continue
    end
    memoryPrecision_trialLevel_resampleMean(tempi) = mean(memoryPrecision_trialLevel_resampleIter{tempi});
end
a1 = sum(~isnan(memoryPrecision_trialLevel_resampleMean));
memoryPrecision_trialLevel_resampleMean;

%% Sub-step B: Compute memory precision on noChoiceError and choiceMemoryError trials
temp_trialIndex_bool = allMemoryErrorBoolIndex;

temp_trialIndex = find(temp_trialIndex_bool==1);
for target_length=1:3
    if target_length == 1
        svm_train_lengthx_outputs = svm_train_length1_outputs;
    elseif target_length == 2
        svm_train_lengthx_outputs = svm_train_length2_outputs;
    elseif target_length == 3
        svm_train_lengthx_outputs = svm_train_length3_outputs;
    end
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    
    temptemp_trialIndex = temp_trialIndex;
    
    if if_memeoryPrecision_stimuli0_response1 == 0
        temptemp_seqIndex = seqIndex(temptemp_trialIndex);
    elseif if_memeoryPrecision_stimuli0_response1 == 1
        temptemp_seqIndex = seqIndex_response(temptemp_trialIndex);
    end
    
    temptempBoolIndex = ismember(temptemp_seqIndex,temp_range);
    %temptempBoolIndex = true(1,length(temptemp_seqIndex));
    
    temptemp_seqIndex_target_length = temptemp_seqIndex(temptempBoolIndex);
    temptemp_trialIndex_target_length = temptemp_trialIndex(temptempBoolIndex);
    
    x = F_dff_currentTime(:,temptemp_trialIndex_target_length);
    %y = boolIndex_location_allTrial(:,temptemp_trialIndex_target_length);
    y = boolIndex_location_seq(:,temptemp_seqIndex_target_length);
    
    x_T = x';
    y_T = y';
    
    temp_svm_resample = svm_train_lengthx_outputs.temp_svm_resample;
    
    %for tempIter=1:resampleIterCount
    for tempIter=1:step_iter:resampleIterCount
        temptemp_svm = temp_svm_resample{tempIter};
        
        if if_testKfold == 1
            temp_Mdl_CV_binary = temptemp_svm.temp_Mdl_CV_binary;
            
            temptemp_Posterior_2d_kfold = zeros(KFold_num,length(temptemp_trialIndex_target_length),numFrames);
            for temploc=1:numFrames
                for tempk=1:KFold_num
                    [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary{temploc}.Trained{tempk},x_T);% Very time-consuming!!!
                    if size(tempPosterior,2) == 1
                        tempPosterior_2 = tempPosterior(:,1);
                    else
                        tempPosterior_2 = tempPosterior(:,2);
                    end
                    temptemp_Posterior_2d_kfold(tempk,:,temploc) = tempPosterior_2;
                end
            end
            temptemp_Posterior_2d_kfoldMean = squeeze(mean(temptemp_Posterior_2d_kfold,1));
            temptemp_Posterior_2d_n11n = temptemp_Posterior_2d_kfoldMean;
            
        elseif if_testKfold == 0           
            temp_Mdl_binary = temptemp_svm.temp_Mdl_binary;
            
            temptemp_Posterior_2d = zeros(length(temptemp_trialIndex_target_length),numFrames);
            for temploc=1:numFrames
                [~,~,~,tempPosterior] = predict(temp_Mdl_binary{temploc},x_T);% Very time-consuming!!!
                if size(tempPosterior,2) == 1
                    tempPosterior_2 = tempPosterior(:,1);
                else
                    tempPosterior_2 = tempPosterior(:,2);
                end
                temptemp_Posterior_2d(:,temploc) = tempPosterior_2;
            end
            temptemp_Posterior_2d_n11n = temptemp_Posterior_2d;
            
            if false
                temp_r = nan(size(a1,1),1);
                for temptempi=1:size(a1,1) %#ok<*UNRCH>
                    [temp_r(temptempi),~] = corr(a2(temptempi,:)',a1(temptempi,:)');
                end
            end
        end
        
        temp_boolIndex_location_seq = boolIndex_location_seq(:,temptemp_seqIndex_target_length)';
        temp_p = temptemp_Posterior_2d_n11n;
        temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
        temp_p_production = prod(temp_p,2);
        
        
            temp_precision = nan(length(temptemp_trialIndex_target_length),1);
            parfor tempi=1:length(temptemp_trialIndex_target_length)
                
                temp_locDistri = temptemp_Posterior_2d_n11n(tempi,:);
                tempSeqIndex = temptemp_seqIndex_target_length(tempi);
                valid_length = 3;
                
                temp_score = score_stimuli_to_response(tempSeqIndex,1:sum(numSeq(1:valid_length)));
                
                temp_options = struct;
                temp_options.boolIndex_location_seq_T = boolIndex_location_seq_T;
                temp_options.numSeq = numSeq;
                temp_options.valid_length = valid_length;
                temp_options.if_precision_meanProb0_sumProb1 = if_precision_meanProb0_sumProb1;
                
                temp_options.temp_locDistri = temp_locDistri;
                temp_options.temp_score = temp_score;
                temp_options.if_entropy = if_entropy;
                temp_options.tempSeqIndex = tempSeqIndex;
                temp_options.tempSeqLength = target_length;
                
                
                temp_precision(tempi) = fun_sigmaBased_singleTrialPrecision(temp_options); %#ok<*PFOUS,*PFBNS>                
            end   
            
            
        temp_r_loc = nan(size(temptemp_Posterior_2d_n11n,1),1);
        for tempi=1:size(temptemp_Posterior_2d_n11n,1)
            [temp_r_loc(tempi),~] = corr(temptemp_Posterior_2d_n11n(tempi,:)',temp_boolIndex_location_seq(tempi,:)');
        end
        for tempi=1:length(temptemp_trialIndex_target_length)
            locCorr_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                [locCorr_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_r_loc(tempi)];      
            
            locCorr_trialLevel_resampleIter_matrix(temptemp_trialIndex_target_length(tempi),tempIter) = temp_r_loc(tempi);
        end        
        
        for tempi=1:length(temptemp_trialIndex_target_length)
            if if_memoryPrecision_accuracy0_sigma1 == 0
                memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_p_production(tempi)];
                
                memoryPrecision_trialLevel_resampleIter_matrix(temptemp_trialIndex_target_length(tempi),tempIter) = temp_p_production(tempi);
            elseif if_memoryPrecision_accuracy0_sigma1 == 1
                memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_precision(tempi)];
                
                memoryPrecision_trialLevel_resampleIter_matrix(temptemp_trialIndex_target_length(tempi),tempIter) = temp_precision(tempi);                
            end
        end

        
    end
end

%% Sub-step C: Compute memory precision on choiceOffload trials

temp_trialIndex_bool = choiceOffloadBoolIndex;
temp_trialIndex = find(temp_trialIndex_bool==1);

temptemp_Posterior_2d_resampleMean_length123 = zeros(length(temp_trialIndex_bool),numFrames);

for target_length=1:3
    if target_length == 1
        svm_train_lengthx_outputs = svm_train_length1_outputs;
    elseif target_length == 2
        svm_train_lengthx_outputs = svm_train_length2_outputs;
    elseif target_length == 3
        svm_train_lengthx_outputs = svm_train_length3_outputs;
    end
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    
    temptemp_trialIndex = temp_trialIndex;
    
    %if if_memeoryPrecision_stimuli0_response1 == 0
    %   temptemp_seqIndex = seqIndex(temptemp_trialIndex);
    %elseif if_memeoryPrecision_stimuli0_response1 == 1
    %   temptemp_seqIndex = seqIndex_response(temptemp_trialIndex);
    %end
    temptemp_seqIndex = seqIndex(temptemp_trialIndex);
    
    temptempBoolIndex = ismember(temptemp_seqIndex,temp_range);
    %temptempBoolIndex = true(1,length(temptemp_seqIndex));
    
    temptemp_seqIndex_target_length = temptemp_seqIndex(temptempBoolIndex);
    temptemp_trialIndex_target_length = temptemp_trialIndex(temptempBoolIndex);
    
    x = F_dff_currentTime(:,temptemp_trialIndex_target_length);
    %y = boolIndex_location_allTrial(:,temptemp_trialIndex_target_length);
    y = boolIndex_location_seq(:,temptemp_seqIndex_target_length);
    
    x_T = x';
    y_T = y';
    
    temp_svm_resample = svm_train_lengthx_outputs.temp_svm_resample;
    
    %temptemp_Posterior_2d_iter = zeros(resampleIterCount,length(temptemp_trialIndex_target_length),numFrames);
    temptemp_Posterior_2d_iter = nan(resampleIterCount,length(temptemp_trialIndex_target_length),numFrames);
    
    %for tempIter=1:resampleIterCount
    for tempIter=1:step_iter:resampleIterCount
        temptemp_svm = temp_svm_resample{tempIter};
        
        if if_testKfold == 1
            temp_Mdl_CV_binary = temptemp_svm.temp_Mdl_CV_binary;
            
            temptemp_Posterior_2d_kfold = zeros(KFold_num,length(temptemp_trialIndex_target_length),numFrames);
            for temploc=1:numFrames
                for tempk=1:KFold_num
                    [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary{temploc}.Trained{tempk},x_T);% Very time-consuming!!!
                    if size(tempPosterior,2) == 1
                        tempPosterior_2 = tempPosterior(:,1);
                    else
                        tempPosterior_2 = tempPosterior(:,2);
                    end
                    temptemp_Posterior_2d_kfold(tempk,:,temploc) = tempPosterior_2;
                end
            end
            temptemp_Posterior_2d_kfoldMean = squeeze(mean(temptemp_Posterior_2d_kfold,1));
            temptemp_Posterior_2d_n11n = temptemp_Posterior_2d_kfoldMean;
            
        elseif if_testKfold == 0
            temp_Mdl_binary = temptemp_svm.temp_Mdl_binary;
            
            temptemp_Posterior_2d = zeros(length(temptemp_trialIndex_target_length),numFrames);
            for temploc=1:numFrames
                [~,~,~,tempPosterior] = predict(temp_Mdl_binary{temploc},x_T);% Very time-consuming!!!
                if size(tempPosterior,2) == 1
                    tempPosterior_2 = tempPosterior(:,1);
                else
                    tempPosterior_2 = tempPosterior(:,2);
                end
                temptemp_Posterior_2d(:,temploc) = tempPosterior_2;                            
            end
            
            temptemp_Posterior_2d_n11n = temptemp_Posterior_2d;            
            
            
            if false
                temp_r = nan(size(a1,1),1);
                for temptempi=1:size(a1,1) %#ok<*UNRCH>
                    [temp_r(temptempi),~] = corr(a2(temptempi,:)',a1(temptempi,:)');
                end
            end
            
        end
        
        %temptemp_Posterior_2d_iter(tempIter,:,:) = temptemp_Posterior_2d_kfoldMean;
        temptemp_Posterior_2d_iter(tempIter,:,:) = temptemp_Posterior_2d_n11n;
        
        if if_memeoryPrecision_stimuli0_response1 == 0
            temp_boolIndex_location_seq = boolIndex_location_seq(:,temptemp_seqIndex_target_length)';
            
        elseif if_memeoryPrecision_stimuli0_response1 == 1
            if if_inferOffloadResponse == 1
                % infer response seq from temptemp_Posterior_2d_n11n
                [M,I] = sort(temptemp_Posterior_2d_n11n,2,'descend');
                temp_boolIndex_location_seq = false(size(temptemp_Posterior_2d_n11n,1),numFrames);
                for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                    temp_boolIndex_location_seq(tempi,I(tempi,1:target_length)) = true;
                end
            elseif if_inferOffloadResponse == 0
                temp_boolIndex_location_seq = boolIndex_location_seq(:,temptemp_seqIndex_target_length)';
            end
        end
        
        
        temp_p = temptemp_Posterior_2d_n11n;
        temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
        temp_p_production = prod(temp_p,2);
        
        
            temp_precision = nan(length(temptemp_trialIndex_target_length),1);
            parfor tempi=1:length(temptemp_trialIndex_target_length)
                
                temp_locDistri = temptemp_Posterior_2d_n11n(tempi,:);
                tempSeqIndex = temptemp_seqIndex_target_length(tempi);
                valid_length = 3;
                
                temp_score = score_stimuli_to_response(tempSeqIndex,1:sum(numSeq(1:valid_length)));
                
                temp_options = struct;
                temp_options.boolIndex_location_seq_T = boolIndex_location_seq_T;
                temp_options.numSeq = numSeq;
                temp_options.valid_length = valid_length;
                temp_options.if_precision_meanProb0_sumProb1 = if_precision_meanProb0_sumProb1;
                
                temp_options.temp_locDistri = temp_locDistri;
                temp_options.temp_score = temp_score;
                temp_options.if_entropy = if_entropy;
                temp_options.tempSeqIndex = tempSeqIndex;
                temp_options.tempSeqLength = target_length;
                
                
                temp_precision(tempi) = fun_sigmaBased_singleTrialPrecision(temp_options); %#ok<*PFOUS,*PFBNS>                
            end   
            
            
        temp_r_loc = nan(size(temptemp_Posterior_2d_n11n,1),1);
        for tempi=1:size(temptemp_Posterior_2d_n11n,1)
            [temp_r_loc(tempi),~] = corr(temptemp_Posterior_2d_n11n(tempi,:)',temp_boolIndex_location_seq(tempi,:)');
        end
        for tempi=1:length(temptemp_trialIndex_target_length)
            locCorr_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                [locCorr_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_r_loc(tempi)];      
            
            locCorr_trialLevel_resampleIter_matrix(temptemp_trialIndex_target_length(tempi),tempIter) = temp_r_loc(tempi);
        end        
        
        for tempi=1:length(temptemp_trialIndex_target_length)
            if if_memoryPrecision_accuracy0_sigma1 == 0
                memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_p_production(tempi)];
                
                memoryPrecision_trialLevel_resampleIter_matrix(temptemp_trialIndex_target_length(tempi),tempIter) = temp_p_production(tempi);
            elseif if_memoryPrecision_accuracy0_sigma1 == 1
                memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_precision(tempi)];
                
                memoryPrecision_trialLevel_resampleIter_matrix(temptemp_trialIndex_target_length(tempi),tempIter) = temp_precision(tempi);                
            end
        end
        
    end
    %temptemp_Posterior_2d_resampleMean = squeeze(mean(temptemp_Posterior_2d_iter,1));
    temptemp_Posterior_2d_resampleMean = squeeze(mean(temptemp_Posterior_2d_iter,1,'omitnan'));
    
    temptemp_Posterior_2d_resampleMean_length123(temptemp_trialIndex_target_length,:) = temptemp_Posterior_2d_resampleMean;
end
%temptemp_Posterior_2d_resampleMean_length123_B = ...
%    temptemp_Posterior_2d_resampleMean_length123(temp_trialIndex_bool&(seqIndex<=41),:);

if if_memeoryPrecision_stimuli0_response1 == 1
    for target_length=1:3
        temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
        
        temptempBoolIndex = ismember(seqIndex,temp_range);
        temptempBoolIndex2 = temptempBoolIndex & choiceOffloadBoolIndex;
        
        if if_inferOffloadResponse == 1
            temptemp_Posterior_2d_n11n = temptemp_Posterior_2d_resampleMean_length123(temptempBoolIndex2,:);
            
            % infer response seq from temptemp_Posterior_2d_n11n
            [M,I] = sort(temptemp_Posterior_2d_n11n,2,'descend');
            temp_boolIndex_location_seq = false(size(temptemp_Posterior_2d_n11n,1),numFrames);
            for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                temp_boolIndex_location_seq(tempi,I(tempi,1:target_length)) = true;
            end
            
            temptempSeqIndex_response = nan(1,size(temp_boolIndex_location_seq,1));
            for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                temptempBoolIndex_seq = temp_boolIndex_location_seq(tempi,:);
                for tempj=1:sum(numSeq(1:3))
                    temptempBoolIndex_seq2 = boolIndex_location_seq_T(tempj,:);
                    if sum(temptempBoolIndex_seq==temptempBoolIndex_seq2) == numFrames
                        temptempSeqIndex_response(tempi) = tempj;
                        break
                    end
                end
            end
            temptempSeqIndex_response;
            seqIndex_response(temptempBoolIndex2) = temptempSeqIndex_response;
        elseif if_inferOffloadResponse == 0
            seqIndex_response(temptempBoolIndex2) = seqIndex(temptempBoolIndex2); %#ok<*SAGROW>
        end
        
    end
    
end


temptempTrialResampleNum = zeros(length(memoryPrecision_trialLevel_resampleIter),1);
for tempi=1:length(memoryPrecision_trialLevel_resampleIter)
    temptempBoolIndex(tempi) = ~isempty(memoryPrecision_trialLevel_resampleIter{tempi});
    temptempTrialResampleNum(tempi) = length(memoryPrecision_trialLevel_resampleIter{tempi});
end
aa1 = sum(seqIndex<=41);
aa2 = sum(temptempTrialResampleNum>0);

a1 = temptempTrialResampleNum(temptempTrialResampleNum>0);
[M,I] = sort(a1);

memoryPrecision_trialLevel_resampleMean = nan(length(memoryPrecision_trialLevel_resampleIter),1);
for tempi=1:length(memoryPrecision_trialLevel_resampleIter)
    if temptempTrialResampleNum(tempi) == 0
        continue
    end
    memoryPrecision_trialLevel_resampleMean(tempi) = mean(memoryPrecision_trialLevel_resampleIter{tempi});
end
% a1 = sum(~isnan(memoryPrecision_trialLevel_resampleMean));
% memoryPrecision_trialLevel_resampleMean;
% 
% a1 = sum(~isnan(memoryPrecision_trialLevel_resampleMean));
% a2 = sum(seqIndex<=41);
% fprintf('valid trial num of memoryPrecision is %d, theretical trial num (stimuli) is %d.\n',a1,a2);

memoryPrecision_trialLevel = memoryPrecision_trialLevel_resampleMean;


locCorr_trialLevel_resampleMean = nan(length(locCorr_trialLevel_resampleIter),1);
for tempi=1:length(locCorr_trialLevel_resampleIter)
    if temptempTrialResampleNum(tempi) == 0
        continue
    end
    locCorr_trialLevel_resampleMean(tempi) = mean(locCorr_trialLevel_resampleIter{tempi});
end
locCorr_trialLevel = locCorr_trialLevel_resampleMean;



locCorr_trialLevel_resampleIter_matrix;
memoryPrecision_trialLevel_resampleIter_matrix;


%%
testMemoryPrecision_output = struct;
testMemoryPrecision_output.memoryPrecision_trialLevel_currentTime = memoryPrecision_trialLevel;
testMemoryPrecision_output.locCorr_trialLevel_currentTime = locCorr_trialLevel;

testMemoryPrecision_output.memoryPrecision_trialLevel_resampleIter_matrix_currentTime = memoryPrecision_trialLevel_resampleIter_matrix;
testMemoryPrecision_output.locCorr_trialLevel_resampleIter_matrix_currentTime = locCorr_trialLevel_resampleIter_matrix;


a = 1;

%% End