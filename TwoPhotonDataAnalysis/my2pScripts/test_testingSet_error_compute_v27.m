% Chuan's 20th script (20251214)
% This script: To test location memory decoder in error trial.
%% Initialization

if_compute_test = 1;

errorTrial_minNum = 1;

if_stimuliError0_responseError1 = 1;
if_decoder_decodingAcc0_pProd1 = 0;%1


if_n11n = 0;%0



%% Compute
if if_compute_test == 1
    t_test = tic;
    for target_length=1:3
        
        if target_length == 1
            svm_train_lengthx_outputs = svm_train_length1_outputs;
        elseif target_length == 2
            svm_train_lengthx_outputs = svm_train_length2_outputs;
        elseif target_length == 3
            svm_train_lengthx_outputs = svm_train_length3_outputs;
        end
        
        % use target length trial
        %temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
        %tempBoolIndex_targetLength = ismember(seqIndex,temp_range);
        %tempBoolIndex = tempBoolIndex_targetLength & trialIndex_bool_memoryError;
        
        % use all length trial
        tempBoolIndex = trialIndex_bool_memoryError;
        
        trialIndex_bool_memoryError_targetLength = tempBoolIndex;
        
        x = F_dff_decisionBin1(:,tempBoolIndex);
        y = boolIndex_location_allTrial(:,tempBoolIndex);
        y2 = boolIndex_location_allTrial_response(:,tempBoolIndex);
        
        x_T = x';
        y_T = y';
        y2_T = y2';
        
        temp_trialNum = size(x,2);
        
        Posterior_2d_kfoldMean_resample = zeros(resampleIterCount,temp_trialNum,numFrames);
        for tempIter=1:resampleIterCount
            %parfor tempIter=1:resampleIterCount
            temptemp_svm = svm_train_lengthx_outputs.temp_svm_resample{tempIter}; %#ok<*PFBNS>
            temp_Mdl_CV_binary = temptemp_svm.temp_Mdl_CV_binary;
            
            Posterior_2d_kfold = zeros(KFold_num,temp_trialNum,numFrames);
            for temploc=1:numFrames
                for tempk=1:KFold_num
                    [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary{temploc}.Trained{tempk},x_T);% Very time-consuming!!!
                    %[~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary{temploc}.Trained{tempk},x_T,'Options',statset('UseParallel',true));% Even more time-consuming!!!
                    if size(tempPosterior,2) == 1
                        tempPosterior_2 = tempPosterior(:,1);
                    else
                        tempPosterior_2 = tempPosterior(:,2);
                    end
                    Posterior_2d_kfold(tempk,:,temploc) = tempPosterior_2;
                end
            end
            Posterior_2d_kfoldMean = squeeze(mean(Posterior_2d_kfold,1));
            Posterior_2d_kfoldMean_resample(tempIter,:,:) = Posterior_2d_kfoldMean;
        end
        a = 1; %#ok<*NASGU>
        Posterior_2d_kfoldMean_resampleMean = squeeze(mean(Posterior_2d_kfoldMean_resample,1));
        Posterior_2d_resampleMean = Posterior_2d_kfoldMean_resampleMean;
        
        Posterior_2d_resample = Posterior_2d_kfoldMean_resample;
        
        test_errorTrial = struct;
        test_errorTrial.Posterior_2d_resampleMean = Posterior_2d_resampleMean;
        test_errorTrial.temp_trialIndex_bool = trialIndex_bool_memoryError_targetLength;
        test_errorTrial.stimuliLabel = y_T;
        test_errorTrial.responseLabel = y2_T;
        
        test_errorTrial.Posterior_2d_resample = Posterior_2d_resample;
        
        svm_train_lengthx_outputs.test_errorTrial = test_errorTrial;
        
        if target_length == 1
            svm_train_length1_outputs = svm_train_lengthx_outputs;
        elseif target_length == 2
            svm_train_length2_outputs = svm_train_lengthx_outputs;
        elseif target_length == 3
            svm_train_length3_outputs = svm_train_lengthx_outputs;
        end
    end
    
    fprintf('t_test = %.1f secs.\n',toc(t_test));
    
end
a = 1;

temp_posterior = svm_train_length1_outputs.test_errorTrial.Posterior_2d_resampleMean;
temp_posterior_n11n = temp_posterior./sum(temp_posterior,2);

a1 = sum(temp_posterior > 0.5,2);
sum(a1>0);

a2 = sum(temp_posterior_n11n > 0.5,2);
sum(a2>0);


errorTrialNum = sum(trialIndex_bool_memoryError);

% %% Further compute
%% Distinguish posterior between stimuli label and response label, because different length use different decoder
seqCount_allLength_memoryError;

posterior_stimuliLength = zeros(errorTrialNum,numFrames);
posterior_responseLength = zeros(errorTrialNum,numFrames);
stimuliLabel_errorTrial = svm_train_length1_outputs.test_errorTrial.stimuliLabel;
responseLabel_errorTrial = svm_train_length1_outputs.test_errorTrial.responseLabel;
for target_length=1:3
    
    if target_length == 1
        svm_train_lengthx_outputs = svm_train_length1_outputs;
    elseif target_length == 2
        svm_train_lengthx_outputs = svm_train_length2_outputs;
    elseif target_length == 3
        svm_train_lengthx_outputs = svm_train_length3_outputs;
    end
    
    temp_posterior = svm_train_lengthx_outputs.test_errorTrial.Posterior_2d_resampleMean;
    temp_seqLength_stimuli = sum(stimuliLabel_errorTrial,2);
    temp_seqLength_response = sum(responseLabel_errorTrial,2);
    
    temptempBoolIndex1 = temp_seqLength_stimuli == target_length;
    temptempBoolIndex2 = temp_seqLength_response == target_length;
    
    posterior_stimuliLength(temptempBoolIndex1,:) = temp_posterior(temptempBoolIndex1,:);
    posterior_responseLength(temptempBoolIndex2,:) = temp_posterior(temptempBoolIndex2,:);
    
end
posterior_stimuliLength;
posterior_responseLength;
stimuliLabel_errorTrial;
responseLabel_errorTrial;

posterior_stimuliLength_raw = posterior_stimuliLength;
posterior_responseLength_raw = posterior_responseLength;

% if_n11n = 1;%0

%% to normalize posterior_stimuliLength and posterior_responseLength
if if_n11n == 1
    %% posterior_stimuliLength
    temp_seqLength_stimuli = sum(stimuliLabel_errorTrial,2);
    for target_length=1:3
        temptempBoolIndex1 = temp_seqLength_stimuli == target_length;
        temp_posterior = posterior_stimuliLength(temptempBoolIndex1,:);
        
        temp_posterior_n11n = (temp_posterior./sum(temp_posterior,2))*target_length;
        a = 0;
        while true
            tempindex = temp_posterior_n11n>1;
            if sum(tempindex,'all') == 0
                break
            end
            temp_posterior_n11n(temp_posterior_n11n>1) = 1;
            temp_posterior_n11n = (temp_posterior_n11n./sum(temp_posterior_n11n,2))*target_length;
            a = a + 1;
        end
        a1=sum(temp_posterior,2);
        a2=sum(temp_posterior_n11n,2);
        posterior_stimuliLength(temptempBoolIndex1,:) = temp_posterior_n11n;
    end
    a3_raw = sum(posterior_stimuliLength_raw,2);
    a3 = sum(posterior_stimuliLength,2);
    a = 1;
    
    %% posterior_responseLength
    temp_seqLength_response = sum(responseLabel_errorTrial,2);
    for target_length=1:3
        temptempBoolIndex1 = temp_seqLength_response == target_length;
        temp_posterior = posterior_responseLength(temptempBoolIndex1,:);
        
        temp_posterior_n11n = (temp_posterior./sum(temp_posterior,2))*target_length;
        a = 0;
        while true
            tempindex = temp_posterior_n11n>1;
            if sum(tempindex,'all') == 0
                break
            end
            temp_posterior_n11n(temp_posterior_n11n>1) = 1;
            temp_posterior_n11n = (temp_posterior_n11n./sum(temp_posterior_n11n,2))*target_length;
            a = a + 1;
        end
        a1=sum(temp_posterior,2);
        a2=sum(temp_posterior_n11n,2);
        posterior_responseLength(temptempBoolIndex1,:) = temp_posterior_n11n;
    end
    a3_raw = sum(posterior_responseLength_raw,2);
    a3 = sum(posterior_responseLength,2);
    a = 1;
end

%% Get seq-level mean posterior from trial-level posterior, and alse get decoding accuracy for error trials
a = 1;

errorTrialNum;

valid_length;

boolIndex_location_seq_T = boolIndex_location_seq';

posterior_stimuliLength_seqMean = zeros(sum(numSeq(1:valid_length)),numFrames);
posterior_responseLength_seqMean = zeros(sum(numSeq(1:valid_length)),numFrames);

svm_posterior_stimuliErrorTrial = zeros(1,sum(numSeq(1:valid_length)));
svm_posterior_responseErrorTrial = zeros(1,sum(numSeq(1:valid_length)));

seqCount_stimuliError = zeros(1,sum(numSeq(1:valid_length)));
seqCount_responseError = zeros(1,sum(numSeq(1:valid_length)));

for tempi=1:sum(numSeq(1:valid_length))
    temp_boolIndex_location_seq_T = boolIndex_location_seq_T(tempi,:);
    
    temptempBoolIndex1 = stimuliLabel_errorTrial==temp_boolIndex_location_seq_T;
    temptempBoolIndex2 = responseLabel_errorTrial==temp_boolIndex_location_seq_T;
    
    temptempIndex1_B = find(sum(temptempBoolIndex1,2)==numFrames);
    temptempIndex2_B = find(sum(temptempBoolIndex2,2)==numFrames);
    
    seqCount_stimuliError(tempi) = length(temptempIndex1_B);
    seqCount_responseError(tempi) = length(temptempIndex2_B);
    
    if length(temptempIndex1_B) >= errorTrial_minNum
        posterior_stimuliLength_seqMean(tempi,:) = mean(posterior_stimuliLength(temptempIndex1_B,:),1);
        
        %temptempBoolIndex1_C = (posterior_stimuliLength(temptempIndex1_B,:)>0.5)==temp_boolIndex_location_seq_T;
        
        temptemp_posterior = posterior_stimuliLength(temptempIndex1_B,:);
        
        if if_decodingAcc_threshold0_sort1 == 0
            predictLabel_boolIndex = temptemp_posterior>0.5;
        elseif if_decodingAcc_threshold0_sort1 == 1
            a = 1;
            [M,I] = sort(temptemp_posterior,2,'descend');
            predictLabel_boolIndex = false(size(temptemp_posterior));
            for temptempi=1:size(predictLabel_boolIndex,1)
                for temptempj=1:sum(temp_boolIndex_location_seq_T)
                    predictLabel_boolIndex(temptempi,I(temptempi,temptempj)) = true;
                end
            end
        end
        
        temptempBoolIndex1_C = (predictLabel_boolIndex)==temp_boolIndex_location_seq_T;
        
        
        svm_posterior_stimuliErrorTrial(tempi) = sum(sum(temptempBoolIndex1_C,2)==numFrames)/length(temptempIndex1_B);
    else
        posterior_stimuliLength_seqMean(tempi,:) = nan;
        svm_posterior_stimuliErrorTrial(tempi) = nan;
    end
    
    if length(temptempIndex2_B) >= errorTrial_minNum
        posterior_responseLength_seqMean(tempi,:) = mean(posterior_responseLength(temptempIndex2_B,:),1);
        
        %temptempBoolIndex2_C = (posterior_responseLength(temptempIndex2_B,:)>0.5)==temp_boolIndex_location_seq_T;
        
        temptemp_posterior = posterior_responseLength(temptempIndex2_B,:);
        
        if if_decodingAcc_threshold0_sort1 == 0
            predictLabel_boolIndex = temptemp_posterior>0.5;
        elseif if_decodingAcc_threshold0_sort1 == 1
            a = 1;
            [M,I] = sort(temptemp_posterior,2,'descend');
            predictLabel_boolIndex = false(size(temptemp_posterior));
            for temptempi=1:size(predictLabel_boolIndex,1)
                for temptempj=1:sum(temp_boolIndex_location_seq_T)
                    predictLabel_boolIndex(temptempi,I(temptempi,temptempj)) = true;
                end
            end
        end
        
        temptempBoolIndex2_C = (predictLabel_boolIndex)==temp_boolIndex_location_seq_T;
        
        
        svm_posterior_responseErrorTrial(tempi) = sum(sum(temptempBoolIndex2_C,2)==numFrames)/length(temptempIndex2_B);
    else
        posterior_responseLength_seqMean(tempi,:) = nan;
        svm_posterior_responseErrorTrial(tempi) = nan;
    end
    
    %if tempi==29
    %    a = 1;
    %end
    
end
posterior_stimuliLength_seqMean;
posterior_responseLength_seqMean;
a = 1;

temp_boolIndex_location_seq = boolIndex_location_seq_T(1:sum(numSeq(1:valid_length)),:);

temp_p = posterior_stimuliLength_seqMean;
temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
temp_p_seq = prod(temp_p,2);
p_seq_prod_stimuli = temp_p_seq';

temp_p = posterior_responseLength_seqMean;
temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
temp_p_seq = prod(temp_p,2);
p_seq_prod_responsei = temp_p_seq';


%% Get test_errorTrial_allLength
test_errorTrial_allLength = struct;
test_errorTrial_allLength.valid_length = valid_length;
test_errorTrial_allLength.trialIndex_bool_memoryError = trialIndex_bool_memoryError;
test_errorTrial_allLength.boolIndex_location_seq_T = boolIndex_location_seq_T;
test_errorTrial_allLength.posterior_stimuliLength = posterior_stimuliLength;
test_errorTrial_allLength.posterior_responseLength = posterior_responseLength;
test_errorTrial_allLength.stimuliLabel_errorTrial = stimuliLabel_errorTrial;
test_errorTrial_allLength.responseLabel_errorTrial = responseLabel_errorTrial;
test_errorTrial_allLength.posterior_stimuliLength_seqMean = posterior_stimuliLength_seqMean;
test_errorTrial_allLength.posterior_responseLength_seqMean = posterior_responseLength_seqMean;
test_errorTrial_allLength.p_seq_prod_stimuli = p_seq_prod_stimuli;
test_errorTrial_allLength.p_seq_prod_responsei = p_seq_prod_responsei;
test_errorTrial_allLength.svm_posterior_stimuliErrorTrial = svm_posterior_stimuliErrorTrial;
test_errorTrial_allLength.svm_posterior_responseErrorTrial = svm_posterior_responseErrorTrial;
test_errorTrial_allLength.seqCount_stimuliError = seqCount_stimuliError;
test_errorTrial_allLength.seqCount_responseError = seqCount_responseError;

test_errorTrial_allLength;

%% Get svm_posterior_lengthx
for target_length=1:3
    
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    
    if if_decoder_decodingAcc0_pProd1 == 0
        if if_stimuliError0_responseError1 == 0
            svm_posterior_lengthx = svm_posterior_stimuliErrorTrial;
        elseif if_stimuliError0_responseError1 == 1
            svm_posterior_lengthx = svm_posterior_responseErrorTrial;
        end
    elseif if_decoder_decodingAcc0_pProd1 == 1
        if if_stimuliError0_responseError1 == 0
            svm_posterior_lengthx = p_seq_prod_stimuli;
        elseif if_stimuliError0_responseError1 == 1
            svm_posterior_lengthx = p_seq_prod_responsei;
        end
    end
    
    if target_length == 1
        svm_posterior_length1 = svm_posterior_lengthx(temp_range);
    elseif target_length == 2
        svm_posterior_length2 = svm_posterior_lengthx(temp_range);
    elseif target_length == 3
        svm_posterior_length3 = svm_posterior_lengthx(temp_range);
    end
end
a = 1;

%% Others
x = [svm_posterior_length1(~isnan(svm_posterior_length1))';svm_posterior_length2(~isnan(svm_posterior_length2))';svm_posterior_length3(~isnan(svm_posterior_length3))'];
y = [gAcc_length1(~isnan(svm_posterior_length1))';gAcc_length2(~isnan(svm_posterior_length2))';gAcc_length3(~isnan(svm_posterior_length3))'];
[r_123,p_123] = corr(x,y);
fprintf('r_length123_decoder_behavior=%.4f, p_length123_decoder_behavior=%.4f.\n',r_123,p_123);


[r_1,p_1] = corr(svm_posterior_length1(~isnan(svm_posterior_length1))',gAcc_length1(~isnan(svm_posterior_length1))');
[r_2,p_2] = corr(svm_posterior_length2(~isnan(svm_posterior_length2))',gAcc_length2(~isnan(svm_posterior_length2))');
[r_3,p_3] = corr(svm_posterior_length3(~isnan(svm_posterior_length3))',gAcc_length3(~isnan(svm_posterior_length3))');
% fprintf('r_1=%.4f, p_1=%.4f, mean=%.4f.\n',r_1,p_1,mean(svm_posterior_length1(~isnan(svm_posterior_length1))')); %#ok<*UDIM>
% fprintf('r_2=%.4f, p_2=%.4f, mean=%.4f.\n',r_2,p_2,mean(svm_posterior_length2(~isnan(svm_posterior_length2))'));
% fprintf('r_3=%.4f, p_3=%.4f, mean=%.4f.\n',r_3,p_3,mean(svm_posterior_length3(~isnan(svm_posterior_length3))'));

%% Temporarily ban this block, but need run it for location distribution of error trials
if false
    if if_stimuliError0_responseError1 == 0 %#ok<*UNRCH>
        Posterior_2d_n11n_mean = test_errorTrial_allLength.posterior_stimuliLength_seqMean;
    elseif if_stimuliError0_responseError1 == 1
        Posterior_2d_n11n_mean = test_errorTrial_allLength.posterior_responseLength_seqMean;
    end
    
    r_n11n = zeros(1,sum(numSeq(1:valid_length)));
    p_n11n = zeros(1,sum(numSeq(1:valid_length)));
    for tempi=1:sum(numSeq(1:valid_length))
        [r_n11n(tempi),p_n11n(tempi)] = corr(Posterior_2d_n11n_mean(tempi,:)',Posterior_2d_model(tempi,:)');
    end
    
    % fprintf('num(p_posterior_seq_n11n<0.05)=%d/%d.\n',sum(p_n11n<0.05),sum(~isnan(p_n11n)));
    a = 1;
    
end



%% Compute cross length generalization
if if_compute_test == 1
    t_test_acrossLength = tic;
    
    if if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2_CMC3 == 0
        temp_trialIndex_bool = trialIndex_bool_memoryCorrect; % use memory correct trial
    elseif if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2_CMC3 == 1
        temp_trialIndex_bool = trial_para_ifSelectOffloading==-1; % use all memory trial
    elseif if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2_CMC3 == 2
        temp_trialIndex_bool = trialIndex_bool_memoryCorrect & (trial_para_choiceCondition_flag==0); % use memory correct trial
    end
    
    resampleIterCount = resampleIterCount; %#ok<*ASGSL>
    
    for target_length=1:3
        
        if target_length == 1
            svm_train_lengthx_outputs = svm_train_length1_outputs;
        elseif target_length == 2
            svm_train_lengthx_outputs = svm_train_length2_outputs;
        elseif target_length == 3
            svm_train_lengthx_outputs = svm_train_length3_outputs;
        end
                
        % use all length trial
        tempBoolIndex = temp_trialIndex_bool;
                
        x = F_dff_decisionBin1(:,tempBoolIndex);
        y = boolIndex_location_allTrial(:,tempBoolIndex);
        
        x_T = x';
        y_T = y';
        
        temp_trialNum = size(x,2);
        
        %Posterior_2d_kfoldMean_resample = zeros(resampleIterCount,temp_trialNum,numFrames);
        Posterior_2d_kfoldMean_resample = nan(resampleIterCount,temp_trialNum,numFrames);        
        for tempIter=1:resampleIterCount
        %for tempIter=1:4:resampleIterCount
            %parfor tempIter=1:resampleIterCount
            temptemp_svm = svm_train_lengthx_outputs.temp_svm_resample{tempIter}; %#ok<*PFBNS>
            temp_Mdl_CV_binary = temptemp_svm.temp_Mdl_CV_binary;
            
            Posterior_2d_kfold = zeros(KFold_num,temp_trialNum,numFrames);
            for temploc=1:numFrames
                for tempk=1:KFold_num
                    [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary{temploc}.Trained{tempk},x_T);% Very time-consuming!!!
                    if size(tempPosterior,2) == 1
                        tempPosterior_2 = tempPosterior(:,1);
                    else
                        tempPosterior_2 = tempPosterior(:,2);
                    end
                    Posterior_2d_kfold(tempk,:,temploc) = tempPosterior_2;
                end
            end
            Posterior_2d_kfoldMean = squeeze(mean(Posterior_2d_kfold,1));
            Posterior_2d_kfoldMean_resample(tempIter,:,:) = Posterior_2d_kfoldMean;
        end
        a = 1; %#ok<*NASGU>
        Posterior_2d_kfoldMean_resampleMean = squeeze(mean(Posterior_2d_kfoldMean_resample,1,'omitnan'));
        Posterior_2d_resampleMean = Posterior_2d_kfoldMean_resampleMean;
        
        Posterior_2d_resample = Posterior_2d_kfoldMean_resample;
        
        
        Posterior_2d_seqMean = nan(sum(numSeq(1:valid_length)),numFrames);        
        for tempi=1:sum(numSeq(1:valid_length))
            temptempBoolIndex = seqIndex(temp_trialIndex_bool) == tempi;            
            temp1 = Posterior_2d_resampleMean(temptempBoolIndex,:);
            Posterior_2d_seqMean(tempi,:) = mean(temp1,1);
        end
        
        test_acrossLength = struct;
        test_acrossLength.Posterior_2d_resampleMean = Posterior_2d_resampleMean;
        test_acrossLength.temp_trialIndex_bool = temp_trialIndex_bool;        
        test_acrossLength.Posterior_2d_resample = Posterior_2d_resample;
        test_acrossLength.Posterior_2d_seqMean = Posterior_2d_seqMean;
        
        svm_train_lengthx_outputs.test_acrossLength = test_acrossLength;
        
        if target_length == 1
            svm_train_length1_outputs = svm_train_lengthx_outputs;
        elseif target_length == 2
            svm_train_length2_outputs = svm_train_lengthx_outputs;
        elseif target_length == 3
            svm_train_length3_outputs = svm_train_lengthx_outputs;
        end
        
    end
    
    fprintf('t_test_acrossLength = %.1f secs.\n',toc(t_test_acrossLength));
    
end

Posterior_2d_seqMean_trainLength1 = svm_train_length1_outputs.test_acrossLength.Posterior_2d_seqMean;
Posterior_2d_seqMean_trainLength2 = svm_train_length2_outputs.test_acrossLength.Posterior_2d_seqMean;
Posterior_2d_seqMean_trainLength3 = svm_train_length3_outputs.test_acrossLength.Posterior_2d_seqMean;


target_length = 1;
temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
Posterior_2d_seqMean_trainLength1(temp_range,:) = Posterior_2d_n11n_mean(temp_range,:);

target_length = 2;
temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
Posterior_2d_seqMean_trainLength2(temp_range,:) = Posterior_2d_n11n_mean(temp_range,:);

target_length = 3;
temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
Posterior_2d_seqMean_trainLength3(temp_range,:) = Posterior_2d_n11n_mean(temp_range,:);

a = 1;

%% End