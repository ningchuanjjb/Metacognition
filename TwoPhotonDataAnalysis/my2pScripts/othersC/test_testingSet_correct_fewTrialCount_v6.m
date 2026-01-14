%% Initialization

if_compute_test = 1;

minTrialCount = 1;

if_decoder_decodingAcc0_pProd1 = 0;


if_n11n = 0;%0



seqCount = zeros(1,sum(numSeq(1:3)));
for target_length=1:3
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    temp_seqCount = zeros(1,numSeq(target_length));
    for tempi=1:numSeq(target_length)
        tempj = temp_range(tempi);
        temp_seqCount(tempi) = sum(seqIndex_valid==tempj,'all');
    end
    seqCount(temp_range) = temp_seqCount;
end

% fewTrialCountSeqBoolIndex = seqCount<resampleTrialCount & seqCount>=minTrialCount;
fewTrialCountSeqBoolIndex = resampleExcludeSeqBoolIndex & seqCount>=minTrialCount;
% fewTrialCountSeqBoolIndex = seqCount>=minTrialCount;
temp_sum = sum(fewTrialCountSeqBoolIndex);

temp1 = 1:sum(numSeq(1:3));
fewTrialCountSeqIndex = temp1(fewTrialCountSeqBoolIndex);

a = 1;

trialIndex_bool_memoryCorrect;
trialIndex_bool_memoryCorrectFew = false(1, trial_para.trial_count);

a = 1;

seqIndex;
for tempi=1:length(fewTrialCountSeqIndex)
    
    temp_seqIndex = fewTrialCountSeqIndex(tempi);
    
    temptempBoolIndex = ismember(seqIndex,temp_seqIndex);
    temptempBoolIndex2 = temptempBoolIndex & trialIndex_bool_memoryCorrect;
    %sum(temptempBoolIndex2)
    trialIndex_bool_memoryCorrectFew(temptempBoolIndex2) = true;
end
temp_sum2 = sum(trialIndex_bool_memoryCorrectFew);


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
        
        % use all length trial
        tempBoolIndex = trialIndex_bool_memoryCorrectFew;
        
        x = F_dff_decisionBin1(:,tempBoolIndex);
        y = boolIndex_location_allTrial(:,tempBoolIndex);
        
        x_T = x';
        y_T = y';
        
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
        
        test_correctTrial_few = struct;
        test_correctTrial_few.Posterior_2d_resampleMean = Posterior_2d_resampleMean;
        test_correctTrial_few.temp_trialIndex_bool = trialIndex_bool_memoryCorrectFew;
        test_correctTrial_few.stimuliLabel = y_T;
        
        test_correctTrial_few.Posterior_2d_resample = Posterior_2d_resample;
        
        svm_train_lengthx_outputs.test_correctTrial_few = test_correctTrial_few;
        
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

temp_posterior = svm_train_length1_outputs.test_correctTrial_few.Posterior_2d_resampleMean;
temp_posterior_n11n = temp_posterior./sum(temp_posterior,2);

correctFewTrialNum = sum(trialIndex_bool_memoryCorrectFew);

% %% Further compute
%% Distinguish posterior between stimuli label and response label (useless now)
posterior_correctFew = zeros(correctFewTrialNum,numFrames);
stimuliLabel_correctTrialFew = svm_train_length1_outputs.test_correctTrial_few.stimuliLabel;
for target_length=1:3
    
    if target_length == 1
        svm_train_lengthx_outputs = svm_train_length1_outputs;
    elseif target_length == 2
        svm_train_lengthx_outputs = svm_train_length2_outputs;
    elseif target_length == 3
        svm_train_lengthx_outputs = svm_train_length3_outputs;
    end
    
    temp_posterior = svm_train_lengthx_outputs.test_correctTrial_few.Posterior_2d_resampleMean;
    temp_seqLength_correctFew = sum(stimuliLabel_correctTrialFew,2);
    
    temptempBoolIndex1 = temp_seqLength_correctFew == target_length;
    
    posterior_correctFew(temptempBoolIndex1,:) = temp_posterior(temptempBoolIndex1,:);
end

posterior_correctFew_raw = posterior_correctFew;

%% to normalize posterior_correctFew
if if_n11n == 1
    %% posterior_correctFew
    temp_seqLength_correctFew = sum(stimuliLabel_correctTrialFew,2);
    for target_length=1:3
        temptempBoolIndex1 = temp_seqLength_correctFew == target_length;
        temp_posterior = posterior_correctFew(temptempBoolIndex1,:);
        
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
        posterior_correctFew(temptempBoolIndex1,:) = temp_posterior_n11n;
    end
    a3_raw = sum(posterior_correctFew_raw,2);
    a3 = sum(posterior_correctFew,2);
    a = 1;
end

%% Get seq-level mean posterior from trial-level posterior, and alse get decoding accuracy
boolIndex_location_seq_T = boolIndex_location_seq';

posterior_correctFew_seqMean = zeros(sum(numSeq(1:valid_length)),numFrames);

svm_posterior_correctTrialFew = zeros(1,sum(numSeq(1:valid_length)));

seqCount_correctFew = zeros(1,sum(numSeq(1:valid_length)));

for tempi=1:sum(numSeq(1:valid_length))
    temp_boolIndex_location_seq_T = boolIndex_location_seq_T(tempi,:);
    
    temptempBoolIndex1 = stimuliLabel_correctTrialFew==temp_boolIndex_location_seq_T;
    
    temptempIndex1_B = find(sum(temptempBoolIndex1,2)==numFrames);
    
    seqCount_correctFew(tempi) = length(temptempIndex1_B);
    
    if length(temptempIndex1_B) >= minTrialCount
        posterior_correctFew_seqMean(tempi,:) = mean(posterior_correctFew(temptempIndex1_B,:),1);

        temptemp_posterior = posterior_correctFew(temptempIndex1_B,:);
        
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

        svm_posterior_correctTrialFew(tempi) = sum(sum(temptempBoolIndex1_C,2)==numFrames)/length(temptempIndex1_B);
        a = 1;
    else
        posterior_correctFew_seqMean(tempi,:) = nan;
        svm_posterior_correctTrialFew(tempi) = nan;
    end
    
end

temp_boolIndex_location_seq = boolIndex_location_seq_T(1:sum(numSeq(1:valid_length)),:);

temp_p = posterior_correctFew_seqMean;
temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
temp_p_seq = prod(temp_p,2);
p_seq_prod_correctFew = temp_p_seq';

%% svm_posterior_lengthx for decodingAcc
% if if_decoder_decodingAcc0_pProd1 == 0
svm_posterior_length1 = svm_train_length1_outputs.svm_posterior_lengthx;
svm_posterior_length2 = svm_train_length2_outputs.svm_posterior_lengthx;
svm_posterior_length3 = svm_train_length3_outputs.svm_posterior_lengthx;


for target_length=1:3
    if target_length == 1
        svm_posterior_lengthx = svm_posterior_length1;
    elseif target_length == 2
        svm_posterior_lengthx = svm_posterior_length2;
    elseif target_length == 3
        svm_posterior_lengthx = svm_posterior_length3;
    end
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    
    svm_posterior_lengthx_raw = svm_posterior_lengthx;
    
    temp_correctFewBoolIndex = fewTrialCountSeqBoolIndex(temp_range);
    
    if if_decoder_decodingAcc0_pProd1 == 0
        temp_posterior_correctTrialFew = svm_posterior_correctTrialFew(temp_range);
        svm_posterior_lengthx(temp_correctFewBoolIndex) = temp_posterior_correctTrialFew(temp_correctFewBoolIndex);
    elseif if_decoder_decodingAcc0_pProd1 == 1
        temp_posterior_correctTrialFew = p_seq_prod_correctFew(temp_range);
        svm_posterior_lengthx(temp_correctFewBoolIndex) = temp_posterior_correctTrialFew(temp_correctFewBoolIndex);
    end
    
    if target_length == 1
        svm_posterior_length1 = svm_posterior_lengthx;
    elseif target_length == 2
        svm_posterior_length2 = svm_posterior_lengthx;
    elseif target_length == 3
        svm_posterior_length3 = svm_posterior_lengthx;
    end
end
svm_posterior_length1_decodingAcc = svm_posterior_length1;
svm_posterior_length2_decodingAcc = svm_posterior_length2;
svm_posterior_length3_decodingAcc = svm_posterior_length3;

%% svm_posterior_lengthx for pProd
% elseif if_decoder_decodingAcc0_pProd1 == 1
svm_posterior_length1 = svm_train_length1_outputs.p_seq_prod;
svm_posterior_length2 = svm_train_length2_outputs.p_seq_prod;
svm_posterior_length3 = svm_train_length3_outputs.p_seq_prod;


for target_length=1:3
    if target_length == 1
        svm_posterior_lengthx = svm_posterior_length1;
    elseif target_length == 2
        svm_posterior_lengthx = svm_posterior_length2;
    elseif target_length == 3
        svm_posterior_lengthx = svm_posterior_length3;
    end
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    
    svm_posterior_lengthx_raw = svm_posterior_lengthx;
    
    temp_correctFewBoolIndex = fewTrialCountSeqBoolIndex(temp_range);
    
    if if_decoder_decodingAcc0_pProd1 == 0
        temp_posterior_correctTrialFew = svm_posterior_correctTrialFew(temp_range);
        svm_posterior_lengthx(temp_correctFewBoolIndex) = temp_posterior_correctTrialFew(temp_correctFewBoolIndex);
    elseif if_decoder_decodingAcc0_pProd1 == 1
        temp_posterior_correctTrialFew = p_seq_prod_correctFew(temp_range);
        svm_posterior_lengthx(temp_correctFewBoolIndex) = temp_posterior_correctTrialFew(temp_correctFewBoolIndex);
    end
    
    if target_length == 1
        svm_posterior_length1 = svm_posterior_lengthx;
    elseif target_length == 2
        svm_posterior_length2 = svm_posterior_lengthx;
    elseif target_length == 3
        svm_posterior_length3 = svm_posterior_lengthx;
    end
end
svm_posterior_length1_pProd = svm_posterior_length1;
svm_posterior_length2_pProd = svm_posterior_length2;
svm_posterior_length3_pProd = svm_posterior_length3;

% end


if if_decoder_decodingAcc0_pProd1 == 0
    svm_posterior_length1 = svm_posterior_length1_decodingAcc;
    svm_posterior_length2 = svm_posterior_length2_decodingAcc;
    svm_posterior_length3 = svm_posterior_length3_decodingAcc;    
elseif if_decoder_decodingAcc0_pProd1 == 1
    svm_posterior_length1 = svm_posterior_length1_pProd;
    svm_posterior_length2 = svm_posterior_length2_pProd;
    svm_posterior_length3 = svm_posterior_length3_pProd;    
end

Posterior_2d_n11n_mean(fewTrialCountSeqBoolIndex,:) = posterior_correctFew_seqMean(fewTrialCountSeqBoolIndex,:);


x = [svm_posterior_length1';svm_posterior_length2(~isnan(svm_posterior_length2))';svm_posterior_length3(~isnan(svm_posterior_length3))'];
y = [gAcc_length1';gAcc_length2(~isnan(svm_posterior_length2))';gAcc_length3(~isnan(svm_posterior_length3))'];
[r_123,p_123] = corr(x,y);
if if_gAcc_behavior0_model1 == 0
    fprintf('r_length123_decoder_behavior=%.4f, p_length123_decoder_behavior=%.4f.\n',r_123,p_123);
elseif if_gAcc_behavior0_model1 == 1
    fprintf('r_length123_decoder_model   =%.4f, p_length123_decoder_model   =%.4f.\n',r_123,p_123);
end

%p_seq_prod_correctFew(fewTrialCountSeqBoolIndex)
%% End
