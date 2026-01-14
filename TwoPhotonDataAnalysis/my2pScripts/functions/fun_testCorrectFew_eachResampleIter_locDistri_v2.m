function [r_n11n,p_n11n] = fun_testCorrectFew_eachResampleIter_locDistri_v2(iterCount,eachIter_options)

%% Initialization
svm_train_length1_outputs = eachIter_options.svm_train_length1_outputs;
svm_train_length2_outputs = eachIter_options.svm_train_length2_outputs;
svm_train_length3_outputs = eachIter_options.svm_train_length3_outputs;

trialIndex_bool_memoryCorrectFew = eachIter_options.trialIndex_bool_memoryCorrectFew;
numFrames = eachIter_options.numFrames;
boolIndex_location_seq = eachIter_options.boolIndex_location_seq;
numSeq = eachIter_options.numSeq;
valid_length = eachIter_options.valid_length;
Posterior_2d_model = eachIter_options.Posterior_2d_model;
minTrialCount = eachIter_options.minTrialCount;
if_n11n = eachIter_options.if_n11n;

% iterCount = 1;
iterCount; %#ok<*VUNUS>

a = 1;

correctFewTrialNum = sum(trialIndex_bool_memoryCorrectFew);

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
        
        temptempBoolIndex1_C = (posterior_correctFew(temptempIndex1_B,:)>0.5)==temp_boolIndex_location_seq_T;
        svm_posterior_correctTrialFew(tempi) = sum(sum(temptempBoolIndex1_C,2)==numFrames)/length(temptempIndex1_B);
    else
        posterior_correctFew_seqMean(tempi,:) = nan;
        svm_posterior_correctTrialFew(tempi) = nan;
    end
    
end

temp_boolIndex_location_seq = boolIndex_location_seq_T(1:sum(numSeq(1:valid_length)),:);

temp_p = posterior_correctFew_seqMean;
temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
temp_p_seq = prod(temp_p,2);
p_seq_prod_correctFew = temp_p_seq'; %#ok<*NASGU>



Posterior_2d_n11n_mean = posterior_correctFew_seqMean;

r_n11n = zeros(1,sum(numSeq(1:valid_length)));
p_n11n = zeros(1,sum(numSeq(1:valid_length)));
for tempi=1:sum(numSeq(1:valid_length))
    [r_n11n(tempi),p_n11n(tempi)] = corr(Posterior_2d_n11n_mean(tempi,:)',Posterior_2d_model(tempi,:)');
end

% fprintf('num(p_posterior_seq_n11n<0.05)=%d/%d.\n',sum(p_n11n<0.05),sum(~isnan(p_n11n)));
r_n11n;
p_n11n;
