function [r_n11n,p_n11n,Posterior_2d_n11n_mean] = fun_testError_eachResampleIter_locDistri_v4(iterCount,eachIter_options)

%% Initialization
svm_train_length1_outputs = eachIter_options.svm_train_length1_outputs;
svm_train_length2_outputs = eachIter_options.svm_train_length2_outputs;
svm_train_length3_outputs = eachIter_options.svm_train_length3_outputs;

trialIndex_bool_memoryError = eachIter_options.trialIndex_bool_memoryError;
numFrames = eachIter_options.numFrames;
boolIndex_location_seq = eachIter_options.boolIndex_location_seq;
numSeq = eachIter_options.numSeq;
valid_length = eachIter_options.valid_length;
Posterior_2d_model = eachIter_options.Posterior_2d_model;
errorTrial_minNum = eachIter_options.errorTrial_minNum;
if_stimuliError0_responseError1 = eachIter_options.if_stimuliError0_responseError1;
if_n11n = eachIter_options.if_n11n;

% iterCount = 1;
iterCount; %#ok<*VUNUS>

a = 1;

errorTrialNum = sum(trialIndex_bool_memoryError);

% %% Further compute
%% Distinguish posterior between stimuli label and response label, because different length use different decoder
% seqCount_allLength_memoryError;

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
   
    a = 1;
    
    %temp_posterior = svm_train_lengthx_outputs.test_errorTrial.Posterior_2d_resampleMean;        
    temp_posterior = squeeze(svm_train_lengthx_outputs.test_errorTrial.Posterior_2d_resample(iterCount,:,:));
    
    temp_seqLength_stimuli = sum(stimuliLabel_errorTrial,2);
    temp_seqLength_response = sum(responseLabel_errorTrial,2);
        
    temptempBoolIndex1 = temp_seqLength_stimuli == target_length;
    temptempBoolIndex2 = temp_seqLength_response == target_length;
    
    posterior_stimuliLength(temptempBoolIndex1,:) = temp_posterior(temptempBoolIndex1,:);
    posterior_responseLength(temptempBoolIndex2,:) = temp_posterior(temptempBoolIndex2,:);
       
end

a = 1;
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
        
        temptempBoolIndex1_C = (posterior_stimuliLength(temptempIndex1_B,:)>0.5)==temp_boolIndex_location_seq_T;
        svm_posterior_stimuliErrorTrial(tempi) = sum(sum(temptempBoolIndex1_C,2)==numFrames)/length(temptempIndex1_B);
    else
        posterior_stimuliLength_seqMean(tempi,:) = nan;
        svm_posterior_stimuliErrorTrial(tempi) = nan;
    end
    
    if length(temptempIndex2_B) >= errorTrial_minNum
        posterior_responseLength_seqMean(tempi,:) = mean(posterior_responseLength(temptempIndex2_B,:),1);

        temptempBoolIndex2_C = (posterior_responseLength(temptempIndex2_B,:)>0.5)==temp_boolIndex_location_seq_T;
        svm_posterior_responseErrorTrial(tempi) = sum(sum(temptempBoolIndex2_C,2)==numFrames)/length(temptempIndex2_B);
    else
        posterior_responseLength_seqMean(tempi,:) = nan;
        svm_posterior_responseErrorTrial(tempi) = nan;
    end
        
end
temp_boolIndex_location_seq = boolIndex_location_seq_T(1:sum(numSeq(1:valid_length)),:);

temp_p = posterior_stimuliLength_seqMean;
temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
temp_p_seq = prod(temp_p,2);
p_seq_prod_stimuli = temp_p_seq'; %#ok<*NASGU>

temp_p = posterior_responseLength_seqMean;
temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
temp_p_seq = prod(temp_p,2);
p_seq_prod_responsei = temp_p_seq';


if if_stimuliError0_responseError1 == 0
    Posterior_2d_n11n_mean = posterior_stimuliLength_seqMean;
elseif if_stimuliError0_responseError1 == 1
    Posterior_2d_n11n_mean = posterior_responseLength_seqMean;
end

r_n11n = zeros(1,sum(numSeq(1:valid_length)));
p_n11n = zeros(1,sum(numSeq(1:valid_length)));
for tempi=1:sum(numSeq(1:valid_length))
    [r_n11n(tempi),p_n11n(tempi)] = corr(Posterior_2d_n11n_mean(tempi,:)',Posterior_2d_model(tempi,:)');
end

% fprintf('num(p_posterior_seq_n11n<0.05)=%d/%d.\n',sum(p_n11n<0.05),sum(~isnan(p_n11n)));

r_n11n;
p_n11n;
Posterior_2d_n11n_mean;
