function [Loss,svm_posterior_lengthx,Posterior_2d_w] = fun_tempModel_v9(params,tempModel_options)

Posterior_2d_n11n = tempModel_options.Posterior_2d_n11n;
temp_svm_Y_valid_T = tempModel_options.temp_svm_Y_valid_T;
numFrames = tempModel_options.numFrames;
boolIndex_location_seq = tempModel_options.boolIndex_location_seq;
gAcc = tempModel_options.gAcc;
seq_range = tempModel_options.seq_range;
target_length = tempModel_options.target_length; %#ok<*NASGU>
if_n11n = tempModel_options.if_n11n;
if_decodingAcc_threshold0_sort1 = tempModel_options.if_decodingAcc_threshold0_sort1;

% if_n11n = 1;

w_powered = params;

temp_Posterior_2d = Posterior_2d_n11n.*w_powered;

% temp_Posterior_2d = Posterior_2d_n11n;


if if_n11n == 1
    
    % temp_Posterior_2d = temp_Posterior_2d./sum(temp_Posterior_2d,2);
    temp_Posterior_2d = (temp_Posterior_2d./sum(temp_Posterior_2d,2))*target_length;
    
    a = 0;
    while true
        tempindex = temp_Posterior_2d>1;
        if sum(tempindex,'all') == 0
            break
        end
        temp_Posterior_2d(temp_Posterior_2d>1) = 1;
        %break
        temp_Posterior_2d = (temp_Posterior_2d./sum(temp_Posterior_2d,2))*target_length; %#ok<*UNRCH>
        a = a + 1;
        %break
        if a == 2
            break
        end
    end
    a;
    
elseif if_n11n == 0
    aa = 1;
end


temp_Posterior_2d(temp_Posterior_2d>1) = 1;

Posterior_2d_w = temp_Posterior_2d;

if if_decodingAcc_threshold0_sort1 == 0
    predictLabel_boolIndex = temp_Posterior_2d > 0.5;%0.5
    
elseif if_decodingAcc_threshold0_sort1 == 1
    temp_Posterior_2d;
    [M,I] = sort(temp_Posterior_2d,2,'descend');
    predictLabel_boolIndex = false(size(temp_Posterior_2d));
    for tempi=1:size(predictLabel_boolIndex,1)
        for tempj=1:target_length
            predictLabel_boolIndex(tempi,I(tempi,tempj)) = true;
        end
    end
end

tempBoolIndex = predictLabel_boolIndex == temp_svm_Y_valid_T;

predict_boolIndex = sum(tempBoolIndex,2) == numFrames;

%% svm_posterior_seqLevel
svm_posterior_lengthx = zeros(1,length(seq_range));
for tempi=1:size(svm_posterior_lengthx,2)
    tempSeq_bool = boolIndex_location_seq(:,seq_range(tempi))';
    temp1_boolIndex = sum(temp_svm_Y_valid_T == tempSeq_bool,2) == numFrames;
    temp_svm_posterior_seqLevel = sum(predict_boolIndex(temp1_boolIndex))/sum(temp1_boolIndex);
    svm_posterior_lengthx(tempi) = temp_svm_posterior_seqLevel;
end
svm_posterior_lengthx; %#ok<*VUNUS>


%% Loss Function

% LSE(Least Square Estimation/Quadratic Loss Function)


p = gAcc(~isnan(svm_posterior_lengthx));
q = svm_posterior_lengthx(~isnan(svm_posterior_lengthx));
Loss = sum((p - q).^2, 'all');