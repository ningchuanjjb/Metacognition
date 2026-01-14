%% Compure beta from GLM
function [glm_beta_output,glm_r2_output] = fun_glm_neuron_par_v14(dff,glm_options)


%% Initialization
order_glm_valid = glm_options.order_glm_valid;
% lasso_repeatNum = glm_options.lasso_repeatNum;
numFrames = glm_options.numFrames;
plot_lengthFlag = glm_options.plot_lengthFlag;
sequence_lengthx_onehot = glm_options.sequence_lengthx_onehot;
sequence_lengthx_onehot_order = glm_options.sequence_lengthx_onehot_order;
fun_lasso_jjb = glm_options.fun_lasso_jjb;
fun_lasso_repeat_jjb = glm_options.fun_lasso_repeat_jjb; %#ok<*NASGU>
designMatrix_allSDT = glm_options.designMatrix_allSDT;
if_seqSubspace_onlyhit0_allSDT1 = glm_options.if_seqSubspace_onlyhit0_allSDT1;
if_strictHit = glm_options.if_strictHit;
if_allSDT_withOnlyHit = glm_options.if_allSDT_withOnlyHit;
if_seqTrialCountWeighted = glm_options.if_seqTrialCountWeighted;
if_glm_12B0_24B1_6B2 = glm_options.if_glm_12B0_24B1_6B2;

a = 1;

if order_glm_valid == 0
    if if_seqSubspace_onlyhit0_allSDT1 == 0
        x = sequence_lengthx_onehot;
    elseif if_seqSubspace_onlyhit0_allSDT1 == 1
        x = designMatrix_allSDT;
    end
    temp_range = 1:size(x,2);
    glm_beta = zeros(size(dff,1),size(x,2));
else
    glm_beta = zeros(size(dff,1),numFrames*plot_lengthFlag);
    x = sequence_lengthx_onehot_order;
    temp_range = (1:numFrames) + numFrames*(order_glm_valid-1);
end

a = 1;



if if_seqSubspace_onlyhit0_allSDT1 == 1
    range_hit_correctRejection = [1:numFrames (1:numFrames)+numFrames*3];
    range_miss_falseAlarm = [(1:numFrames)+numFrames*1 (1:numFrames)+numFrames*2];
    x_hit_correctRejection = x(:,range_hit_correctRejection);
    x_miss_falseAlarm = x(:,range_miss_falseAlarm);
    
    range_hit = 1:numFrames;
    range_miss = (1:numFrames)+numFrames*1;
    range_falseAlarm = (1:numFrames)+numFrames*2;
    range_correctRejection = (1:numFrames)+numFrames*3;
    
    x_hit = x(:,range_hit);
    x_miss = x(:,range_miss);
    x_falseAlarm = x(:,range_falseAlarm);
    x_correctRejection = x(:,range_correctRejection);
end

if if_strictHit == 1
    Hit_boolIndex = sum((x_miss+x_falseAlarm),2)==0;
elseif if_strictHit == 0
    Hit_boolIndex = true(size(x_hit,1),1);  
end
x = x(Hit_boolIndex,:);
temp_sequence_lengthx_onehot = sequence_lengthx_onehot(Hit_boolIndex,:);

if if_seqSubspace_onlyhit0_allSDT1 == 1
    range_hit_correctRejection = [1:numFrames (1:numFrames)+numFrames*3];
    range_miss_falseAlarm = [(1:numFrames)+numFrames*1 (1:numFrames)+numFrames*2];
    x_hit_correctRejection = x(:,range_hit_correctRejection);
    x_miss_falseAlarm = x(:,range_miss_falseAlarm);
    
    range_hit = 1:numFrames;
    range_miss = (1:numFrames)+numFrames*1;
    range_falseAlarm = (1:numFrames)+numFrames*2;
    range_correctRejection = (1:numFrames)+numFrames*3;
    
    x_hit = x(:,range_hit);
    x_miss = x(:,range_miss);
    x_falseAlarm = x(:,range_falseAlarm);
    x_correctRejection = x(:,range_correctRejection);
end

a = 1;

if if_seqTrialCountWeighted == 1    
    trialNum = size(temp_sequence_lengthx_onehot,1);
    sequence_lengthx_onehot_decimal = zeros(trialNum,1);
    for tempi=1:trialNum
        temp_seqBool = temp_sequence_lengthx_onehot(tempi,:);
        temp_str = num2str(temp_seqBool);
        sequence_lengthx_onehot_decimal(tempi) = bin2dec(temp_str);
    end
    [temp_uniqueID,ia,ic] = unique(sequence_lengthx_onehot_decimal);
    
    seqNum = length(temp_uniqueID);
    temp_uniqueIDCount = zeros(seqNum,1);
    for tempi=1:seqNum
        temp_uniqueIDCount(tempi) = sum((temp_uniqueID(tempi) == sequence_lengthx_onehot_decimal));
    end
    % temp_uniqueIDWeight = 1./temp_uniqueIDCount;
    temp_uniqueIDWeight = (1./temp_uniqueIDCount)./seqNum;
    temp_uniqueIDWeight_trialLevel = temp_uniqueIDWeight(ic);
    
elseif if_seqTrialCountWeighted == 0
    trialNum = size(temp_sequence_lengthx_onehot,1);
    n = trialNum;
    temp_uniqueIDWeight_trialLevel = 1/n*ones(n,1);
end


a = 1;
glm_r2 = zeros(size(dff,1),4);
parfor tempi=1:size(dff,1)
    % parfor tempi=1:100
    % for tempi=1:10
    warning('off');
    B = [];
    r2 = [];
    %y0 = dff(tempi,:)';    
    y0 = dff(tempi,Hit_boolIndex)';
    y = (y0-mean(y0))/std(y0); % z-score
    %y = y0-mean(y0); % only substract mean
    %y = y0;
    
    %[B,r2] = fun_lasso_repeat_jjb(x,y,lasso_repeatNum,fun_lasso_jjb);             %#ok<*ASGLU,*PFBNS>
    %[B,r2] = fun_lasso_jjb(x,y,temp_uniqueIDWeight_trialLevel); %#ok<*ASGLU,*PFBNS>
    
    if if_seqSubspace_onlyhit0_allSDT1 == 0
        [B,r2] = fun_lasso_jjb(x,y,temp_uniqueIDWeight_trialLevel); %#ok<*PFTIN,*ASGLU,*PFBNS>
        glm_r2(tempi,:) = r2;
    elseif if_seqSubspace_onlyhit0_allSDT1 == 1
        if if_glm_12B0_24B1_6B2 == 1
            [B,r2] = fun_lasso_jjb(x,y,temp_uniqueIDWeight_trialLevel); %#ok<*PFTIN,*ASGLU,*PFBNS>
            glm_r2(tempi,:) = r2;
        elseif if_glm_12B0_24B1_6B2 == 0
            [B_hit_correctRejection,r2A] = fun_lasso_jjb(x_hit_correctRejection,y,temp_uniqueIDWeight_trialLevel); %#ok<*ASGLU,*PFBNS>
            [B_miss_falseAlarm,r2B] = fun_lasso_jjb(x_miss_falseAlarm,y,temp_uniqueIDWeight_trialLevel); %#ok<*ASGLU,*PFBNS>
            B = zeros(numFrames*4,1);
            B(range_hit_correctRejection) = B_hit_correctRejection;
            B(range_miss_falseAlarm) = B_miss_falseAlarm;
            glm_r2(tempi,:) = [r2A r2B r2B r2A];
        elseif if_glm_12B0_24B1_6B2 == 2
            if if_allSDT_withOnlyHit == 0
                [B_hit,r2A] = fun_lasso_jjb(x_hit,y,temp_uniqueIDWeight_trialLevel); %#ok<*ASGLU,*PFBNS>
                [B_miss,r2B] = fun_lasso_jjb(x_miss,y,temp_uniqueIDWeight_trialLevel); %#ok<*ASGLU,*PFBNS>
                [B_falseAlarm,r2C] = fun_lasso_jjb(x_falseAlarm,y,temp_uniqueIDWeight_trialLevel); %#ok<*ASGLU,*PFBNS>
                [B_correctRejection,r2D] = fun_lasso_jjb(x_correctRejection,y,temp_uniqueIDWeight_trialLevel); %#ok<*ASGLU,*PFBNS>
            elseif if_allSDT_withOnlyHit == 1
                [B_hit,r2A] = fun_lasso_jjb(x_hit,y,temp_uniqueIDWeight_trialLevel); %#ok<*ASGLU,*PFBNS>
                B_miss = zeros(numFrames,1);
                B_falseAlarm = zeros(numFrames,1);
                B_correctRejection = zeros(numFrames,1);
                r2B = 0;
                r2C = 0;
                r2D = 0;
            end
            
            B = zeros(numFrames*4,1);
            B(range_hit) = B_hit;
            B(range_miss) = B_miss;
            B(range_falseAlarm) = B_falseAlarm;
            B(range_correctRejection) = B_correctRejection;
            glm_r2(tempi,:) = [r2A r2B r2C r2D];
        end
    end
    
    glm_beta(tempi,:) = B;
    %glm_r2(tempi) = r2; %#ok<*PFOUS>
    warning('on');
end
glm_beta_output = glm_beta(:,temp_range);
glm_r2_output = glm_r2;



%% End