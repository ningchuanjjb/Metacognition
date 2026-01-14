%% Compure beta from GLM
function [glm_beta_output,glm_r2_output] = fun_glm_neuron_eachLength_par_v3(dff,glm_options)


%% Initialization
order_glm_valid = glm_options.order_glm_valid;
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

x = designMatrix_allSDT;

range_hit = 1:numFrames;
range_miss = (1:numFrames)+numFrames*1;
range_falseAlarm = (1:numFrames)+numFrames*2;
range_correctRejection = (1:numFrames)+numFrames*3;

x_hit = x(:,range_hit);
x_miss = x(:,range_miss);
x_falseAlarm = x(:,range_falseAlarm);
x_correctRejection = x(:,range_correctRejection);

Hit_boolIndex = sum((x_miss+x_falseAlarm),2)==0;



trial_num = size(designMatrix_allSDT,1);

designMatrix_hit = designMatrix_allSDT(:,1:numFrames);

seq_length = sum(designMatrix_hit,2);
seq_length_max = max(sum(designMatrix_hit,2));

trialBoolIndex_hit_length = false(trial_num,seq_length_max);
for tempi=1:seq_length_max
   trialBoolIndex_hit_length(:,tempi) = Hit_boolIndex==true & seq_length==tempi;    
end

designMatrix_hit; %#ok<*VUNUS>
trialBoolIndex_hit_length;


% if_seqTrialCountWeighted = 0
uniqueIDWeight_trialLevel_eachLength = cell(1,seq_length_max);
for tempi=1:seq_length_max
    temp_trialNum = sum(trialBoolIndex_hit_length(:,tempi));
    n = temp_trialNum;
    uniqueIDWeight_trialLevel_eachLength{tempi} = 1/n*ones(n,1);
end

    
    
glm_beta = zeros(size(dff,1),seq_length_max*numFrames);
glm_r2 = zeros(size(dff,1),seq_length_max);
parfor tempi=1:size(dff,1)
%     parfor tempi=1:10
    %for tempi=1:10
    warning('off');
    B = zeros(1,seq_length_max*numFrames);
    r2 = zeros(1,seq_length_max);
    
    for tempj=1:seq_length_max
        temp_range = (tempj-1)*numFrames+1:tempj*numFrames;
        
        tempTrialBoolIndex = trialBoolIndex_hit_length(:,tempj);
        
        y0 = dff(tempi,tempTrialBoolIndex)';
        y = (y0-mean(y0))/std(y0); % z-score
        
        x = designMatrix_hit(tempTrialBoolIndex,:);
        
        [B(temp_range),r2(tempj)] = ...
            fun_lasso_jjb(x,y,uniqueIDWeight_trialLevel_eachLength{tempj}); %#ok<*PFTIN,*ASGLU,*PFBNS>
        
        
        if false
            
        end
        
    end
    
    
    
    glm_r2(tempi,:) = r2;        
    glm_beta(tempi,:) = B;
    warning('on');
end
glm_beta_output = glm_beta;
glm_r2_output = glm_r2;



%% End