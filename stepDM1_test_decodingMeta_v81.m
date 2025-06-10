%% Initialization
% clear
close all


%% Params setting

if_compute_train = 1;%1

if_resample_meta = 1;
resampleIterCount_meta = 16;

if_usePureNeuronDecoding = 0;%0

if_compute_shuffle = 0;
shuffleNum = 112;%32,112

valid_length = 3;%4,3


% if_trainMeta_0baseline_1delay1 = 1;

if_memoryPrecision_accuracy0_sigma1 = if_memoryPrecision_accuracy0_sigma1; %#ok<*ASGSL>

% Only temporary command!!!
% if true
if false
    if if_trainMeta_0baseline_1delay1 == 0
        if_resample_meta = 0;
    end
end

exampleSeq = 35;%17,35
if_plot_meta_trialLevelEvidence = 1;
if_plot_additionalSmooth = 1;

if_trialEvidenceExampleSeq_histogram0_violinPlot1 = 1;
if_trialEvidenceAllSeq_violinplot0_pairline1 = 1;


if_plot = 1;

if_prctile_metaDecoderThreshold = 1;

color_choiceMemory = [1,133+20,113]/255;
color_choiceOffload = [166+40,97,26]/255;


% KFold_num_meta = 10;%10-->5-->15-->10(ROC824)
% KFold_num_meta = 20;%10(ROC799)-->15(ROC811)-->20(ROC818)
KFold_num_meta = 20;%10(ROC812)-->15(ROC816)-->20(ROC819)

rng(5);%5
t_decoder = templateSVM('Standardize',true,'KernelFunction','linear'); % To standardise input data
% t_decoder = templateSVM('Standardize',true,'KernelFunction','polynomial'); % To standardise input data


targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)
output_shortPath = 'D:\twoPhotonData_motionCorrected';

if if_dff_singleSession1_twoSession2_allSession3 == 1
    currentSession_multi = string;
    
    if if_monkey_D0_Z1 == 0
        currentSession_multi = [currentSession_multi; '113Recording_20230426A_Ding_Site16'];%1
        currentSession_multi = [currentSession_multi; '113Recording_20230427A_Ding_Site16_sameFOV0426'];%2
        currentSession_multi = [currentSession_multi; '113Recording_20230502A_Ding_Site13'];%3
        currentSession_multi = [currentSession_multi; '113Recording_20230503A_Ding_Site13_sameFOV0502'];%4
        currentSession_multi = [currentSession_multi; '113Recording_20230504A_Ding_Site02'];%5
        currentSession_multi = [currentSession_multi; '113Recording_20230508A_Ding_Site02_sameFOV0509'];%6
        currentSession_multi = [currentSession_multi; '113Recording_20230509A_Ding_Site02'];%7, 660000 frames, easy to crash
        
        currentSession_multi = [currentSession_multi; '113Recording_20230510A_Ding_Site05_sameFOV0511'];%8
        currentSession_multi = [currentSession_multi; '113Recording_20230510B_Ding_Site05_sameFOV0511'];%9
        currentSession_multi = [currentSession_multi; '113Recording_20230511A_Ding_Site05'];%10
        currentSession_multi = [currentSession_multi; '113Recording_20230512A_Ding_Site09'];%11
        currentSession_multi = [currentSession_multi; '113Recording_20230513A_Ding_Site09_sameFOV0512'];%12
        
        currentSession_multi = [currentSession_multi; '113Recording_20230515A_Ding_Site24_sameFOV0516'];%13
        currentSession_multi = [currentSession_multi; '113Recording_20230516A_Ding_Site24'];%14
        currentSession_multi = [currentSession_multi; '113Recording_20230517A_Ding_Site16B'];%15
        currentSession_multi = [currentSession_multi; '113Recording_20230522A_Ding_Site05B'];%16
        currentSession_multi = [currentSession_multi; '113Recording_20230525A_Ding_Site09B'];%17
        currentSession_multi = [currentSession_multi; '113Recording_20230526A_Ding_Site09B_sameFOV0525'];%18
        
        currentSession_multi = [currentSession_multi; '113Recording_20230527A_Ding_Site02B'];%19
        currentSession_multi = [currentSession_multi; '113Recording_20230528A_Ding_Site02B_sameFOV0527'];%20
        currentSession_multi = [currentSession_multi; '113Recording_20230530A_Ding_Site05C'];%21
        currentSession_multi = [currentSession_multi; '113Recording_20230531A_Ding_Site05C_sameFOV0530'];%22
        currentSession_multi = [currentSession_multi; '113Recording_20230601A_Ding_Site13B'];%23
        currentSession_multi = [currentSession_multi; '113Recording_20230602A_Ding_Site13B_sameFOV0601'];%24
        
        currentSession_multi = [currentSession_multi; '113Recording_20230604A_Ding_Site07'];%25
        currentSession_multi = [currentSession_multi; '113Recording_20230605A_Ding_Site07_sameFOV0604'];%26
        currentSession_multi = [currentSession_multi; '113Recording_20230612A_Ding_Site14'];%27
        currentSession_multi = [currentSession_multi; '113Recording_20230614A_Ding_Site05D'];%28
        currentSession_multi = [currentSession_multi; '113Recording_20230615A_Ding_Site05D_sameFOV0614'];%29
        currentSession_multi = [currentSession_multi; '113Recording_20230619A_Ding_Site02C'];%30
        currentSession_multi = [currentSession_multi; '113Recording_20230620A_Ding_Site05E'];%31
        
    elseif if_monkey_D0_Z1 == 1
        currentSession_multi = [currentSession_multi; '113Recording_20240111A_Zelku_Site09A'];%1
        currentSession_multi = [currentSession_multi; '113Recording_20240112A_Zelku_Site06A'];%2
        currentSession_multi = [currentSession_multi; '113Recording_20240115A_Zelku_Site06A'];%3
        currentSession_multi = [currentSession_multi; '113Recording_20240117A_Zelku_Site14A'];%4
        currentSession_multi = [currentSession_multi; '113Recording_20240118A_Zelku_Site18A'];%5
        currentSession_multi = [currentSession_multi; '113Recording_20240119A_Zelku_Site17A'];%6
        currentSession_multi = [currentSession_multi; '113Recording_20240122A_Zelku_Site09B'];%7
        currentSession_multi = [currentSession_multi; '113Recording_20240123A_Zelku_Site09B_sameFOV0122'];%8
        currentSession_multi = [currentSession_multi; '113Recording_20240124A_Zelku_Site06B'];%9
        currentSession_multi = [currentSession_multi; '113Recording_20240126A_Zelku_Site06B_sameFOV0124'];%10
        
    end
    
    currentSession_multi(1) = [];
    num_FOV = length(currentSession_multi);
    
    currentSession = currentSession_multi{currentSessionIndex};
    
    fprintf('currentSession = %s.\n',currentSession);
    
    temp_currentSession_path = [output_shortPath '\' currentSession];
    temp_if_max0_min1 = 0;
    output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
    
    currentSession_short = currentSession(18:22);
    
elseif if_dff_singleSession1_twoSession2_allSession3 == 2
    
    currentABSession_multi = string;
    
    if if_monkey_D0_Z1 == 0
        currentABSession_multi = [currentABSession_multi; '20230426A_and_20230427A'];%1 few
        currentABSession_multi = [currentABSession_multi; '20230502A_and_20230503A'];%2
        currentABSession_multi = [currentABSession_multi; '20230504A_and_20230508A_and_20230509A'];%3
        
        currentABSession_multi = [currentABSession_multi; '20230510A_and_20230510B_and_20230511A'];%4
        currentABSession_multi = [currentABSession_multi; '20230512A_and_20230513A'];%5 few
        
        currentABSession_multi = [currentABSession_multi; '20230515A_and_20230516A'];%6
        currentABSession_multi = [currentABSession_multi; '20230525A_and_20230526A'];%7
        currentABSession_multi = [currentABSession_multi; '20230527A_and_20230528A'];%8
        currentABSession_multi = [currentABSession_multi; '20230530A_and_20230531A'];%9
        currentABSession_multi = [currentABSession_multi; '20230601A_and_20230602A'];%10
        
        currentABSession_multi = [currentABSession_multi; '20230604A_and_20230605A'];%11
        currentABSession_multi = [currentABSession_multi; '20230614A_and_20230615A'];%12
        
        
    elseif if_monkey_D0_Z1 == 1
        
        currentABSession_multi = [currentABSession_multi; '20240122A_and_20240123A'];%1
        currentABSession_multi = [currentABSession_multi; '20240124A_and_20240126A'];%2
        currentABSession_multi = [currentABSession_multi; '20240129A_and_20240131A'];%3
        currentABSession_multi = [currentABSession_multi; '20240202A_and_20240203A'];%4
        currentABSession_multi = [currentABSession_multi; '20240207A_and_20240208A'];%5
        currentABSession_multi = [currentABSession_multi; '20240210A_and_20240211A'];%6
        currentABSession_multi = [currentABSession_multi; '20240216A_and_20240218A'];%7
        currentABSession_multi = [currentABSession_multi; '20240220A_and_20240221A'];%8
        currentABSession_multi = [currentABSession_multi; '20240226A_and_20240227A'];%9
        currentABSession_multi = [currentABSession_multi; '20240229A_and_20240301A'];%10
        currentABSession_multi = [currentABSession_multi; '20240304A_and_20240305A'];%11
        currentABSession_multi = [currentABSession_multi; '20240307A_and_20240308A'];%12
        currentABSession_multi = [currentABSession_multi; '20240312A_and_20240315A'];%13
        currentABSession_multi = [currentABSession_multi; '20240319A_and_20240320A'];%14
        currentABSession_multi = [currentABSession_multi; '20240322A_and_20240323A'];%15
        currentABSession_multi = [currentABSession_multi; '20240329A_and_20240330A'];%16
        currentABSession_multi = [currentABSession_multi; '20240402A_and_20240403A'];%17
        currentABSession_multi = [currentABSession_multi; '20240410A_and_20240411A'];%18
        
    end
    
    currentABSession_multi(1) = [];
    num_FOV_AB = length(currentABSession_multi);
    
    currentSession = currentABSession_multi{currentSessionIndex_AB};
    currentSession_short = [currentSession(5:9),'+',currentSession(19:23)];
    
    fprintf('currentSession = %s.\n',currentSession);
    
    output_path = 'D:\twoPhotonData_motionCorrected\twoSessions';
elseif if_dff_singleSession1_twoSession2_allSession3 == 3
    
end

%% Load data
t0 = tic;
numFrames = 6;
pointKindsNum = 4;
target_seqSet = get_target_seqSet_v1(numFrames,pointKindsNum);
numSeq = zeros(1,pointKindsNum);
for tempi=1:pointKindsNum
    numSeq(tempi) = length(target_seqSet{tempi});
end

% if false
if true
    if exist('decodingDataSimplified','var') == 0 || if_compute_train == 1 %#ok<*UNRCH>
        if if_dff_singleSession1_twoSession2_allSession3 == 1
            temp_load = load([output_path,'\','decodingDataSimplified.mat'],'decodingDataSimplified');
            decodingDataSimplified = temp_load.decodingDataSimplified;
        elseif if_dff_singleSession1_twoSession2_allSession3 == 2
            temp_load = load([output_path '\' currentSession]);
            decodingDataSimplified = temp_load.decodingDataSimplified_AB;
        elseif if_dff_singleSession1_twoSession2_allSession3 == 3
            
        end
    end
end

temp_fields = fieldnames(decodingDataSimplified);
for tempi=1:length(temp_fields)
    eval([temp_fields{tempi},'=decodingDataSimplified.(temp_fields{tempi});']);
end

if if_monkey_D0_Z1 == 0
    searchName_rProb = 'from23-04-26to23-06-20_D_offloadingProb_1';
elseif if_monkey_D0_Z1 == 1
    %searchName_rProb = 'from01-08to01-23_Z_offloadingProb_1';
    searchName_rProb = 'from01-08to03-30_Z_offloadingProb_1';
end

path_behavior = [output_shortPath '\behavior'];

% Load other processed results
load_rProb = loadMat_single(searchName_rProb, path_behavior);
offloadingProb_collapsed = load_rProb.offloadingProb_all;

offloadingProb_inOne = [];
for tempi=1:pointKindsNum
    offloadingProb_inOne = [offloadingProb_inOne offloadingProb_collapsed{tempi}']; %#ok<*AGROW>
end

cd(targetPATH)

% %% Load trial data
a = 1; %#ok<*NASGU>

trialIndex_bool_memoryCorrect;

trial_para_choiceCondition_flag;

temp_trial_para_ifSelectOffloading = trial_para_ifSelectOffloading==1;

choiceBoolIndex = trial_para_choiceCondition_flag == 2;
allMemoryBoolIndex = ~temp_trial_para_ifSelectOffloading;
choiceMemoryBoolIndex = choiceBoolIndex & allMemoryBoolIndex;
choiceOffloadBoolIndex = choiceBoolIndex & temp_trial_para_ifSelectOffloading;
noChoiceBoolIndex = ~choiceBoolIndex;

allMemoryCorrectBoolIndex = trialIndex_bool_memoryCorrect;
allMemoryErrorBoolIndex = allMemoryBoolIndex & ~trialIndex_bool_memoryCorrect;

choiceMemoryCorrectBoolIndex = choiceMemoryBoolIndex & trialIndex_bool_memoryCorrect;
choiceMemoryErrorBoolIndex = choiceMemoryBoolIndex & ~trialIndex_bool_memoryCorrect;

% temp_seqIndex = seqIndex;
if if_memeoryPrecision_stimuli0_response1 == 0
    temp_seqIndex = seqIndex;
elseif if_memeoryPrecision_stimuli0_response1 == 1
    temp_seqIndex = seqIndex_response;
end


seqIndex_choice = temp_seqIndex(choiceBoolIndex);

temp_range = 1:sum(numSeq(1:valid_length));
temptempBoolIndex_validLength = ismember(temp_seqIndex,temp_range);
seqIndex_choice_validLength = temp_seqIndex(choiceBoolIndex&temptempBoolIndex_validLength);


choiceBoolIndex;
choiceMemoryBoolIndex;
choiceOffloadBoolIndex;
noChoiceBoolIndex;

choiceBoolIndex_validLength = choiceBoolIndex&temptempBoolIndex_validLength;
choiceMemoryBoolIndex_validLength = choiceMemoryBoolIndex&temptempBoolIndex_validLength;
choiceOffloadBoolIndex_validLength = choiceOffloadBoolIndex&temptempBoolIndex_validLength;
noChoiceBoolIndex_validLength = noChoiceBoolIndex&temptempBoolIndex_validLength;


%% if_usePureNeuronDecoding
neuronBoolIndex = true(size(F_dff_decisionPeriodA,1),1); %#ok<*PREALL>

if if_usePureNeuronDecoding == 1
    r2_rProb;
    r2_6loc;
    p_rProb;
    p_6loc;
    
    % two variables
    if false
        temptempBoolIndex1 = (p_rProb < 0.01) & (p_6loc < 0.01);
        temptempBoolIndex2 = ~temptempBoolIndex1;
    end
    
    % three variables
    if true
        temptempBoolIndex1 = (p_rProb < 0.01) & (p_6loc < 0.01) & (p_choiceMemory_baseline < 0.01);
        
        %temptempBoolIndex2 = (p_rProb < 0.01) | (p_6loc < 0.01) | (p_choiceMemory_baseline < 0.01);
        %temptemp2_A = (p_rProb < 0.01) & (p_6loc >= 0.01) & (p_choiceMemory_baseline >= 0.01);
        %temptemp2_B = (p_rProb >= 0.01) & (p_6loc < 0.01) & (p_choiceMemory_baseline >= 0.01);
        %temptemp2_C = (p_rProb >= 0.01) & (p_6loc >= 0.01) & (p_choiceMemory_baseline < 0.01);
        %temptempBoolIndex2 = temptemp2_A | temptemp2_B | temptemp2_C;
        
        temptemp2 = (p_rProb < 0.01) + (p_6loc < 0.01) + (p_choiceMemory_baseline < 0.01);
        temptempBoolIndex2 = temptemp2==1;
    end

    
    temptempBoolIndex_mixed = temptempBoolIndex1;
    temptempBoolIndex_pure = temptempBoolIndex2;
    
    neuronBoolIndex = temptempBoolIndex_pure;
%     neuronBoolIndex = temptempBoolIndex_mixed;% only for test!!!
    
    neuronBoolIndex2 = neuronBoolIndex;
    %neuronBoolIndex = neuronBoolIndex2;
end

%% More preparation
if if_dff_singleSession1_twoSession2_allSession3 < 3
    if if_trainMeta_0baseline_1delay1 == 1
        F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
        
    elseif if_trainMeta_0baseline_1delay1 == 0
        temptemp_range = (baselinePeriod_interval(2)-8):baselinePeriod_interval(3);
        F_dff_decisionBin1 = mean(F_dff_baselinePeriod(:,:,temptemp_range),3);
    end
    
    F_dff_decisionBin1 = double(F_dff_decisionBin1);
    F_dff_decisionBin1 = F_dff_decisionBin1 + eps;
    
    F_dff_decisionBin1 = F_dff_decisionBin1(neuronBoolIndex,:);
    
    roiNum = size(F_dff_decisionBin1,1);
    
elseif if_dff_singleSession1_twoSession2_allSession3 == 3
    
end
a = 1;


%% Get resampleTrialCount_meta
validChoiceMemoryTrialNum = sum(choiceMemoryBoolIndex_validLength);
validChoiceOffloadTrialNum = sum(choiceOffloadBoolIndex_validLength);

resampleTrialCount_meta = min([validChoiceMemoryTrialNum,validChoiceOffloadTrialNum]);

%% Get XY
% x_svm
x_raw = F_dff_decisionBin1';
x_svm = x_raw(choiceBoolIndex_validLength,:);

% y_svm
y_raw = nan((length(choiceMemoryBoolIndex)),1);
y_raw(choiceMemoryBoolIndex_validLength) = true;
y_raw(choiceOffloadBoolIndex_validLength) = false;

y_svm = y_raw(~isnan(y_raw));
y_svm = y_svm==true;


%% Get resampled XY
resampleIterCount_meta;
resampleTrialCount_meta;

validTrialIndex = find((choiceBoolIndex_validLength)==true);
validChoiceMemoryTrialIndex = find((choiceMemoryBoolIndex_validLength)==true);
validChoiceOffloadTrialIndex = find((choiceOffloadBoolIndex_validLength)==true);

validChoiceMemoryTrialNum;
validChoiceOffloadTrialNum;
resampleTrialCount_meta;

validChoiceMemoryTrialIndex_resampled = nan(resampleIterCount_meta,resampleTrialCount_meta);
for tempi=1:resampleIterCount_meta
    a2 = sort(randperm(validChoiceMemoryTrialNum,resampleTrialCount_meta));
    a2_resampled = validChoiceMemoryTrialIndex(a2);
    validChoiceMemoryTrialIndex_resampled(tempi,:) = a2_resampled;
end

validChoiceOffloadTrialIndex_resampled = nan(resampleIterCount_meta,resampleTrialCount_meta);
for tempi=1:resampleIterCount_meta
    a2 = sort(randperm(validChoiceOffloadTrialNum,resampleTrialCount_meta));
    a2_resampled = validChoiceOffloadTrialIndex(a2);
    validChoiceOffloadTrialIndex_resampled(tempi,:) = a2_resampled;
end

validChoiceMemoryTrialIndex_resampled;
validChoiceOffloadTrialIndex_resampled;

y_svm_resampled = false(resampleTrialCount_meta*2,resampleIterCount_meta);
y_svm_resampled(1:resampleTrialCount_meta,:) = true;

x_svm_resampled = nan(resampleTrialCount_meta*2,roiNum,resampleIterCount_meta);
for tempi=1:resampleIterCount_meta
    a1 = x_raw(validChoiceMemoryTrialIndex_resampled(tempi,:),:);
    a2 = x_raw(validChoiceOffloadTrialIndex_resampled(tempi,:),:);
    a12 = [a1;a2];
    x_svm_resampled(:,:,tempi) = a12;
end

x_svm_resampled;
y_svm_resampled;

%% resampled seqIndex
seqIndex_choice_validLength;
temp_seqIndex;

validChoiceMemoryTrialIndex_resampled;
validChoiceOffloadTrialIndex_resampled;

temp_seqIndex_resampled = nan(resampleIterCount_meta,resampleTrialCount_meta*2);
for tempi=1:resampleIterCount_meta
    a1 = temp_seqIndex(validChoiceMemoryTrialIndex_resampled(tempi,:));
    a2 = temp_seqIndex(validChoiceOffloadTrialIndex_resampled(tempi,:));
    a12 = [a1,a2];
    temp_seqIndex_resampled(tempi,:) = a12;
end
seqIndex_choice_validLength_resampled = temp_seqIndex_resampled;

%% Compute svm
if if_compute_train == 1
    fun_metaDecoder_Name_v = autoGetFunName_myScripts('metaDecoder', [targetPATH '\functions']);
    fun_metaDecoder = str2func(fun_metaDecoder_Name_v);
    
    if if_resample_meta == 0
        svm_options = struct;
        svm_options.numSeq = numSeq;
        svm_options.seqIndex_choice = seqIndex_choice_validLength;
        svm_options.KFold_num = KFold_num_meta;
        svm_options.t_decoder = t_decoder;
        svm_options.if_resample = if_resample_meta;
        
        svm_Meta = fun_metaDecoder(x_svm,y_svm,svm_options);
        
        fprintf('Meta-WM decoding accuracy (all trial) = %.3f.\n',svm_Meta.svm_posterior_allTrial);
        
    elseif if_resample_meta == 1
        
        x_svm_resampled;
        y_svm_resampled;
        
        temp_svm_Meta_resampled = cell(resampleIterCount_meta,1);
        parfor tempi=1:resampleIterCount_meta
            %for tempi=1:resampleIterCount_meta
            rng(5);%5
            
            svm_options = struct;
            svm_options.numSeq = numSeq;
            svm_options.seqIndex_choice = seqIndex_choice_validLength_resampled(tempi,:);
            svm_options.KFold_num = KFold_num_meta;
            svm_options.t_decoder = t_decoder;
            svm_options.if_resample = if_resample_meta;
            
            temp_x = x_svm_resampled(:,:,tempi);
            temp_y = y_svm_resampled(:,tempi);
            temp_svm_Meta = fun_metaDecoder(temp_x,temp_y,svm_options); %#ok<*PFBNS>
            
            temp_svm_Meta_resampled{tempi} = temp_svm_Meta;
        end
        
        Posterior_2d = nan(length(choiceBoolIndex),resampleIterCount_meta);
        
        a1 = validChoiceMemoryTrialIndex_resampled;
        a2 = validChoiceOffloadTrialIndex_resampled;
        a12 = [a1,a2];
        
        for tempi=1:resampleIterCount_meta
            temp_posterior = temp_svm_Meta_resampled{tempi}.Posterior_2d_n11n;
            
            temp_trialIndex = a12(tempi,:);
            
            Posterior_2d(temp_trialIndex,tempi) = temp_posterior;
        end
        
        Posterior_2d_mean = nan(length(choiceBoolIndex),1);
        for tempi=1:size(Posterior_2d_mean,1)
            temp1 = Posterior_2d(tempi,:);
            Posterior_2d_mean(tempi) = mean(temp1,'omitnan');
        end
        
        svm_Meta_resampled = struct;
        svm_Meta_resampled.validChoiceMemoryTrialIndex_resampled = validChoiceMemoryTrialIndex_resampled;
        svm_Meta_resampled.validChoiceOffloadTrialIndex_resampled = validChoiceOffloadTrialIndex_resampled;
        svm_Meta_resampled.validChoiceMemoryTrialNum = validChoiceMemoryTrialNum;
        svm_Meta_resampled.validChoiceOffloadTrialNum = validChoiceOffloadTrialNum;
        svm_Meta_resampled.resampleTrialCount = resampleTrialCount_meta;
        
        svm_Meta_resampled.Posterior_2d = Posterior_2d;
        svm_Meta_resampled.Posterior_2d_mean = Posterior_2d_mean;
        
        svm_Meta_resampled.Posterior_2d_mean_valid = Posterior_2d_mean(~isnan(y_raw));
        svm_Meta_resampled.y_svm = y_raw(~isnan(y_raw));
        
        svm_Meta_resampled.temp_svm_Meta_resampled = temp_svm_Meta_resampled;
        
        svm_Meta = svm_Meta_resampled;
    end
end

%%
if if_prctile_metaDecoderThreshold == 1
    temp_thresholdRange = 0.01:0.001:0.99;
    
    hit_minus_falseAlarm_multi = zeros(1,length(temp_thresholdRange));
    svm_choiceMemory_hit_multi = zeros(1,length(temp_thresholdRange));
    svm_choiceMemory_falseAlarm_multi = zeros(1,length(temp_thresholdRange));
        
    for temptempi=1:length(temp_thresholdRange)
        
        temp_threshold = temp_thresholdRange(temptempi);
        
        if if_resample_meta == 0
            predictLabel = svm_Meta.Posterior_2d_n11n > temp_threshold;%0.5
            trueLabel = svm_Meta.temp_svm_Y_valid_T;
        elseif if_resample_meta == 1
            predictLabel = svm_Meta.Posterior_2d_mean_valid > temp_threshold;%0.5
            trueLabel = svm_Meta.y_svm;
        end
        
        svm_choiceMemory_hit = sum(predictLabel(trueLabel==true))/sum(trueLabel==true);
        svm_choiceMemory_falseAlarm = sum(predictLabel(trueLabel==false))/sum(trueLabel==false);  
        
        hit_minus_falseAlarm = svm_choiceMemory_hit - svm_choiceMemory_falseAlarm;
        
        svm_choiceMemory_hit_multi(temptempi) = svm_choiceMemory_hit;
        svm_choiceMemory_falseAlarm_multi(temptempi) = svm_choiceMemory_falseAlarm;   
        
        hit_minus_falseAlarm_multi(temptempi) = hit_minus_falseAlarm;
    end
    
    [M,I] = max(hit_minus_falseAlarm_multi);
    
    metaDecoderThreshold = temp_thresholdRange(I);
    
    svm_choiceMemory_hit = svm_choiceMemory_hit_multi(I);
    svm_choiceMemory_falseAlarm = svm_choiceMemory_falseAlarm_multi(I);
    
elseif if_prctile_metaDecoderThreshold == 0
    metaDecoderThreshold = 0.5;
end

%% shuffle
% if if_compute_train == 1
if true
    if if_compute_shuffle == 1
        
        %% Preparation
        validChoiceMemoryTrialIndex_shuffled = nan(shuffleNum,resampleTrialCount_meta);
        for tempi=1:shuffleNum
            a2 = sort(randperm(validChoiceMemoryTrialNum,resampleTrialCount_meta));
            a2_resampled = validChoiceMemoryTrialIndex(a2);
            validChoiceMemoryTrialIndex_shuffled(tempi,:) = a2_resampled;
        end
        
        validChoiceOffloadTrialIndex_shuffled = nan(shuffleNum,resampleTrialCount_meta);
        for tempi=1:shuffleNum
            a2 = sort(randperm(validChoiceOffloadTrialNum,resampleTrialCount_meta));
            a2_resampled = validChoiceOffloadTrialIndex(a2);
            validChoiceOffloadTrialIndex_shuffled(tempi,:) = a2_resampled;
        end
        
        x_raw = F_dff_decisionBin1';
        
        y_svm_shuffled = false(resampleTrialCount_meta*2,shuffleNum);
        y_svm_shuffled(1:resampleTrialCount_meta,:) = true;
        for tempi=1:shuffleNum
            temp1 = y_svm_shuffled(:,tempi);
            temp1_shuffled = temp1(randperm(length(temp1)));
            y_svm_shuffled(:,tempi) = temp1_shuffled;
        end
        
        
        x_svm_shuffled = nan(resampleTrialCount_meta*2,roiNum,shuffleNum);
        for tempi=1:shuffleNum
            a1 = x_raw(validChoiceMemoryTrialIndex_shuffled(tempi,:),:);
            a2 = x_raw(validChoiceOffloadTrialIndex_shuffled(tempi,:),:);
            a12 = [a1;a2];
            x_svm_shuffled(:,:,tempi) = a12;
        end
        
        x_svm_shuffled;
        y_svm_shuffled;


        temp_seqIndex_shuffled = nan(shuffleNum,resampleTrialCount_meta*2);
        for tempi=1:shuffleNum
            a1 = temp_seqIndex(validChoiceMemoryTrialIndex_shuffled(tempi,:));
            a2 = temp_seqIndex(validChoiceOffloadTrialIndex_shuffled(tempi,:));
            a12 = [a1,a2];
            temp_seqIndex_shuffled(tempi,:) = a12;
        end
        seqIndex_choice_validLength_shuffled = temp_seqIndex_shuffled;

        
        %% Compute
        AUROC_shuffled = zeros(1,shuffleNum);
        
        %for shuffleCount=1:shuffleNum
        parfor shuffleCount=1:shuffleNum
            warning('off');
            predictLabel_shuffled_length1 = [];
            predictLabel_shuffled_length2 = [];
            predictLabel_shuffled_length3 = [];
            trueLabel_shuffled_length1 = [];
            trueLabel_shuffled_length2 = [];
            trueLabel_shuffled_length3 = [];
            
            %% Get choiceMemoryBoolIndex_shuffled
            choiceBoolIndex;
            choiceMemoryBoolIndex;
            
            temptempIndex_choice = find(choiceBoolIndex_validLength==true);
            temptempIndex_choiceMemory = find(choiceMemoryBoolIndex==true);
            % temptempIndex_choiceOffload = find(choiceOffloadBoolIndex==true);
            % a1 = ismember(temptempIndex_choiceMemory,temptempIndex_choice);
            
            tempChoiceNum = sum(choiceBoolIndex_validLength);
            tempChoiceMemoryNum = sum(choiceMemoryBoolIndex);
            
            tempShuffle = sort(randperm(tempChoiceNum,tempChoiceMemoryNum));
            
            choiceMemoryBoolIndex_shuffled = false(1,length(choiceMemoryBoolIndex));
            
            % choiceMemoryBoolIndex_shuffled(temptempIndex_choiceMemory) = true;
            temptempIndex_choiceMemory_shuffled = temptempIndex_choice(tempShuffle);
            choiceMemoryBoolIndex_shuffled(temptempIndex_choiceMemory_shuffled) = true;
            
            a2 = sum(choiceMemoryBoolIndex_shuffled==choiceMemoryBoolIndex);
            
            
            choiceMemoryBoolIndex_shuffled;
            
            svm_options = struct;
            svm_options.numSeq = numSeq;
            svm_options.seqIndex_choice = seqIndex_choice_validLength_shuffled(shuffleCount,:);
            svm_options.KFold_num = KFold_num_meta;
            svm_options.t_decoder = t_decoder;
            svm_options.if_resample = if_resample_meta;
            
            temp_x = x_svm_shuffled(:,:,shuffleCount);
            temp_y = y_svm_shuffled(:,shuffleCount);
            svm_Meta_shuffled = fun_metaDecoder(temp_x,temp_y,svm_options); %#ok<*PFBNS>
            

            %% Get AUROC_shuffled    
            
            temp_thresholdRange_shuffled = 0.01:0.001:0.99;
            
            hit_minus_falseAlarm_multi_shuffled = zeros(1,length(temp_thresholdRange_shuffled));
            svm_choiceMemory_hit_multi_shuffled = zeros(1,length(temp_thresholdRange_shuffled));
            svm_choiceMemory_falseAlarm_multi_shuffled = zeros(1,length(temp_thresholdRange_shuffled));
            
            predictLabel = [];
            trueLabel = [];            
            
            for temptempi=1:length(temp_thresholdRange_shuffled)
                
                temp_threshold = temp_thresholdRange_shuffled(temptempi);
                
                predictLabel = svm_Meta_shuffled.Posterior_2d_n11n > temp_threshold;%0.5
                trueLabel = svm_Meta_shuffled.temp_svm_Y_valid_T;
                
                
                svm_choiceMemory_hit_shuffled = sum(predictLabel(trueLabel==true))/sum(trueLabel==true);
                svm_choiceMemory_falseAlarm_shuffled = sum(predictLabel(trueLabel==false))/sum(trueLabel==false);
                
                
                svm_choiceMemory_hit_multi_shuffled(temptempi) = svm_choiceMemory_hit_shuffled;
                svm_choiceMemory_falseAlarm_multi_shuffled(temptempi) = svm_choiceMemory_falseAlarm_shuffled;
                
            end
            
            x = svm_choiceMemory_falseAlarm_multi_shuffled; %#ok<*PFTUSW>
            y = svm_choiceMemory_hit_multi_shuffled;
            
            x = x(end:-1:1);
            y = y(end:-1:1);
            
            AUROC_shuffled(shuffleCount) = trapz(x,y);
            
            
            warning('on');
        end
        
        %AUROC;
        AUROC_shuffled;
        
        AUROC_shuffled_prctile = prctile(AUROC_shuffled,95);
        
    elseif if_compute_shuffle == 0
        
        if if_trainMeta_0baseline_1delay1 == 1
            if if_monkey_D0_Z1 == 0
                AUROC_shuffled_prctile = 0.530;
            elseif if_monkey_D0_Z1 == 1
                AUROC_shuffled_prctile = 0.530;% FOV9
            end
        elseif if_trainMeta_0baseline_1delay1 == 0
            if if_monkey_D0_Z1 == 0
                AUROC_shuffled_prctile = 0.534;
            elseif if_monkey_D0_Z1 == 1
                AUROC_shuffled_prctile = 0.530;% FOV9
                %AUROC_shuffled_prctile = 0.539;% FOV13
            end            
        end
        
        
    end
end

if if_trainMeta_0baseline_1delay1 == 0
    metaDecoderThreshold_baseline = metaDecoderThreshold;
elseif if_trainMeta_0baseline_1delay1 == 1
    metaDecoderThreshold_delay1 = metaDecoderThreshold;
end


fprintf('t_script = %.1f secs.\n',toc(t0));


%%  Compute meta-memory on training set trials (choice trials)
if if_resample_meta == 0
    meta_trialLevel_choice = svm_Meta.Posterior_2d_n11n;
elseif if_resample_meta == 1
    meta_trialLevel_choice = svm_Meta.Posterior_2d_mean_valid;
end


%%  Compute meta-memory on other trials (noChoice trials)
if if_compute_train == 1
    if if_resample_meta == 0
        x_raw = F_dff_decisionBin1';
        x = x_raw(noChoiceBoolIndex_validLength,:);
        
        temp_Mdl_CV_binary = svm_Meta.temp_Mdl_CV_binary;
        temptemp_Posterior_2d_kfold = zeros(KFold_num_meta,size(x,1));
        for tempk=1:KFold_num_meta
            [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary.Trained{tempk},x);% Very time-consuming!!!
            if size(tempPosterior,2) == 1
                tempPosterior_2 = tempPosterior(:,1);
            else
                tempPosterior_2 = tempPosterior(:,2);
            end
            temptemp_Posterior_2d_kfold(tempk,:) = tempPosterior_2;
        end
        temptemp_Posterior_2d_kfoldMean = squeeze(mean(temptemp_Posterior_2d_kfold,1));
        temptemp_Posterior_2d_n11n = temptemp_Posterior_2d_kfoldMean';
        
        meta_trialLevel_noChoice = temptemp_Posterior_2d_n11n;
        
    elseif if_resample_meta == 1
        
        x_raw = F_dff_decisionBin1';
        x = x_raw(noChoiceBoolIndex_validLength,:);
        
        temptemp_Posterior_2d_n11n_resampled = nan(size(x,1),resampleIterCount_meta);
        
        %for tempIter=1:resampleIterCount_meta
        parfor tempIter=1:resampleIterCount_meta
            temp_svm_Meta = svm_Meta.temp_svm_Meta_resampled{tempIter};
            
            temp_Mdl_CV_binary = temp_svm_Meta.temp_Mdl_CV_binary;
            temptemp_Posterior_2d_kfold = zeros(KFold_num_meta,size(x,1));
            for tempk=1:KFold_num_meta
                [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary.Trained{tempk},x);% Very time-consuming!!!
                if size(tempPosterior,2) == 1
                    tempPosterior_2 = tempPosterior(:,1);
                else
                    tempPosterior_2 = tempPosterior(:,2);
                end
                temptemp_Posterior_2d_kfold(tempk,:) = tempPosterior_2;
            end
            temptemp_Posterior_2d_kfoldMean = squeeze(mean(temptemp_Posterior_2d_kfold,1));
            temptemp_Posterior_2d_n11n = temptemp_Posterior_2d_kfoldMean';
            
            temptemp_Posterior_2d_n11n_resampled(:,tempIter) = temptemp_Posterior_2d_n11n;
        end
        
        temptemp_Posterior_2d_n11n_resampledMean = mean(temptemp_Posterior_2d_n11n_resampled,2);
        
        meta_trialLevel_noChoice = temptemp_Posterior_2d_n11n_resampledMean;
    end
    
end


%% Compute meta-memory on all trials
temp_meta_trialLevel = nan(length(temp_seqIndex),1);
temp_meta_trialLevel(choiceBoolIndex_validLength) = meta_trialLevel_choice;
temp_meta_trialLevel(noChoiceBoolIndex_validLength) = meta_trialLevel_noChoice;


%% Compute meta_seqLevel
meta_seqLevel_choice = zeros(sum(numSeq(1:3)),1);
meta_seqLevel_choiceMemory = zeros(sum(numSeq(1:3)),1);
meta_seqLevel_choiceOffload = zeros(sum(numSeq(1:3)),1);
meta_seqLevel_noChoice = zeros(sum(numSeq(1:3)),1);
meta_seqLevel_all = zeros(sum(numSeq(1:3)),1);

for tempi=1:sum(numSeq(1:3))
    temptempBoolIndex = temp_seqIndex==tempi;
    
    temptempBoolIndex_choice = temptempBoolIndex & choiceBoolIndex;
    temptempBoolIndex_choiceMemory = temptempBoolIndex & choiceMemoryBoolIndex;
    temptempBoolIndex_choiceOffload = temptempBoolIndex & choiceOffloadBoolIndex;
    temptempBoolIndex_noChoice = temptempBoolIndex & (~choiceBoolIndex);
    
    meta_seqLevel_choice(tempi) = mean(temp_meta_trialLevel(temptempBoolIndex_choice),'omitnan');
    meta_seqLevel_choiceMemory(tempi) = mean(temp_meta_trialLevel(temptempBoolIndex_choiceMemory),'omitnan');
    meta_seqLevel_choiceOffload(tempi) = mean(temp_meta_trialLevel(temptempBoolIndex_choiceOffload),'omitnan');
    meta_seqLevel_noChoice(tempi) = mean(temp_meta_trialLevel(temptempBoolIndex_noChoice),'omitnan');
    meta_seqLevel_all(tempi) = mean(temp_meta_trialLevel(temptempBoolIndex),'omitnan');
end
meta_seqLevel_allCondition = [meta_seqLevel_choice,meta_seqLevel_choiceMemory,meta_seqLevel_choiceOffload,meta_seqLevel_noChoice];

temp1 = isnan(meta_seqLevel_allCondition);
temp2 = sum(temp1,2)==0;

% meta_seqLevelMean_allCondition = mean(meta_seqLevel_allCondition,1,'omitnan');
meta_seqLevelMean_allCondition = mean(meta_seqLevel_allCondition(temp2,:),1,'omitnan');


%%
if if_trainMeta_0baseline_1delay1 == 0
    meta_trialLevel_baseline = temp_meta_trialLevel;
    
    
    meta_seqLevel_choice_baseline = meta_seqLevel_choice;
    meta_seqLevel_choiceMemory_baseline = meta_seqLevel_choiceMemory;
    meta_seqLevel_choiceOffload_baseline = meta_seqLevel_choiceOffload;
    meta_seqLevel_noChoice_baseline = meta_seqLevel_noChoice;
    meta_seqLevel_all_baseline = meta_seqLevel_all;
    
    
elseif if_trainMeta_0baseline_1delay1 == 1
    meta_trialLevel_delay1 = temp_meta_trialLevel;
    meta_trialLevel = temp_meta_trialLevel;
    
    
    meta_seqLevel_choice_delay1 = meta_seqLevel_choice;
    meta_seqLevel_choiceMemory_delay1 = meta_seqLevel_choiceMemory;
    meta_seqLevel_choiceOffload_delay1 = meta_seqLevel_choiceOffload;
    meta_seqLevel_noChoice_delay1 = meta_seqLevel_noChoice;
    meta_seqLevel_all_delay1 = meta_seqLevel_all;
    
end


%% meta_trialLevelEvidence
if if_plot_meta_trialLevelEvidence == 1 && if_trainMeta_0baseline_1delay1 == 1
    
    meta_trialLevel;
    
    
    mildSeqBoolIndex_seqLevel = false(1,length(offloadingProb_inOne));
    mildSeqBoolIndex_seqLevel(exampleSeq) = true;
    
    mildSeqBoolIndex = false(size(choiceBoolIndex));
    for tempi=1:length(mildSeqBoolIndex)
        if mildSeqBoolIndex_seqLevel(seqIndex(tempi)) == true
            mildSeqBoolIndex(tempi) = true;
        end
    end
    
    choiceBoolIndex;
    choiceMemoryBoolIndex;
    choiceOffloadBoolIndex;
    
    choiceMemoryBoolIndex_mildSeq = choiceMemoryBoolIndex & mildSeqBoolIndex;
    choiceOffloadBoolIndex_mildSeq = choiceOffloadBoolIndex & mildSeqBoolIndex;
    
    meta_trialLevel_mildSeq_choiceMemory = meta_trialLevel(choiceMemoryBoolIndex_mildSeq);
    meta_trialLevel_mildSeq_choiceOffload = meta_trialLevel(choiceOffloadBoolIndex_mildSeq);
    
    
    meta_trialLevel_allSeq_choiceOffloadMean = nan(sum(numSeq(1:valid_length)),1);
    meta_trialLevel_allSeq_choiceMemoryStimuliMean = nan(sum(numSeq(1:valid_length)),1);
    
    for tempi=1:sum(numSeq(1:valid_length))
        tempSeqBoolIndex = seqIndex==tempi;
        
        meta_trialLevel_allSeq_choiceOffloadMean(tempi) = ...
            mean(meta_trialLevel(choiceOffloadBoolIndex&tempSeqBoolIndex));
        
        meta_trialLevel_allSeq_choiceMemoryStimuliMean(tempi) = ...
            mean(meta_trialLevel(choiceMemoryBoolIndex&tempSeqBoolIndex));
    end
    d1 = meta_trialLevel_allSeq_choiceOffloadMean;
    d2 = meta_trialLevel_allSeq_choiceMemoryStimuliMean;
    
    [~,temp_p_choiceMeta_stimuli] = ttest(d1,d2);
    
    
    
    
    metaDecoderThreshold;
    
    
    trialBoolIndex_metaLow_choiceMemory = (meta_trialLevel<=metaDecoderThreshold) & choiceMemoryBoolIndex';
    trialBoolIndex_metaHigh_choiceMemory = (meta_trialLevel>metaDecoderThreshold) & choiceMemoryBoolIndex';
    
    trialBoolIndex_metaLow_choiceOffload = (meta_trialLevel<=metaDecoderThreshold) & choiceOffloadBoolIndex';
    trialBoolIndex_metaHigh_choiceOffload = (meta_trialLevel>metaDecoderThreshold) & choiceOffloadBoolIndex';
    
    
    seqIndex;
    
    seqRProb_metaLow_choice = nan(sum(numSeq(1:valid_length)),1);
    seqRProb_metaHigh_choice = nan(sum(numSeq(1:valid_length)),1);
    
    for tempi=1:sum(numSeq(1:valid_length))
        % Stimuli-labeled trial
        tempTrialBoolIndex_targetSeq = seqIndex'==tempi;
        
        tempTrialBoolIndex_targetSeq_low_choiceMemory = ...
            tempTrialBoolIndex_targetSeq & trialBoolIndex_metaLow_choiceMemory;
        tempTrialBoolIndex_targetSeq_low_choiceOffload = ...
            tempTrialBoolIndex_targetSeq & trialBoolIndex_metaLow_choiceOffload;
        
        tempRProb_low = sum(tempTrialBoolIndex_targetSeq_low_choiceOffload)/...
            (sum(tempTrialBoolIndex_targetSeq_low_choiceMemory)+sum(tempTrialBoolIndex_targetSeq_low_choiceOffload));
        seqRProb_metaLow_choice(tempi) = tempRProb_low;
        
        tempTrialBoolIndex_targetSeq_high_choiceMemory = ...
            tempTrialBoolIndex_targetSeq & trialBoolIndex_metaHigh_choiceMemory;
        tempTrialBoolIndex_targetSeq_high_choiceOffload = ...
            tempTrialBoolIndex_targetSeq & trialBoolIndex_metaHigh_choiceOffload;
        
        tempRProb_high = sum(tempTrialBoolIndex_targetSeq_high_choiceOffload)/...
            (sum(tempTrialBoolIndex_targetSeq_high_choiceMemory)+sum(tempTrialBoolIndex_targetSeq_high_choiceOffload));
        seqRProb_metaHigh_choice(tempi) = tempRProb_high;
        
        if tempi==41
            a = 1;
        end
    end
    
    b1 = seqRProb_metaLow_choice;
    b2 = seqRProb_metaHigh_choice;
    
    
    b1_mean = mean(b1,'omitnan');
    b2_mean = mean(b2,'omitnan');
    
    
    [~,temp_p_lowHigh] = ttest(b1,b2);
    
    
end


seqPrecision_neuron_choice = nan(sum(numSeq(1:valid_length)),1);
for tempi=1:sum(numSeq(1:valid_length))
    temptempBoolIndex = (seqIndex == tempi) & (~isnan(memoryPrecision_trialLevel))';
    
    temptempBoolIndex2 = temptempBoolIndex & choiceBoolIndex;
    
    temp1 = memoryPrecision_trialLevel(temptempBoolIndex2);
    seqPrecision_neuron_choice(tempi) = mean(temp1,'omitnan');
end
% seqPrecision_neuron_choice = seqPrecision_neuron;

%% Plot
if if_plot == 1
    fig = figure('Name','meta_seqLevel','NumberTitle','off');
    %set(gcf,'Position',[10 50 720*0.92 240*0.9*0.98]);
    %set(gcf,'Position',[10 50 720*0.92*0.85 240*0.9*0.98*0.8*1.06]);
     set(gcf,'Position',[10 50 720*0.92*0.85 240*0.9*0.98*0.8*1.06*0.92]);
   
    t = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
    
    t.Title.String = sprintf('%d seqs (length123), %s',...
        length(meta_seqLevel_choice),currentSession_short);
    t.Title.FontSize = 12;
    t.Title.Interpreter = 'none';
    
    
    %% fig, ROC curve for meta-memory docoder
    nexttile
        
    x = svm_choiceMemory_falseAlarm_multi;
    y = svm_choiceMemory_hit_multi;
    
    x = x(end:-1:1);
    y = y(end:-1:1);
    
    AUROC = trapz(x,y);
    fprintf('AUROC = %.3f\n',AUROC);
    
    if if_trainMeta_0baseline_1delay1 == 0
        AUROC_meta_baseline = AUROC;
    elseif if_trainMeta_0baseline_1delay1 == 1
        AUROC_meta_delay1 = AUROC;
    end
    
    plot(x,y,'color',[0.25 0.25 0.25],'linewidth',1);    
    hold on
    
    plot(0:1,0:1,'--','Color','k');
    
    h = fill([x,1:-1:0],[y,0 0],[0.25 0.25 0.25]);
    
    set(h,'edgealpha',0,'facealpha',0.3)
    
    set(gca,'linewidth',1.5)
    ylim([0 1]);
    xlim([0 1]);
    set(gca, 'FontSize', 8)
    
    xticks([0 1]);
    yticks([0 1]);
    
    set(gca,'box','off');
    xlabel({'False Alarm Rate'}, 'FontSize', 9, 'FontWeight', 'normal');
    ylabel({'Hit Rate'}, 'FontSize', 9, 'FontWeight', 'normal');
    
    text(0.275,0.25,sprintf('AUROC = %.3f',AUROC),'fontsize',7,'FontWeight','normal');
    text(0.275,0.12,sprintf('(95%% = %.3f)',AUROC_shuffled_prctile),'fontsize',7,'FontWeight','normal');
    
    title(sprintf('AUROC'), 'FontSize', 9, 'FontWeight', 'normal');
    
    
    %% fig, meta_seqLevel
    for loopCount=1:2
        nexttile
        if loopCount == 1
            %tempStr = 'Free choice (training set)';
            %tempStr = 'Free-choice';
            tempStr = 'Across seqs';
            x = meta_seqLevel_choice;
        elseif loopCount == 2
            %tempStr = 'Force to test (testing set)';
            tempStr = 'Forced-to-test';
            x = meta_seqLevel_noChoice;
        end
        y = offloadingProb_inOne(1:sum(numSeq(1:3)))';
        
        %[r,p] = corr(x,y);
        [r,p] = corr(x(~isnan(x)),y(~isnan(x)));
        
        mdl = fitglm(x,y);
        
        tempBoolIndex = ~isnan(x);
        
        h = [];
        for tempi=1:3
            temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
            tempIndex = find(tempBoolIndex(temp_range)==true);
            temp_range2 = temp_range(tempIndex);
            
            %temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
            temp_size = 10;
            temp_h = scatter(x(temp_range2), y(temp_range2), ...
                temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
            h = [h temp_h]; %#ok<*AGROW>
        end
        
        [temp_xmin,temp_xmax] = bounds(x);
        [temp_ymin,temp_ymax] = bounds(y);
        
        x_fit = temp_xmin:0.001:temp_xmax;
        y_fit = predict(mdl,x_fit')';
        
        temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        hold on
        h = [h temp_h];
        
        xlim([0 1]);
        ylim([0 1]);
        
        xticks([0 1]);
        yticks([0 1]);
        
        set(gca,'linewidth',1.5)
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
        
        
        tempTxt = 'p > 0.05';
        if p < 0.05
            tempTxt = 'p < 0.05';
        end
        if p < 0.01
            tempTxt = 'p < 0.01';
        end
        if p < 0.001
            tempTxt = 'p < 0.001';
        end
        if if_trainMeta_0baseline_1delay1 == 1
            text(0.55,0.92,sprintf('r = %.3f',r), 'fontsize',7.5,'FontWeight','normal');
            text(0.55,0.80,sprintf('%s',tempTxt), 'fontsize',7.5,'FontWeight','normal');
        elseif if_trainMeta_0baseline_1delay1 == 0
            %text(0.10,0.32,sprintf('r = %.3f',r), 'fontsize',7.5,'FontWeight','normal');
            %text(0.10,0.22,sprintf('%s',tempTxt), 'fontsize',7.5,'FontWeight','normal');
            text(0.58,0.32,sprintf('r = %.3f',r), 'fontsize',7.5,'FontWeight','normal');
            text(0.58,0.20,sprintf('%s',tempTxt), 'fontsize',7.5,'FontWeight','normal');
        end
        
        
        %temp_title = title(sprintf('%s, r=%.3f, p=%.3f',tempStr,r,p), 'FontSize', 12, 'FontWeight', 'bold');
        %temp_title = title(sprintf(' %s \n r=%.3f, p=%.3f ',tempStr,r,p), 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('%s',tempStr), 'FontSize', 9, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
        ylabel(sprintf('Offloading rate'), 'FontSize', 9, 'FontWeight', 'normal');
        
    end
    
    
    %% fig, meta_seqLevel
    nexttile
    tempStr = 'Correlation';
    
    x_raw = seqPrecision_neuron_choice;
    %x = rescale(x_raw,0.01,0.99);
    x = x_raw;
    y = meta_seqLevel_choice;
    
    
    [r,p] = corr(x(~isnan(x)),y(~isnan(x)));
    
    mdl = fitglm(x,y);
    
    tempBoolIndex = ~isnan(x);
    
    h = [];
    for tempi=1:3
        temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
        tempIndex = find(tempBoolIndex(temp_range)==true);
        temp_range2 = temp_range(tempIndex);
        
        temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
        temp_h = scatter(x(temp_range2), y(temp_range2), ...
            temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        h = [h temp_h]; %#ok<*AGROW>
    end
    
    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y);
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    h = [h temp_h];
    
    xlim([temp_xmin-(temp_xmax-temp_xmin)*0.08 temp_xmax+(temp_xmax-temp_xmin)*0.08]);
    ylim([0 1]);
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    text(0.55,0.25,sprintf('r = %.3f',r), 'fontsize',8,'FontWeight','normal');
    text(0.55,0.15,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
    
    if if_monkey_D0_Z1 == 0 && currentSessionIndex_AB == 8
        xticks([floor(temp_xmin*10)/10 floor(temp_xmax*10)/10]);
        yticks([floor(temp_ymin*10)/10 floor(temp_ymax*10)/10]);
    end
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
    
    %temp_title = title(sprintf(' %s \n r=%.3f, p=%.3f ',tempStr,r,p), 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('%s',tempStr), 'FontSize', 9, 'FontWeight', 'normal');
    temp_title.Interpreter = 'none';
    
    
    %% Meta-WM in different conditions 
    if true
        fig = figure('Name','meta_seqLevel','NumberTitle','off');
        set(gcf,'Position',[10 300 720*0.92*0.85 240*0.9*0.98*0.8*1.06*0.92*0.9]);
        
        t = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
                        
        
        %% fig, meta_seqLevel
        for loopCount=1:4
            nexttile
            if loopCount == 1
                tempStr = 'Choice';
                x = meta_seqLevel_choice;
            elseif loopCount == 2
                tempStr = 'Choice-memory';
                x = meta_seqLevel_choiceMemory;            
            elseif loopCount == 3
                tempStr = 'Choice-offload';
                x = meta_seqLevel_choiceOffload;      
            elseif loopCount == 4
                tempStr = 'Forced-to-test';
                x = meta_seqLevel_noChoice;                
            end
            y = offloadingProb_inOne(1:sum(numSeq(1:3)))';
            
            %[r,p] = corr(x,y);
            [r,p] = corr(x(~isnan(x)),y(~isnan(x)));
            
            mdl = fitglm(x,y);
            
            tempBoolIndex = ~isnan(x);
            
            h = [];
            for tempi=1:3
                temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
                tempIndex = find(tempBoolIndex(temp_range)==true);
                temp_range2 = temp_range(tempIndex);
                
                %temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
                temp_size = 10;
                temp_h = scatter(x(temp_range2), y(temp_range2), ...
                    temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
                hold on
                h = [h temp_h]; %#ok<*AGROW>
            end
            
            [temp_xmin,temp_xmax] = bounds(x);
            [temp_ymin,temp_ymax] = bounds(y);
            
            x_fit = temp_xmin:0.001:temp_xmax;
            y_fit = predict(mdl,x_fit')';
            
            temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
            hold on
            h = [h temp_h];
            
            xlim([0 1]);
            ylim([0 1]);
            
            xticks([0 1]);
            yticks([0 1]);
            
            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            xlabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
            
            
            tempTxt = 'p > 0.05';
            if p < 0.05
                tempTxt = 'p < 0.05';
            end
            if p < 0.01
                tempTxt = 'p < 0.01';
            end
            if p < 0.001
                tempTxt = 'p < 0.001';
            end
            if if_trainMeta_0baseline_1delay1 == 1
                text(0.55,0.92,sprintf('r = %.3f',r), 'fontsize',7.5,'FontWeight','normal');
                text(0.55,0.80,sprintf('%s',tempTxt), 'fontsize',7.5,'FontWeight','normal');
            elseif if_trainMeta_0baseline_1delay1 == 0
                text(0.10,0.32,sprintf('r = %.3f',r), 'fontsize',7.5,'FontWeight','normal');
                text(0.10,0.22,sprintf('%s',tempTxt), 'fontsize',7.5,'FontWeight','normal');
            end
            
            
            %temp_title = title(sprintf('%s, r=%.3f, p=%.3f',tempStr,r,p), 'FontSize', 12, 'FontWeight', 'bold');
            %temp_title = title(sprintf(' %s \n r=%.3f, p=%.3f ',tempStr,r,p), 'FontSize', 12, 'FontWeight', 'bold');
            temp_title = title(sprintf('%s',tempStr), 'FontSize', 9, 'FontWeight', 'normal');
            temp_title.Interpreter = 'none';
            
            ylabel(sprintf('Offloading rate'), 'FontSize', 9, 'FontWeight', 'normal');
            
        end
        
    end
    
    
    %% fig, meta_seqLevel
    fig = figure('Name','meta_seqLevel','NumberTitle','off');
    %set(gcf,'Position',[610 390 240*1.3*1.08 240*0.75*0.98]);
    %set(gcf,'Position',[610 390 337 176*1.065*1.2]);
    %     if if_memoryPrecision_accuracy0_sigma1 == 0
    %set(gcf,'Position',[610 390 337*0.8*0.99 176*1.065*1.2*0.85]);
    %set(gcf,'Position',[610 390 337*0.8*0.99*1.15 176*1.065*1.2*0.85]);
    %set(gcf,'Position',[610 390 337*0.8*0.99*1.15*0.57 176*1.065*1.2*0.85*0.80]);
    set(gcf,'Position',[610 390 337*0.8*0.99*1.15*0.57*0.79 176*1.065*1.2*0.85*0.80]);
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    %     elseif if_memoryPrecision_accuracy0_sigma1 == 1
    %         set(gcf,'Position',[610 390 226*0.75 186]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %     end
    
    %tempStr = 'Correlation';
    
    x_raw = seqPrecision_neuron_choice;
    %x = rescale(x_raw,0.01,0.99);
    x = x_raw;
    y = meta_seqLevel_choice;
    
    
    [r,p] = corr(x(~isnan(x)),y(~isnan(x)));
    
    mdl = fitglm(x,y);
    
    tempBoolIndex = ~isnan(x);
    
    h = [];
    for tempi=1:3
        temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
        tempIndex = find(tempBoolIndex(temp_range)==true);
        temp_range2 = temp_range(tempIndex);
        
        %temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
        temp_size = 10;
        temp_h = scatter(x(temp_range2), y(temp_range2), ...
            temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        h = [h temp_h]; %#ok<*AGROW>
    end
    
    %[temp_xmin,temp_xmax] = bounds(x);
    %[temp_ymin,temp_ymax] = bounds(y);
    
    temp_xmin = 0;
    temp_xmax = 1;
    temp_ymin = 0;
    temp_ymax = 1;
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    h = [h temp_h];
    
    xlim([temp_xmin-(temp_xmax-temp_xmin)*0.08 temp_xmax+(temp_xmax-temp_xmin)*0.08]);
    %     if if_memoryPrecision_accuracy0_sigma1 == 0
    %ylim([0 1]);
    ylim([temp_ymin-(temp_ymax-temp_ymin)*0.08 temp_ymax+(temp_ymax-temp_ymin)*0.08]);
    %     elseif if_memoryPrecision_accuracy0_sigma1 == 1
    %         ylim([temp_ymin-(temp_ymax-temp_ymin)*0.08 temp_ymax+(temp_ymax-temp_ymin)*0.08]);
    %     end
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    %if if_memoryPrecision_accuracy0_sigma1 == 0
    %    text(0.76,0.25,sprintf('r=%.3f',r), 'fontsize',9,'FontWeight','normal');
    %    text(0.76,0.12,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
    %elseif if_memoryPrecision_accuracy0_sigma1 == 1
    %    text(temp_xmin+(temp_xmax-temp_xmin)*0.63,temp_ymin+(temp_ymax-temp_ymin)*0.18,sprintf('r=%.3f',r), 'fontsize',9,'FontWeight','normal');
    %    text(temp_xmin+(temp_xmax-temp_xmin)*0.63,temp_ymin+(temp_ymax-temp_ymin)*0.05,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
    %end
    
    %text(0.03,0.98,sprintf('r = %.3f',r), 'fontsize',7.5,'FontWeight','normal');
    %text(0.03,0.85,sprintf('%s',tempTxt), 'fontsize',7.5,'FontWeight','normal');
    text(0.53,0.18,sprintf('r = %.3f',r), 'fontsize',7.5,'FontWeight','normal');
    text(0.53,0.05,sprintf('%s',tempTxt), 'fontsize',7.5,'FontWeight','normal');
    
    
    if false
        if if_monkey_D0_Z1 == 0 && currentSessionIndex_AB == 8
            if if_memoryPrecision_accuracy0_sigma1 == 0
                xticks([floor(temp_xmin*10)/10 floor(temp_xmax*10)/10]);
                yticks([floor(temp_ymin*10)/10 floor(temp_ymax*10)/10]);
            end
        end
    end
    
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
    
    %temp_title = title(sprintf(' %s \n r=%.3f, p=%.3f ',tempStr,r,p), 'FontSize', 12, 'FontWeight', 'bold');
    %temp_title = title(sprintf('r=%.3f, p=%.3f',r,p), 'FontSize', 12, 'FontWeight', 'bold');
    %     if if_memoryPrecision_accuracy0_sigma1 == 0
    %temp_title = title(sprintf('Correlation'), 'FontSize', 10, 'FontWeight', 'normal');
    %temp_title = title(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 8.5, 'FontWeight', 'normal');
    temp_title = title(sprintf('Across seqs'), 'FontSize', 9, 'FontWeight', 'normal');
    %     elseif if_memoryPrecision_accuracy0_sigma1 == 1
    %         temp_title = title(sprintf('Neuronal correlation'), 'FontSize', 9, 'FontWeight', 'normal');
    %     end
    temp_title.Interpreter = 'none';
    
    
    if false
        
        %% fig, meta_seqLevel
        fig = figure('Name','meta_seqLevel','NumberTitle','off'); %#ok<*UNRCH>
        %set(gcf,'Position',[610 390 337*0.8*0.99*1.15 176*1.065*1.2*0.85]);
        set(gcf,'Position',[610 390 337*0.8*0.99*1.15*0.6*0.7*1.1*1.01 176*1.065*1.2*0.85*0.6*0.95*1.1]);
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        
        
        tempStr = 'Correlation';
        
        x_raw = seqPrecision_neuron_choice;
        %x = rescale(x_raw,0.01,0.99);
        x = x_raw;
        y = meta_seqLevel_choice;
        
        
        [r,p] = corr(x(~isnan(x)),y(~isnan(x)));
        
        mdl = fitglm(x,y);
        
        tempBoolIndex = ~isnan(x);
        
        h = [];
        for tempi=1:3
            temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
            tempIndex = find(tempBoolIndex(temp_range)==true);
            temp_range2 = temp_range(tempIndex);
            
            %temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
            temp_size = ((tempi.^3)*0.8 + 3) .* ones(1, length(temp_range2));
            temp_h = scatter(x(temp_range2), y(temp_range2), ...
                temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
            h = [h temp_h]; %#ok<*AGROW>
        end
        
        [temp_xmin,temp_xmax] = bounds(x);
        [temp_ymin,temp_ymax] = bounds(y);
        
        x_fit = temp_xmin:0.001:temp_xmax;
        y_fit = predict(mdl,x_fit')';
        
        temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        hold on
        h = [h temp_h];
        
        xlim([temp_xmin-(temp_xmax-temp_xmin)*0.08 temp_xmax+(temp_xmax-temp_xmin)*0.08]);
        %         if if_memoryPrecision_accuracy0_sigma1 == 0
        ylim([0 1]);
        %         elseif if_memoryPrecision_accuracy0_sigma1 == 1
        %             ylim([temp_ymin-(temp_ymax-temp_ymin)*0.08 temp_ymax+(temp_ymax-temp_ymin)*0.08]);
        %         end
        
        tempTxt = 'p>0.05';
        if p < 0.05
            tempTxt = 'p<0.05';
        end
        if p < 0.01
            tempTxt = 'p<0.01';
        end
        if p < 0.001
            tempTxt = 'p<0.001';
        end
        
        if if_monkey_D0_Z1 == 0 && currentSessionIndex_AB == 8
            %             if if_memoryPrecision_accuracy0_sigma1 == 0
            xticks([floor(temp_xmin*10)/10 floor(temp_xmax*10)/10]);
            yticks([floor(temp_ymin*10)/10 floor(temp_ymax*10)/10]);
            %             end
        end
        
        
        set(gca,'linewidth',1.5)
        set(gca, 'FontSize', 10)
        set(gca,'box','off');% 取消右、上边框
        xlabel(sprintf('WM strength'), 'FontSize', 10, 'FontWeight', 'normal');
        ylabel(sprintf('Meta-WM'), 'FontSize', 10, 'FontWeight', 'normal');
        
        temp_title = title(sprintf('r=%.3f, %s',r,tempTxt), 'FontSize', 7.5, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
        
    end
    
    
    
    if if_plot_meta_trialLevelEvidence == 1 && if_trainMeta_0baseline_1delay1 == 1
        
        %% meta_trialLevelEvidence in choice trials
        fig = figure('Name','meta pdf (mildSeq)','NumberTitle','off');
        
        if if_trialEvidenceExampleSeq_histogram0_violinPlot1 == 0
            %set(gcf,'Position',[10 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355*1.3*1.1 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            set(gcf,'Position',[10 450 360*1.1*0.97 240*0.8*1.03]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            t = tiledlayout(1,2,'TileSpacing','compact','Padding','loose');
            
            
            %nexttile([1 2])
            nexttile
            
            x1 = meta_trialLevel_mildSeq_choiceMemory;
            x2 = meta_trialLevel_mildSeq_choiceOffload;
            
            [~,temp_p] = ttest2(x1,x2);
            
            
            h_NumBins = 8;%10
            
            x = x1;
            h1 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
            hold on
            h1.NumBins = h_NumBins;
            h1.EdgeColor = color_choiceMemory;
            
            x = x2;
            h2 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
            hold on
            h2.NumBins = h_NumBins;
            h2.EdgeColor = color_choiceOffload;
            
            
            y1 = h1.Values;
            y2 = h2.Values;
            
            x_min = 0;
            x_max = 1;
            
            [y1_min,y1_max] = bounds(y1);
            [y2_min,y2_max] = bounds(y2);
            y_min = min([y1_min,y2_min]);
            y_max = max([y1_max,y2_max]);
            
            if if_plot_additionalSmooth == 1
                n=100;
                n=2^ceil(log2(n)); % round up n to the next power of 2;
                
                [pdf1,xmesh1,~] = ksdensity(x1,'NumPoints',n,'Function','pdf');
                plot(xmesh1,pdf1*sum(y1)*h1.BinWidth,':','LineWidth',1.5,'color',color_choiceMemory);
                hold on
                
                [pdf2,xmesh2,~] = ksdensity(x2,'NumPoints',n,'Function','pdf');
                plot(xmesh2,pdf2*sum(y2)*h2.BinWidth,':','LineWidth',1.5,'color',color_choiceOffload);
                hold on
            end
            
            
            plot([metaDecoderThreshold metaDecoderThreshold],[y_min y_max],...
                'LineWidth',1,'color',[0.9290 0.6940 0.1250]);
            hold on
            
            le = legend('Choice-memory','Choice-offload',...
                'Location','northeast','fontsize',7);
            le.ItemTokenSize = ones(1,2)*8;
            
            
            tempTxt = sprintf('');
            if temp_p < 0.001
                tempTxt = sprintf('***');
            elseif temp_p < 0.01
                tempTxt = sprintf('**');
            elseif temp_p < 0.05
                tempTxt = sprintf('*');
            end
            text(0.4,y_min+(y_max-y_min)*0.7,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            
            set(gca,'linewidth',1.5)
            xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
            ylim([y_min y_max+(y_max-y_min)*0.45]);%0.1
            %xticks([0 1]);
            set(gca, 'FontSize', 10)
            set(gca,'box','off');% 取消右、上边框
            xlabel('Meta-WM', 'FontSize', 10, 'FontWeight', 'normal');
            ylabel('Trial count', 'FontSize', 10, 'FontWeight', 'normal');
            
            
            temp_trialNum = length(x1)+length(x2);
            
            %temp_title = title(sprintf('Free-choice trial\nSeq %d (%d trials)',...
            %    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',10,'FontWeight','normal');
            temp_title = title(sprintf('Free-choice\nSeq %d (%d trials)',...
                seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',10,'FontWeight','normal');
            temp_title.Interpreter = 'none';
            
        elseif if_trialEvidenceExampleSeq_histogram0_violinPlot1 == 1
            %set(gcf,'Position',[10 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355 200*2*0.85*0.9*0.98*0.985]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %t = tiledlayout(2,2,'TileSpacing','compact','Padding','loose');
            
            %set(gcf,'Position',[10 450 355*1.5*1.05 295]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355*1.5*1.05*0.97 295]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355*1.5*1.05*0.97*0.84 295*1.06]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355*1.5*1.05*0.97*0.84*0.74 295*1.06]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            set(gcf,'Position',[10 450 355*1.5*1.05*0.97*0.84*0.74 295*1.06*0.88]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            %t = tiledlayout(2,3,'TileSpacing','compact','Padding','loose');
            t = tiledlayout(2,3,'TileSpacing','loose','Padding','loose');
            
            
            %nexttile([1 2])
            nexttile
            
            x1 = meta_trialLevel_mildSeq_choiceOffload;
            x2 = meta_trialLevel_mildSeq_choiceMemory;
            
            x1 = x1(~isnan(x1));
            x2 = x2(~isnan(x2));
            
            [~,temp_p] = ttest2(x1,x2);
            
            temp_1 = x1;
            temp_2 = x2;
            
            if ~isempty(temp_1) && ~isempty(temp_2) == 1
                
                temp_y_min = min([temp_1;temp_2]);
                temp_y_max = max([temp_1;temp_2]);
                
                temp_y_min = min([temp_y_min,0]);
                temp_y_max = max([temp_y_max,1]);
                
                temp_data = [temp_1;temp_2];
                
                g1 = repmat({'A'},length(temp_1),1);
                g2 = repmat({'B'},length(temp_2),1);
                
                temp_label = [g1;g2];
                
                temptemp_color1 = [1 1 1]*0.5;
                temptemp_color2 = repmat(temptemp_color1, 2, 1);
                
                %violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                %    'GroupOrder',[{'A'};{'B'}]);
                
                h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                    'GroupOrder',[{'A'};{'B'}]);
                h(1).ViolinPlot.FaceAlpha = 0.1;
                h(2).ViolinPlot.FaceAlpha = 0.1;
                
                
                tempTxt = sprintf('');
                if temp_p < 0.001
                    tempTxt = sprintf('***');
                elseif temp_p < 0.01
                    tempTxt = sprintf('**');
                elseif temp_p < 0.05
                    tempTxt = sprintf('*');
                end
                text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                    'HorizontalAlignment','center');
                
                set(gca,'linewidth',1.5)
                xlim([0.15 2.65])
                ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
                set(gca, 'FontSize', 8)
                
                %xtl = ["Choice-memory"; "Choice-offload"];
                xtl = ["Choice-offload";"Choice-memory"];
                xt=get(gca,'XTick');
                yt=get(gca,'YTick');
                xtext_xp=xt;
                %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.325;%0.56,0.4
                xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.225;%0.56,0.4
                %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25,10
                text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25,10
                set(gca,'xticklabel','');
                
                yticks([0 1]);
                
                set(gca,'box','off');% 取消右、上边框
                
                ylabel('Meta-WM', 'FontSize', 9, 'FontWeight', 'normal');
                temp_trialNum = length(x1)+length(x2);
                %temp_title = title(sprintf('Free-choice trial\nSeq %d (%d trials)',...
                %    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',10,'FontWeight','normal');
                temp_title = title(sprintf('Free-choice\nSeq %d (%d trials)',...
                    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',9,'FontWeight','normal');
                temp_title.Interpreter = 'none';
                
            end
            
        end
        
        
        % Compare meta-memory of choice-offload and choice-memory trials
        nexttile
        
        temp_p = temp_p_choiceMeta_stimuli;
        temptempBoolIndex = (~isnan(d1)) & (~isnan(d2));
        temp_1 = d1(temptempBoolIndex);
        temp_2 = d2(temptempBoolIndex);
        
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        
        temp_y_min = min([temp_1;temp_2]);
        
        %temp_y_min = min([temp_y_min,0.06]);
        
        temp_y_max = max([temp_1;temp_2]);
        
        temp_y_min = min([temp_y_min,0]);
        temp_y_max = max([temp_y_max,1]);
        
        if if_trialEvidenceAllSeq_violinplot0_pairline1 == 0
            temp_data = [temp_1;temp_2];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            
            temp_label = [g1;g2];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 2, 1);
            
            h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'}]);
            h(1).ViolinPlot.FaceAlpha = 0.1;
            h(2).ViolinPlot.FaceAlpha = 0.1;
            
        elseif if_trialEvidenceAllSeq_violinplot0_pairline1 == 1
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);%30
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);%30
            hold on
            
        end
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.12,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.85])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.20]);
        set(gca, 'FontSize', 8)
        
        xticks([1 2]);
        
        %xtl = ["Choice-offload"; "Choice-memory"];
        xtl = ["Offload"; "Memory"];        
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.345;%0.56,0.4
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.19;%0.56,0.4
        %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%25
        set(gca,'xticklabel','');
        
        
        yticks([0 1]);
        
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Meta-WM', 'FontSize', 9, 'FontWeight', 'normal');
        title(sprintf('Free-choice\nAll seqs'),'fontsize',9);
        
        
        
        % Compare offloading rate of lowMeta & highMeta
        nexttile
        
        temp_p = temp_p_lowHigh;
        temptempBoolIndex = (~isnan(b1)) & (~isnan(b2));
        temp_1 = b1(temptempBoolIndex);
        temp_2 = b2(temptempBoolIndex);
        
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
        
        if if_trialEvidenceAllSeq_violinplot0_pairline1 == 0
            temp_data = [temp_1;temp_2];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            
            temp_label = [g1;g2];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 2, 1);
            
            h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'}]);
            h(1).ViolinPlot.FaceAlpha = 0.1;
            h(2).ViolinPlot.FaceAlpha = 0.1;
            
        elseif if_trialEvidenceAllSeq_violinplot0_pairline1 == 1
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            scatter(1*ones(1,length(temp_1)),temp_1,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            
        end
        
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.65])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
        set(gca, 'FontSize', 8)
        
        xticks([1 2]);
        xtl = ["Low-meta"; "High-meta"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.205;%0.56,0.4
        %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%25
        set(gca,'xticklabel','');
        
        
        yticks([0 1]);
        
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Offloading rate', 'FontSize', 9, 'FontWeight', 'normal');
        %title(sprintf('Free-choice trial\nSeq Offloading rate'),'fontsize',9);
        title(sprintf('Free-choice\nAll seqs'),'fontsize',9);
    end
    
    
    if true
        
        [M,I] = max(hit_minus_falseAlarm_multi);
        metaDecoderThreshold = temp_thresholdRange(I);
        hit_minus_falseAlarm_optimal = M;
        
        fig = figure('Name','asd','NumberTitle','off');
        %set(gcf,'Position',[810 400 240*1.15 240*0.7*0.85]);
        %set(gcf,'Position',[810 600 240 143*1.02]);
        set(gcf,'Position',[810 600 240*0.73 143*1.02*0.89]);
        t = tiledlayout(1,1,'TileSpacing','tight','Padding','Compact');
        
        %% Plot 
        nexttile
        
        x = temp_thresholdRange;
        y = hit_minus_falseAlarm_multi;
        
        plot(x,y,'color',[0.25 0.25 0.25],'linewidth',1);
        hold on
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        y_min = min(y_min,0.5);
        
        plot([1 1].*metaDecoderThreshold,[y_min-(y_max-y_min)*0.1 hit_minus_falseAlarm_optimal],':','color',[0.25 0.25 0.25],'linewidth',1);
        hold on
        
        plot([x_min-(x_max-x_min)*0.05 metaDecoderThreshold],[1 1].*hit_minus_falseAlarm_optimal,':','color',[0.25 0.25 0.25],'linewidth',1);
        hold on
        
        %plot([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05],[1 1].*0.5,':','color',[0.25 0.25 0.25],'linewidth',1);
        %hold on
        
        %set(gca,'yscale','log')
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05]);%
        ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);%0.1
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Decision boundary', 'FontSize', 9);
        ylabel('HR - FAR', 'FontSize', 9);
        temp_title = title(sprintf('Meta-WM decoder'),'FontSize',9);
        temp_title.Interpreter = 'none';
        
    end
    
    if if_trainMeta_0baseline_1delay1 == 0
        fig = figure('Name','meta_seqLevel','NumberTitle','off');
        %set(gcf,'Position',[10 50 720*0.92*0.955 240*0.9*0.98]);
        set(gcf,'Position',[10 50 720*0.92*0.955 240*0.9*0.98*0.65*0.7*1.1*1.05]);
        
        
        t = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
        
        %t.Title.String = sprintf('%d seqs (length123), %s',...
        %    length(meta_seqLevel_choice),currentSession_short);
        %t.Title.FontSize = 12;
        %t.Title.Interpreter = 'none';
        
        
        %% fig, ROC curve for meta-memory docoder
        nexttile
        
    x = svm_choiceMemory_falseAlarm_multi;
    y = svm_choiceMemory_hit_multi;
    
    x = x(end:-1:1);
    y = y(end:-1:1);
    
    AUROC = trapz(x,y);
    %fprintf('AUROC = %.3f\n',AUROC);
    
    if if_trainMeta_0baseline_1delay1 == 0
        AUROC_meta_baseline = AUROC;
    elseif if_trainMeta_0baseline_1delay1 == 1
        AUROC_meta_delay1 = AUROC;
    end
    
    plot(x,y,'color',[0.25 0.25 0.25],'linewidth',1);  
        hold on
        
        plot(0:1,0:1,'--','Color','k');
        
        h = fill([x,1:-1:0],[y,0 0],[0.25 0.25 0.25]);
        
        set(h,'edgealpha',0,'facealpha',0.3)
        
        set(gca,'linewidth',1.5)
        ylim([0 1]);
        xlim([0 1]);
        set(gca, 'FontSize', 8)
        
        xticks([0 1]);
        yticks([0 1]);
        
        set(gca,'box','off');
        xlabel({'False Alarm Rate'}, 'FontSize', 9, 'FontWeight', 'normal');
        ylabel({'Hit Rate'}, 'FontSize', 9, 'FontWeight', 'normal');
        
        text(0.33,0.25,sprintf('AUROC = %.3f',AUROC),'fontsize',7.5,'FontWeight','normal');
        %text(0.290,0.13,sprintf('(95%% = %.3f)',AUROC_shuffled_prctile),'fontsize',7.5,'FontWeight','normal');
        
        %title(sprintf('AUROC'), 'FontSize', 9, 'FontWeight', 'normal');
        %title(sprintf('AUROC = %.3f',AUROC), 'FontSize', 9, 'FontWeight', 'normal');
        
        %% fig, meta_seqLevel
        for loopCount=1:2
            nexttile
            if loopCount == 1
                %tempStr = 'Free choice (training set)';
                tempStr = 'Free-choice';
                x = meta_seqLevel_choice;
            elseif loopCount == 2
                %tempStr = 'Force to test (testing set)';
                tempStr = 'Forced-to-test';
                x = meta_seqLevel_noChoice;
            end
            y = offloadingProb_inOne(1:sum(numSeq(1:3)))';
            
            %[r,p] = corr(x,y);
            [r,p] = corr(x(~isnan(x)),y(~isnan(x)));
            
            mdl = fitglm(x,y);
            
            tempBoolIndex = ~isnan(x);
            
            h = [];
            for tempi=1:3
                temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
                tempIndex = find(tempBoolIndex(temp_range)==true);
                temp_range2 = temp_range(tempIndex);
                
                %temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
                temp_size = 10;
                temp_h = scatter(x(temp_range2), y(temp_range2), ...
                    temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
                hold on
                h = [h temp_h]; %#ok<*AGROW>
            end
            
            [temp_xmin,temp_xmax] = bounds(x);
            [temp_ymin,temp_ymax] = bounds(y);
            
            x_fit = temp_xmin:0.001:temp_xmax;
            y_fit = predict(mdl,x_fit')';
            
            temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
            hold on
            h = [h temp_h];
            
            xlim([0 1]);
            ylim([0 1]);
            
            xticks([0 1]);
            yticks([0 1]);
            
            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            xlabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
            
            
            tempTxt = 'p > 0.05';
            if p < 0.05
                tempTxt = 'p < 0.05';
            end
            if p < 0.01
                tempTxt = 'p < 0.01';
            end
            if p < 0.001
                tempTxt = 'p < 0.001';
            end
            
            %if if_trainMeta_0baseline_1delay1 == 1
            %    text(0.60,0.92,sprintf('r=%.3f',r), 'fontsize',8,'FontWeight','normal');%9
            %    text(0.60,0.79,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
            %elseif if_trainMeta_0baseline_1delay1 == 0
            %    text(0.10,0.32,sprintf('r=%.3f',r), 'fontsize',8,'FontWeight','normal');
            %    text(0.10,0.19,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
            %end
            
            if if_resample_meta == 0
                text(0.10,0.32,sprintf('r = %.3f',r), 'fontsize',7.5,'FontWeight','normal');
                text(0.10,0.22,sprintf('%s',tempTxt), 'fontsize',7.5,'FontWeight','normal');
            elseif if_resample_meta == 1
                text(0.025,0.45,sprintf('r = %.3f',r), 'fontsize',7.5,'FontWeight','normal');
                text(0.025,0.22,sprintf('%s',tempTxt), 'fontsize',7.5,'FontWeight','normal');
            end
            
            
            %temp_title = title(sprintf('%s, r=%.3f, p=%.3f',tempStr,r,p), 'FontSize', 12, 'FontWeight', 'bold');
            %temp_title = title(sprintf(' %s \n r=%.3f, p=%.3f ',tempStr,r,p), 'FontSize', 12, 'FontWeight', 'bold');
            %temp_title = title(sprintf('%s',tempStr), 'FontSize', 9, 'FontWeight', 'normal');
            %temp_title.Interpreter = 'none';
            %title(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 9, 'FontWeight', 'normal');
            
            ylabel(sprintf('Offloading rate'), 'FontSize', 9, 'FontWeight', 'normal');
            
        end
        
        
        %% fig, meta_seqLevel
        nexttile
        tempStr = 'Correlation';
        
        x_raw = seqPrecision_neuron_choice;
        %x = rescale(x_raw,0.01,0.99);
        x = x_raw;
        y = meta_seqLevel_choice;
        
        
        [r,p] = corr(x(~isnan(x)),y(~isnan(x)));
        
        mdl = fitglm(x,y);
        
        tempBoolIndex = ~isnan(x);
        
        h = [];
        for tempi=1:3
            temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
            tempIndex = find(tempBoolIndex(temp_range)==true);
            temp_range2 = temp_range(tempIndex);
            
            %temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
            temp_size = 10;
            temp_h = scatter(x(temp_range2), y(temp_range2), ...
                temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
            h = [h temp_h]; %#ok<*AGROW>
        end
        
        [temp_xmin,temp_xmax] = bounds(x);
        [temp_ymin,temp_ymax] = bounds(y);
        
        x_fit = temp_xmin:0.001:temp_xmax;
        y_fit = predict(mdl,x_fit')';
        
        temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        hold on
        h = [h temp_h];
        
        xlim([temp_xmin-(temp_xmax-temp_xmin)*0.08 temp_xmax+(temp_xmax-temp_xmin)*0.08]);
        ylim([0 1]);
        
        tempTxt = 'p > 0.05';
        if p < 0.05
            tempTxt = 'p < 0.05';
        end
        if p < 0.01
            tempTxt = 'p < 0.01';
        end
        if p < 0.001
            tempTxt = 'p < 0.001';
        end
        text(0.55,0.25,sprintf('r = %.3f',r), 'fontsize',8,'FontWeight','normal');
        text(0.55,0.15,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        
        if if_monkey_D0_Z1 == 0 && currentSessionIndex_AB == 8
            xticks([floor(temp_xmin*10)/10 floor(temp_xmax*10)/10]);
            yticks([floor(temp_ymin*10)/10 floor(temp_ymax*10)/10]);
        end
        
        set(gca,'linewidth',1.5)
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
        ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
        
        %temp_title = title(sprintf(' %s \n r=%.3f, p=%.3f ',tempStr,r,p), 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('%s',tempStr), 'FontSize', 9, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
    end
end



%% Testing for meta-memory ~ memory precision + length + memory precision * length
if true
    meta_seqLevel_choice;
    seqPrecision_neuron_choice;
    
    temp_if_zscore = 1;
    
    temp_seqLength = sum(boolIndex_location_seq_T,2);
    temp_seqLength = temp_seqLength(1:41);
    
    %x = [seqPrecision_neuron_choice,temp_seqLength];
    %x = [seqPrecision_neuron_choice,1./temp_seqLength];
    
    x1 = seqPrecision_neuron_choice;
    %x2 = 1./temp_seqLength;
    x2 = temp_seqLength;
    
    x1_n11n = (x1-mean(x1,'omitnan'))/std(x1,'omitnan'); % z-score
    x2_n11n = (x2-mean(x2,'omitnan'))/std(x2,'omitnan'); % z-score
    
    if temp_if_zscore == 0
        x = [x1,x2];
    elseif temp_if_zscore == 1
        x = [x1_n11n,x2_n11n];
    end
    
    y = meta_seqLevel_choice;
    
    if temp_if_zscore == 1
        y = (y-mean(y,'omitnan'))/std(y,'omitnan'); % z-score
    end
    
    %temp_mdl= fitglm(x,y,'linear','Distribution','normal','Intercept',true);
    temp_mdl= fitglm(x,y,'interactions','Distribution','normal','Intercept',true);
    temp_glm_beta = temp_mdl.Coefficients.Estimate;
    temp_glm_r2 = temp_mdl.Rsquared.Adjusted;
    
end


%% End