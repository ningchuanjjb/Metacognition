%% Initialization
% clear
close all


%% Params setting
if_compute_train = 0;%1

if_dff_singleSession1_twoSession2_allSession3 = 1;

currentSessionIndex_AB = 3;

currentSessionIndex = 20;%20-->10-->28-->5-->16-->21-->19

if_plot = 1;



KFold_num = 5;%10-->5
rng(5);%5
t_decoder = templateSVM('Standardize',true,'KernelFunction','linear'); % To standardise input data



targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)
output_shortPath = 'D:\twoPhotonData_motionCorrected';

if if_dff_singleSession1_twoSession2_allSession3 == 1
    currentSession_multi = string;
    
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
    
    
    %currentABSession_multi = [currentABSession_multi; '20230508A_and_20230509A'];%3 few
    
    currentABSession_multi(1) = [];
    num_FOV_AB = length(currentABSession_multi);
    
    currentSession = currentABSession_multi{currentSessionIndex_AB};
    
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

if exist('decodingDataSimplified','var') == 0 || if_compute_train == 1
    if if_dff_singleSession1_twoSession2_allSession3 == 1
        temp_load = load([output_path,'\','decodingDataSimplified.mat'],'decodingDataSimplified');
        decodingDataSimplified = temp_load.decodingDataSimplified;
    elseif if_dff_singleSession1_twoSession2_allSession3 == 2
        temp_load = load([output_path '\' currentSession]);
        decodingDataSimplified = temp_load.decodingDataSimplified_AB;
    elseif if_dff_singleSession1_twoSession2_allSession3 == 3
        
    end
end

temp_fields = fieldnames(decodingDataSimplified);
for tempi=1:length(temp_fields)
    eval([temp_fields{tempi},'=decodingDataSimplified.(temp_fields{tempi});']);
end

searchName_rProb = 'from23-04-26to23-06-20_D_offloadingProb_1';

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

% load trial_para_choiceCondition_flag, it's ok for single session
temp_load = load([output_path,'\','trial_para.mat']);
trial_para = temp_load.trial_para;
trial_para_choiceCondition_flag = trial_para.choiceCondition_flag;

% trial_para_choiceCondition_flag;

temp_trial_para_ifSelectOffloading = trial_para_ifSelectOffloading==1;

choiceBoolIndex = trial_para_choiceCondition_flag == 2;
allMemoryBoolIndex = ~temp_trial_para_ifSelectOffloading;
choiceMemoryBoolIndex = choiceBoolIndex & allMemoryBoolIndex;
choiceOffloadBoolIndex = choiceBoolIndex & temp_trial_para_ifSelectOffloading;

allMemoryCorrectBoolIndex = trialIndex_bool_memoryCorrect;
allMemoryErrorBoolIndex = allMemoryBoolIndex & ~trialIndex_bool_memoryCorrect;

choiceMemoryCorrectBoolIndex = choiceMemoryBoolIndex & trialIndex_bool_memoryCorrect;
choiceMemoryErrorBoolIndex = choiceMemoryBoolIndex & ~trialIndex_bool_memoryCorrect;

seqIndex_choice = seqIndex(choiceBoolIndex);

%% More preparation
if if_dff_singleSession1_twoSession2_allSession3 < 3
    F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
    % F_dff_decisionBin1 = mean(F_dff_decisionPeriod(:,:,decisionPeriod_interval(2):decisionPeriod_interval(3)),3);
    F_dff_decisionBin1 = double(F_dff_decisionBin1);
    F_dff_decisionBin1 = F_dff_decisionBin1 + eps;
    
    roiNum = size(F_dff_decisionBin1,1);
    
elseif if_dff_singleSession1_twoSession2_allSession3 == 3
    
end
a = 1;

%% Compute svm
if if_compute_train == 1
    
    fun_metaDecoder_Name_v = autoGetFunName_myScripts('metaDecoder', [targetPATH '\functions']);
    fun_metaDecoder = str2func(fun_metaDecoder_Name_v);
    
    svm_options = struct;
    svm_options.numSeq = numSeq;
    svm_options.seqIndex_choice = seqIndex_choice;
    svm_options.KFold_num = KFold_num;
    svm_options.t_decoder = t_decoder;
    
    x_raw = F_dff_decisionBin1';
    x = x_raw(choiceBoolIndex,:);
    y_raw = choiceMemoryBoolIndex';
    y = y_raw(choiceBoolIndex);
    
    
    svm_Meta = fun_metaDecoder(x,y,svm_options);
    
end
svm_Meta;
fprintf('Meta-memory decoding accuracy (all trial) = %.3f.\n',svm_Meta.svm_posterior_allTrial);


fprintf('t_script = %.1f secs.\n',toc(t0));

%%  Compute meta-memory on training set trials (choice trials)
meta_trialLevel_choice = svm_Meta.Posterior_2d_n11n;


%%  Compute meta-memory on other trials (noChoice trials)
if if_compute_train == 1
    x_raw = F_dff_decisionBin1';
    x = x_raw(~choiceBoolIndex,:);
    
    temp_Mdl_CV_binary = svm_Meta.temp_Mdl_CV_binary;
    temptemp_Posterior_2d_kfold = zeros(KFold_num,size(x,1));
    for tempk=1:KFold_num
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
end


%% Compute meta-memory on all trials
meta_trialLevel = zeros(length(seqIndex),1);
meta_trialLevel(choiceBoolIndex) = meta_trialLevel_choice;
meta_trialLevel(~choiceBoolIndex) = meta_trialLevel_noChoice;

%% Compute meta_seqLevel
meta_seqLevel_choice = zeros(sum(numSeq(1:3)),1);
meta_seqLevel_choiceMemory = zeros(sum(numSeq(1:3)),1);
meta_seqLevel_choiceOffload = zeros(sum(numSeq(1:3)),1);
meta_seqLevel_noChoice = zeros(sum(numSeq(1:3)),1);

for tempi=1:sum(numSeq(1:3))
    temptempBoolIndex = seqIndex==tempi;
    
    temptempBoolIndex_choice = temptempBoolIndex & choiceBoolIndex;
    temptempBoolIndex_choiceMemory = temptempBoolIndex & choiceMemoryBoolIndex;
    temptempBoolIndex_choiceOffload = temptempBoolIndex & choiceOffloadBoolIndex;
    temptempBoolIndex_noChoice = temptempBoolIndex & (~choiceBoolIndex);
    
    meta_seqLevel_choice(tempi) = mean(meta_trialLevel(temptempBoolIndex_choice));
    meta_seqLevel_choiceMemory(tempi) = mean(meta_trialLevel(temptempBoolIndex_choiceMemory));
    meta_seqLevel_choiceOffload(tempi) = mean(meta_trialLevel(temptempBoolIndex_choiceOffload));
    meta_seqLevel_noChoice(tempi) = mean(meta_trialLevel(temptempBoolIndex_noChoice));
end
meta_seqLevel_allCondition = [meta_seqLevel_choice,meta_seqLevel_choiceMemory,meta_seqLevel_choiceOffload,meta_seqLevel_noChoice];

temp1 = isnan(meta_seqLevel_allCondition);
temp2 = sum(temp1,2)==0;

% meta_seqLevelMean_allCondition = mean(meta_seqLevel_allCondition,1,'omitnan');
meta_seqLevelMean_allCondition = mean(meta_seqLevel_allCondition(temp2,:),1,'omitnan');



%% Plot
if if_plot == 1        
    fig = figure('Name','meta_seqLevel','NumberTitle','off');
    set(gcf,'Position',[50+0 50+0 1350 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    
    t.Title.String = sprintf('%d seqs (length123), %s',...
        length(y),currentSession_short);
    t.Title.FontSize = 16;
    t.Title.Interpreter = 'none';
    
    
    %% fig, decoding accuracy for meta-memory
    nexttile
    
    y = svm_Meta.svm_posterior_seqLevel';
    
    h = [];
    for tempi=1:3
        temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
        tempIndex = find(tempBoolIndex(temp_range)==true);
        temp_range2 = temp_range(tempIndex);
        
        temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
        temp_h = scatter(temp_range2,y(temp_range2), ...
            temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        h = [h temp_h]; %#ok<*AGROW>
    end
    
    
    xlim([0 sum(numSeq(1:3))+1]);
    ylim([0 1]);
    
    set(gca, 'FontSize', 11)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('Sequence'), 'FontSize', 14, 'FontWeight', 'bold');
    
    temp_title = title(sprintf('Decoding accuracy = %.3f',svm_Meta.svm_posterior_allTrial_length123), 'FontSize', 14, 'FontWeight', 'bold');
    temp_title.Interpreter = 'none';
    
    ylabel(sprintf('Decoding accuracy'), 'FontSize', 14, 'FontWeight', 'bold');
    
    
    %% fig, meta_seqLevel
    for loopCount=1:2
        nexttile
        if loopCount == 1
            tempStr = 'Free choice (training set)';
            x = meta_seqLevel_choice;
        elseif loopCount == 2
            tempStr = 'Force to test (testing set)';
            x = meta_seqLevel_noChoice;
        end
        y = offloadingProb_inOne(1:sum(numSeq(1:3)))';
        
        [r,p] = corr(x,y);
        
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
        
        xlim([0 1]);
        ylim([0 1]);
        
        set(gca, 'FontSize', 11)
        set(gca,'box','off');% 取消右、上边框
        xlabel(sprintf('Meta-memory'), 'FontSize', 14, 'FontWeight', 'bold');
        
        temp_title = title(sprintf('%s, r=%.3f, p=%.3f',tempStr,r,p), 'FontSize', 14, 'FontWeight', 'bold');
        temp_title.Interpreter = 'none';
        
        ylabel(sprintf('Offloading rate'), 'FontSize', 14, 'FontWeight', 'bold');
        
    end
end

%% End