%% Initialization
% clear
close all


% currentSession = '113Recording_20230426A_Ding_Site16';
% currentSession = '113Recording_20230427A_Ding_Site16_sameFOV0426';
% currentSession = '113Recording_20230502A_Ding_Site13';
% currentSession = '113Recording_20230503A_Ding_Site13_sameFOV0502';
% currentSession = '113Recording_20230504A_Ding_Site02';
% currentSession = '113Recording_20230509A_Ding_Site02';
% currentSession = '113Recording_20230508A_Ding_Site02_sameFOV0509';
% currentSession = '113Recording_20230511A_Ding_Site05';
% currentSession = '113Recording_20230510A_Ding_Site05_sameFOV0511';
% currentSession = '113Recording_20230510B_Ding_Site05_sameFOV0511';
% currentSession = '113Recording_20230512A_Ding_Site09';
% currentSession = '113Recording_20230513A_Ding_Site09_sameFOV0512';
% currentSession = '113Recording_20230515A_Ding_Site24_sameFOV0516';
% currentSession = '113Recording_20230516A_Ding_Site24';
% currentSession = '113Recording_20230517A_Ding_Site16B';
% currentSession = '113Recording_20230522A_Ding_Site05B';
% currentSession = '113Recording_20230525A_Ding_Site09B';
% currentSession = '113Recording_20230526A_Ding_Site09B_sameFOV0525';
% currentSession = '113Recording_20230527A_Ding_Site02B';
currentSession = '113Recording_20230528A_Ding_Site02B_sameFOV0527';
% currentSession = '113Recording_20230530A_Ding_Site05C';
% currentSession = '113Recording_20230531A_Ding_Site05C_sameFOV0530';
% currentSession = '113Recording_20230601A_Ding_Site13B';
% currentSession = '113Recording_20230602A_Ding_Site13B_sameFOV0601';
% currentSession = '113Recording_20230604A_Ding_Site07';
% currentSession = '113Recording_20230605A_Ding_Site07_sameFOV0604';
% currentSession = '113Recording_20230612A_Ding_Site14';
% currentSession = '113Recording_20230614A_Ding_Site05D';
% currentSession = '113Recording_20230615A_Ding_Site05D_sameFOV0614';
% currentSession = '113Recording_20230619A_Ding_Site02C';


targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)


output_shortPath = 'D:\twoPhotonData_motionCorrected';
temp_currentSession_path = [output_shortPath '\' currentSession];
temp_if_max0_min1 = 0;
output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);



a = 1;

%% Params setting
if_compute_train = 1;%1
if_compute_test = 1;%1

if_plot = 0;

if_seqLevel0_trialLevel1 = 1;
if_decoder_6location0_56seq1 = 0;%0


numFrames = 6;
rng(1); % For reproducibility
t_decoder = templateSVM('Standardize',true,'KernelFunction','linear'); % To standardise input data




%% Load data
trial_para = []; %#ok<*NASGU>

% load decodingData
% load trial_para
% load markerParse_trialLevel
if exist('decodingData','var') == 0
    load([output_path,'\','decodingData.mat'],'decodingData');
end
load([output_path,'\','trial_para.mat'],'trial_para');
load([output_path,'\','markerParse_trialLevel.mat'],'markerParse_trialLevel');


temp_fields = fieldnames(decodingData);
for tempi=1:length(temp_fields)
    eval([temp_fields{tempi},'=decodingData.(temp_fields{tempi});']);
end

searchName_gAcc = 'from23-04-26to23-06-20_D_gAcc_1';
searchName_rProb = 'from23-04-26to23-06-20_D_offloadingProb_1';

%path_behavior = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\behavior';
path_behavior = [output_shortPath '\behavior'];

% Load other processed results
load_gAcc = loadMat_single(searchName_gAcc, path_behavior);
gAcc_noChoice_collapsed_inOne = load_gAcc.gAcc_noChoice_collapsed_inOne;
gAcc_noChoice_collapsed = load_gAcc.gAcc_noChoice_collapsed;

% Load other processed results
load_rProb = loadMat_single(searchName_rProb, path_behavior);
offloadingProb_collapsed = load_rProb.offloadingProb_all;

offloadingProb_inOne = [];
pointKindsNum = 4;
for tempi=1:pointKindsNum
    offloadingProb_inOne = [offloadingProb_inOne offloadingProb_collapsed{tempi}']; %#ok<*AGROW>
end


cd(targetPATH)

%% Others
a = 1;


numFrames = 6;
pointKindsNum = 4;
target_seqSet = get_target_seqSet_v1(numFrames,pointKindsNum);
numSeq = zeros(1,pointKindsNum);
for tempi=1:pointKindsNum
    numSeq(tempi) = length(target_seqSet{tempi});
end

% F_dff_decisionBin1 = mean(F_dff_decisionPeriod(:,:,decisionPeriod_interval(1):decisionPeriod_interval(2)),3);
% F_dff_decisionBin1 = mean(F_dff_decisionPeriod(:,:,1:27),3);
% F_dff_decisionBin1 = mean(F_dff_decisionPeriod(:,:,decisionPeriod_interval(2):decisionPeriod_interval(3)),3);

F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);

temp_F_dff_decisionBin = F_dff_decisionBin1;

roiNum = size(temp_F_dff_decisionBin,1);

a = 1;
if if_seqLevel0_trialLevel1 == 0
    % abandoned now
elseif if_seqLevel0_trialLevel1 == 1
    
    %     temp_trialIndex_bool = trialIndex_lengthx_bool(2,:);
    %temp_trialIndex_bool = trialIndex_lengthx_bool_memoryCorrect(1,:);
    %     temp_trialIndex_bool = trialIndex_lengthx_bool_memoryCorrect(2,:);
    %temp_trialIndex_bool = trialIndex_lengthx_bool_memoryCorrect(1,:) | trialIndex_lengthx_bool_memoryCorrect(2,:) | trialIndex_lengthx_bool_memoryCorrect(3,:);
    temp_trialIndex_bool = trialIndex_bool_memoryCorrect; % use memory correct trial
    %     temp_trialIndex_bool = true(1,size(F_dff_decisionBin1,2)); % use all trial
    %temp_trialIndex_bool = trialIndex_bool_choiceMemory; % use choice memory trial
    
    temp_F_dff_decisionBin = F_dff_decisionBin1(:,temp_trialIndex_bool);
    
    seqIndex_valid = seqIndex(temp_trialIndex_bool);
    
    boolIndex_location_trial = false(numFrames,length(seqIndex_valid));
    for tempi=1:length(seqIndex_valid)
        currentSequence = target_seqSet_inOne{seqIndex_valid(tempi)};
        boolIndex_location_trial(currentSequence,tempi) = true;
    end
    
    boolIndex_location_seq = false(numFrames,length(target_seqSet_inOne));
    for tempi=1:length(target_seqSet_inOne)
        currentSequence = target_seqSet_inOne{tempi};
        boolIndex_location_seq(currentSequence,tempi) = true;
    end
    
end

%% Train decoder in length 1 trials
if if_compute_train == 1
    t_fun_seqProbSVM_6location = tic;
    if if_seqLevel0_trialLevel1 == 0
        % abandoned now
    elseif if_seqLevel0_trialLevel1 == 1
        if if_decoder_6location0_56seq1 == 0
            seqProbSVM_6location_trialLevel_train_Name_v = autoGetFunName_myScripts('seqProbSVM_6location_trialLevel_train', [targetPATH '\functions']);
            fun_seqProbSVM_6location_trialLevel_train = str2func(seqProbSVM_6location_trialLevel_train_Name_v);
            
            if_resample = 0;
            if if_resample == 0
                a = 1;
                %gAcc = [0.9,0.9,0.9,0.9,0.75,0.8];
                gAcc = gAcc_noChoice_collapsed_inOne(1:numFrames);
                
                svm_options = struct;
                svm_options.boolIndex_location_seq = boolIndex_location_seq;
                svm_options.t_decoder = t_decoder;
                svm_options.numFrames = numFrames;
                svm_options.targetPATH = targetPATH;
                svm_options.gAcc = gAcc;
                svm_options.seq_range = 1:numFrames;
                svm_options.target_length = 1;
                
                % use length 1 trials to train decoder
                temp_range = 1:numFrames;
                tempBoolIndex = ismember(seqIndex_valid,temp_range);
                x = temp_F_dff_decisionBin(:,tempBoolIndex);
                y = boolIndex_location_trial(:,tempBoolIndex);
                
                svm_train_outputs = fun_seqProbSVM_6location_trialLevel_train(x,y,svm_options);
                
                svm_posterior_length1 = svm_train_outputs.svm_posterior_lengthx;
                coeff_w_power = svm_train_outputs.coeff_w_power;
                temp_Mdl_CV_binary = svm_train_outputs.temp_Mdl_CV_binary;
                temp_Mdl_binary = svm_train_outputs.temp_Mdl_binary;
                
                
                a = 1;
            elseif if_resample == 1
                % abandoned now
            end
        elseif if_decoder_6location0_56seq1 == 1
            % abandoned now
        end
    end
    fprintf('t_fun_seqProbSVM_6location = %.1f secs.\n',toc(t_fun_seqProbSVM_6location));
elseif if_compute_train == 0
    
end
a = 1;

%% Test decoder in length 2 3 trials
if if_compute_test == 1
    a = 1;
    
    seqProbSVM_6location_trialLevel_test_Name_v = autoGetFunName_myScripts('seqProbSVM_6location_trialLevel_test', [targetPATH '\functions']);
    fun_seqProbSVM_6location_trialLevel_test = str2func(seqProbSVM_6location_trialLevel_test_Name_v);
    
    svm_train_outputs;
    svm_posterior_length1 = svm_train_outputs.svm_posterior_lengthx;
    coeff_w_power = svm_train_outputs.coeff_w_power;
    temp_Mdl_CV_binary = svm_train_outputs.temp_Mdl_CV_binary;
    temp_Mdl_binary = svm_train_outputs.temp_Mdl_binary;
    
    svm_test_options = struct;
    svm_test_options.boolIndex_location_seq = boolIndex_location_seq;
    svm_test_options.t_decoder = t_decoder;
    svm_test_options.numFrames = numFrames;
    svm_test_options.targetPATH = targetPATH;
    
    svm_test_options.coeff_w_power = coeff_w_power;
    svm_test_options.temp_Mdl_CV_binary = temp_Mdl_CV_binary;
    svm_test_options.temp_Mdl_binary = temp_Mdl_binary;
    
    a = 1;
    %% use length 2 3 trials to test decoder
    for target_length=2:3
        temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
        gAcc = gAcc_noChoice_collapsed_inOne(temp_range);
        
        tempBoolIndex = ismember(seqIndex_valid,temp_range);
        x = temp_F_dff_decisionBin(:,tempBoolIndex);
        y = boolIndex_location_trial(:,tempBoolIndex);
        
        svm_test_options.gAcc = gAcc;
        svm_test_options.seq_range = temp_range;
        svm_test_options.target_length = target_length;
        
        if target_length == 2
            svm_test_length2_outputs = fun_seqProbSVM_6location_trialLevel_test(x,y,svm_test_options);
            svm_posterior_length2 = svm_test_length2_outputs.svm_posterior_lengthx;
        elseif target_length == 3
            svm_test_length3_outputs = fun_seqProbSVM_6location_trialLevel_test(x,y,svm_test_options);
            svm_posterior_length3 = svm_test_length3_outputs.svm_posterior_lengthx;
        end
    end
    a = 1;
    svm_posterior_length1_mean = mean(svm_posterior_length1,'omitnan');
    svm_posterior_length2_mean = mean(svm_posterior_length2,'omitnan');
    svm_posterior_length3_mean = mean(svm_posterior_length3,'omitnan');
    
    gAcc_length1 = gAcc_noChoice_collapsed{1}';
    gAcc_length2 = gAcc_noChoice_collapsed{2}';
    gAcc_length3 = gAcc_noChoice_collapsed{3}';
    
    a = 1;
    
    [r_1,p_1] = corr(svm_posterior_length1',gAcc_length1');
    [r_2,p_2] = corr(svm_posterior_length2(~isnan(svm_posterior_length2))',gAcc_length2(~isnan(svm_posterior_length2))');
    [r_3,p_3] = corr(svm_posterior_length3(~isnan(svm_posterior_length3))',gAcc_length3(~isnan(svm_posterior_length3))');
    fprintf('r_1=%.4f, p_1=%.4f.\n',r_1,p_1);
    fprintf('r_2=%.4f, p_2=%.4f.\n',r_2,p_2);
    fprintf('r_3=%.4f, p_3=%.4f.\n',r_3,p_3);
    
    a = 1;
    
elseif if_compute_test == 0
    
end

a = 1;
svm_test_length2_outputs;
svm_test_length3_outputs;

for target_length=1:3
    if target_length == 2
        temp_svm = svm_test_length2_outputs;
    elseif target_length == 3
        temp_svm = svm_test_length3_outputs;
    elseif target_length == 1
        temp_svm = svm_train_outputs;
    end
    
    temp_trialNum = size(temp_svm.temp_svm_Y_valid_T,1);
    temp_seqIndex_trial = zeros(temp_trialNum,1);
    for tempi=1:temp_trialNum
        temp_seqBoolIndex_trial = temp_svm.temp_svm_Y_valid_T(tempi,:);
        for tempj=1:numSeq(target_length)
            temp_seqIndex = temp_svm.seq_range(tempj);
            temp_seqBoolIndex = temp_svm.boolIndex_location_seq(:,temp_seqIndex)';
            
            if sum(temp_seqBoolIndex==temp_seqBoolIndex_trial) == numFrames
                temp_seqIndex_trial(tempi) = temp_svm.seq_range(tempj);
                break
            end
        end
    end
    temp_svm.seqIndex_trial = temp_seqIndex_trial;
    
    a = 1;
    Posterior_2d_n11n_lengthx_mean = zeros(numSeq(target_length),numFrames);
    Posterior_2d_w_lengthx_mean = zeros(numSeq(target_length),numFrames);
    for tempi=1:numSeq(target_length)
        temp_seqIndex = temp_svm.seq_range(tempi);
        tempBoolIndex = temp_seqIndex_trial == temp_seqIndex;
        temp_Posterior_2d_n11n = temp_svm.Posterior_2d_n11n(tempBoolIndex,:);
        Posterior_2d_n11n_lengthx_mean(tempi,:) = mean(temp_Posterior_2d_n11n,1);
        
        temp_Posterior_2d_w = temp_svm.Posterior_2d_w(tempBoolIndex,:);
        Posterior_2d_w_lengthx_mean(tempi,:) = mean(temp_Posterior_2d_w,1);
    end
    temp_svm.Posterior_2d_n11n_lengthx_mean = Posterior_2d_n11n_lengthx_mean;
    temp_svm.Posterior_2d_w_lengthx_mean = Posterior_2d_w_lengthx_mean;
    
    if target_length == 2
        svm_test_length2_outputs = temp_svm;
    elseif target_length == 3
        svm_test_length3_outputs = temp_svm;
    elseif target_length == 1
        svm_train_outputs = temp_svm;
    end
end

a = 1;

svm_train_outputs;




Posterior_2d_n11n_mean = [svm_train_outputs.Posterior_2d_n11n_lengthx_mean;...
    svm_test_length2_outputs.Posterior_2d_n11n_lengthx_mean;...
    svm_test_length3_outputs.Posterior_2d_n11n_lengthx_mean];

Posterior_2d_w_mean = [svm_train_outputs.Posterior_2d_w_lengthx_mean;...
    svm_test_length2_outputs.Posterior_2d_w_lengthx_mean;...
    svm_test_length3_outputs.Posterior_2d_w_lengthx_mean];


s = load(['C:\ASDROOT\STUDY\temp\data\','FittingResults_CCM_noOrder_distri_kwp0_v11_20230831A']);
% s = load(['C:\ASDROOT\STUDY\temp\data\','FittingResults_CCM_noOrder_distri_kwp0_lp0_v3_20230831A']);
FittingResults = s.FittingResults;

a = 1;

valid_length = 3;
temp_range = 1:sum(numSeq(1:valid_length));
Posterior_2d_model = FittingResults.intermediate.locationDistri_seq(temp_range,:);

a = 1;

Posterior_2d_n11n_mean;
Posterior_2d_w_mean;
Posterior_2d_model;

r_n11n = zeros(1,length(temp_range));
p_n11n = zeros(1,length(temp_range));
for tempi=1:length(temp_range)
    [r_n11n(tempi),p_n11n(tempi)] = corr(Posterior_2d_n11n_mean(tempi,:)',Posterior_2d_model(tempi,:)');
end


r_w = zeros(1,length(temp_range));
p_w = zeros(1,length(temp_range));
for tempi=1:length(temp_range)
    [r_w(tempi),p_w(tempi)] = corr(Posterior_2d_w_mean(tempi,:)',Posterior_2d_model(tempi,:)');
end

fprintf('num(p_n11n<0.05)=%d, sum(p_w<0.05)=%d.\n',sum(p_n11n<0.05),sum(p_w<0.05));

a = 1;

r_loc_n11n = zeros(1,numFrames);
p_loc_n11n = zeros(1,numFrames);
for tempi=1:numFrames
    tempBoolIndex = ~isnan(Posterior_2d_n11n_mean(:,1));
    [r_loc_n11n(tempi),p_loc_n11n(tempi)] = ...
        corr(Posterior_2d_n11n_mean(tempBoolIndex,tempi),Posterior_2d_model(tempBoolIndex,tempi));
end

r_loc_w = zeros(1,numFrames);
p_loc_w = zeros(1,numFrames);
for tempi=1:numFrames
    tempBoolIndex = ~isnan(Posterior_2d_w_mean(:,1));
    [r_loc_w(tempi),p_loc_w(tempi)] = ...
        corr(Posterior_2d_w_mean(tempBoolIndex,tempi),Posterior_2d_model(tempBoolIndex,tempi));
end

a = 1;
r_loc_n11n;
r_loc_w;
fprintf('r_loc_n11n= %s.\n',num2str(round(r_loc_n11n*10000)/10000));
fprintf('r_loc_w=    %s.\n',num2str(round(r_loc_w*10000)/10000));

locMSE_n11n = mean((Posterior_2d_n11n_mean-Posterior_2d_model).^2,1,'omitnan');
locMSE_w = mean((Posterior_2d_w_mean-Posterior_2d_model).^2,1,'omitnan');

locMSE_n11n;
a = 1;
% fprintf('locMSE_n11n=%s.\n',num2str(round(locMSE_n11n*10000)/10000));
% fprintf('locMSE_w=   %s.\n',num2str(round(locMSE_w*10000)/10000));

a = 1;


%% Fig
if if_plot == 1
    close all
    for tempi=length(temp_range):-1:1
        % for tempi=6:-1:1
        
        if sum(isnan(Posterior_2d_n11n_mean(tempi,:)))>0
            continue
        end
        
        fig = figure('Name','','NumberTitle','off');
        set(gcf,'Position',[35+400 35+150 520 320]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        templd1 = Posterior_2d_model(tempi,:);
        templd2 = Posterior_2d_n11n_mean(tempi,:);
        templd3 = Posterior_2d_w_mean(tempi,:);
        
        plot(templd1,'color',[0 0 0],'LineWidth',1.5);
        hold on
        plot(templd2,'color',[0.75 0.25 0.25],'LineWidth',1.5);
        hold on
        plot(templd3,'color',[0.25 0.25 0.75],'LineWidth',1.5);
        hold on
        
        hl = legend('model','n11n','w',...
            'Location','northeast','fontsize',9);
        
        ylim([0 1]);
        set(gca, 'FontSize', 18)
        set(gca,'XLim',[1-0.5 numFrames+0.5]);%X轴的数据显示范围
        set(gca,'XTick',[1:1:numFrames]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        set(gca,'XTickLabel',[1:1:numFrames]);%给坐标加标签
        set(gca,'box','off');% 取消右、上边框
        xlabel('Location', 'FontSize', 16, 'FontWeight', 'bold');
        ylabel('probability', 'FontSize', 16, 'FontWeight', 'bold');
        
        
        temp_title = title(sprintf('%s',num2str(boolIndex_location_seq(:,tempi)')));
        temp_title.FontSize = 16;
        temp_title.Interpreter = 'none';
    end
end
%% End