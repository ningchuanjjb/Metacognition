% Chuan's 28th script (20251214)
% This script: Cross-time location memory decoding.
%% Initialization
% clear
close all



% currentSession = '113Recording_20240111A_Zelku_Site09A';
% currentSession = '113Recording_20240112A_Zelku_Site06A';
% currentSession = '113Recording_20240115A_Zelku_Site06A';
% currentSession = '113Recording_20240117A_Zelku_Site14A';
% currentSession = '113Recording_20240118A_Zelku_Site18A';
% currentSession = '113Recording_20240119A_Zelku_Site17A';
% currentSession = '113Recording_20240122A_Zelku_Site09B';
% currentSession = '113Recording_20240123A_Zelku_Site09B_sameFOV0122';
% currentSession = '113Recording_20240124A_Zelku_Site06B';
% currentSession = '113Recording_20240126A_Zelku_Site06B_sameFOV0124';
% currentSession = '113Recording_20240129A_Zelku_Site07A';
% currentSession = '113Recording_20240131A_Zelku_Site07A_sameFOV0129';
% currentSession = '113Recording_20240202A_Zelku_Site06XA';
% currentSession = '113Recording_20240203A_Zelku_Site06XA_sameFOV0202';
% currentSession = '113Recording_20240207A_Zelku_Site05A';
% currentSession = '113Recording_20240208A_Zelku_Site05A_sameFOV0207';
% currentSession = '113Recording_20240210A_Zelku_Site10A';
% currentSession = '113Recording_20240211A_Zelku_Site10A_sameFOV0210';
% currentSession = '113Recording_20240216A_Zelku_Site09C';
% currentSession = '113Recording_20240218A_Zelku_Site09C_sameFOV0216';
% currentSession = '113Recording_20240220A_Zelku_Site06XB';
% currentSession = '113Recording_20240221A_Zelku_Site06XB_sameFOV0220';
% currentSession = '113Recording_20240226A_Zelku_Site10B';
% currentSession = '113Recording_20240227A_Zelku_Site10B_sameFOV0226';
% currentSession = '113Recording_20240229A_Zelku_Site06C';
% currentSession = '113Recording_20240301A_Zelku_Site06C_sameFOV0229';
% currentSession = '113Recording_20240304A_Zelku_Site09D';
% currentSession = '113Recording_20240305A_Zelku_Site09D_sameFOV0304';
% currentSession = '113Recording_20240307A_Zelku_Site10C';
% currentSession = '113Recording_20240308A_Zelku_Site10C_sameFOV0307';
% currentSession = '113Recording_20240312A_Zelku_Site06RA';
% currentSession = '113Recording_20240315A_Zelku_Site06RA_sameFOV0312';
% currentSession = '113Recording_20240319A_Zelku_Site09E';
% currentSession = '113Recording_20240320A_Zelku_Site09E_sameFOV0319';
% currentSession = '113Recording_20240322A_Zelku_Site07B';
% currentSession = '113Recording_20240323A_Zelku_Site07B_sameFOV0322';
% currentSession = '113Recording_20240329A_Zelku_Site05B';
% currentSession = '113Recording_20240330A_Zelku_Site05B_sameFOV0329';
% currentSession = '113Recording_20240402A_Zelku_Site14B';
% currentSession = '113Recording_20240403A_Zelku_Site14B_sameFOV0402';
% currentSession = '113Recording_20240410A_Zelku_Site17B';
% currentSession = '113Recording_20240411A_Zelku_Site17B_sameFOV0410';


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
% currentSession = '113Recording_20230620A_Ding_Site05E';


targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)


% output_shortPath = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected';
output_shortPath = 'D:\twoPhotonData_motionCorrected';
temp_currentSession_path = [output_shortPath '\' currentSession];
temp_if_max0_min1 = 0;
output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);


%% Params setting
windowSize = 1;%3-->5
if_compute = 1;%1
lengthFlag = 3;
if_plot = 1;
criteria_all1_correct2_new3 = 2;%2
% criteria_allmemory0_nochoicememory1_choicememory2_choiceoffloading3 = 1;
criteria2_trialType = 0; % allmemory0 nochoicememory1 choicememory2 choiceoffloading3

numFrames = 6;
rng(1); % For reproducibility
t_SVM = templateSVM('Standardize',true); % To standardise input data
filter_b = (1/windowSize)*ones(1,windowSize);% param for filter
filter_a = 1;% param for filter

%% Load data
trial_para = [];

% load decodingData
% load trial_para
% load markerParse_trialLevel
load([output_path,'\','decodingData.mat'],'decodingData');
load([output_path,'\','trial_para.mat'],'trial_para');
load([output_path,'\','markerParse_trialLevel.mat'],'markerParse_trialLevel');

temp_fields = fieldnames(decodingData);
for tempi=1:length(temp_fields)
    eval([temp_fields{tempi},'=decodingData.(temp_fields{tempi});']);
end


%% Preparation for length1
% ifvalidTrial_length1 = false(1, length(trialIndex_length1));
ifvalidTrialType_length1 = false(1, length(trialIndex_length1));
ifCorrect_length1 = false(1, length(trialIndex_length1));
responseBool_length1 = false(length(trialIndex_length1) ,numFrames);

ifSelectOffloading = trial_para.ifSelectOffloading(1:trial_para.trial_count);
ifSelectOffloading_bool = ifSelectOffloading==1;

for tempi=1:length(trialIndex_length1)
    trialIndex_marker = trialIndex_length1(tempi);
    trialIndex_PTB = markerParse_trialLevel{trialIndex_marker}.trial_index;
    responseBool_length1(tempi,:) = ~trial_para.isFilled{trialIndex_PTB};
    
    if criteria_all1_correct2_new3 == 1
        ifCorrect_length1(tempi) = true;
    elseif criteria_all1_correct2_new3 == 2
        ifCorrect_length1(tempi) = trial_para.isCorrect(trialIndex_PTB)==1;
    elseif criteria_all1_correct2_new3 == 3
        ifCorrect_length1(tempi) = trial_para.selectFlag_newSeq(trialIndex_PTB)==1;
    end
    %ifCorrect_length1(tempi) = (trial_para.isCorrect(trialIndex_PTB)==1) & (trial_para.selectFlag_newSeq(trialIndex_PTB)==1);
    
    a = 1;
    if criteria2_trialType == 0
        %ifvalidTrialType_length1(tempi) = trial_para.choiceCondition_flag(trialIndex_PTB)==0;
        ifvalidTrialType_length1(tempi) = ifSelectOffloading_bool(trialIndex_PTB)==0;
    end
    
end
ifvalidTrial_length1 = ifCorrect_length1 & ifvalidTrialType_length1;

sequence_length1 = sequence_length1(ifvalidTrial_length1==1);
sequenceBool_length1 = false(length(sequence_length1), numFrames);
for tempi=1:length(sequence_length1)
    sequenceBool_length1(tempi,sequence_length1(:,tempi)') = true;
end
responseBool_length1 = responseBool_length1(ifvalidTrial_length1==1,:);

F_dff_length1_period1_raw = F_dff_length1_period1;
F_dff_length1_period1 = F_dff_length1_period1(:,ifvalidTrial_length1==1,:);

% F_dff_length1_period1 = filter(filter_b,filter_a,F_dff_length1_period1,[],3);
F_dff_length1_period1 = smoothdata(F_dff_length1_period1,3,'gaussian',windowSize);

F_dff_length1_period1 = permute(F_dff_length1_period1,[2 1 3]);
T_length1_period1 = size(F_dff_length1_period1,3);

F_dff_length1_period2_raw = F_dff_length1_period2;
F_dff_length1_period2 = F_dff_length1_period2(:,ifvalidTrial_length1==1,:);

% F_dff_length1_period2 = filter(filter_b,filter_a,F_dff_length1_period2,[],3);
F_dff_length1_period2 = smoothdata(F_dff_length1_period2,3,'gaussian',windowSize);

F_dff_length1_period2 = permute(F_dff_length1_period2,[2 1 3]);
T_length1_period2 = size(F_dff_length1_period2,3);

%% Preparation for length2
ifCorrect_length2 = false(1, length(trialIndex_length2));
responseBool_length2 = false(length(trialIndex_length2) ,numFrames);
for tempi=1:length(trialIndex_length2)
    trialIndex_marker = trialIndex_length2(tempi);
    trialIndex_PTB = markerParse_trialLevel{trialIndex_marker}.trial_index;
    responseBool_length2(tempi,:) = ~trial_para.isFilled{trialIndex_PTB};
    
    if criteria_all1_correct2_new3 == 1
        ifCorrect_length2(tempi) = true;
    elseif criteria_all1_correct2_new3 == 2
        ifCorrect_length2(tempi) = trial_para.isCorrect(trialIndex_PTB)==1;
    elseif criteria_all1_correct2_new3 == 3
        ifCorrect_length2(tempi) = trial_para.selectFlag_newSeq(trialIndex_PTB)==1;
    end
    %ifCorrect_length2(tempi) = (trial_para.isCorrect(trialIndex_PTB)==1) & (trial_para.selectFlag_newSeq(trialIndex_PTB)==1);
end
sequence_length2 = sequence_length2(:,ifCorrect_length2==1);
sequenceBool_length2 = false(length(sequence_length2), numFrames);
for tempi=1:length(sequence_length2)
    sequenceBool_length2(tempi,sequence_length2(:,tempi)') = true;
end
responseBool_length2 = responseBool_length2(ifCorrect_length2==1,:);

F_dff_length2_period1_raw = F_dff_length2_period1;
F_dff_length2_period1 = F_dff_length2_period1(:,ifCorrect_length2==1,:);

% F_dff_length2_period1 = filter(filter_b,filter_a,F_dff_length2_period1,[],3);
F_dff_length2_period1 = smoothdata(F_dff_length2_period1,3,'gaussian',windowSize);

F_dff_length2_period1 = permute(F_dff_length2_period1,[2 1 3]);
T_length2_period1 = size(F_dff_length2_period1,3);

F_dff_length2_period2_raw = F_dff_length2_period2;
F_dff_length2_period2 = F_dff_length2_period2(:,ifCorrect_length2==1,:);

% F_dff_length2_period2 = filter(filter_b,filter_a,F_dff_length2_period2,[],3);
F_dff_length2_period2 = smoothdata(F_dff_length2_period2,3,'gaussian',windowSize);

F_dff_length2_period2 = permute(F_dff_length2_period2,[2 1 3]);
T_length2_period2 = size(F_dff_length2_period2,3);


%% Preparation for length3
ifCorrect_length3 = false(1, length(trialIndex_length3));
responseBool_length3 = false(length(trialIndex_length3) ,numFrames);
for tempi=1:length(trialIndex_length3)
    trialIndex_marker = trialIndex_length3(tempi);
    trialIndex_PTB = markerParse_trialLevel{trialIndex_marker}.trial_index;
    responseBool_length3(tempi,:) = ~trial_para.isFilled{trialIndex_PTB};
    
    if criteria_all1_correct2_new3 == 1
        ifCorrect_length3(tempi) = true;
    elseif criteria_all1_correct2_new3 == 2
        ifCorrect_length3(tempi) = trial_para.isCorrect(trialIndex_PTB)==1;
    elseif criteria_all1_correct2_new3 == 3
        ifCorrect_length3(tempi) = trial_para.selectFlag_newSeq(trialIndex_PTB)==1;
    end
    %ifCorrect_length3(tempi) = (trial_para.isCorrect(trialIndex_PTB)==1) & (trial_para.selectFlag_newSeq(trialIndex_PTB)==1);
end
sequence_length3 = sequence_length3(:,ifCorrect_length3==1);
sequenceBool_length3 = false(length(sequence_length3), numFrames);
for tempi=1:length(sequence_length3)
    sequenceBool_length3(tempi,sequence_length3(:,tempi)') = true;
end
responseBool_length3 = responseBool_length3(ifCorrect_length3==1,:);

F_dff_length3_period1_raw = F_dff_length3_period1;
F_dff_length3_period1 = F_dff_length3_period1(:,ifCorrect_length3==1,:);

% F_dff_length3_period1 = filter(filter_b,filter_a,F_dff_length3_period1,[],3);
F_dff_length3_period1 = smoothdata(F_dff_length3_period1,3,'gaussian',windowSize);

F_dff_length3_period1 = permute(F_dff_length3_period1,[2 1 3]);
T_length3_period1 = size(F_dff_length3_period1,3);

F_dff_length3_period2_raw = F_dff_length3_period2;
F_dff_length3_period2 = F_dff_length3_period2(:,ifCorrect_length3==1,:);

% F_dff_length3_period2 = filter(filter_b,filter_a,F_dff_length3_period2,[],3);
F_dff_length3_period2 = smoothdata(F_dff_length3_period2,3,'gaussian',windowSize);

F_dff_length3_period2 = permute(F_dff_length3_period2,[2 1 3]);
T_length3_period2 = size(F_dff_length3_period2,3);

%% Compute
svm_CCGP_eachLocName = 'svm_CCGP_eachLoc';
svm_CCGP_eachLocName_v = autoGetFunName_myScripts(svm_CCGP_eachLocName, [targetPATH '\functions']);
fprintf('Now runing %s.  ------> ', svm_CCGP_eachLocName_v);
fun_svm_CCGP_eachLoc = str2func(svm_CCGP_eachLocName_v);

a = 1;

if lengthFlag == 1
    temp_svm_Y = sequence_length1;
    responseBool_lengthx = responseBool_length1;
    sequenceBool_lengthx = sequenceBool_length1;
    
    F_dff_lengthx_period1 = F_dff_length1_period1;
    T_lengthx_period1 = T_length1_period1;
    lengthx_period1_interval = length1_period1_interval;
    
    F_dff_lengthx_period2 = F_dff_length1_period2;
    T_lengthx_period2 = T_length1_period2;
    
elseif lengthFlag == 2
    temp_svm_Y = sequence_length2;
    responseBool_lengthx = responseBool_length2;
    sequenceBool_lengthx = sequenceBool_length2;
    
    F_dff_lengthx_period1 = F_dff_length2_period1;
    T_lengthx_period1 = T_length2_period1;
    lengthx_period1_interval = length2_period1_interval;
    
    F_dff_lengthx_period2 = F_dff_length2_period2;
    T_lengthx_period2 = T_length2_period2;
    
elseif lengthFlag == 3
    temp_svm_Y = sequence_length3;
    responseBool_lengthx = responseBool_length3;
    sequenceBool_lengthx = sequenceBool_length3;
    
    F_dff_lengthx_period1 = F_dff_length3_period1;
    T_lengthx_period1 = T_length3_period1;
    lengthx_period1_interval = length3_period1_interval;
    
    F_dff_lengthx_period2 = F_dff_length3_period2;
    T_lengthx_period2 = T_length3_period2;
    
end
if if_compute == 1
    F_dff_lengthx_period12 = cat(3,F_dff_lengthx_period1,F_dff_lengthx_period2);
    T_lengthx_period12 = T_lengthx_period1 + T_lengthx_period2;
    lengthx_period12_interval = [lengthx_period1_interval lengthx_period2_interval+T_lengthx_period1];
    
    %     F_dff_lengthx_period12 = F_dff_lengthx_period12(:,1:200,:);
    %     F_dff_lengthx_period12 = F_dff_lengthx_period12(:,1:100,:);
    
    %     t0 = tic;
    %     [svm_accuracy_CV_lengthx_period1, crossTime_svm_accuracy_lengthx_period1] = ...
    %         fun_svm_CCGP_eachLoc(F_dff_lengthx_period1,temp_svm_Y,t_SVM,T_lengthx_period1,numFrames,responseBool_lengthx,sequenceBool_lengthx);
    %     fprintf('Decoding period1 over, time=%.2f secs.\n',toc(t0));
    %     fprintf('max svm_accuracy_CV_lengthx_period1 = %.3f.\n',max(svm_accuracy_CV_lengthx_period1));
    %     [min_acc_period1,max_acc_period1] = bounds(crossTime_svm_accuracy_lengthx_period1,'all');
    
    t1 = tic;
    
    
    temp_dff = F_dff_lengthx_period12;
    
    %     if(~exist('selectiveCellBoolIndex','var'))
    %         temp_dff = F_dff_lengthx_period12;
    %     else
    %         temp_dff = F_dff_lengthx_period12(:,selectiveCellBoolIndex,:);
    %         fprintf('Only use selective cell to decode. selectiveCellNum = %d.\n',sum(selectiveCellBoolIndex));
    %     end
    
    [svm_accuracy_CV_lengthx_period12, crossTime_svm_accuracy_lengthx_period12] = ...
        fun_svm_CCGP_eachLoc(temp_dff,temp_svm_Y,t_SVM,T_lengthx_period12,numFrames,responseBool_lengthx,sequenceBool_lengthx);
    
    %[svm_accuracy_CV_lengthx_period12, crossTime_svm_accuracy_lengthx_period12] = ...
    %   fun_svm_CCGP_eachLoc(F_dff_lengthx_period1,temp_svm_Y,t_SVM,T_lengthx_period1,numFrames,responseBool_lengthx,sequenceBool_lengthx);
    
    %[svm_accuracy_CV_lengthx_period12, crossTime_svm_accuracy_lengthx_period12] = ...
    %   svm_CCGP_v4(F_dff_lengthx_period12,temp_svm_Y,t_SVM,T_lengthx_period12,numFrames);
    
    fprintf('Decoding period12 over, time=%.2f secs.\n',toc(t1));
    fprintf('max svm_accuracy_CV_lengthx_period12 = %.3f.\n',max(svm_accuracy_CV_lengthx_period12));
    [min_acc_period12,max_acc_period12] = bounds(crossTime_svm_accuracy_lengthx_period12,'all');
end

%% Plot
if if_plot == 1
    %% Plot fig1
    fig1 = figure('Name','Fig1','NumberTitle','off');
    set(gcf,'Position',[35+198 35+325 1300 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
    
    nexttile
    %plot(svm_accuracy_CV_lengthx_period1,'-','LineWidth',1.5,'Color',[0 0.4470 0.7410]);
    plot(svm_accuracy_CV_lengthx_period12(1:T_lengthx_period1),'-','LineWidth',1.5,'Color',[0 0.4470 0.7410]);
    hold on
    for tempi=1:length(lengthx_period1_interval)
        plot(lengthx_period1_interval(tempi)*[1 1],[0 1],...
            '-','LineWidth',1,'Color',[0 0 0]);
        hold on
    end
    plot([1 T_lengthx_period1],[1 1]/nchoosek(numFrames,lengthFlag),'--','LineWidth',1,'Color',[0.4 0.4 0.4]);
    hold on
    
    set(gca,'xticklabel',[])
    set(gca,'xtick',[])
    xticks(lengthx_period1_interval);
    %     if lengthFlag == 1
    %         xticklabels({'Fixation','T1','Delay1','Delay1 end'});
    %     elseif lengthFlag == 2
    %         xticklabels({'Fixation','T1','T2','Delay1','Delay1 end'});
    %     elseif lengthFlag == 3
    %         xticklabels({'Fixation','T1','T2','T3','Delay1','Delay1 end'});
    %     end
    if lengthFlag == 1
        %xticklabels({'PreFixation','Fixation','T1','Delay1','LayoutOff',''});
        xticklabels({'PreFixation','Fixation','T1','Delay1',''});
    elseif lengthFlag == 2
        xticklabels({'PreFixation','Fixation','T1','T2','Delay1',''});
    elseif lengthFlag == 3
        xticklabels({'PreFixation','Fixation','T1','T2','T3','Delay1',''});
    end
    ylim([0 min(max(svm_accuracy_CV_lengthx_period12)+0.2,1)]);
    set(gca, 'FontSize', 11)
    set(gca,'box','off');% 取消右、上边框
    ylabel('Accuray', 'FontSize', 14, 'FontWeight', 'bold');
    
    nexttile
    plot(svm_accuracy_CV_lengthx_period12(T_lengthx_period1+1:end),'-','LineWidth',1.5,'Color',[0 0.4470 0.7410]);
    hold on
    for tempi=1:length(lengthx_period2_interval)
        plot(lengthx_period2_interval(tempi)*[1 1],[0 1],...
            '-','LineWidth',1,'Color',[0 0 0]);
        hold on
    end
    plot([1 T_lengthx_period2],[1 1]/nchoosek(numFrames,lengthFlag),'--','LineWidth',1,'Color',[0.4 0.4 0.4]);
    hold on
    
    set(gca,'xticklabel',[])
    set(gca,'xtick',[])
    xticks(lengthx_period2_interval);
    xticklabels({'-500 ms','Decision','Delay2 end'});
    %ylim([0 1]);
    ylim([0 min(max(svm_accuracy_CV_lengthx_period12)+0.2,1)]);
    set(gca, 'FontSize', 11)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Frames', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Accuray', 'FontSize', 14, 'FontWeight', 'bold');
    
    %% Plot fig2
    fig2 = figure('Name','Fig2','NumberTitle','off');
    %set(gcf,'Position',[450 80 800 700]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[450 80 800*0.4*0.9*1.05 700*0.4*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    % Plot heat map
    %imagesc(crossTime_svm_accuracy_lengthx_period1, [min_acc_period1 max_acc_period1]);
    %imagesc(crossTime_svm_accuracy_lengthx_period12, [min_acc_period12 max_acc_period12]);
    if lengthFlag == 1
        imagesc(crossTime_svm_accuracy_lengthx_period12, [0.16 max_acc_period12]);
    elseif lengthFlag == 2
        imagesc(crossTime_svm_accuracy_lengthx_period12, [0.07 max_acc_period12]);        
    elseif lengthFlag == 3        
        imagesc(crossTime_svm_accuracy_lengthx_period12, [0.05 max_acc_period12]);
    end
    hold on
    axis equal
    
    for tempi=1:length(lengthx_period12_interval)
        plot(lengthx_period12_interval(tempi)*[1 1],[0.5 T_lengthx_period12+0.5],...
            '-','LineWidth',1,'Color',[1 1 1]); %[0 0 0]
        hold on
        plot([0.5 T_lengthx_period12+0.5],lengthx_period12_interval(tempi)*[1 1],...
            '-','LineWidth',1,'Color',[1 1 1]);
        hold on
    end
    
    xlim([0.5 T_lengthx_period12+0.5]);
    ylim([0.5 T_lengthx_period12+0.5]);
    x_lim = xlim;
    y_lim = ylim;
    
    set(gca, 'FontSize', 9) %12
    set(gca,'YDir','normal')
    set(gca,'box','off');% 取消右、上边框
    
    xticks(lengthx_period12_interval);
    %     if lengthFlag == 1
    %         xticklabels({'Fixation','T1','D1','D1END','','Decision','D2END'});
    %     elseif lengthFlag == 2
    %         xticklabels({'Fixation','T1','T2','D1','D1END','','Decision','D2END'});
    %     elseif lengthFlag == 3
    %         xticklabels({'Fixation','T1','T2','T3','D1','D1END','','Decision','D2END'});
    %     end
    
    if lengthFlag == 1
        %xticklabels({'PreFixation','Fixation','T1','Delay1','LayoutOff',' ',' ','Decision','','','Submit',''});
        xticklabels({'Fix','','T1','Delay1',' ',' ','Decision','','','Submit',''});
    elseif lengthFlag == 2
        xticklabels({'Fix','','T1','T2','Delay1',' ',' ','Decision','','','Submit',''});
    elseif lengthFlag == 3
        %xticklabels({'PreFixation','Fixation','T1','T2','T3','Delay1',' ',' ','Decision','','','Submit',''});
        xticklabels({'Fix','','T1','T2','T3','Delay1',' ',' ','Decision','','','Submit',''});
    end
    xtickangle(0);%90
    
    yticks(lengthx_period12_interval(2:end));
    %     if lengthFlag == 1
    %         yticklabels({'T1','D1','D1END','','Decision','D2END'});
    %     elseif lengthFlag == 2
    %         yticklabels({'T1','T2','D1','D1END','','Decision','D2END'});
    %     elseif lengthFlag == 3
    %         yticklabels({'T1','T2','T3','D1','D1END','','Decision','D2END'});
    %     end
    
    if lengthFlag == 1
        %yticklabels({'Fixation','T1','Delay1','LayoutOff',' ',' ','Decision','','','Submit',''});
        yticklabels({'','T1','Delay1',' ',' ','Decision','','','Submit',''});
    elseif lengthFlag == 2
        yticklabels({'','T1','T2','Delay1',' ',' ','Decision','','','Submit',''});
    elseif lengthFlag == 3
        yticklabels({'','T1','T2','T3','Delay1',' ',' ','Decision','','','Submit',''});
    end
    ytickangle(0);%90-->60
    
%     h_xlabel = xlabel('Testing','Position',[mean(x_lim) max(y_lim)+2],'FontSize',11,'FontWeight','bold');
%     h_xlabel.Position(2) = -9;%-2-->-5-->-8
%     h_ylabel = ylabel('Training','Position',[min(x_lim)-2.5 mean(y_lim)],'FontSize',11,'FontWeight','bold');
%     h_ylabel.Position(1) = -16;
    title('Cross time decoding accuracy','FontSize',12,'FontWeight','bold');
    
    if lengthFlag == 1
        c = colorbar('Ticks', 0.16:0.1:1,'FontSize',9);
    elseif lengthFlag == 2
        c = colorbar('Ticks', 0.07:0.1:1,'FontSize',9);
    elseif lengthFlag == 3
        c = colorbar('Ticks', 0.05:0.1:1,'FontSize',9);     
    end
    c.Position(1) = c.Position(1) +0.03 + 0.05;
    
end

a = 1;

%% End
