% Chuan's 34th script (20260114)
%This script: Using input of 'pupilSizeExtraction.m', to conduct pupil size analysis
%% Initialization
% clear
close all

%% Some flags to control the following script
if_monkey_D0_Z1 = 0;% To decide whether dealing with Ding's data or Zelku's data

if_load = 1;%1

if_plot_timeCourse_eachLength_baseline = 0;
if_plot_bar_eachLength_baseline = 0;
if_plot_lengthMerge_baseline = 1;
if_plot_timeCourse_eachLength_allPeriod = 1;

if_plot_errorBar_trial0_session1 = 1;%0
if_plot_errorBar_sem0_std1_prc2 = 0;%no 2, 0
if_trialMean0_trialMedian1 = 0;%0

if_errorBarSession_raw0_subBase1 = 0;%0
if_errorBarSession_raw0_zscore1 = 1;%0

if if_plot_errorBar_trial0_session1 == 0
    if_errorBarSession_raw0_subBase1 = 0;
    if_errorBarSession_raw0_zscore1 = 0;
end
% if if_plot_errorBar_sem0_std1_prc2 == 2
%     if_trialMean0_trialMedian1 = 1;
% end

if_plot_violin0_bar1 = 1;

if_choiceMemory0_choiceMemoryCorrect1 = 0;


Baseline_timeRange = 1:1100;

delay1_timeDuration = 1100;

alpha = 0.05;%0.01



if_filter = 1;

% filter
n_order = 2;%2

band_low = 0.01/500;%0.01/500
band_high = 10/500;%10/500,4/500
[filter_b,filter_a] = butter(n_order,[band_low,band_high]);% band pass, Osorio et al., 2022
% [filter_b,filter_a] = butter(n_order,band_high,'low');% low pass
% [filter_b,filter_a] = butter(n_order,band_low,'high');% high pass


t0 = tic;

targetPATH = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\fixation';
% targetPATH = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\MonkeyDataAnalysis_ScriptsOnly_20260114A\fixation';
cd(targetPATH)

if if_monkey_D0_Z1 == 0
    monkey_name = 'D';
    searchName = 'D2p_edf_';
    
elseif if_monkey_D0_Z1 == 1
    monkey_name = 'Z';
    searchName = 'Z2p_edf_';
end

if if_load == 1
    % Load pupilSize_forAnalysis
    pupilSize_fileName = [searchName,'pupilSize_forAnalysis'];
    temp_load1 = load([pwd,'\data\',pupilSize_fileName]);
    pupilSize_forAnalysis = temp_load1.pupil_extracted;
    
    temp_fields = fieldnames(pupilSize_forAnalysis);
    for tempi=1:length(temp_fields)
        eval([temp_fields{tempi},'=pupilSize_forAnalysis.(temp_fields{tempi});']);
    end
    
    % Load pupilDiameterMeasure
    pupilDiameterMeasure_fileName = 'test_edf2022-09-14_pupilDiameterMeasure_3.5mm.mat';
    temp_load2 = load([pwd,'\data\',pupilDiameterMeasure_fileName]);
    edf0_saved = temp_load2.edf0_saved;
    
    pupilSize_factor = 3.5/mean(edf0_saved.pupilSize);%arbitray unit --> mm
    %pupilSize_factor = 1;
end



color_r_g = [166+40,97,26;1,133+20,113]/255;%bar r,g scatter r,g
color_choiceMemory = [1,133+20,113]/255;
color_choiceOffload = [166+40,97,26]/255;

temp_EdgeColor = 'none';%[0.62,0.62,0.62]-->'none'
temp_FaceAlpha = 0.35;%0.1-->0.05-->0.35


%% Preprocessing of pupil size
eyeSample_pupilSize_collapsed;
eyeDelay1_pupilSize_collapsed;
eyeBaseline_pupilSize_collapsed;
sequence_collapsed;
choice_g_index_collapsed;
choice_r_index_collapsed;
isCorrect_collapsed;
pointKindsNum;
pointShowTime;
pointShowPWM;
trial_num_multi;

pupilSizeBaseline = cell(1,length(sequence_collapsed));
pupilSizeSample = cell(1,length(sequence_collapsed));
pupilSizeDelay1 = cell(1,length(sequence_collapsed));


% % filter
% band_low = 0.01/500;%0.01/500
% band_high = 10/500;%10/500
% % [filter_b,filter_a] = butter(2,[band_low,band_high]);% band pass
% % [filter_b,filter_a] = butter(2,band_high,'low');% low pass
% [filter_b,filter_a] = butter(2,band_low,'high');% high pass


trialNum = length(sequence_collapsed);
length_seq_collapsed = zeros(1, trialNum);

% for tempi = 1:trialNum
parfor tempi = 1:trialNum
    temp_seqLength = length(sequence_collapsed{tempi});
    length_seq_collapsed(tempi) = temp_seqLength;
    
    pupilSizeBaseline{tempi}= eyeBaseline_pupilSize_collapsed{tempi}*pupilSize_factor;%arbitray unit --> mm
    pupilSizeSample{tempi}= eyeSample_pupilSize_collapsed{tempi}*pupilSize_factor;%arbitray unit --> mm
    pupilSizeDelay1{tempi}= eyeDelay1_pupilSize_collapsed{tempi}(1:delay1_timeDuration)*pupilSize_factor;%arbitray unit --> mm
    
    if if_filter == 1
        %pupilSizeBaseline{tempi} = filtfilt(filter_b,filter_a,pupilSizeBaseline{tempi});
        %pupilSizeSample{tempi} = filtfilt(filter_b,filter_a,pupilSizeSample{tempi});
        %pupilSizeDelay1{tempi} = filtfilt(filter_b,filter_a,pupilSizeDelay1{tempi});
        temp1 = pupilSizeBaseline{tempi};
        temp2 = pupilSizeSample{tempi};
        temp3 = pupilSizeDelay1{tempi};
        temp123 = [temp1,temp2,temp3];
        temp123_filter = filtfilt(filter_b,filter_a,temp123);
        temp1_filter = temp123_filter(1:length(temp1));
        temp2_filter = temp123_filter(length(temp1)+(1:length(temp2)));
        temp3_filter = temp123_filter(length(temp1)+length(temp2)+(1:length(temp3)));
        
        pupilSizeBaseline{tempi} = temp1_filter;
        pupilSizeSample{tempi} = temp2_filter;
        pupilSizeDelay1{tempi} = temp3_filter;
        
    end
end

%% sessionIndex
trial_num_multi;
sessionIndex = nan(1,sum(trial_num_multi));
for tempi=1:length(trial_num_multi)
    temp_range = sum(trial_num_multi(1:tempi-1))+1:sum(trial_num_multi(1:tempi));
    sessionIndex(temp_range) = tempi;
end

%% Extract pupil size in memory and offloading respectively
G_pupilSizeBaseline_trial = cell(1,pointKindsNum);
R_pupilSizeBaseline_trial = cell(1,pointKindsNum);
G_pupilSizeSample_trial = cell(1,pointKindsNum);
R_pupilSizeSample_trial = cell(1,pointKindsNum);
G_pupilSizeDelay1_trial = cell(1,pointKindsNum);
R_pupilSizeDelay1_trial = cell(1,pointKindsNum);

G_pupilSizeBaseline_session = cell(length(trial_num_multi),pointKindsNum);
R_pupilSizeBaseline_session = cell(length(trial_num_multi),pointKindsNum);
G_pupilSizeSample_session = cell(length(trial_num_multi),pointKindsNum);
R_pupilSizeSample_session = cell(length(trial_num_multi),pointKindsNum);
G_pupilSizeDelay1_session = cell(length(trial_num_multi),pointKindsNum);
R_pupilSizeDelay1_session = cell(length(trial_num_multi),pointKindsNum);

G_pupilSizeBaseline_session_trialMean = cell(1,pointKindsNum);
R_pupilSizeBaseline_session_trialMean = cell(1,pointKindsNum);
G_pupilSizeSample_session_trialMean = cell(1,pointKindsNum);
R_pupilSizeSample_session_trialMean = cell(1,pointKindsNum);
G_pupilSizeDelay1_session_trialMean = cell(1,pointKindsNum);
R_pupilSizeDelay1_session_trialMean = cell(1,pointKindsNum);

for temp_length = 1:pointKindsNum
    if if_choiceMemory0_choiceMemoryCorrect1 == 1
        G_trialBoolIndex = choice_g_index_collapsed == 1 & length_seq_collapsed == temp_length & isCorrect_collapsed == 1;
    elseif if_choiceMemory0_choiceMemoryCorrect1 == 0
        G_trialBoolIndex = choice_g_index_collapsed == 1 & length_seq_collapsed == temp_length;
    end
    R_trialBoolIndex = choice_r_index_collapsed == 1 & length_seq_collapsed == temp_length;
    
    G_pupilSizeBaseline_trial{temp_length} = cell2mat(pupilSizeBaseline(G_trialBoolIndex)');
    R_pupilSizeBaseline_trial{temp_length} = cell2mat(pupilSizeBaseline(R_trialBoolIndex)');
    G_pupilSizeSample_trial{temp_length} = cell2mat(pupilSizeSample(G_trialBoolIndex)');
    R_pupilSizeSample_trial{temp_length} = cell2mat(pupilSizeSample(R_trialBoolIndex)');
    G_pupilSizeDelay1_trial{temp_length} = cell2mat(pupilSizeDelay1(G_trialBoolIndex)');
    R_pupilSizeDelay1_trial{temp_length} = cell2mat(pupilSizeDelay1(R_trialBoolIndex)');
    
    
    % for each session
    temp_baselineDuration = length(pupilSizeBaseline{1});
    temp_sampleDuration = length(G_pupilSizeSample_trial{temp_length}(1,:));
    temp_delay1Duration = length(pupilSizeDelay1{1});
    
    G_pupilSizeBaseline_session_trialMean{temp_length} = nan(length(trial_num_multi),temp_baselineDuration);
    R_pupilSizeBaseline_session_trialMean{temp_length} = nan(length(trial_num_multi),temp_baselineDuration);
    G_pupilSizeSample_session_trialMean{temp_length} = nan(length(trial_num_multi),temp_sampleDuration);
    R_pupilSizeSample_session_trialMean{temp_length} = nan(length(trial_num_multi),temp_sampleDuration);
    G_pupilSizeDelay1_session_trialMean{temp_length} = nan(length(trial_num_multi),temp_delay1Duration);
    R_pupilSizeDelay1_session_trialMean{temp_length} = nan(length(trial_num_multi),temp_delay1Duration);
    
    for tempi=1:length(trial_num_multi)
        
        temp_sessionBoolIndex = sessionIndex == tempi;
        
        temp1 = cell2mat(pupilSizeBaseline(temp_sessionBoolIndex)');        
        if if_errorBarSession_raw0_subBase1 == 1
            temp1 = temp1 - temp1(:,1);
        end        
        temp2 = mean(temp1,2);
        mean_session = mean(temp2);
        std_session = std(temp2);
        
        % pupilSizeBaseline
        temp1 = cell2mat(pupilSizeBaseline(G_trialBoolIndex & temp_sessionBoolIndex)');
        if if_errorBarSession_raw0_subBase1 == 1
            temp0_G = temp1(:,1);
            temp1 = temp1 - temp0_G;
        end          
        if if_errorBarSession_raw0_zscore1 == 1
            temp1 = (temp1-mean_session)./std_session;
        end
        if if_trialMean0_trialMedian1 == 0
            temp1_trialMean = mean(temp1,1);
        elseif if_trialMean0_trialMedian1 == 1
            temp1_trialMean = median(temp1,1);
        end
        G_pupilSizeBaseline_session_trialMean{temp_length}(tempi,:) = temp1_trialMean;
        G_pupilSizeBaseline_session{tempi,temp_length} = temp1;
        
        temp1 = cell2mat(pupilSizeBaseline(R_trialBoolIndex & temp_sessionBoolIndex)');
        if if_errorBarSession_raw0_subBase1 == 1
            temp0_R = temp1(:,1);            
            temp1 = temp1 - temp0_R;
        end           
        if if_errorBarSession_raw0_zscore1 == 1
            temp1 = (temp1-mean_session)./std_session;
        end
        if if_trialMean0_trialMedian1 == 0
            temp1_trialMean = mean(temp1,1);
        elseif if_trialMean0_trialMedian1 == 1
            temp1_trialMean = median(temp1,1);
        end
        R_pupilSizeBaseline_session_trialMean{temp_length}(tempi,:) = temp1_trialMean;
        R_pupilSizeBaseline_session{tempi,temp_length} = temp1;
        
        % pupilSizeSample
        temp1 = cell2mat(pupilSizeSample(G_trialBoolIndex & temp_sessionBoolIndex)');
        if if_errorBarSession_raw0_subBase1 == 1
            temp1 = temp1 - temp0_G;
        end           
        if if_errorBarSession_raw0_zscore1 == 1
            temp1 = (temp1-mean_session)./std_session;
        end
        if if_trialMean0_trialMedian1 == 0
            temp1_trialMean = mean(temp1,1);
        elseif if_trialMean0_trialMedian1 == 1
            temp1_trialMean = median(temp1,1);
        end
        G_pupilSizeSample_session_trialMean{temp_length}(tempi,:) = temp1_trialMean;
        G_pupilSizeSample_session{tempi,temp_length} = temp1;
        
        temp1 = cell2mat(pupilSizeSample(R_trialBoolIndex & temp_sessionBoolIndex)');
        if if_errorBarSession_raw0_subBase1 == 1
            temp1 = temp1 - temp0_R;
        end           
        if if_errorBarSession_raw0_zscore1 == 1
            temp1 = (temp1-mean_session)./std_session;
        end
        if if_trialMean0_trialMedian1 == 0
            temp1_trialMean = mean(temp1,1);
        elseif if_trialMean0_trialMedian1 == 1
            temp1_trialMean = median(temp1,1);
        end
        R_pupilSizeSample_session_trialMean{temp_length}(tempi,:) = temp1_trialMean;
        R_pupilSizeSample_session{tempi,temp_length} = temp1;
        
        % pupilSizeDelay1
        temp1 = cell2mat(pupilSizeDelay1(G_trialBoolIndex & temp_sessionBoolIndex)');
        if if_errorBarSession_raw0_subBase1 == 1
            temp1 = temp1 - temp0_G;
        end           
        if if_errorBarSession_raw0_zscore1 == 1
            temp1 = (temp1-mean_session)./std_session;
        end
        if if_trialMean0_trialMedian1 == 0
            temp1_trialMean = mean(temp1,1);
        elseif if_trialMean0_trialMedian1 == 1
            temp1_trialMean = median(temp1,1);
        end
        G_pupilSizeDelay1_session_trialMean{temp_length}(tempi,:) = temp1_trialMean;
        G_pupilSizeDelay1_session{tempi,temp_length} = temp1;
        
        temp1 = cell2mat(pupilSizeDelay1(R_trialBoolIndex & temp_sessionBoolIndex)');
        if if_errorBarSession_raw0_subBase1 == 1
            temp1 = temp1 - temp0_R;
        end           
        if if_errorBarSession_raw0_zscore1 == 1
            temp1 = (temp1-mean_session)./std_session;
        end
        if if_trialMean0_trialMedian1 == 0
            temp1_trialMean = mean(temp1,1);
        elseif if_trialMean0_trialMedian1 == 1
            temp1_trialMean = median(temp1,1);
        end
        R_pupilSizeDelay1_session_trialMean{temp_length}(tempi,:) = temp1_trialMean;
        R_pupilSizeDelay1_session{tempi,temp_length} = temp1;
        
    end
    
end

if if_plot_errorBar_trial0_session1 == 0
    G_pupilSizeBaseline = G_pupilSizeBaseline_trial;
    R_pupilSizeBaseline = R_pupilSizeBaseline_trial;
    G_pupilSizeSample = G_pupilSizeSample_trial;
    R_pupilSizeSample = R_pupilSizeSample_trial;
    G_pupilSizeDelay1 = G_pupilSizeDelay1_trial;
    R_pupilSizeDelay1 = R_pupilSizeDelay1_trial;
elseif if_plot_errorBar_trial0_session1 == 1
    G_pupilSizeBaseline = G_pupilSizeBaseline_session_trialMean;
    R_pupilSizeBaseline = R_pupilSizeBaseline_session_trialMean;
    G_pupilSizeSample = G_pupilSizeSample_session_trialMean;
    R_pupilSizeSample = R_pupilSizeSample_session_trialMean;
    G_pupilSizeDelay1 = G_pupilSizeDelay1_session_trialMean;
    R_pupilSizeDelay1 = R_pupilSizeDelay1_session_trialMean;
end

% % FDR
% R_G_p_Baseline = cell(1,pointKindsNum);
% FDR_Baseline = cell(1,pointKindsNum);
% FDR_p_Baseline = cell(1,pointKindsNum);
% for temp_length = 1:pointKindsNum
%     [Rrow2,Rcol2] = size(R_pupilSizeBaseline{temp_length});
%     for i = 1:Rcol2
%         [R_G_sign,p,~,~] = ttest2(R_pupilSizeBaseline{temp_length}(:,i),G_pupilSizeBaseline{temp_length}(:,i),'Alpha',alpha);
%         R_G_p_Baseline{temp_length} = [R_G_p_Baseline{temp_length},p];
%     end
%     [pthr, pcor, FDR_Baseline{temp_length}] = fdr(R_G_p{temp_length});
%     temp_FDR_p_Num = find(FDR_Baseline{temp_length}<alpha);
%     FDR_p_Baseline{temp_length} = zeros(1, length(FDR_Baseline{temp_length}));
%     FDR_p_Baseline{temp_length}(temp_FDR_p_Num) = 1;
% end

% Mean & SEM
meanR_pupilSizeBaseline = cell(1,pointKindsNum);
meanG_pupilSizeBaseline =cell(1,pointKindsNum);
semR_pupilSizeBaseline =cell(1,pointKindsNum);
semG_pupilSizeBaseline = cell(1,pointKindsNum);

meanR_pupilSizeSample = cell(1,pointKindsNum);
meanG_pupilSizeSample =cell(1,pointKindsNum);
semR_pupilSizeSample =cell(1,pointKindsNum);
semG_pupilSizeSample = cell(1,pointKindsNum);

meanR_pupilSizeDelay1 = cell(1,pointKindsNum);
meanG_pupilSizeDelay1 =cell(1,pointKindsNum);
semR_pupilSizeDelay1 =cell(1,pointKindsNum);
semG_pupilSizeDelay1 = cell(1,pointKindsNum);

for temp_length = 1:pointKindsNum
    [temp_r_size,~] = size(R_pupilSizeBaseline{temp_length});
    [temp_g_size,~] = size(G_pupilSizeBaseline{temp_length});
    
    meanR_pupilSizeBaseline{temp_length} = mean(R_pupilSizeBaseline{temp_length},1);
    meanG_pupilSizeBaseline{temp_length} = mean(G_pupilSizeBaseline{temp_length},1);
    if if_plot_errorBar_sem0_std1_prc2 == 0
        semR_pupilSizeBaseline{temp_length} = std(R_pupilSizeBaseline{temp_length})/sqrt(temp_r_size);
        semG_pupilSizeBaseline{temp_length} = std(G_pupilSizeBaseline{temp_length})/sqrt(temp_g_size);
    elseif if_plot_errorBar_sem0_std1_prc2 == 1
        semR_pupilSizeBaseline{temp_length} = std(R_pupilSizeBaseline{temp_length});
        semG_pupilSizeBaseline{temp_length} = std(G_pupilSizeBaseline{temp_length});
        %     elseif if_plot_errorBar_sem0_std1_prc2 == 2
        %         temp_prc1 = prctile(R_pupilSizeBaseline{temp_length},75,1);
        %         temp_prc2 = prctile(R_pupilSizeBaseline{temp_length},25,1);
        %         temp_50prc = prctile(R_pupilSizeBaseline{temp_length},50,1);
        %
        %         temptemp1 = abs(temp_prc1-temp_50prc);
        %         temptemp2 = abs(temp_prc2-temp_50prc);
        %         temptemp12_max = max([temptemp1;temptemp2],[],1);
        %
        %         semR_pupilSizeBaseline{temp_length} = std(R_pupilSizeBaseline{temp_length})/sqrt(temp_r_size);
        %         semG_pupilSizeBaseline{temp_length} = std(G_pupilSizeBaseline{temp_length})/sqrt(temp_g_size);
    end
    
    meanR_pupilSizeSample{temp_length} = mean(R_pupilSizeSample{temp_length},1);
    meanG_pupilSizeSample{temp_length} = mean(G_pupilSizeSample{temp_length},1);
    if if_plot_errorBar_sem0_std1_prc2 == 0
        semR_pupilSizeSample{temp_length} = std(R_pupilSizeSample{temp_length})/sqrt(temp_r_size);
        semG_pupilSizeSample{temp_length} = std(G_pupilSizeSample{temp_length})/sqrt(temp_g_size);
    elseif if_plot_errorBar_sem0_std1_prc2 == 1
        semR_pupilSizeSample{temp_length} = std(R_pupilSizeSample{temp_length});
        semG_pupilSizeSample{temp_length} = std(G_pupilSizeSample{temp_length});
    end
    
    meanR_pupilSizeDelay1{temp_length} = mean(R_pupilSizeDelay1{temp_length},1);
    meanG_pupilSizeDelay1{temp_length} = mean(G_pupilSizeDelay1{temp_length},1);
    if if_plot_errorBar_sem0_std1_prc2 == 0
        semR_pupilSizeDelay1{temp_length} = std(R_pupilSizeDelay1{temp_length})/sqrt(temp_r_size);
        semG_pupilSizeDelay1{temp_length} = std(G_pupilSizeDelay1{temp_length})/sqrt(temp_g_size);
    elseif if_plot_errorBar_sem0_std1_prc2 == 1
        semR_pupilSizeDelay1{temp_length} = std(R_pupilSizeDelay1{temp_length});
        semG_pupilSizeDelay1{temp_length} = std(G_pupilSizeDelay1{temp_length});
    end
end


%% Plot if_plot_timeCourse_eachLength_baseline
if if_plot_timeCourse_eachLength_baseline == 1
    maxy_all = [];
    miny_all = [];
    for temp_length = 1:pointKindsNum
        maxy = max([meanR_pupilSizeBaseline{temp_length}+semR_pupilSizeBaseline{temp_length},meanG_pupilSizeBaseline{temp_length}+semG_pupilSizeBaseline{temp_length}]);%0.1
        miny = min([meanR_pupilSizeBaseline{temp_length}-semR_pupilSizeBaseline{temp_length},meanG_pupilSizeBaseline{temp_length}-semG_pupilSizeBaseline{temp_length}]);%-0.2
        maxy_all = [maxy_all,maxy]; %#ok<*AGROW>
        miny_all = [miny_all,miny];
    end
    % temp_y_min = roundn(min(miny_all),-1);
    % temp_y_max = roundn(max(maxy_all),-1);
    temp_y_min = min(miny_all);
    temp_y_max = max(maxy_all);
    rangeyy = temp_y_max-temp_y_min;
    % breaklength = 100;
    
    %fig1 = figure(1);
    fig = figure('Name',sprintf('fig'),'NumberTitle','off');
    if pointKindsNum == 3
        set(gcf,'Position',[0 50 1650 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,pointKindsNum,'TileSpacing','Compact','Padding','Compact');
    elseif pointKindsNum == 4
        %set(gcf,'Position',[0 50 1100 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[0 50 1100*0.7*0.9 800*0.7*0.4*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
        
        %t.Title.String = sprintf('Monkey %s baseline pupil size',monkey_name);
        t.Title.String = sprintf('Monkey %s baseline pupil size\n%s',monkey_name,pupilSize_fileName);
        
        t.Title.FontSize = 11;
        t.Title.Interpreter = 'none';
    end
    
    for temp_length = 1:pointKindsNum
        nexttile
        
        x = 1:length(meanG_pupilSizeBaseline{temp_length});
        
        [temp_x_min,temp_x_max] = bounds(x);
        
        h_line = [];
        
        % plot baselinePeriod
        x = x; %#ok<*ASGSL>
        y = meanG_pupilSizeBaseline{temp_length};
        h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');];
        hold on
        y_sem = semG_pupilSizeBaseline{temp_length};
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        x = x;
        y = meanR_pupilSizeBaseline{temp_length};
        h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');];
        hold on
        y_sem = semR_pupilSizeBaseline{temp_length};
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        %plot([1 1]*length(Baseline_timeRange),[temp_y_min temp_y_max],'Color','k','LineWidth',0.5)
        
        
        
        %         for i = 1:length(FDR_p_Baseline{temp_length})
        %             if FDR_p_Baseline{temp_length}(i) == 1
        %                 text_signi = '*';
        %                 text(i,temp_y_min-(temp_y_max-temp_y_min)*0.05,text_signi,'FontSize',8,'Color',[0 0 0]);%[0/255 133/255 0/255]
        %             end
        %         end
        %
        %         text_p = sprintf('P < %s(FDR)',num2str(alpha));
        %         text(temp_x_max-(temp_x_max-temp_x_min)*0.5,temp_y_min-(temp_y_max-temp_y_min)*0.10,text_p,'FontSize',7);
        
        
        if temp_length == 1
            if if_choiceMemory0_choiceMemoryCorrect1 == 0
                le = legend(h_line,'Memory','Offload','Location','northeast','fontsize',7);
            elseif if_choiceMemory0_choiceMemoryCorrect1 == 1
                le = legend(h_line,'Memory correct','Offload','Location','northeast','fontsize',7);
            end
            le.ItemTokenSize = ones(1,3)*10;
            legend('boxoff');
        end
        
        
        set(gca,'XTick',[1,length(Baseline_timeRange)])
        
        set(gca,'linewidth',1.5)
        %set(gca,'color',backgounrdColor);
        set(gca, 'FontSize', 8)
        set(gca,'box','off');
        
        %xlim([0,length(Baseline_timeRange)+ceil(pointShowTime*pointShowPWM)])
        
        xlim([temp_x_min,temp_x_max]);
        
        
        set(gca,'xticklabels',{'Fix','T1'})
        
        ytickk = [temp_y_min:0.1:0,0.1:0.1:temp_y_max];
        set(gca,'YTick',ytickk)
        if temp_length == 1
            ylabel('Pupil size (mm)','FontSize',9);
        end
        set(gca,'yticklabels',string(ytickk))
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.18,temp_y_max+(temp_y_max-temp_y_min)*0.1]);
        
        %title(['Monkey ',monkey_name,' baseline pupil size'],'Fontsize',9);
        title(['Length',num2str(temp_length)],'Fontsize',9);
        
    end
    
end


%% Plot if_plot_bar_eachLength_baseline
if if_plot_bar_eachLength_baseline == 1
    %% statistics test
    mean_G= zeros(1,pointKindsNum);
    mean_R = zeros(1,pointKindsNum);
    sem_G = zeros(1,pointKindsNum);
    sem_R = zeros(1,pointKindsNum);
    signi = zeros(1,pointKindsNum);
    temp_p = nan(1,pointKindsNum);
    alpha2 = 0.001;
    for temp_length = 1:pointKindsNum
        temp_1 = mean(G_pupilSizeBaseline{temp_length},2);
        mean_G(temp_length) = mean(temp_1);
        if if_plot_errorBar_sem0_std1_prc2 == 0
            sem_G(temp_length) = std(temp_1)./sqrt(length(temp_1));
        elseif if_plot_errorBar_sem0_std1_prc2 == 1
            sem_G(temp_length) = std(temp_1);
        end
        
        temp_2 = mean(R_pupilSizeBaseline{temp_length},2);
        mean_R(temp_length) = mean(temp_2);
        if if_plot_errorBar_sem0_std1_prc2 == 0
            sem_R(temp_length) = std(temp_2)./sqrt(length(temp_2));
        elseif if_plot_errorBar_sem0_std1_prc2 == 1
            sem_R(temp_length) = std(temp_2);
        end
        
        [~,temp_p(temp_length)] = ttest2(temp_1,temp_2);
    end
    
    temp_y_max = max([mean_R+sem_R,mean_G+sem_G]);
    temp_y_min = min([mean_R-sem_R,mean_G-sem_G]);
    
    %% Plot the mean pupil size as bar graph
    %fig2 = figure(2);
    fig = figure('Name',sprintf('fig'),'NumberTitle','off');
    set(gcf,'Position',[700 50 400 250]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    width = 0.4;
    xlim([1 pointKindsNum*4])
    xtickk = [];
    for temp_length = 1:pointKindsNum
        xtickk = [xtickk,temp_length*3-0.5,temp_length*3+0.5];
    end
    for temp_length = 1:pointKindsNum
        x_seqG = temp_length*3-0.5;
        x_seqR = temp_length*3+0.5;
        bar1 = bar(x_seqG,abs(mean_G(temp_length)),width,'FaceColor',color_choiceMemory);
        hold on
        bar2 = bar(x_seqR,abs(mean_R(temp_length)),width,'FaceColor',color_choiceOffload);
        errorbar1 = errorbar(x_seqG,abs(mean_G(temp_length)),abs(sem_G(temp_length)),'Color','[0 0 0]','LineStyle','none');
        errorbar2 = errorbar(x_seqR,abs(mean_R(temp_length)),abs(sem_R(temp_length)),'Color','[0 0 0]','LineStyle','none');
        
        tempTxt = sprintf('');
        if temp_p(temp_length) < 0.001
            tempTxt = sprintf('***');
        elseif temp_p(temp_length) < 0.01
            tempTxt = sprintf('**');
        elseif temp_p(temp_length) < 0.05
            tempTxt = sprintf('*');
        end
        text(x_seqG+0.49,temp_y_max+(temp_y_max-temp_y_min)*0.03,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        txt_length = ['length',num2str(temp_length)];
        text(x_seqG-0.19,temp_y_min-(temp_y_max-temp_y_min)*0.09,txt_length,'FontSize',9)
    end
    
    if if_choiceMemory0_choiceMemoryCorrect1 == 0
        legend('Memory','Offloading','Location','northeast','fontsize',7)
    elseif if_choiceMemory0_choiceMemoryCorrect1 == 1
        legend('Memory correct','Offloading','Location','northeast','fontsize',7)
    end
    legend('boxoff');
    
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.05,temp_y_max+(temp_y_max-temp_y_min)*0.4]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');
    
    set(gca,'XTick',xtickk)
    xticklabels([ ])
    ylabel('Pupil size (mm)','Fontsize',9)
    title(['Monkey ',monkey_name,' baseline pupil size'],'Fontsize',9);
    
end


%% Plot length-merged results
if if_plot_lengthMerge_baseline == 1
    
    G_pupilSizeBaseline;
    R_pupilSizeBaseline;
    
    G_correct_pupilSize_lengthMerged = [];
    for tempi=1:length(G_pupilSizeBaseline)
        G_correct_pupilSize_lengthMerged = [G_correct_pupilSize_lengthMerged;G_pupilSizeBaseline{tempi}];
    end
    
    R_pupilSize_lengthMerged = [];
    for tempi=1:length(R_pupilSizeBaseline)
        R_pupilSize_lengthMerged = [R_pupilSize_lengthMerged;R_pupilSizeBaseline{tempi}];
    end
    
    G_correct_pupilSize_lengthMerged;
    R_pupilSize_lengthMerged;
    
    G_correct_pupilSize_lengthMerged_trialMean = mean(G_correct_pupilSize_lengthMerged,1);
    if if_plot_errorBar_sem0_std1_prc2 == 0
        G_correct_pupilSize_lengthMerged_trialSEM = std(G_correct_pupilSize_lengthMerged,1,1)./sqrt(size(G_correct_pupilSize_lengthMerged,1));
    elseif if_plot_errorBar_sem0_std1_prc2 == 1
        G_correct_pupilSize_lengthMerged_trialSEM = std(G_correct_pupilSize_lengthMerged,1,1);
    end
    
    R_pupilSize_lengthMerged_trialMean = mean(R_pupilSize_lengthMerged,1);
    if if_plot_errorBar_sem0_std1_prc2 == 0
        R_pupilSize_lengthMerged_trialSEM = std(R_pupilSize_lengthMerged,1,1)./sqrt(size(R_pupilSize_lengthMerged,1));
    elseif if_plot_errorBar_sem0_std1_prc2 == 1
        R_pupilSize_lengthMerged_trialSEM = std(R_pupilSize_lengthMerged,1,1);
    end
    
    
    
    G_correct_pupilSize_lengthMerged_timeMean = mean(G_correct_pupilSize_lengthMerged,2);
    R_pupilSize_lengthMerged_timeMean = mean(R_pupilSize_lengthMerged,2);
    
    %% Plot
    if true
        fig = figure('Name',sprintf('fig'),'NumberTitle','off');
        
        %set(gcf,'Position',[0 50 1100*0.7*0.9*0.7 800*0.7*0.4*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[0 50 1100*0.7*0.9*0.7*0.65*1.04*1.08*0.92*0.985 800*0.7*0.4*0.9*1.09*0.96]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[0 50 1100*0.7*0.9*0.7*0.65*1.04*1.08*0.92*0.985 800*0.7*0.4*0.9*1.09*0.96*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[0 50 1100*0.7*0.9*0.7*0.65*1.04*1.08*0.92*0.985 800*0.7*0.4*0.9*1.09*0.96*0.8*1.09]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[0 50 1100*0.7*0.9*0.7*0.65*1.04*1.08*0.92*0.985 800*0.7*0.4*0.9*1.09*0.96*0.8*1.09*1.02]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        %t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
        t = tiledlayout(1,11,'TileSpacing','Compact','Padding','Compact');
        
        if if_plot_errorBar_trial0_session1 == 0
            t.Title.String = sprintf('Monkey %s baseline pupil size, length-merged\n%s, error bar from trials',monkey_name,pupilSize_fileName);
        elseif if_plot_errorBar_trial0_session1 == 1
            t.Title.String = sprintf('Monkey %s baseline pupil size, length-merged\n%s, error bar from sessions',monkey_name,pupilSize_fileName);
        end
        t.Title.FontSize = 9;
        t.Title.Interpreter = 'none';
        
        temp_y_min = min([G_correct_pupilSize_lengthMerged_trialMean-G_correct_pupilSize_lengthMerged_trialSEM,R_pupilSize_lengthMerged_trialMean-R_pupilSize_lengthMerged_trialSEM]);
        temp_y_max = max([G_correct_pupilSize_lengthMerged_trialMean+G_correct_pupilSize_lengthMerged_trialSEM,R_pupilSize_lengthMerged_trialMean+R_pupilSize_lengthMerged_trialSEM]);
        
        nexttile([1 6])
        
        x = 1:length(G_correct_pupilSize_lengthMerged_trialMean);
        
        [temp_x_min,temp_x_max] = bounds(x);
        
        h_line = [];
        
        % plot baselinePeriod
        x = x; %#ok<*ASGSL>
        y = G_correct_pupilSize_lengthMerged_trialMean;
        h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');];
        hold on
        y_sem = G_correct_pupilSize_lengthMerged_trialSEM;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        x = x;
        y = R_pupilSize_lengthMerged_trialMean;
        h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');];
        hold on
        y_sem = R_pupilSize_lengthMerged_trialSEM;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        if if_choiceMemory0_choiceMemoryCorrect1 == 0
            le = legend(h_line,'Memory','Offload','Location','northeast','fontsize',6.5);
            %le = legend(h_line,'Choice-memory','Choice-offload','Location','northeast','fontsize',6.5);
        elseif if_choiceMemory0_choiceMemoryCorrect1 == 1
            le = legend(h_line,'Memory correct','Offload','Location','northeast','fontsize',7);
        end
        le.ItemTokenSize = ones(1,3)*10;
        legend('boxoff');
        
        
        set(gca,'XTick',[1,length(Baseline_timeRange)])
        
        xlim([temp_x_min,temp_x_max]);
        %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.18,temp_y_max+(temp_y_max-temp_y_min)*0.1]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1,temp_y_max+(temp_y_max-temp_y_min)*0.1]);
        
        set(gca,'xticklabels',{'Fixation','T1'})
        
        if if_monkey_D0_Z1 == 0
            %yticks([5.3:0.1:5.5]);
        end
        
        set(gca,'linewidth',1.5)
        set(gca, 'FontSize', 8)
        set(gca,'box','off');
        %ylabel('Pupil size (mm)','FontSize',9);
        if if_errorBarSession_raw0_zscore1 == 0
            ylabel('Pupil size (mm)','FontSize',9);
        elseif if_errorBarSession_raw0_zscore1 == 1
            ylabel('Pupil size (z-scored)','FontSize',9);
        end
        title(sprintf('Time course'),'Fontsize',9);
        
        
        nexttile([1 5])
        
        
        temp_1 = G_correct_pupilSize_lengthMerged_timeMean;
        temp_2 = R_pupilSize_lengthMerged_timeMean;
        
        [~,temp_p] = ttest2(temp_1,temp_2);
        
        if if_plot_errorBar_sem0_std1_prc2 == 0
            temp1_SEM = std(temp_1)/sqrt(length(temp_1));
            temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        elseif if_plot_errorBar_sem0_std1_prc2 == 1
            temp1_SEM = std(temp_1);
            temp2_SEM = std(temp_2);
        end
        
        if if_plot_violin0_bar1 == 0
            temp_y_min = min([temp_1;temp_2]);
            temp_y_max = max([temp_1;temp_2]);
            
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
            
        elseif if_plot_violin0_bar1 == 1
            temp_y_min = min([mean(temp_1)-temp1_SEM;mean(temp_2)-temp2_SEM]);
            temp_y_max = max([mean(temp_1)+temp1_SEM;mean(temp_2)+temp2_SEM]);
            
            bar([1 2], [mean(temp_1) mean(temp_2)], ...
                'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',0.5);
            hold on
            errorbar([1 2], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 10);%15
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
        xlim([0.30 2.70])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.20]);%0.20
        set(gca, 'FontSize', 8);%10
        
        xticks([1 2]);
        
        if if_choiceMemory0_choiceMemoryCorrect1 == 0
            xtl = ["Memory"; "Offload"];
        elseif if_choiceMemory0_choiceMemoryCorrect1 == 1
            xtl = ["Memory-correct"; "Offload"];
        end
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.22;%0.56,0.4
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',7.5);%9-->7.5
        set(gca,'xticklabel','');
        
        set(gca,'box','off');% 取消右、上边框
        %ylabel('Pupil size (mm)','FontSize',9);
        if if_errorBarSession_raw0_zscore1 == 0
            ylabel('Pupil size (mm)','FontSize',9);
        elseif if_errorBarSession_raw0_zscore1 == 1
            ylabel('Pupil size (z-scored)','FontSize',9);
        end
        title(sprintf('Time averaged'),'FontSize',9,'FontWeight','normal');
        
    end
    
    
    %% Plot
    if true
        fig = figure('Name',sprintf('fig'),'NumberTitle','off');
        
        set(gcf,'Position',[0 50 320.9*0.65*0.6*0.9 187.6*1.8*0.8*0.92*1.01]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');       
        
        temp_y_min = min([G_correct_pupilSize_lengthMerged_trialMean-G_correct_pupilSize_lengthMerged_trialSEM,R_pupilSize_lengthMerged_trialMean-R_pupilSize_lengthMerged_trialSEM]);
        temp_y_max = max([G_correct_pupilSize_lengthMerged_trialMean+G_correct_pupilSize_lengthMerged_trialSEM,R_pupilSize_lengthMerged_trialMean+R_pupilSize_lengthMerged_trialSEM]);
        
        
        nexttile
                
        temp_1 = G_correct_pupilSize_lengthMerged_timeMean;
        temp_2 = R_pupilSize_lengthMerged_timeMean;
        
        [~,temp_p] = ttest2(temp_1,temp_2);
        
        %if_plot_errorBar_sem0_std1_prc2 = if_plot_errorBar_sem0_std1_prc2;
        %if_plot_errorBar_sem0_std1_prc2 = 1;
        
        if if_plot_errorBar_sem0_std1_prc2 == 0
            temp1_SEM = std(temp_1)/sqrt(length(temp_1));
            temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        elseif if_plot_errorBar_sem0_std1_prc2 == 1
            temp1_SEM = std(temp_1);
            temp2_SEM = std(temp_2);
        end
        
        if if_plot_violin0_bar1 == 0
            temp_y_min = min([temp_1;temp_2]);
            temp_y_max = max([temp_1;temp_2]);
            
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
            
        elseif if_plot_violin0_bar1 == 1
            temp_y_min = min([mean(temp_1)-temp1_SEM;mean(temp_2)-temp2_SEM]);
            temp_y_max = max([mean(temp_1)+temp1_SEM;mean(temp_2)+temp2_SEM]);
            
            temp_y_min = -0.5375;
            temp_y_max = 0.5958;
            
            %bar([1 2], [mean(temp_1) mean(temp_2)],'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',0.5);
            bar([1 2], [mean(temp_1) mean(temp_2)],'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0],'LineWidth',0.5);            
            hold on
            %errorbar([1 2], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 10);%15
            %errorbar([1 2], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 1.5, 'CapSize', 8);%15
            errorbar([1 2], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'CapSize', 8);%15
                                   
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
        xlim([0.30 2.70])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.20]);%0.20
        set(gca, 'FontSize', 8);%10
        
        xticks([1 2]);
        
        if if_choiceMemory0_choiceMemoryCorrect1 == 0
            xtl = ["Memory"; "Offload"];
        elseif if_choiceMemory0_choiceMemoryCorrect1 == 1
            xtl = ["Memory-correct"; "Offload"];
        end
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.19;%0.56,0.4,0.17
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->7.5
        set(gca,'xticklabel','');
        
        set(gca,'box','off');% 取消右、上边框
        %ylabel('Pupil size (mm)','FontSize',9);
        if if_errorBarSession_raw0_zscore1 == 0
            ylabel('Pupil size (mm)','FontSize',9);
        elseif if_errorBarSession_raw0_zscore1 == 1
            ylabel('Pupil size (z-scored)','FontSize',9);
        end
        %title(sprintf('Time averaged'),'FontSize',9,'FontWeight','normal');
        
    end
    
    
end


%% Plot if_plot_timeCourse_eachLength_allPeriod
if if_plot_timeCourse_eachLength_allPeriod == 1
    
    maxy_all = [];
    miny_all = [];
    for temp_length = 1:pointKindsNum
        maxy = max([meanR_pupilSizeBaseline{temp_length}+semR_pupilSizeBaseline{temp_length},meanG_pupilSizeBaseline{temp_length}+semG_pupilSizeBaseline{temp_length},...
            meanR_pupilSizeSample{temp_length}+semR_pupilSizeSample{temp_length},meanG_pupilSizeSample{temp_length}+semG_pupilSizeSample{temp_length},...
            meanR_pupilSizeDelay1{temp_length}+semR_pupilSizeDelay1{temp_length},meanG_pupilSizeDelay1{temp_length}+semG_pupilSizeDelay1{temp_length}
            ]);
        miny = min([meanR_pupilSizeBaseline{temp_length}-semR_pupilSizeBaseline{temp_length},meanG_pupilSizeBaseline{temp_length}-semG_pupilSizeBaseline{temp_length},...
            meanR_pupilSizeSample{temp_length}-semR_pupilSizeSample{temp_length},meanG_pupilSizeSample{temp_length}-semG_pupilSizeSample{temp_length},...
            meanR_pupilSizeDelay1{temp_length}-semR_pupilSizeDelay1{temp_length},meanG_pupilSizeDelay1{temp_length}-semG_pupilSizeDelay1{temp_length}
            ]);
        maxy_all = [maxy_all,maxy]; %#ok<*AGROW>
        miny_all = [miny_all,miny];
    end
    temp_y_min = min(miny_all);
    temp_y_max = max(maxy_all);
    rangeyy = temp_y_max-temp_y_min;
    
    if true
        fig = figure('Name',sprintf('fig'),'NumberTitle','off'); %#ok<*NASGU>
        
        %set(gcf,'Position',[500 50 800*0.7*0.4*0.9*2*1.5*0.975*1.10 1100*0.7*0.9*0.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[500 50 800*0.7*0.4*0.9*2*1.5*0.975*1.10 1100*0.7*0.9*0.5*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[500 50 800*0.7*0.4*0.9*2*1.5*0.975*1.10*0.6*0.81*1.03 1100*0.7*0.9*0.5*0.85*1.8*1.03]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        %t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
        t = tiledlayout(4,1,'TileSpacing','Compact','Padding','Compact');
        
        if if_plot_errorBar_trial0_session1 == 0
            t.Title.String = sprintf('Monkey %s pupil size\n%s, error bar from trials',monkey_name,pupilSize_fileName);
        elseif if_plot_errorBar_trial0_session1 == 1
            t.Title.String = sprintf('Monkey %s pupil size\n%s, error bar from sessions',monkey_name,pupilSize_fileName);
        end
        
        t.Title.FontSize = 11;
        t.Title.Interpreter = 'none';
        
        
        for temp_length = 1:pointKindsNum
            nexttile
            
            tempMeanG_pupilSize = [meanG_pupilSizeBaseline{temp_length},meanG_pupilSizeSample{temp_length},meanG_pupilSizeDelay1{temp_length}];
            tempMeanR_pupilSize = [meanR_pupilSizeBaseline{temp_length},meanR_pupilSizeSample{temp_length},meanR_pupilSizeDelay1{temp_length}];
            tempSEMG_pupilSize = [semG_pupilSizeBaseline{temp_length},semG_pupilSizeSample{temp_length},semG_pupilSizeDelay1{temp_length}];
            tempSEMR_pupilSize = [semR_pupilSizeBaseline{temp_length},semR_pupilSizeSample{temp_length},semR_pupilSizeDelay1{temp_length}];
            
            x = 1:length(tempMeanG_pupilSize);
            
            [temp_x_min,temp_x_max] = bounds(x);
            
            h_line = [];
            
            x = x; %#ok<*ASGSL>
            y = tempMeanG_pupilSize;
            h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');];
            hold on
            y_sem = tempSEMG_pupilSize;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
            hold on
            
            
            x = x;
            y = tempMeanR_pupilSize;
            h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');];
            hold on
            y_sem = tempSEMR_pupilSize;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
            hold on
            
            
            for tempi=1:temp_length
                temp_x = [1+length(Baseline_timeRange)+pointShowTime*(tempi-1) temp_y_min-(temp_y_max-temp_y_min)*0.1;...
                    1+length(Baseline_timeRange)+pointShowTime*(tempi-1) temp_y_max+(temp_y_max-temp_y_min)*0.1;...
                    1+length(Baseline_timeRange)+pointShowTime*(tempi-1)+200 temp_y_max+(temp_y_max-temp_y_min)*0.1;...
                    1+length(Baseline_timeRange)+pointShowTime*(tempi-1)+200 temp_y_min-(temp_y_max-temp_y_min)*0.1];
                temp_y = [1 2 3 4];
                patch('Faces',temp_y,'Vertices',temp_x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.25,'EdgeColor','none');%0.1
                hold on
            end
            
            
            if temp_length == 1
                if if_choiceMemory0_choiceMemoryCorrect1 == 0
                    %le = legend(h_line,'Choice-memory','Choice-offload','Location','northeast','fontsize',6.5);
                    le = legend(h_line,'Memory','Offload','Location','northeast','fontsize',6.5);
                elseif if_choiceMemory0_choiceMemoryCorrect1 == 1
                    le = legend(h_line,'Memory correct','Offload','Location','northeast','fontsize',7);
                end
                le.ItemTokenSize = ones(1,3)*10;
                legend('boxoff');
            end
            
            if temp_length == 1
                set(gca,'xticklabels',{'Fixation','T1','Delay-on'})
                set(gca,'XTick',[1,1+length(Baseline_timeRange),...
                    1+length(Baseline_timeRange)+pointShowTime*1,...
                    1+length(Baseline_timeRange)+pointShowTime*1+delay1_timeDuration])
            elseif temp_length == 2
                set(gca,'xticklabels',{'Fixation','T1','T2','Delay-on'})
                set(gca,'XTick',[1,1+length(Baseline_timeRange),...
                    1+length(Baseline_timeRange)+pointShowTime*1,...
                    1+length(Baseline_timeRange)+pointShowTime*2,...
                    1+length(Baseline_timeRange)+pointShowTime*2+delay1_timeDuration])
            elseif temp_length == 3
                set(gca,'xticklabels',{'Fixation','T1','T2','T3','Delay-on'})
                set(gca,'XTick',[1,1+length(Baseline_timeRange),...
                    1+length(Baseline_timeRange)+pointShowTime*1,...
                    1+length(Baseline_timeRange)+pointShowTime*2,...
                    1+length(Baseline_timeRange)+pointShowTime*3,...
                    1+length(Baseline_timeRange)+pointShowTime*3+delay1_timeDuration])
            elseif temp_length == 4
                set(gca,'xticklabels',{'Fixation','T1','T2','T3','T4','Delay-on'})
                set(gca,'XTick',[1,1+length(Baseline_timeRange),...
                    1+length(Baseline_timeRange)+pointShowTime*1,...
                    1+length(Baseline_timeRange)+pointShowTime*2,...
                    1+length(Baseline_timeRange)+pointShowTime*3,...
                    1+length(Baseline_timeRange)+pointShowTime*4,...
                    1+length(Baseline_timeRange)+pointShowTime*4+delay1_timeDuration])
            end
            
            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 8)
            set(gca,'box','off');
            
            xlim([temp_x_min,temp_x_max]);
            
            xtickangle(0);
            
            
            %ytickk = [temp_y_min:0.1:0,0.1:0.1:temp_y_max];
            %set(gca,'YTick',ytickk)
            if temp_length == 1
                if if_errorBarSession_raw0_zscore1 == 0
                    ylabel('Pupil size (mm)','FontSize',9);
                elseif if_errorBarSession_raw0_zscore1 == 1
                    ylabel('Pupil size (z-scored)','FontSize',9);
                end
            end
            %set(gca,'yticklabels',string(ytickk))
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1,temp_y_max+(temp_y_max-temp_y_min)*0.1]);
            
            title(['Length',num2str(temp_length)],'Fontsize',9);
            
        end
        
        
    end
    
    
    
    if true
        fig = figure('Name',sprintf('fig'),'NumberTitle','off');
        
        %set(gcf,'Position',[500 450 389.2*0.65*1.04 176.7*0.94*1.06*0.995]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[500 450 389.2*0.65*1.04 176.7*0.94*1.06*0.995*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        %         if if_plot_errorBar_trial0_session1 == 0
        %             t.Title.String = sprintf('Monkey %s pupil size\n%s, error bar from trials',monkey_name,pupilSize_fileName);
        %         elseif if_plot_errorBar_trial0_session1 == 1
        %             t.Title.String = sprintf('Monkey %s pupil size\n%s, error bar from sessions',monkey_name,pupilSize_fileName);
        %         end
        %
        %         t.Title.FontSize = 11;
        %         t.Title.Interpreter = 'none';
        
        
        %for temp_length = 1:pointKindsNum
        if true
            temp_length = 2;
            
            nexttile
            
            tempMeanG_pupilSize = [meanG_pupilSizeBaseline{temp_length},meanG_pupilSizeSample{temp_length},meanG_pupilSizeDelay1{temp_length}];
            tempMeanR_pupilSize = [meanR_pupilSizeBaseline{temp_length},meanR_pupilSizeSample{temp_length},meanR_pupilSizeDelay1{temp_length}];
            tempSEMG_pupilSize = [semG_pupilSizeBaseline{temp_length},semG_pupilSizeSample{temp_length},semG_pupilSizeDelay1{temp_length}];
            tempSEMR_pupilSize = [semR_pupilSizeBaseline{temp_length},semR_pupilSizeSample{temp_length},semR_pupilSizeDelay1{temp_length}];
            
            x = 1:length(tempMeanG_pupilSize);
            
            [temp_x_min,temp_x_max] = bounds(x);
            
            h_line = [];
            
            x = x; %#ok<*ASGSL>
            y = tempMeanG_pupilSize;
            h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');];
            hold on
            y_sem = tempSEMG_pupilSize;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
            hold on
            
            
            x = x;
            y = tempMeanR_pupilSize;
            h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');];
            hold on
            y_sem = tempSEMR_pupilSize;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
            hold on
            
                        
            for tempi=1:temp_length
                temp_x = [1+length(Baseline_timeRange)+pointShowTime*(tempi-1) temp_y_min-(temp_y_max-temp_y_min)*0.1;...
                    1+length(Baseline_timeRange)+pointShowTime*(tempi-1) temp_y_max+(temp_y_max-temp_y_min)*0.1;...
                    1+length(Baseline_timeRange)+pointShowTime*(tempi-1)+200 temp_y_max+(temp_y_max-temp_y_min)*0.1;...
                    1+length(Baseline_timeRange)+pointShowTime*(tempi-1)+200 temp_y_min-(temp_y_max-temp_y_min)*0.1];
                temp_y = [1 2 3 4];
                patch('Faces',temp_y,'Vertices',temp_x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.25,'EdgeColor','none');%0.1
                hold on
            end
            
            
            temp_x = [1 temp_y_min-(temp_y_max-temp_y_min)*0.1;...
                1 temp_y_max+(temp_y_max-temp_y_min)*0.1;...
                1+length(Baseline_timeRange) temp_y_max+(temp_y_max-temp_y_min)*0.1;...
                1+length(Baseline_timeRange) temp_y_min-(temp_y_max-temp_y_min)*0.1];
            temp_y = [1 2 3 4];
            patch('Faces',temp_y,'Vertices',temp_x,'FaceColor',[217,217,162]/255,'FaceAlpha',0.4,'EdgeColor','none');%0.1
            hold on
            
            
            if temp_length == 2
                if if_choiceMemory0_choiceMemoryCorrect1 == 0
                    %le = legend(h_line,'Choice-memory','Choice-offload','Location','northeast','fontsize',6.5);
                    le = legend(h_line,'Memory','Offload','Location','northeast','fontsize',6.5);
                elseif if_choiceMemory0_choiceMemoryCorrect1 == 1
                    le = legend(h_line,'Memory correct','Offload','Location','northeast','fontsize',7);
                end
                le.ItemTokenSize = ones(1,3)*10;
                legend('boxoff');
            end
            
            if temp_length == 1
                set(gca,'xticklabels',{'Fixation','T1','Delay-on'})
                set(gca,'XTick',[1,1+length(Baseline_timeRange),...
                    1+length(Baseline_timeRange)+pointShowTime*1,...
                    1+length(Baseline_timeRange)+pointShowTime*1+delay1_timeDuration])
            elseif temp_length == 2
                set(gca,'xticklabels',{'Fixation','T1','T2','Delay-on'})
                set(gca,'XTick',[1,1+length(Baseline_timeRange),...
                    1+length(Baseline_timeRange)+pointShowTime*1,...
                    1+length(Baseline_timeRange)+pointShowTime*2,...
                    1+length(Baseline_timeRange)+pointShowTime*2+delay1_timeDuration])
            elseif temp_length == 3
                set(gca,'xticklabels',{'Fixation','T1','T2','T3','Delay-on'})
                set(gca,'XTick',[1,1+length(Baseline_timeRange),...
                    1+length(Baseline_timeRange)+pointShowTime*1,...
                    1+length(Baseline_timeRange)+pointShowTime*2,...
                    1+length(Baseline_timeRange)+pointShowTime*3,...
                    1+length(Baseline_timeRange)+pointShowTime*3+delay1_timeDuration])
            elseif temp_length == 4
                set(gca,'xticklabels',{'Fixation','T1','T2','T3','T4','Delay-on'})
                set(gca,'XTick',[1,1+length(Baseline_timeRange),...
                    1+length(Baseline_timeRange)+pointShowTime*1,...
                    1+length(Baseline_timeRange)+pointShowTime*2,...
                    1+length(Baseline_timeRange)+pointShowTime*3,...
                    1+length(Baseline_timeRange)+pointShowTime*4,...
                    1+length(Baseline_timeRange)+pointShowTime*4+delay1_timeDuration])
            end
            
            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 8)
            set(gca,'box','off');
            
            xlim([temp_x_min,temp_x_max]);
            
            xtickangle(0);
            
            
            %ytickk = [temp_y_min:0.1:0,0.1:0.1:temp_y_max];
            %set(gca,'YTick',ytickk)
            %if temp_length == 1
            if if_errorBarSession_raw0_zscore1 == 0
                ylabel('Pupil size (mm)','FontSize',9);
            elseif if_errorBarSession_raw0_zscore1 == 1
                ylabel('Pupil size (z-scored)','FontSize',9);
            end
            %end
            %set(gca,'yticklabels',string(ytickk))
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1,temp_y_max+(temp_y_max-temp_y_min)*0.1]);
            
            title(['Length',num2str(temp_length)],'Fontsize',9);
            %title(sprintf(''),'Fontsize',9);
            
        end
        
        
    end    
    
    
    
    
    if true
        fig = figure('Name',sprintf('fig'),'NumberTitle','off');
        
        %set(gcf,'Position',[800 450 389.2*0.65*1.04 176.7*0.94*1.06*0.995*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[800 450 389.2*0.65*1.04*0.68 176.7*0.94*1.06*0.995*0.8*1.28]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        if true
            temp_length = 2;
            
            nexttile
            
            tempMeanG_pupilSize = [meanG_pupilSizeBaseline{temp_length},meanG_pupilSizeSample{temp_length},meanG_pupilSizeDelay1{temp_length}];
            tempMeanR_pupilSize = [meanR_pupilSizeBaseline{temp_length},meanR_pupilSizeSample{temp_length},meanR_pupilSizeDelay1{temp_length}];
            tempSEMG_pupilSize = [semG_pupilSizeBaseline{temp_length},semG_pupilSizeSample{temp_length},semG_pupilSizeDelay1{temp_length}];
            tempSEMR_pupilSize = [semR_pupilSizeBaseline{temp_length},semR_pupilSizeSample{temp_length},semR_pupilSizeDelay1{temp_length}];
            
            x = 1:length(tempMeanG_pupilSize);
            
            [temp_x_min,temp_x_max] = bounds(x);
            
            h_line = [];
            
            x = x; %#ok<*ASGSL>
            y = tempMeanG_pupilSize;
            h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');];
            hold on
            y_sem = tempSEMG_pupilSize;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
            hold on
            
            
            x = x;
            y = tempMeanR_pupilSize;
            h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');];
            hold on
            y_sem = tempSEMR_pupilSize;
            patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
            hold on
            
                        
            for tempi=1:temp_length
                temp_x = [1+length(Baseline_timeRange)+pointShowTime*(tempi-1) temp_y_min-(temp_y_max-temp_y_min)*0.1;...
                    1+length(Baseline_timeRange)+pointShowTime*(tempi-1) temp_y_max+(temp_y_max-temp_y_min)*0.1;...
                    1+length(Baseline_timeRange)+pointShowTime*(tempi-1)+200 temp_y_max+(temp_y_max-temp_y_min)*0.1;...
                    1+length(Baseline_timeRange)+pointShowTime*(tempi-1)+200 temp_y_min-(temp_y_max-temp_y_min)*0.1];
                temp_y = [1 2 3 4];
                patch('Faces',temp_y,'Vertices',temp_x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.25,'EdgeColor','none');%0.1
                hold on
            end
            
            
            temp_x = [1 temp_y_min-(temp_y_max-temp_y_min)*0.1;...
                1 temp_y_max+(temp_y_max-temp_y_min)*0.1;...
                1+length(Baseline_timeRange) temp_y_max+(temp_y_max-temp_y_min)*0.1;...
                1+length(Baseline_timeRange) temp_y_min-(temp_y_max-temp_y_min)*0.1];
            temp_y = [1 2 3 4];
            patch('Faces',temp_y,'Vertices',temp_x,'FaceColor',[217,217,162]/255,'FaceAlpha',0.4,'EdgeColor','none');%0.1
            hold on
            
            
            if temp_length == 2
                if if_choiceMemory0_choiceMemoryCorrect1 == 0
                    %le = legend(h_line,'Choice-memory','Choice-offload','Location','northeast','fontsize',6.5);
                    le = legend(h_line,'Memory','Offload','Location','northeast','fontsize',6.5);
                elseif if_choiceMemory0_choiceMemoryCorrect1 == 1
                    le = legend(h_line,'Memory correct','Offload','Location','northeast','fontsize',7);
                end
                le.ItemTokenSize = ones(1,3)*10;
                legend('boxoff');
            end
            
            if temp_length == 1
                set(gca,'xticklabels',{'Fixation','T1','Delay-on'})
                set(gca,'XTick',[1,1+length(Baseline_timeRange),...
                    1+length(Baseline_timeRange)+pointShowTime*1,...
                    1+length(Baseline_timeRange)+pointShowTime*1+delay1_timeDuration])
            elseif temp_length == 2
                set(gca,'xticklabels',{'Fixation','T1','T2','Delay-on'})
                set(gca,'XTick',[1,1+length(Baseline_timeRange),...
                    1+length(Baseline_timeRange)+pointShowTime*1,...
                    1+length(Baseline_timeRange)+pointShowTime*2,...
                    1+length(Baseline_timeRange)+pointShowTime*2+delay1_timeDuration])
            elseif temp_length == 3
                set(gca,'xticklabels',{'Fixation','T1','T2','T3','Delay-on'})
                set(gca,'XTick',[1,1+length(Baseline_timeRange),...
                    1+length(Baseline_timeRange)+pointShowTime*1,...
                    1+length(Baseline_timeRange)+pointShowTime*2,...
                    1+length(Baseline_timeRange)+pointShowTime*3,...
                    1+length(Baseline_timeRange)+pointShowTime*3+delay1_timeDuration])
            elseif temp_length == 4
                set(gca,'xticklabels',{'Fixation','T1','T2','T3','T4','Delay-on'})
                set(gca,'XTick',[1,1+length(Baseline_timeRange),...
                    1+length(Baseline_timeRange)+pointShowTime*1,...
                    1+length(Baseline_timeRange)+pointShowTime*2,...
                    1+length(Baseline_timeRange)+pointShowTime*3,...
                    1+length(Baseline_timeRange)+pointShowTime*4,...
                    1+length(Baseline_timeRange)+pointShowTime*4+delay1_timeDuration])
            end
            
            
            text(575,temp_y_max-(temp_y_max-temp_y_min)*0.08,'***','Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            

            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 8)
            set(gca,'box','off');
            
            xlim([temp_x_min,temp_x_max]);
            
            %xtickangle(0);
            
            
            %ytickk = [temp_y_min:0.1:0,0.1:0.1:temp_y_max];
            %set(gca,'YTick',ytickk)
            %if temp_length == 1
            if if_errorBarSession_raw0_zscore1 == 0
                ylabel('Pupil size (mm)','FontSize',9);
            elseif if_errorBarSession_raw0_zscore1 == 1
                ylabel('Pupil size (z-scored)','FontSize',9);
            end
            %end
            %set(gca,'yticklabels',string(ytickk))
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1,temp_y_max+(temp_y_max-temp_y_min)*0.1]);
            
            title(['Length',num2str(temp_length)],'Fontsize',9);
            %title(sprintf(''),'Fontsize',9);
            
        end
        
        
    end    
    
    
end


fprintf('Pupil size analysis (baseline) time = %.1f secs.\n',toc(t0));


%% End