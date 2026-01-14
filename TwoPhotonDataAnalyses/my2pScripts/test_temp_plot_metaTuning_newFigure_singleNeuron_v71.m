% Chuan's 24th script (20251214)
% This script: To plot single-neuron's meta-memory selectivity, related to figure 3.
%% Initialization

if_plot = 1;

if_plot_multiPanel = 0;%1
if_plot_singlePanel = 0;
if_plot_multiPanelB = 0;
if_plot_multiPanelC = 1;% Baseline
if_plot_multiPanelD = 0;% Delay1
if_plot_multiPanelD2 = 0;% Delay1, plot mean of length 1 2 3, from lastT (for Fig. 1)
if_plot_multiPanelE = 0;% Delay1 (for Fig. 3 & 4)
if_plot_rProb = 0;%1
if_plot_rProb_B = 0;


if_delay1_forward0_backward1 = 0; % new

if if_plot_multiPanelD2 == 1
    if_delay1_forward0_backward1 = 0;
end

if if_plot_multiPanelE == 1
    if_delay1_forward0_backward1 = 1;
end

if if_plot_multiPanelC == 1
    if_delay1_forward0_backward1 = 1;
end

% if currentSessionIndex_AB ~= 8
%     if_delay1_forward0_backward1 = 1;
% end

if_colormap_loadEnhanced = 0;% only in multiPanelB (baseline)

if_sort_rProb0_precision1 = 0;

if_halfTrial = 1;%0

tempBoolIndex_validTrial = true(1,trial_num);
if if_halfTrial == 1
    if if_monkey_D0_Z1 == 0
        tempBoolIndex_validTrial(1:round(trial_num/2)) = false;
        %tempBoolIndex_validTrial(round(trial_num/2)+1:end) = false;
    elseif if_monkey_D0_Z1 == 1
        %tempBoolIndex_validTrial(1:round(trial_num/2)) = false;
        tempBoolIndex_validTrial(round(trial_num/2)+1:end) = false;
    end
    
    %tempBoolIndex_validTrial(416:end) = false;
end


% low location, high meta
% cellIndex_suite2p_temptempRaw = 724;%
% cellIndex_suite2p_temptempRaw = 695;% good
% cellIndex_suite2p_temptempRaw = 1005;%
% cellIndex_suite2p_temptempRaw = 323;%
% cellIndex_suite2p_temptempRaw = 454;%
% cellIndex_suite2p_temptempRaw = 816;% good
% cellIndex_suite2p_temptempRaw = 847;%

% high location, low meta
% cellIndex_suite2p_temptempRaw = 2;% good
% cellIndex_suite2p_temptempRaw = 11;%
% cellIndex_suite2p_temptempRaw = 197;% (0510A)
% cellIndex_suite2p_temptempRaw = 50;% (0614A)
% cellIndex_suite2p_temptempRaw = 17;% (0312A_Z)
% cellIndex_suite2p_temptempRaw = 19;% (0312A_Z)
% cellIndex_suite2p_temptempRaw = 101;% (0312A_Z)
% cellIndex_suite2p_temptempRaw = 29;% (0124A_Z)

% high location, high meta
% cellIndex_suite2p_temptempRaw = 25;%
% cellIndex_suite2p_temptempRaw = 20;% good
% cellIndex_suite2p_temptempRaw = 124;% (0510A) good
% cellIndex_suite2p_temptempRaw = 1444;% (0510A) 

% suite2p cellIndex in sessionB
% offloading rate tuning neuron
% cellIndex_suite2p_temptempRaw = 92;% good
% cellIndex_suite2p_temptempRaw = 791;%
% cellIndex_suite2p_temptempRaw = 20;%
% cellIndex_suite2p_temptempRaw = 736;%
% cellIndex_suite2p_temptempRaw = 126;% good
% cellIndex_suite2p_temptempRaw = 356;%
% cellIndex_suite2p_temptempRaw = 896;%
% cellIndex_suite2p_temptempRaw = 1338;%
% cellIndex_suite2p_temptempRaw = 99;%
% cellIndex_suite2p_temptempRaw = 695;%
% cellIndex_suite2p_temptempRaw = 1207;%
% cellIndex_suite2p_temptempRaw = 1;% (0510A) good
% cellIndex_suite2p_temptempRaw = 78;% (0510A)

% precision tuning neuron
% cellIndex_suite2p_temptempRaw = 410;%
% cellIndex_suite2p_temptempRaw = 656;%(0513A)
% cellIndex_suite2p_temptempRaw = 21;%(0509A)

% baseline tuning neuron
% cellIndex_suite2p_temptempRaw = 1;%
% cellIndex_suite2p_temptempRaw = 18;%
% cellIndex_suite2p_temptempRaw = 5;% Ding, FOV8
% cellIndex_suite2p_temptempRaw = 1216;% Zelku, FOV9
cellIndex_suite2p_temptempRaw = 364;% Zelku, FOV9
% cellIndex_suite2p_temptempRaw = 421;% Zelku, FOV9
% cellIndex_suite2p_temptempRaw = 913;% Zelku, FOV9
% cellIndex_suite2p_temptempRaw = 1075;% Zelku, FOV9
% cellIndex_suite2p_temptempRaw = 1020;% Zelku, FOV9


% pure location tuning neuron (no precision tuning)
% cellIndex_suite2p_temptempRaw = 19;%

% mix tuning for location and precision neuron
% cellIndex_suite2p_temptempRaw = 428;%(0513A)
% cellIndex_suite2p_temptempRaw = 0;%(0513A)
% cellIndex_suite2p_temptempRaw = 27;%(0526A)
% cellIndex_suite2p_temptempRaw = 5;%(0513A)


% pure offloading rate tuning neuron (no precision tuning)
% cellIndex_suite2p_temptempRaw = 78;%(0511A)
% cellIndex_suite2p_temptempRaw = 906;%
% cellIndex_suite2p_temptempRaw = 1;%(0511A)
% cellIndex_suite2p_temptempRaw = 20;%

% pure precision tuning neuron (no offloading rate)
% cellIndex_suite2p_temptempRaw = 265;%

% color_choiceMemory = [141,160,203]/255;%[0.4660 0.6740 0.1880], [179,205,227]/255
% color_choiceOffload = [252,141,98]/255;%[0.6350 0.0780 0.1840], [251,180,174]/255

color_choiceMemory = [1,133+20,113]/255;
color_choiceOffload = [166+40,97,26]/255;


temp_EdgeColor = 'none';%[0.62,0.62,0.62]-->'none'
temp_FaceAlpha = 0.35;%0.1-->0.05-->0.35



% temp_cellIndex_target = find(decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,2)==cellIndex_suite2p_temptempRaw);
temp_cellIndex_target = find(decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,end)==cellIndex_suite2p_temptempRaw);
cellIndex_suite2p_temptemp = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(temp_cellIndex_target,1);

xlim_padSize1 = 2;
periodSkipInterval = 2;%3
periodSkipIntervalX = 11;%10,12,8,11
periodSkipIntervalB = 0;%3
temp_ticks_baselinePeriod = baselinePeriod_interval;
temp_ticks_sampleLastT = temp_ticks_baselinePeriod(end)+periodSkipInterval+length1_sample_interval;
% temp_ticks_decisionPeriodA = temp_ticks_baselinePeriod(end)+periodSkipIntervalX+decisionPeriodA_interval;
temp_ticks_decisionPeriodA = temp_ticks_sampleLastT(end)+periodSkipIntervalB+decisionPeriodA_interval;
temp_ticks_decisionPeriod = temp_ticks_decisionPeriodA(end)+periodSkipInterval+decisionPeriod_interval;

% backgounrdColor = [1 1 1]*0.825;%0.875
backgounrdColor = [1 1 1];

temp_id2 = find(cellIndex_suite2p == cellIndex_suite2p_temptemp);

F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
cellID_F_dff_decisionBin1 = F_dff_decisionBin1(temp_id2,:);
cellID_F_dff_decisionBin1_seqMerged = zeros(1,sum(numSeq));
for tempi=1:sum(numSeq)
    temptempBoolIndex = seqIndex==tempi;
    
    if if_halfTrial == 1
        temptempBoolIndex(1:round(trial_num/2)) = false;
        %temptempBoolIndex(round(trial_num/2):end) = false;
        %temptempBoolIndex(416:end) = false;
    end
    
    temp_dff = cellID_F_dff_decisionBin1(temptempBoolIndex);
    cellID_F_dff_decisionBin1_seqMerged(tempi) = mean(temp_dff);
end

%% offloading rate tuning
% figglm
x = cellID_F_dff_decisionBin1_seqMerged';
y = offloadingProb_inOne';
temp_mdl = fitglm(x,y);
beta_rProb = temp_mdl.Coefficients.Estimate;
r2_rProb = temp_mdl.Rsquared.Adjusted;
p_rProb = temp_mdl.Coefficients.pValue(2);

fprintf('r2_rProb = %.3f, p_rProb = %.3f.\n',r2_rProb,p_rProb);

%% behavior precision tuning
% figglm
x = cellID_F_dff_decisionBin1_seqMerged';
y = seqPrecision_behavior_56';
temp_mdl = fitglm(x,y);
beta_precision = temp_mdl.Coefficients.Estimate;
r2_precision = temp_mdl.Rsquared.Adjusted;
p_precision = temp_mdl.Coefficients.pValue(2);

%%

tempBoolIndex_validLength = seqIndex<=sum(numSeq(1:3));

%%
if true
    temptempIndex_choiceMemory = find((choiceMemoryBoolIndex&tempBoolIndex_validTrial&tempBoolIndex_validLength)==true);
    temptempIndex_choiceOffload = find((choiceOffloadBoolIndex&tempBoolIndex_validTrial&tempBoolIndex_validLength)==true);
 
    [temptempIndex_start,temptempIndex_end] = bounds([temptempIndex_choiceMemory,temptempIndex_choiceOffload]);
    temptempIndex_middle = round((temptempIndex_start+temptempIndex_end)/2);
    
    temptempIndex_middle_choiceMemory = find(temptempIndex_choiceMemory>temptempIndex_middle,1');
    temptempIndex_middle_choiceOffload = find(temptempIndex_choiceOffload>temptempIndex_middle,1');    
end

%% baselinePeriod
cellID_F_dff_baselinePeriod_choiceMemory = squeeze(F_dff_baselinePeriod(temp_id2,choiceMemoryBoolIndex&tempBoolIndex_validTrial&tempBoolIndex_validLength,:));
cellID_F_dff_baselinePeriod_choiceOffload = squeeze(F_dff_baselinePeriod(temp_id2,choiceOffloadBoolIndex&tempBoolIndex_validTrial&tempBoolIndex_validLength,:));

cellID_F_dff_baselinePeriod_choiceMemory_mean = mean(cellID_F_dff_baselinePeriod_choiceMemory,1);
cellID_F_dff_baselinePeriod_choiceOffload_mean = mean(cellID_F_dff_baselinePeriod_choiceOffload,1);

cellID_F_dff_baselinePeriod_choiceMemory_sem = std(cellID_F_dff_baselinePeriod_choiceMemory,1)...
    ./sqrt(size(cellID_F_dff_baselinePeriod_choiceMemory,1));
cellID_F_dff_baselinePeriod_choiceOffload_sem = std(cellID_F_dff_baselinePeriod_choiceOffload,1)...
    ./sqrt(size(cellID_F_dff_baselinePeriod_choiceOffload,1));

%% sampleLastT
cellID_F_dff_length1_sample = squeeze(F_dff_length1_sample(temp_id2,:,:));
cellID_F_dff_length2_sample = squeeze(F_dff_length2_sample(temp_id2,:,:));
cellID_F_dff_length3_sample = squeeze(F_dff_length3_sample(temp_id2,:,:));

cellID_F_dff_lastT_sample = nan(size(cellID_F_dff_length1_sample));

for tempi=1:size(cellID_F_dff_lastT_sample,1)
    temp1 = cellID_F_dff_length1_sample(tempi,:);
    temp2 = cellID_F_dff_length2_sample(tempi,end-17:end);
    temp3 = cellID_F_dff_length3_sample(tempi,end-17:end);
    
    if ~isnan(temp1(1,1))
        cellID_F_dff_lastT_sample(tempi,:) = temp1;
    end
    if ~isnan(temp2(1,1))
        cellID_F_dff_lastT_sample(tempi,:) = temp2;
    end
    if ~isnan(temp3(1,1))
        cellID_F_dff_lastT_sample(tempi,:) = temp3;
    end
    
end

cellID_F_dff_sampleLastT_choiceMemory = cellID_F_dff_lastT_sample(choiceMemoryBoolIndex&tempBoolIndex_validTrial&tempBoolIndex_validLength,:);
cellID_F_dff_sampleLastT_choiceOffload = cellID_F_dff_lastT_sample(choiceOffloadBoolIndex&tempBoolIndex_validTrial&tempBoolIndex_validLength,:);

cellID_F_dff_sampleLastT_choiceMemory_mean = mean(cellID_F_dff_sampleLastT_choiceMemory,1,'omitnan');
cellID_F_dff_sampleLastT_choiceOffload_mean = mean(cellID_F_dff_sampleLastT_choiceOffload,1,'omitnan');

cellID_F_dff_sampleLastT_choiceMemory_sem = std(cellID_F_dff_sampleLastT_choiceMemory,1,'omitnan')...
    ./sqrt(size(cellID_F_dff_sampleLastT_choiceMemory,1));
cellID_F_dff_sampleLastT_choiceOffload_sem = std(cellID_F_dff_sampleLastT_choiceOffload,1,'omitnan')...
    ./sqrt(size(cellID_F_dff_sampleLastT_choiceOffload,1));


%% decisionPeriodA
if if_delay1_forward0_backward1 == 1
    cellID_F_dff_decisionPeriodA_choiceMemory = squeeze(F_dff_decisionPeriodA(temp_id2,choiceMemoryBoolIndex&tempBoolIndex_validTrial&tempBoolIndex_validLength,:));
    cellID_F_dff_decisionPeriodA_choiceOffload = squeeze(F_dff_decisionPeriodA(temp_id2,choiceOffloadBoolIndex&tempBoolIndex_validTrial&tempBoolIndex_validLength,:));
elseif if_delay1_forward0_backward1 == 0
    cellID_F_dff_decisionPeriodA_choiceMemory = squeeze(F_dff_decisionPeriodA(temp_id2,choiceMemoryBoolIndex&tempBoolIndex_validTrial&tempBoolIndex_validLength,:));
    cellID_F_dff_decisionPeriodA_choiceOffload = squeeze(F_dff_decisionPeriodA(temp_id2,choiceOffloadBoolIndex&tempBoolIndex_validTrial&tempBoolIndex_validLength,:));    
    
    cellID_F_dff_decisionPeriodA_choiceMemory(:,1:decisionPeriodB_interval(2)) = squeeze(F_dff_decisionPeriodB(temp_id2,choiceMemoryBoolIndex&tempBoolIndex_validTrial&tempBoolIndex_validLength,:));
    cellID_F_dff_decisionPeriodA_choiceOffload(:,1:decisionPeriodB_interval(2)) = squeeze(F_dff_decisionPeriodB(temp_id2,choiceOffloadBoolIndex&tempBoolIndex_validTrial&tempBoolIndex_validLength,:));    
end

cellID_F_dff_decisionPeriodA_choiceMemory_mean = mean(cellID_F_dff_decisionPeriodA_choiceMemory,1);
cellID_F_dff_decisionPeriodA_choiceOffload_mean = mean(cellID_F_dff_decisionPeriodA_choiceOffload,1);

cellID_F_dff_decisionPeriodA_choiceMemory_sem = std(cellID_F_dff_decisionPeriodA_choiceMemory,1)...
    ./sqrt(size(cellID_F_dff_decisionPeriodA_choiceMemory,1));
cellID_F_dff_decisionPeriodA_choiceOffload_sem = std(cellID_F_dff_decisionPeriodA_choiceOffload,1)...
    ./sqrt(size(cellID_F_dff_decisionPeriodA_choiceOffload,1));


%% decisionPeriod
cellID_F_dff_decisionPeriod_choiceMemory = squeeze(F_dff_decisionPeriod(temp_id2,choiceMemoryBoolIndex&tempBoolIndex_validTrial&tempBoolIndex_validLength,:));
cellID_F_dff_decisionPeriod_choiceOffload = squeeze(F_dff_decisionPeriod(temp_id2,choiceOffloadBoolIndex&tempBoolIndex_validTrial&tempBoolIndex_validLength,:));

cellID_F_dff_decisionPeriod_choiceMemory_mean = mean(cellID_F_dff_decisionPeriod_choiceMemory,1);
cellID_F_dff_decisionPeriod_choiceOffload_mean = mean(cellID_F_dff_decisionPeriod_choiceOffload,1);

cellID_F_dff_decisionPeriod_choiceMemory_sem = std(cellID_F_dff_decisionPeriod_choiceMemory,1)...
    ./sqrt(size(cellID_F_dff_decisionPeriod_choiceMemory,1));
cellID_F_dff_decisionPeriod_choiceOffload_sem = std(cellID_F_dff_decisionPeriod_choiceOffload,1)...
    ./sqrt(size(cellID_F_dff_decisionPeriod_choiceOffload,1));


% %% Others


%% cellID_F_dff_baselinePeriod_sortedSeq_collaped
F_dff_current = F_dff_baselinePeriod;

temp_dff = squeeze(F_dff_current(temp_id2,:,:));
temp_dff_seq = zeros(sum(numSeq),size(F_dff_current,3));
for tempi=1:sum(numSeq)
    temp_range = find(seqIndex == tempi &tempBoolIndex_validTrial);
    temp_dff2 = temp_dff(temp_range,:);
    temp_dff_seq(tempi,:) = mean(temp_dff2,1);
end
if if_sort_rProb0_precision1 == 0
    [temp_B,temp_I] = sort(offloadingProb_inOne);
elseif if_sort_rProb0_precision1 == 1
    [temp_B,temp_I] = sort(seqPrecision_behavior_56);
end
rProbAscendingSortedIndex_seq = temp_I;
temp_dff_sortedSeq = ...
    temp_dff_seq(rProbAscendingSortedIndex_seq,:);

cellID_F_dff_baselinePeriod_sortedSeq_collaped = temp_dff_sortedSeq;

%% cellID_F_dff_baselinePeriod_sortedSeq_collaped
F_dff_current = F_dff_decisionPeriodA;

temp_dff = squeeze(F_dff_current(temp_id2,:,:));
temp_dff_seq = zeros(sum(numSeq),size(F_dff_current,3));
for tempi=1:sum(numSeq)
    temp_range = find(seqIndex == tempi &tempBoolIndex_validTrial);
    temp_dff2 = temp_dff(temp_range,:);
    temp_dff_seq(tempi,:) = mean(temp_dff2,1);
end
if if_sort_rProb0_precision1 == 0
    [temp_B,temp_I] = sort(offloadingProb_inOne);
elseif if_sort_rProb0_precision1 == 1
    [temp_B,temp_I] = sort(seqPrecision_behavior_56);
end
rProbAscendingSortedIndex_seq = temp_I;
temp_dff_sortedSeq = ...
    temp_dff_seq(rProbAscendingSortedIndex_seq,:);

cellID_F_dff_decisionPeriodA_sortedSeq_collaped = temp_dff_sortedSeq;


%% cellID_F_dff_baselinePeriod_sortedSeq_collaped
F_dff_current = F_dff_decisionPeriod;

temp_dff = squeeze(F_dff_current(temp_id2,:,:));
temp_dff_seq = zeros(sum(numSeq),size(F_dff_current,3));
for tempi=1:sum(numSeq)
    temp_range = find(seqIndex == tempi &tempBoolIndex_validTrial);
    temp_dff2 = temp_dff(temp_range,:);
    temp_dff_seq(tempi,:) = mean(temp_dff2,1);
end
if if_sort_rProb0_precision1 == 0
    [temp_B,temp_I] = sort(offloadingProb_inOne);
elseif if_sort_rProb0_precision1 == 1
    [temp_B,temp_I] = sort(seqPrecision_behavior_56);
end
rProbAscendingSortedIndex_seq = temp_I;
temp_dff_sortedSeq = ...
    temp_dff_seq(rProbAscendingSortedIndex_seq,:);

cellID_F_dff_decisionPeriod_sortedSeq_collaped = temp_dff_sortedSeq;


%%

cellID_F_dff_baselinePeriod_choice_collaped = [cellID_F_dff_baselinePeriod_choiceMemory;cellID_F_dff_baselinePeriod_choiceOffload];
cellID_F_dff_sampleLastT_choice_collaped = [cellID_F_dff_sampleLastT_choiceMemory;cellID_F_dff_sampleLastT_choiceOffload];
cellID_F_dff_decisionPeriodA_choice_collaped = [cellID_F_dff_decisionPeriodA_choiceMemory;cellID_F_dff_decisionPeriodA_choiceOffload];
cellID_F_dff_decisionPeriod_choice_collaped = [cellID_F_dff_decisionPeriod_choiceMemory;cellID_F_dff_decisionPeriod_choiceOffload];

temp_choiceIndex = [1, size(cellID_F_dff_decisionPeriodA_choiceMemory,1)+1];



% min_baselineDecisionPeriod_raw = min([...
%     cellID_F_dff_baselinePeriod_choiceMemory_mean-cellID_F_dff_baselinePeriod_choiceMemory_sem,...
%     cellID_F_dff_decisionPeriodA_choiceOffload_mean-cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
% 
% max_baselineDecisionPeriod_raw = max([...
%     cellID_F_dff_baselinePeriod_choiceMemory_mean+cellID_F_dff_baselinePeriod_choiceMemory_sem,...
%     cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
% 
% min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
% %min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.4;
% max_baselineDecisionPeriod = max_baselineDecisionPeriod_raw+(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;


%% Plot
if if_plot == 1
    close all
    
    if if_plot_multiPanel == 1
        %fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptemp,currentSession(18:22)),'NumberTitle','off');
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 530 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 520*0.97 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 504*0.9 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 454*0.6 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 504*0.9 250]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 454*0.6 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 272*0.815 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 222*1.08*1.1 260*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 222*1.08*1.1 260*0.85*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 222*1.08*1.1*1.095 260*0.85*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        set(gcf,'Position',[35+0 42+0 222*1.08*1.1*1.095*0.885 260*0.85*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        
        %t = tiledlayout(2,1,'TileSpacing','compact','Padding','tight'); %#ok<*NASGU>
        %t = tiledlayout(2,1,'TileSpacing','loose','Padding','tight'); %#ok<*NASGU>
        t = tiledlayout(2,1,'TileSpacing','tight','Padding','loose'); %#ok<*NASGU>
        
        %t.Title.String = sprintf('suite2pBID%d,%s',...
        %    cellIndex_suite2p_temptempRaw,FOVName_currentFOV2);
        %t.Title.FontSize = 11;
        %t.Title.Interpreter = 'none';
        
        %% Plot choice tuning
        
        nexttile
        
%         min_baselineDecisionPeriod_raw = min([...
%             cellID_F_dff_baselinePeriod_choiceMemory_mean-cellID_F_dff_baselinePeriod_choiceMemory_sem,...
%             cellID_F_dff_baselinePeriod_choiceOffload_mean-cellID_F_dff_baselinePeriod_choiceOffload_sem,...            
%             cellID_F_dff_decisionPeriodA_choiceMemory_mean-cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
%             cellID_F_dff_decisionPeriodA_choiceOffload_mean-cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
%         
% %         max_baselineDecisionPeriod_raw = max([...
% %             cellID_F_dff_baselinePeriod_choiceMemory_mean+cellID_F_dff_baselinePeriod_choiceMemory_sem,...
% %             cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
%         max_baselineDecisionPeriod_raw = max([...
%             cellID_F_dff_baselinePeriod_choiceMemory_mean+cellID_F_dff_baselinePeriod_choiceMemory_sem,...
%             cellID_F_dff_baselinePeriod_choiceOffload_mean+cellID_F_dff_baselinePeriod_choiceOffload_sem,...            
%             cellID_F_dff_decisionPeriodA_choiceMemory_mean+cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
%             cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        

        min_baselineDecisionPeriod_raw = min([...
            cellID_F_dff_decisionPeriodA_choiceMemory_mean-cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
            cellID_F_dff_decisionPeriodA_choiceOffload_mean-cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        
        max_baselineDecisionPeriod_raw = max([...
            cellID_F_dff_decisionPeriodA_choiceMemory_mean+cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
            cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);

        min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        %min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.4;
        max_baselineDecisionPeriod = max_baselineDecisionPeriod_raw+(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.15;
        
        h_line = [];
        
        % plot baselinePeriod
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceMemory_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceOffload_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        % plot decisionPeriodA
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceMemory_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceOffload_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        % plot decisionPeriod
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceMemory_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceOffload_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        % Others
        
        for tempi=1:length(baselinePeriod_interval)
            plot(temp_ticks_baselinePeriod(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        for tempi=1:length(decisionPeriodA_interval)
            plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        for tempi=1:length(decisionPeriod_interval)
            plot(temp_ticks_decisionPeriod(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        %le = legend(h_line,'memory','offload','Location','northeast','fontsize',8);%10
        %le = legend(h_line,'ChoiceMemory','ChoiceOffload','Location','northeast','fontsize',8);%10        
        %le = legend(h_line,'Choice-memory','Choice-offload','Location','northeast','fontsize',8);%10                
        
        %         le = legend(h_line,'Choice-memory','Choice-offload','Location','northeast','fontsize',8,'NumColumns',2);%10
        %         le.ItemTokenSize = ones(1,3)*10;
        %         le.Color = backgounrdColor;
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriod(end)+xlim_padSize2]);
        % xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)+1]);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]);
        ylim([min_baselineDecisionPeriod,max_baselineDecisionPeriod]);
        
        set(gca,'xticklabel',[])
        set(gca,'xtick',[])
        
        %xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        %xticklabels({'PreFixation','Fixation','','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'','','','Delay1','ChoiceCue','','','',''});
        
        %xtickangle(0);
        
        %text(temp_ticks_baselinePeriod(end),min_baselineDecisionPeriod,'/','fontsize',12);
        %text(temp_ticks_decisionPeriodA(1),min_baselineDecisionPeriod,'/','fontsize',12);
        %
        %     text(temp_ticks_decisionPeriodA(end),min_baselineDecisionPeriod,'/','fontsize',15);
        %     text(temp_ticks_decisionPeriod(1),min_baselineDecisionPeriod,'/','fontsize',15);
        
        yticks([0 0.2]);
        
        set(gca, 'FontSize', 10)
        set(gca,'box','off');% 取消右、上边框
        ylabel(sprintf('dF/F'), 'FontSize', 10, 'FontWeight', 'normal');
        
        
        %% Plot offloading rate tuning
        
        nexttile
        
        cmin = min(cellID_F_dff_decisionPeriodA_sortedSeq_collaped,[],'all');
        cmax = max(cellID_F_dff_decisionPeriodA_sortedSeq_collaped,[],'all');
        
        
        x = [temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end)];
        y = [1,size(cellID_F_dff_baselinePeriod_sortedSeq_collaped,1)];
        C = cellID_F_dff_baselinePeriod_sortedSeq_collaped;
        imagesc(x,y,C,[cmin cmax]);
        hold on
        
        x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)];
        y = [1,size(cellID_F_dff_decisionPeriodA_sortedSeq_collaped,1)];
        C = cellID_F_dff_decisionPeriodA_sortedSeq_collaped;
        imagesc(x,y,C,[cmin cmax]);
        hold on
        
        x = [temp_ticks_decisionPeriod(1) temp_ticks_decisionPeriod(end)];
        y = [1,size(cellID_F_dff_decisionPeriod_sortedSeq_collaped,1)];
        C = cellID_F_dff_decisionPeriod_sortedSeq_collaped;
        imagesc(x,y,C,[cmin cmax]);
        hold on
        
        
        for tempi=1:length(baselinePeriod_interval)
            plot(temp_ticks_baselinePeriod(tempi)*[1 1],[1 size(cellID_F_dff_baselinePeriod_sortedSeq_collaped,1)],...
                '-','LineWidth',1,'Color',[1 1 1])
            hold on
        end
        
        for tempi=1:length(decisionPeriodA_interval)
            plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[1 size(cellID_F_dff_decisionPeriodA_sortedSeq_collaped,1)],...
                '-','LineWidth',1,'Color',[1 1 1]);
            hold on
        end
        
        for tempi=1:length(decisionPeriod_interval)
            plot(temp_ticks_decisionPeriod(tempi)*[1 1],[1 size(cellID_F_dff_decisionPeriod_sortedSeq_collaped,1)],...
                '-','LineWidth',1,'Color',[1 1 1]);
            hold on
        end
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        % xlim([1-xlim_padSize1 temp_ticks_decisionPeriod(end)+xlim_padSize2]);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]);
        ylim([1-1 size(cellID_F_dff_decisionPeriod_sortedSeq_collaped,1)+1]);
        set(gca,'xticklabel',[])
        set(gca,'xtick',[])
        
        %xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        %xticklabels({'PreFixation','Fixation','','','Decision',''});
        %xticklabels({'PreFixation','Fixation','','','ChoiceCue','','','Decision',''});
        %xticklabels('');
        
        %text(temp_ticks_baselinePeriod(end),size(cellID_F_dff_baselinePeriod_sortedSeq_collaped,1),'/','fontsize',12);
        %text(temp_ticks_decisionPeriodA(1),size(cellID_F_dff_baselinePeriod_sortedSeq_collaped,1),'/','fontsize',12);
        
        %text(temp_ticks_decisionPeriodA(end),size(cellID_F_dff_baselinePeriod_sortedSeq_collaped,1),'/','fontsize',12);
        %text(temp_ticks_decisionPeriod(1),size(cellID_F_dff_baselinePeriod_sortedSeq_collaped,1),'/','fontsize',12);
        
        xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        xticklabels({'PreFixation','Fixation','','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'','','','Delay1','ChoiceCue','','','',''});
        
        xtickangle(0);
        
        xtickangle(0);
        set(gca, 'FontSize', 10)
        set(gca,'box','off');% 取消右、上边框
        %ylabel(sprintf('sorted seqs'), 'FontSize', 12, 'FontWeight', 'bold');
        if if_sort_rProb0_precision1 == 0
            ylabel(sprintf('Sort seqs as\nOffloading rate'), 'FontSize', 10, 'FontWeight', 'normal');
        elseif if_sort_rProb0_precision1 == 1
            ylabel(sprintf('Sort seqs as\nPrecision'), 'FontSize', 10, 'FontWeight', 'normal');
        end
        
        c = colorbar('fontsize',9);  
        c.Ticks = [0 0.4];
        
    end
    
    
    if if_plot_multiPanel == 1
        %fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptemp,currentSession(18:22)),'NumberTitle','off');
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 530 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 520*0.97 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 504*0.9 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 454*0.6 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 504*0.9 250]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 454*0.6 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 272*0.815 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 222*1.08*1.1 260*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 222*1.08*1.1 260*0.85*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 222*1.08*1.1*1.095 260*0.85*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 222*1.08*1.1*1.0165*0.874*0.75*0.975 120*1.2*1.15*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 222*1.08*1.1*1.0165*0.874*0.75*0.975*1.83*0.98 120*1.2*1.15*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        set(gcf,'Position',[35+0 42+0 222*1.08*1.1*1.0165*0.874*0.75*0.975*1.83*0.98*0.8*0.9*0.85*0.95*0.975*0.99 120*1.2*1.15*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        
        t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact'); %#ok<*NASGU>

        
        %t.Title.String = sprintf('suite2pBID%d,%s',...
        %    cellIndex_suite2p_temptempRaw,FOVName_currentFOV2);
        %t.Title.FontSize = 11;
        %t.Title.Interpreter = 'none';
        
        %% Plot choice tuning
        
        nexttile
        
%         min_baselineDecisionPeriod_raw = min([...
%             cellID_F_dff_baselinePeriod_choiceMemory_mean-cellID_F_dff_baselinePeriod_choiceMemory_sem,...
%             cellID_F_dff_baselinePeriod_choiceOffload_mean-cellID_F_dff_baselinePeriod_choiceOffload_sem,...            
%             cellID_F_dff_decisionPeriodA_choiceMemory_mean-cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
%             cellID_F_dff_decisionPeriodA_choiceOffload_mean-cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
%         
% %         max_baselineDecisionPeriod_raw = max([...
% %             cellID_F_dff_baselinePeriod_choiceMemory_mean+cellID_F_dff_baselinePeriod_choiceMemory_sem,...
% %             cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
%         max_baselineDecisionPeriod_raw = max([...
%             cellID_F_dff_baselinePeriod_choiceMemory_mean+cellID_F_dff_baselinePeriod_choiceMemory_sem,...
%             cellID_F_dff_baselinePeriod_choiceOffload_mean+cellID_F_dff_baselinePeriod_choiceOffload_sem,...            
%             cellID_F_dff_decisionPeriodA_choiceMemory_mean+cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
%             cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        

        min_baselineDecisionPeriod_raw = min([...
            cellID_F_dff_decisionPeriodA_choiceMemory_mean-cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
            cellID_F_dff_decisionPeriodA_choiceOffload_mean-cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        
        max_baselineDecisionPeriod_raw = max([...
            cellID_F_dff_decisionPeriodA_choiceMemory_mean+cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
            cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);

        min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        %min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.4;
        max_baselineDecisionPeriod = max_baselineDecisionPeriod_raw+(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.15;
        
        h_line = [];
        
        % plot baselinePeriod
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceMemory_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceOffload_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        % plot decisionPeriodA
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceMemory_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceOffload_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        % plot decisionPeriod
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceMemory_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceOffload_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        % Others
        
        for tempi=1:length(baselinePeriod_interval)
            plot(temp_ticks_baselinePeriod(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        %for tempi=1:length(decisionPeriodA_interval)
        for tempi=1:(length(decisionPeriodA_interval)-1)
            plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        for tempi=1:length(decisionPeriod_interval)
            plot(temp_ticks_decisionPeriod(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        %le = legend(h_line,'memory','offload','Location','northeast','fontsize',8);%10
        %le = legend(h_line,'ChoiceMemory','ChoiceOffload','Location','northeast','fontsize',8);%10        
        %le = legend(h_line,'Choice-memory','Choice-offload','Location','northeast','fontsize',8);%10                
        
        le = legend(h_line,'Choice-memory','Choice-offload','Location','northeast','fontsize',8,'NumColumns',1);%10
        le.ItemTokenSize = ones(1,3)*10;
        le.Color = backgounrdColor;
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriod(end)+xlim_padSize2]);
        % xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)+1]);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]);
        ylim([min_baselineDecisionPeriod,max_baselineDecisionPeriod]);
        
        set(gca,'xticklabel',[])
        set(gca,'xtick',[])
        
        %xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        %xticklabels({'PreFixation','Fixation','','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'','','','Delay1','ChoiceCue','','','',''});
        
        xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        xticklabels({'PreFixation','Fixation','','Delay1','ChoiceCue','','','Decision',''});
        
        xtickangle(0);
        
        %text(temp_ticks_baselinePeriod(end),min_baselineDecisionPeriod,'/','fontsize',12);
        %text(temp_ticks_decisionPeriodA(1),min_baselineDecisionPeriod,'/','fontsize',12);
        %
        %     text(temp_ticks_decisionPeriodA(end),min_baselineDecisionPeriod,'/','fontsize',15);
        %     text(temp_ticks_decisionPeriod(1),min_baselineDecisionPeriod,'/','fontsize',15);
        
        
        yticks([0 0.2]);
        
        set(gca, 'FontSize', 8);%10
        set(gca,'box','off');% 取消右、上边框
        ylabel(sprintf('dF/F'), 'FontSize', 9, 'FontWeight', 'normal');%10
                
    end
    
    
    %if if_plot_multiPanel == 1 && if_sort_rProb0_precision1 == 1
    if if_plot_multiPanel == 1
        %fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptemp,currentSession(18:22)),'NumberTitle','off');
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 530 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 520*0.97 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 504*0.9 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 454*0.6 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 504*0.9 250]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 454*0.6 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 272*0.815 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 222*1.08*1.1 260*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 263*0.8*0.90 120]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 263*0.8*0.90 120*1.2*1.15]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 222*1.08*1.1*1.0165 120*1.2*1.15]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 222*1.08*1.1*1.0165*0.874*0.92 120*1.2*1.15]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 222*1.08*1.1*1.0165*0.874*0.92*0.92 120*1.2*1.15]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 222*1.08*1.1*1.0165*0.874*0.92*0.92*0.97*0.98*0.992*0.995 120*1.2*1.15*0.97*0.992]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 222*1.08*1.1*1.0165*0.874*0.92*0.92*0.97*0.98*0.992*0.995 120*1.2*1.15*0.97*0.992*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        set(gcf,'Position',[35+0 42+0 222*1.08*1.1*1.0165*0.874*0.92*0.92*0.97*0.98*0.992*0.995 120*1.2*1.15*0.97*0.992*0.9*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        
        t = tiledlayout(1,4,'TileSpacing','compact','Padding','tight'); %#ok<*NASGU>

        
                
        %% Plot offloading rate tuning        
        nexttile([1 3])
        
        cmin = min(cellID_F_dff_decisionPeriodA_sortedSeq_collaped,[],'all');
        cmax = max(cellID_F_dff_decisionPeriodA_sortedSeq_collaped,[],'all');
        
        cmin = 0;
        
        if cellIndex_suite2p_temptempRaw == 724
            %cmax = 0.6;
        end
        if cellIndex_suite2p_temptempRaw == 695
            cmin = 0;
            cmax = 0.15;
        end
        if cellIndex_suite2p_temptempRaw == 1005
            cmin = 0;
            cmax = 0.4;
        end
        if cellIndex_suite2p_temptempRaw == 92
            cmin = 0;
            cmax = 0.4;
        end    
        if cellIndex_suite2p_temptempRaw == 20
            cmin = 0.1;
            cmax = 1.5;
        end         
        if cellIndex_suite2p_temptempRaw == 2
            cmin = 0;
            cmax = 4;
        end       
        if cellIndex_suite2p_temptempRaw == 816
            cmin = 0.1;
            cmax = 0.3;
        end              
        
        if cellIndex_suite2p_temptempRaw == 124 % 0510A
            cmin = 0.2;
            cmax = 1.2;
        end
        
        if cellIndex_suite2p_temptempRaw == 1 % 0510A
            cmin = 0.025;
            cmax = 0.35;
        end        
        
        
        x = [temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end)];
        %y = [1,size(cellID_F_dff_baselinePeriod_sortedSeq_collaped,1)];
        y = [size(cellID_F_dff_baselinePeriod_sortedSeq_collaped,1),1];        
        C = cellID_F_dff_baselinePeriod_sortedSeq_collaped;
        imagesc(x,y,C,[cmin cmax]);
        hold on
        
        x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)];
        %y = [1,size(cellID_F_dff_decisionPeriodA_sortedSeq_collaped,1)];
        y = [size(cellID_F_dff_decisionPeriodA_sortedSeq_collaped,1),1];        
        C = cellID_F_dff_decisionPeriodA_sortedSeq_collaped;
        imagesc(x,y,C,[cmin cmax]);
        hold on
        
        x = [temp_ticks_decisionPeriod(1) temp_ticks_decisionPeriod(end)];
        %y = [1,size(cellID_F_dff_decisionPeriod_sortedSeq_collaped,1)];
        y = [size(cellID_F_dff_decisionPeriod_sortedSeq_collaped,1),1];        
        C = cellID_F_dff_decisionPeriod_sortedSeq_collaped;
        imagesc(x,y,C,[cmin cmax]);
        hold on
        
        
        for tempi=1:length(baselinePeriod_interval)
            plot(temp_ticks_baselinePeriod(tempi)*[1 1],[1 size(cellID_F_dff_baselinePeriod_sortedSeq_collaped,1)],...
                '-','LineWidth',1,'Color',[1 1 1])
            hold on
        end
        
        for tempi=1:length(decisionPeriodA_interval)
            plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[1 size(cellID_F_dff_decisionPeriodA_sortedSeq_collaped,1)],...
                '-','LineWidth',1,'Color',[1 1 1]);
            hold on
        end
        
        for tempi=1:length(decisionPeriod_interval)
            plot(temp_ticks_decisionPeriod(tempi)*[1 1],[1 size(cellID_F_dff_decisionPeriod_sortedSeq_collaped,1)],...
                '-','LineWidth',1,'Color',[1 1 1]);
            hold on
        end
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        % xlim([1-xlim_padSize1 temp_ticks_decisionPeriod(end)+xlim_padSize2]);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]);
        ylim([1-1 size(cellID_F_dff_decisionPeriod_sortedSeq_collaped,1)+1]);
        set(gca,'xticklabel',[])
        set(gca,'xtick',[])
                
        xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        %xticklabels({'PreFixation','Fixation','','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'PreFixation','Fixation','','Delay1','Decision','','','Decision',''});
        xticklabels({'PreFixation','Fixation','','Delay-on','Decision','','','Decision',''});
        
        xtickangle(0);
        
        xtickangle(0);
                
        set(gca,'ytick',[])
        %set(gca,'ytick',[0 50])
                        
        set(gca, 'FontSize', 8);%10
        set(gca,'box','off');% 取消右、上边框
        if if_sort_rProb0_precision1 == 0
            ylabel(sprintf('Sort seqs as\nOffloading rate'), 'FontSize', 9, 'FontWeight', 'normal');%10
        elseif if_sort_rProb0_precision1 == 1
            %ylabel(sprintf('Sort seqs as\nPrecision'), 'FontSize', 10, 'FontWeight', 'normal');
            ylabel(sprintf('Sort seqs as\nAccuracy'), 'FontSize', 9, 'FontWeight', 'normal');%10        
        end
        
        %colorbar('fontsize',8);        
        %c = colorbar('eastoutside','FontSize',8);
        %c.Position = c.Position+[0.04 0 -0.025 0];
        
        %c = colorbar('fontsize',8);%9
        c = colorbar('fontsize',6.5);%9
        %c.Ticks = [0 0.4];
        if cmax > 1
            c.Ticks = [cmin floor(cmax)];
        else
            c.Ticks = [cmin floor(cmax*10)/10];
        end
        
        %c.Ticks = [cmin cmax];
        
        if cellIndex_suite2p_temptempRaw == 1 % 0510A
            c.Ticks = [0.1 floor(cmax*10)/10];
        end
        
        %c.Position = c.Position+[0.05 0 -0.025 -0.03];
        c.Position = c.Position+[0.03 0 -0.025 -0.03];
        
        %colormap(parula_enhanced);
        %colormap parula
        %colormap(cividis());
        %colormap(viridis());
        %colormap jet
        
        %temp1 = pink;
        %temp1 = gray;
        %temp1 = temp1(end:-1:1,:);
        %colormap(temp1);

        %colormap hot
        %colormap cool
        
        colormap(coolwarm());
        
        %temp1 = coolwarm(300);
        %temp2 = ((300-256)/2)+1;
        %temp3 = temp1(temp2:temp2+255,:);
        %colormap(temp3);
        
    end
    
    
    %% Plot single panel
    if if_plot_singlePanel == 1
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 400*0.94 130*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[35+0 42+0 376*1.15 110*1.275]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
        
        t.Title.String = sprintf('suite2pBID%d,%s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2);
        t.Title.FontSize = 11;
        t.Title.Interpreter = 'none';
        
        nexttile
        
        
        min_baselineDecisionPeriod_raw = min([...
            cellID_F_dff_baselinePeriod_choiceMemory_mean-cellID_F_dff_baselinePeriod_choiceMemory_sem,...
            cellID_F_dff_decisionPeriodA_choiceOffload_mean-cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        
        max_baselineDecisionPeriod_raw = max([...
            cellID_F_dff_baselinePeriod_choiceMemory_mean+cellID_F_dff_baselinePeriod_choiceMemory_sem,...
            cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        
        min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        %min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.4;
        max_baselineDecisionPeriod = max_baselineDecisionPeriod_raw+(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        
        h_line = [];
        
        %% plot baselinePeriod
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceMemory_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceOffload_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        %% plot decisionPeriodA
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceMemory_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceOffload_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        %% plot decisionPeriod
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceMemory_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceOffload_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        %% Others
        
        for tempi=1:length(baselinePeriod_interval)
            plot(temp_ticks_baselinePeriod(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        for tempi=1:length(decisionPeriodA_interval)
            plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        for tempi=1:length(decisionPeriod_interval)
            plot(temp_ticks_decisionPeriod(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        le = legend(h_line,'memory','offload','Location','northeast','fontsize',8);%10
        le.ItemTokenSize = ones(1,3)*10;
        le.Color = backgounrdColor;
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriod(end)+xlim_padSize2]);
        % xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)+1]);
        xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        ylim([min_baselineDecisionPeriod,max_baselineDecisionPeriod]);
        
        set(gca,'xticklabel',[])
        set(gca,'xtick',[])
        
        xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        xticklabels({'PreFixation','Fixation','','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'','','','Delay1','ChoiceCue','','','',''});
        
        xtickangle(0);
        
        text(temp_ticks_baselinePeriod(end),min_baselineDecisionPeriod,'/','fontsize',12);
        text(temp_ticks_decisionPeriodA(1),min_baselineDecisionPeriod,'/','fontsize',12);
        %
        %     text(temp_ticks_decisionPeriodA(end),min_baselineDecisionPeriod,'/','fontsize',15);
        %     text(temp_ticks_decisionPeriod(1),min_baselineDecisionPeriod,'/','fontsize',15);
        
        set(gca, 'FontSize', 11)
        set(gca,'box','off');% 取消右、上边框
        ylabel(sprintf('dF/F'), 'FontSize', 11, 'FontWeight', 'bold');
        
        % temp_title = title(sprintf('suite2pBID%d,%s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'FontSize',11);%11
        
        a = 1;
        %temp_title.Position(2) = temp_title.Position(2) + 0;
    end
    
    if if_plot_rProb == 1
        %% Plot offloading rate tuning
        fig1 = figure('Name','','NumberTitle','off');
        %set(gcf,'Position',[35+0 40+0 273*0.8 240*0.97]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 40+0 218*0.68 233*0.76]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 40+0 148*1.34 177*1.18]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 40+0 148*1.34*0.85 177*1.18]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 40+0 148*1.34*0.85*0.87 177*1.18]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 40+0 148*1.34*0.85*0.87 146*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 40+0 148*1.34*0.85*0.87*1.5*1.1 146*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 40+0 148*1.34*0.85*0.87 146*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 40+0 148*1.34*0.85*0.87 146*1.2*1.068]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 40+0 148*1.34*0.85*0.87 146*1.2*1.068*0.97]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[35+0 40+0 148*1.34*0.85*0.87 146*1.2*1.068*0.97*0.9*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        x = cellID_F_dff_decisionBin1_seqMerged;
        y = offloadingProb_inOne;
        
        tempBoolIndex = ~isnan(x);
        
        h = [];
        for tempi=1:4
            temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
            tempIndex = find(tempBoolIndex(temp_range)==true);
            temp_range2 = temp_range(tempIndex);
            
            %temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
            %temp_size = 10;
            temp_size = (((tempi).^3)*2 + 3) ./ 5;
            temp_h = scatter(x(temp_range2), y(temp_range2), ...
                temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
            h = [h temp_h]; %#ok<*AGROW>
        end
        
        
        
        [temp_min,temp_max] = bounds(x);
        
        n = 1;
        [p_mapping,S] = polyfit(x,y,n);
        x_fit = temp_min:0.01:temp_max;
        y_fit = polyval(p_mapping,x_fit);
        
        plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);%2
        hold on
        
        set(gca,'linewidth',1.5)
        xlim([temp_min temp_max]);
        ylim([0 1]);
        
        xticks([0 0.4]);
        if cellIndex_suite2p_temptempRaw == 2
            xticks([0 4]);
        end
        if cellIndex_suite2p_temptempRaw == 695
            xticks([0 0.1]);
        end
        if cellIndex_suite2p_temptempRaw == 20
            xticks([0.1 1.5]);
        end           
        if cellIndex_suite2p_temptempRaw == 25
            xticks([0.1 0.8]);
        end        
        if cellIndex_suite2p_temptempRaw == 816
            xticks([0.1 0.3]);
        end  
        if cellIndex_suite2p_temptempRaw == 124 % 0510A
            xticks([0.2 1.5]);
        end
        if cellIndex_suite2p_temptempRaw == 1 % 0510A
            xticks([0 0.6]);
        end        
        if cellIndex_suite2p_temptempRaw == 29 && if_monkey_D0_Z1 == 1  % 0124A_Z
            xticks([0.3 1.2]);
        end 
        
        yticks([0 1]);

        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        
        %text(temp_min+(temp_max-temp_min)*0.55,0.10,sprintf('r2=%.3f',r2_rProb), 'fontsize',9,'FontWeight','normal');
        %text(temp_min+(temp_max-temp_min)*0.45,0.10,sprintf('r2 = %.3f',r2_rProb), 'fontsize',7.5,'FontWeight','normal');
        
        %xlabel(sprintf('dF/F, r2=%.3f',r2_rProb), 'FontSize', 11, 'FontWeight', 'bold');
        xlabel(sprintf('dF/F'), 'FontSize', 9, 'FontWeight', 'normal');        
        ylabel(sprintf('cellID=%d\nOffloading rate',cellIndex_suite2p_temptempRaw), 'FontSize', 9, 'FontWeight', 'normal');
        title(sprintf('r2 = %.3f',r2_rProb),'FontSize',9);%11
        
        
        
        %% Plot behavior precision tuning B
        fig1 = figure('Name','','NumberTitle','off');
        %set(gcf,'Position',[35+0 40+0 273*0.8 240*0.97]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 40+0 183 146*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 40+0 148*1.34*0.85 146*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 40+0 148*1.34*0.85*0.87 146*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 40+300 148*1.34*0.85*0.87 146*1.2*1.068]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[35+0 40+300 148*1.34*0.85*0.87 146*1.2*1.068*0.97]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        
        x = cellID_F_dff_decisionBin1_seqMerged;
        y = seqPrecision_behavior_56;
        
        tempBoolIndex = ~isnan(x);
        
        h = [];
        for tempi=1:4
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
        
        
        
        [temp_min,temp_max] = bounds(x);
        
        n = 1;
        [p_mapping,S] = polyfit(x,y,n);
        x_fit = temp_min:0.01:temp_max;
        y_fit = polyval(p_mapping,x_fit);
        
        plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);
        hold on
        
        set(gca,'linewidth',1.5)
        xlim([temp_min temp_max]);
        ylim([0 1]);
        
        xticks([0 0.4]);
        yticks([0 1]);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        
        %text(temp_min+(temp_max-temp_min)*0.55,0.92,sprintf('r2=%.3f',r2_precision), 'fontsize',9,'FontWeight','normal');
        %text(temp_min+(temp_max-temp_min)*0.48,0.92,sprintf('r2 = %.3f',r2_precision), 'fontsize',7.5,'FontWeight','normal');
        
        
        %xlabel(sprintf('dF/F, r2=%.3f',r2_precision), 'FontSize', 11, 'FontWeight', 'bold');
        xlabel(sprintf('dF/F'), 'FontSize', 9, 'FontWeight', 'normal');        
        %ylabel(sprintf('cellID=%d\nBehavior precision',cellIndex_suite2p_temptempRaw), 'FontSize', 12, 'FontWeight', 'bold');
        %ylabel(sprintf('cellID=%d\nPrecision',cellIndex_suite2p_temptempRaw), 'FontSize', 10, 'FontWeight', 'normal');
        ylabel(sprintf('cellID=%d\nAccuracy',cellIndex_suite2p_temptempRaw), 'FontSize', 9, 'FontWeight', 'normal');
        title(sprintf('r2 = %.3f',r2_precision),'FontSize',9);%11
        
    end
    
    
    
    if if_plot_rProb_B == 1
        %% Plot offloading rate tuning
        fig1 = figure('Name','','NumberTitle','off');
        set(gcf,'Position',[35+0 40+0 148*1.34*0.85*0.87*1.12*1.05 146*1.2*1.068*0.97*0.71]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        
        x = cellID_F_dff_decisionBin1_seqMerged;
        y = offloadingProb_inOne;
        
        tempBoolIndex = ~isnan(x);
        
        h = [];
        for tempi=1:4
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
        
        
        
        [temp_min,temp_max] = bounds(x);
        
        n = 1;
        [p_mapping,S] = polyfit(x,y,n);
        x_fit = temp_min:0.01:temp_max;
        y_fit = polyval(p_mapping,x_fit);
        
        plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);
        hold on
        
        set(gca,'linewidth',1.5)
        xlim([temp_min temp_max]);
        ylim([0 1]);
        
        %xticks([0 0.4]);
        yticks([0 1]);

        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        
        %text(temp_min+(temp_max-temp_min)*0.55,0.10,sprintf('r2=%.3f',r2_rProb), 'fontsize',9,'FontWeight','normal');
        %text(temp_min+(temp_max-temp_min)*0.45,0.10,sprintf('r2 = %.3f',r2_rProb), 'fontsize',7.5,'FontWeight','normal');
        
        %xlabel(sprintf('dF/F, r2=%.3f',r2_rProb), 'FontSize', 11, 'FontWeight', 'bold');
        xlabel(sprintf('dF/F'), 'FontSize', 9, 'FontWeight', 'normal');        
        ylabel(sprintf('cellID=%d\nOffloading rate',cellIndex_suite2p_temptempRaw), 'FontSize', 9, 'FontWeight', 'normal');
        title(sprintf('r2 = %.3f',r2_rProb),'FontSize',9);%11
                        
        
        %% Plot behavior precision tuning B
        fig1 = figure('Name','','NumberTitle','off');
        %set(gcf,'Position',[35+0 40+300 148*1.34*0.85*0.87 146*1.2*1.068*0.97]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 40+300 148*1.34*0.85*0.87*1.12*1.05 146*1.2*1.068*0.97*0.71*1.1]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[35+0 40+300 148*1.34*0.85*0.87*1.12*1.05 146*1.2*1.068*0.97*0.71]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        
        
        x = cellID_F_dff_decisionBin1_seqMerged;
        y = seqPrecision_behavior_56;
        
        tempBoolIndex = ~isnan(x);
        
        h = [];
        for tempi=1:4
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
        
        
        
        [temp_min,temp_max] = bounds(x);
        
        n = 1;
        [p_mapping,S] = polyfit(x,y,n);
        x_fit = temp_min:0.01:temp_max;
        y_fit = polyval(p_mapping,x_fit);
        
        plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);
        hold on
        
        set(gca,'linewidth',1.5)
        xlim([temp_min temp_max]);
        ylim([0 1]);
        
        %xticks([0 0.4]);
        yticks([0 1]);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        
        %text(temp_min+(temp_max-temp_min)*0.55,0.92,sprintf('r2=%.3f',r2_precision), 'fontsize',9,'FontWeight','normal');
        %text(temp_min+(temp_max-temp_min)*0.48,0.92,sprintf('r2 = %.3f',r2_precision), 'fontsize',7.5,'FontWeight','normal');
        
        
        %xlabel(sprintf('dF/F, r2=%.3f',r2_precision), 'FontSize', 11, 'FontWeight', 'bold');
        xlabel(sprintf('dF/F'), 'FontSize', 9, 'FontWeight', 'normal');        
        %ylabel(sprintf('cellID=%d\nBehavior precision',cellIndex_suite2p_temptempRaw), 'FontSize', 12, 'FontWeight', 'bold');
        %ylabel(sprintf('cellID=%d\nPrecision',cellIndex_suite2p_temptempRaw), 'FontSize', 10, 'FontWeight', 'normal');
        ylabel(sprintf('cellID=%d\nAccuracy',cellIndex_suite2p_temptempRaw), 'FontSize', 9, 'FontWeight', 'normal');
        title(sprintf('r2 = %.3f',r2_precision),'FontSize',9);%11
        
    end
    
    
    if if_plot_multiPanelB == 1
        %fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptemp,currentSession(18:22)),'NumberTitle','off');
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 530 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 520*0.97 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 504*0.85 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 504*0.85*0.85*0.95 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 504*0.85 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 504*0.85*0.8 260*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[35+0 42+0 504*0.85*0.8 260*0.9*0.85*1.02*1.01]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','tight'); %#ok<*NASGU>
        
        t.Title.String = sprintf('suite2pBID%d,%s',...
            cellIndex_suite2p_temptempRaw,FOVName_currentFOV2);
        t.Title.FontSize = 11;
        t.Title.Interpreter = 'none';
        
        
        min_baselineDecisionPeriod_raw = min([...
            cellID_F_dff_baselinePeriod_choiceMemory_mean-cellID_F_dff_baselinePeriod_choiceMemory_sem,...
            cellID_F_dff_baselinePeriod_choiceOffload_mean-cellID_F_dff_baselinePeriod_choiceOffload_sem,...
            cellID_F_dff_decisionPeriodA_choiceMemory_mean-cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
            cellID_F_dff_decisionPeriodA_choiceOffload_mean-cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        
        max_baselineDecisionPeriod_raw = max([...
            cellID_F_dff_baselinePeriod_choiceMemory_mean+cellID_F_dff_baselinePeriod_choiceMemory_sem,...
            cellID_F_dff_baselinePeriod_choiceOffload_mean+cellID_F_dff_baselinePeriod_choiceOffload_sem,...
            cellID_F_dff_decisionPeriodA_choiceMemory_mean+cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
            cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        
        min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        %min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.4;
        max_baselineDecisionPeriod = max_baselineDecisionPeriod_raw+(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        
        %% Plot trials raster plot        
        nexttile        
        
        %temptempMax = max_baselineDecisionPeriod*6;
        %temptempMax = max_baselineDecisionPeriod*3;
        %temptempMax = max_baselineDecisionPeriod*3;
        temptempMax = max_baselineDecisionPeriod*2;

        %% plot baselinePeriod
        x = [temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end)];
        y = [1,size(cellID_F_dff_baselinePeriod_choice_collaped,1)];
        C = cellID_F_dff_baselinePeriod_choice_collaped;
        imagesc(x,y,C,[0,temptempMax]);
        hold on
        
        %% plot decisionPeriodA
        x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)];
        y = [1,size(cellID_F_dff_decisionPeriodA_choice_collaped,1)];
        C = cellID_F_dff_decisionPeriodA_choice_collaped;
        imagesc(x,y,C,[0,temptempMax]);
        hold on
        
        %% plot decisionPeriod
        x = [temp_ticks_decisionPeriod(1) temp_ticks_decisionPeriod(end)];
        y = [1,size(cellID_F_dff_decisionPeriod_choice_collaped,1)];
        C = cellID_F_dff_decisionPeriod_choice_collaped;
        imagesc(x,y,C,[0,temptempMax]);
        hold on
        
        %% Others
        for tempi=1:length(baselinePeriod_interval)
            if tempi == 2
                continue
            end
            plot(temp_ticks_baselinePeriod(tempi)*[1 1],[1 size(cellID_F_dff_baselinePeriod_choice_collaped,1)],...
                '-','LineWidth',1,'Color',[1 1 1]);
            hold on
        end
        
        for tempi=1:length(decisionPeriodA_interval)
            plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[1 size(cellID_F_dff_baselinePeriod_choice_collaped,1)],...
                '-','LineWidth',1,'Color',[1 1 1]);
            hold on
        end
        
        for tempi=1:length(decisionPeriod_interval)
            plot(temp_ticks_decisionPeriod(tempi)*[1 1],[1 size(cellID_F_dff_baselinePeriod_choice_collaped,1)],...
                '-','LineWidth',1,'Color',[1 1 1]);
            hold on
        end
        

        tempColor = [...
            0.25,0.75,0.25;
            0.75,0.25,0.25;
            0.25,0.25,0.25];
        %for tempi=1:2
%             plot([temp_ticks_decisionPeriodA(end)+0.5, temp_ticks_decisionPeriodA(end)+5]...
%                 ,(temp_choiceIndex(2)-0.5)*[1 1],...
%                 '-','LineWidth',1,'Color',tempColor(2,:));
%             hold on
            
            plot([0.5 temp_ticks_decisionPeriodA(end)+0.5],(temp_choiceIndex(2)-0.5)*[1 1],...
                ':','LineWidth',2,'Color',[1 1 1]);%0.5
            hold on
        %end
        
        plot([temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA(1)],(temptempIndex_middle_choiceMemory-0.5)*[1 1],...
            ':','LineWidth',2,'Color',[0 0 0]);%0.5
        hold on
        
        plot([temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA(1)],(temp_choiceIndex(2)+temptempIndex_middle_choiceOffload-0.5)*[1 1],...
            ':','LineWidth',2,'Color',[0 0 0]);%0.5
        hold on        
        
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        xlim([1 temp_ticks_decisionPeriodA(end)+1]);        
        ylim([1-1 size(cellID_F_dff_decisionPeriod_choice_collaped,1)+1]);
        
        
        set(gca,'xticklabel',[])
        set(gca,'xtick',[])
                
        yticks([]);
        
        xtickangle(0);
        
        text(temp_ticks_baselinePeriod(end),size(cellID_F_dff_baselinePeriod_choice_collaped,1),'/','fontsize',12);
        text(temp_ticks_decisionPeriodA(1),size(cellID_F_dff_baselinePeriod_choice_collaped,1),'/','fontsize',12);
                
        set(gca, 'FontSize', 10)
        set(gca,'box','off');% 取消右、上边框
        ylabel(sprintf('Trials'), 'FontSize', 10, 'FontWeight', 'normal');
        
        
        
        colorbar('fontsize',9);
        
        if if_colormap_loadEnhanced == 1
            load('parula_enhanced');
            colormap(parula_enhanced);
        elseif if_colormap_loadEnhanced == 0
            colormap parula
        end
        
        % %% nexttile
        %% Plot choice tuning
        
        nexttile
        
%         min_baselineDecisionPeriod_raw = min([...
%             cellID_F_dff_baselinePeriod_choiceMemory_mean-cellID_F_dff_baselinePeriod_choiceMemory_sem,...
%             cellID_F_dff_decisionPeriodA_choiceOffload_mean-cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
%         
%         max_baselineDecisionPeriod_raw = max([...
%             cellID_F_dff_baselinePeriod_choiceMemory_mean+cellID_F_dff_baselinePeriod_choiceMemory_sem,...
%             cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
%         
%         min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
%         %min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.4;
%         max_baselineDecisionPeriod = max_baselineDecisionPeriod_raw+(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        
        h_line = [];
        
        %% plot baselinePeriod
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceMemory_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceOffload_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        %% plot decisionPeriodA
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceMemory_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceOffload_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        %% plot decisionPeriod
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceMemory_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceOffload_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        %% Others
        
        for tempi=1:length(baselinePeriod_interval)
            if tempi == 2
                continue
            end
            plot(temp_ticks_baselinePeriod(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        for tempi=1:length(decisionPeriodA_interval)
            plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        for tempi=1:length(decisionPeriod_interval)
            plot(temp_ticks_decisionPeriod(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        %le = legend(h_line,'ChoiceMemory','ChoiceOffload','Location','northeast','fontsize',8);%10
        le = legend(h_line,'Choice-memory','Choice-offload','Location','northeast','fontsize',8);%10                        
        le.ItemTokenSize = ones(1,3)*10;
        le.Color = backgounrdColor;
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriod(end)+xlim_padSize2]);
        % xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)+1]);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        xlim([1 temp_ticks_decisionPeriodA(end)+1]);        
        ylim([min_baselineDecisionPeriod,max_baselineDecisionPeriod]);
        
        set(gca,'xticklabel',[])
        set(gca,'xtick',[])
        
        %xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        xticks([temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        
        %xticklabels({'PreFixation','Fixation','','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'Fixation','','Delay1','ChoiceCue','','','Decision',''});        
        xticklabels({'Fixation','T1','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'','','','Delay1','ChoiceCue','','','',''});
        
        xtickangle(0);
        
        yticks([0 0.8]);
        
        text(temp_ticks_baselinePeriod(end),min_baselineDecisionPeriod,'/','fontsize',12);
        text(temp_ticks_decisionPeriodA(1),min_baselineDecisionPeriod,'/','fontsize',12);
        %
        %     text(temp_ticks_decisionPeriodA(end),min_baselineDecisionPeriod,'/','fontsize',15);
        %     text(temp_ticks_decisionPeriod(1),min_baselineDecisionPeriod,'/','fontsize',15);
        
        set(gca, 'FontSize', 10)
        set(gca,'box','off');% 取消右、上边框
        ylabel(sprintf('dF/F'), 'FontSize', 10, 'FontWeight', 'normal');
        
        
    end
    
    %%
    if if_plot_multiPanelC == 1
        tempPlotIndicator = 2;%0 for small, 1 for big, 2 for smaller
        
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 504*0.85*0.8*0.7 260*0.9*0.85*1.02*1.01]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 504*0.85*0.8*0.7 260*0.9*0.85*1.02*1.01*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        if tempPlotIndicator == 0
            set(gcf,'Position',[35+0 42+0 240*0.8 174.2*0.95*0.95*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        elseif tempPlotIndicator == 1
            set(gcf,'Position',[35+0 42+0 240*0.8 174.2*0.95*0.95*0.95*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        elseif tempPlotIndicator == 2
            set(gcf,'Position',[35+0 42+0 240*0.8 174.2*0.95*0.95*0.95*0.84]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        end
        
        
        %t = tiledlayout(2,1,'TileSpacing','compact','Padding','tight'); %#ok<*NASGU>
        %t = tiledlayout(2,5,'TileSpacing','compact','Padding','tight'); %#ok<*NASGU>
        if tempPlotIndicator == 0 || tempPlotIndicator == 1
            t = tiledlayout(2,5,'TileSpacing','none','Padding','tight'); %#ok<*NASGU>
        elseif tempPlotIndicator == 2
            t = tiledlayout(7,5,'TileSpacing','none','Padding','tight'); %#ok<*NASGU>
        end
        
        
        %t.Title.String = sprintf('suite2pBID%d,%s',...
        %    cellIndex_suite2p_temptempRaw,FOVName_currentFOV2);
        %t.Title.FontSize = 11;
        %t.Title.Interpreter = 'none';
        
        
        min_baselineDecisionPeriod_raw = min([...
            cellID_F_dff_baselinePeriod_choiceMemory_mean-cellID_F_dff_baselinePeriod_choiceMemory_sem,...
            cellID_F_dff_baselinePeriod_choiceOffload_mean-cellID_F_dff_baselinePeriod_choiceOffload_sem,...
            cellID_F_dff_decisionPeriodA_choiceMemory_mean-cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
            cellID_F_dff_decisionPeriodA_choiceOffload_mean-cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        
        max_baselineDecisionPeriod_raw = max([...
            cellID_F_dff_baselinePeriod_choiceMemory_mean+cellID_F_dff_baselinePeriod_choiceMemory_sem,...
            cellID_F_dff_baselinePeriod_choiceOffload_mean+cellID_F_dff_baselinePeriod_choiceOffload_sem,...
            cellID_F_dff_decisionPeriodA_choiceMemory_mean+cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
            cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        
        min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        %min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.4;
        max_baselineDecisionPeriod = max_baselineDecisionPeriod_raw+(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        
        %% Plot trials raster plot        
        if tempPlotIndicator == 0 || tempPlotIndicator == 1        
        nexttile([1 4])  
        elseif tempPlotIndicator == 2
        nexttile([4 4])  
        end
        
        %temptempMax = max_baselineDecisionPeriod*6;
        %temptempMax = max_baselineDecisionPeriod*3;
        %temptempMax = max_baselineDecisionPeriod*3;
        temptempMax = max_baselineDecisionPeriod*2;

        
        if_sortedChoiceTuning = 0;%1
        
        cellID_F_dff_baselinePeriod_choice_collaped;
        temp_choiceIndex;
        
        cellID_F_dff_baselinePeriod_choice_collaped_mean = mean(cellID_F_dff_baselinePeriod_choice_collaped,2);
        cellID_F_dff_baselinePeriod_choice_collaped_mean_choiceMemory = cellID_F_dff_baselinePeriod_choice_collaped_mean(1:temp_choiceIndex(2));
        cellID_F_dff_baselinePeriod_choice_collaped_mean_choiceOffload = cellID_F_dff_baselinePeriod_choice_collaped_mean((temp_choiceIndex(2)+1):end);
        
        [M_choiceMemory,I_choiceMemory] = sort(cellID_F_dff_baselinePeriod_choice_collaped_mean_choiceMemory,'descend');
        [M_choiceOffload,I_choiceOffload] = sort(cellID_F_dff_baselinePeriod_choice_collaped_mean_choiceOffload,'descend');
        I_choiceOffload = I_choiceOffload + temp_choiceIndex(2);
        
        I_sorted = [I_choiceMemory;I_choiceOffload];
        
        if if_sortedChoiceTuning == 1
            temp_I = I_sorted;
        elseif if_sortedChoiceTuning == 0
            temp_I = 1:size(cellID_F_dff_baselinePeriod_choice_collaped,1);
        end

        
        %% plot baselinePeriod
        x = [temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end)];
        y = [1,size(cellID_F_dff_baselinePeriod_choice_collaped,1)];
        C = cellID_F_dff_baselinePeriod_choice_collaped(temp_I,:);
        imagesc(x,y,C,[0,temptempMax]);
        hold on
        
        %% plot decisionPeriodA
        x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)];
        y = [1,size(cellID_F_dff_decisionPeriodA_choice_collaped,1)];
        C = cellID_F_dff_decisionPeriodA_choice_collaped(temp_I,:);
        imagesc(x,y,C,[0,temptempMax]);
        hold on
        
        %% plot decisionPeriod
        x = [temp_ticks_decisionPeriod(1) temp_ticks_decisionPeriod(end)];
        y = [1,size(cellID_F_dff_decisionPeriod_choice_collaped,1)];
        C = cellID_F_dff_decisionPeriod_choice_collaped(temp_I,:);
        imagesc(x,y,C,[0,temptempMax]);
        hold on
        
        %% Others
        for tempi=1:length(baselinePeriod_interval)
            if tempi == 2
                continue
            end
            plot(temp_ticks_baselinePeriod(tempi)*[1 1],[1 size(cellID_F_dff_baselinePeriod_choice_collaped,1)],...
                '-','LineWidth',1,'Color',[1 1 1]);
            hold on
        end
        
        for tempi=1:length(decisionPeriodA_interval)
            plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[1 size(cellID_F_dff_baselinePeriod_choice_collaped,1)],...
                '-','LineWidth',1,'Color',[1 1 1]);
            hold on
        end
        
        for tempi=1:length(decisionPeriod_interval)
            plot(temp_ticks_decisionPeriod(tempi)*[1 1],[1 size(cellID_F_dff_baselinePeriod_choice_collaped,1)],...
                '-','LineWidth',1,'Color',[1 1 1]);
            hold on
        end
        

        tempColor = [...
            0.25,0.75,0.25;
            0.75,0.25,0.25;
            0.25,0.25,0.25];
        %for tempi=1:2
%             plot([temp_ticks_decisionPeriodA(end)+0.5, temp_ticks_decisionPeriodA(end)+5]...
%                 ,(temp_choiceIndex(2)-0.5)*[1 1],...
%                 '-','LineWidth',1,'Color',tempColor(2,:));
%             hold on
            
            plot([0.5 temp_ticks_decisionPeriodA(end)+0.5],(temp_choiceIndex(2)-0.5)*[1 1],...
                ':','LineWidth',2,'Color',[1 1 1]);%0.5
            hold on
        %end
        
        %plot([temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA(1)],(temptempIndex_middle_choiceMemory-0.5)*[1 1],...
        %    ':','LineWidth',2,'Color',[0 0 0]);%0.5
        %hold on
        
        %plot([temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA(1)],(temp_choiceIndex(2)+temptempIndex_middle_choiceOffload-0.5)*[1 1],...
        %    ':','LineWidth',2,'Color',[0 0 0]);%0.5
        %hold on        
        
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        %xlim([1 temp_ticks_decisionPeriodA(end)+1]);  
        xlim([1 temp_ticks_baselinePeriod(end)+4]);  
        ylim([1-1 size(cellID_F_dff_decisionPeriod_choice_collaped,1)+1]);                                
        
        set(gca,'xticklabel',[])
        set(gca,'xtick',[])
                
        yticks([]);
        
            tempLabelStr = string(1:2);
            tempLabelStr(1) = 'Memory';
            tempLabelStr(2) = 'Offload';
            ytl=string(tempLabelStr);
            % 设置ttext的x坐标位置
            %ytext_xp=(xt(1))*ones(1,length(tempLabelStr))-2.4;%-2.4
            ytext_xp=(xt(1))*ones(1,length(tempLabelStr))-6.5+4;%-2.4
            % 设置ttext的y坐标位置
            % ytext_yp=yt;
            ytext_yp = nan(1,2);
            ytext_yp(1) = temp_choiceIndex(2)/2 - 20;
            ytext_yp(2) = temp_choiceIndex(2) + temptempIndex_middle_choiceOffload + 20;
            ytext_yp = ytext_yp + 0;%0.5
            
            if tempPlotIndicator == 0
            text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',90,'fontsize',6.5);%8
            elseif tempPlotIndicator == 1
            text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',90,'fontsize',8);%8                
            elseif tempPlotIndicator == 2
            text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',90,'fontsize',6.5);%8                            
            end
            set(gca,'yticklabel','');
        
        
        xtickangle(0);
        
        %text(temp_ticks_baselinePeriod(end),size(cellID_F_dff_baselinePeriod_choice_collaped,1),'/','fontsize',12);
        %text(temp_ticks_decisionPeriodA(1),size(cellID_F_dff_baselinePeriod_choice_collaped,1),'/','fontsize',12);
                
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        %ylabel(sprintf('Trials'), 'FontSize', 10, 'FontWeight', 'normal');                       
        
        c = colorbar('fontsize',9);
        %c.Position = c.Position+[0.06 0.04 -0.025 -0.05];
        c.Position = c.Position+[0.01 0.055 -0.035 -0.07];
        
        c.Ticks = [0 1];
        
        if if_monkey_D0_Z1 == 1 && cellIndex_suite2p_temptempRaw == 364
            c.Ticks = [0 0.6];
        end
        
            if if_colormap_loadEnhanced == 1
                %load('parula_enhanced');
                %colormap(parula_enhanced);
                %load('gray_enhanced');
                %colormap(gray_enhanced);
                
                %load('gray_reversed_enhanced');
                %colormap(gray_reversed_enhanced);
                
                load('coolwarm_enhanced');
                colormap(coolwarm_enhanced);
                
                %load('coolwarm_n11n_enhanced');
                %colormap(coolwarm_n11n_enhanced);

            elseif if_colormap_loadEnhanced == 0
                %colormap parula
                %colormap gray
                
                %temp1 = gray;
                %temp1 = temp1(end:-1:1,:);
                %colormap(temp1);
                
                colormap(coolwarm());
                
                %temp1 = coolwarm(300);
                %temp2 = ((300-256)/2)+1;
                %temp3 = temp1(temp2:temp2+255,:);
                %colormap(temp3);
            end
            
        
        % %% nexttile
        %% Plot choice tuning
        
        %nexttile([1 4])  
        if tempPlotIndicator == 0 || tempPlotIndicator == 1        
        nexttile([1 4])  
        elseif tempPlotIndicator == 2
        nexttile([3 4])  
        end        
        
%         min_baselineDecisionPeriod_raw = min([...
%             cellID_F_dff_baselinePeriod_choiceMemory_mean-cellID_F_dff_baselinePeriod_choiceMemory_sem,...
%             cellID_F_dff_decisionPeriodA_choiceOffload_mean-cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
%         
%         max_baselineDecisionPeriod_raw = max([...
%             cellID_F_dff_baselinePeriod_choiceMemory_mean+cellID_F_dff_baselinePeriod_choiceMemory_sem,...
%             cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
%         
%         min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
%         %min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.4;
%         max_baselineDecisionPeriod = max_baselineDecisionPeriod_raw+(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        
        h_line = [];
        
        %% plot baselinePeriod
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceMemory_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceOffload_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        %% plot decisionPeriodA
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceMemory_mean;
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceOffload_mean;
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        %% plot decisionPeriod
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceMemory_mean;
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceOffload_mean;
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        %% Others
        
        for tempi=1:length(baselinePeriod_interval)
            if tempi == 2
                continue
            end
            plot(temp_ticks_baselinePeriod(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod+(max_baselineDecisionPeriod-min_baselineDecisionPeriod)*0.35],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        for tempi=1:length(decisionPeriodA_interval)
            plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod+(max_baselineDecisionPeriod-min_baselineDecisionPeriod)*0.35],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        for tempi=1:length(decisionPeriod_interval)
            plot(temp_ticks_decisionPeriod(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod+(max_baselineDecisionPeriod-min_baselineDecisionPeriod)*0.35],...
                '-','LineWidth',1,'Color',[0 0 0]);
            hold on
        end
        
        %le = legend(h_line,'ChoiceMemory','ChoiceOffload','Location','northeast','fontsize',8);%10
        %le = legend(h_line,'Choice-memory','Choice-offload','Location','northeast','fontsize',6.5,'NumColumns',1);%10                        
        if tempPlotIndicator == 0
        le = legend(h_line,'Memory','Offload','Location','northeast','fontsize',6.5,'NumColumns',1);%10
        elseif tempPlotIndicator == 1
        le = legend(h_line,'Memory','Offload','Location','northeast','fontsize',8,'NumColumns',1);%10
        elseif tempPlotIndicator == 2
        le = legend(h_line,'Memory','Offload','Location','northeast','fontsize',6.5,'NumColumns',2);%10        
        end
        le.ItemTokenSize = ones(1,3)*10;
        le.Color = backgounrdColor;
        legend('boxoff');
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriod(end)+xlim_padSize2]);
        % xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)+1]);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        %xlim([1 temp_ticks_decisionPeriodA(end)+1]);        
        xlim([1 temp_ticks_baselinePeriod(end)+4]);
        ylim([min_baselineDecisionPeriod,...
            max_baselineDecisionPeriod+(max_baselineDecisionPeriod-min_baselineDecisionPeriod)*0.35]);
        
        set(gca,'xticklabel',[])
        set(gca,'xtick',[])
        
        %xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        xticks([temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        
        %xticklabels({'PreFixation','Fixation','','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'Fixation','','Delay1','ChoiceCue','','','Decision',''});        
        xticklabels({'Fixation','T1','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'','','','Delay1','ChoiceCue','','','',''});
        
        xtickangle(0);
        
        yticks([0 0.8]);
       
        if if_monkey_D0_Z1 == 1 && cellIndex_suite2p_temptempRaw == 364
            ylim([0,max_baselineDecisionPeriod+(max_baselineDecisionPeriod-min_baselineDecisionPeriod)*0.35]);
            yticks([0 0.3]);
        end
        
        %text(temp_ticks_baselinePeriod(end),min_baselineDecisionPeriod,'/','fontsize',12);
        %text(temp_ticks_decisionPeriodA(1),min_baselineDecisionPeriod,'/','fontsize',12);
        %
        %     text(temp_ticks_decisionPeriodA(end),min_baselineDecisionPeriod,'/','fontsize',15);
        %     text(temp_ticks_decisionPeriod(1),min_baselineDecisionPeriod,'/','fontsize',15);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        ylabel(sprintf('dF/F'), 'FontSize', 9, 'FontWeight', 'normal');
        
        
    end
    
    
    
    
    if if_plot_multiPanelD == 1
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 504*0.85*0.8*0.7 260*0.9*0.85*1.02*1.01]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 504*0.85*0.8*0.7 260*0.9*0.85*1.02*1.01*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        %set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67 127*1.05*0.95*0.93]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67*1.04 127*1.05*0.95*0.93*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
                
        
        %t = tiledlayout(2,1,'TileSpacing','compact','Padding','tight'); %#ok<*NASGU>
        %t = tiledlayout(2,5,'TileSpacing','compact','Padding','tight'); %#ok<*NASGU>
        t = tiledlayout(2,5,'TileSpacing','none','Padding','tight'); %#ok<*NASGU>        
        
        %t.Title.String = sprintf('suite2pBID%d,%s',...
        %    cellIndex_suite2p_temptempRaw,FOVName_currentFOV2);
        %t.Title.FontSize = 11;
        %t.Title.Interpreter = 'none';
        
        
        min_baselineDecisionPeriod_raw = min([...
            cellID_F_dff_baselinePeriod_choiceMemory_mean-cellID_F_dff_baselinePeriod_choiceMemory_sem,...
            cellID_F_dff_baselinePeriod_choiceOffload_mean-cellID_F_dff_baselinePeriod_choiceOffload_sem,...
            cellID_F_dff_decisionPeriodA_choiceMemory_mean-cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
            cellID_F_dff_decisionPeriodA_choiceOffload_mean-cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        
        max_baselineDecisionPeriod_raw = max([...
            cellID_F_dff_baselinePeriod_choiceMemory_mean+cellID_F_dff_baselinePeriod_choiceMemory_sem,...
            cellID_F_dff_baselinePeriod_choiceOffload_mean+cellID_F_dff_baselinePeriod_choiceOffload_sem,...
            cellID_F_dff_decisionPeriodA_choiceMemory_mean+cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
            cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        
        min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        %min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.4;
        max_baselineDecisionPeriod = max_baselineDecisionPeriod_raw+(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        
        %% Plot trials raster plot        
        nexttile([1 4])  
                
        
        %temptempMax = max_baselineDecisionPeriod*6;
        %temptempMax = max_baselineDecisionPeriod*3;
        %temptempMax = max_baselineDecisionPeriod*3;
        %temptempMax = max_baselineDecisionPeriod*2;
        temptempMax = max_baselineDecisionPeriod*1.5;

        if_sortedChoiceTuning = 1;
        
        cellID_F_dff_decisionPeriodA_choice_collaped;
        temp_choiceIndex;
        
        cellID_F_dff_decisionPeriodA_choice_collaped_mean = mean(cellID_F_dff_decisionPeriodA_choice_collaped,2);
        cellID_F_dff_decisionPeriodA_choice_collaped_mean_choiceMemory = cellID_F_dff_decisionPeriodA_choice_collaped_mean(1:temp_choiceIndex(2));
        cellID_F_dff_decisionPeriodA_choice_collaped_mean_choiceOffload = cellID_F_dff_decisionPeriodA_choice_collaped_mean((temp_choiceIndex(2)+1):end);
        
        [M_choiceMemory,I_choiceMemory] = sort(cellID_F_dff_decisionPeriodA_choice_collaped_mean_choiceMemory,'descend');
        [M_choiceOffload,I_choiceOffload] = sort(cellID_F_dff_decisionPeriodA_choice_collaped_mean_choiceOffload,'descend');
        I_choiceOffload = I_choiceOffload + temp_choiceIndex(2);
        
        I_sorted = [I_choiceMemory;I_choiceOffload];
        
        if if_sortedChoiceTuning == 1
            temp_I = I_sorted;
        elseif if_sortedChoiceTuning == 0
            temp_I = 1:size(cellID_F_dff_decisionPeriodA_choice_collaped,1);
        end
        
        
        %% plot baselinePeriod
        x = [temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end)];
        y = [1,size(cellID_F_dff_baselinePeriod_choice_collaped,1)];
        C = cellID_F_dff_baselinePeriod_choice_collaped(temp_I,:);
        imagesc(x,y,C,[0,temptempMax]);
        hold on
        
        %% plot decisionPeriodA
        x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)];
        y = [1,size(cellID_F_dff_decisionPeriodA_choice_collaped,1)];
        C = cellID_F_dff_decisionPeriodA_choice_collaped(temp_I,:);
        imagesc(x,y,C,[0,temptempMax]);
        hold on
        
        %% plot decisionPeriod
        x = [temp_ticks_decisionPeriod(1) temp_ticks_decisionPeriod(end)];
        y = [1,size(cellID_F_dff_decisionPeriod_choice_collaped,1)];
        C = cellID_F_dff_decisionPeriod_choice_collaped(temp_I,:);
        imagesc(x,y,C,[0,temptempMax]);
        hold on
        
        %% Others
        %         for tempi=1:length(baselinePeriod_interval)
        %             if tempi == 2
        %                 continue
        %             end
        %             plot(temp_ticks_baselinePeriod(tempi)*[1 1],[1 size(cellID_F_dff_baselinePeriod_choice_collaped,1)],...
        %                 '-','LineWidth',1,'Color',[1 1 1]);
        %             hold on
        %         end
        %
        %         for tempi=1:length(decisionPeriodA_interval)
        %             plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[1 size(cellID_F_dff_baselinePeriod_choice_collaped,1)],...
        %                 '-','LineWidth',1,'Color',[1 1 1]);
        %             hold on
        %         end
        %
        %         for tempi=1:length(decisionPeriod_interval)
        %             plot(temp_ticks_decisionPeriod(tempi)*[1 1],[1 size(cellID_F_dff_baselinePeriod_choice_collaped,1)],...
        %                 '-','LineWidth',1,'Color',[1 1 1]);
        %             hold on
        %         end
        

        tempColor = [...
            0.25,0.75,0.25;
            0.75,0.25,0.25;
            0.25,0.25,0.25];
        %for tempi=1:2
%             plot([temp_ticks_decisionPeriodA(end)+0.5, temp_ticks_decisionPeriodA(end)+5]...
%                 ,(temp_choiceIndex(2)-0.5)*[1 1],...
%                 '-','LineWidth',1,'Color',tempColor(2,:));
%             hold on
            
            %plot([0.5 temp_ticks_decisionPeriodA(end)+0.5],(temp_choiceIndex(2)-0.5)*[1 1],...
            %   ':','LineWidth',2,'Color',[1 1 1]);%0.5
            %hold on
            
        %end
        
        plot([temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA(1)],(temptempIndex_middle_choiceMemory-0.5)*[1 1],...
            ':','LineWidth',2,'Color',[0 0 0]);%0.5
        hold on
        
        plot([temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA(1)],(temp_choiceIndex(2)+temptempIndex_middle_choiceOffload-0.5)*[1 1],...
            ':','LineWidth',2,'Color',[0 0 0]);%0.5
        hold on        
        
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        %xlim([1 temp_ticks_decisionPeriodA(end)+1]);  
        %xlim([1 temp_ticks_baselinePeriod(end)+4]);  
        xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]);        
        ylim([1-1 size(cellID_F_dff_decisionPeriod_choice_collaped,1)+1]);                                
        
        set(gca,'xticklabel',[])
        %set(gca,'xtick',[])
        xticks([temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        
                
        yticks([]);
        
            tempLabelStr = string(1:2);
            tempLabelStr(1) = 'Memory';
            tempLabelStr(2) = 'Offload';
            ytl=string(tempLabelStr);
            % 设置ttext的x坐标位置
            %ytext_xp=(xt(1))*ones(1,length(tempLabelStr))-2.4;%-2.4
            %ytext_xp=(xt(1))*ones(1,length(tempLabelStr))-6.5+4;%-2.4
            ytext_xp=(xt(1)+temp_ticks_decisionPeriodA(1))*ones(1,length(tempLabelStr))-6.5+0.4;%-2.4            
            % 设置ttext的y坐标位置
            % ytext_yp=yt;
            ytext_yp = nan(1,2);
            ytext_yp(1) = temp_choiceIndex(2)/2 - 20;
            ytext_yp(2) = temp_choiceIndex(2) + temptempIndex_middle_choiceOffload + 20;
            ytext_yp = ytext_yp + 0;%0.5
            
            %text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',90,'fontsize',6.5);%8
            text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',0,'fontsize',6.5);%8            
            set(gca,'yticklabel','');
            
            text(ytext_xp(1),temp_choiceIndex(2),'---------','HorizontalAlignment','center','rotation',0,'fontsize',8);%8            
            
        
        
        xtickangle(0);
        
        %text(temp_ticks_baselinePeriod(end),size(cellID_F_dff_baselinePeriod_choice_collaped,1),'/','fontsize',12);
        %text(temp_ticks_decisionPeriodA(1),size(cellID_F_dff_baselinePeriod_choice_collaped,1),'/','fontsize',12);
                
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        %ylabel(sprintf('Trials'), 'FontSize', 10, 'FontWeight', 'normal');                       
        
        c = colorbar('fontsize',6.5);%9
        %c.Position = c.Position+[0.06 0.04 -0.025 -0.05];
        c.Position = c.Position+[0.01 0.06 -0.038 -0.075];
        
        
        %if if_colormap_loadEnhanced == 1
        %    load('parula_enhanced');
        %    colormap(parula_enhanced);
        %elseif if_colormap_loadEnhanced == 0
        %    colormap parula
        %end
        colormap parula
        
        % %% nexttile
        %% Plot choice tuning
        
        nexttile([1 4])  
        
%         min_baselineDecisionPeriod_raw = min([...
%             cellID_F_dff_baselinePeriod_choiceMemory_mean-cellID_F_dff_baselinePeriod_choiceMemory_sem,...
%             cellID_F_dff_decisionPeriodA_choiceOffload_mean-cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
%         
%         max_baselineDecisionPeriod_raw = max([...
%             cellID_F_dff_baselinePeriod_choiceMemory_mean+cellID_F_dff_baselinePeriod_choiceMemory_sem,...
%             cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
%         
%         min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
%         %min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.4;
%         max_baselineDecisionPeriod = max_baselineDecisionPeriod_raw+(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        
        h_line = [];
        
        %% plot baselinePeriod
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceMemory_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceOffload_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        %% plot decisionPeriodA
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceMemory_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceOffload_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        %% plot decisionPeriod
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceMemory_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceOffload_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        %% Others
        
        %         for tempi=1:length(baselinePeriod_interval)
        %             %if tempi == 2
        %             %    continue
        %             %end
        %             plot(temp_ticks_baselinePeriod(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod+(max_baselineDecisionPeriod-min_baselineDecisionPeriod)*0.35],...
        %                 '-','LineWidth',1,'Color',[0 0 0]);
        %             hold on
        %         end
        %
        %         for tempi=1:length(decisionPeriodA_interval)
        %             plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod+(max_baselineDecisionPeriod-min_baselineDecisionPeriod)*0.35],...
        %                 '-','LineWidth',1,'Color',[0 0 0]);
        %             hold on
        %         end
        %
        %         for tempi=1:length(decisionPeriod_interval)
        %             plot(temp_ticks_decisionPeriod(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod+(max_baselineDecisionPeriod-min_baselineDecisionPeriod)*0.35],...
        %                 '-','LineWidth',1,'Color',[0 0 0]);
        %             hold on
        %         end
        
        %le = legend(h_line,'ChoiceMemory','ChoiceOffload','Location','northeast','fontsize',8);%10
        le = legend(h_line,'Choice-memory','Choice-offload','Location','northwest','fontsize',6.5,'NumColumns',2);%10                        
        le.ItemTokenSize = ones(1,3)*10;
        le.Color = backgounrdColor;
        legend('boxoff');
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriod(end)+xlim_padSize2]);
        xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        %xlim([1 temp_ticks_decisionPeriodA(end)+1]);        
        %xlim([1 temp_ticks_baselinePeriod(end)+4]);
        ylim([min_baselineDecisionPeriod,...
            max_baselineDecisionPeriod+(max_baselineDecisionPeriod-min_baselineDecisionPeriod)*0.35]);
        
        set(gca,'xticklabel',[])
        set(gca,'xtick',[])
        
        %xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        xticks([temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        
        %xticklabels({'PreFixation','Fixation','','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'Fixation','','Delay1','ChoiceCue','','','Decision',''});        
        xticklabels({'Fixation','T1','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'','','','Delay1','ChoiceCue','','','',''});
        
        xtickangle(0);
        
        %yticks([0 0.8]);
        
        %text(temp_ticks_baselinePeriod(end),min_baselineDecisionPeriod,'/','fontsize',12);
        %text(temp_ticks_decisionPeriodA(1),min_baselineDecisionPeriod,'/','fontsize',12);
        %
        %     text(temp_ticks_decisionPeriodA(end),min_baselineDecisionPeriod,'/','fontsize',15);
        %     text(temp_ticks_decisionPeriod(1),min_baselineDecisionPeriod,'/','fontsize',15);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        ylabel(sprintf('dF/F'), 'FontSize', 8, 'FontWeight', 'normal');
        
        
    end
    
    
    if if_plot_multiPanelD2 == 1
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        
        %set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67 127*1.05*0.95*0.93]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        %set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67*1.04 127*1.05*0.95*0.93*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        set(gcf,'Position',[35+0 42+0 233.4*0.93 176.7*0.86]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
                
        
        t = tiledlayout(2,5,'TileSpacing','none','Padding','tight'); %#ok<*NASGU>        
                
        
        min_baselineDecisionPeriod_raw = min([...
            cellID_F_dff_baselinePeriod_choiceMemory_mean-cellID_F_dff_baselinePeriod_choiceMemory_sem,...
            cellID_F_dff_baselinePeriod_choiceOffload_mean-cellID_F_dff_baselinePeriod_choiceOffload_sem,...
            cellID_F_dff_decisionPeriodA_choiceMemory_mean-cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
            cellID_F_dff_decisionPeriodA_choiceOffload_mean-cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        
        max_baselineDecisionPeriod_raw = max([...
            cellID_F_dff_baselinePeriod_choiceMemory_mean+cellID_F_dff_baselinePeriod_choiceMemory_sem,...
            cellID_F_dff_baselinePeriod_choiceOffload_mean+cellID_F_dff_baselinePeriod_choiceOffload_sem,...
            cellID_F_dff_decisionPeriodA_choiceMemory_mean+cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
            cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        
        min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        %min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.4;
        max_baselineDecisionPeriod = max_baselineDecisionPeriod_raw+(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        
        %% Plot trials raster plot        
        nexttile([1 4])  
                
        
        %temptempMax = max_baselineDecisionPeriod*6;
        %temptempMax = max_baselineDecisionPeriod*3;
        %temptempMax = max_baselineDecisionPeriod*3;
        %temptempMax = max_baselineDecisionPeriod*2;
        temptempMax = max_baselineDecisionPeriod*1.5;

        temptempMin = 0;

        if cellIndex_suite2p_temptempRaw == 126
            %temptempMax = 0.4;
            
            %temptempMin = 0.10;
        end
        
        
        if_sortedChoiceTuning = 1;%1
        
        cellID_F_dff_decisionPeriodA_choice_collaped;
        temp_choiceIndex;
        
        cellID_F_dff_decisionPeriodA_choice_collaped_mean = mean(cellID_F_dff_decisionPeriodA_choice_collaped,2);
        cellID_F_dff_decisionPeriodA_choice_collaped_mean_choiceMemory = cellID_F_dff_decisionPeriodA_choice_collaped_mean(1:temp_choiceIndex(2));
        cellID_F_dff_decisionPeriodA_choice_collaped_mean_choiceOffload = cellID_F_dff_decisionPeriodA_choice_collaped_mean((temp_choiceIndex(2)+1):end);
        
        [M_choiceMemory,I_choiceMemory] = sort(cellID_F_dff_decisionPeriodA_choice_collaped_mean_choiceMemory,'descend');
        [M_choiceOffload,I_choiceOffload] = sort(cellID_F_dff_decisionPeriodA_choice_collaped_mean_choiceOffload,'descend');
        I_choiceOffload = I_choiceOffload + temp_choiceIndex(2);
        
        I_sorted = [I_choiceMemory;I_choiceOffload];
        
        if if_sortedChoiceTuning == 1
            temp_I = I_sorted;
        elseif if_sortedChoiceTuning == 0
            temp_I = 1:size(cellID_F_dff_decisionPeriodA_choice_collaped,1);
        end
        
        
        %% plot baselinePeriod
        x = [temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end)];
        y = [1,size(cellID_F_dff_baselinePeriod_choice_collaped,1)];
        C = cellID_F_dff_baselinePeriod_choice_collaped(temp_I,:);
        imagesc(x,y,C,[temptempMin,temptempMax]);
        hold on
        
        %% plot sampleLastT
        x = [temp_ticks_sampleLastT(1) temp_ticks_sampleLastT(end)];
        y = [1,size(cellID_F_dff_sampleLastT_choice_collaped,1)];
        C = cellID_F_dff_sampleLastT_choice_collaped(temp_I,:);
        imagesc(x,y,C,[temptempMin,temptempMax]);
        hold on
        
        %% plot decisionPeriodA
        x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)];
        y = [1,size(cellID_F_dff_decisionPeriodA_choice_collaped,1)];
        C = cellID_F_dff_decisionPeriodA_choice_collaped(temp_I,:);
        imagesc(x,y,C,[temptempMin,temptempMax]);
        hold on
        
        %% plot decisionPeriod
        x = [temp_ticks_decisionPeriod(1) temp_ticks_decisionPeriod(end)];
        y = [1,size(cellID_F_dff_decisionPeriod_choice_collaped,1)];
        C = cellID_F_dff_decisionPeriod_choice_collaped(temp_I,:);
        imagesc(x,y,C,[temptempMin,temptempMax]);
        hold on

        
        %%
        
        plot((temp_ticks_decisionPeriodA(1)-0.25)*[1 1],[1,size(cellID_F_dff_decisionPeriodA_choice_collaped,1)],...
            '-','LineWidth',1,'Color',[1 1 1]);%3,2
        hold on
        
        plot((temp_ticks_decisionPeriodA(2)+0.25)*[1 1],[1,size(cellID_F_dff_decisionPeriodA_choice_collaped,1)],...
            '-','LineWidth',1,'Color',[1 1 1]);%3,2
        hold on
        
        
        
        
        tempColor = [...
            0.25,0.75,0.25;
            0.75,0.25,0.25;
            0.25,0.25,0.25];                

        
        %plot([temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA(1)],(temptempIndex_middle_choiceMemory-0.5)*[1 1],...
        %    ':','LineWidth',2,'Color',[0 0 0]);%0.5
        %hold on
        %
        %plot([temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA(1)],(temp_choiceIndex(2)+temptempIndex_middle_choiceOffload-0.5)*[1 1],...
        %    ':','LineWidth',2,'Color',[0 0 0]);%0.5
        %hold on        
        
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        %xlim([1 temp_ticks_decisionPeriodA(end)+1]);  
        %xlim([1 temp_ticks_baselinePeriod(end)+4]);  
        %xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]);        
        xlim([temp_ticks_sampleLastT(1) temp_ticks_decisionPeriodA(end)]);                
        ylim([1-1 size(cellID_F_dff_decisionPeriod_choice_collaped,1)+1]);                                
        
        set(gca,'xticklabel',[])
        %set(gca,'xtick',[])
        xticks([temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        
                
        yticks([]);
        
            tempLabelStr = string(1:2);
            tempLabelStr(1) = 'Memory';
            tempLabelStr(2) = 'Offload';
            ytl=string(tempLabelStr);
            % 设置ttext的x坐标位置
            %ytext_xp=(xt(1))*ones(1,length(tempLabelStr))-2.4;%-2.4
            %ytext_xp=(xt(1))*ones(1,length(tempLabelStr))-6.5+4;%-2.4
            %ytext_xp=(xt(1)+temp_ticks_sampleLastT(1))*ones(1,length(tempLabelStr))-12.5;%-2.4            
            ytext_xp=(xt(1)+temp_ticks_sampleLastT(1))*ones(1,length(tempLabelStr))-8.5;%-2.4
            % 设置ttext的y坐标位置
            % ytext_yp=yt;
            ytext_yp = nan(1,2);
            ytext_yp(1) = temp_choiceIndex(2)/2 - 20;
            ytext_yp(2) = temp_choiceIndex(2) + temptempIndex_middle_choiceOffload + 20;
            ytext_yp = ytext_yp + 0;%0.5
            
            %text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',90,'fontsize',6.5);%8
            text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',0,'fontsize',6.5);%8            
            set(gca,'yticklabel','');
            
            text(ytext_xp(1),temp_choiceIndex(2),'---------','HorizontalAlignment','center','rotation',0,'fontsize',8);%8            
            
        
        
        xtickangle(0);
        
        %text(temp_ticks_baselinePeriod(end),size(cellID_F_dff_baselinePeriod_choice_collaped,1),'/','fontsize',12);
        %text(temp_ticks_decisionPeriodA(1),size(cellID_F_dff_baselinePeriod_choice_collaped,1),'/','fontsize',12);
                
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        %ylabel(sprintf('Trials'), 'FontSize', 10, 'FontWeight', 'normal');                       
        
        c = colorbar('fontsize',6.5);%9
        %c.Position = c.Position+[0.06 0.04 -0.025 -0.05];
        c.Position = c.Position+[0.01 0.06 -0.038 -0.075];
        
        
        %if if_colormap_loadEnhanced == 1
        %    load('parula_enhanced');
        %    colormap(parula_enhanced);
        %elseif if_colormap_loadEnhanced == 0
        %    colormap parula
        %end
        %colormap parula
        
        %temp1 = gray;
        %temp1 = temp1(end:-1:1,:);
        %colormap(temp1);
        
        colormap(coolwarm());
        
        %temp1 = coolwarm(300);
        %temp2 = ((300-256)/2)+1;
        %temp3 = temp1(temp2:temp2+255,:);
        %colormap(temp3);
        
        % %% nexttile
        %% Plot choice tuning
        
        nexttile([1 4])  
        
        h_line = [];
        
        %% plot baselinePeriod
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceMemory_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceOffload_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        %% plot sampleLastT
        x = temp_ticks_sampleLastT(1):temp_ticks_sampleLastT(end);
        y = cellID_F_dff_sampleLastT_choiceMemory_mean;
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_sampleLastT_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_sampleLastT(1):temp_ticks_sampleLastT(end);
        y = cellID_F_dff_sampleLastT_choiceOffload_mean;
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_sampleLastT_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on

        
        %% plot decisionPeriodA
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(2);
        y = cellID_F_dff_decisionPeriodA_choiceMemory_mean(decisionPeriodA_interval(1):decisionPeriodA_interval(2));
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceMemory_sem(decisionPeriodA_interval(1):decisionPeriodA_interval(2));
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriodA(2)+1:temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceMemory_mean(decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceMemory_sem(decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on        
        
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(2);
        y = cellID_F_dff_decisionPeriodA_choiceOffload_mean(decisionPeriodA_interval(1):decisionPeriodA_interval(2));
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceOffload_sem(decisionPeriodA_interval(1):decisionPeriodA_interval(2));
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriodA(2)+1:temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceOffload_mean(decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceOffload_sem(decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on        
        
        
        %% plot decisionPeriod
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceMemory_mean;
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceOffload_mean;
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        %% Others
        
        %le = legend(h_line,'ChoiceMemory','ChoiceOffload','Location','northeast','fontsize',8);%10
        %le = legend(h_line,'Choice-memory','Choice-offload','Location','northwest','fontsize',6.5,'NumColumns',2);%10                        
        le = legend(h_line,'Memory','Offload','Location','northwest','fontsize',6.5,'NumColumns',1);%10
        le.ItemTokenSize = ones(1,3)*10;
        le.Color = backgounrdColor;
        legend('boxoff');
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriod(end)+xlim_padSize2]);
        %xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]);
        xlim([temp_ticks_sampleLastT(1) temp_ticks_decisionPeriodA(end)]);                        
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        %xlim([1 temp_ticks_decisionPeriodA(end)+1]);        
        %xlim([1 temp_ticks_baselinePeriod(end)+4]);
        ylim([min_baselineDecisionPeriod,...
            max_baselineDecisionPeriod+(max_baselineDecisionPeriod-min_baselineDecisionPeriod)*0.35]);
        
        set(gca,'xticklabel',[])
        set(gca,'xtick',[])
        
        %xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        %xticks([temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        xticks([temp_ticks_sampleLastT(1) temp_ticks_decisionPeriodA]);

        %xticklabels({'PreFixation','Fixation','','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'Fixation','','Delay1','ChoiceCue','','','Decision',''});        
        %xticklabels({'Fixation','T1','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'Fixation','LastT','Delay1','Decision','','','Decision',''});        
        %xticklabels({'','','','Delay1','ChoiceCue','','','',''});
        %xticklabels({'LastT','Delay1','Decision'});        
        xticklabels({'LastT','Delay-on','Decision'});        
        
        xtickangle(0);
        
        %yticks([0 0.8]);
        
        %text(temp_ticks_baselinePeriod(end),min_baselineDecisionPeriod,'/','fontsize',12);
        %text(temp_ticks_decisionPeriodA(1),min_baselineDecisionPeriod,'/','fontsize',12);
        %
        %     text(temp_ticks_decisionPeriodA(end),min_baselineDecisionPeriod,'/','fontsize',15);
        %     text(temp_ticks_decisionPeriod(1),min_baselineDecisionPeriod,'/','fontsize',15);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        ylabel(sprintf('dF/F'), 'FontSize', 8, 'FontWeight', 'normal');
        
        
    end    
    
    if if_plot_multiPanelE == 1
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67*1.04 127*1.05*0.95*0.93*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        %set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67*1.04*0.8*1.01 127*1.05*0.95*0.93*1.5*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        %set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67*1.04*0.8*1.01*0.92 127*1.05*0.95*0.93*1.5*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67*1.04*0.8*1.01*0.92 127*1.05*0.95*0.93*1.5*0.9*0.9*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
                
        
        t = tiledlayout(2,5,'TileSpacing','none','Padding','tight'); %#ok<*NASGU>        
                
        min_baselineDecisionPeriod_raw = min([...
            cellID_F_dff_baselinePeriod_choiceMemory_mean-cellID_F_dff_baselinePeriod_choiceMemory_sem,...
            cellID_F_dff_baselinePeriod_choiceOffload_mean-cellID_F_dff_baselinePeriod_choiceOffload_sem,...
            cellID_F_dff_decisionPeriodA_choiceMemory_mean-cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
            cellID_F_dff_decisionPeriodA_choiceOffload_mean-cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        
        max_baselineDecisionPeriod_raw = max([...
            cellID_F_dff_baselinePeriod_choiceMemory_mean+cellID_F_dff_baselinePeriod_choiceMemory_sem,...
            cellID_F_dff_baselinePeriod_choiceOffload_mean+cellID_F_dff_baselinePeriod_choiceOffload_sem,...
            cellID_F_dff_decisionPeriodA_choiceMemory_mean+cellID_F_dff_decisionPeriodA_choiceMemory_sem,...            
            cellID_F_dff_decisionPeriodA_choiceOffload_mean+cellID_F_dff_decisionPeriodA_choiceOffload_sem]);
        
        min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        %min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.4;
        max_baselineDecisionPeriod = max_baselineDecisionPeriod_raw+(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        
        if cellIndex_suite2p_temptempRaw == 124 % 0510A
            max_baselineDecisionPeriod = 1.45;
            min_baselineDecisionPeriod = -0.1;
        end
        
        if cellIndex_suite2p_temptempRaw == 1 % 0510A
            max_baselineDecisionPeriod = 0.6;
        end        
        
        %% Plot trials raster plot        
        nexttile([1 4])                  
        
        temptempMax = max_baselineDecisionPeriod*1.5;

        temptempMin = 0;
        
        if cellIndex_suite2p_temptempRaw == 92
            %temptempMax = 0.4;
            
            temptempMin = 0;
        end
        
        if cellIndex_suite2p_temptempRaw == 20
            %temptempMax = 1;
            %temptempMin = 0.4;
        end        
        
        if cellIndex_suite2p_temptempRaw == 124 % 0510A
            temptempMax = 1.05;
        end
        if cellIndex_suite2p_temptempRaw == 1 % 0510A
            temptempMax = 0.42;
        end        
            
        if_sortedChoiceTuning = 1;
        
        cellID_F_dff_decisionPeriodA_choice_collaped;
        temp_choiceIndex;
        
        cellID_F_dff_decisionPeriodA_choice_collaped_mean = mean(cellID_F_dff_decisionPeriodA_choice_collaped,2);
        cellID_F_dff_decisionPeriodA_choice_collaped_mean_choiceMemory = cellID_F_dff_decisionPeriodA_choice_collaped_mean(1:temp_choiceIndex(2));
        cellID_F_dff_decisionPeriodA_choice_collaped_mean_choiceOffload = cellID_F_dff_decisionPeriodA_choice_collaped_mean((temp_choiceIndex(2)+1):end);
        
        [M_choiceMemory,I_choiceMemory] = sort(cellID_F_dff_decisionPeriodA_choice_collaped_mean_choiceMemory,'descend');
        [M_choiceOffload,I_choiceOffload] = sort(cellID_F_dff_decisionPeriodA_choice_collaped_mean_choiceOffload,'descend');
        I_choiceOffload = I_choiceOffload + temp_choiceIndex(2);
        
        I_sorted = [I_choiceMemory;I_choiceOffload];
        
        if if_sortedChoiceTuning == 1
            temp_I = I_sorted;
        elseif if_sortedChoiceTuning == 0
            temp_I = 1:size(cellID_F_dff_decisionPeriodA_choice_collaped,1);
        end
        
        
        %% plot baselinePeriod
        x = [temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end)];
        y = [1,size(cellID_F_dff_baselinePeriod_choice_collaped,1)];
        C = cellID_F_dff_baselinePeriod_choice_collaped(temp_I,:);
        imagesc(x,y,C,[temptempMin,temptempMax]);
        hold on
        
        %% plot decisionPeriodA
        x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)];
        y = [1,size(cellID_F_dff_decisionPeriodA_choice_collaped,1)];
        C = cellID_F_dff_decisionPeriodA_choice_collaped(temp_I,:);
        imagesc(x,y,C,[temptempMin,temptempMax]);
        hold on
        
        %% plot decisionPeriod
        x = [temp_ticks_decisionPeriod(1) temp_ticks_decisionPeriod(end)];
        y = [1,size(cellID_F_dff_decisionPeriod_choice_collaped,1)];
        C = cellID_F_dff_decisionPeriod_choice_collaped(temp_I,:);
        imagesc(x,y,C,[temptempMin,temptempMax]);
        hold on
                
        
        
        plot((temp_ticks_decisionPeriodA(2)+0.25)*[1 1],[1,size(cellID_F_dff_decisionPeriodA_choice_collaped,1)],...
            '-','LineWidth',1,'Color',[1 1 1]);%3,2
        hold on
           
        
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        %xlim([1 temp_ticks_decisionPeriodA(end)+1]);  
        %xlim([1 temp_ticks_baselinePeriod(end)+4]);  
        xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]);        
        ylim([1-1 size(cellID_F_dff_decisionPeriod_choice_collaped,1)+1]);                                
        
        set(gca,'xticklabel',[])
        %set(gca,'xtick',[])
        xticks([temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        
                
        yticks([]);
        
            tempLabelStr = string(1:2);
            tempLabelStr(1) = 'Memory';
            tempLabelStr(2) = 'Offload';
            ytl=string(tempLabelStr);
            % 设置ttext的x坐标位置
            %ytext_xp=(xt(1))*ones(1,length(tempLabelStr))-2.4;%-2.4
            %ytext_xp=(xt(1))*ones(1,length(tempLabelStr))-6.5+4;%-2.4
            ytext_xp=(xt(1)+temp_ticks_decisionPeriodA(1))*ones(1,length(tempLabelStr))-6.5-1.0-0.5;%-2.4            
            % 设置ttext的y坐标位置
            % ytext_yp=yt;
            ytext_yp = nan(1,2);
            ytext_yp(1) = temp_choiceIndex(2)/2 - 20;
            ytext_yp(2) = temp_choiceIndex(2) + temptempIndex_middle_choiceOffload + 20;
            ytext_yp = ytext_yp + 0;%0.5
            
            %text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',90,'fontsize',6.5);%8
            text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',0,'fontsize',6.5);%8            
            set(gca,'yticklabel','');
            
            text(ytext_xp(1),temp_choiceIndex(2),'---------','HorizontalAlignment','center','rotation',0,'fontsize',8);%8            
            
        
        
        xtickangle(0);
        
        %text(temp_ticks_baselinePeriod(end),size(cellID_F_dff_baselinePeriod_choice_collaped,1),'/','fontsize',12);
        %text(temp_ticks_decisionPeriodA(1),size(cellID_F_dff_baselinePeriod_choice_collaped,1),'/','fontsize',12);
                
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        %ylabel(sprintf('Trials'), 'FontSize', 10, 'FontWeight', 'normal');                       
        
        c = colorbar('fontsize',6.5);%9
        %c.Position = c.Position+[0.06 0.04 -0.025 -0.05];
        %c.Position = c.Position+[0.01 0.06 -0.038 -0.075];
        c.Position = c.Position+[0.00 0.06 -0.046 -0.075];
        
        if temptempMax > 1
            c.Ticks = [0 floor(temptempMax)];
        else
            c.Ticks = [0 floor(temptempMax*10)/10];
        end
            
        %c.Ticks = [temptempMin temptempMax];
        
        %if if_colormap_loadEnhanced == 1
        %    load('parula_enhanced');
        %    colormap(parula_enhanced);
        %elseif if_colormap_loadEnhanced == 0
        %    colormap parula
        %end
        %colormap parula
        
        colormap(coolwarm());
        
        %temp1 = coolwarm(300);
        %temp2 = ((300-256)/2)+1;
        %temp3 = temp1(temp2:temp2+255,:);
        %colormap(temp3);
        
        
        % %% nexttile
        %% Plot choice tuning
        
        nexttile([1 4])  
                
        h_line = [];
        
        %% plot baselinePeriod
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceMemory_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceOffload_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        %% plot decisionPeriodA
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(2);
        y = cellID_F_dff_decisionPeriodA_choiceMemory_mean(decisionPeriodA_interval(1):decisionPeriodA_interval(2));
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceMemory_sem(decisionPeriodA_interval(1):decisionPeriodA_interval(2));
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriodA(2)+1:temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceMemory_mean(decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceMemory_sem(decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on        
        
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(2);
        y = cellID_F_dff_decisionPeriodA_choiceOffload_mean(decisionPeriodA_interval(1):decisionPeriodA_interval(2));
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceOffload_sem(decisionPeriodA_interval(1):decisionPeriodA_interval(2));
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriodA(2)+1:temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceOffload_mean(decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceOffload_sem(decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on           
        
        
        %% plot decisionPeriod
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceMemory_mean;
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceMemory,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceMemory_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceOffload_mean;
        plot(x,y,'-','LineWidth',0.5,'Color',color_choiceOffload,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceOffload_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
        hold on
        
        
        %le = legend(h_line,'ChoiceMemory','ChoiceOffload','Location','northeast','fontsize',8);%10
        %le = legend(h_line,'Choice-memory','Choice-offload','Location','northwest','fontsize',6.5,'NumColumns',1);%10
        le = legend(h_line,'Memory','Offload','Location','northwest','fontsize',6.5,'NumColumns',2);%10        
        le.ItemTokenSize = ones(1,3)*10;
        le.Color = backgounrdColor;
        legend('boxoff');
        
        set(gca,'linewidth',1.5)
        set(gca,'color',backgounrdColor);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriod(end)+xlim_padSize2]);
        xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]);
        %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+1]);
        %xlim([1 temp_ticks_decisionPeriodA(end)+1]);        
        %xlim([1 temp_ticks_baselinePeriod(end)+4]);
        ylim([min_baselineDecisionPeriod,...
            max_baselineDecisionPeriod+(max_baselineDecisionPeriod-min_baselineDecisionPeriod)*0.35]);
        
        set(gca,'xticklabel',[])
        set(gca,'xtick',[])
        
        %xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        %xticks([temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end) temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
        xticks([temp_ticks_decisionPeriodA]);
        
        %xticklabels({'PreFixation','Fixation','','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'Fixation','','Delay1','ChoiceCue','','','Decision',''});        
        %xticklabels({'Fixation','T1','Delay1','ChoiceCue','','','Decision',''});
        %xticklabels({'','','','Delay1','ChoiceCue','','','',''});
        %xticklabels({'Delay1','Decision'});
        xticklabels({'Delay-on','Decision'});
        
        xtickangle(0);
        
        %yticks([0 0.8]);
        
        %text(temp_ticks_baselinePeriod(end),min_baselineDecisionPeriod,'/','fontsize',12);
        %text(temp_ticks_decisionPeriodA(1),min_baselineDecisionPeriod,'/','fontsize',12);
        %
        %     text(temp_ticks_decisionPeriodA(end),min_baselineDecisionPeriod,'/','fontsize',15);
        %     text(temp_ticks_decisionPeriod(1),min_baselineDecisionPeriod,'/','fontsize',15);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        ylabel(sprintf('dF/F'), 'FontSize', 8, 'FontWeight', 'normal');
        
        
    end
    
    
end

%% End