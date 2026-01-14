%% Initialization
% clear
close all
home;% To scroll down in command window

currentSession = '329Recording_20231227A_YRK';

if_plotTuning_fixation1_sample2_delay3 = 3;

if_plotSelectiveCell = 0;
temp_id_suite2p = 464;%83-->417-->127-->464

plotRoiNum = 30;

if_load = 0;
if_compute = 0;

significantThreshold = 0.05;
significantThreshold_lowBound = -1;


if_plot = 1;


fixationMarkerValue = 16;
sampleMarkerValue = 20;
delayMarkerValue = 40;
rewardMarkerValue = 90;


fixationPeriodDuration = 16;%frames
samplePeriodDuration = 7;%frames
delayPeriodDuration = 48;%frames
rewardPeriodDuration = 90;%frames


numFrames = 6;


multi_rgbColor = ...
    [228,26,28;
    55,126,184;
    77,175,74;
    152,78,163;
    255,127,0;
    [255,255,51].*0.8]/255;
backgounrdColor = [1 1 1];


t0 = tic;

targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)
targetPATH_norm = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';

rawData_path = 'R:\twoPhotonRawData\ToBeProcessed';
markerExtractTestName_2p='fun_markerExtract_2p';
markerExtractTestName_PTB='fun_markerExtract_PTB';
markerDownSamplingTestName='fun_markerDownSampling';
markerParse_trialLevelTestName='fun_markerParse_trialLevel';


currentSession_path = [rawData_path '\' currentSession];


monkeyLogic_fullName = [rawData_path,'\',currentSession,'\','231227Len1.mat'];

output_shortPath = 'D:\twoPhotonData_motionCorrected';
temp_currentSession_path = [output_shortPath '\' currentSession];
temp_if_max0_min1 = 0;


MAT_file=dir([temp_currentSession_path,'\Result*']);
if isempty(MAT_file) == 1
    datestr_now30 = datestr(now,30);
    output_path = [temp_currentSession_path,'\Result',datestr_now30];
    mkdir(output_path);
else
    output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
end

%% Load dff
if if_load == 1
    cd(targetPATH)
    
    path_plane = [output_path,'\plane0'];
    fileName_Fall = 'Fall.mat';
    fileName_iscell = 'iscell.npy';
    fullFileName_Fall = [path_plane,'\',fileName_Fall];
    fullFileName_iscell = [path_plane,'\',fileName_iscell];
    
    load(fullFileName_Fall,'F_dff_raw');
    iscell = readNPY(fullFileName_iscell);
    
    cellIndex = find(iscell(:,1)==1);
    cellIndex_suite2p = cellIndex - 1;
    
    
    F_dff = F_dff_raw(iscell(:,1)==1, :);
    clear F_dff_raw
    
    s = load(fullFileName_Fall,'stat');
    roi_stats_raw = s.stat;
    temp_cellIndex = find(iscell(:,1)==1);
    roi_stats = roi_stats_raw(temp_cellIndex);
    
    
    F_dff_copy2 = nan(size(F_dff,1),size(F_dff,2)*2);
    for tempi=1:size(F_dff,1)
        temp_dff = F_dff(tempi,:);
        
        temp_dff_copy2 = repmat(temp_dff,2,1);
        temp_dff_copy2 = reshape(temp_dff_copy2,1,[]);
        
        F_dff_copy2(tempi,:) = temp_dff_copy2;
    end
    F_dff = F_dff_copy2;
    
    
    windowSize = 5;%3-->5
    x = F_dff;
    F_dff = smoothdata(x,2,'gaussian',windowSize);
    
    compensateImageDelay = 3;%1-->8-->3, to compensate imaging delay
    F_dff(:,1:end-compensateImageDelay) = F_dff(:,1+compensateImageDelay:end); % to compensate imaging delay
    
    
    temp_if_max0_min1 = 0;
    template_path = autoGetFileName_general('maxProjection*.tif', output_path,temp_if_max0_min1);
    template = double(loadtiff(template_path));
end

%% Load marker and behavior
temp_load = load(monkeyLogic_fullName);
load_monkeyLogic = temp_load.data;

trialNum = length(load_monkeyLogic);

behaviorData_cell = cell(trialNum,1);
behaviorMarker_cell = cell(trialNum,1);
twoPhotonScanVoltage_cell = cell(trialNum,1);
trialErrorData = nan(trialNum,1);
sequence = nan(trialNum,1);
blockType = nan(trialNum,1);
for tempi=1:length(load_monkeyLogic)
    temp_struct = load_monkeyLogic(tempi);
    behaviorData_cell{tempi} = temp_struct.UserVars;
    behaviorMarker_cell{tempi} = temp_struct.BehavioralCodes;
    twoPhotonScanVoltage_cell{tempi} = temp_struct.AnalogData.General.Gen2;
    trialErrorData(tempi) = temp_struct.TrialError;
    sequence(tempi) = temp_struct.UserVars.SequencePos(1);
    %sequence(tempi) = temp_struct.UserVars.Sequence(1);
    blockType(tempi) = temp_struct.Block;
end
isCorrect = trialErrorData == 0;
locationBlockBoolIndex = blockType==2;
% locationBlockBoolIndex = blockType==1;

locationBlockCorrectBoolIndex = locationBlockBoolIndex & isCorrect;


%% markerTime_merged_NICard
markerTime_fixation = nan(trialNum,1);
markerTime_sample = nan(trialNum,1);
markerTime_delay = nan(trialNum,1);
markerTime_reward = nan(trialNum,1);
for tempi=1:trialNum
    %temp_behaviorData = behaviorData_cell{tempi};
    temp_behaviorMarker = behaviorMarker_cell{tempi};
    %temp_twoPhotonScanVoltage = twoPhotonScanVoltage_cell{tempi};
    
    temptempIndex = find(temp_behaviorMarker.CodeNumbers==fixationMarkerValue);
    if ~isempty(temptempIndex)
        markerTime_fixation(tempi) = temp_behaviorMarker.CodeTimes(temptempIndex);
    end
    
    temptempIndex = find(temp_behaviorMarker.CodeNumbers==sampleMarkerValue);
    if ~isempty(temptempIndex)
        markerTime_sample(tempi) = temp_behaviorMarker.CodeTimes(temptempIndex);
    end
    
    temptempIndex = find(temp_behaviorMarker.CodeNumbers==delayMarkerValue);
    if ~isempty(temptempIndex)
        markerTime_delay(tempi) = temp_behaviorMarker.CodeTimes(temptempIndex);
    end
    
    temptempIndex = find(temp_behaviorMarker.CodeNumbers==rewardMarkerValue);
    if ~isempty(temptempIndex)
        markerTime_reward(tempi) = temp_behaviorMarker.CodeTimes(temptempIndex);
    end    
end
markerTime_merged = [markerTime_fixation,markerTime_sample,markerTime_delay,markerTime_reward];
markerTime_merged = round(markerTime_merged);

markerTime_merged_NICard = markerTime_merged;

%% twoPhotonMarker_cell
twoPhotonMarker_cell = cell(trialNum,1);
for tempi=1:trialNum
    temp_twoPhotonScanVoltage = twoPhotonScanVoltage_cell{tempi};
    
    temptempBoolIndex = false(length(temp_twoPhotonScanVoltage),1);
    for tempj=2:length(temp_twoPhotonScanVoltage)
        if temp_twoPhotonScanVoltage(tempj-1)<0 && temp_twoPhotonScanVoltage(tempj)>0
            if abs(temp_twoPhotonScanVoltage(tempj-1)) < abs(temp_twoPhotonScanVoltage(tempj))
                temptempBoolIndex(tempj-1) = true;
            else
                temptempBoolIndex(tempj) = true;
            end
        end
    end
    temp_twoPhotonMarker_raw = find(temptempBoolIndex==true);
    %if round(length(temp_twoPhotonMarker_raw)/2)*2 ~= length(temp_twoPhotonMarker_raw)
    %    temp_twoPhotonMarker_raw = [temp_twoPhotonMarker_raw;temp_twoPhotonMarker_raw(end)]; %#ok<*AGROW>
    %end
    %temp1 = reshape(temp_twoPhotonMarker_raw,2,[]);
    %temp2 = round(mean(temp1,1));
    %temp_twoPhotonMarker = temp2';
    %twoPhotonMarker_cell{tempi} = temp_twoPhotonMarker;
    
    twoPhotonMarker_cell{tempi} = temp_twoPhotonMarker_raw;
end

%% markerTime_merged_frame
twoPhotonMarker_cell;
markerTime_merged_frame = nan(size(markerTime_merged_NICard));
for tempi=1:trialNum
    temp_markerTime_merged_NICard = markerTime_merged_NICard(tempi,:);
    temp_twoPhotonMarker = twoPhotonMarker_cell{tempi};
    
    for tempj=1:length(temp_markerTime_merged_NICard)
        temptempIndex = find(temp_markerTime_merged_NICard(tempj)<temp_twoPhotonMarker,1);
        if isempty(temptempIndex)
            continue
        end
        if temptempIndex >= 2
            temp1 = abs(temp_twoPhotonMarker(temptempIndex)-temp_markerTime_merged_NICard(tempj));
            temp2 = abs(temp_twoPhotonMarker(temptempIndex-1)-temp_markerTime_merged_NICard(tempj));
            if temp1 > temp2
                temptempIndex = temptempIndex - 1;
            end
        end
        
        markerTime_merged_frame(tempi,tempj) = temptempIndex;
    end
end

twoPhotonMarker_cell;

cumulativeFrame = zeros(trialNum,1);
for tempi=2:trialNum
    cumulativeFrame(tempi) = cumulativeFrame(tempi-1) + length(twoPhotonMarker_cell{tempi-1});    
end

%%
F_dff;
markerTime_merged_frame;
sequence;
locationBlockCorrectBoolIndex;

fixationPeriodDuration;
samplePeriodDuration;
delayPeriodDuration;
rewardPeriodDuration;

roiNum = size(F_dff,1);

F_dff_fixationPeriod = nan(roiNum,trialNum,fixationPeriodDuration);
F_dff_samplePeriod = nan(roiNum,trialNum,samplePeriodDuration);
F_dff_delayPeriod = nan(roiNum,trialNum,delayPeriodDuration);
F_dff_rewardPeriod = nan(roiNum,trialNum,rewardPeriodDuration);

for tempi=1:trialNum
    
    if tempi >= 2
        
    end
    
    temp_markerTime = markerTime_merged_frame(tempi,:) + cumulativeFrame(tempi);
    
    if isnan(temp_markerTime(end))
        continue
    end
    
    temp_dff_fixation = F_dff(:,temp_markerTime(1):(temp_markerTime(1)+fixationPeriodDuration-1));
    temp_dff_sample = F_dff(:,temp_markerTime(2):(temp_markerTime(2)+samplePeriodDuration-1));
    temp_dff_delay = F_dff(:,temp_markerTime(3):(temp_markerTime(3)+delayPeriodDuration-1));
    
    rewardPeriodDurationA = round(rewardPeriodDuration/2);
    rewardPeriodDurationB = rewardPeriodDuration - rewardPeriodDurationA;
    
    temp_dff_reward = F_dff(:,(temp_markerTime(4)-rewardPeriodDurationA):(temp_markerTime(4)+rewardPeriodDurationB-1));
    
    F_dff_fixationPeriod(:,tempi,:) = temp_dff_fixation;
    F_dff_samplePeriod(:,tempi,:) = temp_dff_sample;
    F_dff_delayPeriod(:,tempi,:) = temp_dff_delay;
    F_dff_rewardPeriod(:,tempi,:) = temp_dff_reward;
    
end
F_dff_fixationPeriod;
F_dff_samplePeriod;
F_dff_delayPeriod;
F_dff_rewardPeriod;

%% Anova
if if_compute == 1
    frameAnova1Name_v = autoGetFunName_myScripts('frameAnova1', [targetPATH '\functions']);
    fun_frameAnova1 = str2func(frameAnova1Name_v);
    
    target_seqSet = {1,2,3,4,5,6};
    if_selective_seq0_loc1 = 1;
    
    % selectiveCellNum_fixationBin
    temp_dff = F_dff_fixationPeriod(:,locationBlockCorrectBoolIndex,:);
    trialLabel = sequence(locationBlockCorrectBoolIndex);
    [selectiveCellBoolIndex_fixationBin,selectiveCellNum_fixationBin] = ...
        fun_frameAnova1(temp_dff,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet,if_selective_seq0_loc1);
    selectiveCellIndex_suite2p_fixationBin = cellIndex_suite2p(selectiveCellBoolIndex_fixationBin);
    fprintf('selectiveCellNum_fixationBin = %d.\n',selectiveCellNum_fixationBin);
    
    % selectiveCellNum_sampleBin
    temp_dff = F_dff_samplePeriod(:,locationBlockCorrectBoolIndex,:);
    trialLabel = sequence(locationBlockCorrectBoolIndex);
    [selectiveCellBoolIndex_sampleBin,selectiveCellNum_sampleBin] = ...
        fun_frameAnova1(temp_dff,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet,if_selective_seq0_loc1);
    selectiveCellIndex_suite2p_sampleBin = cellIndex_suite2p(selectiveCellBoolIndex_sampleBin);
    fprintf('selectiveCellNum_sampleBin = %d.\n',selectiveCellNum_sampleBin);
    
    % selectiveCellNum_delayBin
    temp_dff = F_dff_delayPeriod(:,locationBlockCorrectBoolIndex,:);
    trialLabel = sequence(locationBlockCorrectBoolIndex);
    [selectiveCellBoolIndex_delayBin,selectiveCellNum_delayBin] = ...
        fun_frameAnova1(temp_dff,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet,if_selective_seq0_loc1);
    selectiveCellIndex_suite2p_delayBin = cellIndex_suite2p(selectiveCellBoolIndex_delayBin);
    fprintf('selectiveCellNum_delayBin = %d.\n',selectiveCellNum_delayBin);
    
end
selectiveCellIndex_suite2p_fixationBin;
selectiveCellIndex_suite2p_sampleBin;
selectiveCellIndex_suite2p_delayBin;

if if_plotSelectiveCell == 1
    if if_plotTuning_fixation1_sample2_delay3 == 1
        selectiveCellIndex_suite2p = selectiveCellIndex_suite2p_fixationBin;
    elseif if_plotTuning_fixation1_sample2_delay3 == 2
        selectiveCellIndex_suite2p = selectiveCellIndex_suite2p_sampleBin;
    elseif if_plotTuning_fixation1_sample2_delay3 == 3
        selectiveCellIndex_suite2p = selectiveCellIndex_suite2p_delayBin;
    end
    selectiveCellIndex_suite2p = selectiveCellIndex_suite2p(selectiveCellIndex_suite2p>=temp_id_suite2p);
elseif if_plotSelectiveCell == 0
    selectiveCellIndex_suite2p = temp_id_suite2p;
end

plotRoiNum = min(plotRoiNum,length(selectiveCellIndex_suite2p));
fprintf('plotRoiNum=%d.\n',plotRoiNum);



%% Plot
for loopCount=1:plotRoiNum
    
     
    cellIndex_suite2p_temptemp = selectiveCellIndex_suite2p(loopCount);
    
    fig1 = figure('Name',sprintf('Cell id %d',cellIndex_suite2p_temptemp),'NumberTitle','off');
    set(gcf,'Position',[450+0 202+0 512*1.2 200*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,1,'TileSpacing','tight','Padding','loose'); %#ok<*NASGU>
    
    t.Title.String = sprintf('suite2pID%d',cellIndex_suite2p_temptemp);
    t.Title.FontSize = 11;
    t.Title.Interpreter = 'none';
    
    
    nexttile
    
    
    temp_id2 = find(cellIndex_suite2p == cellIndex_suite2p_temptemp);
    
    cellID_F_dff_fixation = squeeze(F_dff_fixationPeriod(temp_id2,:,:));
    cellID_F_dff_sample = squeeze(F_dff_samplePeriod(temp_id2,:,:));
    cellID_F_dff_delay = squeeze(F_dff_delayPeriod(temp_id2,:,:));
    cellID_F_dff_reward = squeeze(F_dff_rewardPeriod(temp_id2,:,:));
    
    cellID_F_dff_fixation_location = cell(numFrames,1);
    cellID_F_dff_sample_location = cell(numFrames,1);
    cellID_F_dff_delay_location = cell(numFrames,1);
    cellID_F_dff_reward_location = cell(numFrames,1);    
    for tempi=1:numFrames
        temptempBoolIndex = locationBlockCorrectBoolIndex & (sequence==tempi);
        
        cellID_F_dff_fixation_location{tempi} = cellID_F_dff_fixation(temptempBoolIndex,:);
        cellID_F_dff_sample_location{tempi} = cellID_F_dff_sample(temptempBoolIndex,:);
        cellID_F_dff_delay_location{tempi} = cellID_F_dff_delay(temptempBoolIndex,:);
        cellID_F_dff_reward_location{tempi} = cellID_F_dff_reward(temptempBoolIndex,:);        
    end
    
    cellID_F_dff_fixation_location_mean = nan(numFrames,fixationPeriodDuration);
    cellID_F_dff_fixation_location_sem = nan(numFrames,fixationPeriodDuration);
    for tempi=1:numFrames
        temp_dff = cellID_F_dff_fixation_location{tempi};
        cellID_F_dff_fixation_location_mean(tempi,:) = mean(temp_dff,1,'omitnan');
        cellID_F_dff_fixation_location_sem(tempi,:) = std(temp_dff,1,1,'omitnan')./sqrt(size(temp_dff,1));
    end
    
    cellID_F_dff_sample_location_mean = nan(numFrames,samplePeriodDuration);
    cellID_F_dff_sample_location_sem = nan(numFrames,samplePeriodDuration);
    for tempi=1:numFrames
        temp_dff = cellID_F_dff_sample_location{tempi};
        cellID_F_dff_sample_location_mean(tempi,:) = mean(temp_dff,1,'omitnan');
        cellID_F_dff_sample_location_sem(tempi,:) = std(temp_dff,1,1,'omitnan')./sqrt(size(temp_dff,1));
    end
    
    cellID_F_dff_delay_location_mean = nan(numFrames,delayPeriodDuration);
    cellID_F_dff_delay_location_sem = nan(numFrames,delayPeriodDuration);
    for tempi=1:numFrames
        temp_dff = cellID_F_dff_delay_location{tempi};
        cellID_F_dff_delay_location_mean(tempi,:) = mean(temp_dff,1,'omitnan');
        cellID_F_dff_delay_location_sem(tempi,:) = std(temp_dff,1,1,'omitnan')./sqrt(size(temp_dff,1));
    end
    
    cellID_F_dff_reward_location_mean = nan(numFrames,rewardPeriodDuration);
    cellID_F_dff_reward_location_sem = nan(numFrames,rewardPeriodDuration);
    for tempi=1:numFrames
        temp_dff = cellID_F_dff_reward_location{tempi};
        cellID_F_dff_reward_location_mean(tempi,:) = mean(temp_dff,1,'omitnan');
        cellID_F_dff_reward_location_sem(tempi,:) = std(temp_dff,1,1,'omitnan')./sqrt(size(temp_dff,1));
    end
    
    
    
    temp_ticks_fixation = [1,fixationPeriodDuration];    
    temp_ticks_sample = temp_ticks_fixation(end)+[1,samplePeriodDuration];
    temp_ticks_delay = temp_ticks_sample(end)+[1,delayPeriodDuration];
    temp_ticks_reward = temp_ticks_delay(end)+[1,rewardPeriodDurationA,rewardPeriodDuration];
    
    %% Plot fixation
    h_line = [];
    for tempi=1:numFrames
        x = temp_ticks_fixation(1):temp_ticks_fixation(end);
        y = cellID_F_dff_fixation_location_mean(tempi,:);
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
        hold on
        y_sem = cellID_F_dff_fixation_location_sem(tempi,:);
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
    end
    
    %% Plot sample
    h_line = [];
    for tempi=1:numFrames
        x = temp_ticks_sample(1):temp_ticks_sample(end);
        y = cellID_F_dff_sample_location_mean(tempi,:);
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
        hold on
        y_sem = cellID_F_dff_sample_location_sem(tempi,:);
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
    end

    
    %% Plot delay
    for tempi=1:numFrames
        x = temp_ticks_delay(1):temp_ticks_delay(2);
        y = cellID_F_dff_delay_location_mean(tempi,1:length(x));
        plot(x,y,'-','LineWidth',1,'Color',multi_rgbColor(tempi,:));
        hold on
        y_sem = cellID_F_dff_delay_location_sem(tempi,1:length(x));
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
    end
    
    %% Plot reward
    for tempi=1:numFrames
        x = temp_ticks_reward(1):temp_ticks_reward(end);
        y = cellID_F_dff_reward_location_mean(tempi,1:length(x));
        plot(x,y,'-','LineWidth',1,'Color',multi_rgbColor(tempi,:));
        hold on
        y_sem = cellID_F_dff_reward_location_sem(tempi,1:length(x));
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
    end    
    
    
    le = legend(h_line,'1','2','3',...
        '4','5','6','Location','northwest','fontsize',9,'NumColumns',6);%10--.8
    le.ItemTokenSize = ones(1,6)*12;%10
    le.Color = backgounrdColor;
    a = 1;
    
    set(gca,'linewidth',1.5);
    set(gca,'color',backgounrdColor);
    set(gca,'xticklabel',[])
    set(gca,'xtick',[])
    xticks([temp_ticks_fixation(1:end-1) temp_ticks_sample(1:end-1) temp_ticks_delay(1:end-1) temp_ticks_reward(1:end-1)]);
    xticklabels({'Fixation','T1','Delay1','','Reward'});

    xtickangle(30);
    
    set(gca, 'FontSize', 10)
    set(gca,'box','off');% 取消右、上边框
    ylabel(sprintf('dF / F'), 'FontSize', 12);

    
end

fprintf('Time of the marker_extract is %.1f secs.\n',toc(t0));

cd(targetPATH)
%% End