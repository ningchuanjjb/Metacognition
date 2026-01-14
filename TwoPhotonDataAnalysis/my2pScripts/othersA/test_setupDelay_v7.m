clear all
close all

plotRoiNum = 1;
plot_lengthFlag = 1;

if_load_rawData = 0;

currentSession = 'test0608';

rawData_path = 'R:\twoPhotonRawData\ToBeProcessed';
currentSession_path = [rawData_path '\' currentSession];

raw_fileName = 'Image_001_001.raw';

raw_fileFullName = [currentSession_path '\' raw_fileName];

%% Load raw data
if if_load_rawData == 1
    fid = fopen(raw_fileFullName);
    FOV = [512,512];
    bitsize = 2;
    imsize = FOV(1)*FOV(2)*bitsize;                                                   % Bit size of single frame
    current_seek = ftell(fid);
    fseek(fid, 0, 1);
    file_length = ftell(fid);
    fseek(fid, current_seek, -1);
    T = file_length/imsize;
    T = round(T);
    sizY = [FOV,T];
    Y = fread(fid,512*512*T,'uint16=>single',0,'l')';
    Y = reshape(Y,[512 512 T]);
    Y = pagetranspose(Y);
    
    Y_mean = squeeze(mean(Y,[1,2]))';
    save([currentSession_path,'\','Y_mean.mat'],'Y_mean');
    fclose(fid); 
else
    load([currentSession_path,'\','Y_mean.mat'],'Y_mean');
    Y_mean;
end

Y_mean;
% plot(Y_mean);
F_dff = [Y_mean; Y_mean];


windowSize = 1;%3-->5-->7

% F_dff = smoothdata(F_dff,2,'gaussian',windowSize);
% F_dff = smoothdata(F_dff,2,'sgolay',windowSize,'Degree',5);
% F_dff = smoothdata(F_dff,2,'rloess',windowSize);
% F_dff = halfGaussianFilter_v1(F_dff,windowSize);

%% Load behavior and marker
output_shortPath = 'D:\twoPhotonData_motionCorrected';
temp_currentSession_path = [output_shortPath '\' currentSession];
temp_if_max0_min1 = 0;
output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);

load([output_path,'\','markerParse_trialLevel.mat'],'markerParse_trialLevel');

PTB_name = '.mat';
PTB_fullName = autoGetMarkerFileName(['2*',PTB_name], currentSession_path);

basic_para = []; %#ok<*NASGU>
trial_para = [];
load(PTB_fullName,'basic_para','trial_para');

trialIndex_length1 = [];
for tempi=1:length(markerParse_trialLevel)
    temp_seq_length = length(markerParse_trialLevel{tempi}.currentSequence);
    if temp_seq_length == 1
        trialIndex_length1 = [trialIndex_length1 tempi]; %#ok<*AGROW>
    end
end

sequence_length1 = zeros(1,length(trialIndex_length1));
for tempi=1:length(trialIndex_length1)
    temp_trialIndex = trialIndex_length1(tempi);
    sequence_length1(tempi) = markerParse_trialLevel{temp_trialIndex}.currentSequence;
end

numFrames = 6;
pointKindsNum = 4;
target_seqSet = get_target_seqSet_v1(numFrames,pointKindsNum);
numSeq = zeros(1,pointKindsNum);
for tempi=1:pointKindsNum
    numSeq(tempi) = length(target_seqSet{tempi});    
end

trialIndex_length1_seq = cell(numFrames,1);
for tempi=1:numFrames
    temp_boolIndex = (sequence_length1 == tempi);
    temp_index = trialIndex_length1(temp_boolIndex);
    trialIndex_length1_seq{tempi} = temp_index;
end
a = 1;
trialIndex_length1_location = trialIndex_length1_seq;

%% Extract length1
%DELAY 1, DELAY 2 = 1500~1700
% Period 1: prefixation --> TRIALID --> TARGET 1 ITEM x ON --> DELAY 1 ON --> layout disappear--> DELAY 1 ON + 1500 ms
% Period 2: SELECTING AND DELAY 2 ON - 500 ms--> SELECTING AND DELAY 2 ON --> SELECTING AND DELAY 2 ON + 1500 ms
% Period 3: SUBMIT - 600 ms --> SUBMIT --> SUBMIT + 600 ms
%F_dff
%F_dff_length1_period1
%F_dff_length1_period2

length1_period1_frameIndex = zeros(length(trialIndex_length1),6);%4
length1_period2_frameIndex = zeros(length(trialIndex_length1),3);
length1_period3_frameIndex = zeros(length(trialIndex_length1),3);
for tempi=1:length(trialIndex_length1)
    temp_trialIndex = trialIndex_length1(tempi);
    
    temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(1:3);
    %temp_frameIndex = [temp_frameIndex(1:2) temp_frameIndex(2)+6 temp_frameIndex(end)];
    %length1_period1_frameIndex(tempi,:) = [temp_frameIndex temp_frameIndex(end)+45];
    %length1_period1_frameIndex(tempi,:) = [temp_frameIndex(1)-21 temp_frameIndex temp_frameIndex(end)+45];
    length1_period1_frameIndex(tempi,:) = [temp_frameIndex(1)-21 temp_frameIndex temp_frameIndex(end)+12 temp_frameIndex(end)+45];
    
    temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(5);
    length1_period2_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+45];
    
    temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(end);    
    length1_period3_frameIndex(tempi,:) = [temp_frameIndex-18 temp_frameIndex temp_frameIndex+18];
end
% length1_period1_interval = length1_period1_frameIndex(1,:)-length1_period1_frameIndex(1,1)+1;
% length1_period2_interval = length1_period2_frameIndex(1,:)-length1_period2_frameIndex(1,1)+1;
% length1_period3_interval = length1_period3_frameIndex(1,:)-length1_period3_frameIndex(1,1)+1;

% length1_period1_interval = median(length1_period1_frameIndex-length1_period1_frameIndex(:,1)+1,1);
% length1_period2_interval = median(length1_period2_frameIndex-length1_period2_frameIndex(:,1)+1,1);
% length1_period3_interval = median(length1_period3_frameIndex-length1_period3_frameIndex(:,1)+1,1);

length1_period1_interval = max(length1_period1_frameIndex-length1_period1_frameIndex(:,1)+1,[],1);
length1_period2_interval = max(length1_period2_frameIndex-length1_period2_frameIndex(:,1)+1,[],1);
length1_period3_interval = max(length1_period3_frameIndex-length1_period3_frameIndex(:,1)+1,[],1);

% align to fixation
% [length1_period1_frameRangeMin,~] = bounds(length1_period1_frameIndex,2);
% temp_frameRange = repmat(length1_period1_frameRangeMin,1,length1_period1_interval(end)-length1_period1_interval(1)+1);
% length1_period1_frameRange = temp_frameRange + ((length1_period1_interval(1):length1_period1_interval(end))-1);

% align to T1
length1_period1_frameRangeMin = length1_period1_frameIndex(:,3);
temp_frameRange = repmat(length1_period1_frameRangeMin,1,length1_period1_interval(end)-length1_period1_interval(1)+1);
length1_period1_frameRange = temp_frameRange + ((length1_period1_interval(1):length1_period1_interval(end))-length1_period1_interval(3));


% b1 = reshape(F_dff(2,length1_period1_frameRange(3,:)'),size(length1_period1_frameRange(3,:)'));
temp_dff = reshape(F_dff(:,length1_period1_frameRange'),[size(F_dff,1),size(length1_period1_frameRange')]);
F_dff_length1_period1 = permute(temp_dff,[1 3 2]);

a = 1;

[length1_period2_frameRangeMin,~] = bounds(length1_period2_frameIndex,2);
temp_frameRange = repmat(length1_period2_frameRangeMin,1,length1_period2_interval(end)-length1_period2_interval(1)+1);
length1_period2_frameRange = temp_frameRange + ((length1_period2_interval(1):length1_period2_interval(end))-1);

temp_dff = reshape(F_dff(:,length1_period2_frameRange'),[size(F_dff,1),size(length1_period2_frameRange')]);
F_dff_length1_period2 = permute(temp_dff,[1 3 2]);

[length1_period3_frameRangeMin,~] = bounds(length1_period3_frameIndex,2);
temp_frameRange = repmat(length1_period3_frameRangeMin,1,length1_period3_interval(end)-length1_period3_interval(1)+1);
length1_period3_frameRange = temp_frameRange + ((length1_period3_interval(1):length1_period3_interval(end))-1);

temp_dff = reshape(F_dff(:,length1_period3_frameRange'),[size(F_dff,1),size(length1_period3_frameRange')]);
F_dff_length1_period3 = permute(temp_dff,[1 3 2]);

%% Others
sequence_length1;
Y_mean;
markerParse_trialLevel;

lengthx_period1_interval = length1_period1_interval;
lengthx_period2_interval = length1_period2_interval;
lengthx_period3_interval = length1_period3_interval;
F_dff_lengthx_period1 = F_dff_length1_period1;
F_dff_lengthx_period2 = F_dff_length1_period2;
F_dff_lengthx_period3 = F_dff_length1_period3;
trialIndex_lengthx = trialIndex_length1;
trialIndex_lengthx_location = trialIndex_length1_location;
trialIndex_lengthx_memoryCorrect = trialIndex_length1;

%% Plot single cell
for tempIndex=1:plotRoiNum    
    %% compute cellID_F_dff
    temp_cellIndex_suite2p = tempIndex;
    
    
    cellID_F_dff_lengthx_period1 = squeeze(F_dff_lengthx_period1(tempIndex,:,:));
    cellID_F_dff_lengthx_period2 = squeeze(F_dff_lengthx_period2(tempIndex,:,:));
    cellID_F_dff_lengthx_period3 = squeeze(F_dff_lengthx_period3(tempIndex,:,:));
    
    cellID_F_dff_lengthx_period1_location = cell(1,numFrames);
    cellID_F_dff_lengthx_period2_location = cell(1,numFrames);
    cellID_F_dff_lengthx_period3_location = cell(1,numFrames);
    cellID_F_dff_lengthx_period1_location_mean = zeros(numFrames,size(cellID_F_dff_lengthx_period1,2));
    cellID_F_dff_lengthx_period2_location_mean = zeros(numFrames,size(cellID_F_dff_lengthx_period2,2));
    cellID_F_dff_lengthx_period3_location_mean = zeros(numFrames,size(cellID_F_dff_lengthx_period3,2));
    cellID_F_dff_lengthx_period1_location_sem = zeros(numFrames,size(cellID_F_dff_lengthx_period1,2));
    cellID_F_dff_lengthx_period2_location_sem = zeros(numFrames,size(cellID_F_dff_lengthx_period2,2));
    cellID_F_dff_lengthx_period3_location_sem = zeros(numFrames,size(cellID_F_dff_lengthx_period3,2));
    
    trialIndex_lengthx;
    %trialIndex_lengthx_seq;
    trialIndex_lengthx_location;
    
    for tempi=1:numFrames
        %[~,Locb] = ismember(trialIndex_lengthx_seq{tempi},trialIndex_lengthx_memoryCorrect);
        [~,Locb] = ismember(trialIndex_lengthx_location{tempi},trialIndex_lengthx_memoryCorrect);
        
        Locb(Locb==0) = [];
        cellID_F_dff_lengthx_period1_location{tempi} = cellID_F_dff_lengthx_period1(Locb,:);
        cellID_F_dff_lengthx_period2_location{tempi} = cellID_F_dff_lengthx_period2(Locb,:);
        cellID_F_dff_lengthx_period3_location{tempi} = cellID_F_dff_lengthx_period3(Locb,:);
        
        cellID_F_dff_lengthx_period1_location_mean(tempi,:) = mean(cellID_F_dff_lengthx_period1_location{tempi},1);
        cellID_F_dff_lengthx_period2_location_mean(tempi,:) = mean(cellID_F_dff_lengthx_period2_location{tempi},1);
        cellID_F_dff_lengthx_period3_location_mean(tempi,:) = mean(cellID_F_dff_lengthx_period3_location{tempi},1);
        
        cellID_F_dff_lengthx_period1_location_sem(tempi,:) = std(cellID_F_dff_lengthx_period1_location{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period1_location{tempi},1));
        cellID_F_dff_lengthx_period2_location_sem(tempi,:) = std(cellID_F_dff_lengthx_period2_location{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period2_location{tempi},1));
        cellID_F_dff_lengthx_period3_location_sem(tempi,:) = std(cellID_F_dff_lengthx_period3_location{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period3_location{tempi},1));
    end
    
    cellID_F_dff_lengthx_period1_location_collapsed = [];
    temp_lengthx_locationIndex = zeros(1, numFrames);
    for tempi=1:numFrames
        temp_dff = cellID_F_dff_lengthx_period1_location{tempi};
        temp_dff_mean = mean(temp_dff,2);
        %[~,I] = sort(temp_dff_mean,'descend');
        I = 1:size(temp_dff_mean,1);
        temp_dff_sorted = temp_dff(I,:);
        
        cellID_F_dff_lengthx_period1_location_collapsed = ...
            [cellID_F_dff_lengthx_period1_location_collapsed;temp_dff_sorted]; %#ok<*AGROW>
        temp_lengthx_locationIndex(tempi) = ...
            size(cellID_F_dff_lengthx_period1_location_collapsed,1) - ...
            size(cellID_F_dff_lengthx_period1_location{tempi},1) + 1;
    end
    
    cellID_F_dff_lengthx_period2_location_collapsed = [];
    for tempi=1:numFrames
        temp_dff1 = cellID_F_dff_lengthx_period1_location{tempi};
        temp_dff1_mean = mean(temp_dff1,2);
        %[~,I] = sort(temp_dff1_mean,'descend');
        I = 1:size(temp_dff1_mean,1);
        
        temp_dff2 = cellID_F_dff_lengthx_period2_location{tempi};
        %temp_dff2_mean = mean(temp_dff2,2);
        %[~,I] = sort(temp_dff2_mean,'descend');
        %temp_dff_sorted = temp_dff2(I,:);
        
        temp_dff2_sorted = temp_dff2(I,:);
        
        cellID_F_dff_lengthx_period2_location_collapsed = ...
            [cellID_F_dff_lengthx_period2_location_collapsed;temp_dff2_sorted];
    end
    
    cellID_F_dff_lengthx_period3_location_collapsed = [];
    for tempi=1:numFrames
        temp_dff1 = cellID_F_dff_lengthx_period1_location{tempi};
        temp_dff1_mean = mean(temp_dff1,2);
        %[~,I] = sort(temp_dff1_mean,'descend');
        I = 1:size(temp_dff1_mean,1);
        
        temp_dff2 = cellID_F_dff_lengthx_period3_location{tempi};
        temp_dff2_sorted = temp_dff2(I,:);
        cellID_F_dff_lengthx_period3_location_collapsed = ...
            [cellID_F_dff_lengthx_period3_location_collapsed;temp_dff2_sorted];
    end
    
    
    periodSkipInterval = 3;
    temp_ticks_period1 = lengthx_period1_interval;
    temp_ticks_period2 = lengthx_period1_interval(end)+periodSkipInterval+lengthx_period2_interval;
    temp_ticks_period3 = lengthx_period1_interval(end)+periodSkipInterval+lengthx_period2_interval(end)+periodSkipInterval+lengthx_period3_interval;
    
    xlim_padSize1 = 2;
    xlim_padSize2 = 20;%9-->15-->18-->20
    
    multi_rgbColor = ...
        [228,26,28;
        55,126,184;
        77,175,74;
        152,78,163;
        255,127,0;
        255,255,51]/255;
    
    backgounrdColor = [1 1 1]*0.825;%0.875
    
    [min_period1,max_period1] = bounds(cellID_F_dff_lengthx_period1_location_mean,'all');
    [min_period2,max_period2] = bounds(cellID_F_dff_lengthx_period2_location_mean,'all');
    [min_period3,max_period3] = bounds(cellID_F_dff_lengthx_period3_location_mean,'all');
    
    min_period12 = min(min_period1,min_period2) - max([cellID_F_dff_lengthx_period1_location_sem cellID_F_dff_lengthx_period2_location_sem],[],'all');
    max_period12 = max(max_period1,max_period2) + max([cellID_F_dff_lengthx_period1_location_sem cellID_F_dff_lengthx_period2_location_sem],[],'all');
    
    fig1 = figure('Name','Fig1','NumberTitle','off');
    %set(gcf,'Position',[35+260 35+0 1100 850]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[35+260 35+0 1100 1000]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[35+0 35+0 1630 1000]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
    
    %set(gcf,'color',backgounrdColor);
    
    %% subplot 1, raster plot in period1
    nexttile
    x = [temp_ticks_period1(1) temp_ticks_period1(end)];
    y = [1,size(cellID_F_dff_lengthx_period1_location_collapsed,1)];
    C = cellID_F_dff_lengthx_period1_location_collapsed;
    imagesc(x,y,C);
    hold on
    
    x = [temp_ticks_period2(1) temp_ticks_period2(end)];
    y = [1,size(cellID_F_dff_lengthx_period1_location_collapsed,1)];
    C = cellID_F_dff_lengthx_period2_location_collapsed;
    imagesc(x,y,C);
    hold on
    
    x = [temp_ticks_period3(1) temp_ticks_period3(end)];
    y = [1,size(cellID_F_dff_lengthx_period1_location_collapsed,1)];
    C = cellID_F_dff_lengthx_period3_location_collapsed;
    imagesc(x,y,C);
    hold on
    
    
    for tempi=1:length(lengthx_period1_interval)
        plot(temp_ticks_period1(tempi)*[1 1],[1 size(cellID_F_dff_lengthx_period1_location_collapsed,1)],...
            '-','LineWidth',1,'Color',[1 1 1])
        hold on
    end
    
    for tempi=1:length(lengthx_period2_interval)
        plot(temp_ticks_period2(tempi)*[1 1],[1 size(cellID_F_dff_lengthx_period2_location_collapsed,1)],...
            '-','LineWidth',1,'Color',[1 1 1]);
        hold on
    end
    
    for tempi=1:length(lengthx_period3_interval)
        plot(temp_ticks_period3(tempi)*[1 1],[1 size(cellID_F_dff_lengthx_period3_location_collapsed,1)],...
            '-','LineWidth',1,'Color',[1 1 1]);
        hold on
    end
    
    for tempi=1:numFrames
        plot(temp_ticks_period3(end)+[0.5 5],(temp_lengthx_locationIndex(tempi)-0.5)*[1 1],...
            '-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
        hold on
    end
    
    set(gca,'color',backgounrdColor);
    xlim([1-xlim_padSize1 temp_ticks_period3(end)+xlim_padSize2]);
    ylim([1-1 size(cellID_F_dff_lengthx_period1_location_collapsed,1)+1]);
    set(gca,'xticklabel',[])
    set(gca,'xtick',[])
    
    xticks([temp_ticks_period1 temp_ticks_period2 temp_ticks_period3]);
    if plot_lengthFlag == 1
        xticklabels({'PreFixation','Fixation','T1','Delay1','LayoutOff',' ',' ','Decision','','','Submit',''});
    elseif plot_lengthFlag == 2
        xticklabels({'PreFixation','Fixation','T1','T2','Delay1','LayoutOff',' ',' ','Decision','','','Submit',''});
    elseif plot_lengthFlag == 3
        xticklabels({'PreFixation','Fixation','T1','T2','T3','Delay1','LayoutOff',' ',' ','Decision','','','Submit',''});
    end
    
    text(temp_ticks_period1(end),size(cellID_F_dff_lengthx_period1_location_collapsed,1),'/','fontsize',15);
    text(temp_ticks_period2(1),size(cellID_F_dff_lengthx_period1_location_collapsed,1),'/','fontsize',15);
    
    text(temp_ticks_period2(end),size(cellID_F_dff_lengthx_period1_location_collapsed,1),'/','fontsize',15);
    text(temp_ticks_period3(1),size(cellID_F_dff_lengthx_period1_location_collapsed,1),'/','fontsize',15);
    
    set(gca, 'FontSize', 11)
    set(gca,'box','off');% 取消右、上边框
    ylabel(sprintf('Trials\ncellID=%d',temp_cellIndex_suite2p), 'FontSize', 14, 'FontWeight', 'bold');
    
    
    %% subplot 2, trial average in period1
    nexttile
    h_line = [];
    for tempi=1:numFrames
        x = temp_ticks_period1(1):temp_ticks_period1(end);
        y = cellID_F_dff_lengthx_period1_location_mean(tempi,:);
        h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:))];
        hold on
        y_sem = cellID_F_dff_lengthx_period1_location_sem(tempi,:);
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
        hold on
    end
    
    
    for tempi=1:numFrames
        x = temp_ticks_period2(1):temp_ticks_period2(end);
        y = cellID_F_dff_lengthx_period2_location_mean(tempi,:);
        plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
        hold on
        y_sem = cellID_F_dff_lengthx_period2_location_sem(tempi,:);
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
        hold on
    end
    
    for tempi=1:numFrames
        x = temp_ticks_period3(1):temp_ticks_period3(end);
        y = cellID_F_dff_lengthx_period3_location_mean(tempi,:);
        plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
        hold on
        y_sem = cellID_F_dff_lengthx_period3_location_sem(tempi,:);
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
        hold on
    end
    
    for tempi=1:length(lengthx_period1_interval)
        plot(temp_ticks_period1(tempi)*[1 1],[min_period12 max_period12],...
            '-','LineWidth',1,'Color',[0 0 0]);
        hold on
    end
    
    for tempi=1:length(lengthx_period2_interval)
        plot(temp_ticks_period2(tempi)*[1 1],[min_period12 max_period12],...
            '-','LineWidth',1,'Color',[0 0 0]);
        
        hold on
    end
    
    for tempi=1:length(lengthx_period3_interval)
        plot(temp_ticks_period3(tempi)*[1 1],[min_period12 max_period12],...
            '-','LineWidth',1,'Color',[0 0 0]);
        hold on
    end
    
    legend(h_line,'location1','location2','location3','location4','location5','location6','Location','northeast','fontsize',10)
    
    set(gca,'color',backgounrdColor);
    xlim([1-xlim_padSize1 temp_ticks_period3(end)+xlim_padSize2]);
    ylim([min_period12,max_period12]);
    set(gca,'xticklabel',[])
    set(gca,'xtick',[])
    
    xticks([temp_ticks_period1 temp_ticks_period2 temp_ticks_period3]);
    if plot_lengthFlag == 1
        xticklabels({'PreFixation','Fixation','T1','Delay1','LayoutOff',' ',' ','Decision','','','Submit',''});
    elseif plot_lengthFlag == 2
        xticklabels({'PreFixation','Fixation','T1','T2','Delay1','LayoutOff',' ',' ','Decision','','','Submit',''});
    elseif plot_lengthFlag == 3
        xticklabels({'PreFixation','Fixation','T1','T2','T3','Delay1','LayoutOff',' ',' ','Decision','','','Submit',''});
    end
    
    text(temp_ticks_period1(end),min_period12,'/','fontsize',15);
    text(temp_ticks_period2(1),min_period12,'/','fontsize',15);
    
    text(temp_ticks_period2(end),min_period12,'/','fontsize',15);
    text(temp_ticks_period3(1),min_period12,'/','fontsize',15);
    
    set(gca, 'FontSize', 11)
    set(gca,'box','off');% 取消右、上边框
    ylabel(sprintf('dF / F\ncellID=%d',temp_cellIndex_suite2p), 'FontSize', 14, 'FontWeight', 'bold');
    
end



%% End