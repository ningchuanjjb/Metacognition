% Chuan's 25th script (20251214)
if_plot_memoryMetaMismatch = 1;

% temp_id2;

% cellIndex_suite2p_temptemp = 16;%
% cellIndex_suite2p_temptemp = 17;%


% cellIndex_suite2p_temptempRaw = 92;%

% memoryPrecisionHigh_choiceMemoryHigh
cellIndex_suite2p_temptempRaw = 20;%
% cellIndex_suite2p_temptempRaw = 25;%25
% cellIndex_suite2p_temptempRaw = 354;%
% cellIndex_suite2p_temptempRaw = 290;%

% memoryPrecisionLow_choiceMemoryHigh
% cellIndex_suite2p_temptempRaw = 18;%18
% cellIndex_suite2p_temptempRaw = 107;%107,19,404
% cellIndex_suite2p_temptempRaw = 404;%
% cellIndex_suite2p_temptempRaw = 107;%19

% memoryPrecisionHigh_choiceMemoryLow
% cellIndex_suite2p_temptempRaw = 673;%673,193
% cellIndex_suite2p_temptempRaw = 193;%673,193
% cellIndex_suite2p_temptempRaw = 745;%745
% cellIndex_suite2p_temptempRaw = 9;%9

temp_cellIndex_target = find(decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,2)==cellIndex_suite2p_temptempRaw);
cellIndex_suite2p_temptemp = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(temp_cellIndex_target,1);

% if_p_linear0_each1 = 1;
% 
% if if_p_linear0_each1 == 0
%     cellIndex_suite2p_temptemp = 16;%17
% elseif if_p_linear0_each1 == 1
%     cellIndex_suite2p_temptemp = 17;%17
% end

% color_choiceMemoryHigh = [127,201,127]/255;%[0.5,1,0.5]
% color_choiceMemoryLow = [0,0.5,0];%[0,0.5,0]
% color_choiceOffloadHigh = [0.5,0,0];%[1,0.5,0.5]
% color_choiceOffloadLow = [1,0.5,0.5];%[0.5,0,0]

color_choiceMemoryHigh = [128,205,193]/255; %[146,197,222]/255
color_choiceMemoryLow = [1,133,113]/255; %[5,113,176]/255
color_choiceOffloadHigh = [166,97,26]/255; %[202,0,32]/255
color_choiceOffloadLow = [223,194,125]/255; %[244,165,130]/255


periodSkipInterval = 2;%3
temp_ticks_baselinePeriod = baselinePeriod_interval;
temp_ticks_decisionPeriodA = temp_ticks_baselinePeriod(end)+periodSkipInterval+decisionPeriodA_interval;
temp_ticks_decisionPeriod = temp_ticks_decisionPeriodA(end)+periodSkipInterval+decisionPeriod_interval;

% backgounrdColor = [1 1 1]*0.825;%0.875
backgounrdColor = [1 1 1];

temp_id2 = find(cellIndex_suite2p == cellIndex_suite2p_temptemp);

if if_plot_memoryMetaMismatch == 1
    close all
    
    %fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptemp,currentSession(18:22)),'NumberTitle','off');
    fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptemp,FOVName_currentFOV2),'NumberTitle','off');    
    %set(gcf,'Position',[35+0 42+0 530 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[35+0 42+0 340 130]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point    
    %set(gcf,'Position',[35+0 42+0 340 130*1.15]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
    set(gcf,'Position',[35+0 42+0 340 130*1.15*0.9*0.95*1.01]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point            
    t = tiledlayout(1,2,'TileSpacing','tight','Padding','tight'); %#ok<*NASGU>
    
    %t.Title.String = sprintf('suite2pID%d,%s',...
    %    cellIndex_suite2p_temptemp,FOVName_currentFOV2);
    t.Title.String = sprintf('suite2pBID%d,%s',...
        cellIndex_suite2p_temptempRaw,FOVName_currentFOV2);
    t.Title.FontSize = 10;
    t.Title.Interpreter = 'none';
       
    %nexttile([1 2])
    nexttile
    
    %% memoryMetaMismatch, trial average in decision-making
    if if_plot_memoryMetaMismatch == 1
        memoryMetaMismatch;
        
        tempTrialBoolIndex_lowMatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionLow_choiceOffloadLow;
        tempTrialBoolIndex_highMatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh;
        tempTrialBoolIndex_overMismatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh;
        tempTrialBoolIndex_underMismatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow;
%         tempTrialBoolIndex_lowMatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionLow_choiceOffloadLow(1036:end);
%         tempTrialBoolIndex_highMatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh(1036:end);
%         tempTrialBoolIndex_overMismatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh(1036:end);
%         tempTrialBoolIndex_underMismatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow(1036:end);
        
        
        
        %% baselinePeriod
        cellID_F_dff_baselinePeriod_choiceMemoryLow = squeeze(F_dff_baselinePeriod(temp_id2,tempTrialBoolIndex_overMismatch,:));
        cellID_F_dff_baselinePeriod_choiceMemoryHigh = squeeze(F_dff_baselinePeriod(temp_id2,tempTrialBoolIndex_highMatch,:));
        cellID_F_dff_baselinePeriod_choiceOffloadLow = squeeze(F_dff_baselinePeriod(temp_id2,tempTrialBoolIndex_lowMatch,:));
        cellID_F_dff_baselinePeriod_choiceOffloadHigh = squeeze(F_dff_baselinePeriod(temp_id2,tempTrialBoolIndex_underMismatch,:));
        
        cellID_F_dff_baselinePeriod_choiceMemoryLow_mean = mean(cellID_F_dff_baselinePeriod_choiceMemoryLow,1);
        cellID_F_dff_baselinePeriod_choiceMemoryHigh_mean = mean(cellID_F_dff_baselinePeriod_choiceMemoryHigh,1);
        cellID_F_dff_baselinePeriod_choiceOffloadLow_mean = mean(cellID_F_dff_baselinePeriod_choiceOffloadLow,1);
        cellID_F_dff_baselinePeriod_choiceOffloadHigh_mean = mean(cellID_F_dff_baselinePeriod_choiceOffloadHigh,1);
        
        cellID_F_dff_baselinePeriod_choiceMemoryLow_sem = std(cellID_F_dff_baselinePeriod_choiceMemoryLow,1)...
            ./sqrt(size(cellID_F_dff_baselinePeriod_choiceMemoryLow,1));
        cellID_F_dff_baselinePeriod_choiceMemoryHigh_sem = std(cellID_F_dff_baselinePeriod_choiceMemoryHigh,1)...
            ./sqrt(size(cellID_F_dff_baselinePeriod_choiceMemoryHigh,1));
        cellID_F_dff_baselinePeriod_choiceOffloadLow_sem = std(cellID_F_dff_baselinePeriod_choiceOffloadLow,1)...
            ./sqrt(size(cellID_F_dff_baselinePeriod_choiceOffloadLow,1));
        cellID_F_dff_baselinePeriod_choiceOffloadHigh_sem = std(cellID_F_dff_baselinePeriod_choiceOffloadHigh,1)...
            ./sqrt(size(cellID_F_dff_baselinePeriod_choiceOffloadHigh,1));
        
        
        %% decisionPeriodA
        cellID_F_dff_decisionPeriodA_choiceMemoryLow = squeeze(F_dff_decisionPeriodA(temp_id2,tempTrialBoolIndex_overMismatch,:));
        cellID_F_dff_decisionPeriodA_choiceMemoryHigh = squeeze(F_dff_decisionPeriodA(temp_id2,tempTrialBoolIndex_highMatch,:));
        cellID_F_dff_decisionPeriodA_choiceOffloadLow = squeeze(F_dff_decisionPeriodA(temp_id2,tempTrialBoolIndex_lowMatch,:));
        cellID_F_dff_decisionPeriodA_choiceOffloadHigh = squeeze(F_dff_decisionPeriodA(temp_id2,tempTrialBoolIndex_underMismatch,:));
        
        cellID_F_dff_decisionPeriodA_choiceMemoryLow_mean = mean(cellID_F_dff_decisionPeriodA_choiceMemoryLow,1);
        cellID_F_dff_decisionPeriodA_choiceMemoryHigh_mean = mean(cellID_F_dff_decisionPeriodA_choiceMemoryHigh,1);
        cellID_F_dff_decisionPeriodA_choiceOffloadLow_mean = mean(cellID_F_dff_decisionPeriodA_choiceOffloadLow,1);
        cellID_F_dff_decisionPeriodA_choiceOffloadHigh_mean = mean(cellID_F_dff_decisionPeriodA_choiceOffloadHigh,1);
        
        cellID_F_dff_decisionPeriodA_choiceMemoryLow_sem = std(cellID_F_dff_decisionPeriodA_choiceMemoryLow,1)...
            ./sqrt(size(cellID_F_dff_decisionPeriodA_choiceMemoryLow,1));
        cellID_F_dff_decisionPeriodA_choiceMemoryHigh_sem = std(cellID_F_dff_decisionPeriodA_choiceMemoryHigh,1)...
            ./sqrt(size(cellID_F_dff_decisionPeriodA_choiceMemoryHigh,1));
        cellID_F_dff_decisionPeriodA_choiceOffloadLow_sem = std(cellID_F_dff_decisionPeriodA_choiceOffloadLow,1)...
            ./sqrt(size(cellID_F_dff_decisionPeriodA_choiceOffloadLow,1));
        cellID_F_dff_decisionPeriodA_choiceOffloadHigh_sem = std(cellID_F_dff_decisionPeriodA_choiceOffloadHigh,1)...
            ./sqrt(size(cellID_F_dff_decisionPeriodA_choiceOffloadHigh,1));
        
        
        %% decisionPeriod
        cellID_F_dff_decisionPeriod_choiceMemoryLow = squeeze(F_dff_decisionPeriod(temp_id2,tempTrialBoolIndex_overMismatch,:));
        cellID_F_dff_decisionPeriod_choiceMemoryHigh = squeeze(F_dff_decisionPeriod(temp_id2,tempTrialBoolIndex_highMatch,:));
        cellID_F_dff_decisionPeriod_choiceOffloadLow = squeeze(F_dff_decisionPeriod(temp_id2,tempTrialBoolIndex_lowMatch,:));
        cellID_F_dff_decisionPeriod_choiceOffloadHigh = squeeze(F_dff_decisionPeriod(temp_id2,tempTrialBoolIndex_underMismatch,:));
        
        cellID_F_dff_decisionPeriod_choiceMemoryLow_mean = mean(cellID_F_dff_decisionPeriod_choiceMemoryLow,1);
        cellID_F_dff_decisionPeriod_choiceMemoryHigh_mean = mean(cellID_F_dff_decisionPeriod_choiceMemoryHigh,1);
        cellID_F_dff_decisionPeriod_choiceOffloadLow_mean = mean(cellID_F_dff_decisionPeriod_choiceOffloadLow,1);
        cellID_F_dff_decisionPeriod_choiceOffloadHigh_mean = mean(cellID_F_dff_decisionPeriod_choiceOffloadHigh,1);
        
        cellID_F_dff_decisionPeriod_choiceMemoryLow_sem = std(cellID_F_dff_decisionPeriod_choiceMemoryLow,1)...
            ./sqrt(size(cellID_F_dff_decisionPeriod_choiceMemoryLow,1));
        cellID_F_dff_decisionPeriod_choiceMemoryHigh_sem = std(cellID_F_dff_decisionPeriod_choiceMemoryHigh,1)...
            ./sqrt(size(cellID_F_dff_decisionPeriod_choiceMemoryHigh,1));
        cellID_F_dff_decisionPeriod_choiceOffloadLow_sem = std(cellID_F_dff_decisionPeriod_choiceOffloadLow,1)...
            ./sqrt(size(cellID_F_dff_decisionPeriod_choiceOffloadLow,1));
        cellID_F_dff_decisionPeriod_choiceOffloadHigh_sem = std(cellID_F_dff_decisionPeriod_choiceOffloadHigh,1)...
            ./sqrt(size(cellID_F_dff_decisionPeriod_choiceOffloadHigh,1));
        
        
        %% Others        
        %min_baselineDecisionPeriod;
        %max_baselineDecisionPeriod;
        
%         min_baselineDecisionPeriod = min([min_baselineDecisionPeriod ...
%             cellID_F_dff_baselinePeriod_choiceMemoryLow_mean,cellID_F_dff_baselinePeriod_choiceMemoryHigh_mean,...
%             cellID_F_dff_baselinePeriod_choiceOffloadLow_mean,cellID_F_dff_baselinePeriod_choiceOffloadHigh_mean,...
%             cellID_F_dff_decisionPeriodA_choiceMemoryLow_mean,cellID_F_dff_decisionPeriodA_choiceMemoryHigh_mean,...
%             cellID_F_dff_decisionPeriodA_choiceOffloadLow_mean,cellID_F_dff_decisionPeriodA_choiceOffloadHigh_mean,...
%             cellID_F_dff_decisionPeriod_choiceMemoryLow_mean,cellID_F_dff_decisionPeriod_choiceMemoryHigh_mean,...
%             cellID_F_dff_decisionPeriod_choiceOffloadLow_mean,cellID_F_dff_decisionPeriod_choiceOffloadHigh_mean]);
        min_baselineDecisionPeriod_raw = min([...
            cellID_F_dff_decisionPeriodA_choiceMemoryLow_mean-cellID_F_dff_decisionPeriodA_choiceMemoryLow_sem,...
            cellID_F_dff_decisionPeriodA_choiceMemoryHigh_mean-cellID_F_dff_decisionPeriodA_choiceMemoryHigh_sem,...
            cellID_F_dff_decisionPeriodA_choiceOffloadLow_mean-cellID_F_dff_decisionPeriodA_choiceOffloadLow_sem,...
            cellID_F_dff_decisionPeriodA_choiceOffloadHigh_mean-cellID_F_dff_decisionPeriodA_choiceOffloadHigh_sem]);
        
%         max_baselineDecisionPeriod = max([max_baselineDecisionPeriod ...
%             cellID_F_dff_baselinePeriod_choiceMemoryLow_mean,cellID_F_dff_baselinePeriod_choiceMemoryHigh_mean,...
%             cellID_F_dff_baselinePeriod_choiceOffloadLow_mean,cellID_F_dff_baselinePeriod_choiceOffloadHigh_mean,...
%             cellID_F_dff_decisionPeriodA_choiceMemoryLow_mean,cellID_F_dff_decisionPeriodA_choiceMemoryHigh_mean,...
%             cellID_F_dff_decisionPeriodA_choiceOffloadLow_mean,cellID_F_dff_decisionPeriodA_choiceOffloadHigh_mean,...
%             cellID_F_dff_decisionPeriod_choiceMemoryLow_mean,cellID_F_dff_decisionPeriod_choiceMemoryHigh_mean,...
%             cellID_F_dff_decisionPeriod_choiceOffloadLow_mean,cellID_F_dff_decisionPeriod_choiceOffloadHigh_mean]);
        max_baselineDecisionPeriod_raw = max([...
            cellID_F_dff_decisionPeriodA_choiceMemoryLow_mean+cellID_F_dff_decisionPeriodA_choiceMemoryLow_sem,...
            cellID_F_dff_decisionPeriodA_choiceMemoryHigh_mean+cellID_F_dff_decisionPeriodA_choiceMemoryHigh_sem,...
            cellID_F_dff_decisionPeriodA_choiceOffloadLow_mean+cellID_F_dff_decisionPeriodA_choiceOffloadLow_sem,...
            cellID_F_dff_decisionPeriodA_choiceOffloadHigh_mean+cellID_F_dff_decisionPeriodA_choiceOffloadHigh_sem]);
        
        %min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        %max_baselineDecisionPeriod = max_baselineDecisionPeriod_raw+(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.1;
        min_baselineDecisionPeriod = min_baselineDecisionPeriod_raw-(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.02;
        max_baselineDecisionPeriod = max_baselineDecisionPeriod_raw+(max_baselineDecisionPeriod_raw-min_baselineDecisionPeriod_raw)*0.02;
        
        h_line = [];
        
        %% plot baselinePeriod
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceMemoryHigh_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',color_choiceMemoryHigh,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceMemoryHigh_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryHigh,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceOffloadLow_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',color_choiceOffloadLow,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceOffloadLow_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
        
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceMemoryLow_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',color_choiceMemoryLow,'LineStyle','-');];
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceMemoryLow_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
        
        
        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
        y = cellID_F_dff_baselinePeriod_choiceOffloadHigh_mean;
        h_line = [h_line plot(x,y,'-','LineWidth',1,'Color',color_choiceOffloadHigh,'LineStyle','-');];                
        hold on
        y_sem = cellID_F_dff_baselinePeriod_choiceOffloadHigh_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadHigh,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
        
        %% plot decisionPeriodA
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceMemoryHigh_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceMemoryHigh,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceMemoryHigh_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryHigh,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceOffloadLow_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceOffloadLow,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceOffloadLow_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceMemoryLow_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceMemoryLow,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceMemoryLow_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
                
        
        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
        y = cellID_F_dff_decisionPeriodA_choiceOffloadHigh_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceOffloadHigh,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriodA_choiceOffloadHigh_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadHigh,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
        
        %% plot decisionPeriod
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceMemoryHigh_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceMemoryHigh,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceMemoryHigh_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryHigh,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceOffloadLow_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceOffloadLow,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceOffloadLow_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceMemoryLow_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceMemoryLow,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceMemoryLow_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on                
        
        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
        y = cellID_F_dff_decisionPeriod_choiceOffloadHigh_mean;
        plot(x,y,'-','LineWidth',1,'Color',color_choiceOffloadHigh,'LineStyle','-');
        hold on
        y_sem = cellID_F_dff_decisionPeriod_choiceOffloadHigh_sem;
        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadHigh,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
    end
    
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
    
%     le = legend(h_line,'highMatch','lowMatch','overMismatch','underMismatch','Location','northeast','fontsize',8);%10
%     le.ItemTokenSize = ones(1,3)*10;
%     le.Color = backgounrdColor;
    
    set(gca,'linewidth',1.5)
    set(gca,'color',backgounrdColor);
    %xlim([1-xlim_padSize1 temp_ticks_decisionPeriod(end)+xlim_padSize2]);
    xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)+1]);    
    ylim([min_baselineDecisionPeriod,max_baselineDecisionPeriod]);
    
    set(gca,'xticklabel',[])
    set(gca,'xtick',[])
    
    xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
    xticklabels({'PreFixation','Fixation','','Delay1','ChoiceCue','','','Decision',''});
    %xticklabels({'','','','Delay1','ChoiceCue','','','',''});
    
    xtickangle(0);
    
%     text(temp_ticks_baselinePeriod(end),min_baselineDecisionPeriod,'/','fontsize',15);
%     text(temp_ticks_decisionPeriodA(1),min_baselineDecisionPeriod,'/','fontsize',15);
%     
%     text(temp_ticks_decisionPeriodA(end),min_baselineDecisionPeriod,'/','fontsize',15);
%     text(temp_ticks_decisionPeriod(1),min_baselineDecisionPeriod,'/','fontsize',15);
    
    set(gca, 'FontSize', 10)
    set(gca,'box','off');% 取消右、上边框
    ylabel(sprintf('dF/F'), 'FontSize', 10, 'FontWeight', 'normal');
    
    
    %% bar plot
    nexttile           
    
    cellID_F_dff_decisionPeriodA_choiceMemoryHigh;
    cellID_F_dff_decisionPeriodA_choiceMemoryLow;    
    cellID_F_dff_decisionPeriodA_choiceOffloadHigh;
    cellID_F_dff_decisionPeriodA_choiceOffloadLow;    
    
    cellID_F_dff_decisionABin_choiceMemoryHigh = ...
        mean(cellID_F_dff_decisionPeriodA_choiceMemoryHigh(:,1:decisionPeriodA_interval(2)),2);    
    cellID_F_dff_decisionABin_choiceMemoryLow = ...
        mean(cellID_F_dff_decisionPeriodA_choiceMemoryLow(:,1:decisionPeriodA_interval(2)),2);
    cellID_F_dff_decisionABin_choiceOffloadHigh = ...
        mean(cellID_F_dff_decisionPeriodA_choiceOffloadHigh(:,1:decisionPeriodA_interval(2)),2);
    cellID_F_dff_decisionABin_choiceOffloadLow = ...
        mean(cellID_F_dff_decisionPeriodA_choiceOffloadLow(:,1:decisionPeriodA_interval(2)),2);
    
    cellID_F_dff_decisionABin_choiceMemoryHigh;
    cellID_F_dff_decisionABin_choiceMemoryLow;    
    cellID_F_dff_decisionABin_choiceOffloadHigh;
    cellID_F_dff_decisionABin_choiceOffloadLow;    
    
%     [~,temp_p12,~,~] = ttest2(cellID_F_dff_decisionABin_choiceMemoryLow,cellID_F_dff_decisionABin_choiceMemoryHigh);
%     [~,temp_p13,~,~] = ttest2(cellID_F_dff_decisionABin_choiceMemoryLow,cellID_F_dff_decisionABin_choiceOffloadLow);
%     [~,temp_p14,~,~] = ttest2(cellID_F_dff_decisionABin_choiceMemoryLow,cellID_F_dff_decisionABin_choiceOffloadHigh);
%     [~,temp_p23,~,~] = ttest2(cellID_F_dff_decisionABin_choiceMemoryHigh,cellID_F_dff_decisionABin_choiceOffloadLow);
%     [~,temp_p24,~,~] = ttest2(cellID_F_dff_decisionABin_choiceMemoryHigh,cellID_F_dff_decisionABin_choiceOffloadHigh);
%     [~,temp_p34,~,~] = ttest2(cellID_F_dff_decisionABin_choiceOffloadLow,cellID_F_dff_decisionABin_choiceOffloadHigh);
         
    
x = [];
y = [];
for tempi=1:4
    if tempi==1
        temp_dff = cellID_F_dff_decisionABin_choiceMemoryHigh;
    elseif tempi==2
        temp_dff = cellID_F_dff_decisionABin_choiceMemoryLow;
    elseif tempi==3
        temp_dff = cellID_F_dff_decisionABin_choiceOffloadHigh;
    elseif tempi==4        
        temp_dff = cellID_F_dff_decisionABin_choiceOffloadLow;
    end
    x = [x;tempi*ones(size(temp_dff,1),1)];%#ok<*AGROW>    
    y = [y;temp_dff];
end
temp_mdl = fitglm(x,y,'linear');
temp_p_linear = temp_mdl.Coefficients.pValue(2);
    
    
%     temp_p12 = p12;
%     temp_p13 = p13;
%     temp_p14 = p14;
%     temp_p23 = p23;
%     temp_p24 = p24;
%     temp_p34 = p34;
    
    temp_1 = cellID_F_dff_decisionABin_choiceMemoryHigh;
    temp_2 = cellID_F_dff_decisionABin_choiceOffloadHigh;
    temp_3 = cellID_F_dff_decisionABin_choiceOffloadLow;
    temp_4 = cellID_F_dff_decisionABin_choiceMemoryLow;

    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    temp3_SEM = std(temp_3)/sqrt(length(temp_3));
    temp4_SEM = std(temp_4)/sqrt(length(temp_4));
    
    temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM,mean(temp_3)-temp3_SEM,mean(temp_4)-temp4_SEM]);
    temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM,mean(temp_4)+temp4_SEM]);
    
    temp_y_max12 = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
    temp_y_max13 = max([mean(temp_1)+temp1_SEM,mean(temp_3)+temp3_SEM]);
    temp_y_max14 = max([mean(temp_1)+temp1_SEM,mean(temp_4)+temp4_SEM]);    
    temp_y_max23 = max([mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM]);
    temp_y_max24 = max([mean(temp_2)+temp2_SEM,mean(temp_4)+temp4_SEM]);
    temp_y_max34 = max([mean(temp_3)+temp3_SEM,mean(temp_4)+temp4_SEM]);
    
    
    %temp_bar = bar([0 1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3) mean(temp_4)], ...
    %    'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1);
    temp_bar = bar([0 1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3) mean(temp_4)], ...
       'FaceColor','flat');
        
   temp_bar.CData(1,:) = color_choiceMemoryHigh;
   temp_bar.CData(2,:) = color_choiceOffloadHigh;
   temp_bar.CData(3,:) = color_choiceOffloadLow;
   temp_bar.CData(4,:) = color_choiceMemoryLow;

    
    hold on
    
    %errorbar([0 1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3) mean(temp_4)],...
    %    [temp1_SEM temp2_SEM temp3_SEM temp4_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 6);
    
    errorbar([0], [mean(temp_1)],[temp1_SEM],'.','Color',color_choiceMemoryHigh*0.65,...
        'LineWidth',2,'CapSize',6);
    hold on
    errorbar([1], [mean(temp_2)],[temp2_SEM],'.','Color',color_choiceOffloadHigh*0.65,...
        'LineWidth',2,'CapSize',6);
    hold on
    errorbar([2], [mean(temp_3)],[temp3_SEM],'.','Color',color_choiceOffloadLow*0.65,...
        'LineWidth',2,'CapSize',6);
    hold on
    errorbar([3], [mean(temp_4)],[temp4_SEM],'.','Color',color_choiceMemoryLow*0.65,...
        'LineWidth',2,'CapSize',6); %#ok<*NBRAK>
    hold on

%     if if_p_linear0_each1 == 0
%         tempTxt = sprintf('');
%         if temp_p_linear < 0.001
%             tempTxt = sprintf('***');
%         elseif temp_p_linear < 0.01
%             tempTxt = sprintf('**');
%         elseif temp_p_linear < 0.05
%             tempTxt = sprintf('*');
%         end
%         
%         text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.1,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
%             'HorizontalAlignment','center');
%         hold on
%     end
    
    
%     tempTxt = sprintf('');
%     if temp_p12 < 0.001
%         tempTxt = sprintf('***');
%     elseif temp_p12 < 0.01
%         tempTxt = sprintf('**');
%     elseif temp_p12 < 0.05
%         tempTxt = sprintf('*');
%     end
% 
%     text(0.5,(temp_y_max12-temp_y_min)*1.11,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
%         'HorizontalAlignment','center');
%     plot([0 1],(temp_y_max12-temp_y_min)*1.10*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
%     hold on
%     
%     if temp_p13 < 0.001
%         tempTxt = sprintf('***');
%     elseif temp_p13 < 0.01
%         tempTxt = sprintf('**');
%     elseif temp_p13 < 0.05
%         tempTxt = sprintf('*');
%     end
%     % 0.4660 0.6740 0.1880
%     %text(1,(temp_y_max13-temp_y_min)*1.05,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',20,'FontWeight','bold',...
%     %    'HorizontalAlignment','center');
%     %plot([0 2],(temp_y_max13-temp_y_min)*1.04*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
%     text(1,(temp_y_max13-temp_y_min)*1.04,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
%        'HorizontalAlignment','center');
%     plot([0 2],(temp_y_max13-temp_y_min)*1.03*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);    
%     hold on
%     
%     if temp_p23 < 0.001
%         tempTxt = sprintf('***');
%     elseif temp_p23 < 0.01
%         tempTxt = sprintf('**');
%     elseif temp_p23 < 0.05
%         tempTxt = sprintf('*');
%     end
%     % 0 0.4470 0.7410
%     text(1.5,(temp_y_max23-temp_y_min)*1.05,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
%         'HorizontalAlignment','center');
%     plot([1 2],(temp_y_max23-temp_y_min)*1.04*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
%     hold on
    
    
    set(gca,'linewidth',1.5)
    
    %ylim([0 1.1]);
%     if if_p_linear0_each1 == 0
%         ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
%     elseif if_p_linear0_each1 == 1
%         ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.5]);
%     end    
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 10) %14
    %set(gca,'XTickLabel', ["overMismatch";"highMatch";"lowMatch";"underMismatch"],'FontSize', 12);%给坐标加标签
    set(gca,'XTickLabel','');
    %set(gca,'YTick',0:0.2:1,'FontSize', 12);%给坐标加标签    
    set(gca,'box','off');% 取消右、上边框
    ylabel('dF/F', 'FontSize', 10, 'FontWeight', 'normal');    
    %temp_title = title(sprintf(' Location distribution \n correlation (with model) '),'fontsize',12);        
    %temp_title.Interpreter = 'none';    
    
end