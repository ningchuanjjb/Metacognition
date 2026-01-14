
if if_plot_memoryMetaMismatch == 1
    memoryMetaMismatch;
    
%     tempTrialBoolIndex_lowMatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionLow_choiceOffloadLow;
%     tempTrialBoolIndex_highMatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh;
%     tempTrialBoolIndex_overMismatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh;
%     tempTrialBoolIndex_underMismatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow;
    tempTrialBoolIndex_lowMatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionLow_choiceOffloadLow(1036:end);
    tempTrialBoolIndex_highMatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh(1036:end);
    tempTrialBoolIndex_overMismatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh(1036:end);
    tempTrialBoolIndex_underMismatch = memoryMetaMismatch.trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow(1036:end);
    
    
    
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
    
    color_choiceMemoryLow = [0,0.5,0];
    color_choiceMemoryHigh = [0.5,1,0.5];
    color_choiceOffloadLow = [0.5,0,0];
    color_choiceOffloadHigh = [1,0.5,0.5];
    
    
    min_baselineDecisionPeriod;
    max_baselineDecisionPeriod;
    
    min_baselineDecisionPeriod = min([min_baselineDecisionPeriod ...
        cellID_F_dff_baselinePeriod_choiceMemoryLow_mean,cellID_F_dff_baselinePeriod_choiceMemoryHigh_mean,...
        cellID_F_dff_baselinePeriod_choiceOffloadLow_mean,cellID_F_dff_baselinePeriod_choiceOffloadHigh_mean,...
        cellID_F_dff_decisionPeriodA_choiceMemoryLow_mean,cellID_F_dff_decisionPeriodA_choiceMemoryHigh_mean,...
        cellID_F_dff_decisionPeriodA_choiceOffloadLow_mean,cellID_F_dff_decisionPeriodA_choiceOffloadHigh_mean,...
        cellID_F_dff_decisionPeriod_choiceMemoryLow_mean,cellID_F_dff_decisionPeriod_choiceMemoryHigh_mean,...
        cellID_F_dff_decisionPeriod_choiceOffloadLow_mean,cellID_F_dff_decisionPeriod_choiceOffloadHigh_mean]);
    
    max_baselineDecisionPeriod = max([max_baselineDecisionPeriod ...
        cellID_F_dff_baselinePeriod_choiceMemoryLow_mean,cellID_F_dff_baselinePeriod_choiceMemoryHigh_mean,...
        cellID_F_dff_baselinePeriod_choiceOffloadLow_mean,cellID_F_dff_baselinePeriod_choiceOffloadHigh_mean,...
        cellID_F_dff_decisionPeriodA_choiceMemoryLow_mean,cellID_F_dff_decisionPeriodA_choiceMemoryHigh_mean,...
        cellID_F_dff_decisionPeriodA_choiceOffloadLow_mean,cellID_F_dff_decisionPeriodA_choiceOffloadHigh_mean,...
        cellID_F_dff_decisionPeriod_choiceMemoryLow_mean,cellID_F_dff_decisionPeriod_choiceMemoryHigh_mean,...
        cellID_F_dff_decisionPeriod_choiceOffloadLow_mean,cellID_F_dff_decisionPeriod_choiceOffloadHigh_mean]);
    
    
    %% plot baselinePeriod
    x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
    y = cellID_F_dff_baselinePeriod_choiceMemoryLow_mean;
    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceMemoryLow,'LineStyle','--');
    hold on
%     y_sem = cellID_F_dff_baselinePeriod_choiceMemoryLow_sem;
%     patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryLow,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
%     hold on
    
    x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
    y = cellID_F_dff_baselinePeriod_choiceMemoryHigh_mean;
    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceMemoryHigh,'LineStyle','--');
    hold on
%     y_sem = cellID_F_dff_baselinePeriod_choiceMemoryHigh_sem;
%     patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryHigh,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
%     hold on
    
    x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
    y = cellID_F_dff_baselinePeriod_choiceOffloadLow_mean;
    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceOffloadLow,'LineStyle','--');
    hold on
%     y_sem = cellID_F_dff_baselinePeriod_choiceOffloadLow_sem;
%     patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadLow,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
%     hold on
    
    x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
    y = cellID_F_dff_baselinePeriod_choiceOffloadHigh_mean;
    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceOffloadHigh,'LineStyle','--');
    hold on
%     y_sem = cellID_F_dff_baselinePeriod_choiceOffloadHigh_sem;
%     patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadHigh,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
%     hold on
    
    
    %% plot decisionPeriodA
    x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
    y = cellID_F_dff_decisionPeriodA_choiceMemoryLow_mean;
    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceMemoryLow,'LineStyle','--');
    hold on
%     y_sem = cellID_F_dff_decisionPeriodA_choiceMemoryLow_sem;
%     patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryLow,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
%     hold on
    
    x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
    y = cellID_F_dff_decisionPeriodA_choiceMemoryHigh_mean;
    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceMemoryHigh,'LineStyle','--');
    hold on
%     y_sem = cellID_F_dff_decisionPeriodA_choiceMemoryHigh_sem;
%     patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryHigh,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
%     hold on
    
    x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
    y = cellID_F_dff_decisionPeriodA_choiceOffloadLow_mean;
    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceOffloadLow,'LineStyle','--');
    hold on
%     y_sem = cellID_F_dff_decisionPeriodA_choiceOffloadLow_sem;
%     patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadLow,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
%     hold on
    
    x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
    y = cellID_F_dff_decisionPeriodA_choiceOffloadHigh_mean;
    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceOffloadHigh,'LineStyle','--');
    hold on
%     y_sem = cellID_F_dff_decisionPeriodA_choiceOffloadHigh_sem;
%     patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadHigh,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
%     hold on
    
    
    %% plot decisionPeriod
    x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
    y = cellID_F_dff_decisionPeriod_choiceMemoryLow_mean;
    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceMemoryLow,'LineStyle','--');
    hold on
%     y_sem = cellID_F_dff_decisionPeriod_choiceMemoryLow_sem;
%     patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryLow,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
%     hold on
    
    x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
    y = cellID_F_dff_decisionPeriod_choiceMemoryHigh_mean;
    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceMemoryHigh,'LineStyle','--');
    hold on
%     y_sem = cellID_F_dff_decisionPeriod_choiceMemoryHigh_sem;
%     patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryHigh,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
%     hold on
    
    x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
    y = cellID_F_dff_decisionPeriod_choiceOffloadLow_mean;
    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceOffloadLow,'LineStyle','--');
    hold on
%     y_sem = cellID_F_dff_decisionPeriod_choiceOffloadLow_sem;
%     patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadLow,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
%     hold on
    
    x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
    y = cellID_F_dff_decisionPeriod_choiceOffloadHigh_mean;
    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceOffloadHigh,'LineStyle','--');
    hold on
%     y_sem = cellID_F_dff_decisionPeriod_choiceOffloadHigh_sem;
%     patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadHigh,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
%     hold on
    
end


