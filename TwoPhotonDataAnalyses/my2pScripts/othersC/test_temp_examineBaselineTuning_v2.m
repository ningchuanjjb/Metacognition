close all

shuffleNum = 20;

if_compute = 1;

if_plot = 1;

if_baseline0_meta1 = 0;

%% Examine tuning for decision (baseline period)
if if_compute == 1
    temptemp_range = (baselinePeriod_interval(2)-8):baselinePeriod_interval(3);
    F_dff_baselineBin = mean(F_dff_baselinePeriod(:,:,temptemp_range),3);
    
    if if_baseline0_meta1 == 0
        x_raw = F_dff_baselineBin';
        
    elseif if_baseline0_meta1 == 1
        x_raw = F_dff_decisionBin1';
    end
    
    x = x_raw(choiceBoolIndex,:);
    y_raw = choiceMemoryBoolIndex';
    y = y_raw(choiceBoolIndex);
    
    
    x;
    y;
    
    temptempTrialNum = length(x);
    temptempTrialNumA = round(temptempTrialNum/2);
    temptempTrialNumB = temptempTrialNum - temptempTrialNumA;
    
    temptempBoolIndexA = false(temptempTrialNum,1);
    temptempBoolIndexA(1:temptempTrialNumA) = true;
    
    temptempBoolIndexB = false(temptempTrialNum,1);
    temptempBoolIndexB((temptempTrialNumA+1):end) = true;
    
    beta_A_shuffled = zeros(roiNum,shuffleNum);
    beta_B_shuffled = zeros(roiNum,shuffleNum);
    
    for tempShuffleIndex=1:shuffleNum
        
        temptempTrialIndex_shuffled = randperm(temptempTrialNum);
        
        temptempTrialIndexA_shuffled = temptempTrialIndex_shuffled(temptempBoolIndexA);
        temptempTrialIndexB_shuffled = temptempTrialIndex_shuffled(temptempBoolIndexB);
        
        x_A_shuffled = x(temptempTrialIndexA_shuffled,:);
        y_A_shuffled = y(temptempTrialIndexA_shuffled);
        
        x_B_shuffled = x(temptempTrialIndexB_shuffled,:);
        y_B_shuffled = y(temptempTrialIndexB_shuffled);
        
        parfor tempi=1:roiNum
            % for tempi=1:roiNum
            warning('off');
            temp_mdl = fitglm(x_A_shuffled(:,tempi),y_A_shuffled);
            beta_A_shuffled(tempi,tempShuffleIndex) = temp_mdl.Coefficients.Estimate(2);
            
            temp_mdl = fitglm(x_B_shuffled(:,tempi),y_B_shuffled);
            beta_B_shuffled(tempi,tempShuffleIndex) = temp_mdl.Coefficients.Estimate(2);
            warning('on');
        end
        
    end
    
    beta_A_shuffled;
    beta_B_shuffled;
    
    % beta_A_shuffled_median = median(beta_A_shuffled,2);
    % beta_B_shuffled_median = median(beta_B_shuffled,2);
    % [r_betaAB_shuffled,p_betaAB_shuffled] = corr(beta_A_shuffled_median,beta_B_shuffled_median);
    %
    % r_betaAB_shuffled
    
    r_betaAB_shuffled = nan(1,shuffleNum);
    p_betaAB_shuffled = nan(1,shuffleNum);
    for tempi=1:shuffleNum
        [r_betaAB_shuffled(tempi),p_betaAB_shuffled(tempi)] = corr(beta_A_shuffled(:,tempi),beta_B_shuffled(:,tempi));
    end
    
    r_betaAB_shuffled;
    r_betaAB_shuffled_median = median(r_betaAB_shuffled);
    
end

%% Plot
if if_plot == 1
    
    if if_baseline0_meta1 == 0
        tempStr = 'Baseline';
    elseif if_baseline0_meta1 == 1
        tempStr = 'Meta (delay1)';
    end
    
    
    fig = figure('Name','tuning property (seq-level)','NumberTitle','off');
    
    set(gcf,'Position',[400 50 445*0.8*1.02*0.965*0.88*1.5 379*0.7*0.9*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,2,'TileSpacing','Loose','Padding','Compact');
    
    
    %%
    nexttile
    
    x = beta_A_shuffled(:,1);
    y = beta_B_shuffled(:,1);
    
    [temp_r, temp_p] = corr(x,y);
    
    [x_min,x_max] = bounds(x);
    [y_min,y_max] = bounds(y);
    
    temp_MarkerFaceAlpha = 1;%0.5
    
    scatter(x,y,10,...
        'filled','MarkerFaceColor',[1 1 1]*0.3,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    
    
    tempTxt = 'p>0.05';
    if temp_p < 0.05
        tempTxt = 'p<0.05';
    end
    if temp_p < 0.01
        tempTxt = 'p<0.01';
    end
    if temp_p < 0.001
        tempTxt = 'p<0.001';
    end
    text(x_min+(x_max-x_min)*0.03,y_min+(y_max-y_min)*0.98,sprintf('r=%.3f',temp_r), 'fontsize',8,'FontWeight','normal');
    text(x_min+(x_max-x_min)*0.03,y_min+(y_max-y_min)*0.85,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
    
    
    set(gca,'linewidth',1.5)
    xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
    ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
    set(gca, 'FontSize', 10)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Beta A', 'FontSize', 10, 'FontWeight', 'normal');
    ylabel('Beta B', 'FontSize', 10, 'FontWeight', 'normal');
    temp_title = title(sprintf('%s, one shuffle',tempStr),'FontSize',10);
    temp_title.Interpreter = 'none';
    
    
    %% 
    nexttile
    
    temp_1 = r_betaAB_shuffled;
    temp_1_chanceLevel = 0;
    
    [~,temp_p] = ttest(temp_1,temp_1_chanceLevel);

    
    temp_y_min = min([temp_1 temp_1_chanceLevel]);
    temp_y_max = max(temp_1);
    
    
    temp_data = temp_1';
    temp_label = repmat({'A'},length(temp_1),1);
    
    temptemp_color1 = [1 1 1]*0.5;
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',{'A'});
    h(1).ViolinPlot.FaceAlpha = 0.1;
    hold on
    
    
    plot([0 2],temp_1_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    
    xlim([0 2]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    set(gca, 'FontSize', 10) %14,12
    %set(gca,'YTick',0:1,'FontSize', 10);%12    
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""]; %#ok<*NBRAK>
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Correlation', 'FontSize', 10, 'FontWeight', 'normal');    
    temp_title = title(sprintf('%s, %d shuffle',tempStr,shuffleNum),'fontsize',10);
    temp_title.Interpreter = 'none';
    
end