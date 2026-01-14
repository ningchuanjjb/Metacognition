
topX = 10;
topX2 = 50;%50

threshold_r2 = 0.01;

temp_rProb_glm_output = rProb_glm_output_all.delay1;

r2 = temp_rProb_glm_output.r2_rProb_glm;

selectiveCellBoolIndex_rProb_glm = temp_rProb_glm_output.selectiveCellBoolIndex_rProb_glm;

%selectiveCellBoolIndex_rProb_glm = selectiveCellBoolIndex_rProb_glm & tempSelectiveCellBoolIndex_length1_currentFOV;

%selectiveCellBoolIndex_rProb_glm = selectiveCellBoolIndex_rProb_glm & (~tempSelectiveCellBoolIndex_length1_currentFOV);

r2(~selectiveCellBoolIndex_rProb_glm) = -1;

[r2_sorted,I_sorted] = sort(r2,'descend');

highCorrIndex_rProb = I_sorted(1:topX);
highCorrIndex_suite2p_rProb = cellIndex(highCorrIndex_rProb) - 1;

highCorrIndex_rProb;
highCorrIndex_suite2p_rProb;
r2;

temp_F_dff_decisionBin = temp_rProb_glm_output.temp_F_dff_mean;


temp_topX2 = min(topX2,sum(r2_sorted>threshold_r2));
temp_I = I_sorted(1:temp_topX2);
x = temp_F_dff_decisionBin(temp_I,:)';
y = offloadingProb_inOne';
mdl_all = fitglm(x,y);

r2_all = mdl_all.Rsquared.Adjusted;

offloadingProb_inOne_hat = predict(mdl_all,x);

%% Plot
close all
fig1 = figure('Name','Fig1','NumberTitle','off');
% set(gcf,'Position',[35+0 35+0 1630 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
% set(gcf,'Position',[35+0 35+0 1630 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set(gcf,'Position',[35+0 35+0 1630*0.54 500*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(2,5,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>

for tempi=1:length(highCorrIndex_rProb)
    
    if tempi > 10
        break
    end
    
    nexttile
    if tempi <= 9
        x = temp_F_dff_decisionBin(highCorrIndex_rProb(tempi),:);
    elseif tempi == 10
        %x = offloadingProb_inOne;
        x = offloadingProb_inOne_hat;
    end
    y = offloadingProb_inOne;
    scatter(x, y, 5, 'filled',...
        'MarkerFaceColor', [0 0 0], 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4);
    hold on
    
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
    set(gca, 'FontSize', 12)
    set(gca,'box','off');% 取消右、上边框
    
    if tempi <= 9
        xlabel(sprintf('dF / F, r2=%.3f',r2(highCorrIndex_rProb(tempi))), 'FontSize', 12, 'FontWeight', 'bold');
        %ylabel(sprintf('offloading rate\ncellID=%d',highCorrIndex_suite2p_rProb(tempi)), 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(sprintf('cellID=%d\noffloading rate',highCorrIndex_suite2p_rProb(tempi)), 'FontSize', 12, 'FontWeight', 'bold');        
    elseif tempi == 10
        xlabel(sprintf('regression, r2=%.3f',r2_all), 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(sprintf('offloading rate\nPopulation (top%d)',topX2), 'FontSize', 12, 'FontWeight', 'bold');        
    end
    
end

%% End

