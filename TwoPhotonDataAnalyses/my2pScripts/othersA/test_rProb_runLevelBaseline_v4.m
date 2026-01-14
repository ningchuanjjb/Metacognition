close all

runSize = 207;%207-->104-->52-->42

threshold_r2 = 0.01;%0.01-->0.15
topX = 1;%10
topX2 = 1;%50

if_shuffle = 0;
shuffleNum = 1000;


runNum = ceil(trial_para.trial_count/runSize);

F_dff_baselineBin = mean(F_dff_baselinePeriod(:,:,baselinePeriod_interval(1):baselinePeriod_interval(end)),3);
temp_F_dff_baselineBin = F_dff_baselineBin;

offloadingProb_bin_runLevel = zeros(1,runNum);
temp_F_dff_baselineBin_runLevel = zeros(size(F_dff_baselinePeriod,1),runNum);
for tempi=1:runNum
    temp_range = (1:runSize) + (tempi-1)*runSize;
    temp_range(temp_range>trial_para.trial_count) = [];
    
    temp_choiceOffloadNum = sum(trialIndex_bool_choiceOffload(temp_range));
    temp_choiceMemoryNum = sum(trialIndex_bool_choiceMemory(temp_range));
    
    offloadingProb_bin_runLevel(tempi) = temp_choiceOffloadNum/(temp_choiceMemoryNum+temp_choiceOffloadNum);
    temp_F_dff_baselineBin_runLevel(:,tempi) = mean(temp_F_dff_baselineBin(:,temp_range),2);
end
offloadingProb_bin_runLevel;
temp_F_dff_baselineBin_runLevel;

temp_offloadingProb = offloadingProb_bin_runLevel;
temp_F_dff = temp_F_dff_baselineBin_runLevel;

roiNum = size(temp_F_dff,1);



x = temp_F_dff';
y = temp_offloadingProb';


mdl = cell(roiNum,1);
beta = zeros(roiNum,2);
r2 = zeros(roiNum,1);
warning('off');
for tempi=1:roiNum
    temp_mdl = fitglm(x(:,tempi),y);
    beta(tempi,:) = temp_mdl.Coefficients.Estimate; %#ok<*NASGU>
    r2(tempi) = temp_mdl.Rsquared.Adjusted;
    mdl{tempi} = temp_mdl;
    a = 1;
end

[r2_sorted,I_sorted] = sort(r2,'descend');

highCorrIndex_rProb = I_sorted(1:topX);
highCorrIndex_suite2p_rProb = cellIndex(highCorrIndex_rProb) - 1;


temp_topX2 = min(topX2,sum(r2_sorted>threshold_r2));
temp_I = I_sorted(1:temp_topX2);
x = temp_F_dff(temp_I,:)';

y = temp_offloadingProb';
mdl_all = fitglm(x,y);
r2_all = mdl_all.Rsquared.Adjusted;

a = 1;
if if_shuffle == 1
    r2_all_shuffle = zeros(1,shuffleNum);
    
    t0 = tic;
    
    parfor tempi=1:shuffleNum
        
        warning('off');
        x = temp_F_dff';
        temp_offloadingProb_shuffle = temp_offloadingProb(randperm(length(temp_offloadingProb))); %#ok<*PFBNS>
        y = temp_offloadingProb_shuffle';
        
        mdl_shuffle = cell(roiNum,1);
        beta_shuffle = zeros(roiNum,2);
        r2_shuffle = zeros(roiNum,1);
        for tempj=1:roiNum
            temp_mdl = []; %#ok<*PFTUSW>
            temp_mdl = fitglm(x(:,tempj),y);
            beta_shuffle(tempj,:) = temp_mdl.Coefficients.Estimate; %#ok<*NASGU>
            r2_shuffle(tempj) = temp_mdl.Rsquared.Adjusted;
            mdl_shuffle{tempj} = temp_mdl;
            a = 1;
        end
        [r2_sorted_shuffle,I_sorted_shuffle] = sort(r2_shuffle,'descend');
        
        temp_topX2_shuffle = min(topX2,sum(r2_sorted_shuffle>threshold_r2));
        temp_I = I_sorted_shuffle(1:temp_topX2_shuffle);
        x = temp_F_dff(temp_I,:)';
        
        mdl_all_shuffle = fitglm(x,y);
        r2_all_shuffle(tempi) = mdl_all_shuffle.Rsquared.Adjusted;
        a = 1;
    end
    r2_all_shuffle_median = median(r2_all_shuffle);
    h = histogram(r2_all_shuffle);
    hold on
    plot(r2_all*[1 1],[min(h.Values) max(h.Values)],'Color',[0.75 0.25 0.25],'LineWidth',2);
    hold on
    tempPercentile = sum(r2_all_shuffle>r2_all)/shuffleNum;
    
    text((max(r2_all_shuffle)-min(r2_all_shuffle))*0.05+min(r2_all_shuffle),(max(h.Values)-min(h.Values))*0.9+min(h.Values),...
        sprintf('percentile = %.4f',tempPercentile),'FontSize',16);
    
    fprintf('Shuffle done, time = %.0f secs, tempPercentile = %.4f.\n',toc(t0),tempPercentile);
    a = 1;
end

x = temp_F_dff(highCorrIndex_rProb,:)';
y = temp_offloadingProb';
mdl_topX = fitglm(x,y);
r2_topX = mdl_topX.Rsquared.Adjusted;

warning('on');
a = 1;


%% Plot
fig1 = figure('Name','Fig1','NumberTitle','off');
set(gcf,'Position',[35+0 35+0 1630 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(2,5,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
title(t,sprintf('r2(top%d) = %.3f, r2(top%d) = %.3f',temp_topX2,r2_all,topX,r2_topX),'FontSize',20,'FontWeight','bold');

for tempi=1:length(highCorrIndex_rProb)
    
    if tempi > 10
        break
    end
    
    nexttile
    x = temp_F_dff(highCorrIndex_rProb(tempi),:);
    y = temp_offloadingProb;
    scatter(x, y, 5, 'filled',...
        'MarkerFaceColor', [0 0 0], 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4);
    hold on
    
    [temp_min,temp_max] = bounds(x);
    

        x_fit = temp_min:0.001:temp_max;
        y_fit = predict(mdl{highCorrIndex_rProb(tempi)},x_fit')';
        a = 1;

    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    
    xlim([temp_min temp_max]);
    ylim([0 1]);
    set(gca, 'FontSize', 11)
    set(gca,'box','off');% 取消右、上边框
    if if_linear0_nonlinear1 == 0
        xlabel(sprintf('dF / F, r=%.3f',r(highCorrIndex_rProb(tempi))), 'FontSize', 14, 'FontWeight', 'bold');
    elseif if_linear0_nonlinear1 == 1
        xlabel(sprintf('dF / F, r2=%.3f',r2(highCorrIndex_rProb(tempi))), 'FontSize', 14, 'FontWeight', 'bold');
    end
    ylabel(sprintf('offloading rate\ncellID=%d',highCorrIndex_suite2p_rProb(tempi)), 'FontSize', 14, 'FontWeight', 'bold');
end

a = 1;
