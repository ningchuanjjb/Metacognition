close all

if_loopSave = 0;
ifSave_fig = 0;

if_regress_trial0_seq1 = 1;
if_glm_linear0_log1 = 0;


topX = 10;



if_linear0_nonlinear1 = 1;



temp_gAcc_noChoice_inOne = gAcc_noChoice_collapsed_inOne;

F_dff_decisionBin1 = mean(F_dff_decisionPeriod(:,:,decisionPeriod_interval(1):decisionPeriod_interval(2)),3);
% F_dff_decisionBin2 = mean(F_dff_decisionPeriod(:,:,17:46),3);
% F_dff_decisionBin3 = mean(F_dff_decisionPeriod(:,:,1:30),3);
% F_dff_decisionBin4 = mean(F_dff_decisionPeriod(:,:,45:75),3);
% F_dff_decisionBin5 = mean(F_dff_decisionPeriod(:,:,25:46),3);

temp_F_dff_decisionBin = F_dff_decisionBin1;

seqIndex = zeros(1,trial_para.trial_count);
target_seqSet;
for tempi=1:trial_para.trial_count
    currentSequence = trial_para.currentSequence{tempi};
    temp_seq_length = length(currentSequence);
    for tempj=1:numSeq(temp_seq_length)
        if sum(ismember(currentSequence,target_seqSet{temp_seq_length}{tempj})) == temp_seq_length
            break
        end
    end
    seqIndex(tempi) = sum(numSeq(1:temp_seq_length-1)) + tempj;
end

gAcc_inTrial = temp_gAcc_noChoice_inOne(seqIndex);
gAcc_inSeq = gAcc_noChoice_collapsed_inOne;


roiNum = size(temp_F_dff_decisionBin,1);

a = 1;

temp_F_dff_decisionBin_seqMerged = zeros(roiNum,sum(numSeq));
temp_F_dff_decisionBin;

for tempi=1:sum(numSeq)
    temp_dff = temp_F_dff_decisionBin(:,seqIndex==tempi);
    temp_F_dff_decisionBin_seqMerged(:,tempi) = mean(temp_dff,2);
end
temp_F_dff_decisionBin_raw = temp_F_dff_decisionBin;

if if_loopSave == 1
    loopNum = 4;
else
    loopNum = 1;
end
for loopCount=1:loopNum
    if if_loopSave == 1
        if loopCount == 1
            if_regress_trial0_seq1 = 0;
            if_glm_linear0_log1 = 0;
        elseif loopCount == 2
            if_regress_trial0_seq1 = 0;
            if_glm_linear0_log1 = 1;
        elseif loopCount == 3
            if_regress_trial0_seq1 = 1;
            if_glm_linear0_log1 = 0;
        elseif loopCount == 4
            if_regress_trial0_seq1 = 1;
            if_glm_linear0_log1 = 1;
        end
    end
    
    if if_regress_trial0_seq1 == 1
        temp_F_dff_decisionBin = temp_F_dff_decisionBin_seqMerged;
        temp_gAcc = gAcc_inSeq;
    elseif if_regress_trial0_seq1 == 0
        temp_gAcc = gAcc_inTrial;
    end
    
    
    x = temp_F_dff_decisionBin';
    y = temp_gAcc';
    if if_linear0_nonlinear1 == 0
        [r,p] = corr(x,y);
        [r_sorted,I_sorted] = sort(abs(r),'descend');
        p_corr_sorted = p(I_sorted);
        
        %x = temp_F_dff_decisionBin(highCorrIndex_gAcc,:)';
        %y = gAcc_inTrial';
        %[temp_r,temp_p_corr] = corr(x,y);
        
        highCorrIndex_gAcc = I_sorted(1:10);
        highCorrIndex_suite2p_gAcc = highCorrIndex_gAcc - 1;
        
    elseif if_linear0_nonlinear1 == 1
        mdl = cell(roiNum,1);
        beta = zeros(roiNum,2);
        r2 = zeros(roiNum,1);
        warning('off');
        for tempi=1:roiNum
            if if_glm_linear0_log1 == 0
                temp_mdl = fitglm(x(:,tempi),y);
            elseif if_glm_linear0_log1 == 1
                temp_mdl = fitglm(x(:,tempi),y,'Link','log');
            end
            
            %temp_mdl = fitglm(x(:,tempi),y,'Distribution','poisson');
            %temp_mdl = fitglm(x(:,tempi),y,'Distribution','Gamma');
            
            beta(tempi,:) = temp_mdl.Coefficients.Estimate; %#ok<*NASGU>
            r2(tempi) = temp_mdl.Rsquared.Adjusted;
            mdl{tempi} = temp_mdl;
            a = 1;
        end
        
        [r2_sorted,I_sorted] = sort(r2,'descend');
        
        highCorrIndex_gAcc = I_sorted(1:topX);
        highCorrIndex_suite2p_gAcc = cellIndex(highCorrIndex_gAcc) - 1;
        
        if if_regress_trial0_seq1 == 0
            x = temp_F_dff_decisionBin';
        elseif if_regress_trial0_seq1 == 1
            temp_I = I_sorted(1:50);
            x = temp_F_dff_decisionBin(temp_I,:)';
        end
        y = temp_gAcc';
        if if_glm_linear0_log1 == 0
            mdl_all = fitglm(x,y);
        elseif if_glm_linear0_log1 == 1
            mdl_all = fitglm(x,y,'Link','log');
        end
        %mdl_all = fitglm(x,y,'Distribution','poisson');
        %mdl_all = fitglm(x,y,'Distribution','Gamma');
        
        r2_all = mdl_all.Rsquared.Adjusted;
        
        x = temp_F_dff_decisionBin(highCorrIndex_gAcc,:)';
        y = temp_gAcc';
        if if_glm_linear0_log1 == 0
            mdl_topX = fitglm(x,y);
        elseif if_glm_linear0_log1 == 1
            mdl_topX = fitglm(x,y,'Link','log');
        end
        
        
        %mdl_topX = fitglm(x,y,'Distribution','poisson');
        %mdl_topX = fitglm(x,y,'Distribution','Gamma');
        
        r2_topX = mdl_topX.Rsquared.Adjusted;
        
        a = 1;
        warning('on');
        
        
        
        %
        %     mdl = cell(roiNum,1);
        %     beta = zeros(roiNum,2);
        %     r2 = zeros(roiNum,1);
        %     warning('off');
        %     for tempi=1:roiNum
        % %         temp_mdl = fitglm(x(:,tempi),y);
        % %         temp_mdl = fitglm(x(:,tempi),y,'Distribution','poisson');
        % %         temp_mdl = fitglm(x(:,tempi),y,'Distribution','Gamma');
        %         temp_mdl = fitglm(x(:,tempi),y,'Link','log');
        %         beta(tempi,:) = temp_mdl.Coefficients.Estimate; %#ok<*NASGU>
        %         r2(tempi) = temp_mdl.Rsquared.Adjusted;
        %         mdl{tempi} = temp_mdl;
        %         a = 1;
        %     end
        %
        %     [r2_sorted,I_sorted] = sort(r2,'descend');
        %
        %     highCorrIndex_gAcc = I_sorted(1:10);
        %     highCorrIndex_suite2p_gAcc = highCorrIndex_gAcc - 1;
        %
        %     x = temp_F_dff_decisionBin';
        %     y = gAcc_inTrial';
        % %     mdl_all = fitglm(x,y);
        % %     mdl_all = fitglm(x,y,'Distribution','poisson');
        % %     mdl_all = fitglm(x,y,'Distribution','Gamma');
        %     mdl_all = fitglm(x,y,'Link','log');
        %     r2_all = mdl_all.Rsquared.Adjusted;
        %
        %     x = temp_F_dff_decisionBin(highCorrIndex_gAcc,:)';
        %     y = gAcc_inTrial';
        % %     mdl_topX = fitglm(x,y);
        % %     mdl_topX = fitglm(x,y,'Distribution','poisson');
        % %     mdl_topX = fitglm(x,y,'Distribution','Gamma');
        %     mdl_topX = fitglm(x,y,'Link','log');
        %     r2_topX = mdl_topX.Rsquared.Adjusted;
        %
        %     a = 1;
        %     warning('on');
        
        
        
        
        %     n = 2;
        %     r2 = zeros(roiNum,1);
        %     coeff = zeros(roiNum,n+1);
        %     for tempi=1:roiNum
        %         [coeff(tempi,:),S] = polyfit(x(:,tempi),y,n);
        %         r2(tempi) = 1 - (S.normr/norm(y - mean(y)))^2;
        %         a = 1;
        %     end
        %     [r2_sorted,I_sorted] = sort(r2,'descend');
        %     highCorrIndex_gAcc = I_sorted(1:10);
        %     highCorrIndex_suite2p_gAcc = highCorrIndex_gAcc - 1;
        
        
    end
    
    
    
    %% Plot
    fig1 = figure('Name','Fig1','NumberTitle','off');
    set(gcf,'Position',[35+0 35+0 1630 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,5,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
    if if_linear0_nonlinear1 == 1
        if if_regress_trial0_seq1 == 0
            title(t,sprintf('r2(all) = %.3f, r2(top%d) = %.3f',r2_all,topX,r2_topX),'FontSize',20,'FontWeight','bold');
        elseif if_regress_trial0_seq1 == 1
            title(t,sprintf('r2(top50) = %.3f, r2(top%d) = %.3f',r2_all,topX,r2_topX),'FontSize',20,'FontWeight','bold');
        end
    end
    
    for tempi=1:length(highCorrIndex_gAcc)
        nexttile
        x = temp_F_dff_decisionBin(highCorrIndex_gAcc(tempi),:);
        y = temp_gAcc;
        scatter(x, y, 5, 'filled',...
            'MarkerFaceColor', [0 0 0], 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4);
        hold on
        
        [temp_min,temp_max] = bounds(x);
        
        if if_linear0_nonlinear1 == 0
            n = 1;
            [p_mapping,S] = polyfit(x,y,n);
            x_fit = temp_min:0.1:temp_max;
            y_fit = polyval(p_mapping,x_fit);
            
        elseif if_linear0_nonlinear1 == 1
            %n = 2;
            %[p_mapping,S] = polyfit(x,y,n);
            %x_fit = temp_min:0.1:temp_max;
            %y_fit = polyval(p_mapping,x_fit);
            
            x_fit = temp_min:0.01:temp_max;
            y_fit = predict(mdl{highCorrIndex_gAcc(tempi)},x_fit')';
            a = 1;
        end
        plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);
        hold on
        
        xlim([temp_min temp_max]);
        ylim([0 1]);
        set(gca, 'FontSize', 11)
        set(gca,'box','off');% 取消右、上边框
        if if_linear0_nonlinear1 == 0
            xlabel(sprintf('dF / F, r=%.3f',r(highCorrIndex_gAcc(tempi))), 'FontSize', 14, 'FontWeight', 'bold');
        elseif if_linear0_nonlinear1 == 1
            xlabel(sprintf('dF / F, r2=%.3f',r2(highCorrIndex_gAcc(tempi))), 'FontSize', 14, 'FontWeight', 'bold');
        end
        ylabel(sprintf('Accracy\ncellID=%d',highCorrIndex_suite2p_gAcc(tempi)), 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    
    % Save fig
    if if_regress_trial0_seq1 == 0
        if if_glm_linear0_log1 == 0
            currentFigName = sprintf('regress_trial_gAcc_linear');
        elseif if_glm_linear0_log1 == 1
            currentFigName = sprintf('regress_trial_gAcc_log');
        end
        
    elseif if_regress_trial0_seq1 == 1
        if if_glm_linear0_log1 == 0
            currentFigName = sprintf('regress_seq_gAcc_linear');
        elseif if_glm_linear0_log1 == 1
            currentFigName = sprintf('regress_seq_gAcc_log');
        end
    end
    % currentFigName = sprintf('regress_gAcc');
    fileName_fig = [output_path_singleCell '\' currentFigName '.png'];%emf,pdf,png,jpg,tiff
    % to judge whether save figure or not
    if ifSave_fig == 1
        %exportgraphics(fig1,fileName_fig,'BackgroundColor','none','ContentType','vector');
        exportgraphics(fig1,fileName_fig,'Resolution',300);
        close all
    end
    
end

%% End

