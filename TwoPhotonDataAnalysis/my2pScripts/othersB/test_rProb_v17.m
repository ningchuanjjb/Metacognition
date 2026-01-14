close all


if_loopSave = 0;
ifSave_fig = 0;

if_regress_trial0_seq1 = 1;
if_glm_linear0_log1 = 0;


if_shuffle = 0;
shuffleNum = 1000;

threshold_r2 = 0.01;%0.01-->0.15

topX = 10;%10
topX2 = 50;%50

temp_lengthFlag_all0_1_2_3_4 = 0;

if_linear0_nonlinear1 = 1;


temp_offloadingProb_inOne = offloadingProb_inOne;

temp_F_dff_decisionBin = F_dff_delay1Bin;
% temp_F_dff_decisionBin = F_dff_delay2Bin;
% temp_F_dff_decisionBin = F_dff_baselineBin;


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
offloadingProb_inTrial = temp_offloadingProb_inOne(seqIndex);
offloadingProb_inSeq = offloadingProb_inOne;
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
        temp_offloadingProb = offloadingProb_inSeq;
    elseif if_regress_trial0_seq1 == 0
        temp_offloadingProb = offloadingProb_inTrial;
    end
    
    if if_regress_trial0_seq1 == 1
        tempi = temp_lengthFlag_all0_1_2_3_4;
        if tempi == 0
            
        else
            temp_F_dff_decisionBin_raw2 = temp_F_dff_decisionBin;
            temp_offloadingProb_raw2 = temp_offloadingProb;
            if tempi == 1
                temp_range = 1:numSeq(1);
            elseif tempi > 1
                temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
            end
            temp_F_dff_decisionBin = temp_F_dff_decisionBin_raw2(:,temp_range);
            temp_offloadingProb = temp_offloadingProb_raw2(temp_range);
            
            topX = min(topX,numSeq(tempi)-2);
            
        end
    end
    
    x = temp_F_dff_decisionBin';
    y = temp_offloadingProb';
    if if_linear0_nonlinear1 == 0
        [r,p] = corr(x,y);
        [r_sorted,I_sorted] = sort(abs(r),'descend');
        p_corr_sorted = p(I_sorted);
                
        highCorrIndex_rProb = I_sorted(1:topX);
        highCorrIndex_suite2p_rProb = highCorrIndex_rProb - 1;
        
    elseif if_linear0_nonlinear1 == 1
        mdl = cell(roiNum,1);
        beta = zeros(roiNum,2);
        r2 = zeros(roiNum,1);
        p = zeros(roiNum,1);
        warning('off');
        for tempi=1:roiNum
            if if_glm_linear0_log1 == 0
                temp_mdl = fitglm(x(:,tempi),y);
            elseif if_glm_linear0_log1 == 1
                temp_mdl = fitglm(x(:,tempi),y,'Link','log');
            end
                        
            beta(tempi,:) = temp_mdl.Coefficients.Estimate; %#ok<*NASGU>
            r2(tempi) = temp_mdl.Rsquared.Adjusted;
            p(tempi) = temp_mdl.Coefficients.pValue(2);
            mdl{tempi} = temp_mdl;
            a = 1;
        end
        
        [r2_sorted,I_sorted] = sort(r2,'descend');
        
        highCorrIndex_rProb = I_sorted(1:topX);
        highCorrIndex_suite2p_rProb = cellIndex(highCorrIndex_rProb) - 1;
        
        if if_regress_trial0_seq1 == 0
            x = temp_F_dff_decisionBin';
        elseif if_regress_trial0_seq1 == 1
            temp_topX2 = min(topX2,sum(r2_sorted>threshold_r2));
            temp_I = I_sorted(1:temp_topX2);
            x = temp_F_dff_decisionBin(temp_I,:)';
        end
        y = temp_offloadingProb';
        if if_glm_linear0_log1 == 0
            mdl_all = fitglm(x,y);
        elseif if_glm_linear0_log1 == 1
            mdl_all = fitglm(x,y,'Link','log');
        end
        
        r2_all = mdl_all.Rsquared.Adjusted;
        
        a = 1;
        
        if if_shuffle == 1 && if_regress_trial0_seq1 == 1
            %shuffleNum = 10000;
            r2_all_shuffle = zeros(1,shuffleNum);
            
            t0 = tic;
            
%             for tempi=1:shuffleNum
            parfor tempi=1:shuffleNum

                warning('off');
                x = temp_F_dff_decisionBin';
                temp_offloadingProb_shuffle = temp_offloadingProb(randperm(length(temp_offloadingProb))); %#ok<*PFBNS>
                y = temp_offloadingProb_shuffle';                
    
                mdl_shuffle = cell(roiNum,1);
                beta_shuffle = zeros(roiNum,2);
                r2_shuffle = zeros(roiNum,1);
                for tempj=1:roiNum
                %parfor tempj=1:roiNum    
                    temp_mdl = []; %#ok<*PFTUSW>
                    if if_glm_linear0_log1 == 0
                        temp_mdl = fitglm(x(:,tempj),y);
                    elseif if_glm_linear0_log1 == 1
                        temp_mdl = fitglm(x(:,tempj),y,'Link','log');
                    end
                    beta_shuffle(tempj,:) = temp_mdl.Coefficients.Estimate; %#ok<*NASGU>
                    r2_shuffle(tempj) = temp_mdl.Rsquared.Adjusted;
                    mdl_shuffle{tempj} = temp_mdl;
                    a = 1;
                end
                [r2_sorted_shuffle,I_sorted_shuffle] = sort(r2_shuffle,'descend');
                
                temp_topX2_shuffle = min(topX2,sum(r2_sorted_shuffle>threshold_r2)); 
                temp_I = I_sorted_shuffle(1:temp_topX2_shuffle);
                x = temp_F_dff_decisionBin(temp_I,:)';
                
                if if_glm_linear0_log1 == 0
                    mdl_all_shuffle = fitglm(x,y);
                elseif if_glm_linear0_log1 == 1
                    mdl_all_shuffle = fitglm(x,y,'Link','log');
                end
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
        
        x = temp_F_dff_decisionBin(highCorrIndex_rProb,:)';
        y = temp_offloadingProb';
        if if_glm_linear0_log1 == 0
            mdl_topX = fitglm(x,y);
        elseif if_glm_linear0_log1 == 1
            mdl_topX = fitglm(x,y,'Link','log');
        end
        
        r2_topX = mdl_topX.Rsquared.Adjusted;
        
        warning('on');
    end
    
    
    
    %% Plot
    fig1 = figure('Name','Fig1','NumberTitle','off');
    set(gcf,'Position',[35+0 35+0 1630 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,5,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
    if if_linear0_nonlinear1 == 1
        if if_regress_trial0_seq1 == 0
            title(t,sprintf('r2(all) = %.3f, r2(top%d) = %.3f',r2_all,topX,r2_topX),'FontSize',20,'FontWeight','bold');
        elseif if_regress_trial0_seq1 == 1
            title(t,sprintf('r2(top%d) = %.3f, r2(top%d) = %.3f',temp_topX2,r2_all,topX,r2_topX),'FontSize',20,'FontWeight','bold');
        end
    end
    
    for tempi=1:length(highCorrIndex_rProb)
        
        if tempi > 10
            break
        end
        
        nexttile
        x = temp_F_dff_decisionBin(highCorrIndex_rProb(tempi),:);
        y = temp_offloadingProb;
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
            
            x_fit = temp_min:0.001:temp_max;
            y_fit = predict(mdl{highCorrIndex_rProb(tempi)},x_fit')';
            a = 1;
        end
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
    
    % Save fig
    if if_regress_trial0_seq1 == 0
        if if_glm_linear0_log1 == 0
            currentFigName = sprintf('regress_trial_rProb_linear');
        elseif if_glm_linear0_log1 == 1
            currentFigName = sprintf('regress_trial_rProb_log');
        end
        
    elseif if_regress_trial0_seq1 == 1
        if if_glm_linear0_log1 == 0
            currentFigName = sprintf('regress_seq_rProb_linear');
        elseif if_glm_linear0_log1 == 1
            currentFigName = sprintf('regress_seq_rProb_log');
        end
    end
    % currentFigName = sprintf('regress_rProb');
    fileName_fig = [output_path_singleCell '\' currentFigName '.png'];%emf,pdf,png,jpg,tiff
    % to judge whether save figure or not
    if ifSave_fig == 1
        %exportgraphics(fig1,fileName_fig,'BackgroundColor','none','ContentType','vector');
        exportgraphics(fig1,fileName_fig,'Resolution',300);
        close all
    end
    
end


%% End

