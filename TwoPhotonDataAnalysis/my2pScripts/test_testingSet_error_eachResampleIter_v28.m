% Chuan's 21th script (20251214)
% This script: To test location memory decoder in error trial.
if_plot = 1;

if_correlateWith_theoretical0_behavior1 = 0;


%% Get signiSeqNum_resample_stimuliError and signiSeqNum_resample_responseError
eachIter_options = struct;
eachIter_options.svm_train_length1_outputs = svm_train_length1_outputs;
eachIter_options.svm_train_length2_outputs = svm_train_length2_outputs;
eachIter_options.svm_train_length3_outputs = svm_train_length3_outputs;

eachIter_options.trialIndex_bool_memoryError = trialIndex_bool_memoryError;
eachIter_options.numFrames = numFrames;
eachIter_options.boolIndex_location_seq = boolIndex_location_seq;
eachIter_options.numSeq = numSeq;
eachIter_options.valid_length = valid_length;
if if_correlateWith_theoretical0_behavior1 == 0
    C = double(boolIndex_location_seq_T(1:sum(numSeq(1:3)),:));
elseif if_correlateWith_theoretical0_behavior1 == 1
    C = Posterior_2d_model;
end
eachIter_options.Posterior_2d_model = C;
eachIter_options.errorTrial_minNum = errorTrial_minNum;
eachIter_options.if_stimuliError0_responseError1 = if_stimuliError0_responseError1;
eachIter_options.if_n11n = if_n11n;
% eachIter_options.if_n11n = 0;

fun_testError_eachResampleIter_locDistri_Name_v = autoGetFunName_myScripts('fun_testError_eachResampleIter_locDistri', [targetPATH '\functions']);
fun_testError_eachResampleIter_locDistri = str2func(fun_testError_eachResampleIter_locDistri_Name_v);

r_n11n_resample_stimuliError = zeros(resampleIterCount,sum(numSeq(1:valid_length)));
p_n11n_resample_stimuliError = zeros(resampleIterCount,sum(numSeq(1:valid_length)));
Posterior_2d_n11n_resample_stimuliError = zeros(resampleIterCount,sum(numSeq(1:valid_length)),numFrames);
eachIter_options.if_stimuliError0_responseError1 = 0;
for iterCount=1:resampleIterCount
    [r_n11n_resample_stimuliError(iterCount,:),p_n11n_resample_stimuliError(iterCount,:),Posterior_2d_n11n_resample_stimuliError(iterCount,:,:)] = ...
        fun_testError_eachResampleIter_locDistri(iterCount,eachIter_options);
end
r_n11n_resampleMean_stimuliError = mean(r_n11n_resample_stimuliError,1);
p_n11n_resampleMean_stimuliError = mean(p_n11n_resample_stimuliError,1);
Posterior_2d_n11n_resampleMean_stimuliError = squeeze(mean(Posterior_2d_n11n_resample_stimuliError,1));
temp1 = p_n11n_resample_stimuliError<0.05;
signiSeqNum_resample_stimuliError = sum(temp1,2);
signiSeqProportion_resample_stimuliError = signiSeqNum_resample_stimuliError./sum(~isnan(p_n11n_resampleMean_stimuliError));
% fprintf('num(p_posterior_seq_n11n<0.05, stimuliError)=%.1f/%d.\n',mean(signiSeqNum_resample_stimuliError),sum(~isnan(p_n11n_resampleMean_stimuliError)));

r_n11n_resample_responseError = zeros(resampleIterCount,sum(numSeq(1:valid_length)));
p_n11n_resample_responseError = zeros(resampleIterCount,sum(numSeq(1:valid_length)));
Posterior_2d_n11n_resample_responseError = zeros(resampleIterCount,sum(numSeq(1:valid_length)),numFrames);
eachIter_options.if_stimuliError0_responseError1 = 1;
for iterCount=1:resampleIterCount
    [r_n11n_resample_responseError(iterCount,:),p_n11n_resample_responseError(iterCount,:),Posterior_2d_n11n_resample_responseError(iterCount,:,:)] = ...
        fun_testError_eachResampleIter_locDistri(iterCount,eachIter_options);
end
r_n11n_resampleMean_responseError = mean(r_n11n_resample_responseError,1);
p_n11n_resampleMean_responseError = mean(p_n11n_resample_responseError,1);
Posterior_2d_n11n_resampleMean_responseError = squeeze(mean(Posterior_2d_n11n_resample_responseError,1));
temp1 = p_n11n_resample_responseError<0.05;
signiSeqNum_resample_responseError = sum(temp1,2);
signiSeqProportion_resample_responseError = signiSeqNum_resample_responseError./sum(~isnan(p_n11n_resampleMean_responseError));
% fprintf('num(p_posterior_seq_n11n<0.05, responseError)=%.1f/%d.\n',mean(signiSeqNum_resample_responseError),sum(~isnan(p_n11n_resampleMean_responseError)));

%% Get signiSeqNum_resample_correctFew
if exist('fewTrialCountSeqBoolIndex','var') == 1
    
    fun_testCorrectFew_eachResampleIter_locDistri_Name_v = autoGetFunName_myScripts('fun_testCorrectFew_eachResampleIter_locDistri', [targetPATH '\functions']);
    fun_testCorrectFew_eachResampleIter_locDistri = str2func(fun_testCorrectFew_eachResampleIter_locDistri_Name_v);
    
    eachIter_options.svm_train_length1_outputs = svm_train_length1_outputs;
    eachIter_options.svm_train_length2_outputs = svm_train_length2_outputs;
    eachIter_options.svm_train_length3_outputs = svm_train_length3_outputs;
    
    eachIter_options.trialIndex_bool_memoryCorrectFew = trialIndex_bool_memoryCorrectFew;
    eachIter_options.numFrames = numFrames;
    eachIter_options.boolIndex_location_seq = boolIndex_location_seq;
    eachIter_options.numSeq = numSeq;
    eachIter_options.valid_length = valid_length;
    eachIter_options.Posterior_2d_model = Posterior_2d_model;
    eachIter_options.minTrialCount = minTrialCount;
    eachIter_options.if_n11n = if_n11n;
    
    r_n11n_resample_correctFew = zeros(resampleIterCount,sum(numSeq(1:valid_length)));
    p_n11n_resample_correctFew = zeros(resampleIterCount,sum(numSeq(1:valid_length)));
    eachIter_options.if_stimuliError0_responseError1 = 0;
    for iterCount=1:resampleIterCount
        [r_n11n_resample_correctFew(iterCount,:),p_n11n_resample_correctFew(iterCount,:)] = ...
            fun_testCorrectFew_eachResampleIter_locDistri(iterCount,eachIter_options);
    end
    r_n11n_resampleMean_correctFew = mean(r_n11n_resample_correctFew,1);
    p_n11n_resampleMean_correctFew = mean(p_n11n_resample_correctFew,1);
end

%% Get signiSeqNum_resample_correct
Posterior_2d_n11n_resample = zeros(resampleIterCount,size(Posterior_2d_n11n_mean,1),size(Posterior_2d_n11n_mean,2));
temp_r_n11n_resample_correct = zeros(resampleIterCount,sum(numSeq(1:valid_length)));
temp_p_n11n_resample_correct = zeros(resampleIterCount,sum(numSeq(1:valid_length)));
temp_r_n11n_resample_correct_chanceLevel = zeros(resampleIterCount,sum(numSeq(1:valid_length)));
for tempi=1:resampleIterCount
    temp_posterior = [];
    temp_posterior = [temp_posterior; svm_train_length1_outputs.temp_svm_resample{tempi}.Posterior_2d_n11n_lengthx_mean]; %#ok<*AGROW>
    temp_posterior = [temp_posterior; svm_train_length2_outputs.temp_svm_resample{tempi}.Posterior_2d_n11n_lengthx_mean];
    temp_posterior = [temp_posterior; svm_train_length3_outputs.temp_svm_resample{tempi}.Posterior_2d_n11n_lengthx_mean];
    % sum(~isnan(sum(temp_posterior,2)))
    
    Posterior_2d_n11n_resample(tempi,:,:) = temp_posterior;
    
    temp_r_n11n = zeros(1,sum(numSeq(1:valid_length)));
    temp_p_n11n = zeros(1,sum(numSeq(1:valid_length)));
    temp_r_n11n_chanceLevel = zeros(1,sum(numSeq(1:valid_length)));
    
    if if_correlateWith_theoretical0_behavior1 == 0
        C = double(boolIndex_location_seq_T(1:sum(numSeq(1:3)),:));
    elseif if_correlateWith_theoretical0_behavior1 == 1
        C = Posterior_2d_model;
    end
    
    for tempj=1:sum(numSeq(1:valid_length))
        [temp_r_n11n(tempj),temp_p_n11n(tempj)] = corr(temp_posterior(tempj,:)',C(tempj,:)');
        %[temp_r_n11n_chanceLevel(tempj),~] = corr(temp_posterior(tempj,:)',rand(numFrames,1));
        [temp_r_n11n_chanceLevel(tempj),~] = corr(rand(numFrames,1),C(tempj,:)');
    end
    temp_r_n11n_resample_correct(tempi,:) = temp_r_n11n;
    temp_p_n11n_resample_correct(tempi,:) = temp_p_n11n;
    temp_r_n11n_resample_correct_chanceLevel(tempi,:) = temp_r_n11n_chanceLevel;
end
temp_r_n11n_resampleMean_correct = mean(temp_r_n11n_resample_correct,1);
temp_p_n11n_resampleMean_correct = mean(temp_p_n11n_resample_correct,1);
temp_r_n11n_resampleMean_correct_chanceLevel = mean(temp_r_n11n_resample_correct_chanceLevel,1);

if exist('fewTrialCountSeqBoolIndex','var') == 1
    r_n11n_resampleMean_correct_raw = temp_r_n11n_resampleMean_correct;
    temp_r_n11n_resampleMean_correct(fewTrialCountSeqBoolIndex) = r_n11n_resampleMean_correctFew(fewTrialCountSeqBoolIndex);
    p_n11n_resampleMean_correct_raw = temp_p_n11n_resampleMean_correct;
    temp_p_n11n_resampleMean_correct(fewTrialCountSeqBoolIndex) = p_n11n_resampleMean_correctFew(fewTrialCountSeqBoolIndex);
end

temp1 = temp_p_n11n_resample_correct<0.05;
signiSeqNum_resample_correct = sum(temp1,2);
signiSeqProportion_resample_correct = signiSeqNum_resample_correct./sum(~isnan(temp_p_n11n_resampleMean_correct));
% fprintf('num(p_posterior_seq_n11n<0.05, correct)=%.1f/%d.\n',mean(signiSeqNum_resample_correct),sum(~isnan(temp_p_n11n_resampleMean_correct)));

a = 1;

signiSeqNum_resample_correct;
signiSeqNum_resample_stimuliError;
signiSeqNum_resample_responseError;

signiSeqProportion_resample_correct;
signiSeqProportion_resample_stimuliError;
signiSeqProportion_resample_responseError;

% [~,p12,~,~] = ttest(signiSeqProportion_resample_correct,signiSeqProportion_resample_stimuliError);
% [~,p13,~,~] = ttest(signiSeqProportion_resample_correct,signiSeqProportion_resample_responseError);
% [~,p23,~,~] = ttest(signiSeqProportion_resample_stimuliError,signiSeqProportion_resample_responseError);


[~,p12,~,~] = ttest(temp_r_n11n_resampleMean_correct,r_n11n_resampleMean_stimuliError);
[~,p13,~,~] = ttest(temp_r_n11n_resampleMean_correct,r_n11n_resampleMean_responseError);
[~,p23,~,~] = ttest(r_n11n_resampleMean_stimuliError,r_n11n_resampleMean_responseError);
% [~,p12,~,~] = ttest2(temp_r_n11n_resampleMean_correct,r_n11n_resampleMean_stimuliError);
% [~,p13,~,~] = ttest2(temp_r_n11n_resampleMean_correct,r_n11n_resampleMean_responseError);
% [~,p23,~,~] = ttest2(r_n11n_resampleMean_stimuliError,r_n11n_resampleMean_responseError);

% temp_1 = temp_r_n11n_resampleMean_correct(~isnan(temp_r_n11n_resampleMean_correct));
% temp_2 = r_n11n_resampleMean_stimuliError(~isnan(r_n11n_resampleMean_stimuliError));
% temp_3 = r_n11n_resampleMean_responseError(~isnan(r_n11n_resampleMean_responseError));


a = 1;

% temp_p_n11n = zeros(1,sum(numSeq(1:valid_length)));
% temp_posterior = zeros(sum(numSeq(1:valid_length)),numFrames);
% for tempj=1:sum(numSeq(1:valid_length))
%     [~,temp_p_n11n(tempj)] = corr(temp_posterior(tempj,:)',Posterior_2d_model(tempj,:)');
% end

%% r_resample_x
temp_r_n11n_resampleMean_correct;
r_n11n_resampleMean_stimuliError;
r_n11n_resampleMean_responseError;

a1 = mean(temp_r_n11n_resampleMean_correct,'omitnan');
a2 = mean(r_n11n_resampleMean_stimuliError,'omitnan');
a3 = mean(r_n11n_resampleMean_responseError,'omitnan');

b1 = mean(temp_r_n11n_resampleMean_correct_chanceLevel,'omitnan');

%% Plot
if if_plot == 1
    close all
    
    if true
        fig = figure('Name','locDistri','NumberTitle','off');
        %set(gcf,'Position',[50+0 50+0 240 252*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 240*1.1 336*0.96]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 240*0.80 336*1.11]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 240*0.80*0.78 336*1.11*0.78]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[50+0 50+0 240*0.80*0.78 336*1.11*0.78*0.9*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_p12 = p12;
        temp_p13 = p13;
        temp_p23 = p23;
        
        temp_1 = temp_r_n11n_resampleMean_correct(~isnan(temp_r_n11n_resampleMean_correct));
        temp_2 = r_n11n_resampleMean_stimuliError(~isnan(r_n11n_resampleMean_stimuliError));
        temp_3 = r_n11n_resampleMean_responseError(~isnan(r_n11n_resampleMean_responseError));
        temp_1_chanceLevel = temp_r_n11n_resampleMean_correct_chanceLevel(~isnan(temp_r_n11n_resampleMean_correct_chanceLevel));
        
        
        temp_y_min = mean(temp_1_chanceLevel);
        temp_y_max = max([temp_1,temp_2,temp_3]);
        
        temp_data = [temp_1';temp_2';temp_3'];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'C'},length(temp_3),1);
        
        temp_label = [g1;g2;g3];
        
        temptemp_color1 = [1 1 1]*0.5;%0.5
        temptemp_color2 = repmat(temptemp_color1, 3, 1);
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        h(3).ViolinPlot.FaceAlpha = 0.1;
        
        plot([0 4],mean(temp_1_chanceLevel)*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        tempTxt = sprintf('');
        if temp_p12 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p12 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p12 < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.13,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1.1 1.9],temp_y_max+(temp_y_max-temp_y_min)*0.11*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        %     tempTxt = sprintf('');
        %     if temp_p13 < 0.001
        %         tempTxt = sprintf('***');
        %     elseif temp_p13 < 0.01
        %         tempTxt = sprintf('**');
        %     elseif temp_p13 < 0.05
        %         tempTxt = sprintf('*');
        %     end
        %     text(2,temp_y_max+(temp_y_max-temp_y_min)*0.25,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        %        'HorizontalAlignment','center');
        %     plot([1.1 2.9],temp_y_max+(temp_y_max-temp_y_min)*0.23*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        %     hold on
        
        tempTxt = sprintf('');
        if temp_p23 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p23 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p23 < 0.05
            tempTxt = sprintf('*');
        end
        text(2.5,temp_y_max+(temp_y_max-temp_y_min)*0.13,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([2.1 2.9],temp_y_max+(temp_y_max-temp_y_min)*0.11*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        
        set(gca,'linewidth',1.5)
        
        %xlim([0 4]);
        xlim([0.4 4]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.3]);
        set(gca, 'FontSize', 8) %14
        %set(gca,'XTickLabel', ["correct"; "stimuliError"; "responseError"],'FontSize', 12);%给坐标加标签
        %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
        set(gca,'box','off');% 取消右、上边框
        
        
        %xtl = ["correct", "stimuliError", "responseError"];
        xtl = ["Correct", "Stimuli-error", "Response-error"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-0.20;%-0.26
        %xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.40;
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12-->11
        
        set(gca,'xticklabel','');
        
        yticks([0 1]);
        
        %ylabel('Proportion', 'FontSize', 15, 'FontWeight', 'bold');
        ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf('Location distribution correlated sequence (with model)\nerrorTrial_minNum=%d',errorTrial_minNum),'fontsize',16);
        %temp_title = title(sprintf('Location distribution correlated sequence (with model)'),'fontsize',16);
        %temp_title = title(sprintf(' Location distribution \n correlation (with model) '),'fontsize',12);
        %temp_title = title(sprintf('Location correlation'),'fontsize',10);
        temp_title = title(sprintf('Location-level'),'fontsize',9);
        temp_title.Interpreter = 'none';
        
    end
    
    
    %% Posterior_2d_n11n_mean
    if true
        fig = figure('Name','Posterior_2d_n11n_mean','NumberTitle','off'); %#ok<*NASGU>
        set(gcf,'Position',[310 50 137*0.95*1.3*1.1 396*0.9*1.02*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,7,'TileSpacing','loose','Padding','loose');
        
        temp_locDistri_confusionMatrix = Posterior_2d_n11n_mean;
        
        nexttile
        set(gca,'Visible','off');
        
        nexttile([1 5])
        
        C = temp_locDistri_confusionMatrix;
        C_max = 1;
        
        imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
        hold on
        
        for tempi=1:size(C,1)
            if ~isnan(C(tempi,1))
                continue
            end
            
            plot([1-0.4,size(C,2)+0.4],[tempi-0.3 tempi+0.3],...
                '-','LineWidth',1,'Color',[1 1 1]*0.7); %#ok<*UNRCH>
            hold on
        end
        
        
        set(gca,'linewidth',1.5)
        set(gca,'box','off');% 取消右、上边框
        set(gca, 'FontSize', 12) %12
        
        % Plot length bound in response
        for temptempi=1:size(temp_locDistri_confusionMatrix,1)
            if temptempi > 1
                if length(target_seqSet_inOne{temptempi}) > ...
                        length(target_seqSet_inOne{temptempi-1})
                    plot([temptempi temptempi]-0.5, [0.5 size(temp_locDistri_confusionMatrix,1)+0.5], 'color', [0.25 0.25 0.25]);
                    hold on
                end
            end
        end
        seqSet_inOne_inOne = [];
        for temptempi=1:length(seqSet_inOne)
            seqSet_inOne_inOne = [seqSet_inOne_inOne  seqSet_inOne{temptempi}];
        end
        
        set(gca,'YTick',1:size(temp_locDistri_confusionMatrix,1));
        
        ytl=string(seqSet_inOne_inOne(1:size(temp_locDistri_confusionMatrix,1)));
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        % 设置xtext的x坐标位置
        xtext_xp=xt;
        % 设置xtext的y坐标位置
        xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
        % 设置ttext的x坐标位置
        ytext_xp=(xt(1))*ones(1,length(yt))-1.6;%-0.75
        % 设置ttext的y坐标位置
        ytext_yp=yt;
        text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',6.5);%8,7
        set(gca,'yticklabel','');
        
        
        set(gca,'XTick',[1:1:numFrames]);
        %xtickangle(0);
        
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        
        xtl=string(1:1:numFrames);
        xtext_xp=xt;
        xtext_yp=(yt(end))*ones(1,length(xt))+1.35;
        
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',9);%8
        set(gca,'xticklabel','');
        
        
        set(gca,'TickLength',[0 0]);
        
        if if_colormap_loadEnhanced == 1
            load('parula_enhanced');
            colormap(parula_enhanced);
        elseif if_colormap_loadEnhanced == 0
            colormap parula
        end
        
        x_lim = xlim;
        y_lim = ylim;
        
        xlabel('Location','Position',[mean(x_lim) max(y_lim)+1.4],'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Sequence','Position',[min(x_lim)-1.15 mean(y_lim)],'FontSize',10,'FontWeight','bold');
        
        
        temp_title = title(sprintf(' Population decoder\n Correct'),'FontSize',9);
        temp_title.Interpreter = 'none';
        
    end
    
    %% Posterior_2d_n11n_resampleMean_stimuliError
    if true
        fig = figure('Name','Posterior_2d_n11n_resampleMean_stimuliError','NumberTitle','off');
        set(gcf,'Position',[510 50 137*0.95*1.3*1.1 396*0.9*1.02*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,7,'TileSpacing','loose','Padding','loose');
        
        temp_locDistri_confusionMatrix = Posterior_2d_n11n_resampleMean_stimuliError;
        
        nexttile
        set(gca,'Visible','off');
        
        nexttile([1 5])
        
        C = temp_locDistri_confusionMatrix;
        C_max = 1;
        
        imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
        hold on
        
        for tempi=1:size(C,1)
            if ~isnan(C(tempi,1))
                continue
            end
            
            plot([1-0.4,size(C,2)+0.4],[tempi-0.3 tempi+0.3],...
                '-','LineWidth',1,'Color',[1 1 1]*0.7); %#ok<*UNRCH>
            hold on
        end
        
        
        set(gca,'linewidth',1.5)
        set(gca,'box','off');% 取消右、上边框
        set(gca, 'FontSize', 12) %12
        
        % Plot length bound in response
        for temptempi=1:size(temp_locDistri_confusionMatrix,1)
            if temptempi > 1
                if length(target_seqSet_inOne{temptempi}) > ...
                        length(target_seqSet_inOne{temptempi-1})
                    plot([temptempi temptempi]-0.5, [0.5 size(temp_locDistri_confusionMatrix,1)+0.5], 'color', [0.25 0.25 0.25]);
                    hold on
                end
            end
        end
        seqSet_inOne_inOne = [];
        for temptempi=1:length(seqSet_inOne)
            seqSet_inOne_inOne = [seqSet_inOne_inOne  seqSet_inOne{temptempi}];
        end
        
        set(gca,'YTick',1:size(temp_locDistri_confusionMatrix,1));
        
        ytl=string(seqSet_inOne_inOne(1:size(temp_locDistri_confusionMatrix,1)));
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        % 设置xtext的x坐标位置
        xtext_xp=xt;
        % 设置xtext的y坐标位置
        xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
        % 设置ttext的x坐标位置
        ytext_xp=(xt(1))*ones(1,length(yt))-1.6;%-0.75
        % 设置ttext的y坐标位置
        ytext_yp=yt;
        text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',6.5);%8,7
        set(gca,'yticklabel','');
        
        
        set(gca,'XTick',[1:1:numFrames]);
        %xtickangle(0);
        
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        
        xtl=string(1:1:numFrames);
        xtext_xp=xt;
        xtext_yp=(yt(end))*ones(1,length(xt))+1.35;
        
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',9);%8
        set(gca,'xticklabel','');
        
        
        set(gca,'TickLength',[0 0]);
        
        if if_colormap_loadEnhanced == 1
            load('parula_enhanced');
            colormap(parula_enhanced);
        elseif if_colormap_loadEnhanced == 0
            colormap parula
        end
        
        x_lim = xlim;
        y_lim = ylim;
        
        xlabel('Location','Position',[mean(x_lim) max(y_lim)+1.4],'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Sequence','Position',[min(x_lim)-1.15 mean(y_lim)],'FontSize',10,'FontWeight','bold');
        
        
        temp_title = title(sprintf(' Population decoder\n Stimuli-error'),'FontSize',9);
        temp_title.Interpreter = 'none';
        
    end
    
    %% Posterior_2d_n11n_resampleMean_responseError
    if true
        fig = figure('Name','Posterior_2d_n11n_resampleMean_responseError','NumberTitle','off');
        set(gcf,'Position',[710 50 137*0.95*1.3*1.1 396*0.9*1.02*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,7,'TileSpacing','loose','Padding','loose');
        
        temp_locDistri_confusionMatrix = Posterior_2d_n11n_resampleMean_responseError;
        
        nexttile
        set(gca,'Visible','off');
        
        nexttile([1 5])
        
        C = temp_locDistri_confusionMatrix;
        C_max = 1;
        
        imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
        hold on
        
        for tempi=1:size(C,1)
            if ~isnan(C(tempi,1))
                continue
            end
            
            plot([1-0.4,size(C,2)+0.4],[tempi-0.3 tempi+0.3],...
                '-','LineWidth',1,'Color',[1 1 1]*0.7); %#ok<*UNRCH>
            hold on
        end
        
        
        set(gca,'linewidth',1.5)
        set(gca,'box','off');% 取消右、上边框
        set(gca, 'FontSize', 12) %12
        
        % Plot length bound in response
        for temptempi=1:size(temp_locDistri_confusionMatrix,1)
            if temptempi > 1
                if length(target_seqSet_inOne{temptempi}) > ...
                        length(target_seqSet_inOne{temptempi-1})
                    plot([temptempi temptempi]-0.5, [0.5 size(temp_locDistri_confusionMatrix,1)+0.5], 'color', [0.25 0.25 0.25]);
                    hold on
                end
            end
        end
        seqSet_inOne_inOne = [];
        for temptempi=1:length(seqSet_inOne)
            seqSet_inOne_inOne = [seqSet_inOne_inOne  seqSet_inOne{temptempi}];
        end
        
        set(gca,'YTick',1:size(temp_locDistri_confusionMatrix,1));
        
        ytl=string(seqSet_inOne_inOne(1:size(temp_locDistri_confusionMatrix,1)));
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        % 设置xtext的x坐标位置
        xtext_xp=xt;
        % 设置xtext的y坐标位置
        xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
        % 设置ttext的x坐标位置
        ytext_xp=(xt(1))*ones(1,length(yt))-1.6;%-0.75
        % 设置ttext的y坐标位置
        ytext_yp=yt;
        text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',6.5);%8,7
        set(gca,'yticklabel','');
        
        
        set(gca,'XTick',[1:1:numFrames]);
        %xtickangle(0);
        
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        
        xtl=string(1:1:numFrames);
        xtext_xp=xt;
        xtext_yp=(yt(end))*ones(1,length(xt))+1.35;
        
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',9);%8
        set(gca,'xticklabel','');
        
        
        set(gca,'TickLength',[0 0]);
        
        if if_colormap_loadEnhanced == 1
            load('parula_enhanced');
            colormap(parula_enhanced);
        elseif if_colormap_loadEnhanced == 0
            colormap parula
        end
        
        x_lim = xlim;
        y_lim = ylim;
        
        xlabel('Location','Position',[mean(x_lim) max(y_lim)+1.4],'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Sequence','Position',[min(x_lim)-1.15 mean(y_lim)],'FontSize',10,'FontWeight','bold');
        
        
        temp_title = title(sprintf(' Population decoder\n Response-error'),'FontSize',9);
        temp_title.Interpreter = 'none';
        
    end
    
    
    %% Posterior_2d_seqMean_trainLength1
    if true
        fig = figure('Name','Posterior_2d_seqMean_trainLength1','NumberTitle','off');
        set(gcf,'Position',[910 50 137*0.95*1.3*1.1 396*0.9*1.02*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,7,'TileSpacing','loose','Padding','loose');
        
        temp_locDistri_confusionMatrix = Posterior_2d_seqMean_trainLength1;
        
        nexttile
        set(gca,'Visible','off');
        
        nexttile([1 5])
        
        C = temp_locDistri_confusionMatrix;
        C_max = 1;
        
        imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
        hold on
        
        for tempi=1:size(C,1)
            if ~isnan(C(tempi,1))
                continue
            end
            
            plot([1-0.4,size(C,2)+0.4],[tempi-0.3 tempi+0.3],...
                '-','LineWidth',1,'Color',[1 1 1]*0.7); %#ok<*UNRCH>
            hold on
        end
        
        
        set(gca,'linewidth',1.5)
        set(gca,'box','off');% 取消右、上边框
        set(gca, 'FontSize', 12) %12
        
        % Plot length bound in response
        for temptempi=1:size(temp_locDistri_confusionMatrix,1)
            if temptempi > 1
                if length(target_seqSet_inOne{temptempi}) > ...
                        length(target_seqSet_inOne{temptempi-1})
                    plot([temptempi temptempi]-0.5, [0.5 size(temp_locDistri_confusionMatrix,1)+0.5], 'color', [0.25 0.25 0.25]);
                    hold on
                end
            end
        end
        seqSet_inOne_inOne = [];
        for temptempi=1:length(seqSet_inOne)
            seqSet_inOne_inOne = [seqSet_inOne_inOne  seqSet_inOne{temptempi}];
        end
        
        set(gca,'YTick',1:size(temp_locDistri_confusionMatrix,1));
        
        ytl=string(seqSet_inOne_inOne(1:size(temp_locDistri_confusionMatrix,1)));
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        % 设置xtext的x坐标位置
        xtext_xp=xt;
        % 设置xtext的y坐标位置
        xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
        % 设置ttext的x坐标位置
        ytext_xp=(xt(1))*ones(1,length(yt))-1.6;%-0.75
        % 设置ttext的y坐标位置
        ytext_yp=yt;
        text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',6.5);%8,7
        set(gca,'yticklabel','');
        
        
        set(gca,'XTick',[1:1:numFrames]);
        %xtickangle(0);
        
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        
        xtl=string(1:1:numFrames);
        xtext_xp=xt;
        xtext_yp=(yt(end))*ones(1,length(xt))+1.35;
        
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',9);%8
        set(gca,'xticklabel','');
        
        
        set(gca,'TickLength',[0 0]);
        
        if if_colormap_loadEnhanced == 1
            load('parula_enhanced');
            colormap(parula_enhanced);
        elseif if_colormap_loadEnhanced == 0
            colormap parula
        end
        
        x_lim = xlim;
        y_lim = ylim;
        
        xlabel('Location','Position',[mean(x_lim) max(y_lim)+1.4],'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Sequence','Position',[min(x_lim)-1.15 mean(y_lim)],'FontSize',10,'FontWeight','bold');
        
        
        temp_title = title(sprintf(' Population decoder\n Train length1'),'FontSize',9);
        temp_title.Interpreter = 'none';
        
    end
    
    %% Posterior_2d_seqMean_trainLength2
    if true
        fig = figure('Name','Posterior_2d_seqMean_trainLength2','NumberTitle','off');
        set(gcf,'Position',[1110 50 137*0.95*1.3*1.1 396*0.9*1.02*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,7,'TileSpacing','loose','Padding','loose');
        
        temp_locDistri_confusionMatrix = Posterior_2d_seqMean_trainLength2;
        
        nexttile
        set(gca,'Visible','off');
        
        nexttile([1 5])
        
        C = temp_locDistri_confusionMatrix;
        C_max = 1;
        
        imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
        hold on
        
        for tempi=1:size(C,1)
            if ~isnan(C(tempi,1))
                continue
            end
            
            plot([1-0.4,size(C,2)+0.4],[tempi-0.3 tempi+0.3],...
                '-','LineWidth',1,'Color',[1 1 1]*0.7); %#ok<*UNRCH>
            hold on
        end
        
        
        set(gca,'linewidth',1.5)
        set(gca,'box','off');% 取消右、上边框
        set(gca, 'FontSize', 12) %12
        
        % Plot length bound in response
        for temptempi=1:size(temp_locDistri_confusionMatrix,1)
            if temptempi > 1
                if length(target_seqSet_inOne{temptempi}) > ...
                        length(target_seqSet_inOne{temptempi-1})
                    plot([temptempi temptempi]-0.5, [0.5 size(temp_locDistri_confusionMatrix,1)+0.5], 'color', [0.25 0.25 0.25]);
                    hold on
                end
            end
        end
        seqSet_inOne_inOne = [];
        for temptempi=1:length(seqSet_inOne)
            seqSet_inOne_inOne = [seqSet_inOne_inOne  seqSet_inOne{temptempi}];
        end
        
        set(gca,'YTick',1:size(temp_locDistri_confusionMatrix,1));
        
        ytl=string(seqSet_inOne_inOne(1:size(temp_locDistri_confusionMatrix,1)));
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        % 设置xtext的x坐标位置
        xtext_xp=xt;
        % 设置xtext的y坐标位置
        xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
        % 设置ttext的x坐标位置
        ytext_xp=(xt(1))*ones(1,length(yt))-1.6;%-0.75
        % 设置ttext的y坐标位置
        ytext_yp=yt;
        text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',6.5);%8,7
        set(gca,'yticklabel','');
        
        
        set(gca,'XTick',[1:1:numFrames]);
        %xtickangle(0);
        
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        
        xtl=string(1:1:numFrames);
        xtext_xp=xt;
        xtext_yp=(yt(end))*ones(1,length(xt))+1.35;
        
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',9);%8
        set(gca,'xticklabel','');
        
        
        set(gca,'TickLength',[0 0]);
        
        if if_colormap_loadEnhanced == 1
            load('parula_enhanced');
            colormap(parula_enhanced);
        elseif if_colormap_loadEnhanced == 0
            colormap parula
        end
        
        x_lim = xlim;
        y_lim = ylim;
        
        xlabel('Location','Position',[mean(x_lim) max(y_lim)+1.4],'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Sequence','Position',[min(x_lim)-1.15 mean(y_lim)],'FontSize',10,'FontWeight','bold');
        
        
        temp_title = title(sprintf(' Population decoder\n Train length2'),'FontSize',9);
        temp_title.Interpreter = 'none';
        
    end
    
    
    %% Posterior_2d_seqMean_trainLength3
    if true
        fig = figure('Name','Posterior_2d_seqMean_trainLength3','NumberTitle','off');
        set(gcf,'Position',[1310 50 137*0.95*1.3*1.1 396*0.9*1.02*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,7,'TileSpacing','loose','Padding','loose');
        
        temp_locDistri_confusionMatrix = Posterior_2d_seqMean_trainLength3;
        
        nexttile
        set(gca,'Visible','off');
        
        nexttile([1 5])
        
        C = temp_locDistri_confusionMatrix;
        C_max = 1;
        
        imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
        hold on
        
        for tempi=1:size(C,1)
            if ~isnan(C(tempi,1))
                continue
            end
            
            plot([1-0.4,size(C,2)+0.4],[tempi-0.3 tempi+0.3],...
                '-','LineWidth',1,'Color',[1 1 1]*0.7); %#ok<*UNRCH>
            hold on
        end
        
        
        set(gca,'linewidth',1.5)
        set(gca,'box','off');% 取消右、上边框
        set(gca, 'FontSize', 12) %12
        
        % Plot length bound in response
        for temptempi=1:size(temp_locDistri_confusionMatrix,1)
            if temptempi > 1
                if length(target_seqSet_inOne{temptempi}) > ...
                        length(target_seqSet_inOne{temptempi-1})
                    plot([temptempi temptempi]-0.5, [0.5 size(temp_locDistri_confusionMatrix,1)+0.5], 'color', [0.25 0.25 0.25]);
                    hold on
                end
            end
        end
        seqSet_inOne_inOne = [];
        for temptempi=1:length(seqSet_inOne)
            seqSet_inOne_inOne = [seqSet_inOne_inOne  seqSet_inOne{temptempi}];
        end
        
        set(gca,'YTick',1:size(temp_locDistri_confusionMatrix,1));
        
        ytl=string(seqSet_inOne_inOne(1:size(temp_locDistri_confusionMatrix,1)));
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        % 设置xtext的x坐标位置
        xtext_xp=xt;
        % 设置xtext的y坐标位置
        xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
        % 设置ttext的x坐标位置
        ytext_xp=(xt(1))*ones(1,length(yt))-1.6;%-0.75
        % 设置ttext的y坐标位置
        ytext_yp=yt;
        text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',6.5);%8,7
        set(gca,'yticklabel','');
        
        
        set(gca,'XTick',[1:1:numFrames]);
        %xtickangle(0);
        
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        
        xtl=string(1:1:numFrames);
        xtext_xp=xt;
        xtext_yp=(yt(end))*ones(1,length(xt))+1.35;
        
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',9);%8
        set(gca,'xticklabel','');
        
        
        set(gca,'TickLength',[0 0]);
        
        if if_colormap_loadEnhanced == 1
            load('parula_enhanced');
            colormap(parula_enhanced);
        elseif if_colormap_loadEnhanced == 0
            colormap parula
        end
        
        x_lim = xlim;
        y_lim = ylim;
        
        xlabel('Location','Position',[mean(x_lim) max(y_lim)+1.4],'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Sequence','Position',[min(x_lim)-1.15 mean(y_lim)],'FontSize',10,'FontWeight','bold');
        
        
        temp_title = title(sprintf(' Population decoder\n Train length3'),'FontSize',9);
        temp_title.Interpreter = 'none';
        
    end
    
    
    
    %% Length 1 generalization
    if if_correlateWith_theoretical0_behavior1 == 0
        C = double(boolIndex_location_seq_T(1:sum(numSeq(1:3)),:));
    elseif if_correlateWith_theoretical0_behavior1 == 1
        C = Posterior_2d_model;
    end
    
    Posterior_2d_seqMean_trainLength1;
    Posterior_2d_seqMean_trainLength2;
    Posterior_2d_seqMean_trainLength3;
    
    
    temp_posterior = Posterior_2d_seqMean_trainLength1;
    temp_r_n11n = zeros(1,sum(numSeq(1:valid_length)));
    for tempi=1:sum(numSeq(1:valid_length))
        [temp_r_n11n(tempi),~] = corr(temp_posterior(tempi,:)',C(tempi,:)');
    end
    
    target_length = 1;
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    temp_1 = temp_r_n11n(temp_range);
    
    target_length = 2;
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    temp_2 = temp_r_n11n(temp_range);
    
    target_length = 3;
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    temp_3 = temp_r_n11n(temp_range);
    
    temp_1_chanceLevel = temp_r_n11n_resampleMean_correct_chanceLevel(~isnan(temp_r_n11n_resampleMean_correct_chanceLevel));
    
    [~,temp_p1] = ttest2(temp_1,temp_1_chanceLevel);
    [~,temp_p2] = ttest2(temp_2,temp_1_chanceLevel);
    [~,temp_p3] = ttest2(temp_3,temp_1_chanceLevel);
    
    
    if true
        fig = figure('Name','Length 1 generalization','NumberTitle','off');
        set(gcf,'Position',[350+0 550+0 240*0.80*0.78 336*1.11*0.78*0.9*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        
        temp_y_min = mean(temp_1_chanceLevel);
        temp_y_max = max([temp_1,temp_2,temp_3]);
        
        temp_data = [temp_1';temp_2';temp_3'];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'C'},length(temp_3),1);
        
        temp_label = [g1;g2;g3];
        
        temptemp_color1 = [1 1 1]*0.5;%0.5
        temptemp_color2 = repmat(temptemp_color1, 3, 1);
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        h(3).ViolinPlot.FaceAlpha = 0.1;
        
        plot([0 4],mean(temp_1_chanceLevel)*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        tempTxt = sprintf('');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y_max+(temp_y_max-temp_y_min)*0.13,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        
        tempTxt = sprintf('');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(2,temp_y_max+(temp_y_max-temp_y_min)*0.13,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        
        tempTxt = sprintf('');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        text(3,temp_y_max+(temp_y_max-temp_y_min)*0.13,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        
        set(gca,'linewidth',1.5)
        
        %xlim([0.4 4]);
        xlim([0.5 3.5]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.3]);
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        xtl = ["Len1-Len1", "Len1-Len2", "Len1-Len3"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.19;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7);%9-->12-->11
        
        set(gca,'xticklabel','');
        
        yticks([0 1]);
        
        ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('Location-level'),'fontsize',9);
        subtitle(sprintf('Length 1 generalization'),'fontsize',7.5);
        
    end
    
    
    %% Length 2 generalization
    if if_correlateWith_theoretical0_behavior1 == 0
        C = double(boolIndex_location_seq_T(1:sum(numSeq(1:3)),:));
    elseif if_correlateWith_theoretical0_behavior1 == 1
        C = Posterior_2d_model;
    end
    
    Posterior_2d_seqMean_trainLength1;
    Posterior_2d_seqMean_trainLength2;
    Posterior_2d_seqMean_trainLength3;
    
    
    temp_posterior = Posterior_2d_seqMean_trainLength2;
    temp_r_n11n = zeros(1,sum(numSeq(1:valid_length)));
    for tempi=1:sum(numSeq(1:valid_length))
        [temp_r_n11n(tempi),~] = corr(temp_posterior(tempi,:)',C(tempi,:)');
    end
    
    target_length = 1;
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    temp_1 = temp_r_n11n(temp_range);
    
    target_length = 2;
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    temp_2 = temp_r_n11n(temp_range);
    
    target_length = 3;
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    temp_3 = temp_r_n11n(temp_range);
    
    temp_1_chanceLevel = temp_r_n11n_resampleMean_correct_chanceLevel(~isnan(temp_r_n11n_resampleMean_correct_chanceLevel));
    
    [~,temp_p1] = ttest2(temp_1,temp_1_chanceLevel);
    [~,temp_p2] = ttest2(temp_2,temp_1_chanceLevel);
    [~,temp_p3] = ttest2(temp_3,temp_1_chanceLevel);
    
    
    if true
        fig = figure('Name','Length 2 generalization','NumberTitle','off');
        set(gcf,'Position',[550+0 550+0 240*0.80*0.78 336*1.11*0.78*0.9*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        
        temp_y_min = mean(temp_1_chanceLevel);
        temp_y_max = max([temp_1,temp_2,temp_3]);
        
        temp_data = [temp_1';temp_2';temp_3'];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'C'},length(temp_3),1);
        
        temp_label = [g1;g2;g3];
        
        temptemp_color1 = [1 1 1]*0.5;%0.5
        temptemp_color2 = repmat(temptemp_color1, 3, 1);
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        h(3).ViolinPlot.FaceAlpha = 0.1;
        
        plot([0 4],mean(temp_1_chanceLevel)*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        tempTxt = sprintf('');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y_max+(temp_y_max-temp_y_min)*0.13,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        
        tempTxt = sprintf('');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(2,temp_y_max+(temp_y_max-temp_y_min)*0.13,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        
        tempTxt = sprintf('');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        text(3,temp_y_max+(temp_y_max-temp_y_min)*0.13,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        
        set(gca,'linewidth',1.5)
        
        %xlim([0.4 4]);
        xlim([0.5 3.5]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.3]);
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        xtl = ["Len2-Len1", "Len2-Len2", "Len2-Len3"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.19;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7);%9-->12-->11
        
        set(gca,'xticklabel','');
        
        yticks([0 1]);
        
        ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('Location-level'),'fontsize',9);
        subtitle(sprintf('Length 2 generalization'),'fontsize',7.5);
        
    end
    
    
    %% Length 3 generalization
    if if_correlateWith_theoretical0_behavior1 == 0
        C = double(boolIndex_location_seq_T(1:sum(numSeq(1:3)),:));
    elseif if_correlateWith_theoretical0_behavior1 == 1
        C = Posterior_2d_model;
    end
    
    Posterior_2d_seqMean_trainLength1;
    Posterior_2d_seqMean_trainLength2;
    Posterior_2d_seqMean_trainLength3;
    
    
    temp_posterior = Posterior_2d_seqMean_trainLength3;
    temp_r_n11n = zeros(1,sum(numSeq(1:valid_length)));
    for tempi=1:sum(numSeq(1:valid_length))
        [temp_r_n11n(tempi),~] = corr(temp_posterior(tempi,:)',C(tempi,:)');
    end
    
    target_length = 1;
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    temp_1 = temp_r_n11n(temp_range);
    
    target_length = 2;
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    temp_2 = temp_r_n11n(temp_range);
    
    target_length = 3;
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
    temp_3 = temp_r_n11n(temp_range);
    
    temp_1_chanceLevel = temp_r_n11n_resampleMean_correct_chanceLevel(~isnan(temp_r_n11n_resampleMean_correct_chanceLevel));
    
    [~,temp_p1] = ttest2(temp_1,temp_1_chanceLevel);
    [~,temp_p2] = ttest2(temp_2,temp_1_chanceLevel);
    [~,temp_p3] = ttest2(temp_3,temp_1_chanceLevel);
    
    
    if true
        fig = figure('Name','Length 3 generalization','NumberTitle','off');
        set(gcf,'Position',[750+0 550+0 240*0.80*0.78 336*1.11*0.78*0.9*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        
        temp_y_min = mean(temp_1_chanceLevel);
        temp_y_max = max([temp_1,temp_2,temp_3]);
        
        temp_data = [temp_1';temp_2';temp_3'];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'C'},length(temp_3),1);
        
        temp_label = [g1;g2;g3];
        
        temptemp_color1 = [1 1 1]*0.5;%0.5
        temptemp_color2 = repmat(temptemp_color1, 3, 1);
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        h(3).ViolinPlot.FaceAlpha = 0.1;
        
        plot([0 4],mean(temp_1_chanceLevel)*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        tempTxt = sprintf('');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y_max+(temp_y_max-temp_y_min)*0.13,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        
        tempTxt = sprintf('');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(2,temp_y_max+(temp_y_max-temp_y_min)*0.13,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        
        tempTxt = sprintf('');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        text(3,temp_y_max+(temp_y_max-temp_y_min)*0.13,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        
        set(gca,'linewidth',1.5)
        
        %xlim([0.4 4]);
        xlim([0.5 3.5]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.3]);
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        xtl = ["Len3-Len1", "Len3-Len2", "Len3-Len3"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.19;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7);%9-->12-->11
        
        set(gca,'xticklabel','');
        
        yticks([0 1]);
        
        ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('Location-level'),'fontsize',9);
        subtitle(sprintf('Length 3 generalization'),'fontsize',7.5);
        
    end
    
    
end