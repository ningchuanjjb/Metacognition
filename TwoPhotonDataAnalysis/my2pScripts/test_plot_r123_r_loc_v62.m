% Chuan's 5th script (20251214)
% This script: One part of 'stepD3_train123_test123.m'.
close all

if_plot = 1;

valid_length = 3;

if_plot_locCorr = 0;
if_plot_seqCorr = 0;
if_plot_eachSeq = 0;
if_plot_eachSeq_confusionMatrixB = 0;
if_plot_eachSeq_confusionMatrixC = 0;
if_plot_eachSeq_confusionMatrixD = 0;
if_plot_eachSeq_confusionMatrixE = 1;

if_plot_seqDecode_errorRate0_trialCount1 = 0;

if_plot_eachSeq_singleFig0_multiFig1 = 0;

if_colormap_loadEnhanced = 0;

if_gAcc_noChoice0_allMemory1_choiceMemory2 = if_gAcc_noChoice0_allMemory1_choiceMemory2;

if_locationDistri_behavior0_model1 = 0;%1


Posterior_2d_behavior = nan(size(Posterior_2d_model));

% temp_gAcc = gAcc_target_collapsed_inOne(1:sum(numSeq(1:valid_length)));

responseMatrix_noChoice;
responseMatrix_noChoice_valid = responseMatrix_noChoice(1:sum(numSeq(1:4)),1:sum(numSeq(1:4)));
responseMatrix_noChoice_valid_n11n = responseMatrix_noChoice_valid ./ sum(responseMatrix_noChoice_valid,2);

responseMatrix_allMemory;
responseMatrix_allMemory_valid = responseMatrix_allMemory(1:sum(numSeq(1:4)),1:sum(numSeq(1:4)));
responseMatrix_allMemory_valid_n11n = responseMatrix_allMemory_valid ./ sum(responseMatrix_allMemory_valid,2);

responseMatrix_choice;
responseMatrix_choice_valid = responseMatrix_choice(1:sum(numSeq(1:4)),1:sum(numSeq(1:4)));
responseMatrix_choice_valid_n11n = responseMatrix_choice_valid ./ sum(responseMatrix_choice_valid,2);


if if_gAcc_noChoice0_allMemory1_choiceMemory2 == 0
    temp_responseMatrix_valid_n11n = responseMatrix_noChoice_valid_n11n(1:sum(numSeq(1:valid_length)),:);
elseif if_gAcc_noChoice0_allMemory1_choiceMemory2 == 1
    temp_responseMatrix_valid_n11n = responseMatrix_allMemory_valid_n11n(1:sum(numSeq(1:valid_length)),:);
    
    %temp_responseMatrix_valid_n11n = responseMatrix_choice_valid_n11n(1:sum(numSeq(1:valid_length)),:);
elseif if_gAcc_noChoice0_allMemory1_choiceMemory2 == 2
    %temp_responseMatrix_valid_n11n = responseMatrix_choice_valid_n11n(1:sum(numSeq(1:valid_length)),:);
    temp_responseMatrix_valid_n11n = responseMatrix_noChoice_valid_n11n(1:sum(numSeq(1:valid_length)),:);
end


temp_responseMatrix_valid_n11n;

boolIndex_location_seq_T;

for tempi=1:size(Posterior_2d_behavior,1)
    temp_seqDistri = temp_responseMatrix_valid_n11n(tempi,:);
    
    temp_locDistri = temp_seqDistri * boolIndex_location_seq_T;
    
    %temp_locDistri_n11n = temp_locDistri ./ sum(temp_locDistri) .* sum(boolIndex_location_seq_T(tempi,:));
    %Posterior_2d_behavior(tempi,:) = temp_locDistri_n11n;
    
    Posterior_2d_behavior(tempi,:) = temp_locDistri;
end
% Posterior_2d_behavior = Posterior_2d_behavior ./ sum(Posterior_2d_behavior,2);

if if_locationDistri_behavior0_model1 == 0
    Posterior_2d_model = Posterior_2d_behavior; %#ok<*ASGSL>
elseif if_locationDistri_behavior0_model1 == 1
    Posterior_2d_model = Posterior_2d_model;
end


%% Others
x = [svm_posterior_length1(~isnan(svm_posterior_length1))';svm_posterior_length2(~isnan(svm_posterior_length2))';svm_posterior_length3(~isnan(svm_posterior_length3))'];
y = [gAcc_length1(~isnan(svm_posterior_length1))';gAcc_length2(~isnan(svm_posterior_length2))';gAcc_length3(~isnan(svm_posterior_length3))'];
[r_123,p_123] = corr(x,y);

r_loc_n11n = zeros(1,numFrames);
p_loc_n11n = zeros(1,numFrames);
for tempi=1:numFrames
    tempBoolIndex = ~isnan(Posterior_2d_n11n_mean(:,1));
    [r_loc_n11n(tempi),p_loc_n11n(tempi)] = ...
        corr(Posterior_2d_n11n_mean(tempBoolIndex,tempi),Posterior_2d_model(tempBoolIndex,tempi));
end

% [temp_r,temp_p]=corr(gAcc_noChoice_inOne_model',gAcc_noChoice_collapsed_inOne');


r_n11n = zeros(1,sum(numSeq(1:valid_length)));
p_n11n = zeros(1,sum(numSeq(1:valid_length)));
for tempi=1:sum(numSeq(1:valid_length))
    [r_n11n(tempi),p_n11n(tempi)] = corr(Posterior_2d_n11n_mean(tempi,:)',Posterior_2d_model(tempi,:)');
end

% fprintf('num(p_posterior_seq_n11n<0.05)=%d/%d.\n',sum(p_n11n<0.05),sum(~isnan(p_n11n)));


seqTrialCount = zeros(sum(numSeq(1:valid_length)),1);
for tempi=1:size(seqTrialCount,1)
    seqTrialCount(tempi) = sum(seqIndex_valid==tempi);
end

seqTrialCount_valid = seqTrialCount(seqTrialCount>0);


%% temptempIndex = y>0.99 & x<0.2;
x = [svm_posterior_length1(~isnan(svm_posterior_length1))';svm_posterior_length2(~isnan(svm_posterior_length2))';svm_posterior_length3(~isnan(svm_posterior_length3))'];
y = [gAcc_length1(~isnan(svm_posterior_length1))';gAcc_length2(~isnan(svm_posterior_length2))';gAcc_length3(~isnan(svm_posterior_length3))'];

temptempIndex = y>0.99 & x<0.2;
x = x(~temptempIndex);
y = y(~temptempIndex);

[r,p] = corr(x,y);

a = 1;



%% locDistri_confusionMatrix_A
r_n11n;
locDistri_confusionMatrix_A = zeros(length(r_n11n)*numFrames,length(r_n11n)*numFrames);
% locDistri_confusionMatrix_A = locDistri_confusionMatrix_A + rand(size(locDistri_confusionMatrix_A));
%
% locDistri_confusionMatrix_A(1:10,1:100) = 0.5; % x is neuron decoder, y is behavior model


Posterior_2d_n11n_mean;
Posterior_2d_model;

for tempi=1:sum(numSeq(1:valid_length))
    temp_range = ((tempi-1)*numFrames+1):(tempi*numFrames);
    
    temp_posterior_neuron = Posterior_2d_n11n_mean(tempi,:);
    temp_posterior_model = Posterior_2d_model(tempi,:);
    
    temp_subMatrix = nan(numFrames,numFrames); % x is neuron decoder, y is behavior model
    for tempj=1:numFrames
        %temp_subMatrix(tempj,:) = 1-abs(temp_posterior_neuron(tempj) - temp_posterior_model);
        temp_subMatrix(tempj,:) = temp_posterior_neuron(tempj) .* temp_posterior_model;
    end
    
    locDistri_confusionMatrix_A(temp_range,temp_range) = temp_subMatrix;
end




%% locDistri_confusionMatrix_B
r_n11n;
% locDistri_confusionMatrix_B = zeros(length(r_n11n)*numFrames,length(r_n11n)*numFrames);
% locDistri_confusionMatrix_B = locDistri_confusionMatrix_B + rand(size(locDistri_confusionMatrix_B));
%
% locDistri_confusionMatrix_B(1:10,1:100) = 0.5; % x is neuron decoder, y is behavior model

locDistri_confusionMatrix_B = zeros(length(r_n11n),numFrames,2);

locDistri_confusionMatrix_B(:,:,1) = Posterior_2d_model;
locDistri_confusionMatrix_B(:,:,2) = Posterior_2d_n11n_mean;

% Posterior_2d_n11n_mean;
% Posterior_2d_model;
%
% for tempi=1:sum(numSeq(1:valid_length))
%     temp_range = ((tempi-1)*numFrames+1):(tempi*numFrames);
%
%     temp_posterior_neuron = Posterior_2d_n11n_mean(tempi,:);
%     temp_posterior_model = Posterior_2d_model(tempi,:);
%
%     temp_subMatrix = nan(numFrames,numFrames); % x is neuron decoder, y is behavior model
%     for tempj=1:numFrames
%         %temp_subMatrix(tempj,:) = 1-abs(temp_posterior_neuron(tempj) - temp_posterior_model);
%         temp_subMatrix(tempj,:) = temp_posterior_neuron(tempj) .* temp_posterior_model;
%     end
%
%     locDistri_confusionMatrix_A(temp_range,temp_range) = temp_subMatrix;
% end


%% Get r_n11n_resampleMean_correct
Posterior_2d_n11n_resample = zeros(resampleIterCount,size(Posterior_2d_n11n_mean,1),size(Posterior_2d_n11n_mean,2));
r_n11n_resample_correct = zeros(resampleIterCount,sum(numSeq(1:valid_length)));
p_n11n_resample_correct = zeros(resampleIterCount,sum(numSeq(1:valid_length)));
r_n11n_resample_correct_chanceLevel = zeros(resampleIterCount,sum(numSeq(1:valid_length)));
for tempi=1:resampleIterCount
    temp_posterior = [];
    temp_posterior = [temp_posterior; svm_train_length1_outputs.temp_svm_resample{tempi}.Posterior_2d_n11n_lengthx_mean]; %#ok<*AGROW>
    temp_posterior = [temp_posterior; svm_train_length2_outputs.temp_svm_resample{tempi}.Posterior_2d_n11n_lengthx_mean];
    temp_posterior = [temp_posterior; svm_train_length3_outputs.temp_svm_resample{tempi}.Posterior_2d_n11n_lengthx_mean];
    
    Posterior_2d_n11n_resample(tempi,:,:) = temp_posterior;
    
    temp_r_n11n = zeros(1,sum(numSeq(1:valid_length)));
    temp_p_n11n = zeros(1,sum(numSeq(1:valid_length)));
    temp_r_n11n_chanceLevel = zeros(1,sum(numSeq(1:valid_length)));
    for tempj=1:sum(numSeq(1:valid_length))
        [temp_r_n11n(tempj),temp_p_n11n(tempj)] = corr(temp_posterior(tempj,:)',Posterior_2d_model(tempj,:)');
        [temp_r_n11n_chanceLevel(tempj),~] = corr(rand(numFrames,1),Posterior_2d_model(tempj,:)');
    end
    r_n11n_resample_correct(tempi,:) = temp_r_n11n;
    p_n11n_resample_correct(tempi,:) = temp_p_n11n;
    r_n11n_resample_correct_chanceLevel(tempi,:) = temp_r_n11n_chanceLevel;
end
r_n11n_resampleMean_correct = mean(r_n11n_resample_correct,1);
p_n11n_resampleMean_correct = mean(p_n11n_resample_correct,1);
r_n11n_resampleMean_correct_chanceLevel = mean(r_n11n_resample_correct_chanceLevel,1);


%% Get r_n11n_resampleMean_correct
% Posterior_2d_n11n_resample = zeros(resampleIterCount,size(Posterior_2d_n11n_mean,1),size(Posterior_2d_n11n_mean,2));
r_n11n_resample_correct_theoretical = zeros(resampleIterCount,sum(numSeq(1:valid_length)));
p_n11n_resample_correct_theoretical = zeros(resampleIterCount,sum(numSeq(1:valid_length)));
r_n11n_resample_correct_theoretical_chanceLevel = zeros(resampleIterCount,sum(numSeq(1:valid_length)));

Posterior_2d_theoretical = double(boolIndex_location_seq_T(1:sum(numSeq(1:3)),:));

for tempi=1:resampleIterCount
    temp_posterior = [];
    temp_posterior = [temp_posterior; svm_train_length1_outputs.temp_svm_resample{tempi}.Posterior_2d_n11n_lengthx_mean]; %#ok<*AGROW>
    temp_posterior = [temp_posterior; svm_train_length2_outputs.temp_svm_resample{tempi}.Posterior_2d_n11n_lengthx_mean];
    temp_posterior = [temp_posterior; svm_train_length3_outputs.temp_svm_resample{tempi}.Posterior_2d_n11n_lengthx_mean];
    
    %Posterior_2d_n11n_resample(tempi,:,:) = temp_posterior;
    
    temp_r_n11n = zeros(1,sum(numSeq(1:valid_length)));
    temp_p_n11n = zeros(1,sum(numSeq(1:valid_length)));
    temp_r_n11n_chanceLevel = zeros(1,sum(numSeq(1:valid_length)));
    for tempj=1:sum(numSeq(1:valid_length))
        [temp_r_n11n(tempj),temp_p_n11n(tempj)] = corr(temp_posterior(tempj,:)',Posterior_2d_theoretical(tempj,:)');
        [temp_r_n11n_chanceLevel(tempj),~] = corr(rand(numFrames,1),Posterior_2d_theoretical(tempj,:)');
    end
    r_n11n_resample_correct_theoretical(tempi,:) = temp_r_n11n;
    p_n11n_resample_correct_theoretical(tempi,:) = temp_p_n11n;
    r_n11n_resample_correct_theoretical_chanceLevel(tempi,:) = temp_r_n11n_chanceLevel;
end
r_n11n_resampleMean_correct_theoretical = mean(r_n11n_resample_correct_theoretical,1);
p_n11n_resampleMean_correct_theoretical = mean(p_n11n_resample_correct_theoretical,1);
r_n11n_resampleMean_correct_theoretical_chanceLevel = mean(r_n11n_resample_correct_theoretical_chanceLevel,1);


%% Plot
if if_plot == 1
    close all
    
    if if_plot_locCorr == 1
        %% Fig, location distribution
        fig = figure('Name','location distribution','NumberTitle','off');
        set(gcf,'Position',[10 50 1100 700]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
        
        tempBoolIndex = ~isnan(Posterior_2d_n11n_mean(:,1));
        
        t.Title.String = sprintf('%d seqs (length123), size represents length',sum(tempBoolIndex));
        t.Title.FontSize = 20;
        t.Title.Interpreter = 'none';
        
        for tempi=1:numFrames
            nexttile
            
            %tempBoolIndex = ~isnan(Posterior_2d_n11n_mean(:,1));
            
            x = Posterior_2d_n11n_mean(tempBoolIndex,tempi);
            y = Posterior_2d_model(tempBoolIndex,tempi);
            [r,p] = corr(x,y);
            
            mdl = fitglm(x,y);
            %r2 = mdl.Rsquared.Adjusted;
            
            temp_x = Posterior_2d_n11n_mean(:,tempi);
            temp_y = Posterior_2d_model(:,tempi);
            
            h = [];
            for tempj=1:valid_length
                temp_range = (sum(numSeq(1:tempj-1))+1):sum(numSeq(1:tempj));
                tempIndex = find(tempBoolIndex(temp_range)==true);
                temp_range2 = temp_range(tempIndex);
                
                temp_size = ((tempj.^3)*2 + 3) .* ones(1, length(temp_range2));
                temp_h = scatter(temp_x(temp_range2), temp_y(temp_range2), ...
                    temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
                hold on
                h = [h temp_h]; %#ok<*AGROW>
            end
            
            [temp_xmin,temp_xmax] = bounds(x);
            [temp_ymin,temp_ymax] = bounds(y);
            
            x_fit = temp_xmin:0.001:temp_xmax;
            y_fit = predict(mdl,x_fit')';
            
            temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
            hold on
            h = [h temp_h];
            
            xlim([0 1]);
            ylim([0 1]);
            
            set(gca, 'FontSize', 11)
            set(gca,'box','off');% 取消右、上边框
            xlabel(sprintf('neuron decoder'), 'FontSize', 14, 'FontWeight', 'bold');
            if if_locationDistri_behavior0_model1 == 0
                ylabel(sprintf('behavior'), 'FontSize', 14, 'FontWeight', 'bold');
            elseif if_locationDistri_behavior0_model1 == 1
                ylabel(sprintf('behavior model'), 'FontSize', 14, 'FontWeight', 'bold');
            end
            title(sprintf('location %d, r=%.3f, p=%.3f',tempi,r,p), 'FontSize', 14, 'FontWeight', 'bold');
            
        end
    end
    
    if if_plot_seqCorr == 1
        %% fig, Sequence decoding
        fig = figure('Name','Sequence decoding','NumberTitle','off'); %#ok<*NASGU>
        %set(gcf,'Position',[50+0 50+0 600 450]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 600 320]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 240 252]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[50+0 50+0 240 252*0.925]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        x = [svm_posterior_length1(~isnan(svm_posterior_length1))';svm_posterior_length2(~isnan(svm_posterior_length2))';svm_posterior_length3(~isnan(svm_posterior_length3))'];
        if if_plot_seqDecode_errorRate0_trialCount1 == 0
            %y = 1-[gAcc_length1';gAcc_length2(~isnan(svm_posterior_length2))';gAcc_length3(~isnan(svm_posterior_length3))'];
            y = [gAcc_length1(~isnan(svm_posterior_length1))';gAcc_length2(~isnan(svm_posterior_length2))';gAcc_length3(~isnan(svm_posterior_length3))'];
        elseif if_plot_seqDecode_errorRate0_trialCount1 == 1
            y =  seqTrialCount_valid;
            %y(1:numFrames) = 20;
        end
        
        %         temptempIndex = y<0.99;
        %         x = x(temptempIndex);
        %         y = y(temptempIndex);
        
        [r,p] = corr(x,y);
        
        mdl = fitglm(x,y);
        %r2 = mdl.Rsquared.Adjusted;
        
        %scatter(x, y, 15, 'filled','MarkerFaceColor', [0 0 0], 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.8);%0.4
        %hold on
        
        temp_x = [svm_posterior_length1';svm_posterior_length2';svm_posterior_length3'];
        %temp_x = seqTrialCount;
        if if_plot_seqDecode_errorRate0_trialCount1 == 0
            %temp_y = 1-[gAcc_length1';gAcc_length2';gAcc_length3'];
            temp_y = [gAcc_length1';gAcc_length2';gAcc_length3'];
        elseif if_plot_seqDecode_errorRate0_trialCount1 == 1
            temp_y = seqTrialCount;
            %temp_y(1:numFrames) = 20;
        end
        tempBoolIndex = ~isnan(temp_x);
        
        h = [];
        for tempi=1:valid_length
            %for tempi=1:2
            temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
            tempIndex = find(tempBoolIndex(temp_range)==true);
            temp_range2 = temp_range(tempIndex);
            
            temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
            temp_h = scatter(temp_x(temp_range2), temp_y(temp_range2), ...
                temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
            h = [h temp_h];
        end
        
        [temp_xmin,temp_xmax] = bounds(x);
        [temp_ymin,temp_ymax] = bounds(y);
        
        x_fit = temp_xmin:0.001:temp_xmax;
        y_fit = predict(mdl,x_fit')';
        
        temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        hold on
        h = [h temp_h];
        
        %legend(h,'length 1','length 2','length 3','fit','Location','northeast','fontsize',12)
        
        xlim([0 1.1]);
        if if_plot_seqDecode_errorRate0_trialCount1 == 0
            ylim([0 1.1]);
        end
        
        set(gca,'linewidth',1.5)
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        set(gca,'XTick',0:0.5:1,'FontSize', 12);%给坐标加标签
        set(gca,'YTick',0:0.5:1,'FontSize', 12);%给坐标加标签
        
        if if_decoder_decodingAcc0_pProd1 == 0
            xlabel(sprintf('decoding accuracy'), 'FontSize', 12, 'FontWeight', 'bold');
        elseif if_decoder_decodingAcc0_pProd1== 1
            xlabel(sprintf('decoding p production'), 'FontSize', 12, 'FontWeight', 'bold');
        end
        %title(sprintf('%d seqs (length123), r=%.3f, p=%.3f\nsize represents length',...
        %    length(y),r,p), 'FontSize', 14, 'FontWeight', 'bold');
        %temp_title = title(sprintf('%d seqs (length123), r=%.3f, p=%.3f\n%s',length(y),r,p,currentSession), 'FontSize', 12, 'FontWeight', 'bold');
        %temp_title = title(sprintf('%d seqs, r=%.3f, p=%.3f\n%s',length(y),r,p,currentSession), 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('%s\n%d seqs, r=%.3f, p=%.3f',currentSession,length(y),r,p), 'FontSize', 12, 'FontWeight', 'bold');
        temp_title.Interpreter = 'none';
        
        if if_plot_seqDecode_errorRate0_trialCount1 == 0
            %ylabel(sprintf('Error rate'), 'FontSize', 14, 'FontWeight', 'bold');
            ylabel(sprintf('Accuracy'), 'FontSize', 12, 'FontWeight', 'bold');
        elseif if_plot_seqDecode_errorRate0_trialCount1 == 1
            ylabel(sprintf('Trial count'), 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
    
end

%% Fig, location distribution of each sequence
%color_choiceOffload = [0.6350 0.0780 0.1840];

if if_plot_eachSeq == 1
    if if_plot_eachSeq_singleFig0_multiFig1 == 1
        % Abandon now
    elseif if_plot_eachSeq_singleFig0_multiFig1 == 0
        fig = figure('Name','locDistri','NumberTitle','off');
        %set(gcf,'Position',[10 50 1300 600]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 50 1300 700]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 50 1100 600]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout('flow','TileSpacing','Compact','Padding','Compact');
        %t = tiledlayout(6,6,'TileSpacing','Compact','Padding','Compact');
        
        t.Title.String = sprintf('num(p_posterior_seq_n11n<0.05)=%d/%d.\n',sum(p_n11n<0.05),sum(~isnan(p_n11n)));
        t.Title.FontSize = 16;
        t.Title.Interpreter = 'none';
        
        temp_range = 1:sum(numSeq(1:valid_length));
        for tempi=1:length(temp_range)
            if sum(isnan(Posterior_2d_n11n_mean(tempi,:)))>0
                continue
            end
            nexttile
            
            templd1 = Posterior_2d_model(tempi,:);
            templd2 = Posterior_2d_n11n_mean(tempi,:);
            
            %plot(templd1,'color',[0 0 0],'LineWidth',1.5);
            %plot(templd1,'color',[0.7 0.7 0.7],'LineWidth',1.5);
            plot(templd1,'color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
            hold on
            %plot(templd2,'color',[0.75 0.25 0.25],'LineWidth',1.5);
            %plot(templd2,'color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
            plot(templd2,'color',[0 0 0],'LineWidth',1.5);
            hold on
            
            if tempi == 1
                if if_locationDistri_behavior0_model1 == 0
                    hl = legend('behavior','neuron decoder',...
                        'Location','northeast','fontsize',9); %5
                elseif if_locationDistri_behavior0_model1 == 1
                    hl = legend('behavior model','neuron decoder',...
                        'Location','northeast','fontsize',9); %5
                end
                hl.ItemTokenSize = ones(1,2)*16;
            end
            ylim([0 1]);
            set(gca, 'FontSize', 10)
            set(gca,'XLim',[1-0.5 numFrames+0.5]);%X轴的数据显示范围
            set(gca,'XTick',[1:1:numFrames]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
            %set(gca,'XTickLabel',[1:1:numFrames]);%给坐标加标签
            
            tempSeqBoolIndex = boolIndex_location_seq_T(tempi,:);
            
            ticklabels = cell(numFrames,1);
            for temptempi = 1:numFrames
                
                if tempSeqBoolIndex(temptempi) == true
                    ticklabels{temptempi} = ['\color[rgb]{0.6350 0.0780 0.1840}' num2str(temptempi)];
                else
                    ticklabels{temptempi} = ['\color[rgb]{0.25 0.25 0.25}' num2str(temptempi)];
                end
            end
            set(gca,'XTickLabel',ticklabels);
            
            set(gca,'box','off');% 取消右、上边框
            if tempi == 1
                %xlabel('Location', 'FontSize', 9, 'FontWeight', 'bold');
                ylabel('probability', 'FontSize', 9, 'FontWeight', 'bold');
            end
            
            temp_p = p_n11n(tempi);
            temp_str = 'ns';
            if temp_p < 0.001
                temp_str = sprintf('***');
            elseif temp_p < 0.01
                temp_str = sprintf('**');
            elseif temp_p < 0.05
                temp_str = sprintf('*');
            end
            
            
            %temp_title = title(sprintf('%s, %s',num2str(boolIndex_location_seq(:,tempi)'),temp_str));
            temp_title = title(sprintf('%s',temp_str));
            temp_title.FontSize = 8;%10
            temp_title.Interpreter = 'none';
            
            
        end
        
        
    end
end


%% Fig, location distribution of each sequence, confusion matrix
if if_plot_eachSeq_confusionMatrixB == 1
    %% locDistri_confusionMatrix_A
    %     fig = figure('Name','locDistri_confusionMatrix_A','NumberTitle','off');
    %     set(gcf,'Position',[10 50 650 600]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %
    %     locDistri_confusionMatrix_A_T = locDistri_confusionMatrix_A';
    %
    %     C = locDistri_confusionMatrix_A_T;
    %
    %     C_max = max(C,[],'all') * 1.0;
    %
    %     imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
    %     hold on
    %
    %     axis equal
    %     set(gca, 'YDir','normal');%normal, reverse
    %     set(gca,'box','off');% 取消右、上边框
    %     set(gca, 'FontSize', 11) %12
    %     c = colorbar('FontSize',14);
    %     xlabel('Population decoder', 'FontSize', 14, 'FontWeight', 'bold');
    %     ylabel('Behavior model', 'FontSize', 14, 'FontWeight', 'bold');
    %     temp_title = title(sprintf('Location distribution similarity'),'FontSize',14,'FontWeight','bold');
    %     temp_title.Interpreter = 'none';
    
    
    
    %% locDistri_confusionMatrix_B
    fig = figure('Name','locDistri_confusionMatrix_B','NumberTitle','off');
    %set(gcf,'Position',[10 50 600*0.8 650*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 50 480 520*0.828]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 50 288*1.05 431]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 50 288*1.00 431]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    %t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Loose');
    t = tiledlayout(1,11,'TileSpacing','loose','Padding','loose');
    
    t.Title.String = sprintf('Location distribution probability');
    %t.Title.String = sprintf('      Location distribution probability');
    t.Title.FontSize = 12;
    t.Title.Interpreter = 'none';
    t.Title.FontWeight = 'bold';
    
    locDistri_confusionMatrix_B;
    
    %% locDistri_confusionMatrix_B(:,:,1)
    %nexttile
    nexttile([1 5])
    
    C = locDistri_confusionMatrix_B(:,:,1);
    C_max = 1;
    
    imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
    hold on
    
    set(gca,'linewidth',1.5)
    set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 12) %12
    
    % Plot length bound in response
    for temptempi=1:size(locDistri_confusionMatrix_B,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(locDistri_confusionMatrix_B,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    seqSet_inOne_inOne = [];
    for temptempi=1:length(seqSet_inOne)
        seqSet_inOne_inOne = [seqSet_inOne_inOne  seqSet_inOne{temptempi}];
    end
    
    set(gca,'YTick',1:size(locDistri_confusionMatrix_B,1));
    
    ytl=string(seqSet_inOne_inOne(1:size(locDistri_confusionMatrix_B,1)));
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
    text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',7);%8
    set(gca,'yticklabel','');
    
    set(gca,'XTick',[1:1:numFrames]);
    
    xtickangle(0);
    
    set(gca,'TickLength',[0 0]);
    
    xlabel('Location', 'FontSize', 12, 'FontWeight', 'bold');
    
    x_lim = xlim;
    y_lim = ylim;
    %ylabel('Sequence','Position',[min(x_lim)-0.8 mean(y_lim)],'FontSize',12,'FontWeight','bold');
    ylabel('Sequence','Position',[min(x_lim)-1.1 mean(y_lim)],'FontSize',12,'FontWeight','bold');
    
    
    if if_locationDistri_behavior0_model1 == 0
        temp_title = title(sprintf('Behavior'),'FontSize',11);
    elseif if_locationDistri_behavior0_model1 == 1
        temp_title = title(sprintf('Behavior model'),'FontSize',11);
    end
    temp_title.Interpreter = 'none';
    
    
    
    %% locDistri_confusionMatrix_B(:,:,2)
    %nexttile
    nexttile([1 5])
    
    C = locDistri_confusionMatrix_B(:,:,2);
    C_max = max(C,[],'all') * 1.0;
    
    imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
    hold on
    
    a = 1;
    
    
    set(gca,'linewidth',1.5)
    set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 12) %12
    %c = colorbar('FontSize',8);
    c = colorbar('westoutside','FontSize',8);
    %c.Layout.Tile = 'south';
    %c.Position(4) = c.Position(4)/2;
    c.Position = c.Position+[0+0.11+0.26+0.03 0 0-0.012 0];
    
    yticklabels('');
    
    set(gca,'TickLength',[0 0]);
    
    xtickangle(0);
    
    set(gca,'XTick',[1:1:numFrames]);
    xlabel('Location', 'FontSize', 12, 'FontWeight', 'bold');
    %temp_title = title(sprintf('Population decoder'),'FontSize',11);
    temp_title = title(sprintf('  Population decoder'),'FontSize',11);
    temp_title.Interpreter = 'none';
    
    
end

if if_plot_eachSeq_confusionMatrixC == 1
    %% locDistri_confusionMatrix_C
    fig = figure('Name','locDistri_confusionMatrix_C','NumberTitle','off');
    %set(gcf,'Position',[10 50 288 431]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 50 533 224]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 50 533*0.85*0.95 224*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    %t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Loose');
    t = tiledlayout(11,1,'TileSpacing','tight','Padding','Loose');
    
    t.Title.String = sprintf('Location distribution probability');
    t.Title.FontSize = 9;%12
    t.Title.Interpreter = 'none';
    t.Title.FontWeight = 'bold';
    
    locDistri_confusionMatrix_B;
    
    %% locDistri_confusionMatrix_C(:,:,1)
    %nexttile
    nexttile([5 1])
    
    C = locDistri_confusionMatrix_B(:,:,1)';
    C_max = 1;
    
    imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
    hold on
    
    set(gca,'linewidth',1.5)
    set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 8) %12
    
    
    set(gca,'xticklabel','');
    
    %yticklabels('');
    
    set(gca,'YTick',[1:1:numFrames]);
    
    xtickangle(0);
    
    set(gca,'TickLength',[0 0]);
    
    ylabel('Location', 'FontSize', 9, 'FontWeight', 'bold');%12
    
    
    if if_locationDistri_behavior0_model1 == 0
        temp_title = title(sprintf('Behavior'),'FontSize',9);%11
    elseif if_locationDistri_behavior0_model1 == 1
        temp_title = title(sprintf('Behavior model'),'FontSize',9);%11
    end
    temp_title.Interpreter = 'none';
    
    
    
    %% locDistri_confusionMatrix_C(:,:,2)
    %nexttile
    nexttile([5 1])
    
    C = locDistri_confusionMatrix_B(:,:,2)';
    %C_max = max(C,[],'all') * 1.0;
    C_max = 1;
    
    
    imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
    hold on
    
    for tempi=1:size(C,2)
        if ~isnan(C(1,tempi))
            continue
        end
        %plot(tempi*[1 1],[1-0.4,size(C,1)+0.4],...
        %    '-','LineWidth',6,'Color',[1 1 1]*0.5); %#ok<*UNRCH>
        %hold on
        %plot([tempi-0.3 tempi+0.3],[1-0.4,size(C,1)+0.4],...
        %    '-','LineWidth',1,'Color',[1 1 1]*0.25); %#ok<*UNRCH>
        %hold on
        
        plot([tempi-0.3 tempi+0.3],[1-0.4,size(C,1)+0.4],...
            '-','LineWidth',1,'Color',[1 1 1]*0.7); %#ok<*UNRCH>
        hold on
        
    end
    
    set(gca,'linewidth',1.5)
    set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 8) %12
    
    
    % Plot length bound in response
    for temptempi=1:size(locDistri_confusionMatrix_B,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(locDistri_confusionMatrix_B,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    seqSet_inOne_inOne = [];
    for temptempi=1:length(seqSet_inOne)
        seqSet_inOne_inOne = [seqSet_inOne_inOne  seqSet_inOne{temptempi}];
    end
    
    set(gca,'XTick',1:size(locDistri_confusionMatrix_B,1));
    
    xtl=string(seqSet_inOne_inOne(1:size(locDistri_confusionMatrix_B,1)));
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
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',90,'fontsize',6.5);%8
    set(gca,'xticklabel','');
    
    
    %yticklabels('');
    
    set(gca,'TickLength',[0 0]);
    
    xtickangle(0);
    
    set(gca,'YTick',[1:1:numFrames]);
    
    xlabel('Sequence', 'FontSize', 9, 'FontWeight', 'bold');%12
    ylabel('Location', 'FontSize', 9, 'FontWeight', 'bold');%12
    %ylabel(sprintf('Population decoder\nLocation'), 'FontSize', 10, 'FontWeight', 'bold');
    
    x_lim = xlim;
    y_lim = ylim;
    xlabel('Sequence','Position',[mean(x_lim) max(y_lim)+1.8],'FontSize',9,'FontWeight','bold');%12
    
    
    temp_title = title(sprintf('Population decoder'),'FontSize',9);%11
    temp_title.Interpreter = 'none';
    
    %c = colorbar('eastoutside','FontSize',8);
    %c.Position = c.Position+[0+0.11+0.26+0.03 0 0-0.012 0];
    c = colorbar('FontSize',7);
    %c = colorbar('east','FontSize',8);
    c.Layout.Tile = 'east';
    %c.Layout.Tile = 'eastoutside';
    c.Position = c.Position+[-0.05 0 -0.02 0];%[-0.15 0 0.0 0]
    
    if if_colormap_loadEnhanced == 1
        load('parula_enhanced');
        colormap(parula_enhanced);
    elseif if_colormap_loadEnhanced == 0
        colormap parula
    end
    
end



%% Fig, location distribution of each sequence, confusion matrix
if if_plot_eachSeq_confusionMatrixD == 1
    
    %% locDistri_confusionMatrix_D1
    fig = figure('Name','locDistri_confusionMatrix_D1','NumberTitle','off');
    %set(gcf,'Position',[10 50 137*0.95 393+3]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 50 137*0.95*1.3*1.1 396*0.9*1.02]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,7,'TileSpacing','loose','Padding','loose');
    
    
    locDistri_confusionMatrix_B;
    
    nexttile
    set(gca,'Visible','off');
    
    %% locDistri_confusionMatrix_D(:,:,1)
    nexttile([1 5])
    
    %C = locDistri_confusionMatrix_B(:,:,1);
    C = double(boolIndex_location_seq_T(1:sum(numSeq(1:3)),:));
    C_max = 1;
    
    imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
    hold on
    
    set(gca,'linewidth',1.5)
    set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 12) %12
    
    % Plot length bound in response
    for temptempi=1:size(locDistri_confusionMatrix_B,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(locDistri_confusionMatrix_B,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    seqSet_inOne_inOne = [];
    for temptempi=1:length(seqSet_inOne)
        seqSet_inOne_inOne = [seqSet_inOne_inOne  seqSet_inOne{temptempi}];
    end
    
    set(gca,'YTick',1:size(locDistri_confusionMatrix_B,1));
    
    ytl=string(seqSet_inOne_inOne(1:size(locDistri_confusionMatrix_B,1)));
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
    
    c = colorbar('westoutside','FontSize',7);
    %c.Position = c.Position+[0.4+0.24 0.01 -0.032 -0.03];
    %c.Position = c.Position+[0.4+0.22 0.01 -0.072 -0.03];
    %c.Position = c.Position+[0.4+0.23 0.01 -0.077 -0.03];
    %c.Position = c.Position+[0.4+0.23 0.00 -0.05 -0.5];
    c.Position = c.Position+[0.4+0.24 0.00 -0.04 -0.5];
    
    %c.Ticks = [0:0.1:1];
    c.Ticks = [0 1];
    
    if if_colormap_loadEnhanced == 1
        load('parula_enhanced');
        colormap(parula_enhanced);
    elseif if_colormap_loadEnhanced == 0
        colormap parula
    end
    
    x_lim = xlim;
    y_lim = ylim;
    
    %xlabel('Location', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Location','Position',[mean(x_lim) max(y_lim)+1.4],'FontSize', 10, 'FontWeight', 'bold');
    
    %ylabel('Sequence','FontSize',12,'FontWeight','bold');
    ylabel('Sequence','Position',[min(x_lim)-1.15 mean(y_lim)],'FontSize',10,'FontWeight','bold');%-1.5
    
    
    temp_title = title(sprintf('Theoretical'),'FontSize',9);%11
    temp_title.Interpreter = 'none';
    
    
    %% locDistri_confusionMatrix_D1
    fig = figure('Name','locDistri_confusionMatrix_D1','NumberTitle','off');
    %set(gcf,'Position',[10 50 137*0.95 393+3]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 50 137*0.95*1.3*1.1 396*0.9*1.02]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,7,'TileSpacing','loose','Padding','loose');
    
    
    locDistri_confusionMatrix_B;
    
    nexttile
    set(gca,'Visible','off');
    
    %% locDistri_confusionMatrix_D(:,:,1)
    nexttile([1 5])
    
    C = locDistri_confusionMatrix_B(:,:,1);
    C_max = 1;
    
    imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
    hold on
    
    set(gca,'linewidth',1.5)
    set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 12) %12
    
    % Plot length bound in response
    for temptempi=1:size(locDistri_confusionMatrix_B,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(locDistri_confusionMatrix_B,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    seqSet_inOne_inOne = [];
    for temptempi=1:length(seqSet_inOne)
        seqSet_inOne_inOne = [seqSet_inOne_inOne  seqSet_inOne{temptempi}];
    end
    
    set(gca,'YTick',1:size(locDistri_confusionMatrix_B,1));
    
    ytl=string(seqSet_inOne_inOne(1:size(locDistri_confusionMatrix_B,1)));
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
    
    c = colorbar('westoutside','FontSize',7);
    %c.Position = c.Position+[0.4+0.24 0.01 -0.032 -0.03];
    %c.Position = c.Position+[0.4+0.22 0.01 -0.072 -0.03];
    %c.Position = c.Position+[0.4+0.23 0.01 -0.077 -0.03];
    %c.Position = c.Position+[0.4+0.23 0.00 -0.05 -0.5];
    c.Position = c.Position+[0.4+0.24 0.00 -0.04 -0.5];
    
    %c.Ticks = [0:0.1:1];
    c.Ticks = [0 1];
    
    if if_colormap_loadEnhanced == 1
        load('parula_enhanced');
        colormap(parula_enhanced);
    elseif if_colormap_loadEnhanced == 0
        colormap parula
    end
    
    x_lim = xlim;
    y_lim = ylim;
    
    %xlabel('Location', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Location','Position',[mean(x_lim) max(y_lim)+1.4],'FontSize', 10, 'FontWeight', 'bold');
    
    %ylabel('Sequence','FontSize',12,'FontWeight','bold');
    ylabel('Sequence','Position',[min(x_lim)-1.15 mean(y_lim)],'FontSize',10,'FontWeight','bold');
    
    
    if if_locationDistri_behavior0_model1 == 0
        temp_title = title(sprintf('Behavior'),'FontSize',9);%11
    elseif if_locationDistri_behavior0_model1 == 1
        temp_title = title(sprintf('Behavior model'),'FontSize',9);%11
    end
    temp_title.Interpreter = 'none';
    
    
    
    %% locDistri_confusionMatrix_D2
    fig = figure('Name','locDistri_confusionMatrix_D2','NumberTitle','off');
    %set(gcf,'Position',[160 50 137*0.95 393+3]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 50 137*0.95*1.3*1.1 396*0.9*1.02]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,7,'TileSpacing','loose','Padding','loose');
    
    
    locDistri_confusionMatrix_B;
    
    nexttile
    set(gca,'Visible','off');
    
    %% locDistri_confusionMatrix_D(:,:,1)
    nexttile([1 5])
    
    C = locDistri_confusionMatrix_B(:,:,2);
    C_max = 1;
    
    imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
    hold on
    
    for tempi=1:size(C,1)
        if ~isnan(C(tempi,1))
            continue
        end
        %plot([1-0.4,size(C,2)+0.4],tempi*[1 1],...
        %    '-','LineWidth',6,'Color',[1 1 1]*0.5); %#ok<*UNRCH>
        %hold on
        %plot([1-0.4,size(C,2)+0.4],[tempi-0.3 tempi+0.3],...
        %    '-','LineWidth',1,'Color',[1 1 1]*0.25); %#ok<*UNRCH>
        %hold on
        
        plot([1-0.4,size(C,2)+0.4],[tempi-0.3 tempi+0.3],...
            '-','LineWidth',1,'Color',[1 1 1]*0.7); %#ok<*UNRCH>
        hold on
    end
    
    
    set(gca,'linewidth',1.5)
    set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 12) %12
    
    % Plot length bound in response
    for temptempi=1:size(locDistri_confusionMatrix_B,1)
        if temptempi > 1
            if length(target_seqSet_inOne{temptempi}) > ...
                    length(target_seqSet_inOne{temptempi-1})
                plot([temptempi temptempi]-0.5, [0.5 size(locDistri_confusionMatrix_B,1)+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    seqSet_inOne_inOne = [];
    for temptempi=1:length(seqSet_inOne)
        seqSet_inOne_inOne = [seqSet_inOne_inOne  seqSet_inOne{temptempi}];
    end
    
    set(gca,'YTick',1:size(locDistri_confusionMatrix_B,1));
    
    ytl=string(seqSet_inOne_inOne(1:size(locDistri_confusionMatrix_B,1)));
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
    
    %c = colorbar('westoutside','FontSize',7);
    %c.Position = c.Position+[0.4+0.24 0.01 -0.032 -0.03];
    %c.Position = c.Position+[0.4+0.22 0.01 -0.072 -0.03];
    %c.Position = c.Position+[0.4+0.23 0.01 -0.077 -0.03];
    
    if if_colormap_loadEnhanced == 1
        load('parula_enhanced');
        colormap(parula_enhanced);
    elseif if_colormap_loadEnhanced == 0
        colormap parula
    end
    
    
    x_lim = xlim;
    y_lim = ylim;
    
    %xlabel('Location', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Location','Position',[mean(x_lim) max(y_lim)+1.4],'FontSize', 10, 'FontWeight', 'bold');
    
    %ylabel('Sequence','FontSize',12,'FontWeight','bold');
    ylabel('Sequence','Position',[min(x_lim)-1.15 mean(y_lim)],'FontSize',10,'FontWeight','bold');
    
    
    temp_title = title(sprintf(' Population decoder'),'FontSize',9);
    temp_title.Interpreter = 'none';
    
    
    %     %% locDistri_confusionMatrix_D(:,:,2)
    %     nexttile([1 5])
    %
    %     C = locDistri_confusionMatrix_B(:,:,2);
    %     C_max = max(C,[],'all') * 1.0;
    %
    %     imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
    %     hold on
    %
    %     set(gca,'linewidth',1.5)
    %     set(gca,'box','off');% 取消右、上边框
    %     set(gca, 'FontSize', 12) %12
    %
    %     yticklabels('');
    %
    %     set(gca,'TickLength',[0 0]);
    %
    %     xtickangle(0);
    %
    %     set(gca,'XTick',[1:1:numFrames]);
    %
    %     c = colorbar('westoutside','FontSize',8);
    %     c.Position = c.Position+[0+0.11+0.26+0.03 0 0-0.012 0];
    %
    %     xlabel('Location', 'FontSize', 12, 'FontWeight', 'bold');
    %
    %     x_lim = xlim;
    %     y_lim = ylim;
    %     %ylabel('Sequence','Position',[min(x_lim)-1.1 mean(y_lim)],'FontSize',12,'FontWeight','bold');
    %     ylabel('Sequence','FontSize',12,'FontWeight','bold');
    %
    %     %temp_title = title(sprintf('Population decoder'),'FontSize',11);
    %         temp_title = title(sprintf('  Population decoder'),'FontSize',11);
    %     temp_title.Interpreter = 'none';
    
    
end



%% Fig, location distribution of each sequence, confusion matrix
if if_plot_eachSeq_confusionMatrixE == 1
    
    seqSet_inOne_inOne = [];
    for temptempi=1:length(seqSet_inOne)
        seqSet_inOne_inOne = [seqSet_inOne_inOne  seqSet_inOne{temptempi}];
    end
    locDistri_confusionMatrix_B;
    
    
    %% locDistri_confusionMatrix_E1 (theoretical)
    if true
        fig = figure('Name','locDistri_confusionMatrix_E1','NumberTitle','off');
        
        set(gcf,'Position',[10 50 390 103]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(3,31,'TileSpacing','none','Padding','compact');
        
        nexttile([2 30])
        
        C = double(boolIndex_location_seq(:,1:sum(numSeq(1:3))));
        C_max = 1;
        
        imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
        hold on
        
        set(gca,'linewidth',1.5)
        set(gca,'box','off');% 取消右、上边框
        set(gca, 'FontSize', 8) %12
        
        
        set(gca,'XTick',1:size(C,2));
        
        xtl=string(seqSet_inOne_inOne(1:size(C,2)));
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
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',90,'fontsize',6.5);%8
        set(gca,'xticklabel','');
        xtickangle(0);
        
        set(gca,'YTick',[1:1:numFrames]);
        set(gca,'TickLength',[0 0]);
        
        x_lim = xlim;
        y_lim = ylim;
        xlabel('Sequence','Position',[mean(x_lim) max(y_lim)+1.35],'FontSize',9,'FontWeight','bold');%12
        ylabel('Location', 'FontSize', 9, 'FontWeight', 'bold');%12
        temp_title = title(sprintf('Theoretical'),'FontSize',9);%11
        temp_title.Interpreter = 'none';
        
        c = colorbar('westoutside','FontSize',7);
        c.Position = c.Position+[0.847 -0.01 -0.02 -0.00];
        c.Ticks = [0 1];
        
        if if_colormap_loadEnhanced == 1
            %load('parula_enhanced');
            %colormap(parula_enhanced);
            %load('gray_enhanced');
            %colormap(gray_enhanced);
            
            %load('gray_reversed_enhanced');
            %colormap(gray_reversed_enhanced);
            
            load('coolwarm_enhanced');
            colormap(coolwarm_enhanced);
            
            %load('coolwarm_n11n_enhanced');
            %colormap(coolwarm_n11n_enhanced);
            
        elseif if_colormap_loadEnhanced == 0
            %colormap parula
            %colormap gray
            
            %temp1 = gray;
            %temp1 = temp1(end:-1:1,:);
            %colormap(temp1);
            
            colormap(coolwarm());
            
            %temp1 = coolwarm(300);
            %temp2 = ((300-256)/2)+1;
            %temp3 = temp1(temp2:temp2+255,:);
            %colormap(temp3);
        end
        
        
    end
    
    
    %% locDistri_confusionMatrix_E2 (population decoder)
    if true
        fig = figure('Name','locDistri_confusionMatrix_E2','NumberTitle','off');
        
        set(gcf,'Position',[10 50 390 103]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(3,31,'TileSpacing','none','Padding','compact');
        
        nexttile([2 30])
        
        C = locDistri_confusionMatrix_B(:,:,2)';
        C_max = 1;
        
        imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
        hold on
        
        for tempi=1:size(C,2)
            if ~isnan(C(1,tempi))
                continue
            end
            plot([tempi-0.3 tempi+0.3],[1-0.4,size(C,1)+0.4],...
                '-','LineWidth',1,'Color',[1 1 1]*0.7); %#ok<*UNRCH>
            hold on
            
        end
        
        set(gca,'linewidth',1.5)
        set(gca,'box','off');% 取消右、上边框
        set(gca, 'FontSize', 8) %12
        
        
        set(gca,'XTick',1:size(C,2));
        
        xtl=string(seqSet_inOne_inOne(1:size(C,2)));
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
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',90,'fontsize',6.5);%8
        set(gca,'xticklabel','');
        xtickangle(0);
        
        set(gca,'YTick',[1:1:numFrames]);
        set(gca,'TickLength',[0 0]);
        
        x_lim = xlim;
        y_lim = ylim;
        xlabel('Sequence','Position',[mean(x_lim) max(y_lim)+1.3],'FontSize',9,'FontWeight','bold');%12
        ylabel('Location', 'FontSize', 9, 'FontWeight', 'bold');%12
        temp_title = title(sprintf('Population decoder'),'FontSize',9);%11
        temp_title.Interpreter = 'none';
        
        if if_colormap_loadEnhanced == 1
            %load('parula_enhanced');
            %colormap(parula_enhanced);
            %load('gray_enhanced');
            %colormap(gray_enhanced);
            
            %load('gray_reversed_enhanced');
            %colormap(gray_reversed_enhanced);
            
            load('coolwarm_enhanced');
            colormap(coolwarm_enhanced);
            
            %load('coolwarm_n11n_enhanced');
            %colormap(coolwarm_n11n_enhanced);
            
        elseif if_colormap_loadEnhanced == 0
            %colormap parula
            %colormap gray
            
            %temp1 = gray;
            %temp1 = temp1(end:-1:1,:);
            %colormap(temp1);
            
            colormap(coolwarm());
            
            %temp1 = coolwarm(300);
            %temp2 = ((300-256)/2)+1;
            %temp3 = temp1(temp2:temp2+255,:);
            %colormap(temp3);
        end
        
    end
    
    
    
    
    %% locDistri_confusionMatrix_E3 (behavior)
    if true
        fig = figure('Name','locDistri_confusionMatrix_E3','NumberTitle','off');
        
        set(gcf,'Position',[10 50 390 103]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(3,31,'TileSpacing','none','Padding','compact');
        
        nexttile([2 30])
        
        C = locDistri_confusionMatrix_B(:,:,1)';
        C_max = 1;
        
        imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
        hold on
        
        set(gca,'linewidth',1.5)
        set(gca,'box','off');% 取消右、上边框
        set(gca, 'FontSize', 8) %12
        
        
        set(gca,'XTick',1:size(C,2));
        
        xtl=string(seqSet_inOne_inOne(1:size(C,2)));
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
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',90,'fontsize',6.5);%8
        set(gca,'xticklabel','');
        xtickangle(0);
        
        set(gca,'YTick',[1:1:numFrames]);
        set(gca,'TickLength',[0 0]);
        
        x_lim = xlim;
        y_lim = ylim;
        xlabel('Sequence','Position',[mean(x_lim) max(y_lim)+1.3],'FontSize',9,'FontWeight','bold');%12
        ylabel('Location', 'FontSize', 9, 'FontWeight', 'bold');%12
        if if_locationDistri_behavior0_model1 == 0
            temp_title = title(sprintf('Behavior'),'FontSize',9);%11
        elseif if_locationDistri_behavior0_model1 == 1
            temp_title = title(sprintf('Behavior model'),'FontSize',9);%11
        end
        temp_title.Interpreter = 'none';
        
        if if_colormap_loadEnhanced == 1
            %load('parula_enhanced');
            %colormap(parula_enhanced);
            %load('gray_enhanced');
            %colormap(gray_enhanced);
            
            %load('gray_reversed_enhanced');
            %colormap(gray_reversed_enhanced);
            
            load('coolwarm_enhanced');
            colormap(coolwarm_enhanced);
            
            %load('coolwarm_n11n_enhanced');
            %colormap(coolwarm_n11n_enhanced);
            
        elseif if_colormap_loadEnhanced == 0
            %colormap parula
            %colormap gray
            
            %temp1 = gray;
            %temp1 = temp1(end:-1:1,:);
            %colormap(temp1);
            
            colormap(coolwarm());
            
            %temp1 = coolwarm(300);
            %temp2 = ((300-256)/2)+1;
            %temp3 = temp1(temp2:temp2+255,:);
            %colormap(temp3);
        end
        
    end
    
end


%% Location correlation \n Correct trials
temp_1 = r_n11n_resampleMean_correct;
temp_1_chanceLevel = mean(r_n11n_resampleMean_correct_chanceLevel);

[~,temp_p] = ttest(temp_1,temp_1_chanceLevel);

if false
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.9*0.85 336*1.11*0.5*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.9*0.85*0.9*0.9 336*1.11*0.5*0.9*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    %temp_y_min = min([temp_1 temp_1_chanceLevel]);
    %temp_y_max = max(temp_1);
    
    if if_dff_singleSession1_twoSession2_allSession3 == 3
        temp_y_min = -0.4;
        temp_y_max = 1;
    end
    
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
    set(gca, 'FontSize', 8) %14,12
    %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
    set(gca,'YTick',0:1,'FontSize', 9);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""]; %#ok<*NBRAK>
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    
    set(gca,'xticklabel','');
    
    
    %ylabel('Correlation', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
    %temp_title = title(sprintf('Location correlation \n Correct trials'),'fontsize',11);
    temp_title = title(sprintf('Location-level'),'fontsize',9);
    temp_title.Interpreter = 'none';
    
end


%% Location correlation Correct trials
temp_1 = r_n11n_resampleMean_correct_theoretical;
temp_1_chanceLevel = mean(r_n11n_resampleMean_correct_theoretical_chanceLevel);

% [~,temp_p] = ttest(temp_1,temp_1_chanceLevel);
[~,temp_p] = ttest2(temp_1,r_n11n_resampleMean_correct_theoretical_chanceLevel);
% temp_p

if false
    fig = figure('Name','','NumberTitle','off');
    set(gcf,'Position',[50+0 50+0 240*0.80*1.1 336*1.11*0.5*0.7]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
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
    set(gca, 'FontSize', 9) %14
    set(gca,'YTick',0:0.2:1,'FontSize', 9);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""]; %#ok<*NBRAK>
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Correlation', 'FontSize', 10, 'FontWeight', 'normal');
    %temp_title = title(sprintf('Location correlation\ntheoretical'),'fontsize',10);
    %temp_title = title(sprintf('Location-level\ntheoretical VS. correct'),'fontsize',10);
    temp_title = title(sprintf('Location-level\ntheoretical & decoder'),'fontsize',10);
    temp_title.Interpreter = 'none';
    
end


if false
    fig = figure('Name','','NumberTitle','off');
    set(gcf,'Position',[50+0 50+0 240*0.80*1.1 336*1.11*0.5*0.7*0.85*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    
    x1 = temp_1;
    
    
    h_NumBins = 8;%10
    
    x = x1;
    h1 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
    hold on
    h1.NumBins = h_NumBins;
    h1.EdgeColor = [0 0 0];
    
    y1 = h1.Values;
    
    x_min = 0;
    x_max = 1;
    
    [y_min,y_max] = bounds(y1);
    
    plot(temp_1_chanceLevel*[1 1],[y_min y_max],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(0.2,y_min+(y_max-y_min)*0.7,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    
    
    set(gca,'linewidth',1.5)
    xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
    ylim([y_min y_max+(y_max-y_min)*0.1]);%0.1
    xticks([0 1]);
    set(gca, 'FontSize', 9)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Correlation', 'FontSize', 10, 'FontWeight', 'normal','Position',[0.5 y_min-(y_max-y_min)*0.095]);
    ylabel('Seq count', 'FontSize', 10, 'FontWeight', 'normal');
    
    temp_title = title(sprintf('Location-level\ntheoretical & decoder'),'fontsize',8.5);
    temp_title.Interpreter = 'none';
    
    
    %     temp_y_min = min([temp_1 temp_1_chanceLevel]);
    %     temp_y_max = max(temp_1);
    %
    %
    %     temp_data = temp_1';
    %     temp_label = repmat({'A'},length(temp_1),1);
    %
    %     temptemp_color1 = [1 1 1]*0.5;
    %
    %     h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
    %         'GroupOrder',{'A'},'Width',0.8);
    %     h(1).ViolinPlot.FaceAlpha = 0.1;
    %     hold on
    %
    %
    %     plot([0 2],temp_1_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    %     hold on
    %
    %
    %     tempTxt = sprintf('');
    %     if temp_p < 0.001
    %         tempTxt = sprintf('***');
    %     elseif temp_p < 0.01
    %         tempTxt = sprintf('**');
    %     elseif temp_p < 0.05
    %         tempTxt = sprintf('*');
    %     end
    %     text(1.8,temp_y_min+(temp_y_max-temp_y_min)*0.3,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %         'HorizontalAlignment','center');
    %
    %     set(gca,'linewidth',1.5)
    %
    %     xlim([0 2]);
    %     ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    %     set(gca, 'FontSize', 9) %14
    %     %set(gca,'YTick',0:0.2:1,'FontSize', 9);%12
    %     set(gca,'YTick',0:1,'FontSize', 9);%12
    %     set(gca,'box','off');% 取消右、上边框
    %
    %
    %     xtl = [""]; %#ok<*NBRAK>
    %     xt=get(gca,'XTick');
    %     yt=get(gca,'YTick');
    %     xtext_xp=xt;
    %     xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    %     text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    %
    %     set(gca,'xticklabel','');
    %
    %     view([90 -90]);
    %
    %     ylabel('Correlation', 'FontSize', 10, 'FontWeight', 'normal');
    %     %temp_title = title(sprintf('Location correlation\ntheoretical'),'fontsize',10);
    %     %temp_title = title(sprintf('Location-level\ntheoretical VS. correct'),'fontsize',10);
    %     temp_title = title(sprintf('Location-level\ntheoretical & decoder'),'fontsize',10);
    %     temp_title.Interpreter = 'none';
    
end




if true
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[50+0 250+0 240*0.80*1.1*0.7*1.1 336*1.11*0.5*0.7*0.85*0.95*1.15*1.03]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[450+0 50+0 240*0.80*1.1*0.7*1.1 336*1.11*0.5*0.7*0.85*0.95*1.15*1.03*1.17]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[650+0 50+0 158.4*0.94*1.03 167.5*0.94*0.95*1.2*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    
    x1 = temp_1;
    
    
    h_NumBins = 8;%10
    
    x = x1;
    h1 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
    hold on
    h1.NumBins = h_NumBins;
    h1.EdgeColor = [0 0 0];
    
    y1 = h1.Values;
    
    x_min = 0;
    x_max = 1;
    
    [y_min,y_max] = bounds(y1);
    y_max = max(y_max,10);
    
    plot(temp_1_chanceLevel*[1 1],[y_min y_max],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    %text(0.2,y_max-(y_max-y_min)*0.03,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %    'HorizontalAlignment','center');
    text(0.2,y_max+(y_max-y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    
    
    set(gca,'linewidth',1.5)
    xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
    ylim([y_min y_max+(y_max-y_min)*0.15]);%0.1
    xticks([0 1]);
    yticks([y_min y_max]);
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    %xlabel('Correlation', 'FontSize', 10, 'FontWeight', 'normal','Position',[0.5 y_min-(y_max-y_min)*0.095]);
    xlabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
    ylabel('Seq number', 'FontSize', 9, 'FontWeight', 'normal');
    
    temp_title = title(sprintf('Location-level'),'fontsize',8.5);
    temp_title.Interpreter = 'none';
    temp_subtitle = subtitle(sprintf('theoretical & decoder'),'fontsize',7.5);
    
end

%% Location-level, Behavior model & decoder
temp_1 = r_n11n_resampleMean_correct;
temp_1_chanceLevel = mean(r_n11n_resampleMean_correct_chanceLevel);

[~,temp_p] = ttest(temp_1,temp_1_chanceLevel);

if true
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[650+0 50+0 158.4*0.94 167.5*0.94*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[650+0 50+0 158.4*0.94*1.03 167.5*0.94*0.95*1.2*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    
    x1 = temp_1;
    
    
    h_NumBins = 8;%10
    
    x = x1;
    h1 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
    hold on
    h1.NumBins = h_NumBins;
    h1.EdgeColor = [0 0 0];
    
    y1 = h1.Values;
    
    x_min = 0;
    x_max = 1;
    
    %if if_dff_singleSession1_twoSession2_allSession3 == 3
    if if_monkey_D0_Z1 == 1
        x_min = min(x);
    end
    %end
    
    [y_min,y_max] = bounds(y1); %#ok<*ASGLU>
    
    y_min = 0;
    
    plot(temp_1_chanceLevel*[1 1],[y_min y_max],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    %text(0.2,y_min+(y_max-y_min)*0.7,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %    'HorizontalAlignment','center');
    text(0.2,y_max+(y_max-y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    
    
    set(gca,'linewidth',1.5)
    xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
    %xlim([-0.1 1.1]);
    ylim([y_min y_max+(y_max-y_min)*0.15]);%0.1
    xticks([0 1]);
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    %xlabel('Correlation', 'FontSize', 10, 'FontWeight', 'normal','Position',[0.5 y_min-(y_max-y_min)*0.095]);
    xlabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
    ylabel('Seq number', 'FontSize', 9, 'FontWeight', 'normal');
    
    %temp_title = title(sprintf('Location-level\nBehavior model & decoder'),'fontsize',8.5);
    temp_title = title(sprintf('Location-level'),'fontsize',9);
    temp_title.Interpreter = 'none';
    %temp_subtitle = subtitle(sprintf('Behavior model & decoder'),'fontsize',7.5);
    if if_locationDistri_behavior0_model1 == 0
        %temp_subtitle = subtitle(sprintf('decoder & behavior'),'fontsize',7.5);
        temp_subtitle = subtitle(sprintf('behavior & decoder'),'fontsize',7.5);
    elseif if_locationDistri_behavior0_model1 == 1
        %temp_subtitle = subtitle(sprintf('decoder & behavior model'),'fontsize',7.5);
        temp_subtitle = subtitle(sprintf('behavior model & decoder'),'fontsize',7.5);
    end
    
    
end


%% Compare r_n11n_resampleMean_correct_theoretical & r_n11n_resampleMean_correct
if true
    [~,temp_p1] = ttest(r_n11n_resampleMean_correct_theoretical,r_n11n_resampleMean_correct);
    [~,temp_p2] = ttest2(r_n11n_resampleMean_correct_theoretical,r_n11n_resampleMean_correct);
    
    if_temp_violinplot0_pairline1 = 1;
    
    [~,temp_p] = ttest(r_n11n_resampleMean_correct_theoretical,r_n11n_resampleMean_correct);
    
    fig = figure('Name','theoretical & correct','NumberTitle','off');
    
    set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35 295*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    
    
    % Compare memory precision of error and correct trials
    nexttile
    
    temp_p = temp_p;
    temptempBoolIndex = (~isnan(r_n11n_resampleMean_correct_theoretical)) & (~isnan(r_n11n_resampleMean_correct));
    temp_1 = r_n11n_resampleMean_correct_theoretical(temptempBoolIndex)';
    temp_2 = r_n11n_resampleMean_correct(temptempBoolIndex)';
    
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    
    temp_y_min = min([temp_1;temp_2]);
    temp_y_max = max([temp_1;temp_2]);
    
    
    if if_temp_violinplot0_pairline1 == 0
        temp_data = [temp_1;temp_2];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        
        temp_label = [g1;g2];
        
        temptemp_color1 = [1 1 1]*0.5;
        temptemp_color2 = repmat(temptemp_color1, 2, 1);
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        
    elseif if_temp_violinplot0_pairline1 == 1
        
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);%30
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);%30
        hold on
        
    end
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.12,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 2.85])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.20]);
    set(gca, 'FontSize', 8);%10
    
    xticks([1 2]);
    
    xtl = ["Theoretical"; "Behavior"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.325;%0.56,0.4
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.22;%0.56,0.4
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9
    set(gca,'xticklabel','');
    
    
    yticks([0 1]);
    
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');%10
    
    title(sprintf('Location correlation'),'FontSize',9,'FontWeight','normal');
    subtitle(sprintf('With decoder'),'FontSize',7,'FontWeight','normal');
    
    
end


%% Location accuracy
Posterior_2d_model;
locAcc_mean = mean(Posterior_2d_model,1);
locAcc_sem = std(Posterior_2d_model,0,1)./sqrt(size(Posterior_2d_model,1));

if true
    fig = figure('Name','locAcc','NumberTitle','off');
    set(gcf,'Position',[0 50 200 200*0.7]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');

    nexttile  

    temp_y_min = min(locAcc_mean-locAcc_sem);
    temp_y_max = max(locAcc_mean+locAcc_sem);
    
    %temp_y_min = 0.17;%0.2,0.25,0.21
    %temp_y_max = 1-temp_y_min;%0.6,0.55,0.51,0.79
    
    temp_y_min = 0.17;%0.25,0.21,0.17
    temp_y_max = 1-temp_y_min;%
    
    plot(1:6,locAcc_mean,'Color',[0.3 0.3 0.3],'LineWidth',1);
    hold on

    errorbar(1:6,locAcc_mean,locAcc_sem,'-o','Color',[0.3 0.3 0.3],'LineWidth',1,'CapSize',6,'MarkerSize',1);
    hold on
    

    set(gca,'linewidth',1.5)
    set(gca,'box','off');% 取消右、上边框
    set(gca, 'FontSize', 8)

    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.5 temp_y_max+(temp_y_max-temp_y_min)*0.5]);
    set(gca,'XLim',[1-0.75 6+0.75]);%X轴的数据显示范围
    set(gca,'XTick',[1:1:6]);%设置要显示坐标刻度的范围
    
    xlabel('Location', 'FontSize', 9);
    ylabel('Recall accuracy', 'FontSize', 9);
    
end

%% End