%This script: To conduct behavioral analyses and plot them, using input of 'A_kMeans_acc_monkey***.m' and 'A_rProbAccRoughAnalysis_monkey***.m'
%% Initialization
close all

targetPATH = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis';
cd(targetPATH)

targetPATH2 = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';

if_title = 1;

if_plot_precision = 1;
if_entropy = 1;
if_entropy_behaviorInverted = 0;

if_singlePanel0_multiPanel1 = 1;

if_precision_meanProb0_sumProb1 = 0;%0

if_plot_rProb_scale = 1;
if_plot_rProb_fit = 0;

if_plot_multiPanelB = 1;
if_plot_multiPanelC = 0;

if_plot_submissionRT = 0;

if_monkey_D0_Z1 = 0;% To decide whether dealing with Ding's data or Zelku's data
if_compute = 1;

% if if_monkey_D0_Z1 == 1
%     if_title = 0;
% end

% spatialInterferenceWeight = 0.5;

exampleSeq = 17;%25,32,14,12

% figure_xPixel = 173*1.2*0.9*0.8;%173
% figure_yPixel = 147*0.8;%147

figure_xPixel = 150*1.2;
figure_yPixel = 118*1.2;


% color_choiceMemory = [102,194,165]/255;%[0.4660 0.6740 0.1880], [179,205,227]/255, [141,160,203]/255
% color_choiceOffload = [252,141,98]/255;%[0.6350 0.0780 0.1840], [251,180,174]/255


color_choiceMemory = [1,133+20,113]/255;
color_choiceOffload = [166+40,97,26]/255;


color_noChoice = [0.3 0.3 0.3];%[0.25 0.25 0.25]
color_rProb = [0.5 0.5 0.5];


if if_compute == 1
    A_kMeans_acc_monkey_isCorrectBugDone_Name_v = autoGetFunName_myScripts('A_kMeans_acc_monkey_isCorrectBugDone', [targetPATH]);
    A_kMeans_acc_monkey= str2func(A_kMeans_acc_monkey_isCorrectBugDone_Name_v);
    A_kMeans_acc_monkey();
    close all
    A_rProbAccRoughAnalysis_monkey_Name_v = autoGetFunName_myScripts('A_rProbAccRoughAnalysis_monkey', [targetPATH]);
    A_rProbAccRoughAnalysis_monkey= str2func(A_rProbAccRoughAnalysis_monkey_Name_v);
    A_rProbAccRoughAnalysis_monkey();
    close all
end




valid_length = 4;
if if_plot_precision == 1
    % if if_plot_precision == 1 && if_compute == 1
    valid_length = 4;
    
    boolIndex_location_seq = false(numFrames,length(target_seqSet_inOne));
    for tempi=1:length(target_seqSet_inOne)
        currentSequence = target_seqSet_inOne{tempi};
        boolIndex_location_seq(currentSequence,tempi) = true;
    end
    boolIndex_location_seq_T = boolIndex_location_seq';
    
    responseMatrix_noChoice;
    responseMatrix_noChoice_valid = responseMatrix_noChoice(1:sum(numSeq(1:valid_length)),1:sum(numSeq(1:valid_length)));
    responseMatrix_noChoice_valid_n11n = responseMatrix_noChoice_valid ./ sum(responseMatrix_noChoice_valid,2);
    
    if if_entropy == 0
        
        %% To get optimal spatialInterferenceWeight
        fun_seqSimilarity_Name_v = autoGetFunName_myScripts('fun_seqSimilarity', [targetPATH2 '\functions']);
        fun_seqSimilarity = str2func(fun_seqSimilarity_Name_v);
        
        temp_spatialInterferenceWeight = 0:0.005:1;
        temp_r_mean_iter = nan(length(temp_spatialInterferenceWeight),1);
        for tempIter=1:length(temp_spatialInterferenceWeight)
            temp_score_stimuli_to_response = nan(sum(numSeq(1:4)),sum(numSeq(1:4)));
            for tempi=1:sum(numSeq(1:4))
                boolIndex_location_A = boolIndex_location_seq_T(tempi,:);
                
                for tempj=1:sum(numSeq(1:4))
                    boolIndex_location_B = boolIndex_location_seq_T(tempj,:);
                    
                    temp_score_stimuli_to_response(tempi,tempj) = ...
                        fun_seqSimilarity(boolIndex_location_A,boolIndex_location_B,temp_spatialInterferenceWeight(tempIter));
                    
                    if tempi == 24 && tempj == 23
                        a = 1; %#ok<*NASGU>
                    end
                end
            end
            
            temp_r = nan(size(responseMatrix_noChoice_valid_n11n,1),1);
            for tempi=1:length(temp_r)
                [temp_r(tempi),~] = corr(responseMatrix_noChoice_valid_n11n(tempi,:)',temp_score_stimuli_to_response(tempi,:)');
            end
            temp_r_mean_iter(tempIter) = mean(temp_r);
        end
        [M,I] = max(temp_r_mean_iter);
        temp_r_max_iter = M;
        spatialInterferenceWeight = temp_spatialInterferenceWeight(I);
        fprintf('spatialInterferenceWeight=%.4f,r=%.4f.\n',spatialInterferenceWeight,temp_r_max_iter);
        
        %% score_stimuli_to_response
        score_stimuli_to_response = nan(sum(numSeq(1:4)),sum(numSeq(1:4)));
        for tempi=1:sum(numSeq(1:4))
            boolIndex_location_A = boolIndex_location_seq_T(tempi,:);
            for tempj=1:sum(numSeq(1:4))
                %if tempi == 24 && tempj == 23
                if tempi == 1 && tempj == 6
                    a = 1; %#ok<*NASGU>
                end
                boolIndex_location_B = boolIndex_location_seq_T(tempj,:);
                score_stimuli_to_response(tempi,tempj) = ...
                    fun_seqSimilarity(boolIndex_location_A,boolIndex_location_B,spatialInterferenceWeight);
                
            end
        end
        
        
        %% seqPrecision_behavior
        score_stimuli_to_response;
        
        %responseMatrix_noChoice;
        %responseMatrix_noChoice_valid = responseMatrix_noChoice(1:sum(numSeq(1:valid_length)),1:sum(numSeq(1:valid_length)));
        %responseMatrix_noChoice_valid_n11n = responseMatrix_noChoice_valid ./ sum(responseMatrix_noChoice_valid,2);
        
        
        % only for test, Precision that accuracy is controlled
        %responseMatrix_noChoice_valid_n11n;
        %for tempi=1:size(responseMatrix_noChoice_valid_n11n,1)
        %    responseMatrix_noChoice_valid_n11n(tempi,:) = responseMatrix_noChoice_valid_n11n(tempi,:) ./ (1-responseMatrix_noChoice_valid_n11n(tempi,tempi));
        %    responseMatrix_noChoice_valid_n11n(tempi,tempi) = 0;
        %    score_stimuli_to_response(tempi,tempi) = 1;
        %end
        
        
        seqSigma = nan(sum(numSeq(1:valid_length)),1);
        for tempi=1:sum(numSeq(1:valid_length))
            temp_score = score_stimuli_to_response(tempi,1:sum(numSeq(1:valid_length)));
            temp_seqProb = responseMatrix_noChoice_valid_n11n(tempi,:);
            
            if isnan(temp_seqProb(1))
                continue
            end
            
            if if_precision_meanProb0_sumProb1 == 0
                x = temp_score';
                y = temp_seqProb';
                % cftool(x,y);
                
            elseif if_precision_meanProb0_sumProb1 == 1
                temp_score_unique = unique(temp_score);
                temp_seqProb_uniqueSum = nan(size(temp_score_unique));
                for tempj=1:length(temp_score_unique)
                    temptempBoolIndex = temp_score == temp_score_unique(tempj);
                    
                    temp_seqProb_uniqueSum(tempj) = sum(temp_seqProb(temptempBoolIndex));
                end
                x = temp_score_unique';
                y = temp_seqProb_uniqueSum';
                %cftool(x,y);
                
            end
            
            
            %mu = max(x);
            mu = 6;
            
            [xData, yData] = prepareCurveData(x,double(y));
            
            temp_str = sprintf('(1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-%.2f)/sigma)^2))',mu);
            
            %ft = fittype('(1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-6)/sigma)^2))',...
            %   'dependent',{'y'},'independent',{'x'},...
            %   'coefficients',{'sigma'});
            ft = fittype(temp_str,...
                'dependent',{'y'},'independent',{'x'},...
                'coefficients',{'sigma'});
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.StartPoint = 0.5;
            [fitresult, gof] = fit(xData,yData,ft,opts);
            
            temp_r2 = gof.rsquare;
            temp_sigma = fitresult.sigma;
            
            seqSigma(tempi) = temp_sigma;
            
            %if tempi == 25
            if tempi == 17
                a = 1;
            end
        end
        
        seqPrecision_behavior = 1./seqSigma;
        %seqPrecision_behavior = 1./(seqSigma.^2);
        
        
        
        %exampleSeq = 25;
        
        temp_score = score_stimuli_to_response(exampleSeq,1:sum(numSeq(1:valid_length)));
        temp_seqProb = responseMatrix_noChoice_valid_n11n(exampleSeq,:);
        
        temp_score_unique = unique(temp_score);
        temp_seqProb_uniqueSum = nan(size(temp_score_unique));
        for tempj=1:length(temp_score_unique)
            temptempBoolIndex = temp_score == temp_score_unique(tempj);
            
            temp_seqProb_uniqueSum(tempj) = sum(temp_seqProb(temptempBoolIndex));
        end
        
        if if_precision_meanProb0_sumProb1 == 0
            x = temp_score';
            y = temp_seqProb';
            
        elseif if_precision_meanProb0_sumProb1 == 1
            x = temp_score_unique';
            y = temp_seqProb_uniqueSum';
            
        end
        
        %mu = max(x);
        mu = 6;
        
        [xData, yData] = prepareCurveData(x,double(y));
        temp_str = sprintf('(1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-%.2f)/sigma)^2))',mu);
        ft = fittype(temp_str,...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'sigma'});
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = 0.5;
        [fitresult, gof] = fit(xData,yData,ft,opts);
        temp_r2 = gof.rsquare;
        temp_sigma = fitresult.sigma;
        
        x_fit = max(x):-0.1:min(x);
        y_fit = (1/(temp_sigma*sqrt(2*pi)))*exp(-0.5*(((x_fit-mu)/temp_sigma).^2));
        
        
        score_exampleSeq = temp_score';
        seqProb_exampleSeq = temp_seqProb';
        
        score_exampleSeq_fit = x_fit;
        seqProb_exampleSeq_fit = y_fit;
        
        score_exampleSeq_unique = temp_score_unique';
        seqProb_exampleSeq_uniqueSum = temp_seqProb_uniqueSum';
        
        score_exampleSeq_unique_fit = x_fit;
        seqProb_exampleSeq_uniqueSum_fit = y_fit;
        
        sigma_exampleSeq = temp_sigma;
        
        
        
        %%
        %responseMatrix_noChoice_valid_n11n;
        %seqProb_raw_exampleSeq = responseMatrix_noChoice_valid_n11n(exampleSeq,:);
        %seqProb_exampleSeq;
        
        [M,I] = sort(seqProb_exampleSeq,'descend');
        
        topX = sum(M>eps);
        I_topX = I(1:topX);
        
        seqProbIndex_exampleSeq_sortedTopX = I_topX;
        seqProb_exampleSeq_sortedTopX = M(1:topX);
        score_exampleSeq_sortedTopX = score_exampleSeq(seqProbIndex_exampleSeq_sortedTopX);
        
        
        [MB,IB] = sort(score_exampleSeq_sortedTopX,'descend');
        
        seqProbIndex_exampleSeq_sortedTopXB = seqProbIndex_exampleSeq_sortedTopX(IB);
        seqProb_exampleSeq_sortedTopXB = seqProb_exampleSeq_sortedTopX(IB);
        score_exampleSeq_sortedTopXB = MB;
        
        
        
    elseif if_entropy == 1
        
        p = ones(1,sum(numSeq(1:valid_length)))./sum(numSeq(1:valid_length));
        p = p + eps;
        p = p./sum(p);
        entropy_max = -sum(p.*log2(p));
        
        seqPrecision_behavior = nan(sum(numSeq(1:valid_length)),1);
        for tempi=1:sum(numSeq(1:valid_length))
            temp_seqProb = responseMatrix_noChoice_valid_n11n(tempi,:);
            
            %seqPrecision_neuron(tempi) = temp_seqProb(tempi);
            
            %p = [1 0 0 0 0 0];
            %p = [1 1 1 1 1 1]./6;
            p = temp_seqProb;
            p = p + eps;
            p = p./sum(p);
            entropy = -sum(p.*log2(p));
            
            %seqPrecision_behavior(tempi) = 1./entropy;
            %seqPrecision_behavior(tempi) = entropy_max-entropy;
            if if_entropy_behaviorInverted == 1
                seqPrecision_behavior(tempi) = (entropy_max-entropy)./entropy_max;
            elseif if_entropy_behaviorInverted == 0
                seqPrecision_behavior(tempi) = entropy./entropy_max;
            end
            
        end
        
        
    end
    
end



%% Plot
if if_singlePanel0_multiPanel1 == 1
    
    %% Trial-based length level
    %if true
    %if if_plot_submissionRT == 0
    if if_plot_multiPanelB == 1
        if true
            %fig1 = figure(1);
            fig = figure('Name','Trial-based length level','NumberTitle','off');
            %set(gcf,'Position',[0 50 figure_xPixel*4 figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[0 50 figure_xPixel*4*1.1 figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[0 50 figure_xPixel*4*1.1 figure_yPixel*0.9*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[0 50 figure_xPixel*4*1.1*0.87*0.99 figure_yPixel*0.9*0.95*0.975*0.975*0.99*0.975]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[0 50 figure_xPixel*4*1.1*0.87*0.99*0.84*1.03 figure_yPixel*0.9*0.95*0.975*0.975*0.99*0.975]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            set(gcf,'Position',[0 50 figure_xPixel*4*1.1*0.87*0.99*0.84*1.03 figure_yPixel*0.8975*0.97]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            %t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
            t = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
            %t = tiledlayout(1,4,'TileSpacing','tight','Padding','compact');
            
            %nexttile
            %set(gca, 'visible', 'off')
            
            nexttile
            %% Plot accuracy
            
            if true
                temp1 = internalChoiceAccuracy;
                temp2 = internalNoChoiceAccuracy;
                
                temptemp_p_Acc = nan(1,4);
                for tempi=1:4                    
                    [~,temptemp_p_Acc(tempi)] = ttest(temp1(:,tempi),temp2(:,tempi));                    
                end
                
            end
            
            plot(seq_length_rangeHead:seq_length_rangeTail,mean(offloadChoiceAccuracy),'Color',color_choiceOffload,'LineWidth',1);
            hold on
            plot(seq_length_rangeHead:seq_length_rangeTail,mean(internalChoiceAccuracy),'Color',color_choiceMemory,'LineWidth',1);
            hold on
            plot(seq_length_rangeHead:seq_length_rangeTail,mean(internalNoChoiceAccuracy),'Color',color_noChoice,'LineWidth',1);
            hold on
            
            % errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(internalNoChoiceAccuracy),std(internalNoChoiceAccuracy)./sqrt(size(internalNoChoiceAccuracy, 1)),'-o','Color',color_noChoice,'LineWidth',1,'CapSize',12,'MarkerSize',4);
            errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(offloadChoiceAccuracy),std(offloadChoiceAccuracy)./sqrt(size(offloadChoiceAccuracy, 1)),'-o','Color',color_choiceOffload,'LineWidth',1,'CapSize',6,'MarkerSize',1);
            hold on
            errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(internalChoiceAccuracy),std(internalChoiceAccuracy)./sqrt(size(internalChoiceAccuracy, 1)),'-o','Color',color_choiceMemory,'LineWidth',1,'CapSize',6,'MarkerSize',1);
            hold on
            errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(internalNoChoiceAccuracy),std(internalNoChoiceAccuracy)./sqrt(size(internalNoChoiceAccuracy, 1)),'-o','Color',color_noChoice,'LineWidth',1,'CapSize',6,'MarkerSize',1);
            hold on
            
            % errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(ProbOffloadingInSeqLength),std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)),'-o','Color',color_rProb,'LineWidth',1,'CapSize',7,'MarkerSize',2);
            % hold on
            
            if false
                [~,temp_p1] = ttest(internalChoiceAccuracy,internalNoChoiceAccuracy);  %#ok<*UNRCH>
                [~,temp_p2] = ttest2(internalChoiceAccuracy,internalNoChoiceAccuracy);  %#ok<*UNRCH>
            end
            
            if if_monkey_D0_Z1 == 0
                %le = legend('Choice-offload',...
                %    'Choice-memory','Forced-to-test',...
                %    'Location','southwest','fontsize',7);%9
                le = legend('Offload',...
                    'Memory','Forced-to-test',...
                    'Location','southwest','fontsize',6.5);%9
                le.ItemTokenSize = ones(1,2)*8;%10
                
                %le.Location = 'westoutside';
                %le.Position(1) = le.Position(1) - 0.2;
                %le.Position(2) = le.Position(2) + 0.14;
                
                legend('boxoff');
            end
            
            set(gca,'linewidth',1.5)
            ylim([0 1]);
            set(gca, 'FontSize', 8)
            %set(gca,'YTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
            set(gca,'YTick',[0:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
            % set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
            set(gca,'XLim',[seq_length_rangeHead-0.75 seq_length_rangeTail+0.75]);%X轴的数据显示范围
            set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
            % set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
            % set(gca,'XTickLabel',{'length1','length2','length3','length4'});%给坐标加标签
            %set(gca,'XTickLabel',{'1','2','3','4'});%给坐标加标签
            %text(4.15,-0.14,sprintf('length'),'fontsize',9);
            %ax = gca;
            %ax.XAxis.FontSize = 9;
            %ax.YAxis.FontSize = 10;%11
            
            %set(gca,'XTickLabel',{'Len1','Len2','Len3','Len4'});%给坐标加标签
            %xtickangle(0);
            %xtickangle(25);
            
            temp_y_min = 0;
            temp_y_max = 1;
            
%             xtl = ["1","2","3","4"];%Len1
%             xt=get(gca,'XTick');
%             yt=get(gca,'YTick');
%             xtext_xp=xt;
%             xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.28;
%             text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',7);%9-->12-->11
%             set(gca,'xticklabel','');
            
            
            set(gca,'box','off');% 取消右、上边框
            
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.1]);%0.3
            
            xlabel('Length', 'FontSize', 8);
            ylabel('Recall accuracy', 'FontSize', 8);
            if if_title == 1
                %title(sprintf('Recall accuracy'), 'FontSize', 9);
            end
            
        end
        
        %% Plot rProb
        if true
            %fig2 = figure(2);
            %set(gcf,'Position',[200 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            nexttile
            
            x = [];
            y = [];
            for tempi=1:4
                ProbOffloadingInSeqLength;
                if if_plot_rProb_fit == 1
                    temp_size = 10;
                    scatter(tempi,ProbOffloadingInSeqLength(:,tempi),...
                        temp_size,'filled','MarkerFaceColor',[1 1 1]*0.5,...
                        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
                    hold on
                end
                x = [x ones(1,size(ProbOffloadingInSeqLength,1))*tempi]; %#ok<*AGROW>
                y = [y ProbOffloadingInSeqLength(:,tempi)'];
            end
            if if_plot_rProb_fit == 1
                %y = a1;
                [r, p_corr] = corr(x',y'); %#ok<*ASGLU>
                
                n = 1;
                [p_mapping,S] = polyfit(x,y,n);
                r2 = 1 - (S.normr/norm(y - mean(y)))^2;
                
                temp_mdl = fitglm(x(~isnan(x)),y(~isnan(x)));
                r2 = temp_mdl.Rsquared.Adjusted;
                beta = temp_mdl.Coefficients.Estimate(2);
                p = temp_mdl.Coefficients.pValue(2);
                
                
                x_fit = 1:1:4;
                y_fit = polyval(p_mapping,x_fit);
                
                plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
            end
            
            if if_plot_rProb_fit == 0
                errorbar(seq_length_rangeHead: seq_length_rangeTail,...
                    mean(ProbOffloadingInSeqLength),std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)),...
                    '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'CapSize', 6, 'MarkerSize', 1);
            end
            
            set(gca,'box','off');% 取消右、上边框
            set(gca,'linewidth',1.5)
            
            temp_y_min = 0;
            temp_y_max = 1;
            
%             xtl = ["1","2","3","4"];%Len1
%             xt=get(gca,'XTick');
%             yt=get(gca,'YTick');
%             xtext_xp=xt;
%             xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.28;
%             text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',7);%9-->12-->11
%             set(gca,'xticklabel','');
            
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.1]);%0.3
            
            %ylim([0 1]);
            set(gca, 'FontSize', 8)
            %set(gca,'YTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
            % set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
            set(gca,'XLim',[seq_length_rangeHead-0.75 seq_length_rangeTail+0.75]);%X轴的数据显示范围
            set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
            % set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
            % set(gca,'XTickLabel',{'length1','length2','length3','length4'});%给坐标加标签
            %set(gca,'XTickLabel',{'1','2','3','4'});%给坐标加标签
            %text(4.15,-0.14,sprintf('length'),'fontsize',9);
            %xtickangle(0);
            
            %xticks([]);
            %xticklabels([]);
            
            %ax = gca;
            %ax.XAxis.FontSize = 9;
            %ax.YAxis.FontSize = 10;%11
            
            if if_plot_rProb_scale == 1
                y_max = max(mean(ProbOffloadingInSeqLength)+std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)));
                y_min = min(mean(ProbOffloadingInSeqLength)-std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)));
                
                %ylim([y_min-(y_max-y_min)*0.6 y_max+(y_max-y_min)*0.6]);%0.3
                %set(gca,'YTick',[0:0.2:1]);
                set(gca,'YTick',[0,1]);
            else
                set(gca,'YTick',[0:0.2:1]);
            end
            
            
            
            
            xlabel('Length', 'FontSize', 8);
            ylabel('Offloading rate', 'FontSize', 8);
            if if_title == 1
                %title(sprintf('Offloading rate'), 'FontSize', 9);
            end
            
        end
        
        %% Plot correlation
        if true
            %fig3 = figure(3);
            %set(gcf,'Position',[400 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            nexttile
            
            % isvalid = (~isnan(gError_noChoice_collapsed_inOne)) & (~isnan(offloadingProb_all_inOne));
            isvalid = (~isnan(gAcc_noChoice_collapsed_inOne)) & (~isnan(offloadingProb_all_inOne));
            
            % x = gError_noChoice_collapsed_inOne(isvalid);
            x = gAcc_noChoice_collapsed_inOne(isvalid);
            y = offloadingProb_all_inOne(isvalid);
            n = 1;
            [p_mapping,S] = polyfit(x,y,n);
            r2 = 1 - (S.normr/norm(y - mean(y)))^2;
            [r, p_corr] = corr(x',y');
            tempTxt = 'p > 0.05';
            if p_corr < 0.05
                tempTxt = 'p < 0.05';
            end
            if p_corr < 0.01
                tempTxt = 'p < 0.01';
            end
            if p_corr < 0.001
                tempTxt = 'p < 0.001';
            end
            
            [p_mapping_poly,S_poly] = polyfit(x,y,3);
            r2_poly = 1 - (S_poly.normr/norm(y - mean(y)))^2;
            
            
            
            
            x_fit = 0:0.01:1;
            y_fit = polyval(p_mapping,x_fit);
            
            temp_color = [[166,97,26]/255;
                [223,194,125]/255;
                [128,205,193]/255;
                [1,133,113]/255];
            
            
            for target_seqLength=1:pointKindsNum
                %temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(gError_noChoice_collapsed{target_seqLength}));
                temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(gAcc_noChoice_collapsed{target_seqLength}));
                
                temp_size = temp_size ./ 5;
                
                %scatter(gError_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
                scatter(gAcc_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
                    temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
                hold on
            end
            plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
            
            % text(0.1,0.9,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',10,'FontWeight','bold');
            % text(0.2,0.1,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',10,'FontWeight','bold');
            %text(0.62,0.97,sprintf('r=%.3f',r), 'fontsize',10,'FontWeight','normal');
            %text(0.62,0.78,sprintf('%s',tempTxt), 'fontsize',10,'FontWeight','normal');
            
            set(gca,'linewidth',1.5)
            % xlim([0 1]);
            x_min = 0;
            x_max = 1;
            xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05]);
            
            %         y_min = 0;
            %         y_max = 1;
            %         %ylim([0 1]);
            %         ylim([y_min-(y_max-y_min)*0.05 y_max+(y_max-y_min)*0.05]);
            
            temp_y_min = 0;
            temp_y_max = 1;
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.1]);%0.3
            
            set(gca, 'FontSize', 8)
            %set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
            set(gca,'XTick',[0,0.5,1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
            %set(gca,'XTick',[0,0.5]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
            
            xtickangle(0);
            
            % set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
            %set(gca,'YTick',[0:0.5:1]);%设置要显示坐标刻度的范围
            set(gca,'YTick',[0:1]);%设置要显示坐标刻度的范围
            set(gca,'box','off');% 取消右、上边框
            % xlabel('Error rate', 'FontSize', 11);
            xlabel('Recall accuracy', 'FontSize', 8);
            %text(0.5,-0.27,sprintf('Recall accuracy'),'fontsize',8,'HorizontalAlignment','center');
            
            temp_ylabel = ylabel('Offloading rate', 'FontSize', 8);
            % temp_Position = temp_ylabel.Position;
            % temp_Position(1) = temp_Position(1) - 0.03 + 1;
            % temp_Position(2) = temp_Position(2) - 0.01 - 0.2;
            % temp_ylabel.Position = temp_Position;
            if if_title == 1
                %title('Correlation', 'FontSize', 10);
            end
            title(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 7.5);
            
        end
        
        
        
        %% Plot correlation with precision
        if if_plot_precision == 1
            nexttile
            isvalid = (~isnan(seqPrecision_behavior')) & (~isnan(offloadingProb_all_inOne));
            
            x = seqPrecision_behavior(isvalid)';
            y = offloadingProb_all_inOne(isvalid);
            n = 1;
            [p_mapping,S] = polyfit(x,y,n);
            r2 = 1 - (S.normr/norm(y - mean(y)))^2;
            [r, p_corr] = corr(x',y');
            tempTxt = 'p > 0.05';
            if p_corr < 0.05
                tempTxt = 'p < 0.05';
            end
            if p_corr < 0.01
                tempTxt = 'p < 0.01';
            end
            if p_corr < 0.001
                tempTxt = 'p < 0.001';
            end
            
            [p_mapping_poly,S_poly] = polyfit(x,y,3);
            r2_poly = 1 - (S_poly.normr/norm(y - mean(y)))^2;
            
            
            
            if if_entropy_behaviorInverted == 1
                x_min = 0.3;
                x_max = 1;
            elseif if_entropy_behaviorInverted == 0
                x_min = 0;
                x_max = 0.7;
            end
            
            %x_fit = 0:0.01:1;
            x_fit = x_min:0.01:x_max;
            y_fit = polyval(p_mapping,x_fit);
            
            temp_color = [[166,97,26]/255;
                [223,194,125]/255;
                [128,205,193]/255;
                [1,133,113]/255];
            
            
            for target_seqLength=1:pointKindsNum
                temp_range = (sum(numSeq(1:target_seqLength-1))+1):sum(numSeq(1:target_seqLength));
                
                temptempBoolIndex = false(sum(numSeq),1);
                
                temptempBoolIndex(temp_range) = true;
                temptempBoolIndex(~isvalid) = false;
                
                temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(gAcc_noChoice_collapsed{target_seqLength}));
                
                temp_size = temp_size ./ 5;
                
                scatter(seqPrecision_behavior(temptempBoolIndex)', offloadingProb_all_inOne(temptempBoolIndex), ...
                    temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
                hold on
            end
            plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
            
            % text(0.1,0.9,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',10,'FontWeight','bold');
            % text(0.2,0.1,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',10,'FontWeight','bold');
            %text(0.62,0.97,sprintf('r=%.3f',r), 'fontsize',10,'FontWeight','normal');
            %text(0.62,0.78,sprintf('%s',tempTxt), 'fontsize',10,'FontWeight','normal');
            
            set(gca,'linewidth',1.5)
            % xlim([0 1]);
            %x_min = 0.4;

            xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05]);
            
            %         y_min = 0;
            %         y_max = 1;
            %         %ylim([0 1]);
            %         ylim([y_min-(y_max-y_min)*0.05 y_max+(y_max-y_min)*0.05]);
            
            temp_y_min = 0;
            temp_y_max = 1;
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.1]);%0.3
            
            set(gca, 'FontSize', 8)
            %set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
            %set(gca,'XTick',[0:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
            set(gca,'XTick',[x_min,(x_min+x_max)/2,x_max]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
            
            xtickangle(0);
            
            % set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
            %set(gca,'YTick',[0:0.5:1]);%设置要显示坐标刻度的范围
            set(gca,'YTick',[0:1]);%设置要显示坐标刻度的范围
            set(gca,'box','off');% 取消右、上边框
            % xlabel('Error rate', 'FontSize', 11);
            xlabel('Recall variability', 'FontSize', 8);
            %text((x_min+x_max)/2,-0.27,sprintf('Recall variability'),'fontsize',8,'HorizontalAlignment','center');
            
            
            temp_ylabel = ylabel('Offloading rate', 'FontSize', 8);
            % temp_Position = temp_ylabel.Position;
            % temp_Position(1) = temp_Position(1) - 0.03 + 1;
            % temp_Position(2) = temp_Position(2) - 0.01 - 0.2;
            % temp_ylabel.Position = temp_Position;
            if if_title == 1
                %title('Correlation', 'FontSize', 10);
            end
            title(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 7.5);
            
        end
        
    end
    
    
    if if_plot_multiPanelC == 1
        %fig1 = figure(1);
        fig = figure('Name','Trial-based length level','NumberTitle','off');
        %set(gcf,'Position',[0 50 figure_xPixel*4 figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[0 50 figure_xPixel*4*1.1 figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[0 50 figure_xPixel*4*1.1 figure_yPixel*0.9*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[0 50 figure_xPixel*3*1.1 figure_yPixel*0.9*0.95*2*0.8*0.9*0.96]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[0 50 figure_xPixel*3*1.1*1.15 figure_yPixel*0.9*0.95*2*0.8*0.9*0.96*1.1*1.1]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[0 50 figure_xPixel*3*1.1*1.15 figure_yPixel*0.9*0.95*2*0.8*0.9*0.96*1.1*1.1]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
        %t = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
        %t = tiledlayout(1,4,'TileSpacing','tight','Padding','compact');
        t = tiledlayout(2,4,'TileSpacing','tight','Padding','compact');
        
        %nexttile
        %set(gca, 'visible', 'off')
        
        nexttile
        %% Plot accuracy
        
        % errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(internalNoChoiceAccuracy),std(internalNoChoiceAccuracy)./sqrt(size(internalNoChoiceAccuracy, 1)),'-o','Color',color_noChoice,'LineWidth',1,'CapSize',12,'MarkerSize',4);
        errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(offloadChoiceAccuracy),std(offloadChoiceAccuracy)./sqrt(size(offloadChoiceAccuracy, 1)),'-o','Color',color_choiceOffload,'LineWidth',1,'CapSize',6,'MarkerSize',1);
        hold on
        errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(internalChoiceAccuracy),std(internalChoiceAccuracy)./sqrt(size(internalChoiceAccuracy, 1)),'-o','Color',color_choiceMemory,'LineWidth',1,'CapSize',6,'MarkerSize',1);
        hold on
        errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(internalNoChoiceAccuracy),std(internalNoChoiceAccuracy)./sqrt(size(internalNoChoiceAccuracy, 1)),'-o','Color',color_noChoice,'LineWidth',1,'CapSize',6,'MarkerSize',1);
        hold on
        
        % errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(ProbOffloadingInSeqLength),std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)),'-o','Color',color_rProb,'LineWidth',1,'CapSize',7,'MarkerSize',2);
        % hold on
        
        if false
            [~,temp_p1] = ttest(internalChoiceAccuracy,internalNoChoiceAccuracy);  %#ok<*UNRCH>
            [~,temp_p2] = ttest2(internalChoiceAccuracy,internalNoChoiceAccuracy);  %#ok<*UNRCH>
        end
        
        le = legend('Choice-offload',...
            'Choice-memory','Forced-to-test',...
            'Location','southwest','fontsize',6.5);%9
        le.ItemTokenSize = ones(1,2)*10;
        
        legend('boxoff');
        
        %%le.Location = 'westoutside';
        %le.Position(1) = le.Position(1) - 0.2;
        %le.Position(2) = le.Position(2) + 0.14;
        
        set(gca,'linewidth',1.5)
        ylim([0 1]);
        set(gca, 'FontSize', 8)
        %set(gca,'YTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        set(gca,'YTick',[0:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        % set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XLim',[seq_length_rangeHead-0.75 seq_length_rangeTail+0.75]);%X轴的数据显示范围
        set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
        % set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
        % set(gca,'XTickLabel',{'length1','length2','length3','length4'});%给坐标加标签
        %set(gca,'XTickLabel',{'1','2','3','4'});%给坐标加标签
        %text(4.15,-0.14,sprintf('length'),'fontsize',9);
        ax = gca;
        ax.XAxis.FontSize = 9;
        ax.YAxis.FontSize = 10;%11
        
        %set(gca,'XTickLabel',{'Len1','Len2','Len3','Len4'});%给坐标加标签
        %xtickangle(0);
        %xtickangle(25);
        
        temp_y_min = 0;
        temp_y_max = 1;
        
        xtl = ["Len1","Len2","Len3","Len4"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.16;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',8);%9-->12-->11
        set(gca,'xticklabel','');
        
        
        set(gca,'box','off');% 取消右、上边框
        %xlabel('Length', 'FontSize', 10);
        % ylabel('Recall accuracy', 'FontSize', 11);
        if if_title == 1
            title(sprintf('Recall accuracy'), 'FontSize', 9);
        end
        
        %% Plot seq-level gAcc_choice VS. gAcc_noChoice
        nexttile
        
        
        q1 = gAcc_choice_collapsed_inOne';
        q2 = gAcc_noChoice_collapsed_inOne';
        
        q1_mean = mean(q1,'omitnan');
        q2_mean = mean(q2,'omitnan');
        
        temp_threshold = 0;
        
        temptempBoolIndex1 = (q1>temp_threshold) & (q2>temp_threshold);
        [~,temp_p_choiceNoChoice] = ttest(q1(temptempBoolIndex1),q2(temptempBoolIndex1));
        
        
        temp_p = temp_p_choiceNoChoice;
        temptempBoolIndex = (~isnan(q1)) & (~isnan(q2)) & temptempBoolIndex1;
        temp_1 = q1(temptempBoolIndex);
        temp_2 = q2(temptempBoolIndex);
        
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        %temp_y_min = 0;
        %temp_y_max = 1;
        
        % paired line plot
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.65])
        %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
        %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
        ylim([-0.15 1.15]);
        %ylim([temp_y_min temp_y_max]);
        set(gca, 'FontSize', 10)
        
        xticks([1 2]);
        
        %xtl = ["Choice-memory"; "Forced-to-test"];
        xtl = ["Memory"; "Forced-to-test"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.435;%0.56,0.4
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.375;%0.56,0.4
        
        if if_monkey_D0_Z1 == 0
            %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.375;%0.56,0.4
            xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.21;%0.56,0.4
        elseif if_monkey_D0_Z1 == 1
            %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.215;%0.56,0.4
            xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.27;%0.56,0.4
        end
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
        %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',7.5);%25
        set(gca,'xticklabel','');
        
        
        yticks([0 1]);
        
        
        set(gca,'box','off');% 取消右、上边框
        %ylabel('Recall accuracy', 'FontSize', 9, 'FontWeight', 'normal');
        %temp1 = title(sprintf('Seq accuracy'),'fontsize',9);
        
        %% Plot rProb
        %fig2 = figure(2);
        %set(gcf,'Position',[200 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        nexttile
        
        x = [];
        y = [];
        for tempi=1:4
            ProbOffloadingInSeqLength;
            if if_plot_rProb_fit == 1
                temp_size = 10;
                scatter(tempi,ProbOffloadingInSeqLength(:,tempi),...
                    temp_size,'filled','MarkerFaceColor',[1 1 1]*0.5,...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
                hold on
            end
            x = [x ones(1,size(ProbOffloadingInSeqLength,1))*tempi]; %#ok<*AGROW>
            y = [y ProbOffloadingInSeqLength(:,tempi)'];
        end
        if if_plot_rProb_fit == 1
            %y = a1;
            [r, p_corr] = corr(x',y'); %#ok<*ASGLU>
            
            n = 1;
            [p_mapping,S] = polyfit(x,y,n);
            r2 = 1 - (S.normr/norm(y - mean(y)))^2;
            
            temp_mdl = fitglm(x(~isnan(x)),y(~isnan(x)));
            r2 = temp_mdl.Rsquared.Adjusted;
            beta = temp_mdl.Coefficients.Estimate(2);
            p = temp_mdl.Coefficients.pValue(2);
            
            
            x_fit = 1:1:4;
            y_fit = polyval(p_mapping,x_fit);
            
            plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        end
        
        if if_plot_rProb_fit == 0
            errorbar(seq_length_rangeHead: seq_length_rangeTail,...
                mean(ProbOffloadingInSeqLength),std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)),...
                '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'CapSize', 6, 'MarkerSize', 1);
        end
        
        set(gca,'linewidth',1.5)
        ylim([0 1]);
        set(gca, 'FontSize', 8)
        %set(gca,'YTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        % set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XLim',[seq_length_rangeHead-0.75 seq_length_rangeTail+0.75]);%X轴的数据显示范围
        set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
        % set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
        % set(gca,'XTickLabel',{'length1','length2','length3','length4'});%给坐标加标签
        %set(gca,'XTickLabel',{'1','2','3','4'});%给坐标加标签
        %text(4.15,-0.14,sprintf('length'),'fontsize',9);
        %xtickangle(0);
        
        %xticks([]);
        xticklabels([]);
        
        ax = gca;
        ax.XAxis.FontSize = 9;
        ax.YAxis.FontSize = 10;%11
        
        if if_plot_rProb_scale == 1
            y_max = max(mean(ProbOffloadingInSeqLength)+std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)));
            y_min = min(mean(ProbOffloadingInSeqLength)-std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)));
            
            ylim([y_min-(y_max-y_min)*0.6 y_max+(y_max-y_min)*0.6]);%0.3
            set(gca,'YTick',[0:0.2:1]);
        else
            set(gca,'YTick',[0:0.2:1]);
        end
        
        set(gca,'box','off');% 取消右、上边框
        %xlabel('Length', 'FontSize', 11);
        % ylabel('Offloading rate', 'FontSize', 11);
        if if_title == 1
            title(sprintf('Offloading rate'), 'FontSize', 9);
        end
        
        %% Plot correlation
        %fig3 = figure(3);
        %set(gcf,'Position',[400 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        nexttile
        
        % isvalid = (~isnan(gError_noChoice_collapsed_inOne)) & (~isnan(offloadingProb_all_inOne));
        isvalid = (~isnan(gAcc_noChoice_collapsed_inOne)) & (~isnan(offloadingProb_all_inOne));
        
        % x = gError_noChoice_collapsed_inOne(isvalid);
        x = gAcc_noChoice_collapsed_inOne(isvalid);
        y = offloadingProb_all_inOne(isvalid);
        n = 1;
        [p_mapping,S] = polyfit(x,y,n);
        r2 = 1 - (S.normr/norm(y - mean(y)))^2;
        [r, p_corr] = corr(x',y');
        tempTxt = 'p > 0.05';
        if p_corr < 0.05
            tempTxt = 'p < 0.05';
        end
        if p_corr < 0.01
            tempTxt = 'p < 0.01';
        end
        if p_corr < 0.001
            tempTxt = 'p < 0.001';
        end
        
        [p_mapping_poly,S_poly] = polyfit(x,y,3);
        r2_poly = 1 - (S_poly.normr/norm(y - mean(y)))^2;
        
        
        
        
        x_fit = 0:0.01:1;
        y_fit = polyval(p_mapping,x_fit);
        
        temp_color = [[166,97,26]/255;
            [223,194,125]/255;
            [128,205,193]/255;
            [1,133,113]/255];
        
        
        for target_seqLength=1:pointKindsNum
            %temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(gError_noChoice_collapsed{target_seqLength}));
            temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(gAcc_noChoice_collapsed{target_seqLength}));
            
            temp_size = temp_size ./ 5;
            
            %scatter(gError_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
            scatter(gAcc_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
                temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
        end
        plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        
        % text(0.1,0.9,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',10,'FontWeight','bold');
        % text(0.2,0.1,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',10,'FontWeight','bold');
        %text(0.62,0.97,sprintf('r=%.3f',r), 'fontsize',10,'FontWeight','normal');
        %text(0.62,0.78,sprintf('%s',tempTxt), 'fontsize',10,'FontWeight','normal');
        
        set(gca,'linewidth',1.5)
        % xlim([0 1]);
        x_min = 0;
        x_max = 1;
        xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05]);
        
        y_min = 0;
        y_max = 1;
        %ylim([0 1]);
        ylim([y_min-(y_max-y_min)*0.05 y_max+(y_max-y_min)*0.05]);
        
        set(gca, 'FontSize', 8)
        %set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        set(gca,'XTick',[0:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        %set(gca,'XTick',[0,0.5]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        
        xtickangle(0);
        
        % set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        %set(gca,'YTick',[0:0.5:1]);%设置要显示坐标刻度的范围
        set(gca,'YTick',[0:1]);%设置要显示坐标刻度的范围
        set(gca,'box','off');% 取消右、上边框
        % xlabel('Error rate', 'FontSize', 11);
        %xlabel('Recall accuracy', 'FontSize', 10);
        %text(0.7,-0.17,sprintf('Recall accuracy'),'fontsize',9);
        
        temp_ylabel = ylabel('Offloading rate', 'FontSize', 9);
        % temp_Position = temp_ylabel.Position;
        % temp_Position(1) = temp_Position(1) - 0.03 + 1;
        % temp_Position(2) = temp_Position(2) - 0.01 - 0.2;
        % temp_ylabel.Position = temp_Position;
        if if_title == 1
            %title('Correlation', 'FontSize', 10);
        end
        title(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 9);
        
    end
    
    
    %% Seq-based length level
    if false
        offloadingAcc_all_inOne;
        gAcc_choice_collapsed_inOne;
        gAcc_noChoice_collapsed_inOne;
        offloadingProb_all_inOne;
        
        
        
        %fig2 = figure(2);
        fig = figure('Name','Seq-based length level','NumberTitle','off');
        set(gcf,'Position',[0 300 figure_xPixel*4 figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
        t = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
        
        nexttile
        set(gca, 'visible', 'off')
        
        nexttile
        %% Plot accuracy
        
        a1 = offloadingAcc_all_inOne;
        a2 = gAcc_choice_collapsed_inOne;
        a3 = gAcc_noChoice_collapsed_inOne;
        
        b1_mean = nan(1,4);
        b2_mean = nan(1,4);
        b3_mean = nan(1,4);
        
        b1_sem = nan(1,4);
        b2_sem = nan(1,4);
        b3_sem = nan(1,4);
        
        for tempi=1:4
            temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
            
            temp_b1 = a1(temp_range);
            b1_mean(tempi) = mean(temp_b1);
            b1_sem(tempi) = std(temp_b1)./sqrt(length(temp_b1));
            
            temp_b2 = a2(temp_range);
            b2_mean(tempi) = mean(temp_b2);
            b2_sem(tempi) = std(temp_b2)./sqrt(length(temp_b2));
            
            temp_b3 = a3(temp_range);
            b3_mean(tempi) = mean(temp_b3);
            b3_sem(tempi) = std(temp_b3)./sqrt(length(temp_b3));
        end
        
        errorbar(seq_length_rangeHead:seq_length_rangeTail,b1_mean,b1_sem,'-o','Color',color_choiceOffload,'LineWidth',1,'CapSize',6,'MarkerSize',1);
        hold on
        errorbar(seq_length_rangeHead:seq_length_rangeTail,b2_mean,b2_sem,'-o','Color',color_choiceMemory,'LineWidth',1,'CapSize',6,'MarkerSize',1);
        hold on
        errorbar(seq_length_rangeHead:seq_length_rangeTail,b3_mean,b3_sem,'-o','Color',color_noChoice,'LineWidth',1,'CapSize',6,'MarkerSize',1);
        hold on
        
        
        if false
            [~,temp_p1] = ttest(gAcc_choice_collapsed_inOne,gAcc_noChoice_collapsed_inOne);  %#ok<*UNRCH>
            [~,temp_p2] = ttest2(gAcc_choice_collapsed_inOne,gAcc_noChoice_collapsed_inOne);  %#ok<*UNRCH>
            
            temp_p1 = [0 0 0 0];
            temp_p2 = [0 0 0 0];
            
            for tempi=1:4
                temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
                
                temp_b2 = a2(temp_range);
                temp_b3 = a3(temp_range);
                
                [~,temp_p1(tempi)] = ttest(temp_b2,temp_b3);  %#ok<*UNRCH>
                [~,temp_p2(tempi)] = ttest2(temp_b2,temp_b3);  %#ok<*UNRCH>
            end
            
        end
        
        
        le = legend('Choice-offload',...
            'Choice-memory','Forced-to-test',...
            'Location','southwest','fontsize',9);
        le.ItemTokenSize = ones(1,2)*10;
        
        %le.Location = 'westoutside';
        le.Position(1) = le.Position(1) - 0.2;
        le.Position(2) = le.Position(2) + 0.14;
        
        set(gca,'linewidth',1.5)
        ylim([0 1]);
        set(gca, 'FontSize', 11)
        set(gca,'YTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        % set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XLim',[seq_length_rangeHead-0.75 seq_length_rangeTail+0.75]);%X轴的数据显示范围
        set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
        % set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
        % set(gca,'XTickLabel',{'length1','length2','length3','length4'});%给坐标加标签
        set(gca,'XTickLabel',{'1','2','3','4'});%给坐标加标签
        %text(4.15,-0.14,sprintf('length'),'fontsize',9);
        ax = gca;
        ax.XAxis.FontSize = 9;
        ax.YAxis.FontSize = 11;
        xtickangle(0);
        
        set(gca,'box','off');% 取消右、上边框
        xlabel('Length', 'FontSize', 11);
        % ylabel('Recall accuracy', 'FontSize', 11);
        if if_title == 1
            title(sprintf('Recall accuracy'), 'FontSize', 11);
        end
        
        
        %% Plot rProb
        %fig2 = figure(2);
        %set(gcf,'Position',[200 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        nexttile
        
        a1 = offloadingProb_all_inOne;
        b1_mean = nan(1,4);
        b1_sem = nan(1,4);
        x = [];
        for tempi=1:4
            temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
            
            temp_b1 = a1(temp_range);
            b1_mean(tempi) = mean(temp_b1);
            b1_sem(tempi) = std(temp_b1)./sqrt(length(temp_b1));
            
            if if_plot_rProb_fit == 1
                temp_size = 10;
                scatter(tempi,temp_b1,...
                    temp_size,'filled','MarkerFaceColor',[1 1 1]*0.5,...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
                hold on
            end
            x = [x ones(1,numSeq(tempi))*tempi];
        end
        if if_plot_rProb_fit == 1
            y = a1;
            [r, p_corr] = corr(x',y'); %#ok<*ASGLU>
            
            n = 1;
            [p_mapping,S] = polyfit(x,y,n);
            r2 = 1 - (S.normr/norm(y - mean(y)))^2;
            
            temp_mdl = fitglm(x(~isnan(x)),y(~isnan(x)));
            r2 = temp_mdl.Rsquared.Adjusted;
            beta = temp_mdl.Coefficients.Estimate(2);
            p = temp_mdl.Coefficients.pValue(2);
            
            
            x_fit = 1:1:4;
            y_fit = polyval(p_mapping,x_fit);
            
            plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        end
        
        if if_plot_rProb_fit == 0
            errorbar(seq_length_rangeHead: seq_length_rangeTail,b1_mean,b1_sem,...
                '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'CapSize', 6, 'MarkerSize', 1);
        end
        
        set(gca,'linewidth',1.5)
        ylim([0 1]);
        set(gca, 'FontSize', 11)
        %set(gca,'YTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        % set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XLim',[seq_length_rangeHead-0.75 seq_length_rangeTail+0.75]);%X轴的数据显示范围
        set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
        % set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
        % set(gca,'XTickLabel',{'length1','length2','length3','length4'});%给坐标加标签
        set(gca,'XTickLabel',{'1','2','3','4'});%给坐标加标签
        %text(4.15,-0.14,sprintf('length'),'fontsize',9);
        ax = gca;
        ax.XAxis.FontSize = 9;
        ax.YAxis.FontSize = 11;
        xtickangle(0);
        
        
        if if_plot_rProb_scale == 1
            y_max = max(b1_mean+b1_sem);
            y_min = min(b1_mean-b1_sem);
            
            ylim([y_min-(y_max-y_min)*0.8 y_max+(y_max-y_min)*0.8]);%0.3
            set(gca,'YTick',[0:0.2:1]);
        else
            set(gca,'YTick',[0:0.2:1]);
        end
        
        
        set(gca,'box','off');% 取消右、上边框
        xlabel('Length', 'FontSize', 11);
        % ylabel('Offloading rate', 'FontSize', 11);
        if if_title == 1
            title(sprintf('Offloading rate'), 'FontSize', 11);
        end
        
        %% Plot correlation
        %fig3 = figure(3);
        %set(gcf,'Position',[400 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        nexttile
        
        % isvalid = (~isnan(gError_noChoice_collapsed_inOne)) & (~isnan(offloadingProb_all_inOne));
        isvalid = (~isnan(gAcc_noChoice_collapsed_inOne)) & (~isnan(offloadingProb_all_inOne));
        
        % x = gError_noChoice_collapsed_inOne(isvalid);
        x = gAcc_noChoice_collapsed_inOne(isvalid);
        y = offloadingProb_all_inOne(isvalid);
        n = 1;
        [p_mapping,S] = polyfit(x,y,n);
        r2 = 1 - (S.normr/norm(y - mean(y)))^2;
        [r, p_corr] = corr(x',y');
        tempTxt = 'p>0.05';
        if p_corr < 0.05
            tempTxt = 'p<0.05';
        end
        if p_corr < 0.01
            tempTxt = 'p<0.01';
        end
        if p_corr < 0.001
            tempTxt = 'p<0.001';
        end
        
        [p_mapping_poly,S_poly] = polyfit(x,y,3);
        r2_poly = 1 - (S_poly.normr/norm(y - mean(y)))^2;
        
        
        
        
        x_fit = 0:0.01:1;
        y_fit = polyval(p_mapping,x_fit);
        
        temp_color = [[166,97,26]/255;
            [223,194,125]/255;
            [128,205,193]/255;
            [1,133,113]/255];
        
        
        for target_seqLength=1:pointKindsNum
            %temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(gError_noChoice_collapsed{target_seqLength}));
            temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(gAcc_noChoice_collapsed{target_seqLength}));
            
            temp_size = temp_size ./ 5;
            
            %scatter(gError_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
            scatter(gAcc_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
                temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
        end
        plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        
        % text(0.1,0.9,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',10,'FontWeight','bold');
        % text(0.2,0.1,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',10,'FontWeight','bold');
        text(0.66,0.99,sprintf('r=%.3f',r), 'fontsize',10,'FontWeight','bold');
        text(0.66,0.80,sprintf('%s',tempTxt), 'fontsize',10,'FontWeight','bold');
        
        set(gca,'linewidth',1.5)
        % xlim([0 1]);
        x_min = 0;
        x_max = 1;
        xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
        
        ylim([0 1]);
        set(gca, 'FontSize', 11)
        set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        %set(gca,'XTick',[0,0.5]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        
        xtickangle(0);
        
        % set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        set(gca,'YTick',[0:0.5:1]);%设置要显示坐标刻度的范围
        set(gca,'box','off');% 取消右、上边框
        % xlabel('Error rate', 'FontSize', 11);
        xlabel('Recall accuracy', 'FontSize', 11);
        %text(0.7,-0.17,sprintf('Recall accuracy'),'fontsize',9);
        
        temp_ylabel = ylabel('Offloading rate', 'FontSize', 11);
        % temp_Position = temp_ylabel.Position;
        % temp_Position(1) = temp_Position(1) - 0.03 + 1;
        % temp_Position(2) = temp_Position(2) - 0.01 - 0.2;
        % temp_ylabel.Position = temp_Position;
        if if_title == 1
            title('Correlation', 'FontSize', 11);
        end
        
    end
    
    
    %     %% Seq-based X Session-based length level
    %     if true
    %         offloadingAcc_all_inOne;
    %         gAcc_choice_collapsed_inOne;
    %         gAcc_noChoice_collapsed_inOne;
    %         offloadingProb_all_inOne;
    %
    %         %offloadingAcc;
    %
    %         offloadingCorrect_trial_count_newSeq;
    %         offloading_trial_count_newSeq;
    %         gAcc_noChoice;
    %         gAcc_choice;
    %         offloadingProb;
    %
    %         fileSize = size(offloadingProb,1);
    %         c1 = cell(1,4);
    %         c2 = cell(1,4);
    %         c3 = cell(1,4);
    %         c4 = cell(1,4);
    %         c1_mean = nan(1,4);
    %         c2_mean = nan(1,4);
    %         c3_mean = nan(1,4);
    %         c4_mean = nan(1,4);
    %         c1_sem = nan(1,4);
    %         c2_sem = nan(1,4);
    %         c3_sem = nan(1,4);
    %         c4_sem = nan(1,4);
    %         for tempj=1:4
    %             for tempi=1:fileSize
    %                 temp1 = offloadingCorrect_trial_count_newSeq{tempi,tempj} ...
    %                     ./ offloading_trial_count_newSeq{tempi,tempj};
    %                 c1{tempj} = [c1{tempj} temp1'];
    %
    %                 temp2 = gAcc_noChoice{tempi,tempj};
    %                 c2{tempj} = [c2{tempj} temp2'];
    %
    %                 temp3 = gAcc_choice{tempi,tempj};
    %                 c3{tempj} = [c3{tempj} temp3'];
    %
    %                 temp4 = offloadingProb{tempi,tempj};
    %                 c4{tempj} = [c4{tempj} temp4'];
    %             end
    %             c1_mean(tempj) = mean(c1{tempj},'omitnan');
    %             c2_mean(tempj) = mean(c2{tempj},'omitnan');
    %             c3_mean(tempj) = mean(c3{tempj},'omitnan');
    %             c4_mean(tempj) = mean(c4{tempj},'omitnan');
    %
    %             c1_sem(tempj) = std(c1{tempj},1,'omitnan')./sqrt(length(c1{tempj}));
    %             c2_sem(tempj) = std(c2{tempj},1,'omitnan')./sqrt(length(c2{tempj}));
    %             c3_sem(tempj) = std(c3{tempj},1,'omitnan')./sqrt(length(c3{tempj}));
    %             c4_sem(tempj) = std(c4{tempj},1,'omitnan')./sqrt(length(c4{tempj}));
    %         end
    %
    %         %fig2 = figure(2);
    %         fig = figure('Name','Seq-based X Session-based length level','NumberTitle','off');
    %         set(gcf,'Position',[0 550 figure_xPixel*4 figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %         %t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    %         t = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
    %
    %         nexttile
    %         set(gca, 'visible', 'off')
    %
    %         nexttile
    %         %% Plot accuracy
    %
    %         errorbar(seq_length_rangeHead:seq_length_rangeTail,c1_mean,c1_sem,'-o','Color',color_choiceOffload,'LineWidth',1,'CapSize',6,'MarkerSize',1);
    %         hold on
    %         errorbar(seq_length_rangeHead:seq_length_rangeTail,c2_mean,c2_sem,'-o','Color',color_choiceMemory,'LineWidth',1,'CapSize',6,'MarkerSize',1);
    %         hold on
    %         errorbar(seq_length_rangeHead:seq_length_rangeTail,c3_mean,c3_sem,'-o','Color',color_noChoice,'LineWidth',1,'CapSize',6,'MarkerSize',1);
    %         hold on
    %
    %
    %         %if false
    %         %   [~,temp_p1] = ttest(internalChoiceAccuracy,internalNoChoiceAccuracy);  %#ok<*UNRCH>
    %         %   [~,temp_p2] = ttest2(internalChoiceAccuracy,internalNoChoiceAccuracy);  %#ok<*UNRCH>
    %         %end
    %
    %         le = legend('Choice-offload',...
    %             'Choice-memory','Forced-to-test',...
    %             'Location','southwest','fontsize',9);
    %         le.ItemTokenSize = ones(1,2)*10;
    %
    %         %le.Location = 'westoutside';
    %         le.Position(1) = le.Position(1) - 0.2;
    %         le.Position(2) = le.Position(2) + 0.14;
    %
    %         set(gca,'linewidth',1.5)
    %         ylim([0 1]);
    %         set(gca, 'FontSize', 11)
    %         set(gca,'YTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    %         % set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    %         set(gca,'XLim',[seq_length_rangeHead-0.75 seq_length_rangeTail+0.75]);%X轴的数据显示范围
    %         set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    %         % set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    %         % set(gca,'XTickLabel',{'length1','length2','length3','length4'});%给坐标加标签
    %         set(gca,'XTickLabel',{'1','2','3','4'});%给坐标加标签
    %         %text(4.15,-0.14,sprintf('length'),'fontsize',9);
    %         ax = gca;
    %         ax.XAxis.FontSize = 9;
    %         ax.YAxis.FontSize = 11;
    %         xtickangle(0);
    %
    %         set(gca,'box','off');% 取消右、上边框
    %         xlabel('Length', 'FontSize', 11);
    %         % ylabel('Recall accuracy', 'FontSize', 11);
    %         if if_title == 1
    %             title(sprintf('Recall accuracy'), 'FontSize', 11);
    %         end
    %
    %
    %         %% Plot rProb
    %         %fig2 = figure(2);
    %         %set(gcf,'Position',[200 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %
    %         nexttile
    %
    %         errorbar(seq_length_rangeHead: seq_length_rangeTail,...
    %             c4_mean,c4_sem,...
    %             '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'CapSize', 6, 'MarkerSize', 1);
    %
    %         set(gca,'linewidth',1.5)
    %         ylim([0 1]);
    %         set(gca, 'FontSize', 11)
    %         %set(gca,'YTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    %         % set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    %         set(gca,'XLim',[seq_length_rangeHead-0.75 seq_length_rangeTail+0.75]);%X轴的数据显示范围
    %         set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    %         % set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    %         % set(gca,'XTickLabel',{'length1','length2','length3','length4'});%给坐标加标签
    %         set(gca,'XTickLabel',{'1','2','3','4'});%给坐标加标签
    %         %text(4.15,-0.14,sprintf('length'),'fontsize',9);
    %         ax = gca;
    %         ax.XAxis.FontSize = 9;
    %         ax.YAxis.FontSize = 11;
    %         xtickangle(0);
    %
    %         if if_plot_rProb_scale == 1
    %             y_max = max(c4_mean+c4_sem);
    %             y_min = min(c4_mean-c4_sem);
    %
    %             ylim([y_min-(y_max-y_min)*0.3 y_max+(y_max-y_min)*0.3]);
    %             set(gca,'YTick',[0:0.1:1]);
    %         else
    %             set(gca,'YTick',[0:0.2:1]);
    %         end
    %
    %         set(gca,'box','off');% 取消右、上边框
    %         xlabel('Length', 'FontSize', 11);
    %         % ylabel('Offloading rate', 'FontSize', 11);
    %         if if_title == 1
    %             title(sprintf('Offloading rate'), 'FontSize', 11);
    %         end
    %
    %         %% Plot correlation
    %         %fig3 = figure(3);
    %         %set(gcf,'Position',[400 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %
    %         nexttile
    %
    %         % isvalid = (~isnan(gError_noChoice_collapsed_inOne)) & (~isnan(offloadingProb_all_inOne));
    %         isvalid = (~isnan(gAcc_noChoice_collapsed_inOne)) & (~isnan(offloadingProb_all_inOne));
    %
    %         % x = gError_noChoice_collapsed_inOne(isvalid);
    %         x = gAcc_noChoice_collapsed_inOne(isvalid);
    %         y = offloadingProb_all_inOne(isvalid);
    %         n = 1;
    %         [p_mapping,S] = polyfit(x,y,n);
    %         r2 = 1 - (S.normr/norm(y - mean(y)))^2;
    %         [r, p_corr] = corr(x',y');
    %         tempTxt = 'p>0.05';
    %         if p_corr < 0.05
    %             tempTxt = 'p<0.05';
    %         end
    %         if p_corr < 0.01
    %             tempTxt = 'p<0.01';
    %         end
    %         if p_corr < 0.001
    %             tempTxt = 'p<0.001';
    %         end
    %
    %         [p_mapping_poly,S_poly] = polyfit(x,y,3);
    %         r2_poly = 1 - (S_poly.normr/norm(y - mean(y)))^2;
    %
    %
    %
    %
    %         x_fit = 0:0.1:1;
    %         y_fit = polyval(p_mapping,x_fit);
    %
    %         temp_color = [[166,97,26]/255;
    %             [223,194,125]/255;
    %             [128,205,193]/255;
    %             [1,133,113]/255];
    %
    %
    %         for target_seqLength=1:pointKindsNum
    %             %temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(gError_noChoice_collapsed{target_seqLength}));
    %             temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(gAcc_noChoice_collapsed{target_seqLength}));
    %
    %             temp_size = temp_size ./ 5;
    %
    %             %scatter(gError_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
    %             scatter(gAcc_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
    %                 temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
    %                 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    %             hold on
    %         end
    %         plot(x_fit, y_fit, '-', 'LineWidth', 1, 'Color', [0.35 0.35 0.35 0.7]);
    %
    %         % text(0.1,0.9,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',10,'FontWeight','bold');
    %         % text(0.2,0.1,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',10,'FontWeight','bold');
    %         text(0.66,0.99,sprintf('r=%.3f',r), 'fontsize',10,'FontWeight','bold');
    %         text(0.66,0.80,sprintf('%s',tempTxt), 'fontsize',10,'FontWeight','bold');
    %
    %         set(gca,'linewidth',1.5)
    %         % xlim([0 1]);
    %         x_min = 0;
    %         x_max = 1;
    %         xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
    %
    %         ylim([0 1]);
    %         set(gca, 'FontSize', 11)
    %         set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    %         %set(gca,'XTick',[0,0.5]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    %
    %         xtickangle(0);
    %
    %         % set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    %         set(gca,'YTick',[0:0.5:1]);%设置要显示坐标刻度的范围
    %         set(gca,'box','off');% 取消右、上边框
    %         % xlabel('Error rate', 'FontSize', 11);
    %         xlabel('Recall accuracy', 'FontSize', 11);
    %         %text(0.7,-0.17,sprintf('Recall accuracy'),'fontsize',9);
    %
    %         temp_ylabel = ylabel('Offloading rate', 'FontSize', 11);
    %         % temp_Position = temp_ylabel.Position;
    %         % temp_Position(1) = temp_Position(1) - 0.03 + 1;
    %         % temp_Position(2) = temp_Position(2) - 0.01 - 0.2;
    %         % temp_ylabel.Position = temp_Position;
    %         if if_title == 1
    %             title('Correlation', 'FontSize', 11);
    %         end
    
    %     end
    
    
    
end

if if_singlePanel0_multiPanel1 == 0
    %% Plot accuracy
    fig1 = figure(1);
    set(gcf,'Position',[0 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    
    % errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(internalNoChoiceAccuracy),std(internalNoChoiceAccuracy)./sqrt(size(internalNoChoiceAccuracy, 1)),'-o','Color',color_noChoice,'LineWidth',1,'CapSize',12,'MarkerSize',4);
    errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(offloadChoiceAccuracy),std(offloadChoiceAccuracy)./sqrt(size(offloadChoiceAccuracy, 1)),'-o','Color',color_choiceOffload,'LineWidth',1,'CapSize',6,'MarkerSize',1);
    hold on
    errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(internalChoiceAccuracy),std(internalChoiceAccuracy)./sqrt(size(internalChoiceAccuracy, 1)),'-o','Color',color_choiceMemory,'LineWidth',1,'CapSize',6,'MarkerSize',1);
    hold on
    errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(internalNoChoiceAccuracy),std(internalNoChoiceAccuracy)./sqrt(size(internalNoChoiceAccuracy, 1)),'-o','Color',color_noChoice,'LineWidth',1,'CapSize',6,'MarkerSize',1);
    hold on
    
    % errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(ProbOffloadingInSeqLength),std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)),'-o','Color',color_rProb,'LineWidth',1,'CapSize',7,'MarkerSize',2);
    % hold on
    
    % le = legend('ChoiceOffload',...
    %     'ChoiceMemory','ForcedToTest',...
    %     'Location','southwest','fontsize',9);
    % le.ItemTokenSize = ones(1,2)*10;
    
    
    
    set(gca,'linewidth',1.5)
    ylim([0 1]);
    set(gca, 'FontSize', 11)
    set(gca,'YTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    % set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    % set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    % set(gca,'XTickLabel',{'length1','length2','length3','length4'});%给坐标加标签
    set(gca,'XTickLabel',{'1','2','3','4'});%给坐标加标签
    text(4.15,-0.14,sprintf('length'),'fontsize',9);
    ax = gca;
    ax.XAxis.FontSize = 9;
    ax.YAxis.FontSize = 11;
    xtickangle(0);
    
    set(gca,'box','off');% 取消右、上边框
    % xlabel('Length', 'FontSize', 11);
    % ylabel('Recall accuracy', 'FontSize', 11);
    if if_title == 1
        title(sprintf('Recall accuracy'), 'FontSize', 11);
    end
    
    
    %% Plot rProb
    fig2 = figure(2);
    set(gcf,'Position',[200 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    errorbar(seq_length_rangeHead: seq_length_rangeTail,...
        mean(ProbOffloadingInSeqLength),std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)),...
        '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'CapSize', 6, 'MarkerSize', 1);
    
    set(gca,'linewidth',1.5)
    ylim([0 1]);
    set(gca, 'FontSize', 11)
    set(gca,'YTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    % set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    % set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    % set(gca,'XTickLabel',{'length1','length2','length3','length4'});%给坐标加标签
    set(gca,'XTickLabel',{'1','2','3','4'});%给坐标加标签
    text(4.15,-0.14,sprintf('length'),'fontsize',9);
    ax = gca;
    ax.XAxis.FontSize = 9;
    ax.YAxis.FontSize = 11;
    xtickangle(0);
    
    set(gca,'box','off');% 取消右、上边框
    % xlabel('Length', 'FontSize', 11);
    % ylabel('Offloading rate', 'FontSize', 11);
    if if_title == 1
        title(sprintf('Offloading rate'), 'FontSize', 11);
    end
    
    %% Plot correlation
    fig3 = figure(3);
    set(gcf,'Position',[400 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    % isvalid = (~isnan(gError_noChoice_collapsed_inOne)) & (~isnan(offloadingProb_all_inOne));
    isvalid = (~isnan(gAcc_noChoice_collapsed_inOne)) & (~isnan(offloadingProb_all_inOne));
    
    % x = gError_noChoice_collapsed_inOne(isvalid);
    x = gAcc_noChoice_collapsed_inOne(isvalid);
    y = offloadingProb_all_inOne(isvalid);
    n = 1;
    [p_mapping,S] = polyfit(x,y,n);
    r2 = 1 - (S.normr/norm(y - mean(y)))^2;
    [r, p_corr] = corr(x',y');
    tempTxt = 'p>0.05';
    if p_corr < 0.05
        tempTxt = 'p<0.05';
    end
    if p_corr < 0.01
        tempTxt = 'p<0.01';
    end
    if p_corr < 0.001
        tempTxt = 'p<0.001';
    end
    
    [p_mapping_poly,S_poly] = polyfit(x,y,3);
    r2_poly = 1 - (S_poly.normr/norm(y - mean(y)))^2;
    
    
    
    
    x_fit = 0:0.01:1;
    y_fit = polyval(p_mapping,x_fit);
    
    temp_color = [[166,97,26]/255;
        [223,194,125]/255;
        [128,205,193]/255;
        [1,133,113]/255];
    
    
    for target_seqLength=1:pointKindsNum
        %temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(gError_noChoice_collapsed{target_seqLength}));
        temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(gAcc_noChoice_collapsed{target_seqLength}));
        
        temp_size = temp_size ./ 5;
        
        %scatter(gError_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
        scatter(gAcc_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
            temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
    end
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    
    % text(0.1,0.9,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',10,'FontWeight','bold');
    % text(0.2,0.1,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',10,'FontWeight','bold');
    text(0.66,0.99,sprintf('r=%.3f',r), 'fontsize',10,'FontWeight','bold');
    text(0.66,0.80,sprintf('%s',tempTxt), 'fontsize',10,'FontWeight','bold');
    
    set(gca,'linewidth',1.5)
    % xlim([0 1]);
    x_min = 0;
    x_max = 1;
    xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
    
    ylim([0 1]);
    set(gca, 'FontSize', 11)
    % set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    set(gca,'XTick',[0,0.5]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    
    xtickangle(0);
    
    % set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    set(gca,'YTick',[0:0.5:1]);%设置要显示坐标刻度的范围
    set(gca,'box','off');% 取消右、上边框
    % xlabel('Error rate', 'FontSize', 11);
    % xlabel('Recall accuracy', 'FontSize', 11);
    text(0.7,-0.17,sprintf('Recall accuracy'),'fontsize',9);
    
    % temp_ylabel = ylabel('Offloading rate', 'FontSize', 11);
    % temp_Position = temp_ylabel.Position;
    % temp_Position(1) = temp_Position(1) - 0.03 + 1;
    % temp_Position(2) = temp_Position(2) - 0.01 - 0.2;
    % temp_ylabel.Position = temp_Position;
    if if_title == 1
        %     title('Correlation', 'FontSize', 11);
    end
    
end


%% Plot correlationB
if if_plot_precision == 1
    if true
        fig41 = figure(41);
        %set(gcf,'Position',[800 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[800 50 figure_xPixel*1.2 figure_yPixel*1.56]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[800 50 226*0.75 186*0.9*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
        nexttile
        
        x_raw = seqPrecision_behavior';
        %y_raw = offloadingProb_all_inOne;
        y_raw = gAcc_noChoice_collapsed_inOne;
        
        
        isvalid = (~isnan(x_raw)) & (~isnan(y_raw));
        
        x = x_raw(isvalid);
        y = y_raw(isvalid);
        
        [x_min,x_max] = bounds(x);
        
        n = 1;
        [p_mapping,S] = polyfit(x,y,n);
        r2 = 1 - (S.normr/norm(y - mean(y)))^2;
        [r, p_corr] = corr(x',y');
        tempTxt = 'p > 0.05';
        if p_corr < 0.05
            tempTxt = 'p < 0.05';
        end
        if p_corr < 0.01
            tempTxt = 'p < 0.01';
        end
        if p_corr < 0.001
            tempTxt = 'p < 0.001';
        end
        
        [p_mapping_poly,S_poly] = polyfit(x,y,3);
        r2_poly = 1 - (S_poly.normr/norm(y - mean(y)))^2;
        
        
        if if_entropy_behaviorInverted == 1
            x_min = 0.3;
            x_max = 1;
        elseif if_entropy_behaviorInverted == 0
            x_min = 0;
            x_max = 0.7;
        end        
        
        
        x_fit = x_min:0.01:x_max;
        y_fit = polyval(p_mapping,x_fit);
        
        temp_color = [[166,97,26]/255;
            [223,194,125]/255;
            [128,205,193]/255;
            [1,133,113]/255];
        
        
        for target_seqLength=1:pointKindsNum
            temp_range = (sum(numSeq(1:target_seqLength-1))+1):sum(numSeq(1:target_seqLength));
            temp_size = (((target_seqLength.^3)*2 + 3) .* ones(1, numSeq(target_seqLength))) ./ 5;
            
            scatter(x(temp_range),y(temp_range), ...
                temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
        end
        plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        
        %x_min = 0.4;
        
        if if_entropy_behaviorInverted == 1
            text(x_min+(x_max-x_min)*0.63,0.22,sprintf('r = %.3f',r), 'fontsize',9,'FontWeight','normal');
            text(x_min+(x_max-x_min)*0.63,0.10,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
        elseif if_entropy_behaviorInverted == 0
            text(x_min+(x_max-x_min)*0.0,0.22,sprintf('r = %.3f',r), 'fontsize',9,'FontWeight','normal');
            text(x_min+(x_max-x_min)*0.0,0.10,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
        end
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
        
        %x_min = 0.4;
        %x_max = 1;
        xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05]);
        temp_y_min = 0;
        temp_y_max = 1;
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.1]);%0.3
        
        %ylim([0 1]);
        set(gca, 'FontSize', 8)
        %set(gca,'XTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        %set(gca,'XTick',[1,1.6]);
        
        xtickangle(0);
        
        %set(gca,'XTick',[0.4,1]);%设置要显示坐标刻度的范围
        set(gca,'XTick',[x_min,x_max]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        % set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        set(gca,'YTick',[0:1]);%设置要显示坐标刻度的范围
        set(gca,'box','off');% 取消右、上边框
        % xlabel('Error rate', 'FontSize', 11);
        %xlabel('Recall variability', 'FontSize', 9);
        text((x_min+x_max)/2,-0.205,sprintf('Recall variability'),'fontsize',9,'HorizontalAlignment','center');
        %text(1.90,-0.17,sprintf('Precision'),'fontsize',9);
        
        temp_ylabel = ylabel('Recall accuracy', 'FontSize', 9);
        % temp_Position = temp_ylabel.Position;
        % temp_Position(1) = temp_Position(1) - 0.03;
        % temp_Position(2) = temp_Position(2) - 0.01;
        % temp_ylabel.Position = temp_Position;
        if if_title == 1
            %title('Precision (1/sigma)', 'FontSize', 11);
            title('Behavioral correlation', 'FontSize', 9);
        end
        
    end
    
    
    %if true
    if false
        fig42 = figure(42);
        %set(gcf,'Position',[800 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[1050 50 figure_xPixel*1.2 figure_yPixel*1.56]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[1050 50 226*0.75 186]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
        %nexttile
        
        %x_raw = gAcc_noChoice_collapsed_inOne;
        x_raw = seqPrecision_behavior';
        y_raw = offloadingProb_all_inOne;
        
        
        isvalid = (~isnan(x_raw)) & (~isnan(y_raw));
        
        x = x_raw(isvalid);
        y = y_raw(isvalid);
        
        [x_min,x_max] = bounds(x);
        
        n = 1;
        [p_mapping,S] = polyfit(x,y,n);
        r2 = 1 - (S.normr/norm(y - mean(y)))^2;
        [r, p_corr] = corr(x',y');
        tempTxt = 'p>0.05';
        if p_corr < 0.05
            tempTxt = 'p<0.05';
        end
        if p_corr < 0.01
            tempTxt = 'p<0.01';
        end
        if p_corr < 0.001
            tempTxt = 'p<0.001';
        end
        
        [p_mapping_poly,S_poly] = polyfit(x,y,3);
        r2_poly = 1 - (S_poly.normr/norm(y - mean(y)))^2;
        
        
        
        
        x_fit = x_min:0.01:x_max;
        y_fit = polyval(p_mapping,x_fit);
        
        temp_color = [[166,97,26]/255;
            [223,194,125]/255;
            [128,205,193]/255;
            [1,133,113]/255];
        
        
        for target_seqLength=1:pointKindsNum
            temp_range = (sum(numSeq(1:target_seqLength-1))+1):sum(numSeq(1:target_seqLength));
            temp_size = (((target_seqLength.^3)*2 + 3) .* ones(1, numSeq(target_seqLength))) ./ 5;
            
            scatter(x(temp_range),y(temp_range), ...
                temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
        end
        plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        
        text(x_min+(x_max-x_min)*0.00,0.23,sprintf('r=%.3f',r), 'fontsize',10,'FontWeight','normal');
        text(x_min+(x_max-x_min)*0.00,0.10,sprintf('%s',tempTxt), 'fontsize',10,'FontWeight','normal');
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
        ylim([0 1]);
        set(gca, 'FontSize', 8)
        %set(gca,'XTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        %set(gca,'XTick',[1,1.6]);
        
        xtickangle(0);
        
        % set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        set(gca,'YTick',[0:1]);%设置要显示坐标刻度的范围
        set(gca,'box','off');% 取消右、上边框
        % xlabel('Error rate', 'FontSize', 11);
        xlabel('Recall variability', 'FontSize', 9);
        %text(1.90,-0.17,sprintf('Precision'),'fontsize',9);
        
        temp_ylabel = ylabel('Offloading rate', 'FontSize', 9);
        % temp_Position = temp_ylabel.Position;
        % temp_Position(1) = temp_Position(1) - 0.03;
        % temp_Position(2) = temp_Position(2) - 0.01;
        % temp_ylabel.Position = temp_Position;
        if if_title == 1
            %title('Precision (1/sigma)', 'FontSize', 11);
            title('Behavioral correlation', 'FontSize', 9);
        end
        
    end
    
end

%% Plot exampleSeqProb
if if_plot_precision == 1 && if_entropy == 0
    fig5 = figure(5);
    %set(gcf,'Position',[800 300 figure_xPixel*2.7 figure_yPixel*2*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[800 300 figure_xPixel*1.5*0.875 figure_yPixel*1.56]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    
    tempSeq = seqSet_inOne_inOne(exampleSeq);
    %t.Title.String = sprintf('Example sequence (%d)',tempSeq);
    t.Title.String = sprintf('Example stimulus sequence: %d',tempSeq);
    t.Title.FontSize = 11;
    t.Title.Interpreter = 'none';
    t.Title.FontWeight = 'normal';%bold
    %t.Position(1) = t.Position(1) - 0.01;
    
    nexttile
    
    %% exampleSeq probability distribution
    seqProbIndex_exampleSeq_sortedTopXB;
    
    x = 1:(length(seqProbIndex_exampleSeq_sortedTopXB)+1);
    y = [seqProb_exampleSeq_sortedTopXB;0];
    
    [x_min,x_max] = bounds(x);
    
    scatter(x,y,20, 'filled', 'MarkerFaceColor', [1 1 1]*0.15, ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    set(gca,'linewidth',1.5)
    xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05]);
    ylim([0 1]);
    set(gca, 'FontSize', 10)
    
    xtickangle(0);
    %temp_labels = [string(seqSet_inOne_inOne(seqProbIndex_exampleSeq_sortedTopX));'others'];
    %xticklabels(temp_labels);
    xticklabels('');
    
    set(gca,'XTick',x);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    set(gca,'YTick',[0,1]);%设置要显示坐标刻度的范围
    set(gca,'box','off');% 取消右、上边框
    %xlabel('Sequence', 'FontSize', 11);
    ylabel('Probability', 'FontSize', 10);
    
    
    %% exampleSeq score distribution
    nexttile
    
    seqProbIndex_exampleSeq_sortedTopXB;
    
    x_raw = seqProbIndex_exampleSeq_sortedTopXB;
    x = 1:(length(seqProbIndex_exampleSeq_sortedTopXB));
    y = score_exampleSeq(seqProbIndex_exampleSeq_sortedTopXB);
    
    [x_min,x_max] = bounds(x);
    [y_min,y_max] = bounds(y);
    
    scatter(x,y,20, 'filled', 'MarkerFaceColor', [1 1 1]*0.15, ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    set(gca,'linewidth',1.5)
    xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05+1]);
    ylim([y_min-(y_max-y_min)*0.05 y_max+(y_max-y_min)*0.05]);
    set(gca, 'FontSize', 10)
    
    xtickangle(90);%20
    temp_labels = [string(seqSet_inOne_inOne(x_raw));'others'];
    xticklabels(temp_labels);
    
    set(gca,'XTick',[x,max(x)+1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    %set(gca,'YTick',[0:0.5:1]);%设置要显示坐标刻度的范围
    set(gca,'box','off');% 取消右、上边框
    %xlabel('Sequence', 'FontSize', 11);
    xlabel('Response sequence', 'FontSize', 10);
    ylabel('Similarity', 'FontSize', 10);
    
    a = 1;
    
    fig6 = figure(6);
    %set(gcf,'Position',[800 700 figure_xPixel*1.9 figure_yPixel*1.15]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[800 700 figure_xPixel*1 figure_yPixel*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[800 600 figure_xPixel*1.2 figure_yPixel*1.56]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
    
    tempSeq = seqSet_inOne_inOne(exampleSeq);
    %t.Title.String = sprintf('Example sequence (%d)',tempSeq);
    %t.Title.String = sprintf('precision (1/sigma) = %.2f',1/sigma_exampleSeq);
    t.Title.String = sprintf('Precision = 1/sigma');
    %t.Title.String = sprintf('Precision = 1/sigma^2');
    t.Title.FontSize = 11;
    t.Title.Interpreter = 'none';
    t.Title.FontWeight = 'bold';
    
    nexttile
    
    
    %% seqProb_exampleSeq (each score bin)
    
    %score_exampleSeq_unique;
    %seqProb_exampleSeq_uniqueSum;
    %score_exampleSeq_unique_fit;
    %seqProb_exampleSeq_uniqueSum_fit;
    
    if if_precision_meanProb0_sumProb1 == 0
        x = score_exampleSeq;
        y = seqProb_exampleSeq;
        x_fit = score_exampleSeq_fit;
        y_fit = seqProb_exampleSeq_fit;
    elseif if_precision_meanProb0_sumProb1 == 1
        x = score_exampleSeq_unique;
        y = seqProb_exampleSeq_uniqueSum;
        x_fit = score_exampleSeq_unique_fit;
        y_fit = seqProb_exampleSeq_uniqueSum_fit;
    end
    
    
    %x_fit = x_fit(x_fit>=x_min);
    
    [x_min,x_max] = bounds(x);
    
    scatter(x,y,20, 'filled', 'MarkerFaceColor', [1 1 1]*0.15, ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    plot(x_fit,y_fit,'Color',[0.5 0.5 0.5],'LineWidth',1.5);
    hold on
    
    set(gca,'linewidth',1.5)
    xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05]);
    ylim([0 1]);
    set(gca, 'FontSize', 11)
    
    %text(5,0.5,sprintf('sigma=%.2f',sigma_exampleSeq),'fontsize',11);
    %text(5,0.5,sprintf('precision (1/sigma) = %.2f',1/sigma_exampleSeq),'fontsize',11);
    
    xtickangle(0);
    set(gca,'XDir','reverse');
    % set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    set(gca,'YTick',[0,1]);%设置要显示坐标刻度的范围
    set(gca,'box','off');% 取消右、上边框
    %xlabel('Sequence similarity', 'FontSize', 11);
    xlabel('Similarity', 'FontSize', 11);
    
    temp_ylabel = ylabel('Probability', 'FontSize', 11);
    %     if if_title == 1
    %         %tempBoolIndex_seq = boolIndex_location_seq_T(exampleSeq,:);
    %         %tempSeq = target_seqSet_inOne{exampleSeq};
    %         tempSeq = seqSet_inOne_inOne(exampleSeq);
    %         title(sprintf('Example sequence probability (%d)',tempSeq), 'FontSize', 11);
    %     end
    
    
    
    
    %% seqProb_exampleSeq (each score bin) no axis
    
    fig7 = figure(7);
    set(gcf,'Position',[500 700 figure_xPixel*1.5 figure_yPixel*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
    
    nexttile
    
    if if_precision_meanProb0_sumProb1 == 0
        x = score_exampleSeq;
        y = seqProb_exampleSeq;
        x_fit = score_exampleSeq_fit;
        y_fit = seqProb_exampleSeq_fit;
    elseif if_precision_meanProb0_sumProb1 == 1
        x = score_exampleSeq_unique;
        y = seqProb_exampleSeq_uniqueSum;
        x_fit = score_exampleSeq_unique_fit;
        y_fit = seqProb_exampleSeq_uniqueSum_fit;
    end
    
    [x_min,x_max] = bounds(x);
    
    scatter(x,y,20, 'filled', 'MarkerFaceColor', [1 1 1]*0.15, ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    plot(x_fit,y_fit,'Color',[1 1 1]*0,'LineWidth',1.5);
    hold on
    
    set(gca,'linewidth',1.5)
    xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05]);
    ylim([0 1]);
    set(gca,'XDir','reverse');
    set(gca, 'visible', 'off');
    
    
    %% spatialInterferenceWeight
    fig = figure;
    %set(gcf,'Position',[1040 50 figure_xPixel*1.2 figure_yPixel*1.56]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[1300 50 figure_xPixel*1.4 figure_yPixel*1.56]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
    nexttile
    
    x = temp_spatialInterferenceWeight';
    y = temp_r_mean_iter;
    
    [x_min,x_max] = bounds(x);
    [y_min,y_max] = bounds(y);
    
    plot(x, y, '-', 'LineWidth', 1, 'Color', [0.15 0.15 0.15]);
    hold on
    plot([1 1].*spatialInterferenceWeight, [y_min-(y_max-y_min)*0.1 y_max],...
        '--', 'LineWidth', 1, 'Color', [0.55 0.55 0.55]);
    hold on
    plot([x_min-(x_max-x_min)*0.1 spatialInterferenceWeight], [1 1].*temp_r_max_iter,...
        '--', 'LineWidth', 1, 'Color', [0.55 0.55 0.55]);
    hold on
    
    
    text(spatialInterferenceWeight,temp_r_max_iter+(y_max-y_min)*0.125,sprintf('optimal'),'fontsize',11,'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
    ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.3]);
    set(gca, 'FontSize', 11)
    
    xtickangle(0);
    
    %set(gca,'YTick',[0:0.5:1]);%设置要显示坐标刻度的范围
    set(gca,'box','off');% 取消右、上边框
    xlabel('Spatial Interference Weight', 'FontSize', 11);
    
    temp_ylabel = ylabel('r', 'FontSize', 11);
    
    if if_title == 1
        title(sprintf('Correlation between\nProbability and Similarity'), 'FontSize', 11);
    end
end



if if_plot_submissionRT == 1
    endingHold.noChoiceMemory.all_mean_inOne;
    endingHold.noChoiceMemory.correct;
    endingHold.noChoiceMemory.error;
    
    endingHold_noCoiceMemory_seq = endingHold.noChoiceMemory.all_mean_inOne;
    
    endingHold_noCoiceMemoryCorrect_trial = [];
    endingHold_noCoiceMemoryError_trial = [];
    for tempi=1:length(endingHold.noChoiceMemory.correct)
        tempA_correct = endingHold.noChoiceMemory.correct{tempi};
        tempA_error = endingHold.noChoiceMemory.error{tempi};
        
        for tempj=1:length(tempA_correct)
            tempB_correct = tempA_correct{tempj};
            tempB_error = tempA_error{tempj};
            
            endingHold_noCoiceMemoryCorrect_trial = [endingHold_noCoiceMemoryCorrect_trial tempB_correct];
            endingHold_noCoiceMemoryError_trial = [endingHold_noCoiceMemoryError_trial tempB_error];
        end
    end
    
    endingHold_noCoiceMemory_seq;
    endingHold_noCoiceMemoryCorrect_trial;
    endingHold_noCoiceMemoryError_trial;
    
    lowThreshold_endingHold = median([endingHold_noCoiceMemoryCorrect_trial endingHold_noCoiceMemoryError_trial]);
    
    endingHold;
    
    seqAccuracy_endingHoldLow_noCoice = nan(sum(numSeq(1:valid_length)),1);
    seqAccuracy_endingHoldHigh_noCoice = nan(sum(numSeq(1:valid_length)),1);
    
    for tempi=1:valid_length
        tempA_correct = endingHold.noChoiceMemory.correct{tempi};
        tempA_error = endingHold.noChoiceMemory.error{tempi};
        
        
        for tempj=1:length(tempA_correct)
            tempB_correct = tempA_correct{tempj};
            tempB_error = tempA_error{tempj};
            
            tempTrialNum_correct_low = sum(tempB_correct<lowThreshold_endingHold);
            tempTrialNum_correct_high = sum(tempB_correct>=lowThreshold_endingHold);
            
            tempTrialNum_error_low = sum(tempB_error<lowThreshold_endingHold);
            tempTrialNum_error_high = sum(tempB_error>=lowThreshold_endingHold);
            
            tempAcc_low = tempTrialNum_correct_low./(tempTrialNum_correct_low+tempTrialNum_error_low);
            tempAcc_high = tempTrialNum_correct_high./(tempTrialNum_correct_high+tempTrialNum_error_high);
            
            temp1 = sum(numSeq(1:tempi-1))+tempj;
            
            seqAccuracy_endingHoldLow_noCoice(temp1) = tempAcc_low;
            seqAccuracy_endingHoldHigh_noCoice(temp1) = tempAcc_high;
            
        end
        
    end
    
    
    
    b1 = seqAccuracy_endingHoldLow_noCoice;
    b2 = seqAccuracy_endingHoldHigh_noCoice;
    
    b1_mean = mean(b1,'omitnan');
    b2_mean = mean(b2,'omitnan');
    
    [~,temp_p_lowHighEndingHold] = ttest(b1,b2);
    
    
    if true
        endingHold_allSeq_noChoiceMemoryCorrectMean = nan(sum(numSeq(1:valid_length)),1);
        endingHold_allSeq_noChoiceMemoryErrorMean = nan(sum(numSeq(1:valid_length)),1);
        
        for tempi=1:valid_length
            tempA_correct = endingHold.noChoiceMemory.correct{tempi};
            tempA_error = endingHold.noChoiceMemory.error{tempi};
            
            for tempj=1:length(tempA_correct)
                tempB_correct = tempA_correct{tempj};
                tempB_error = tempA_error{tempj};
                
                temp1 = sum(numSeq(1:tempi-1))+tempj;
                
                endingHold_allSeq_noChoiceMemoryCorrectMean(temp1) = mean(tempB_correct);
                endingHold_allSeq_noChoiceMemoryErrorMean(temp1) = mean(tempB_error);
            end
        end
        d1 = endingHold_allSeq_noChoiceMemoryCorrectMean;
        d2 = endingHold_allSeq_noChoiceMemoryErrorMean;
        [~,temp_p_correctErrorEndingHold] = ttest(d1,d2);
    end
    
    
    
    
    temp_1 = endingHold_noCoiceMemoryError_trial;
    temp_2 = endingHold_noCoiceMemoryCorrect_trial;
    [~,temp_p12] = ttest2(temp_1,temp_2);
    
    %% Plot
    fig = figure('Name','locDistri','NumberTitle','off');
    %set(gcf,'Position',[50+0 50+0 240*0.80*2 336*1.11*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 350+0 240*0.80*3 336*1.11*0.9*0.6]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 350+0 240*0.80*3 336*1.11*0.9*0.6*2*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 350+0 240*0.80*3*1.3 336*1.11*0.9*0.6*2*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[50+0 350+0 240*0.80*3*1.3 336*1.11*0.9*0.6*1.3]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,4,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    
    
    temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM]);
    temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
    
    temp_bar = bar([1 2], [mean(temp_1) mean(temp_2)], ...
        'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1);
    hold on
    a = 1;
    temp_bar.CData(1,:) = [1 1 1]*0.5;%0.3
    temp_bar.CData(2,:) = [1 1 1]*0.5;
    
    errorbar([1 2], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 8);
    
    %
    %     temp_y_min = min([temp_1 temp_2]);
    %     temp_y_max = max([temp_1 temp_2]);
    %
    %     temp_data = [temp_1';temp_2'];
    %
    %     g1 = repmat({'A'},length(temp_1),1);
    %     g2 = repmat({'B'},length(temp_2),1);
    %
    %     temp_label = [g1;g2];
    %
    %     temptemp_color1 = [1 1 1]*0.5;%0.5
    %     temptemp_color2 = repmat(temptemp_color1, 2, 1);
    %
    %     h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
    %         'GroupOrder',[{'A'};{'B'}]);
    %     h(1).ViolinPlot.FaceAlpha = 0.1;
    %     h(2).ViolinPlot.FaceAlpha = 0.1;
    
    
    tempTxt = sprintf('');
    if temp_p12 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p12 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p12 < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.06,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([1.1 1.9],temp_y_max+(temp_y_max-temp_y_min)*0.05*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    
    set(gca,'linewidth',1.5)
    
    xlim([0.35 2.65]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.25]);%0.3
    set(gca, 'FontSize', 10) %14
    %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = ["Error", "Correct"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    if if_monkey_D0_Z1 == 0
        %xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*1.30;
        xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*1.19;
    elseif if_monkey_D0_Z1 == 1
        xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*1.68;
    end
    
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',9);%9-->12
    
    set(gca,'xticklabel','');
    
    %xlabel(sprintf('Submission RT'), 'FontSize', 11, 'FontWeight', 'normal');
    ylabel('Reaction time (s)', 'FontSize', 10, 'FontWeight', 'normal');
    temp_title = title(sprintf('Submission RT'),'fontsize',10);
    temp_title.Interpreter = 'none';
    
    
    nexttile
    
    temp_x = endingHold_noCoiceMemory_seq';
    temp_y = offloadingProb_all_inOne';
    tempBoolIndex = ~isnan(temp_x);
    
    x = temp_x(~isnan(temp_x));
    y = temp_y(1:sum(numSeq(1:valid_length)));
    y = y(~isnan(temp_x));
    
    [r,p] = corr(x,y);
    
    mdl = fitglm(x,y);
    
    for tempi=1:valid_length
        temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
        tempIndex = find(tempBoolIndex(temp_range)==true);
        temp_range2 = temp_range(tempIndex);
        
        temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
        
        temp_size = temp_size ./ 5;
        
        temp_h = scatter(temp_x(temp_range2), temp_y(temp_range2), ...
            temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
    end
    
    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y);
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    
    tempTxt = 'p>0.05';
    if p < 0.05
        tempTxt = 'p<0.05';
    end
    if p < 0.01
        tempTxt = 'p<0.01';
    end
    if p < 0.001
        tempTxt = 'p<0.001';
    end
    
    text(temp_xmin+(temp_xmax-temp_xmin)*0.66,temp_ymin+(temp_ymax-temp_ymin)*0.2,sprintf('r=%.3f',r), 'fontsize',10,'FontWeight','bold');
    text(temp_xmin+(temp_xmax-temp_xmin)*0.66,temp_ymin+(temp_ymax-temp_ymin)*0.05,sprintf('%s',tempTxt), 'fontsize',10,'FontWeight','bold');
    
    
    xlim([temp_xmin-(temp_xmax-temp_xmin)*0.1 temp_xmax+(temp_xmax-temp_xmin)*0.1]);
    ylim([temp_ymin-(temp_ymax-temp_ymin)*0.1 temp_ymax+(temp_ymax-temp_ymin)*0.1]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 11)
    set(gca,'box','off');% 取消右、上边框
    %set(gca,'XTick',0:0.5:1,'FontSize', 12);%给坐标加标签
    set(gca,'YTick',0:1,'FontSize', 11);%给坐标加标签
    
    xlabel(sprintf('Submission RT'), 'FontSize', 11, 'FontWeight', 'normal');
    %temp_title = title(sprintf('r=%.3f, p=%.3f',r,p), 'FontSize', 11, 'FontWeight', 'bold');
    %temp_title.Interpreter = 'none';
    
    ylabel(sprintf('Offloading rate'), 'FontSize', 11, 'FontWeight', 'normal');
    
    
    %% Low RT VS. high RT
    nexttile
    %b1 = seqAccuracy_endingHoldLow_noCoice;
    %b2 = seqAccuracy_endingHoldHigh_noCoice;
    %[~,temp_p_lowHighEndingHold] = ttest(b1,b2);
    
    temp_p = temp_p_lowHighEndingHold;
    temptempBoolIndex = (~isnan(b1)) & (~isnan(b2));
    temp_1 = b1(temptempBoolIndex);
    temp_2 = b2(temptempBoolIndex);
    
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    
    temp_y_min = min([temp_1;temp_2]);
    temp_y_max = max([temp_1;temp_2]);
    
    % violin plot
    %     temp_data = [temp_1;temp_2];
    %
    %     g1 = repmat({'A'},length(temp_1),1);
    %     g2 = repmat({'B'},length(temp_2),1);
    %
    %     temp_label = [g1;g2];
    %
    %     temptemp_color1 = [1 1 1]*0.5;
    %     temptemp_color2 = repmat(temptemp_color1, 2, 1);
    %
    %     h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
    %         'GroupOrder',[{'A'};{'B'}]);
    %     h(1).ViolinPlot.FaceAlpha = 0.1;
    %     h(2).ViolinPlot.FaceAlpha = 0.1;
    
    % paired line plot
    plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
    hold on
    plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
    hold on
    
    for tempi=1:length(temp_1)
        plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
        hold on
    end
    
    scatter(1*ones(1,length(temp_1)),temp_1,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
        'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
    hold on
    scatter(2*ones(1,length(temp_2)),temp_2,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
        'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 2.65])
    %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 10)
    
    xticks([1 2]);
    
    xtl = ["Low-RT"; "High-RT"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.435;%0.56,0.4
    %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.375;%0.56,0.4
    
    if if_monkey_D0_Z1 == 0
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.375;%0.56,0.4
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.265;%0.56,0.4
    elseif if_monkey_D0_Z1 == 1
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.285;%0.56,0.4
    end
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%25
    set(gca,'xticklabel','');
    
    
    %yticks([0.7 0.8]);
    
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Recall accuracy', 'FontSize', 11, 'FontWeight', 'bold');
    temp1 = title(sprintf('Seq accuracy'),'fontsize',10);
    
    
    %% Correct VS. Error
    nexttile
    %b1 = seqAccuracy_endingHoldLow_noCoice;
    %b2 = seqAccuracy_endingHoldHigh_noCoice;
    %[~,temp_p_lowHighEndingHold] = ttest(b1,b2);
    
    temp_p = temp_p_correctErrorEndingHold;
    temptempBoolIndex = (~isnan(d1)) & (~isnan(d2));
    temp_1 = d2(temptempBoolIndex);
    temp_2 = d1(temptempBoolIndex);
    
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    
    temp_y_min = min([temp_1;temp_2]);
    temp_y_max = max([temp_1;temp_2]);
    
    % violin plot
    %     temp_data = [temp_1;temp_2];
    %
    %     g1 = repmat({'A'},length(temp_1),1);
    %     g2 = repmat({'B'},length(temp_2),1);
    %
    %     temp_label = [g1;g2];
    %
    %     temptemp_color1 = [1 1 1]*0.5;
    %     temptemp_color2 = repmat(temptemp_color1, 2, 1);
    %
    %     h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
    %         'GroupOrder',[{'A'};{'B'}]);
    %     h(1).ViolinPlot.FaceAlpha = 0.1;
    %     h(2).ViolinPlot.FaceAlpha = 0.1;
    
    % paired line plot
    plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
    hold on
    plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
    hold on
    
    for tempi=1:length(temp_1)
        plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
        hold on
    end
    
    scatter(1*ones(1,length(temp_1)),temp_1,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
        'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
    hold on
    scatter(2*ones(1,length(temp_2)),temp_2,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
        'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.35 2.65])
    %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 10)
    
    xticks([1 2]);
    
    xtl = ["Error"; "Correct"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.435;%0.56,0.4
    %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.375;%0.56,0.4
    
    if if_monkey_D0_Z1 == 0
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.375;%0.56,0.4
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.125;%0.56,0.4
    elseif if_monkey_D0_Z1 == 1
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.285;%0.56,0.4
    end
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',9);%25
    set(gca,'xticklabel','');
    
    
    %yticks([0.7 0.8]);
    
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Reaction time (s)', 'FontSize', 10, 'FontWeight', 'normal');
    temp1 = title(sprintf('All Seqs'),'fontsize',10);
    
    
    %% Trial-based length level
    %if true
    if if_plot_submissionRT == 1
        %fig1 = figure(1);
        fig = figure('Name','Trial-based length level','NumberTitle','off');
        %set(gcf,'Position',[0 50 figure_xPixel*4 figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[0 50 figure_xPixel*4*1.1 figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[0 50 figure_xPixel*4*1.1 figure_yPixel*0.9*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[0 50 figure_xPixel*4*1.1*0.8*0.97 figure_yPixel*0.9*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[0 50 figure_xPixel*4*1.1*0.8 figure_yPixel*0.9*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        %t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
        %t = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
        t = tiledlayout(1,4,'TileSpacing','tight','Padding','compact');
        
        %nexttile
        %set(gca, 'visible', 'off')
        
        nexttile
        %% Plot accuracy
        
        % errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(internalNoChoiceAccuracy),std(internalNoChoiceAccuracy)./sqrt(size(internalNoChoiceAccuracy, 1)),'-o','Color',color_noChoice,'LineWidth',1,'CapSize',12,'MarkerSize',4);
        errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(offloadChoiceAccuracy),std(offloadChoiceAccuracy)./sqrt(size(offloadChoiceAccuracy, 1)),'-o','Color',color_choiceOffload,'LineWidth',1,'CapSize',6,'MarkerSize',1);
        hold on
        errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(internalChoiceAccuracy),std(internalChoiceAccuracy)./sqrt(size(internalChoiceAccuracy, 1)),'-o','Color',color_choiceMemory,'LineWidth',1,'CapSize',6,'MarkerSize',1);
        hold on
        errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(internalNoChoiceAccuracy),std(internalNoChoiceAccuracy)./sqrt(size(internalNoChoiceAccuracy, 1)),'-o','Color',color_noChoice,'LineWidth',1,'CapSize',6,'MarkerSize',1);
        hold on
        
        % errorbar(seq_length_rangeHead:seq_length_rangeTail,mean(ProbOffloadingInSeqLength),std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)),'-o','Color',color_rProb,'LineWidth',1,'CapSize',7,'MarkerSize',2);
        % hold on
        
        if false
            [~,temp_p1] = ttest(internalChoiceAccuracy,internalNoChoiceAccuracy);  %#ok<*UNRCH>
            [~,temp_p2] = ttest2(internalChoiceAccuracy,internalNoChoiceAccuracy);  %#ok<*UNRCH>
        end
        
        %le = legend('Choice-offload',...
        %    'Choice-memory','Forced-to-test',...
        %    'Location','southwest','fontsize',7);%9
        %le.ItemTokenSize = ones(1,2)*10;
        %
        %%le.Location = 'westoutside';
        %le.Position(1) = le.Position(1) - 0.2;
        %le.Position(2) = le.Position(2) + 0.14;
        
        set(gca,'linewidth',1.5)
        ylim([0 1]);
        set(gca, 'FontSize', 10)
        %set(gca,'YTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        set(gca,'YTick',[0:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        % set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XLim',[seq_length_rangeHead-0.75 seq_length_rangeTail+0.75]);%X轴的数据显示范围
        set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
        % set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
        % set(gca,'XTickLabel',{'length1','length2','length3','length4'});%给坐标加标签
        %set(gca,'XTickLabel',{'1','2','3','4'});%给坐标加标签
        %text(4.15,-0.14,sprintf('length'),'fontsize',9);
        ax = gca;
        ax.XAxis.FontSize = 9;
        ax.YAxis.FontSize = 10;%11
        
        %set(gca,'XTickLabel',{'Len1','Len2','Len3','Len4'});%给坐标加标签
        %xtickangle(0);
        %xtickangle(25);
        
        temp_y_min = 0;
        temp_y_max = 1;
        
        xtl = ["Len1","Len2","Len3","Len4"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.18;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%9
        set(gca,'xticklabel','');
        
        
        set(gca,'box','off');% 取消右、上边框
        %xlabel('Length', 'FontSize', 10);
        % ylabel('Recall accuracy', 'FontSize', 11);
        if if_title == 1
            title(sprintf('Recall accuracy'), 'FontSize', 10);
        end
        
        
        %% Plot rProb
        %fig2 = figure(2);
        %set(gcf,'Position',[200 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        nexttile
        
        x = [];
        y = [];
        for tempi=1:4
            ProbOffloadingInSeqLength;
            if if_plot_rProb_fit == 1
                temp_size = 10;
                scatter(tempi,ProbOffloadingInSeqLength(:,tempi),...
                    temp_size,'filled','MarkerFaceColor',[1 1 1]*0.5,...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
                hold on
            end
            x = [x ones(1,size(ProbOffloadingInSeqLength,1))*tempi]; %#ok<*AGROW>
            y = [y ProbOffloadingInSeqLength(:,tempi)'];
        end
        if if_plot_rProb_fit == 1
            %y = a1;
            [r, p_corr] = corr(x',y'); %#ok<*ASGLU>
            
            n = 1;
            [p_mapping,S] = polyfit(x,y,n);
            r2 = 1 - (S.normr/norm(y - mean(y)))^2;
            
            temp_mdl = fitglm(x(~isnan(x)),y(~isnan(x)));
            r2 = temp_mdl.Rsquared.Adjusted;
            beta = temp_mdl.Coefficients.Estimate(2);
            p = temp_mdl.Coefficients.pValue(2);
            
            
            x_fit = 1:1:4;
            y_fit = polyval(p_mapping,x_fit);
            
            plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        end
        
        if if_plot_rProb_fit == 0
            errorbar(seq_length_rangeHead: seq_length_rangeTail,...
                mean(ProbOffloadingInSeqLength),std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)),...
                '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'CapSize', 6, 'MarkerSize', 1);
        end
        
        set(gca,'linewidth',1.5)
        ylim([0 1]);
        set(gca, 'FontSize', 10)
        %set(gca,'YTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        % set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XLim',[seq_length_rangeHead-0.75 seq_length_rangeTail+0.75]);%X轴的数据显示范围
        set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
        % set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
        % set(gca,'XTickLabel',{'length1','length2','length3','length4'});%给坐标加标签
        %set(gca,'XTickLabel',{'1','2','3','4'});%给坐标加标签
        %text(4.15,-0.14,sprintf('length'),'fontsize',9);
        %xtickangle(0);
        
        %xticks([]);
        xticklabels([]);
        
        ax = gca;
        ax.XAxis.FontSize = 9;
        ax.YAxis.FontSize = 10;%11
        
        if if_plot_rProb_scale == 1
            y_max = max(mean(ProbOffloadingInSeqLength)+std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)));
            y_min = min(mean(ProbOffloadingInSeqLength)-std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)));
            
            ylim([y_min-(y_max-y_min)*0.6 y_max+(y_max-y_min)*0.6]);%0.3
            set(gca,'YTick',[0:0.2:1]);
        else
            set(gca,'YTick',[0:0.2:1]);
        end
        
        set(gca,'box','off');% 取消右、上边框
        %xlabel('Length', 'FontSize', 11);
        % ylabel('Offloading rate', 'FontSize', 11);
        if if_title == 1
            title(sprintf('Offloading rate'), 'FontSize', 10);
        end
        
        %% Plot correlation
        %fig3 = figure(3);
        %set(gcf,'Position',[400 50 figure_xPixel figure_yPixel]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        nexttile
        
        % isvalid = (~isnan(gError_noChoice_collapsed_inOne)) & (~isnan(offloadingProb_all_inOne));
        isvalid = (~isnan(gAcc_noChoice_collapsed_inOne)) & (~isnan(offloadingProb_all_inOne));
        
        % x = gError_noChoice_collapsed_inOne(isvalid);
        x = gAcc_noChoice_collapsed_inOne(isvalid);
        y = offloadingProb_all_inOne(isvalid);
        n = 1;
        [p_mapping,S] = polyfit(x,y,n);
        r2 = 1 - (S.normr/norm(y - mean(y)))^2;
        [r, p_corr] = corr(x',y');
        tempTxt = 'p>0.05';
        if p_corr < 0.05
            tempTxt = 'p<0.05';
        end
        if p_corr < 0.01
            tempTxt = 'p<0.01';
        end
        if p_corr < 0.001
            tempTxt = 'p<0.001';
        end
        
        [p_mapping_poly,S_poly] = polyfit(x,y,3);
        r2_poly = 1 - (S_poly.normr/norm(y - mean(y)))^2;
        
        
        
        
        x_fit = 0:0.01:1;
        y_fit = polyval(p_mapping,x_fit);
        
        temp_color = [[166,97,26]/255;
            [223,194,125]/255;
            [128,205,193]/255;
            [1,133,113]/255];
        
        
        for target_seqLength=1:pointKindsNum
            %temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(gError_noChoice_collapsed{target_seqLength}));
            temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(gAcc_noChoice_collapsed{target_seqLength}));
            
            temp_size = temp_size ./ 5;
            
            %scatter(gError_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
            scatter(gAcc_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
                temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
        end
        plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        
        % text(0.1,0.9,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',10,'FontWeight','bold');
        % text(0.2,0.1,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',10,'FontWeight','bold');
        %         if if_monkey_D0_Z1 == 0
        %             text(0,0.30,sprintf('r=%.3f',r), 'fontsize',10,'FontWeight','normal');
        %             text(0,0.10,sprintf('%s',tempTxt), 'fontsize',10,'FontWeight','normal');
        %         elseif if_monkey_D0_Z1 == 1
        %             text(0.2,1,sprintf('r=%.3f',r), 'fontsize',10,'FontWeight','normal');
        %             text(0,0.10,sprintf('%s',tempTxt), 'fontsize',10,'FontWeight','normal');
        %         end
        
        set(gca,'linewidth',1.5)
        % xlim([0 1]);
        x_min = 0;
        x_max = 1;
        xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05]);
        
        y_min = 0;
        y_max = 1;
        %ylim([0 1]);
        ylim([y_min-(y_max-y_min)*0.05 y_max+(y_max-y_min)*0.05]);
        
        set(gca, 'FontSize', 10)
        %set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        set(gca,'XTick',[0:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        %set(gca,'XTick',[0,0.5]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        
        xtickangle(0);
        
        % set(gca,'XTick',[0:0.5:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        %set(gca,'YTick',[0:0.5:1]);%设置要显示坐标刻度的范围
        set(gca,'YTick',[0:1]);%设置要显示坐标刻度的范围
        set(gca,'box','off');% 取消右、上边框
        % xlabel('Error rate', 'FontSize', 11);
        %xlabel('Recall accuracy', 'FontSize', 10);
        %text(0.7,-0.17,sprintf('Recall accuracy'),'fontsize',9);
        
        temp_ylabel = ylabel('Offloading rate', 'FontSize', 10);
        % temp_Position = temp_ylabel.Position;
        % temp_Position(1) = temp_Position(1) - 0.03 + 1;
        % temp_Position(2) = temp_Position(2) - 0.01 - 0.2;
        % temp_ylabel.Position = temp_Position;
        if if_title == 1
            %title('Correlation', 'FontSize', 10);
        end
        title(sprintf('r=%.3f, %s',r,tempTxt), 'FontSize', 10);
        
        
        
        %% Plot correlation of endingHold
        nexttile
        
        temp_x = endingHold_noCoiceMemory_seq';
        temp_y = offloadingProb_all_inOne';
        tempBoolIndex = ~isnan(temp_x);
        
        x = temp_x(~isnan(temp_x));
        y = temp_y(1:sum(numSeq(1:valid_length)));
        y = y(~isnan(temp_x));
        
        [r,p] = corr(x,y);
        
        mdl = fitglm(x,y);
        
        for tempi=1:valid_length
            temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
            tempIndex = find(tempBoolIndex(temp_range)==true);
            temp_range2 = temp_range(tempIndex);
            
            temp_size = ((tempi.^3)*2 + 3) .* ones(1, length(temp_range2));
            
            temp_size = temp_size ./ 5;
            
            temp_h = scatter(temp_x(temp_range2), temp_y(temp_range2), ...
                temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
            
        end
        
        [temp_xmin,temp_xmax] = bounds(x);
        [temp_ymin,temp_ymax] = bounds(y);
        
        temp_ymin = 0;
        temp_ymax = 1;
        
        x_fit = temp_xmin:0.001:temp_xmax;
        y_fit = predict(mdl,x_fit')';
        
        plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        hold on
        
        tempTxt = 'p>0.05';
        if p < 0.05
            tempTxt = 'p<0.05';
        end
        if p < 0.01
            tempTxt = 'p<0.01';
        end
        if p < 0.001
            tempTxt = 'p<0.001';
        end
        
        %text(temp_xmin+(temp_xmax-temp_xmin)*0.57,temp_ymin+(temp_ymax-temp_ymin)*0.33,sprintf('r=%.3f',r), 'fontsize',10,'FontWeight','normal');
        %text(temp_xmin+(temp_xmax-temp_xmin)*0.57,temp_ymin+(temp_ymax-temp_ymin)*0.10,sprintf('%s',tempTxt), 'fontsize',10,'FontWeight','normal');
        
        
        xlim([temp_xmin-(temp_xmax-temp_xmin)*0.05 temp_xmax+(temp_xmax-temp_xmin)*0.05]);
        ylim([temp_ymin-(temp_ymax-temp_ymin)*0.05 temp_ymax+(temp_ymax-temp_ymin)*0.05]);
        
        set(gca,'linewidth',1.5)
        set(gca, 'FontSize', 10)
        set(gca,'box','off');% 取消右、上边框
        set(gca,'XTick',0.5:1.5,'FontSize', 12);%给坐标加标签
        set(gca,'YTick',0:1,'FontSize', 10);%给坐标加标签
        
        %xlabel(sprintf('Submission RT'), 'FontSize', 11, 'FontWeight', 'normal');
        %temp_title = title(sprintf('r=%.3f, p=%.3f',r,p), 'FontSize', 11, 'FontWeight', 'bold');
        %temp_title.Interpreter = 'none';
        
        ylabel(sprintf('Offloading rate'), 'FontSize', 10, 'FontWeight', 'normal');
        
        title(sprintf('r=%.3f, %s',r,tempTxt), 'FontSize', 10);
        
    end
    
    
end


%% End
% cd 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts'
cd 'C:\ASDROOT\STUDY\MonkeyDataAnalysis'