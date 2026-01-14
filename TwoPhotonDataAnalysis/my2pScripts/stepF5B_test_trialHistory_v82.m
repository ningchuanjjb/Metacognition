% Chuan's 12th script (20251214)
% This script: To conduct trial reward history analysis, related to figure 5.
%% Initialization
close all

if_compute = 1;%1
if_plot = 1;%1


% if_fitEnd = 1;
if_fitNormal0_fitEnd1_fitBeta2_fitExp3 = 3;%3

if_plot_fineTuningMismatch = 0;%0

if_useFittedW = 1;%1

% range_trialHistory = 1:50;%99-->79
% range_trialHistory = 1:50;%99-->79
range_trialHistory = 1:20;%50

range_trialHistory_min = 11;%11

if_xAxia_normal0_log1 = 0;%0

if_scoreType_cumulative0_meanS1_singleTrial2 = 1;%1


if_noFit0_fitScore1_fitHistroyW2 = 2;%2

if if_noFit0_fitScore1_fitHistroyW2 == 1
    if_scoreType_cumulative0_meanS1_singleTrial2_raw = if_scoreType_cumulative0_meanS1_singleTrial2;
    if_scoreType_cumulative0_meanS1_singleTrial2 = 2;
end
if if_noFit0_fitScore1_fitHistroyW2 == 2
    if_scoreType_cumulative0_meanS1_singleTrial2_raw = if_scoreType_cumulative0_meanS1_singleTrial2;
    if_scoreType_cumulative0_meanS1_singleTrial2 = 2;
end

% score_choiceMemoryCorrect = 2;
% score_choiceMemoryError = 0;
% score_choiceOffloadCorrect = 1;
% score_choiceOffloadError = 0;
% score_noChoiceCorrect = 2;
% score_noChoiceError = 0;

if if_monkey_D0_Z1 == 0
    score_choiceMemoryCorrect = 3;
    score_choiceMemoryError = 0;
    score_choiceOffloadCorrect = 2;
    score_choiceOffloadError = 0;
    score_noChoiceCorrect = 3;
    score_noChoiceError = 0;
    
elseif if_monkey_D0_Z1 == 1
    %     score_choiceMemoryCorrect = 3;
    %     score_choiceMemoryError = 0;
    %     score_choiceOffloadCorrect = 2;
    %     score_choiceOffloadError = 0;
    %     score_noChoiceCorrect = 3;
    %     score_noChoiceError = 0;
    
    score_choiceMemoryCorrect = 2;
    score_choiceMemoryError = 0;
    score_choiceOffloadCorrect = 1;
    score_choiceOffloadError = 0;
    score_noChoiceCorrect = 2;
    score_noChoiceError = 0;
    
end


% color_choiceMemoryHigh = [127,201,127]/255;%[0.5,1,0.5]
% color_choiceMemoryLow = [0,0.5,0];%[0,0.5,0]
% color_choiceOffloadHigh = [0.5,0,0];%[1,0.5,0.5]
% color_choiceOffloadLow = [1,0.5,0.5];%[0.5,0,0]
% color_highMeta = [0.4660 0.6740 0.1880];
% color_lowMeta = [0.6350 0.0780 0.1840];
% color_choiceMemory = [0.4660 0.6740 0.1880];
% color_choiceOffload = [0.6350 0.0780 0.1840];


color_choiceMemoryHigh = [128,205,193]/255; %[146,197,222]/255
color_choiceMemoryLow = [1,133,113]/255; %[5,113,176]/255
color_choiceOffloadHigh = [166,97,26]/255; %[202,0,32]/255
color_choiceOffloadLow = [223,194,125]/255; %[244,165,130]/255
color_choiceMemory = [1,133+20,113]/255;
color_choiceOffload = [166+40,97,26]/255;

color_highMeta = color_choiceMemory;
color_lowMeta = color_choiceOffload;


weight_trialHistory = ones(1,length(range_trialHistory));
weight_trialHistory_raw = weight_trialHistory;

% normStd = 10;
% x = range_trialHistory;
% y = normpdf(x,0,normStd);
% plot(x,y)
% weight_trialHistory = y;


% temp_load = load([output_path,'\','trial_para.mat']);
% trial_para = temp_load.trial_para;
% trial_para_isCorrect = trial_para.isCorrect;
trial_para_isCorrect;

choiceMemoryCorrectBoolIndex;
choiceMemoryErrorBoolIndex;

choiceOffloadCorrectBoolIndex = choiceOffloadBoolIndex & (trial_para_isCorrect==1);
choiceOffloadErrorBoolIndex = choiceOffloadBoolIndex & (trial_para_isCorrect==-1);

noChoiceCorrectBoolIndex = (~choiceBoolIndex) & (trial_para_isCorrect==1);
noChoiceErrorBoolIndex = (~choiceBoolIndex) & (trial_para_isCorrect==-1);

% sum(choiceMemoryCorrectBoolIndex)+sum(choiceMemoryErrorBoolIndex)+sum(choiceOffloadCorrectBoolIndex)+sum(choiceOffloadErrorBoolIndex)+sum(noChoiceCorrectBoolIndex)+sum(noChoiceErrorBoolIndex)

%%

for tempi=1:4
    temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
    
    temptempBoolIndex = ismember(seqIndex,temp_range)';
    
    if tempi==1
        trialBoolIndex_length1 = temptempBoolIndex;
    elseif tempi==2
        trialBoolIndex_length2 = temptempBoolIndex;
    elseif tempi==3
        trialBoolIndex_length3 = temptempBoolIndex;
    elseif tempi==4
        trialBoolIndex_length4 = temptempBoolIndex;
    end
end
trialBoolIndex_length1;
trialBoolIndex_length2;
trialBoolIndex_length3;
trialBoolIndex_length4;

% trialNum_A = decodingDataSimplified.extraForMerged.trialNum_A;
trialNum_B = decodingDataSimplified.extraForMerged.trialNum_B;
trialNum_A = length(trialBoolIndex_length1) - trialNum_B;

%% fit history W
if if_compute == 1
    if if_noFit0_fitScore1_fitHistroyW2 == 2
        
        %     normStd = 10;
        %     x = range_trialHistory;
        %     y = normpdf(x,0,normStd);
        %     % plot(x,y)
        %     weight_trialHistory = y;
        
        
        tempTrialTypeHistory = nan(trial_num,length(range_trialHistory));
        for tempi=1:trial_num
            for  tempj=1:length(range_trialHistory)
                temp_bufferSize = range_trialHistory(tempj);
                if tempi<=temp_bufferSize
                    continue
                end
                
                temp_trialRange = (tempi-temp_bufferSize):(tempi-1);
                
                if tempi == 880
                    a = 1;
                end
                
                % Correct tempTrialTypeHistory in terms of two session translation
                if ismember(trialNum_A,temp_trialRange)
                    continue
                end
                
                if if_scoreType_cumulative0_meanS1_singleTrial2 == 0
                    %tempWeight_trialRange = weight_trialHistory(temp_bufferSize:-1:1);
                elseif if_scoreType_cumulative0_meanS1_singleTrial2 == 1
                    %tempWeight_trialRange = weight_trialHistory(temp_bufferSize:-1:1);
                    %tempWeight_trialRange = tempWeight_trialRange ./ sum(tempWeight_trialRange);
                elseif if_scoreType_cumulative0_meanS1_singleTrial2 == 2
                    tempWeight_trialRange = weight_trialHistory(temp_bufferSize:-1:1);
                    
                    tempWeight_trialRange(2:end) = 0;
                    
                    tempWeight_trialRange = tempWeight_trialRange ./ sum(tempWeight_trialRange);
                end
                
                tempTrialCount_choiceMemoryCorrect = sum(choiceMemoryCorrectBoolIndex(temp_trialRange).*tempWeight_trialRange);
                tempTrialCount_choiceMemoryError = sum(choiceMemoryErrorBoolIndex(temp_trialRange).*tempWeight_trialRange);
                tempTrialCount_choiceOffloadCorrect = sum(choiceOffloadCorrectBoolIndex(temp_trialRange).*tempWeight_trialRange);
                tempTrialCount_choiceOffloadError = sum(choiceOffloadErrorBoolIndex(temp_trialRange).*tempWeight_trialRange);
                tempTrialCount_noChoiceCorrect = sum(noChoiceCorrectBoolIndex(temp_trialRange).*tempWeight_trialRange);
                tempTrialCount_noChoiceError = sum(noChoiceErrorBoolIndex(temp_trialRange).*tempWeight_trialRange);
                
                temp_score = tempTrialCount_choiceMemoryCorrect*score_choiceMemoryCorrect + ...
                    tempTrialCount_choiceMemoryError*score_choiceMemoryError + ...
                    tempTrialCount_choiceOffloadCorrect*score_choiceOffloadCorrect + ...
                    tempTrialCount_choiceOffloadError*score_choiceOffloadError + ...
                    tempTrialCount_noChoiceCorrect*score_noChoiceCorrect+ ...
                    tempTrialCount_noChoiceError*score_noChoiceError;
                
                tempTrialTypeHistory(tempi,tempj) = temp_score;
                
            end
        end
        tempTrialTypeHistory;
        
        
        
        tempChoiceTypeBoolIndex = nan(1,length(choiceMemoryCorrectBoolIndex));
        for tempi=1:trial_num
            if choiceMemoryBoolIndex(tempi) == true
                tempChoiceTypeBoolIndex(tempi) = true;
            elseif choiceOffloadBoolIndex(tempi) == true
                tempChoiceTypeBoolIndex(tempi) = false;
            end
        end
        
        % remove length 4 trials
        %tempChoiceTypeBoolIndex(trialBoolIndex_length4) = nan;
        
        tempChoiceTypeBoolIndex_T = tempChoiceTypeBoolIndex';
        
        
        %% normMean_optimal
        normMean_optimal = 1;%1
        
        
        %% normStd_optimal
        if if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 1
            endRange = 1:length(range_trialHistory);
            
            normMean = normMean_optimal;%0
            normStdRange = 1:0.2:length(range_trialHistory);%0.02,0.2
            
            temp_r2_fitting_2d = nan(length(range_trialHistory),length(normStdRange));
            AUROC_fitting_2d = nan(length(range_trialHistory),length(normStdRange));
            
            for tempIndex = 1:length(endRange)
                temp_r2_fitting_1d = nan(size(normStdRange));
                AUROC_fitting_1d = nan(size(normStdRange));
                
                %for tempi=1:length(normStdRange)
                parfor tempi=1:length(normStdRange)
                    temp_mdl = [];
                    
                    normStd = normStdRange(tempi);
                    %x = range_trialHistory;
                    x = range_trialHistory(1:endRange(tempIndex)); %#ok<*PFBNS>
                    y = normpdf(x,normMean,normStd);
                    %plot(x,y)
                    weight_trialHistory = y./sum(y); %#ok<*PFTUSW>
                    
                    temp1 = tempTrialTypeHistory(:,x) .* weight_trialHistory;
                    %temp2 = sum(temp1,2,'omitnan');
                    temp2 = sum(temp1,2);
                    
                    tempHistorySum = temp2;
                    
                    tempHistorySum;
                    tempChoiceTypeBoolIndex_T;
                    
                    x = tempHistorySum;
                    y = tempChoiceTypeBoolIndex_T;
                    
                    %temptempBoolIndex = (~isnan(x)) & (~isnan(y));
                    %x = x(temptempBoolIndex);
                    %y = y(temptempBoolIndex);
                    
                    temp_if_intercept = 1;
                    
                    if temp_if_intercept == 1
                        %temp_mdl = fitglm(x,y,'linear');
                        %temp_mdl = fitglm(x,y,'quadratic');
                        temp_mdl = fitglm(x,y,'linear','Distribution','binomial');
                        beta = temp_mdl.Coefficients.Estimate;
                        beta(1) = [];
                    elseif temp_if_intercept == 0
                        temp_mdl = fitglm(x,y,'linear','intercept',false);
                        beta = temp_mdl.Coefficients.Estimate;
                    end
                    temp_r2_fitting_1d(tempi) = temp_mdl.Rsquared.Adjusted;
                    
                    
                    % auroc
                    y_pred = predict(temp_mdl,x);
                    
                    temptempBoolIndex = (~isnan(x)) & (~isnan(y));
                    
                    y_valid = y(temptempBoolIndex);
                    y_pred_valid = y_pred(temptempBoolIndex);
                    
                    temp_thresholdRange = 0.01:0.001:0.99;
                    temp_AUROC_multi = zeros(1,length(temp_thresholdRange));
                    for temptempi=1:length(temp_thresholdRange)
                        
                        temp_threshold = temp_thresholdRange(temptempi);
                        predictLabel = y_pred_valid > temp_threshold;%0.5
                        trueLabel = y_valid;
                        
                        temp_Hit = sum(predictLabel(trueLabel==true))/sum(trueLabel==true);
                        temp_FA = sum(predictLabel(trueLabel==false))/sum(trueLabel==false);
                        
                        x = [0 temp_FA 1];
                        y = [0 temp_Hit 1];
                        temp_AUROC = trapz(x,y);
                        temp_AUROC_multi(temptempi) = temp_AUROC;
                    end
                    [M,I] = max(temp_AUROC_multi); %#ok<*ASGLU>
                    AUROC_fitting_1d(tempi) = M;
                    
                end
                %[M,I] = max(temp_r2_fitting);
                %temp_r2_optimal = M;
                
                %[M,I] = max(AUROC_fitting);
                %temp_r2_optimal = temp_r2_fitting(I);
                %AUROC_optimal = M;
                %normStd_optimal = normStdRange(I);
                
                AUROC_fitting_1d;
                temp_r2_fitting_1d;
                
                AUROC_fitting_2d(tempIndex,:) = AUROC_fitting_1d;
                temp_r2_fitting_2d(tempIndex,:) = temp_r2_fitting_1d;
                
            end
            
            AUROC_fitting_2d;
            temp_r2_fitting_2d;
            
            [M,I] = max(AUROC_fitting_2d,[],'all','linear');
            %[M,I] = max(temp_r2_fitting_2d,[],'all','linear');
            [row,col] = ind2sub(size(AUROC_fitting_2d),I);
            
            temp_r2_optimal = temp_r2_fitting_2d(row,col);
            AUROC_optimal = AUROC_fitting_2d(row,col)
            end_optimal = endRange(row);
            normStd_optimal = normStdRange(col);
            
            
        elseif if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 0
            
            
            if false
                x_step = 0.01;
                x = 0:x_step:100;
                y = exppdf(x,2);
                %cftool(x,y);
                %mean(y*0.1)
                %mean(y)
                
                sum(x.*x_step.*y)
                
            end
            
            %% normStd_optimal
            %normStdRange = 0.1:0.02:length(range_trialHistory);
            normStdRange = 1:0.02:length(range_trialHistory);
            normMean = normMean_optimal;%0
            temp_r2_fitting = nan(size(normStdRange));
            AUROC_fitting = nan(size(normStdRange));
            %for tempi=1:length(normStdRange)
            parfor tempi=1:length(normStdRange)
                temp_mdl = [];
                
                normStd = normStdRange(tempi);
                x = range_trialHistory;
                y = normpdf(x,normMean,normStd);
                %plot(x,y)
                weight_trialHistory = y./sum(y); %#ok<*PFTUSW>
                
                temp1 = tempTrialTypeHistory .* weight_trialHistory;
                %temp2 = sum(temp1,2,'omitnan');
                temp2 = sum(temp1,2);
                
                tempHistorySum = temp2;
                
                tempHistorySum;
                tempChoiceTypeBoolIndex_T;
                
                x = tempHistorySum;
                y = tempChoiceTypeBoolIndex_T;
                
                %temptempBoolIndex = (~isnan(x)) & (~isnan(y));
                %x = x(temptempBoolIndex);
                %y = y(temptempBoolIndex);
                
                temp_if_intercept = 1;
                
                if temp_if_intercept == 1
                    %temp_mdl = fitglm(x,y,'linear');
                    %temp_mdl = fitglm(x,y,'quadratic');
                    temp_mdl = fitglm(x,y,'linear','Distribution','binomial');
                    beta = temp_mdl.Coefficients.Estimate;
                    beta(1) = [];
                elseif temp_if_intercept == 0
                    temp_mdl = fitglm(x,y,'linear','intercept',false);
                    beta = temp_mdl.Coefficients.Estimate;
                end
                temp_r2_fitting(tempi) = temp_mdl.Rsquared.Adjusted;
                
                
                % auroc
                y_pred = predict(temp_mdl,x);
                
                temptempBoolIndex = (~isnan(x)) & (~isnan(y));
                
                y_valid = y(temptempBoolIndex);
                y_pred_valid = y_pred(temptempBoolIndex);
                
                temp_thresholdRange = 0.01:0.001:0.99;
                temp_AUROC_multi = zeros(1,length(temp_thresholdRange));
                for temptempi=1:length(temp_thresholdRange)
                    
                    temp_threshold = temp_thresholdRange(temptempi);
                    predictLabel = y_pred_valid > temp_threshold;%0.5
                    trueLabel = y_valid;
                    
                    temp_Hit = sum(predictLabel(trueLabel==true))/sum(trueLabel==true);
                    temp_FA = sum(predictLabel(trueLabel==false))/sum(trueLabel==false);
                    
                    x = [0 temp_FA 1];
                    y = [0 temp_Hit 1];
                    temp_AUROC = trapz(x,y);
                    temp_AUROC_multi(temptempi) = temp_AUROC;
                end
                [M,I] = max(temp_AUROC_multi); %#ok<*ASGLU>
                AUROC_fitting(tempi) = M;
                
            end
            %[M,I] = max(temp_r2_fitting);
            %temp_r2_optimal = M;
            
            [M,I] = max(AUROC_fitting);
            temp_r2_optimal = temp_r2_fitting(I);
            AUROC_optimal = M;
            normStd_optimal = normStdRange(I);
            
        end
        
        if if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 3
            
            if false
                x_step = 1;
                x = 1:x_step:length(range_trialHistory);
                
                temp_expMu = expMu_optimal;%2
                
                y = exppdf(x,temp_expMu);
                %cftool(x,y);
                %mean(y*0.1)
                %mean(y)
                
                %sum(x.*x_step.*y)
                
            end
            
            %% normStd_optimal
            expMuRange = 0.01:0.02:length(range_trialHistory);
            temp_r2_fitting = nan(size(expMuRange));
            AUROC_fitting = nan(size(expMuRange));
            %for tempi=1:length(expMuRange)
            parfor tempi=1:length(expMuRange)
                temp_mdl = [];
                
                temp_expMu = expMuRange(tempi);
                x = range_trialHistory;
                y = exppdf(x,temp_expMu);
                %plot(x,y)
                weight_trialHistory = y./sum(y); %#ok<*PFTUSW>
                
                temp1 = tempTrialTypeHistory .* weight_trialHistory;
                %temp2 = sum(temp1,2,'omitnan');
                temp2 = sum(temp1,2);
                
                tempHistorySum = temp2;
                
                tempHistorySum;
                tempChoiceTypeBoolIndex_T;
                
                x = tempHistorySum;
                y = tempChoiceTypeBoolIndex_T;
                
                %temptempBoolIndex = (~isnan(x)) & (~isnan(y));
                %x = x(temptempBoolIndex);
                %y = y(temptempBoolIndex);
                
                temp_if_intercept = 1;
                
                if temp_if_intercept == 1
                    %temp_mdl = fitglm(x,y,'linear');
                    %temp_mdl = fitglm(x,y,'quadratic');
                    temp_mdl = fitglm(x,y,'linear','Distribution','binomial');
                    beta = temp_mdl.Coefficients.Estimate;
                    beta(1) = [];
                elseif temp_if_intercept == 0
                    temp_mdl = fitglm(x,y,'linear','intercept',false);
                    beta = temp_mdl.Coefficients.Estimate;
                end
                temp_r2_fitting(tempi) = temp_mdl.Rsquared.Adjusted;
                
                
                % auroc
                y_pred = predict(temp_mdl,x);
                
                temptempBoolIndex = (~isnan(x)) & (~isnan(y));
                
                y_valid = y(temptempBoolIndex);
                y_pred_valid = y_pred(temptempBoolIndex);
                
                
                temp_thresholdRange = 0.01:0.001:0.99;
                
                temp_choiceMemory_hit_multi = zeros(1,length(temp_thresholdRange));
                temp_choiceMemory_falseAlarm_multi = zeros(1,length(temp_thresholdRange));
                
                for temptempi=1:length(temp_thresholdRange)
                    
                    temp_threshold = temp_thresholdRange(temptempi);
                    
                    predictLabel = y_pred_valid > temp_threshold;
                    trueLabel = y_valid;
                    
                    temp_choiceMemory_hit = sum(predictLabel(trueLabel==true))/sum(trueLabel==true);
                    temp_choiceMemory_falseAlarm = sum(predictLabel(trueLabel==false))/sum(trueLabel==false);
                    
                    temp_choiceMemory_hit_multi(temptempi) = temp_choiceMemory_hit;
                    temp_choiceMemory_falseAlarm_multi(temptempi) = temp_choiceMemory_falseAlarm;
                end
                x = temp_choiceMemory_falseAlarm_multi;
                y = temp_choiceMemory_hit_multi;
                
                x = x(end:-1:1);
                y = y(end:-1:1);
                
                temp_AUROC = trapz(x,y);
                AUROC_fitting(tempi) = temp_AUROC;
                
            end
            %[M,I] = max(temp_r2_fitting);
            %temp_r2_optimal = M;
            
            [M,I] = max(AUROC_fitting);
            temp_r2_optimal = temp_r2_fitting(I);
            AUROC_optimal = M;
            expMu_optimal = expMuRange(I);
            
        end
        
        %% glm based fit
        %if false
        if if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 2
            a = 1;
            
            %close all
            
            x = tempTrialTypeHistory(:,1:length(range_trialHistory));
            y = tempChoiceTypeBoolIndex_T;
            
            %x = tempTrialTypeHistory(1:1035,1:20);
            %y = tempChoiceTypeBoolIndex_T(1:1035);
            
            %x = tempTrialTypeHistory(1085:end,1:20);
            %y = tempChoiceTypeBoolIndex_T(1085:end);
            
            temp_if_intercept = 1;
            
            if temp_if_intercept == 1
                %temp_mdl = fitglm(x,y,'linear');
                temp_mdl = fitglm(x,y,'linear','Distribution','binomial');
                beta = temp_mdl.Coefficients.Estimate;
                beta(1) = [];
            elseif temp_if_intercept == 0
                temp_mdl = fitglm(x,y,'linear','intercept',false);
                beta = temp_mdl.Coefficients.Estimate;
            end
            temp_r2_fitting_C = temp_mdl.Rsquared.Adjusted;
            
            %plot(beta)
            
            % auroc
            y_pred = predict(temp_mdl,x);
            
            temptempBoolIndex = (~isnan(x(:,end))) & (~isnan(y));
            %temptempBoolIndex = ~isnan(y);
            
            y_valid = y(temptempBoolIndex);
            y_pred_valid = y_pred(temptempBoolIndex);
            
            temp_thresholdRange = 0.01:0.001:0.99;
            temp_AUROC_multi = zeros(1,length(temp_thresholdRange));
            for temptempi=1:length(temp_thresholdRange)
                
                temp_threshold = temp_thresholdRange(temptempi);
                predictLabel = y_pred_valid > temp_threshold;%0.5
                trueLabel = y_valid;
                
                temp_Hit = sum(predictLabel(trueLabel==true))/sum(trueLabel==true);
                temp_FA = sum(predictLabel(trueLabel==false))/sum(trueLabel==false);
                
                x = [0 temp_FA 1];
                y = [0 temp_Hit 1];
                temp_AUROC = trapz(x,y);
                temp_AUROC_multi(temptempi) = temp_AUROC;
            end
            [M,I] = max(temp_AUROC_multi); %#ok<*ASGLU>
            AUROC_fitting_C = M;
            
            temp_r2_optimal = temp_r2_fitting_C;
            AUROC_optimal = M
            
            weight_beta = beta';
            
        end
        
        if false
            
            close all
            
            temp_asd = 50;
            
            temp_r2_fitting_C = nan(1,temp_asd);
            AUROC_fitting_C = nan(1,temp_asd);
            
            for tempi=1:temp_asd
                
                x = tempTrialTypeHistory(:,tempi);
                y = tempChoiceTypeBoolIndex_T;
                
                temp_if_intercept = 1;
                
                if temp_if_intercept == 1
                    %temp_mdl = fitglm(x,y,'linear');
                    temp_mdl = fitglm(x,y,'linear','Distribution','binomial');
                    beta = temp_mdl.Coefficients.Estimate;
                    beta(1) = [];
                elseif temp_if_intercept == 0
                    temp_mdl = fitglm(x,y,'linear','intercept',false);
                    beta = temp_mdl.Coefficients.Estimate;
                end
                temp_r2_fitting_C(tempi) = temp_mdl.Rsquared.Adjusted;
                
                % auroc
                y_pred = predict(temp_mdl,x);
                
                temptempBoolIndex = (~isnan(x(:,end))) & (~isnan(y));
                
                y_valid = y(temptempBoolIndex);
                y_pred_valid = y_pred(temptempBoolIndex);
                
                temp_thresholdRange = 0.01:0.001:0.99;
                temp_AUROC_multi = zeros(1,length(temp_thresholdRange));
                for temptempi=1:length(temp_thresholdRange)
                    
                    temp_threshold = temp_thresholdRange(temptempi);
                    predictLabel = y_pred_valid > temp_threshold;%0.5
                    trueLabel = y_valid;
                    
                    temp_Hit = sum(predictLabel(trueLabel==true))/sum(trueLabel==true);
                    temp_FA = sum(predictLabel(trueLabel==false))/sum(trueLabel==false);
                    
                    x = [0 temp_FA 1];
                    y = [0 temp_Hit 1];
                    temp_AUROC = trapz(x,y);
                    temp_AUROC_multi(temptempi) = temp_AUROC;
                end
                [M,I] = max(temp_AUROC_multi); %#ok<*ASGLU>
                AUROC_fitting_C(tempi) = M;
                
            end
            
            %plot(AUROC_fitting_C);
            %plot(temp_r2_fitting_C);
            
            
        end
        
        
        if false
            
            close all
            
            temp_asd = 15;
            
            temp_range = 1:1035;
            %temp_range = 1085:2070;
            %temp_range = 1:2070;
            
            temp_r2_fitting_CPD = nan(1,temp_asd);
            AUROC_fitting_CPD  = nan(1,temp_asd);
            
            for tempi=1:temp_asd
                
                temptempBoolIndexB = false(1,size(tempTrialTypeHistory,2));
                temptempBoolIndexB(1:temp_asd) = true;
                temptempBoolIndexB(tempi) = false;
                
                
                x = tempTrialTypeHistory(temp_range,temptempBoolIndexB);
                y = tempChoiceTypeBoolIndex_T(temp_range);
                
                temp_if_intercept = 1;
                
                if temp_if_intercept == 1
                    %temp_mdl = fitglm(x,y,'linear');
                    temp_mdl = fitglm(x,y,'linear','Distribution','binomial');
                    beta = temp_mdl.Coefficients.Estimate;
                    beta(1) = [];
                elseif temp_if_intercept == 0
                    temp_mdl = fitglm(x,y,'linear','intercept',false);
                    beta = temp_mdl.Coefficients.Estimate;
                end
                temp_r2_fitting_CPD(tempi) = temp_mdl.Rsquared.Adjusted;
                
                % auroc
                y_pred = predict(temp_mdl,x);
                
                temptempBoolIndex = (~isnan(x(:,end))) & (~isnan(y));
                
                y_valid = y(temptempBoolIndex);
                y_pred_valid = y_pred(temptempBoolIndex);
                
                temp_thresholdRange = 0.01:0.001:0.99;
                temp_AUROC_multi = zeros(1,length(temp_thresholdRange));
                for temptempi=1:length(temp_thresholdRange)
                    
                    temp_threshold = temp_thresholdRange(temptempi);
                    predictLabel = y_pred_valid > temp_threshold;%0.5
                    trueLabel = y_valid;
                    
                    temp_Hit = sum(predictLabel(trueLabel==true))/sum(trueLabel==true);
                    temp_FA = sum(predictLabel(trueLabel==false))/sum(trueLabel==false);
                    
                    x = [0 temp_FA 1];
                    y = [0 temp_Hit 1];
                    temp_AUROC = trapz(x,y);
                    temp_AUROC_multi(temptempi) = temp_AUROC;
                end
                [M,I] = max(temp_AUROC_multi); %#ok<*ASGLU>
                AUROC_fitting_CPD(tempi) = M;
                
            end
            
            %plot(AUROC_fitting_C);
            %plot(temp_r2_fitting_C);
            
            
            x = tempTrialTypeHistory(temp_range,1:temp_asd);
            y = tempChoiceTypeBoolIndex_T(temp_range);
            
            temp_if_intercept = 1;
            
            if temp_if_intercept == 1
                %temp_mdl = fitglm(x,y,'linear');
                temp_mdl = fitglm(x,y,'linear','Distribution','binomial');
                beta = temp_mdl.Coefficients.Estimate;
                beta(1) = [];
            elseif temp_if_intercept == 0
                temp_mdl = fitglm(x,y,'linear','intercept',false);
                beta = temp_mdl.Coefficients.Estimate;
            end
            temp_r2_fitting_C = temp_mdl.Rsquared.Adjusted;
            
            % auroc
            y_pred = predict(temp_mdl,x);
            
            temptempBoolIndex = (~isnan(x(:,end))) & (~isnan(y));
            %temptempBoolIndex = ~isnan(y);
            
            y_valid = y(temptempBoolIndex);
            y_pred_valid = y_pred(temptempBoolIndex);
            
            temp_thresholdRange = 0.01:0.001:0.99;
            temp_AUROC_multi = zeros(1,length(temp_thresholdRange));
            for temptempi=1:length(temp_thresholdRange)
                
                temp_threshold = temp_thresholdRange(temptempi);
                predictLabel = y_pred_valid > temp_threshold;%0.5
                trueLabel = y_valid;
                
                temp_Hit = sum(predictLabel(trueLabel==true))/sum(trueLabel==true);
                temp_FA = sum(predictLabel(trueLabel==false))/sum(trueLabel==false);
                
                x = [0 temp_FA 1];
                y = [0 temp_Hit 1];
                temp_AUROC = trapz(x,y);
                temp_AUROC_multi(temptempi) = temp_AUROC;
            end
            [M,I] = max(temp_AUROC_multi); %#ok<*ASGLU>
            AUROC_fitting_C = M;
            
            
            AUROC_fitting_CPD;
            AUROC_fitting_C;
            temp_r2_fitting_CPD;
            temp_r2_fitting_C;
            
            diff_AUROC_fitting_CPD = (AUROC_fitting_C-AUROC_fitting_CPD)./AUROC_fitting_CPD;
            diff_temp_r2_fitting_CPD = (temp_r2_fitting_C-temp_r2_fitting_CPD)./temp_r2_fitting_CPD;
            
            figure(1)
            plot(diff_AUROC_fitting_CPD);
            
            figure(2)
            plot(diff_temp_r2_fitting_CPD);
            
            figure(3)
            plot(beta)
            
            
        end
        
    end
end


%%
% normStd = normStd_optimal;
% normMean = normMean_optimal;

range_trialHistory_raw = range_trialHistory;

% range_trialHistory = 1:ceil(normMean+3*normStd);

if if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 0
    normStd = normStd_optimal;
    normMean = normMean_optimal;
    
    temp1 = min(ceil(normMean+3*normStd),max(range_trialHistory_raw));
    temp1 = max(temp1,range_trialHistory_min);
    range_trialHistory = 1:temp1;
    
    x = range_trialHistory;
    y = normpdf(x,normMean,normStd);
    % plot(x,y)
    weight_trialHistory = y./sum(y);
    
elseif if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 1
    normStd = normStd_optimal;
    normMean = normMean_optimal;
    
    range_trialHistory = 1:end_optimal;
    
    %x = range_trialHistory_raw;
    x = range_trialHistory;
    y = normpdf(x,normMean,normStd);
    % plot(x,y)
    weight_trialHistory = y./sum(y);
    weight_trialHistory = weight_trialHistory(range_trialHistory);
    
elseif if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 2
    weight_trialHistory = weight_beta;
    
    if false
        x = 1:length(weight_trialHistory);
        y = weight_trialHistory;
        cftool(x,y);
    end
elseif if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 3
    expMu = expMu_optimal;
    
    x = range_trialHistory;
    y = exppdf(x,expMu);
    % plot(x,y)
    weight_trialHistory = y./sum(y);
end





if_scoreType_cumulative0_meanS1_singleTrial2 = if_scoreType_cumulative0_meanS1_singleTrial2_raw;

if if_useFittedW == 1
    weight_trialHistory_raw = weight_trialHistory;
end


%% fit score
if if_noFit0_fitScore1_fitHistroyW2 == 1
    score_choiceMemoryCorrect_raw = score_choiceMemoryCorrect;
    score_choiceMemoryError_raw = score_choiceMemoryError;
    score_choiceOffloadCorrect_raw = score_choiceOffloadCorrect;
    score_choiceOffloadError_raw = score_choiceOffloadError;
    score_noChoiceCorrect_raw = score_noChoiceCorrect;
    score_noChoiceError_raw = score_noChoiceError;
    
    choiceMemoryCorrectBoolIndex;
    choiceMemoryErrorBoolIndex;
    choiceOffloadCorrectBoolIndex;
    choiceOffloadErrorBoolIndex;
    noChoiceCorrectBoolIndex;
    noChoiceErrorBoolIndex;
    
    tempTrialTypeIndex = nan(1,length(choiceMemoryCorrectBoolIndex));
    for tempi=1:trial_num
        if choiceMemoryCorrectBoolIndex(tempi) == true
            tempTrialTypeIndex(tempi) = 1;
        elseif choiceMemoryErrorBoolIndex(tempi) == true
            tempTrialTypeIndex(tempi) = 2;
        elseif choiceOffloadCorrectBoolIndex(tempi) == true
            tempTrialTypeIndex(tempi) = 3;
        elseif choiceOffloadErrorBoolIndex(tempi) == true
            tempTrialTypeIndex(tempi) = 4;
        elseif noChoiceCorrectBoolIndex(tempi) == true
            tempTrialTypeIndex(tempi) = 5;
        elseif noChoiceErrorBoolIndex(tempi) == true
            tempTrialTypeIndex(tempi) = 6;
        end
    end
    
    tempTrialTypeIndex;
    choiceMemoryBoolIndex;
    choiceOffloadBoolIndex;
    
    tempChoiceTypeBoolIndex = nan(1,length(choiceMemoryCorrectBoolIndex));
    for tempi=1:trial_num
        if choiceMemoryBoolIndex(tempi) == true
            tempChoiceTypeBoolIndex(tempi) = true;
        elseif choiceOffloadBoolIndex(tempi) == true
            tempChoiceTypeBoolIndex(tempi) = false;
        end
    end
    
    tempTrialTypeIndex;
    tempChoiceTypeBoolIndex;
    range_trialHistory;
    
    temp_range_trialHistory = 1:10;
    %temp_range_trialHistory = 1:99;
    
    %temp_range_trialHistory = range_trialHistory;
    
    tempTrialTypeHistory = nan(trial_num,length(temp_range_trialHistory));
    for tempi=1:trial_num
        for  tempj=1:length(temp_range_trialHistory)
            temp_bufferSize = temp_range_trialHistory(tempj);
            if tempi<=temp_bufferSize
                continue
            end
            temp_trialRange = (tempi-temp_bufferSize):(tempi-1);
            temp_trialRange = temp_trialRange(1);
            
            tempType = tempTrialTypeIndex(temp_trialRange);
            tempTrialTypeHistory(tempi,tempj) = tempType;
        end
    end
    tempTrialTypeHistory;
    tempChoiceTypeBoolIndex;
    
    tempChoiceTypeBoolIndex_multi = repmat(tempChoiceTypeBoolIndex,length(temp_range_trialHistory),1)';
    tempChoiceTypeBoolIndex_multi_1d = reshape(tempChoiceTypeBoolIndex_multi',[],1);
    
    
    tempTrialTypeHistory_1d = reshape(tempTrialTypeHistory',[],1);
    
    tempTrialTypeBoolHistory_1d = nan(length(tempTrialTypeHistory_1d),6);
    for tempi=1:length(tempTrialTypeHistory_1d)
        if isnan(tempTrialTypeHistory_1d(tempi))
            continue
        end
        tempTrialTypeBoolHistory_1d(tempi,:) = false;
        tempTrialTypeBoolHistory_1d(tempi,tempTrialTypeHistory_1d(tempi)) = true;
    end
    
    x = tempTrialTypeBoolHistory_1d;
    y = tempChoiceTypeBoolIndex_multi_1d;
    
    temp_if_intercept = 0;
    
    if temp_if_intercept == 1
        temp_mdl = fitglm(x,y,'linear')
        beta = temp_mdl.Coefficients.Estimate;
        beta(1) = [];
    elseif temp_if_intercept == 0
        temp_mdl = fitglm(x,y,'linear','intercept',false)
        beta = temp_mdl.Coefficients.Estimate;
    end
    r2 = temp_mdl.Rsquared.Adjusted;
    
    score_choiceMemoryCorrect = beta(1)*1;
    score_choiceMemoryError = beta(2)*1;
    score_choiceOffloadCorrect = beta(3)*1;
    score_choiceOffloadError = beta(4)*1;
    score_noChoiceCorrect = beta(5)*1;
    score_noChoiceError = beta(6)*1;
    
    if_scoreType_cumulative0_meanS1_singleTrial2 = if_scoreType_cumulative0_meanS1_singleTrial2_raw;
end


%% score_trialHistory
score_trialHistory = nan(trial_num,length(range_trialHistory));
for tempi=1:trial_num
    for  tempj=1:length(range_trialHistory)
        temp_bufferSize = range_trialHistory(tempj);
        if tempi<=temp_bufferSize
            continue
        end
        
        %temp_trialRange = (tempi-temp_bufferSize):(tempi-1);
        temp_trialRange = (tempi-temp_bufferSize):(tempi-1);
        
        
        % Correct tempTrialTypeHistory in terms of two session translation
        if ismember(trialNum_A,temp_trialRange)
            continue
        end
        
        if if_scoreType_cumulative0_meanS1_singleTrial2 == 0
            tempWeight_trialRange = weight_trialHistory_raw(temp_bufferSize:-1:1);
        elseif if_scoreType_cumulative0_meanS1_singleTrial2 == 1
            tempWeight_trialRange = weight_trialHistory_raw(temp_bufferSize:-1:1);
            %if if_fitNormal0_fitEnd1_fitBeta2_fitExp3 < 2
            if if_fitNormal0_fitEnd1_fitBeta2_fitExp3 < 2 || if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 3
                tempWeight_trialRange = tempWeight_trialRange ./ sum(tempWeight_trialRange);
            end
        elseif if_scoreType_cumulative0_meanS1_singleTrial2 == 2
            tempWeight_trialRange = weight_trialHistory_raw(temp_bufferSize:-1:1);
            
            tempWeight_trialRange(2:end) = 0;
            
            tempWeight_trialRange = tempWeight_trialRange ./ sum(tempWeight_trialRange);
        end
        
        tempTrialCount_choiceMemoryCorrect = sum(choiceMemoryCorrectBoolIndex(temp_trialRange).*tempWeight_trialRange);
        tempTrialCount_choiceMemoryError = sum(choiceMemoryErrorBoolIndex(temp_trialRange).*tempWeight_trialRange);
        tempTrialCount_choiceOffloadCorrect = sum(choiceOffloadCorrectBoolIndex(temp_trialRange).*tempWeight_trialRange);
        tempTrialCount_choiceOffloadError = sum(choiceOffloadErrorBoolIndex(temp_trialRange).*tempWeight_trialRange);
        tempTrialCount_noChoiceCorrect = sum(noChoiceCorrectBoolIndex(temp_trialRange).*tempWeight_trialRange);
        tempTrialCount_noChoiceError = sum(noChoiceErrorBoolIndex(temp_trialRange).*tempWeight_trialRange);
        
        temp_score = tempTrialCount_choiceMemoryCorrect*score_choiceMemoryCorrect + ...
            tempTrialCount_choiceMemoryError*score_choiceMemoryError + ...
            tempTrialCount_choiceOffloadCorrect*score_choiceOffloadCorrect + ...
            tempTrialCount_choiceOffloadError*score_choiceOffloadError + ...
            tempTrialCount_noChoiceCorrect*score_noChoiceCorrect+ ...
            tempTrialCount_noChoiceError*score_noChoiceError;
        
        score_trialHistory(tempi,tempj) = temp_score;
        
        
        if tempi == 3
            if tempj == 2
                a = 1;
            end
        end
        
        if tempi == 17
            if tempj == 1
                a = 1;
            end
            if tempj == 15
                a = 1;
            end
        end
        
        if tempi == 18
            if tempj == 1
                a = 1;
            end
        end
        
    end
end
% score_trialHistory_raw = score_trialHistory;
% if if_scoreType_cumulative0_meanS1_singleTrial2 == 1
%     score_trialHistory = score_trialHistory_raw ./ range_trialHistory;
% end


%% match and mismatch trials
% match trials
trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory;
trialBoolIndex_memoryPrecisionLow_choiceOffload;
trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh;
trialBoolIndex_memoryPrecisionLow_choiceOffloadLow;

% mismatch trials
trialBoolIndex_memoryPrecisionLowError_choiceMemory;
trialBoolIndex_memoryPrecisionHigh_choiceOffload;
trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh;
trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow;

%
% trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh
% trialBoolIndex_memoryPrecisionLow_choiceOffloadLow
%
% trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh
% trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow


% temp_lowThreshold_meta = metaDecoderThreshold_delay1;
%
% highThreshold_meta = lowThreshold_meta;

% temp_lowThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel(choiceBoolIndex),40);%40
% temp_highThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel(choiceBoolIndex),60);%60
% % temp_lowThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel(choiceBoolIndex),50);%40
% % temp_highThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel(choiceBoolIndex),50);%60
% % temp_lowThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel(choiceBoolIndex),45);%40-->35
% % temp_highThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel(choiceBoolIndex),55);%60-->65
%
% % temp_lowThreshold_meta = prctile(meta_trialLevel(choiceBoolIndex),50);%40
% % temp_highThreshold_meta = prctile(meta_trialLevel(choiceBoolIndex),50);%60
% temp_lowThreshold_meta = prctile(meta_trialLevel(choiceBoolIndex),45);%40
% temp_highThreshold_meta = prctile(meta_trialLevel(choiceBoolIndex),55);%60


% temp_lowThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel(choiceBoolIndex),40);%40
% temp_highThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel(choiceBoolIndex),60);%60
% temp_lowThreshold_meta = prctile(meta_trialLevel(choiceBoolIndex),45);%45
% temp_highThreshold_meta = prctile(meta_trialLevel(choiceBoolIndex),55);%55


temp_lowThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel(choiceBoolIndex),40);%40
temp_highThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel(choiceBoolIndex),60);%60
temp_lowThreshold_meta = prctile(meta_trialLevel(choiceBoolIndex),45);%45
temp_highThreshold_meta = prctile(meta_trialLevel(choiceBoolIndex),55);%55



% if if_plot_fineTuningMismatch == 1
% neuron label + behavior label
%     temp_trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh = ...
%         trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh & (memoryPrecision_trialLevel>temp_highThreshold_memoryPrecision) & (meta_trialLevel>temp_highThreshold_meta);
%     temp_trialBoolIndex_memoryPrecisionLow_choiceOffloadLow = ...
%         trialBoolIndex_memoryPrecisionLow_choiceOffloadLow & (memoryPrecision_trialLevel<temp_lowThreshold_memoryPrecision) & (meta_trialLevel<temp_lowThreshold_meta);
%     temp_trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh = ...
%         trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh & (memoryPrecision_trialLevel<temp_lowThreshold_memoryPrecision) & (meta_trialLevel>temp_highThreshold_meta);
%     temp_trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow = ...
%         trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow & (memoryPrecision_trialLevel>temp_highThreshold_memoryPrecision) & (meta_trialLevel<temp_lowThreshold_meta);

% elseif if_plot_fineTuningMismatch == 0

%     % only neuron label
temp_trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh = ...
    trialBoolIndex_memoryPrecisionHigh_choice & trialBoolIndex_metaHigh_choice & (memoryPrecision_trialLevel>temp_highThreshold_memoryPrecision) & (meta_trialLevel>temp_highThreshold_meta);
temp_trialBoolIndex_memoryPrecisionLow_choiceOffloadLow = ...
    trialBoolIndex_memoryPrecisionLow_choice & trialBoolIndex_metaLow_choice & (memoryPrecision_trialLevel<temp_lowThreshold_memoryPrecision) & (meta_trialLevel<temp_lowThreshold_meta);
temp_trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh = ...
    trialBoolIndex_memoryPrecisionLow_choice & trialBoolIndex_metaHigh_choice & (memoryPrecision_trialLevel<temp_lowThreshold_memoryPrecision) & (meta_trialLevel>temp_highThreshold_meta);
temp_trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow = ...
    trialBoolIndex_memoryPrecisionHigh_choice & trialBoolIndex_metaLow_choice & (memoryPrecision_trialLevel>temp_highThreshold_memoryPrecision) & (meta_trialLevel<temp_lowThreshold_meta);
% temp_trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh = ...
%     trialBoolIndex_memoryPrecisionHigh_choice & trialBoolIndex_metaHigh_choice;
% temp_trialBoolIndex_memoryPrecisionLow_choiceOffloadLow = ...
%     trialBoolIndex_memoryPrecisionLow_choice & trialBoolIndex_metaLow_choice;
% temp_trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh = ...
%     trialBoolIndex_memoryPrecisionLow_choice & trialBoolIndex_metaHigh_choice;
% temp_trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow = ...
%     trialBoolIndex_memoryPrecisionHigh_choice & trialBoolIndex_metaLow_choice;

%neuron label + partial behavior label (correct/error)
%     temp_trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh = ...
%         choiceMemoryCorrectBoolIndex' & trialBoolIndex_memoryPrecisionHigh_choice & trialBoolIndex_metaHigh_choice & (memoryPrecision_trialLevel>temp_highThreshold_memoryPrecision);
%     temp_trialBoolIndex_memoryPrecisionLow_choiceOffloadLow = ...
%         trialBoolIndex_memoryPrecisionLow_choice & trialBoolIndex_metaLow_choice & (memoryPrecision_trialLevel<temp_lowThreshold_memoryPrecision);
%     temp_trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh = ...
%         choiceMemoryErrorBoolIndex' & trialBoolIndex_memoryPrecisionLow_choice & trialBoolIndex_metaHigh_choice & (memoryPrecision_trialLevel<temp_lowThreshold_memoryPrecision);
%     temp_trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow = ...
%         trialBoolIndex_memoryPrecisionHigh_choice & trialBoolIndex_metaLow_choice & (memoryPrecision_trialLevel>temp_highThreshold_memoryPrecision);

%     % neuron label + partial behavior label (Choice-memory/Choice-offload), bad
%     temp_trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh = ...
%         trialBoolIndex_memoryPrecisionHigh_choiceMemory & (memoryPrecision_trialLevel>temp_highThreshold_memoryPrecision);
%     temp_trialBoolIndex_memoryPrecisionLow_choiceOffloadLow = ...
%         trialBoolIndex_memoryPrecisionLow_choiceOffload & (memoryPrecision_trialLevel<temp_highThreshold_memoryPrecision);
%     temp_trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh = ...
%         trialBoolIndex_memoryPrecisionLow_choiceMemory & (memoryPrecision_trialLevel<temp_lowThreshold_memoryPrecision);
%     temp_trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow = ...
%         trialBoolIndex_memoryPrecisionHigh_choiceOffload  & (memoryPrecision_trialLevel>temp_highThreshold_memoryPrecision);

% neuron label + partial behavior label (correct/error + Choice-memory/Choice-offload)
%     temp_trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh = ...
%         choiceMemoryCorrectBoolIndex' & trialBoolIndex_memoryPrecisionHigh_choiceMemory & (memoryPrecision_trialLevel>temp_highThreshold_memoryPrecision);
%     temp_trialBoolIndex_memoryPrecisionLow_choiceOffloadLow = ...
%         trialBoolIndex_memoryPrecisionLow_choiceOffload & (memoryPrecision_trialLevel<temp_highThreshold_memoryPrecision);
%     temp_trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh = ...
%         choiceMemoryErrorBoolIndex' & trialBoolIndex_memoryPrecisionLow_choiceMemory & (memoryPrecision_trialLevel<temp_lowThreshold_memoryPrecision);
%     temp_trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow = ...
%         trialBoolIndex_memoryPrecisionHigh_choiceOffload  & (memoryPrecision_trialLevel>temp_highThreshold_memoryPrecision);


% end

% % a1 = score_trialHistory(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory,:);
% % a2 = score_trialHistory(trialBoolIndex_memoryPrecisionLow_choiceOffload,:);
% a1 = score_trialHistory(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh,:);
% a2 = score_trialHistory(trialBoolIndex_memoryPrecisionLow_choiceOffloadLow,:);
%
% % a3 = score_trialHistory(trialBoolIndex_memoryPrecisionLowError_choiceMemory,:);
% % a4 = score_trialHistory(trialBoolIndex_memoryPrecisionHigh_choiceOffload,:);
% a3 = score_trialHistory(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh,:);
% a4 = score_trialHistory(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow,:);


a1 = score_trialHistory(temp_trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh,:);
a2 = score_trialHistory(temp_trialBoolIndex_memoryPrecisionLow_choiceOffloadLow,:);

a3 = score_trialHistory(temp_trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh,:);
a4 = score_trialHistory(temp_trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow,:);

a12 = score_trialHistory(temp_trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh|temp_trialBoolIndex_memoryPrecisionLow_choiceOffloadLow,:);

%
% [~,p_a3_a1] = ttest2(a3,a1);
% [~,p_a1_a2] = ttest2(a1,a2);
% [~,p_a2_a4] = ttest2(a2,a4);

% [~,p_a2_a3] = ttest2(a2,a3);
% [~,p_a1_a4] = ttest2(a1,a4);
[~,p_a2_a3] = ttest2(a2,a3,'Tail','left');
[~,p_a1_a4] = ttest2(a1,a4,'Tail','right');

if false
    [~,p_a2_a3] = ttest2(a2,a3,'Tail','left'); %#ok<*UNRCH>
    [~,p_a1_a4] = ttest2(a1,a4,'Tail','right');
end


x = []; %#ok<*NASGU>
y = [];
% for tempi=1:seq_length_max
%     x = [x;tempi*ones(size(temp_glm_std_length1234,1),1)];    %#ok<*AGROW>
%     y = [y;temp_glm_std_length1234(:,tempi)];
% end

% overMismatch-->highMatch-->lowMatch-->underMismatch
x = [1*ones(size(a3,1),1);2*ones(size(a1,1),1);3*ones(size(a2,1),1);4*ones(size(a4,1),1)];
y = [a3;a1;a2;a4];

% % overMismatch-->highMatch-->underMismatch-->lowMatch
% x = [1*ones(size(a3,1),1);2*ones(size(a1,1),1);3*ones(size(a4,1),1);4*ones(size(a2,1),1)];
% y = [a3;a1;a4;a2];


% % overMismatch-->highMatch-->lowMatch
% x = [1*ones(size(a3,1),1);2*ones(size(a1,1),1);3*ones(size(a2,1),1)];
% y = [a3;a1;a2];

%cftool(x,y);
p_linear = nan(1,length(range_trialHistory));
for tempi=1:length(range_trialHistory)
    temp_mdl = fitglm(x,y(:,tempi),'linear');
    p_linear(tempi) = temp_mdl.Coefficients.pValue(2);
end

[M,I]=min(p_linear);




% overMismatch-->Match-->underMismatch
x = [1*ones(size(a3,1),1);2*ones(size(a12,1),1);3*ones(size(a4,1),1)];
y = [a3;a12;a4];

p_linear_B = nan(1,length(range_trialHistory));
for tempi=1:length(range_trialHistory)
    temp_mdl = fitglm(x,y(:,tempi),'linear');
    p_linear_B(tempi) = temp_mdl.Coefficients.pValue(2);
end

[M_B,I_B]=min(p_linear_B);



a1_mean = mean(a1,1,'omitnan');
a2_mean = mean(a2,1,'omitnan');
a3_mean = mean(a3,1,'omitnan');
a4_mean = mean(a4,1,'omitnan');

a1_sem = std(a1,1,1,'omitnan')./sqrt(size(a1,1));
a2_sem = std(a2,1,1,'omitnan')./sqrt(size(a2,1));
a3_sem = std(a3,1,1,'omitnan')./sqrt(size(a3,1));
a4_sem = std(a4,1,1,'omitnan')./sqrt(size(a4,1));


a1_median = median(a1,1,'omitnan');
a2_median = median(a2,1,'omitnan');
a3_median = median(a3,1,'omitnan');
a4_median = median(a4,1,'omitnan');


%% all trials (baseline)
trialBoolIndex_metaHigh_choice_baseline;
trialBoolIndex_metaLow_choice_baseline;

a5 = score_trialHistory(trialBoolIndex_metaHigh_choice_baseline,:);
a6 = score_trialHistory(trialBoolIndex_metaLow_choice_baseline,:);

[~,p_a5_a6] = ttest2(a5,a6);

a5_mean = mean(a5,1,'omitnan');
a6_mean = mean(a6,1,'omitnan');

a5_sem = std(a5,1,1,'omitnan')./sqrt(size(a5,1));
a6_sem = std(a6,1,1,'omitnan')./sqrt(size(a6,1));

%% all trials (behavior evidence)
choiceMemoryBoolIndex;
choiceOffloadBoolIndex;

a7 = score_trialHistory(choiceMemoryBoolIndex,:);
a8 = score_trialHistory(choiceOffloadBoolIndex,:);

[~,p_a7_a8] = ttest2(a7,a8);

a7_mean = mean(a7,1,'omitnan');
a8_mean = mean(a8,1,'omitnan');

a7_sem = std(a7,1,1,'omitnan')./sqrt(size(a7,1));
a8_sem = std(a8,1,1,'omitnan')./sqrt(size(a8,1));


%% all trials (delay1)
trialBoolIndex_metaHigh_choice;
trialBoolIndex_metaLow_choice;

a9 = score_trialHistory(trialBoolIndex_metaHigh_choice,:);
a10 = score_trialHistory(trialBoolIndex_metaLow_choice,:);

[~,p_a9_a10] = ttest2(a9,a10);

a9_mean = mean(a9,1,'omitnan');
a10_mean = mean(a10,1,'omitnan');

a9_sem = std(a9,1,1,'omitnan')./sqrt(size(a9,1));
a10_sem = std(a10,1,1,'omitnan')./sqrt(size(a10,1));


%% meta_baseline and meta_delay1 GLM
tempNANBoolIndex1 = isnan(meta_trialLevel_baseline);
tempNANBoolIndex2 = isnan(memoryPrecision_trialLevel);
tempNANBoolIndex3 = isnan(meta_trialLevel);

tempNANBoolIndex123 = tempNANBoolIndex1 | tempNANBoolIndex2 | tempNANBoolIndex3;

tempNONNANBoolIndex123 = ~tempNANBoolIndex123;

tempNONNAN_meta_trialLevel_baseline = meta_trialLevel_baseline(tempNONNANBoolIndex123);
tempNONNAN_memoryPrecision_trialLevel = memoryPrecision_trialLevel(tempNONNANBoolIndex123);
tempNONNAN_meta_trialLevel = meta_trialLevel(tempNONNANBoolIndex123);

% Case A
x = tempNONNAN_memoryPrecision_trialLevel;
y = tempNONNAN_meta_trialLevel;
temp_mdl_caseA = fitglm(x,y,'linear');
r2_caseA = temp_mdl_caseA.Rsquared.Adjusted;
beta0_caseA = temp_mdl_caseA.Coefficients.Estimate(1);
beta1_caseA = temp_mdl_caseA.Coefficients.Estimate(2);

% Case B
x = tempNONNAN_meta_trialLevel_baseline;
y = tempNONNAN_meta_trialLevel;
temp_mdl_caseB = fitglm(x,y,'linear');
r2_caseB = temp_mdl_caseB.Rsquared.Adjusted;
beta0_caseB = temp_mdl_caseB.Coefficients.Estimate(1);
beta1_caseB = temp_mdl_caseB.Coefficients.Estimate(2);

% Case C
x = [tempNONNAN_memoryPrecision_trialLevel,tempNONNAN_meta_trialLevel_baseline];
y = tempNONNAN_meta_trialLevel;
temp_mdl_caseC = fitglm(x,y,'linear');
r2_caseC = temp_mdl_caseC.Rsquared.Adjusted;
beta0_caseC = temp_mdl_caseC.Coefficients.Estimate(1);
beta1_caseC= temp_mdl_caseC.Coefficients.Estimate(2);
beta2_caseC= temp_mdl_caseC.Coefficients.Estimate(3);

% Seq level (control)
memoryPrecision_trialLevel;
meta_trialLevel;
seqIndex;

memoryPrecision_seqLevel = zeros(sum(numSeq(1:3)),1);
meta_seqLevel = zeros(sum(numSeq(1:3)),1);
for tempi=1:sum(numSeq(1:3))
    %     temptempBoolIndex = seqIndex==tempi;
    temptempBoolIndex = seqIndex_response==tempi;
    
    temptempBoolIndex2 = temptempBoolIndex;
    %     temptempBoolIndex2 = temptempBoolIndex & choiceBoolIndex;
    %     temptempBoolIndex2 = temptempBoolIndex & ~choiceBoolIndex;
    
    memoryPrecision_seqLevel(tempi) = mean(memoryPrecision_trialLevel(temptempBoolIndex2),'omitnan');
    meta_seqLevel(tempi) = mean(meta_trialLevel(temptempBoolIndex2),'omitnan');
end
memoryPrecision_seqLevel;
meta_seqLevel;

x = memoryPrecision_seqLevel;
y = meta_seqLevel;
temp_mdl_caseA_seq = fitglm(x,y,'linear');
r2_caseA_seq = temp_mdl_caseA_seq.Rsquared.Adjusted;
beta0_caseA_seq = temp_mdl_caseA_seq.Coefficients.Estimate(1);
beta1_caseA_seq = temp_mdl_caseA_seq.Coefficients.Estimate(2);
% r2_caseA_seq



%% Baseline meta-memory & 4 conditions
meta_trialLevel_baseline;

temp_trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh; % High-match
temp_trialBoolIndex_memoryPrecisionLow_choiceOffloadLow;% Low-match
temp_trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh;% Over-mismatch
temp_trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow;% Under-mismatch


z1 = meta_trialLevel_baseline(temp_trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh,:);
z2 = meta_trialLevel_baseline(temp_trialBoolIndex_memoryPrecisionLow_choiceOffloadLow,:);
z3 = meta_trialLevel_baseline(temp_trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh,:);
z4 = meta_trialLevel_baseline(temp_trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow,:);

z12 = [z1;z2];

[~,p_z2_z3] = ttest2(z2,z3);
[~,p_z1_z4] = ttest2(z1,z4);


% overMismatch-->highMatch-->lowMatch-->underMismatch
x = [1*ones(size(z3,1),1);2*ones(size(z1,1),1);3*ones(size(z2,1),1);4*ones(size(z4,1),1)];
y = [z3;z1;z2;z4];
temp_mdl = fitglm(x,y,'linear');
p_linear_baseline = temp_mdl.Coefficients.pValue(2);


% overMismatch-->match-->underMismatch
x = [1*ones(size(z3,1),1);2*ones(size(z12,1),1);3*ones(size(z4,1),1)];
y = [z3;z12;z4];
temp_mdl = fitglm(x,y,'linear');
p_linear_baseline3 = temp_mdl.Coefficients.pValue(2);




z5 = meta_trialLevel_baseline(trialBoolIndex_metaHigh_choice,:);
z6 = meta_trialLevel_baseline(trialBoolIndex_metaLow_choice,:);

[~,p_z5_z6] = ttest2(z5,z6);

z5_mean = mean(z5,1,'omitnan');
z6_mean = mean(z6,1,'omitnan');

z5_sem = std(z5,1,1,'omitnan')./sqrt(size(z5,1));
z6_sem = std(z6,1,1,'omitnan')./sqrt(size(z6,1));



z7 = meta_trialLevel_baseline(trialBoolIndex_memoryPrecisionHigh_choice,:);
z8 = meta_trialLevel_baseline(trialBoolIndex_memoryPrecisionLow_choice,:);

[~,p_z7_z8] = ttest2(z7,z8);

z7_mean = mean(z7,1,'omitnan');
z8_mean = mean(z8,1,'omitnan');

z7_sem = std(z7,1,1,'omitnan')./sqrt(size(z7,1));
z8_sem = std(z8,1,1,'omitnan')./sqrt(size(z8,1));


%% Offload trial VS. Low-strength mismatch trials
za5 = meta_trialLevel_baseline(trialBoolIndex_memoryPrecisionLow_metaHigh_choiceMemoryError,:);
za6 = meta_trialLevel_baseline(choiceOffloadBoolIndex,:);
% [~,p_za5_za6] = ttest2(za5,za6);

zb5 = meta_trialLevel(trialBoolIndex_memoryPrecisionLow_metaHigh_choiceMemoryError,:);
zb6 = meta_trialLevel(choiceOffloadBoolIndex,:);
% [~,p_zb5_zb6] = ttest2(zb5,zb6);

zc5 = memoryPrecision_trialLevel(trialBoolIndex_memoryPrecisionLow_metaHigh_choiceMemoryError,:);
zc6 = memoryPrecision_trialLevel(choiceOffloadBoolIndex,:);
% [~,p_zc5_zc6] = ttest2(zc5,zc6);


% [~,p_za5_za6] = ttest2(za5,za6);
% [~,p_zb5_zb6] = ttest2(zb5,zb6);
% [~,p_zc5_zc6] = ttest2(zc5,zc6);

[~,p_za5_za6] = ttest2(za5,za6,'tail','right');
[~,p_zb5_zb6] = ttest2(zb5,zb6,'tail','right');
[~,p_zc5_zc6] = ttest2(zc5,zc6,'tail','right');

% p_za5_za6 = kruskalwallis([za5;za6],[zeros(size(za5));ones(size(za6))],'off');
% p_zb5_zb6 = kruskalwallis([zb5;zb6],[zeros(size(zb5));ones(size(zb6))],'off');
% p_zc5_zc6 = kruskalwallis([zc5;zc6],[zeros(size(zc5));ones(size(zc6))],'off');

% p_zc5_zc6 = ranksum(zc5,zc6);


% za5_minus_za6 = za5-za6;
% zb5_minus_zb6 = zb5-zb6;
% zc5_minus_zc6 = zc5-zb6;

zb5_minus_za5 = zb5-za5;
zb6_minus_za6 = zb6-za6;
[~,p_zab5_zab6] = ttest2(zb5_minus_za5,zb6_minus_za6);
fprintf('zb5_minus_za5 = %.3f, zb6_minus_za6 = %.3f.\n',mean(zb5_minus_za5,'omitnan'),mean(zb6_minus_za6,'omitnan'));

za5_minus_zc5 = za5-zc5;
za6_minus_zc6 = za6-zc6;
[~,p_zac5_zac6] = ttest2(za5_minus_zc5,za6_minus_zc6);
fprintf('za5_minus_zc5 = %.3f, za6_minus_zc6 = %.3f.\n',mean(za5_minus_zc5,'omitnan'),mean(za6_minus_zc6,'omitnan'));


%% Plot
% close all
if if_plot == 1
    fig = figure('Name','asd','NumberTitle','off');
    % set(gcf,'Position',[10 50 1450 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 50 720 240]);
    % t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    %t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    t = tiledlayout(1,4,'TileSpacing','tight','Padding','Compact');
    t.Title.String = sprintf('Trial history, %s',FOVName_currentFOV2);
    t.Title.FontSize = 12;
    t.Title.Interpreter = 'none';
    
    
    %% Plot behavior evidence
    nexttile
    
    h_line = [];
    x = range_trialHistory;
    
    y = a7_mean;
    h_line = [h_line plot(x,y,'color',color_highMeta,'linewidth',1)]; %#ok<*AGROW>
    hold on
    y_sem = a7_sem;
    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_highMeta,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
    hold on
    
    y = a8_mean;
    h_line = [h_line plot(x,y,'color',color_lowMeta,'linewidth',1)]; %#ok<*AGROW>
    hold on
    y_sem = a8_sem;
    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_lowMeta,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
    hold on
    
    
    
    % [y_min,y_max] = bounds([a7_mean,a8_mean]);
    [y_min,y_max] = bounds([a7_mean-a7_sem,a7_mean+a7_sem,a8_mean-a8_sem,a8_mean+a8_sem]);
    
    scatter(range_trialHistory(p_a7_a8<0.05),y_max+(y_max-y_min)*0.05,20,[0 0 0],'*');
    hold on
    
    le = legend(h_line,sprintf('Choice-memory'),...
        sprintf('Choice-offload'),...
        'Location','southeast','fontsize',6);
    le.ItemTokenSize = ones(1,2)*10;
    
    if if_xAxia_normal0_log1 == 1
        set(gca,'xscale','log')
    end
    % xticks([1:8,10,15,20,30,range_trialHistory(end)]);
    %xticks([1:8,10,range_trialHistory(end)]);
    xticks([1:range_trialHistory(end)]);
    % xticks([range_trialHistory(1)-1:5:range_trialHistory(end)]);
    xtickangle(0);
    
    set(gca,'linewidth',1.5)
    xlim([range_trialHistory(1)-1 range_trialHistory(end)+1]);
    ylim([y_min-(y_max-y_min)*0.2 y_max+(y_max-y_min)*0.1]);%0.1
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    % xticks([range_trialHistory(1)-1:5:range_trialHistory(end)]);
    xlabel('Trial history', 'FontSize', 12, 'FontWeight', 'bold');
    if if_scoreType_cumulative0_meanS1_singleTrial2 == 0
        ylabel('Cumulative score', 'FontSize', 12, 'FontWeight', 'bold');
    elseif if_scoreType_cumulative0_meanS1_singleTrial2 == 1
        %ylabel('Cumulative score (mean)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Mean reward', 'FontSize', 12, 'FontWeight', 'bold');
    elseif if_scoreType_cumulative0_meanS1_singleTrial2 == 2
        ylabel('Reward', 'FontSize', 12, 'FontWeight', 'bold');
    end
    % temp_title = title(sprintf('Trial history, %s',FOVName_currentFOV2),...
    %     'FontSize',14,'FontWeight','bold');
    temp_title = title(sprintf('Behavior evidence'),...
        'FontSize',12,'FontWeight','bold');
    temp_title.Interpreter = 'none';
    
    
    %% Plot baseline
    nexttile
    
    h_line = [];
    x = range_trialHistory;
    
    y = a5_mean;
    h_line = [h_line plot(x,y,'color',color_highMeta,'linewidth',1)]; %#ok<*AGROW>
    hold on
    y_sem = a5_sem;
    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_highMeta,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
    hold on
    
    y = a6_mean;
    h_line = [h_line plot(x,y,'color',color_lowMeta,'linewidth',1)]; %#ok<*AGROW>
    hold on
    y_sem = a6_sem;
    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_lowMeta,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
    hold on
    
    % [y_min,y_max] = bounds([a5_mean,a6_mean]);
    [y_min,y_max] = bounds([a5_mean-a5_sem,a5_mean+a5_sem,a6_mean-a6_sem,a6_mean+a6_sem]);
    
    scatter(range_trialHistory(p_a5_a6<0.05),y_max+(y_max-y_min)*0.05,20,[0 0 0],'*');
    hold on
    
    le = legend(h_line,sprintf('High-meta'),...
        sprintf('Low-meta'),...
        'Location','southeast','fontsize',6);
    le.ItemTokenSize = ones(1,2)*10;
    
    
    set(gca,'linewidth',1.5)
    if if_xAxia_normal0_log1 == 1
        set(gca,'xscale','log')
    end
    % xticks([1:8,10,15,20,30,range_trialHistory(end)]);
    %xticks([1:8,10,range_trialHistory(end)]);
    xticks([1:range_trialHistory(end)]);
    % xticks([range_trialHistory(1)-1:5:range_trialHistory(end)]);
    xtickangle(0);
    
    xlim([range_trialHistory(1)-1 range_trialHistory(end)+1]);
    ylim([y_min-(y_max-y_min)*0.2 y_max+(y_max-y_min)*0.1]);%0.1
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Trial history', 'FontSize', 12, 'FontWeight', 'bold');
    % if if_scoreType_cumulative0_meanS1_singleTrial2 == 0
    %     ylabel('Cumulative score', 'FontSize', 12, 'FontWeight', 'bold');
    % elseif if_scoreType_cumulative0_meanS1_singleTrial2 == 1
    %     ylabel('Cumulative score (mean)', 'FontSize', 12, 'FontWeight', 'bold');
    % end
    % temp_title = title(sprintf('Trial history, %s',FOVName_currentFOV2),...
    %     'FontSize',14,'FontWeight','bold');
    temp_title = title(sprintf('Baseline-based'),...
        'FontSize',12,'FontWeight','bold');
    temp_title.Interpreter = 'none';
    
    %% Plot delay1
    nexttile
    
    h_line = [];
    x = range_trialHistory;
    
    y = a3_mean;
    h_line = [h_line plot(x,y,'color',color_choiceMemoryLow,'linewidth',1)]; %#ok<*AGROW>
    hold on
    % y_sem = a3_sem;
    % patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
    % hold on
    
    y = a1_mean;
    h_line = [h_line plot(x,y,'color',color_choiceMemoryHigh,'linewidth',1)]; %#ok<*AGROW>
    hold on
    % y_sem = a1_sem;
    % patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemoryHigh,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
    % hold on
    
    y = a2_mean;
    h_line = [h_line plot(x,y,'color',color_choiceOffloadLow,'linewidth',1)]; %#ok<*AGROW>
    hold on
    % y_sem = a2_sem;
    % patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadLow,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
    % hold on
    
    y = a4_mean;
    h_line = [h_line plot(x,y,'color',color_choiceOffloadHigh,'linewidth',1)]; %#ok<*AGROW>
    hold on
    % y_sem = a4_sem;
    % patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffloadHigh,'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
    % hold on
    
    
    [y_min,y_max] = bounds([a1_mean,a2_mean,a3_mean,a4_mean]);
    
    
    scatter(range_trialHistory(p_linear<0.05),y_max+(y_max-y_min)*0.05,20,[0 0 0],'*');
    hold on
    
    le = legend(h_line,sprintf('Over-mismatch'),...
        sprintf('High-match'),...
        sprintf('Low-match'),...
        sprintf('Under-mismatch'),...
        'Location','southeast','fontsize',6);%9
    le.ItemTokenSize = ones(1,2)*10;
    
    set(gca,'linewidth',1.5)
    if if_xAxia_normal0_log1 == 1
        set(gca,'xscale','log')
    end
    % xticks([1:8,10,15,20,30,range_trialHistory(end)]);
    %xticks([1:8,10,range_trialHistory(end)]);
    xticks([1:range_trialHistory(end)]);
    % xticks([range_trialHistory(1)-1:5:range_trialHistory(end)]);
    xtickangle(0);
    
    xlim([range_trialHistory(1)-1 range_trialHistory(end)+1]);
    ylim([y_min-(y_max-y_min)*0.2 y_max+(y_max-y_min)*0.1]);%0.1
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    % xticks([range_trialHistory(1)-1:5:range_trialHistory(end)]);
    xlabel('Trial history', 'FontSize', 12, 'FontWeight', 'bold');
    % if if_scoreType_cumulative0_meanS1_singleTrial2 == 0
    %     ylabel('Cumulative score', 'FontSize', 12, 'FontWeight', 'bold');
    % elseif if_scoreType_cumulative0_meanS1_singleTrial2 == 1
    %     ylabel('Cumulative score (mean)', 'FontSize', 12, 'FontWeight', 'bold');
    % end
    % temp_title = title(sprintf('Trial history, %s',FOVName_currentFOV2),...
    %     'FontSize',14,'FontWeight','bold');
    temp_title = title(sprintf('Delay1-based'),...
        'FontSize',12,'FontWeight','bold');
    temp_title.Interpreter = 'none';
    
    
    
    %% Plot history weight
    nexttile
    
    x = range_trialHistory;
    
    y = weight_trialHistory;
    plot(x,y,'color',[0.25 0.25 0.25],'linewidth',1);
    hold on
    
    
    [y_min,y_max] = bounds(y);
    %[y_min,y_max] = bounds([a7_mean-a7_sem,a7_mean+a7_sem,a8_mean-a8_sem,a8_mean+a8_sem]);
    
    if if_xAxia_normal0_log1 == 1
        set(gca,'xscale','log')
    end
    % xticks([1:8,10,15,20,30,range_trialHistory(end)]);
    %xticks([1:8,10,range_trialHistory(end)]);
    xticks([1:range_trialHistory(end)]);
    % xticks([range_trialHistory(1)-1:5:range_trialHistory(end)]);
    xtickangle(0);
    
    set(gca,'linewidth',1.5)
    xlim([range_trialHistory(1)-1 range_trialHistory(end)+1]);
    ylim([y_min y_max+(y_max-y_min)*0.1]);%0.1
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    % xticks([range_trialHistory(1)-1:5:range_trialHistory(end)]);
    xlabel('Trial history', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Weight', 'FontSize', 12, 'FontWeight', 'bold');
    % temp_title = title(sprintf('Trial history, %s',FOVName_currentFOV2),...
    %     'FontSize',14,'FontWeight','bold');
    temp_title = title(sprintf('History weight'),...
        'FontSize',12,'FontWeight','bold');
    temp_title.Interpreter = 'none';
    
    
    
    
    
    
    
    fig = figure('Name','asd','NumberTitle','off');
    set(gcf,'Position',[10 400 720 240]);
    t = tiledlayout(1,4,'TileSpacing','tight','Padding','Compact');
    t.Title.String = sprintf('Trial history, %s',FOVName_currentFOV2);
    t.Title.FontSize = 12;
    t.Title.Interpreter = 'none';
    
    
    %% Plot history weight
    nexttile
    
    x = range_trialHistory;
    
    y = weight_trialHistory;
    plot(x,y,'color',[0.25 0.25 0.25],'linewidth',1);
    hold on
    
    
    [y_min,y_max] = bounds(y);
    %[y_min,y_max] = bounds([a7_mean-a7_sem,a7_mean+a7_sem,a8_mean-a8_sem,a8_mean+a8_sem]);
    
    if if_xAxia_normal0_log1 == 1
        set(gca,'xscale','log')
    end
    % xticks([1:8,10,15,20,30,range_trialHistory(end)]);
    %xticks([1:8,10,range_trialHistory(end)]);
    xticks([1:range_trialHistory(end)]);
    % xticks([range_trialHistory(1)-1:5:range_trialHistory(end)]);
    xtickangle(0);
    
    set(gca,'linewidth',1.5)
    xlim([range_trialHistory(1)-1 range_trialHistory(end)+1]);
    ylim([y_min y_max+(y_max-y_min)*0.1]);%0.1
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    % xticks([range_trialHistory(1)-1:5:range_trialHistory(end)]);
    xlabel('Trial history', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Weight', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('History weight'),...
        'FontSize',12,'FontWeight','bold');
    temp_title.Interpreter = 'none';
    
    
    %% Plot behavior evidence
    nexttile
    
    
    temp_p = p_a7_a8(end);
    temp_1_raw = a7(:,end);
    temp_2_raw = a8(:,end);
    
    temp_1 = temp_1_raw(~isnan(temp_1_raw));
    temp_2 = temp_2_raw(~isnan(temp_2_raw));
    
    temp1_SEM = a7_sem(end);
    temp2_SEM = a8_sem(end);
    
    
    temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM]);
    temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
    
    temp_bar = bar([1 2], [mean(temp_1) mean(temp_2)], ...
        'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1);
    hold on
    a = 1;
    temp_bar.CData(1,:) = [1 1 1]*0.5;%0.3
    temp_bar.CData(2,:) = [1 1 1]*0.5;
    
    errorbar([1 2], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 8);
    
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
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
    set(gca, 'FontSize', 10)
    
    xtl = ["Choice-memory"; "Choice-offload"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.375;%0.56,0.4
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',20,'fontsize',10);%25
    set(gca,'xticklabel','');
    
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Mean reward', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('Behavior evidence'),'FontSize',12,'FontWeight','bold');
    temp_title.Interpreter = 'none';
    
    
    %% Plot baseline
    nexttile
    
    temp_p = p_a5_a6(end);
    temp_1_raw = a5(:,end);
    temp_2_raw = a6(:,end);
    
    temp_1 = temp_1_raw(~isnan(temp_1_raw));
    temp_2 = temp_2_raw(~isnan(temp_2_raw));
    
    temp1_SEM = a5_sem(end);
    temp2_SEM = a6_sem(end);
    
    
    temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM]);
    temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
    
    temp_bar = bar([1 2], [mean(temp_1) mean(temp_2)], ...
        'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1);
    hold on
    a = 1;
    temp_bar.CData(1,:) = [1 1 1]*0.5;%0.3
    temp_bar.CData(2,:) = [1 1 1]*0.5;
    
    errorbar([1 2], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 8);
    
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
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
    set(gca, 'FontSize', 10)
    
    xtl = ["High-meta"; "Low-meta"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.425;%0.56,0.4
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',20,'fontsize',10);%25
    set(gca,'xticklabel','');
    
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Mean reward', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('Baseline-based'),'FontSize',12,'FontWeight','bold');
    temp_title.Interpreter = 'none';
    
    
    %% Plot delay1
    nexttile
    
    
    temp_p = p_linear(end);
    temp_1_raw = a3(:,end);
    temp_2_raw = a1(:,end);
    temp_3_raw = a2(:,end);
    temp_4_raw = a4(:,end);
    
    temp_1 = temp_1_raw(~isnan(temp_1_raw));
    temp_2 = temp_2_raw(~isnan(temp_2_raw));
    temp_3 = temp_3_raw(~isnan(temp_3_raw));
    temp_4 = temp_4_raw(~isnan(temp_4_raw));
    
    temp1_SEM = a3_sem(end);
    temp2_SEM = a1_sem(end);
    temp3_SEM = a2_sem(end);
    temp4_SEM = a4_sem(end);
    
    temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM,mean(temp_3)-temp3_SEM,mean(temp_4)-temp4_SEM]);
    temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM,mean(temp_4)+temp4_SEM]);
    
    temp_bar = bar([1 2 3 4], [mean(temp_1) mean(temp_2) mean(temp_3) mean(temp_4)], ...
        'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1);
    hold on
    a = 1;
    temp_bar.CData(1,:) = [1 1 1]*0.5;%0.3
    temp_bar.CData(2,:) = [1 1 1]*0.5;
    temp_bar.CData(3,:) = [1 1 1]*0.5;%0.3
    temp_bar.CData(4,:) = [1 1 1]*0.5;
    
    errorbar([1 2 3 4], [mean(temp_1) mean(temp_2) mean(temp_3) mean(temp_4)],...
        [temp1_SEM temp2_SEM temp3_SEM temp4_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 8);
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(2.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 4.65])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
    set(gca, 'FontSize', 10)
    
    xtl = ["Over-mismatch"; "High-match";"Low-match"; "Under-mismatch"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;%0.56,0.4
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
    set(gca,'xticklabel','');
    
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Mean reward', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('Delay1-based'),'FontSize',12,'FontWeight','bold');
    temp_title.Interpreter = 'none';
    
    
    
    
    
    
    
    if if_noFit0_fitScore1_fitHistroyW2 == 2
        
        %         fig = figure('Name','asd','NumberTitle','off');
        %         %set(gcf,'Position',[810 400 240*1.15 240*0.7*0.85]);
        %         set(gcf,'Position',[810 400 240 143*0.97]);
        %         t = tiledlayout(1,1,'TileSpacing','tight','Padding','Compact');
        %
        %         %% Plot history weight fitting
        %         nexttile
        %
        %         x = normStdRange;
        %         y = temp_r2_fitting;
        %
        %         plot(x,y,'color',[0.25 0.25 0.25],'linewidth',1);
        %         hold on
        %
        %         [x_min,x_max] = bounds(x);
        %         [y_min,y_max] = bounds(y);
        %
        %         plot([1 1].*normStd_optimal,[y_min-(y_max-y_min)*0.1 temp_r2_optimal],':','color',[0.25 0.25 0.25],'linewidth',1);
        %         hold on
        %
        %         plot([x_min-(x_max-x_min)*0.05 normStd_optimal],[1 1].*temp_r2_optimal,':','color',[0.25 0.25 0.25],'linewidth',1);
        %         hold on
        %
        %
        %         %set(gca,'yscale','log')
        %
        %         set(gca,'linewidth',1.5)
        %         xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05]);%
        %         ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);%0.1
        %         set(gca, 'FontSize', 8)
        %         set(gca,'box','off');% 取消右、上边框
        %         xlabel('std of norm distribution', 'FontSize', 12, 'FontWeight', 'bold');
        %         ylabel('r2', 'FontSize', 12, 'FontWeight', 'bold');
        %         temp_title = title(sprintf('History weight fitting'),...
        %             'FontSize',12,'FontWeight','bold');
        %         temp_title.Interpreter = 'none';
        %         temp_title.Position(1) = temp_title.Position(1) + 5;
        
        if if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 1
            temp_r2_optimal;
            AUROC_optimal;
            end_optimal;
            normStd_optimal;
            temp_r2_fitting_2d;
            AUROC_fitting_2d;
            row;
            col;
            
            
            fig = figure('Name','asd','NumberTitle','off');
            %set(gcf,'Position',[810 200 370 300]);
            set(gcf,'Position',[810 200 370*0.8*1.05*1.05 300*0.8*0.95*1.05]);
            
            t = tiledlayout(1,1,'TileSpacing','tight','Padding','Compact');
            
            C = AUROC_fitting_2d;
            %imagesc([1 size(C,2)],[1 size(C,1)],C);
            imagesc([normStdRange(1) normStdRange(end)],[1 size(C,1)],C);
            hold on
            
            plot([normStd_optimal normStd_optimal],[1 size(C,1)],':',...
                'LineWidth',1,'color',[1 1 1]*0);
            hold on
            plot([1 size(C,2)],[end_optimal end_optimal],':',...
                'LineWidth',1,'color',[1 1 1]*0);
            hold on
            
            
            set(gca,'Ydir','normal') %reverse
            xtickangle(0);
            set(gca, 'FontSize', 9) %12
            
            c = colorbar('FontSize',8);
            
            colormap parula
            
            xlabel('Std of norm distribution', 'FontSize', 10, 'FontWeight', 'normal');
            ylabel('Trial range', 'FontSize', 10, 'FontWeight', 'normal');
            
            temp_title = title(sprintf('AUROC'),'FontSize',10,'FontWeight','normal');
            temp_title.Interpreter = 'none';
            
            
        elseif if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 0
            
            fig = figure('Name','asd','NumberTitle','off');
            %set(gcf,'Position',[810 400 240*1.15 240*0.7*0.85]);
            %set(gcf,'Position',[810 600 240 143*1.02]);
            %set(gcf,'Position',[810 600 240 143*1.02*0.88]);
            %set(gcf,'Position',[810 600 240*1.26 143*1.02*0.88]);
            set(gcf,'Position',[810 600 240*1.26*0.7 143*1.02*0.88]);
            t = tiledlayout(1,1,'TileSpacing','tight','Padding','Compact');
            
            %% Plot history weight fitting
            nexttile
            
            x = normStdRange;
            y = AUROC_fitting;
            
            plot(x,y,'color',[0.25 0.25 0.25],'linewidth',1);
            hold on
            
            [x_min,x_max] = bounds(x);
            [y_min,y_max] = bounds(y);
            y_min = min(y_min,0.5);
            
            plot([1 1].*normStd_optimal,[y_min-(y_max-y_min)*0.1 AUROC_optimal],':','color',[0.25 0.25 0.25],'linewidth',1);
            hold on
            
            plot([x_min-(x_max-x_min)*0.05 normStd_optimal],[1 1].*AUROC_optimal,':','color',[0.25 0.25 0.25],'linewidth',1);
            hold on
            
            plot([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05],[1 1].*0.5,':','color',[0.25 0.25 0.25],'linewidth',1);
            hold on
            
            %set(gca,'yscale','log')
            
            yticks([0.5 0.56]);
            
            set(gca,'linewidth',1.5)
            xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05]);%
            ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);%0.1
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            xlabel('Std of norm distribution', 'FontSize', 10, 'FontWeight', 'normal');
            ylabel('AUROC', 'FontSize', 10, 'FontWeight', 'normal');
            temp_title = title(sprintf('History weight fitting'),...
                'FontSize',10,'FontWeight','normal');
            temp_title.Interpreter = 'none';
            %temp_title.Position(1) = temp_title.Position(1) + 5;
            
        end
        
        fig = figure('Name','asd','NumberTitle','off');
        %set(gcf,'Position',[810 400 240*1.15 240*0.7*0.85]);
        %set(gcf,'Position',[1110 400 240*0.8*0.85 143*1.02]);
        %set(gcf,'Position',[1110 400 240*0.8*0.85 143*1.02*0.88]);
        %set(gcf,'Position',[1110 400 240*0.8*0.85*1.85 143*1.02*0.88]);
        %set(gcf,'Position',[1110 400 240*0.8*0.85*1.85*0.7 143*1.02*0.88]);
        set(gcf,'Position',[1110 400 240*0.8*0.85*1.85*0.7*1.5*0.95 143*1.02*0.88*1.2*1.2]);
        t = tiledlayout(1,1,'TileSpacing','tight','Padding','Compact');
        
        %% Plot history weight
        nexttile
        
        x = range_trialHistory;
        y = weight_trialHistory;
        
        plot(x,y,'color',[0.25 0.25 0.25],'linewidth',1);
        hold on
        
        
        [y_min,y_max] = bounds(y);
        
        if if_xAxia_normal0_log1 == 1
            set(gca,'xscale','log')
        end
        if if_fitNormal0_fitEnd1_fitBeta2_fitExp3 ~= 1
            xticks([1:8,10,range_trialHistory(end)]);
        elseif if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 1
            xticks([1:range_trialHistory(end)]);
        end
        xtickangle(0);
        
        %yticks([0 0.1]);
        
        set(gca,'linewidth',1.5)
        xlim([range_trialHistory(1)-1 range_trialHistory(end)+1]);
        %ylim([y_min y_max+(y_max-y_min)*0.1]);%0.1
        %ylim([0 y_max+(y_max-y_min)*0.1]);%0.1
        ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);%0.1
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Trial history', 'FontSize', 10, 'FontWeight', 'normal');
        ylabel('Weight', 'FontSize', 10, 'FontWeight', 'normal');
        temp_title = title(sprintf('History weight'),...
            'FontSize',10,'FontWeight','normal');
        temp_title.Interpreter = 'none';
        
        if if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 1
            subtitle(sprintf('AUROC=%.2f, %d trials, std=%.1f',AUROC_optimal,end_optimal,normStd_optimal),...
                'FontSize',8,'FontWeight','normal');
        elseif if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 2
            subtitle(sprintf('AUROC=%.2f',AUROC_optimal),...
                'FontSize',8,'FontWeight','normal');
        end
        
        
    end
    
    
    
    %     fig = figure('Name','asd','NumberTitle','off');
    %     set(gcf,'Position',[10 400 360*0.85*1.05 500*0.8*0.9*1.05*1.03]);
    %     %t = tiledlayout(2,5,'TileSpacing','tight','Padding','Compact');
    %     t = tiledlayout(2,5,'TileSpacing','tight','Padding','Loose');
    %
    %
    %     %t.Title.String = sprintf('Trial history, %s',FOVName_currentFOV2);
    %     %t.Title.FontSize = 12;
    %     %t.Title.Interpreter = 'none';
    %
    %
    %
    %     %% Plot baseline
    %     nexttile([1 2])
    %
    %     temp_p = p_a5_a6(end);
    %     temp_1_raw = a5(:,end);
    %     temp_2_raw = a6(:,end);
    %
    %     temp_1 = temp_1_raw(~isnan(temp_1_raw));
    %     temp_2 = temp_2_raw(~isnan(temp_2_raw));
    %
    %     temp1_SEM = a5_sem(end);
    %     temp2_SEM = a6_sem(end);
    %
    %
    %     %temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM]);
    %     %temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
    %     %
    %     %temp_bar = bar([1 2], [mean(temp_1) mean(temp_2)], ...
    %     %    'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1);
    %     %hold on
    %     %a = 1;
    %     %temp_bar.CData(1,:) = [1 1 1]*0.5;%0.3
    %     %temp_bar.CData(2,:) = [1 1 1]*0.5;
    %     %
    %     %errorbar([1 2], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 8);
    %
    %
    %     temp_y_min = min([temp_1;temp_2]);
    %     temp_y_max = max([temp_1;temp_2]);
    %
    %     %temp_category = ones(size(temp_1,1)+size(temp_2,1),1);
    %     %temp_category(size(temp_1,1)+1:end) = 2;
    %     %violinplot([temp_1;temp_2],temp_category,'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
    %
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
    %     %text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %     %   'HorizontalAlignment','center');
    %     text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %         'HorizontalAlignment','center');
    %
    %     set(gca,'linewidth',1.5)
    %     xlim([0.15 2.65])
    %     %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
    %     ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
    %     set(gca, 'FontSize', 10)
    %
    %     %xtl = ["High-meta"; "Low-meta"];
    %     xtl = ["High-baseline"; "Low-baseline"];
    %     xt=get(gca,'XTick');
    %     yt=get(gca,'YTick');
    %     xtext_xp=xt;
    %     %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.425;%0.56,0.4
    %     if if_monkey_D0_Z1 == 0
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.315;%0.56,0.4
    %     elseif if_monkey_D0_Z1 == 1
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.365;%0.56,0.4
    %     end
    %     text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',20,'fontsize',10);%25
    %     set(gca,'xticklabel','');
    %
    %
    %     set(gca,'box','off');% 取消右、上边框
    %     ylabel('Mean reward', 'FontSize', 12, 'FontWeight', 'bold');
    %     %temp_title = title(sprintf('Baseline'),'FontSize',12,'FontWeight','bold');
    %     temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
    %     temp_title.Interpreter = 'none';
    %
    %
    %     %% Plot delay1
    %     nexttile([1 3])
    %
    %
    %     temp_p = p_linear(end);
    %     temp_1_raw = a3(:,end);
    %     temp_2_raw = a1(:,end);
    %     temp_3_raw = a2(:,end);
    %     temp_4_raw = a4(:,end);
    %
    %     temp_1 = temp_1_raw(~isnan(temp_1_raw));
    %     temp_2 = temp_2_raw(~isnan(temp_2_raw));
    %     temp_3 = temp_3_raw(~isnan(temp_3_raw));
    %     temp_4 = temp_4_raw(~isnan(temp_4_raw));
    %
    %     temp1_SEM = a3_sem(end);
    %     temp2_SEM = a1_sem(end);
    %     temp3_SEM = a2_sem(end);
    %     temp4_SEM = a4_sem(end);
    %
    %     %temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM,mean(temp_3)-temp3_SEM,mean(temp_4)-temp4_SEM]);
    %     %temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM,mean(temp_4)+temp4_SEM]);
    %     %
    %     %temp_bar = bar([1 2 3 4], [mean(temp_1) mean(temp_2) mean(temp_3) mean(temp_4)], ...
    %     %    'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1);
    %     %hold on
    %     %a = 1;
    %     %temp_bar.CData(1,:) = [1 1 1]*0.5;%0.3
    %     %temp_bar.CData(2,:) = [1 1 1]*0.5;
    %     %temp_bar.CData(3,:) = [1 1 1]*0.5;%0.3
    %     %temp_bar.CData(4,:) = [1 1 1]*0.5;
    %     %
    %     %errorbar([1 2 3 4], [mean(temp_1) mean(temp_2) mean(temp_3) mean(temp_4)],...
    %     %    [temp1_SEM temp2_SEM temp3_SEM temp4_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 8);
    %
    %
    %
    %     temp_y_min = min([temp_1;temp_2;temp_3;temp_4]);
    %     temp_y_max = max([temp_1;temp_2;temp_3;temp_4]);
    %
    %     %temp_category = ones(size(temp_1,1)+size(temp_2,1)+size(temp_3,1)+size(temp_4,1),1);
    %     %temp_category(size(temp_1,1)+(1:size(temp_2,1))) = 2;
    %     %temp_category(size(temp_1,1)+size(temp_2,1)+(1:size(temp_3,1))) = 3;
    %     %temp_category(size(temp_1,1)+size(temp_2,1)+size(temp_3,1)+(1:size(temp_4,1))) = 4;
    %     %
    %     %violinplot([temp_1;temp_2;temp_3;temp_4],temp_category,'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
    %
    %
    %     temp_data = [temp_1;temp_2;temp_3;temp_4];
    %
    %     g1 = repmat({'A'},length(temp_1),1);
    %     g2 = repmat({'B'},length(temp_2),1);
    %     g3 = repmat({'C'},length(temp_3),1);
    %     g4 = repmat({'D'},length(temp_4),1);
    %
    %     temp_label = [g1;g2;g3;g4];
    %
    %     temptemp_color1 = [1 1 1]*0.5;
    %     temptemp_color2 = repmat(temptemp_color1, 4, 1);
    %
    %     h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
    %         'GroupOrder',[{'A'};{'B'};{'C'};{'D'}]);
    %     h(1).ViolinPlot.FaceAlpha = 0.1;
    %     h(2).ViolinPlot.FaceAlpha = 0.1;
    %     h(3).ViolinPlot.FaceAlpha = 0.1;
    %     h(4).ViolinPlot.FaceAlpha = 0.1;
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
    %     %text(2.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %     %    'HorizontalAlignment','center');
    %     text(2.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %         'HorizontalAlignment','center');
    %
    %     set(gca,'linewidth',1.5)
    %     xlim([0.15 4.65])
    %     %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
    %     ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
    %     set(gca, 'FontSize', 10)
    %
    %     xtl = ["Over-mismatch"; "High-match";"Low-match"; "Under-mismatch"];
    %     xt=get(gca,'XTick');
    %     yt=get(gca,'YTick');
    %     xtext_xp=xt;
    %     %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;%0.56,0.4
    %     if if_monkey_D0_Z1 == 0
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.30;%0.56,0.4
    %     elseif if_monkey_D0_Z1 == 1
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.40;%0.56,0.4
    %     end
    %     text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
    %     set(gca,'xticklabel','');
    %
    %
    %     set(gca,'box','off');% 取消右、上边框
    %     %ylabel('Mean reward', 'FontSize', 12, 'FontWeight', 'bold');
    %     %temp_title = title(sprintf('Delay1'),'FontSize',12,'FontWeight','bold');
    %     temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
    %     temp_title.Interpreter = 'none';
    
    
    
    
    %     fig = figure('Name','asd','NumberTitle','off');
    %     set(gcf,'Position',[10 400 360*0.85*1.05 500*0.8*0.9*1.05*1.03]);
    %     %t = tiledlayout(2,5,'TileSpacing','tight','Padding','Compact');
    %     t = tiledlayout(2,2,'TileSpacing','tight','Padding','Loose');
    %
    %
    %     %t.Title.String = sprintf('Trial history, %s',FOVName_currentFOV2);
    %     %t.Title.FontSize = 12;
    %     %t.Title.Interpreter = 'none';
    %
    %
    %
    %     %% Plot baseline
    %     nexttile
    %     %nexttile([1 2])
    %
    %     temp_p = p_a5_a6(end);
    %     temp_1_raw = a5(:,end);
    %     temp_2_raw = a6(:,end);
    %
    %     temp_1 = temp_1_raw(~isnan(temp_1_raw));
    %     temp_2 = temp_2_raw(~isnan(temp_2_raw));
    %
    %     temp1_SEM = a5_sem(end);
    %     temp2_SEM = a6_sem(end);
    %
    %
    %     %temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM]);
    %     %temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
    %     %
    %     %temp_bar = bar([1 2], [mean(temp_1) mean(temp_2)], ...
    %     %    'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1);
    %     %hold on
    %     %a = 1;
    %     %temp_bar.CData(1,:) = [1 1 1]*0.5;%0.3
    %     %temp_bar.CData(2,:) = [1 1 1]*0.5;
    %     %
    %     %errorbar([1 2], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 8);
    %
    %
    %     temp_y_min = min([temp_1;temp_2]);
    %     temp_y_max = max([temp_1;temp_2]);
    %
    %     %temp_category = ones(size(temp_1,1)+size(temp_2,1),1);
    %     %temp_category(size(temp_1,1)+1:end) = 2;
    %     %violinplot([temp_1;temp_2],temp_category,'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
    %
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
    %     %text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %     %   'HorizontalAlignment','center');
    %     text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %         'HorizontalAlignment','center');
    %
    %     set(gca,'linewidth',1.5)
    %     xlim([0.15 2.65])
    %     %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
    %     ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
    %     set(gca, 'FontSize', 10)
    %
    %     xtl = ["High-meta"; "Low-meta"];
    %     xt=get(gca,'XTick');
    %     yt=get(gca,'YTick');
    %     xtext_xp=xt;
    %     %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.425;%0.56,0.4
    %     xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.315;%0.56,0.4
    %     text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',20,'fontsize',10);%25
    %     set(gca,'xticklabel','');
    %
    %
    %     set(gca,'box','off');% 取消右、上边框
    %     ylabel('Mean reward', 'FontSize', 12, 'FontWeight', 'bold');
    %     temp_title = title(sprintf('Baseline'),'FontSize',12,'FontWeight','bold');
    %     temp_title.Interpreter = 'none';
    %
    %
    %     %% Plot delay1
    %     nexttile
    %     %nexttile([1 2])
    %
    %     temp_p = p_a9_a10(end);
    %     temp_1_raw = a9(:,end);
    %     temp_2_raw = a10(:,end);
    %
    %     temp_1 = temp_1_raw(~isnan(temp_1_raw));
    %     temp_2 = temp_2_raw(~isnan(temp_2_raw));
    %
    %     temp1_SEM = a9_sem(end);
    %     temp2_SEM = a10_sem(end);
    %
    %
    %     %temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM]);
    %     %temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
    %     %
    %     %temp_bar = bar([1 2], [mean(temp_1) mean(temp_2)], ...
    %     %    'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1);
    %     %hold on
    %     %a = 1;
    %     %temp_bar.CData(1,:) = [1 1 1]*0.5;%0.3
    %     %temp_bar.CData(2,:) = [1 1 1]*0.5;
    %     %
    %     %errorbar([1 2], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 8);
    %
    %
    %     temp_y_min = min([temp_1;temp_2]);
    %     temp_y_max = max([temp_1;temp_2]);
    %
    %     %temp_category = ones(size(temp_1,1)+size(temp_2,1),1);
    %     %temp_category(size(temp_1,1)+1:end) = 2;
    %     %violinplot([temp_1;temp_2],temp_category,'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
    %
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
    %     %text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %     %   'HorizontalAlignment','center');
    %     text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %         'HorizontalAlignment','center');
    %
    %     set(gca,'linewidth',1.5)
    %     xlim([0.15 2.65])
    %     %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
    %     ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
    %     set(gca, 'FontSize', 10)
    %
    %     xtl = ["High-meta"; "Low-meta"];
    %     xt=get(gca,'XTick');
    %     yt=get(gca,'YTick');
    %     xtext_xp=xt;
    %     %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.425;%0.56,0.4
    %     xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.315;%0.56,0.4
    %     text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',20,'fontsize',10);%25
    %     set(gca,'xticklabel','');
    %
    %
    %     set(gca,'box','off');% 取消右、上边框
    %     ylabel('Mean reward', 'FontSize', 12, 'FontWeight', 'bold');
    %     temp_title = title(sprintf('Meta-memory'),'FontSize',12,'FontWeight','bold');
    %     temp_title.Interpreter = 'none';
    
    
    
    %     fig = figure('Name','asd','NumberTitle','off');
    %     %set(gcf,'Position',[810 400 240*1.15 240*0.7*0.85]);
    %     set(gcf,'Position',[810 400 240 143*0.97]);
    %     t = tiledlayout(1,1,'TileSpacing','tight','Padding','Compact');
    %
    %     %% Plot history weight fitting
    %     nexttile
    %
    %     x = normStdRange;
    %     y = temp_r2_fitting;
    %
    %     plot(x,y,'color',[0.25 0.25 0.25],'linewidth',1);
    %     hold on
    %
    %     [x_min,x_max] = bounds(x);
    %     [y_min,y_max] = bounds(y);
    %
    %     plot([1 1].*normStd_optimal,[y_min-(y_max-y_min)*0.1 temp_r2_optimal],':','color',[0.25 0.25 0.25],'linewidth',1);
    %     hold on
    %
    %     plot([x_min-(x_max-x_min)*0.05 normStd_optimal],[1 1].*temp_r2_optimal,':','color',[0.25 0.25 0.25],'linewidth',1);
    %     hold on
    %
    %
    %     %set(gca,'yscale','log')
    %
    %     set(gca,'linewidth',1.5)
    %     xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05]);%
    %     ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);%0.1
    %     set(gca, 'FontSize', 8)
    %     set(gca,'box','off');% 取消右、上边框
    %     xlabel('std of norm distribution', 'FontSize', 12, 'FontWeight', 'bold');
    %     ylabel('r2', 'FontSize', 12, 'FontWeight', 'bold');
    %     temp_title = title(sprintf('History weight fitting'),...
    %         'FontSize',12,'FontWeight','bold');
    %     temp_title.Interpreter = 'none';
    %     temp_title.Position(1) = temp_title.Position(1) + 5;
    %
    %
    %
    %     fig = figure('Name','asd','NumberTitle','off');
    %     %set(gcf,'Position',[810 400 240*1.15 240*0.7*0.85]);
    %     set(gcf,'Position',[810 600 240 143*1.02]);
    %     t = tiledlayout(1,1,'TileSpacing','tight','Padding','Compact');
    %
    %     %% Plot history weight fitting
    %     nexttile
    %
    %     x = normStdRange;
    %     y = AUROC_fitting;
    %
    %     plot(x,y,'color',[0.25 0.25 0.25],'linewidth',1);
    %     hold on
    %
    %     [x_min,x_max] = bounds(x);
    %     [y_min,y_max] = bounds(y);
    %     y_min = min(y_min,0.5);
    %
    %     plot([1 1].*normStd_optimal,[y_min-(y_max-y_min)*0.1 AUROC_optimal],':','color',[0.25 0.25 0.25],'linewidth',1);
    %     hold on
    %
    %     plot([x_min-(x_max-x_min)*0.05 normStd_optimal],[1 1].*AUROC_optimal,':','color',[0.25 0.25 0.25],'linewidth',1);
    %     hold on
    %
    %     plot([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05],[1 1].*0.5,':','color',[0.25 0.25 0.25],'linewidth',1);
    %     hold on
    %
    %     %set(gca,'yscale','log')
    %
    %     set(gca,'linewidth',1.5)
    %     xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05]);%
    %     ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);%0.1
    %     set(gca, 'FontSize', 8)
    %     set(gca,'box','off');% 取消右、上边框
    %     xlabel('std of norm distribution', 'FontSize', 12, 'FontWeight', 'bold');
    %     ylabel('AUROC', 'FontSize', 12, 'FontWeight', 'bold');
    %     temp_title = title(sprintf('History weight fitting'),...
    %         'FontSize',12,'FontWeight','bold');
    %     temp_title.Interpreter = 'none';
    %     %temp_title.Position(1) = temp_title.Position(1) + 5;
    
    
    
    
    
    %% Plot reward history
    fig = figure('Name','asd','NumberTitle','off');
    %set(gcf,'Position',[810 400 240*1.15 240*0.7*0.85]);
    set(gcf,'Position',[1110 600 240*0.8*0.85 143*1.02]);
    t = tiledlayout(1,1,'TileSpacing','tight','Padding','Compact');
    
    %% Plot history weight
    nexttile
    
    temp_trialID = 150;%100
    
    x = range_trialHistory;
    y = score_trialHistory(temp_trialID-1+range_trialHistory,1);
    
    plot(x,y,'color',[0.25 0.25 0.25],'linewidth',1);
    hold on
    
    
    [y_min,y_max] = bounds(y);
    
    if if_xAxia_normal0_log1 == 1
        set(gca,'xscale','log')
    end
    %xticks([1:8,10,range_trialHistory(end)]);
    xticks([1:range_trialHistory(end)]);     %#ok<*NBRAK>
    xtickangle(0);
    
    set(gca,'linewidth',1.5)
    xlim([range_trialHistory(1)-1 range_trialHistory(end)+1]);
    ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);%0.1
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Trial history', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Reward', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('Reward History'),...
        'FontSize',12,'FontWeight','bold');
    temp_title.Interpreter = 'none';
    
    
    
    
    %     %% Baseline & 4 conditions
    %     if true
    %         fig = figure('Name','asd','NumberTitle','off');
    %         set(gcf,'Position',[10 400 360*0.85*1.05*0.7 500*0.8*0.9*1.05*1.03]);
    %         t = tiledlayout(2,1,'TileSpacing','tight','Padding','Loose');
    %
    %         nexttile
    %
    %         temp_p = p_linear_baseline(end);
    %         temp_1_raw = z3;
    %         temp_2_raw = z1;
    %         temp_3_raw = z2;
    %         temp_4_raw = z4;
    %
    %         temp_1 = temp_1_raw(~isnan(temp_1_raw));
    %         temp_2 = temp_2_raw(~isnan(temp_2_raw));
    %         temp_3 = temp_3_raw(~isnan(temp_3_raw));
    %         temp_4 = temp_4_raw(~isnan(temp_4_raw));
    %
    %
    %         temp_y_min = min([temp_1;temp_2;temp_3;temp_4]);
    %         temp_y_max = max([temp_1;temp_2;temp_3;temp_4]);
    %
    %         temp_data = [temp_1;temp_2;temp_3;temp_4];
    %
    %         g1 = repmat({'A'},length(temp_1),1);
    %         g2 = repmat({'B'},length(temp_2),1);
    %         g3 = repmat({'C'},length(temp_3),1);
    %         g4 = repmat({'D'},length(temp_4),1);
    %
    %         temp_label = [g1;g2;g3;g4];
    %
    %         temptemp_color1 = [1 1 1]*0.5;
    %         temptemp_color2 = repmat(temptemp_color1, 4, 1);
    %
    %         h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
    %             'GroupOrder',[{'A'};{'B'};{'C'};{'D'}]);
    %         h(1).ViolinPlot.FaceAlpha = 0.1;
    %         h(2).ViolinPlot.FaceAlpha = 0.1;
    %         h(3).ViolinPlot.FaceAlpha = 0.1;
    %         h(4).ViolinPlot.FaceAlpha = 0.1;
    %
    %
    %         tempTxt = sprintf('');
    %         if temp_p < 0.001
    %             tempTxt = sprintf('***');
    %         elseif temp_p < 0.01
    %             tempTxt = sprintf('**');
    %         elseif temp_p < 0.05
    %             tempTxt = sprintf('*');
    %         end
    %         %text(2.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %         %    'HorizontalAlignment','center');
    %         text(2.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %             'HorizontalAlignment','center');
    %
    %         set(gca,'linewidth',1.5)
    %         xlim([0.15 4.65])
    %         %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
    %         ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
    %         set(gca, 'FontSize', 10)
    %
    %         xtl = ["Over-mismatch"; "High-match";"Low-match"; "Under-mismatch"];
    %         xt=get(gca,'XTick');
    %         yt=get(gca,'YTick');
    %         xtext_xp=xt;
    %         %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;%0.56,0.4
    %         if if_monkey_D0_Z1 == 0
    %             xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;%0.56,0.4
    %         elseif if_monkey_D0_Z1 == 1
    %             xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.40;%0.56,0.4
    %         end
    %         text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
    %         set(gca,'xticklabel','');
    %
    %
    %         set(gca,'box','off');% 取消右、上边框
    %         ylabel('Baseline meta', 'FontSize', 12, 'FontWeight', 'bold');
    %         temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
    %         temp_title.Interpreter = 'none';
    %
    %     end
    %
    %     %% Baseline & 3 conditions
    %     if true
    %         fig = figure('Name','asd','NumberTitle','off');
    %         set(gcf,'Position',[10 400 360*0.85*1.05*0.7 500*0.8*0.9*1.05*1.03]);
    %         t = tiledlayout(2,1,'TileSpacing','tight','Padding','Loose');
    %
    %         nexttile
    %
    %         temp_p = p_linear_baseline3(end);
    %         temp_1_raw = z3;
    %         temp_2_raw = z12;
    %         temp_3_raw = z4;
    %
    %         temp_1 = temp_1_raw(~isnan(temp_1_raw));
    %         temp_2 = temp_2_raw(~isnan(temp_2_raw));
    %         temp_3 = temp_3_raw(~isnan(temp_3_raw));
    %
    %
    %         temp_y_min = min([temp_1;temp_2;temp_3]);
    %         temp_y_max = max([temp_1;temp_2;temp_3]);
    %
    %         temp_data = [temp_1;temp_2;temp_3];
    %
    %         g1 = repmat({'A'},length(temp_1),1);
    %         g2 = repmat({'B'},length(temp_2),1);
    %         g3 = repmat({'C'},length(temp_3),1);
    %
    %         temp_label = [g1;g2;g3];
    %
    %         temptemp_color1 = [1 1 1]*0.5;
    %         temptemp_color2 = repmat(temptemp_color1, 3, 1);
    %
    %         h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
    %             'GroupOrder',[{'A'};{'B'};{'C'}]);
    %         h(1).ViolinPlot.FaceAlpha = 0.1;
    %         h(2).ViolinPlot.FaceAlpha = 0.1;
    %         h(3).ViolinPlot.FaceAlpha = 0.1;
    %
    %
    %         tempTxt = sprintf('');
    %         if temp_p < 0.001
    %             tempTxt = sprintf('***');
    %         elseif temp_p < 0.01
    %             tempTxt = sprintf('**');
    %         elseif temp_p < 0.05
    %             tempTxt = sprintf('*');
    %         end
    %         text(2,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %             'HorizontalAlignment','center');
    %
    %         set(gca,'linewidth',1.5)
    %         xlim([0.15 3.65])
    %         ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.13]);
    %         set(gca, 'FontSize', 10)
    %
    %         xtl = ["Over-mismatch"; "Match"; "Under-mismatch"];
    %         xt=get(gca,'XTick');
    %         yt=get(gca,'YTick');
    %         xtext_xp=xt;
    %         if if_monkey_D0_Z1 == 0
    %             xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;%0.56,0.4
    %         elseif if_monkey_D0_Z1 == 1
    %             xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.40;%0.56,0.4
    %         end
    %         text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
    %         set(gca,'xticklabel','');
    %
    %
    %         set(gca,'box','off');% 取消右、上边框
    %         ylabel('Baseline meta', 'FontSize', 12, 'FontWeight', 'bold');
    %         temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
    %         temp_title.Interpreter = 'none';
    %
    %     end
    
    
    %     %%Baseline & meta
    %     if true
    %         fig = figure('Name','asd','NumberTitle','off');
    %         set(gcf,'Position',[10 400 360*0.85*1.05*0.6 500*0.8*0.9*1.05*1.03]);
    %         %t = tiledlayout(2,5,'TileSpacing','tight','Padding','Compact');
    %         t = tiledlayout(2,1,'TileSpacing','tight','Padding','Loose');
    %
    %
    %         %% Plot baseline & meta
    %         nexttile
    %         %nexttile([1 2])
    %
    %         temp_p = p_z5_z6;
    %         temp_1_raw = z5;
    %         temp_2_raw = z6;
    %
    %         temp_1 = temp_1_raw(~isnan(temp_1_raw));
    %         temp_2 = temp_2_raw(~isnan(temp_2_raw));
    %
    %         temp1_SEM = z5_sem(end);
    %         temp2_SEM = z6_sem(end);
    %
    %
    %         temp_y_min = min([temp_1;temp_2]);
    %         temp_y_max = max([temp_1;temp_2]);
    %
    %         temp_data = [temp_1;temp_2];
    %
    %         g1 = repmat({'A'},length(temp_1),1);
    %         g2 = repmat({'B'},length(temp_2),1);
    %
    %         temp_label = [g1;g2];
    %
    %         temptemp_color1 = [1 1 1]*0.5;
    %         temptemp_color2 = repmat(temptemp_color1, 2, 1);
    %
    %         h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
    %             'GroupOrder',[{'A'};{'B'}]);
    %         h(1).ViolinPlot.FaceAlpha = 0.1;
    %         h(2).ViolinPlot.FaceAlpha = 0.1;
    %
    %
    %         tempTxt = sprintf('');
    %         if temp_p < 0.001
    %             tempTxt = sprintf('***');
    %         elseif temp_p < 0.01
    %             tempTxt = sprintf('**');
    %         elseif temp_p < 0.05
    %             tempTxt = sprintf('*');
    %         end
    %         text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %             'HorizontalAlignment','center');
    %
    %         set(gca,'linewidth',1.5)
    %         xlim([0.15 2.65])
    %         ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
    %         set(gca, 'FontSize', 10)
    %
    %         xtl = ["High-meta"; "Low-meta"];
    %         xt=get(gca,'XTick');
    %         yt=get(gca,'YTick');
    %         xtext_xp=xt;
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.15;%0.56,0.4
    %         text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',20,'fontsize',10);%25
    %         set(gca,'xticklabel','');
    %
    %
    %         set(gca,'box','off');% 取消右、上边框
    %         ylabel('Baseline meta', 'FontSize', 12, 'FontWeight', 'bold');
    %         temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
    %         temp_title.Interpreter = 'none';
    %
    %     end
    %
    %
    %     %% Baseline & 4 conditions (B)
    %     if true
    %         fig = figure('Name','asd','NumberTitle','off');
    %         set(gcf,'Position',[10 400 360*0.85*1.05*0.7 500*0.8*0.9*1.05*1.03]);
    %         t = tiledlayout(2,1,'TileSpacing','tight','Padding','Loose');
    %
    %         nexttile
    %
    %         %temp_p = p_linear_baseline(end);
    %
    %         temp_p12 = p_z1_z4;
    %         temp_p34 = p_z2_z3;
    %
    %         temp_1_raw = z1;
    %         temp_2_raw = z4;
    %         temp_3_raw = z2;
    %         temp_4_raw = z3;
    %
    %         temp_1 = temp_1_raw(~isnan(temp_1_raw));
    %         temp_2 = temp_2_raw(~isnan(temp_2_raw));
    %         temp_3 = temp_3_raw(~isnan(temp_3_raw));
    %         temp_4 = temp_4_raw(~isnan(temp_4_raw));
    %
    %
    %         temp_y_min = min([temp_1;temp_2;temp_3;temp_4]);
    %         temp_y_max = max([temp_1;temp_2;temp_3;temp_4]);
    %
    %         temp_data = [temp_1;temp_2;temp_3;temp_4];
    %
    %         g1 = repmat({'A'},length(temp_1),1);
    %         g2 = repmat({'B'},length(temp_2),1);
    %         g3 = repmat({'C'},length(temp_3),1);
    %         g4 = repmat({'D'},length(temp_4),1);
    %
    %         temp_label = [g1;g2;g3;g4];
    %
    %         temptemp_color1 = [1 1 1]*0.5;
    %         temptemp_color2 = repmat(temptemp_color1, 4, 1);
    %
    %         h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
    %             'GroupOrder',[{'A'};{'B'};{'C'};{'D'}]);
    %         h(1).ViolinPlot.FaceAlpha = 0.1;
    %         h(2).ViolinPlot.FaceAlpha = 0.1;
    %         h(3).ViolinPlot.FaceAlpha = 0.1;
    %         h(4).ViolinPlot.FaceAlpha = 0.1;
    %
    %
    %         tempTxt = sprintf('');
    %         if temp_p12 < 0.001
    %             tempTxt = sprintf('***');
    %         elseif temp_p12 < 0.01
    %             tempTxt = sprintf('**');
    %         elseif temp_p12 < 0.05
    %             tempTxt = sprintf('*');
    %         end
    %         text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %             'HorizontalAlignment','center');
    %
    %         tempTxt = sprintf('');
    %         if temp_p34 < 0.001
    %             tempTxt = sprintf('***');
    %         elseif temp_p34 < 0.01
    %             tempTxt = sprintf('**');
    %         elseif temp_p34 < 0.05
    %             tempTxt = sprintf('*');
    %         end
    %         text(3.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %             'HorizontalAlignment','center');
    %
    %         set(gca,'linewidth',1.5)
    %         xlim([0.15 4.65])
    %         %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
    %         ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
    %         set(gca, 'FontSize', 10)
    %
    %         %xtl = ["Over-mismatch"; "High-match";"Low-match"; "Under-mismatch"];
    %         xtl = ["High-match"; "Under-mismatch";"Low-match"; "Over-mismatch"];
    %         xt=get(gca,'XTick');
    %         yt=get(gca,'YTick');
    %         xtext_xp=xt;
    %         %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;%0.56,0.4
    %         if if_monkey_D0_Z1 == 0
    %             xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;%0.56,0.4
    %         elseif if_monkey_D0_Z1 == 1
    %             xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.40;%0.56,0.4
    %         end
    %         text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
    %         set(gca,'xticklabel','');
    %
    %
    %         set(gca,'box','off');% 取消右、上边框
    %         ylabel('Baseline meta', 'FontSize', 12, 'FontWeight', 'bold');
    %         temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
    %         temp_title.Interpreter = 'none';
    %
    %     end
    
    
    %% Baseline & meta, Baseline & mismatch
    
    if true
        
        
        
    end
    
    if true
        fig = figure('Name','asd','NumberTitle','off');
        %set(gcf,'Position',[10 100 360*0.85*1.05*0.4*1.2*1.8*1.5 500*0.8*0.9*1.05*1.03*0.8*0.9]);
        set(gcf,'Position',[10 100 360*0.85*1.05*0.4*1.2*1.8*1.5*0.94 500*0.8*0.9*1.05*1.03*0.8*0.9*1.14]);
        t = tiledlayout(2,3,'TileSpacing','Loose','Padding','Loose');
        
        
        %% Plot baseline
        nexttile
        
        temp_p = p_za5_za6;
        temp_1_raw = za6;
        temp_2_raw = za5;
        
        temp_1 = temp_1_raw(~isnan(temp_1_raw));
        temp_2 = temp_2_raw(~isnan(temp_2_raw));
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
        temp_y_min = 0;
        temp_y_max = 1;
        
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
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1.1 1.9],temp_y_max+(temp_y_max-temp_y_min)*0.04*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        set(gca,'linewidth',1.5)
        xlim([0.35 2.65])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
        set(gca, 'FontSize', 8)
        
        xtl = ["Offload"; "Low-strength mismatch"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.125;%0.56,0.4
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',20,'fontsize',7);%25
        set(gca,'xticklabel','');
        
        %yticks([0 1]);
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Baseline meta', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
        %temp_title.Interpreter = 'none';
        
        %% Plot memory
        nexttile
        
        temp_p = p_zc5_zc6;
        temp_1_raw = zc6;
        temp_2_raw = zc5;
        
        temp_1 = temp_1_raw(~isnan(temp_1_raw));
        temp_2 = temp_2_raw(~isnan(temp_2_raw));
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
        temp_y_min = 0;
        temp_y_max = 1;
        
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
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1.1 1.9],temp_y_max+(temp_y_max-temp_y_min)*0.04*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        set(gca,'linewidth',1.5)
        xlim([0.35 2.65])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
        set(gca, 'FontSize', 8)
        
        xtl = ["Offload"; "Low-strength mismatch"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.125;%0.56,0.4
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',20,'fontsize',7);%25
        set(gca,'xticklabel','');
        
        %yticks([0 1]);
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('WM strength', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
        %temp_title.Interpreter = 'none';
        
        
        %% Plot meta
        nexttile
        
        temp_p = p_zb5_zb6;
        temp_1_raw = zb6;
        temp_2_raw = zb5;
        
        temp_1 = temp_1_raw(~isnan(temp_1_raw));
        temp_2 = temp_2_raw(~isnan(temp_2_raw));
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
        temp_y_min = 0;
        temp_y_max = 1;
        
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
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1.1 1.9],temp_y_max+(temp_y_max-temp_y_min)*0.04*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        set(gca,'linewidth',1.5)
        xlim([0.35 2.65])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
        set(gca, 'FontSize', 8)
        
        xtl = ["Offload"; "Low-strength mismatch"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.125;%0.56,0.4
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',20,'fontsize',7);%25
        set(gca,'xticklabel','');
        
        %yticks([0 1]);
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Meta-WM', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
        %temp_title.Interpreter = 'none';
        
        
        a = 1;
    end
    
    %% Baseline & meta, Baseline & mismatch
    if true
        fig = figure('Name','asd','NumberTitle','off');
        set(gcf,'Position',[10 400 360*0.85*1.05 500*0.8*0.9*1.05*1.03]);
        %t = tiledlayout(2,5,'TileSpacing','tight','Padding','Compact');
        %t = tiledlayout(2,5,'TileSpacing','tight','Padding','Loose');
        t = tiledlayout(2,5,'TileSpacing','Loose','Padding','Loose');
        
        %% Plot baseline & meta
        nexttile([1 2])
        
        temp_p = p_z5_z6;
        temp_1_raw = z5;
        temp_2_raw = z6;
        
        temp_1 = temp_1_raw(~isnan(temp_1_raw));
        temp_2 = temp_2_raw(~isnan(temp_2_raw));
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
        
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
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.85])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
        set(gca, 'FontSize', 10)
        
        xtl = ["High-meta"; "Low-meta"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %         if if_monkey_D0_Z1 == 0
        %             xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.135;%0.56,0.4
        %         elseif if_monkey_D0_Z1 == 1
        %             xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.365;%0.56,0.4
        %         end
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.209;%0.56,0.4
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',20,'fontsize',9);%25
        set(gca,'xticklabel','');
        
        %yticks([0 1]);
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Baseline meta', 'FontSize', 10, 'FontWeight', 'normal');
        %temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
        %temp_title.Interpreter = 'none';
        
        a = 1;
        
        %% Plot baseline & mismatch
        nexttile([1 3])
        
        temp_p12 = p_z1_z4;
        temp_p34 = p_z2_z3;
        
        temp_1_raw = z1;
        temp_2_raw = z4;
        temp_3_raw = z2;
        temp_4_raw = z3;
        
        temp_1 = temp_1_raw(~isnan(temp_1_raw));
        temp_2 = temp_2_raw(~isnan(temp_2_raw));
        temp_3 = temp_3_raw(~isnan(temp_3_raw));
        temp_4 = temp_4_raw(~isnan(temp_4_raw));
        
        
        temp_y_min = min([temp_1;temp_2;temp_3;temp_4]);
        temp_y_max = max([temp_1;temp_2;temp_3;temp_4]);
        
        temp_data = [temp_1;temp_2;temp_3;temp_4];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'C'},length(temp_3),1);
        g4 = repmat({'D'},length(temp_4),1);
        
        temp_label = [g1;g2;g3;g4];
        
        %temptemp_color1 = [1 1 1]*0.5;
        %temptemp_color2 = repmat(temptemp_color1, 4, 1);
        
        temptemp_color2 = [...
            color_choiceMemoryHigh;
            color_choiceOffloadHigh;
            color_choiceOffloadLow;
            color_choiceMemoryLow
            ];
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'};{'D'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        h(3).ViolinPlot.FaceAlpha = 0.1;
        h(4).ViolinPlot.FaceAlpha = 0.1;
        
        
        tempTxt = sprintf('');
        if temp_p12 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p12 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p12 < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p34 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p34 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p34 < 0.05
            tempTxt = sprintf('*');
        end
        text(3.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 4.85])
        %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
        set(gca, 'FontSize', 10)
        
        %xtl = ["Over-mismatch"; "High-match";"Low-match"; "Under-mismatch"];
        xtl = ["High-match"; "Under-mismatch";"Low-match"; "Over-mismatch"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;%0.56,0.4
        if if_monkey_D0_Z1 == 0
            xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.20;%0.56,0.4
        elseif if_monkey_D0_Z1 == 1
            xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.40;%0.56,0.4
        end
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
        set(gca,'xticklabel','');
        
        yticks([0 1]);
        set(gca,'yticklabel','');
        
        set(gca,'box','off');% 取消右、上边框
        %ylabel('Baseline meta', 'FontSize', 12, 'FontWeight', 'bold');
        %temp_title = title(sprintf(''),'FontSize',10,'FontWeight','bold');
        %temp_title.Interpreter = 'none';
    end
    
    
    %% Baseline & precision
    if true
        fig = figure('Name','asd','NumberTitle','off');
        set(gcf,'Position',[10 400 360*0.85*1.05*0.4*1.1*1.1*0.975 500*0.8*0.9*1.05*1.03]);
        t = tiledlayout(2,1,'TileSpacing','Loose','Padding','Loose');
        
        %% Plot baseline & meta
        nexttile
        
        temp_p = p_z7_z8;
        temp_1_raw = z7;
        temp_2_raw = z8;
        
        temp_1 = temp_1_raw(~isnan(temp_1_raw));
        temp_2 = temp_2_raw(~isnan(temp_2_raw));
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
        
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
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.85])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
        set(gca, 'FontSize', 10)
        
        xtl = ["High-precision"; "Low-precision"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %         if if_monkey_D0_Z1 == 0
        %             xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.175;%0.56,0.4
        %         elseif if_monkey_D0_Z1 == 1
        %             xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.365;%0.56,0.4
        %         end
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.209;%0.56,0.4
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',20,'fontsize',9);%25
        set(gca,'xticklabel','');
        
        %yticks([0 1]);
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Baseline meta', 'FontSize', 10, 'FontWeight', 'normal');
        %temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
        %temp_title.Interpreter = 'none';
    end
    
    
    
    %% History & baseline, History & mismatch
    if true
        fig = figure('Name','asd','NumberTitle','off');
        %set(gcf,'Position',[410 400 360*0.85*1.05 500*0.8*0.9*1.05*1.03]);
        %set(gcf,'Position',[410 400 360*0.85*1.05 500*0.8*0.9*1.05*1.03*1.2]);
        set(gcf,'Position',[10 400 360*0.85*1.05 500*0.8*0.9*1.05*1.03]);
        %t = tiledlayout(2,5,'TileSpacing','tight','Padding','Loose');
        t = tiledlayout(2,5,'TileSpacing','Loose','Padding','Loose');
        
        
        %% Plot History & baseline
        nexttile([1 2])
        
        temp_p = p_a5_a6(end);
        temp_1_raw = a5(:,end);
        temp_2_raw = a6(:,end);
        
        temp_1 = temp_1_raw(~isnan(temp_1_raw));
        temp_2 = temp_2_raw(~isnan(temp_2_raw));
        
        temp1_SEM = a5_sem(end);
        temp2_SEM = a6_sem(end);
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
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
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.85])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
        set(gca, 'FontSize', 10)
        
        %xtl = ["High-meta"; "Low-meta"];
        xtl = ["High-baseline"; "Low-baseline"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %         if if_monkey_D0_Z1 == 0
        %             xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.315;%0.56,0.4
        %         elseif if_monkey_D0_Z1 == 1
        %             xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.365;%0.56,0.4
        %         end
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.15;%0.56,0.4
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',20,'fontsize',9);%25
        set(gca,'xticklabel','');
        
        yticks([1 2 3]);
        
        set(gca,'box','off');% 取消右、上边框
        %ylabel('History reward', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Weighted reward', 'FontSize', 10, 'FontWeight', 'normal');
        %temp_title = title(sprintf('Baseline'),'FontSize',12,'FontWeight','bold');
        %temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
        %temp_title.Interpreter = 'none';
        
        
        %% Plot History & mismatch
        nexttile([1 3])
        
        temp_p12 = p_a1_a4(end);
        temp_p34 = p_a2_a3(end);
        temp_1_raw = a1(:,end);
        temp_2_raw = a4(:,end);
        temp_3_raw = a2(:,end);
        temp_4_raw = a3(:,end);
        
        temp_1 = temp_1_raw(~isnan(temp_1_raw));
        temp_2 = temp_2_raw(~isnan(temp_2_raw));
        temp_3 = temp_3_raw(~isnan(temp_3_raw));
        temp_4 = temp_4_raw(~isnan(temp_4_raw));
        
        
        temp_y_min = min([temp_1;temp_2;temp_3;temp_4]);
        temp_y_max = max([temp_1;temp_2;temp_3;temp_4]);
        
        temp_data = [temp_1;temp_2;temp_3;temp_4];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'C'},length(temp_3),1);
        g4 = repmat({'D'},length(temp_4),1);
        
        temp_label = [g1;g2;g3;g4];
        
        %temptemp_color1 = [1 1 1]*0.5;
        %temptemp_color2 = repmat(temptemp_color1, 4, 1);
        
        temptemp_color2 = [...
            color_choiceMemoryHigh;
            color_choiceOffloadHigh;
            color_choiceOffloadLow;
            color_choiceMemoryLow
            ];
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'};{'D'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        h(3).ViolinPlot.FaceAlpha = 0.1;
        h(4).ViolinPlot.FaceAlpha = 0.1;
        
        
        tempTxt = sprintf('');
        if temp_p12 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p12 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p12 < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p34 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p34 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p34 < 0.05
            tempTxt = sprintf('*');
        end
        text(3.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 4.85])
        %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
        set(gca, 'FontSize', 10)
        
        %xtl = ["Over-mismatch"; "High-match";"Low-match"; "Under-mismatch"];
        xtl = ["High-match"; "Under-mismatch";"Low-match"; "Over-mismatch"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;%0.56,0.4
        if if_monkey_D0_Z1 == 0
            xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.33;%0.56,0.4
        elseif if_monkey_D0_Z1 == 1
            xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.40;%0.56,0.4
        end
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
        set(gca,'xticklabel','');
        
        yticks([1 2 3]);
        set(gca,'yticklabel','');
        
        set(gca,'box','off');% 取消右、上边框
        %ylabel('Mean reward', 'FontSize', 12, 'FontWeight', 'bold');
        %temp_title = title(sprintf('Delay1'),'FontSize',12,'FontWeight','bold');
        %temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
        %temp_title.Interpreter = 'none';
        
        
    end
    
    
    
    %% History & baseline, History & meta, History & mismatch
    if true
        fig = figure('Name','asd','NumberTitle','off');
        %set(gcf,'Position',[410 400 360*0.85*1.05 500*0.8*0.9*1.05*1.03]);
        set(gcf,'Position',[410 400 360*0.85*1.05*1.4 500*0.8*0.9*1.05*1.03*1.2]);
        %t = tiledlayout(2,5,'TileSpacing','tight','Padding','Loose');
        t = tiledlayout(2,7,'TileSpacing','Loose','Padding','Loose');
        
        
        %% Plot History & baseline
        nexttile([1 2])
        
        temp_p = p_a5_a6(end);
        temp_1_raw = a5(:,end);
        temp_2_raw = a6(:,end);
        
        temp_1 = temp_1_raw(~isnan(temp_1_raw));
        temp_2 = temp_2_raw(~isnan(temp_2_raw));
        
        temp1_SEM = a5_sem(end);
        temp2_SEM = a6_sem(end);
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
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
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.85])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
        set(gca, 'FontSize', 10)
        
        %xtl = ["High-meta"; "Low-meta"];
        xtl = ["High-baseline"; "Low-baseline"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        if if_monkey_D0_Z1 == 0
            xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.315;%0.56,0.4
        elseif if_monkey_D0_Z1 == 1
            xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.365;%0.56,0.4
        end
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',20,'fontsize',9);%25
        set(gca,'xticklabel','');
        
        yticks([1 2 3]);
        
        set(gca,'box','off');% 取消右、上边框
        %ylabel('History reward', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Weighted reward', 'FontSize', 10, 'FontWeight', 'normal');
        %temp_title = title(sprintf('Baseline'),'FontSize',12,'FontWeight','bold');
        %temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
        %temp_title.Interpreter = 'none';
        
        
        %% Plot History & meta
        nexttile([1 2])
        
        temp_p = p_a9_a10(end);
        temp_1_raw = a9(:,end);
        temp_2_raw = a10(:,end);
        
        temp_1 = temp_1_raw(~isnan(temp_1_raw));
        temp_2 = temp_2_raw(~isnan(temp_2_raw));
        
        temp1_SEM = a9_sem(end);
        temp2_SEM = a10_sem(end);
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
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
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.85])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
        set(gca, 'FontSize', 10)
        
        xtl = ["High-meta"; "Low-meta"];
        %xtl = ["High-baseline"; "Low-baseline"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        if if_monkey_D0_Z1 == 0
            xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.315;%0.56,0.4
        elseif if_monkey_D0_Z1 == 1
            xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.365;%0.56,0.4
        end
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',20,'fontsize',9);%25
        set(gca,'xticklabel','');
        
        yticks([1 2 3]);
        set(gca,'yticklabel','');
        
        set(gca,'box','off');% 取消右、上边框
        %ylabel('History reward', 'FontSize', 12, 'FontWeight', 'bold');
        %ylabel('Weighted reward', 'FontSize', 10, 'FontWeight', 'normal');
        %temp_title = title(sprintf('Baseline'),'FontSize',12,'FontWeight','bold');
        %temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
        %temp_title.Interpreter = 'none';
        
        
        %% Plot History & mismatch
        nexttile([1 3])
        
        temp_p12 = p_a1_a4(end);
        temp_p34 = p_a2_a3(end);
        temp_1_raw = a1(:,end);
        temp_2_raw = a4(:,end);
        temp_3_raw = a2(:,end);
        temp_4_raw = a3(:,end);
        
        temp_1 = temp_1_raw(~isnan(temp_1_raw));
        temp_2 = temp_2_raw(~isnan(temp_2_raw));
        temp_3 = temp_3_raw(~isnan(temp_3_raw));
        temp_4 = temp_4_raw(~isnan(temp_4_raw));
        
        
        temp_y_min = min([temp_1;temp_2;temp_3;temp_4]);
        temp_y_max = max([temp_1;temp_2;temp_3;temp_4]);
        
        temp_data = [temp_1;temp_2;temp_3;temp_4];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'C'},length(temp_3),1);
        g4 = repmat({'D'},length(temp_4),1);
        
        temp_label = [g1;g2;g3;g4];
        
        %temptemp_color1 = [1 1 1]*0.5;
        %temptemp_color2 = repmat(temptemp_color1, 4, 1);
        
        temptemp_color2 = [...
            color_choiceMemoryHigh;
            color_choiceOffloadHigh;
            color_choiceOffloadLow;
            color_choiceMemoryLow
            ];
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'};{'D'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        h(3).ViolinPlot.FaceAlpha = 0.1;
        h(4).ViolinPlot.FaceAlpha = 0.1;
        
        
        tempTxt = sprintf('');
        if temp_p12 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p12 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p12 < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p34 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p34 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p34 < 0.05
            tempTxt = sprintf('*');
        end
        text(3.5,temp_y_max+(temp_y_max-temp_y_min)*0.04,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 4.85])
        %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.08 temp_y_max+(temp_y_max-temp_y_min)*0.08]);
        set(gca, 'FontSize', 10)
        
        %xtl = ["Over-mismatch"; "High-match";"Low-match"; "Under-mismatch"];
        xtl = ["High-match"; "Under-mismatch";"Low-match"; "Over-mismatch"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;%0.56,0.4
        if if_monkey_D0_Z1 == 0
            xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.33;%0.56,0.4
        elseif if_monkey_D0_Z1 == 1
            xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.40;%0.56,0.4
        end
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
        set(gca,'xticklabel','');
        
        yticks([1 2 3]);
        set(gca,'yticklabel','');
        
        set(gca,'box','off');% 取消右、上边框
        %ylabel('Mean reward', 'FontSize', 12, 'FontWeight', 'bold');
        %temp_title = title(sprintf('Delay1'),'FontSize',12,'FontWeight','bold');
        %temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
        %temp_title.Interpreter = 'none';
        
        
    end
    
    
end

%% history weight
if true
    fig = figure('Name','asd','NumberTitle','off');
    %set(gcf,'Position',[810 400 240*1.15 240*0.7*0.85]);
    %set(gcf,'Position',[1110 400 240*0.8*0.85 143*1.02]);
    %set(gcf,'Position',[1110 400 240*0.8*0.85 143*1.02*0.88]);
    %set(gcf,'Position',[1110 400 240*0.8*0.85*1.85 143*1.02*0.88]);
    %set(gcf,'Position',[1110 400 240*0.8*0.85*1.85*0.7 143*1.02*0.88]);
    set(gcf,'Position',[1110 400 240*0.8*0.85*1.85*0.7*1.5*0.95 143*1.02*0.88*1.2*1.2]);
    t = tiledlayout(1,1,'TileSpacing','tight','Padding','Compact');
    
    %% Plot history weight
    nexttile
    
    x = range_trialHistory;
    y = weight_trialHistory;
    
    plot(x,y,'color',[0.25 0.25 0.25],'linewidth',1);
    hold on
    
    
    [y_min,y_max] = bounds(y);
    
    if if_xAxia_normal0_log1 == 1
        set(gca,'xscale','log')
    end
    if if_fitNormal0_fitEnd1_fitBeta2_fitExp3 ~= 1
        xticks([1:8,10,range_trialHistory(end)]);
    elseif if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 1
        xticks([1:range_trialHistory(end)]);
    end
    xtickangle(0);
    
    %yticks([0 0.1]);
    
    set(gca,'linewidth',1.5)
    xlim([range_trialHistory(1)-1 range_trialHistory(end)+1]);
    %ylim([y_min y_max+(y_max-y_min)*0.1]);%0.1
    %ylim([0 y_max+(y_max-y_min)*0.1]);%0.1
    ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);%0.1
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Trial history', 'FontSize', 10, 'FontWeight', 'normal');
    ylabel('Weight', 'FontSize', 10, 'FontWeight', 'normal');
    temp_title = title(sprintf('History weight'),...
        'FontSize',10,'FontWeight','normal');
    temp_title.Interpreter = 'none';
    
    if if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 1
        subtitle(sprintf('AUROC=%.3f, %d trials, std=%.1f',AUROC_optimal,end_optimal,normStd_optimal),...
            'FontSize',8,'FontWeight','normal');
        %elseif if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 2 || if_fitNormal0_fitEnd1_fitBeta2_fitExp3 == 3
    else
        subtitle(sprintf('AUROC=%.3f',AUROC_optimal),...
            'FontSize',8,'FontWeight','normal');
    end
    
end

AUROC_optimal_trialHistory = AUROC_optimal;



score_trialHistory;
historyWeightedReward = score_trialHistory(:,end);

%% End