% Chuan's 26th script (20251214)
% This script:
% (1) Meta-memory linear regression using baseline meta-memory and (delay) memory strength as predictors.
% (2) Baseline meta linear regression using pupil size, reward history, and their interaction as predictors.
% (3) Comparison of Akaike information criterion (AIC) across three meta-WM regression models.
% (1) (2) (3) related to figure 5.
% (4) Relationship between three subspaces, related to figure 6.
if_usePupilSize = 0;%0

if_onlyGLM = 0;%0;

if_compute = 1;%1

if_resample = 1;
resampleIterCount_GLM = 1000;%1000

if_BIC0_AIC1 = 1;%0,1

if_allTrial0_freeChoice1 = 1;

if exist('if_compute_summary','var') == 1
    if_allTrial0_freeChoice1 = 0;% why?
end

if_plot = 1;%1

if_plot_violinplot0_pairline1 = 1;

%% Seq level
memoryPrecision_trialLevel;
meta_trialLevel;

if if_memeoryPrecision_stimuli0_response1 == 0
    temp_seqIndex = seqIndex;
elseif if_memeoryPrecision_stimuli0_response1 == 1
    temp_seqIndex = seqIndex_response;
end

meta_seqLevel_all = zeros(sum(numSeq(1:3)),1);
memoryPrecision_seqLevel_all = zeros(sum(numSeq(1:3)),1);
for tempi=1:sum(numSeq(1:3))
    temptempBoolIndex = temp_seqIndex==tempi;
    meta_seqLevel_all(tempi) = mean(meta_trialLevel(temptempBoolIndex));
    memoryPrecision_seqLevel_all(tempi) = mean(memoryPrecision_trialLevel(temptempBoolIndex));
end

x = memoryPrecision_seqLevel_all;
y = meta_seqLevel_all;
temp_mdl_seqLevel = fitglm(x,y,'linear');
r2_seqLevel = temp_mdl_seqLevel.Rsquared.Adjusted;
beta0_seqLevel = temp_mdl_seqLevel.Coefficients.Estimate(1);
beta1_seqLevel = temp_mdl_seqLevel.Coefficients.Estimate(2);



%%
tempNANBoolIndex1 = isnan(meta_trialLevel_baseline);
tempNANBoolIndex2 = isnan(memoryPrecision_trialLevel);
tempNANBoolIndex3 = isnan(meta_trialLevel);

tempNANBoolIndex123 = tempNANBoolIndex1 | tempNANBoolIndex2 | tempNANBoolIndex3;

tempNONNANBoolIndex123 = ~tempNANBoolIndex123;

if if_allTrial0_freeChoice1 == 1
    tempNONNANBoolIndex123 = tempNONNANBoolIndex123 & choiceBoolIndex';
    %tempNONNANBoolIndex123 = tempNONNANBoolIndex123 & choiceMemoryBoolIndex';
    %tempNONNANBoolIndex123 = tempNONNANBoolIndex123 & choiceOffloadBoolIndex';
    %tempNONNANBoolIndex123 = tempNONNANBoolIndex123 & (~choiceBoolIndex)';
end

tempNONNAN_meta_trialLevel_baseline = meta_trialLevel_baseline(tempNONNANBoolIndex123);
tempNONNAN_memoryPrecision_trialLevel = memoryPrecision_trialLevel(tempNONNANBoolIndex123);
tempNONNAN_meta_trialLevel = meta_trialLevel(tempNONNANBoolIndex123);


temptempBoolIndex = trialBoolIndex_memoryPrecisionLow_metaHigh_choiceMemoryError | choiceOffloadBoolIndex';
tempNONNAN_meta_CME_CF_trialLevel_baseline = meta_trialLevel_baseline(tempNONNANBoolIndex123&temptempBoolIndex);
tempNONNAN_memoryPrecision_CME_CF_trialLevel = memoryPrecision_trialLevel(tempNONNANBoolIndex123&temptempBoolIndex);
tempNONNAN_meta_CME_CF_trialLevel = meta_trialLevel(tempNONNANBoolIndex123&temptempBoolIndex);
tempNONNAN_boolIndex_CME = ~choiceOffloadBoolIndex(tempNONNANBoolIndex123&temptempBoolIndex)';


temptempBoolIndex = choiceMemoryCorrectBoolIndex';
tempNONNAN_meta_CMC_trialLevel_baseline = meta_trialLevel_baseline(tempNONNANBoolIndex123&temptempBoolIndex);
tempNONNAN_memoryPrecision_CMC_trialLevel = memoryPrecision_trialLevel(tempNONNANBoolIndex123&temptempBoolIndex);
tempNONNAN_meta_CMC_trialLevel = meta_trialLevel(tempNONNANBoolIndex123&temptempBoolIndex);

temptempBoolIndex = choiceMemoryBoolIndex';
tempNONNAN_meta_CM_trialLevel_baseline = meta_trialLevel_baseline(tempNONNANBoolIndex123&temptempBoolIndex);
tempNONNAN_memoryPrecision_CM_trialLevel = memoryPrecision_trialLevel(tempNONNANBoolIndex123&temptempBoolIndex);
tempNONNAN_meta_CM_trialLevel = meta_trialLevel(tempNONNANBoolIndex123&temptempBoolIndex);

temptempBoolIndex = choiceOffloadBoolIndex';
tempNONNAN_meta_CF_trialLevel_baseline = meta_trialLevel_baseline(tempNONNANBoolIndex123&temptempBoolIndex);
tempNONNAN_memoryPrecision_CF_trialLevel = memoryPrecision_trialLevel(tempNONNANBoolIndex123&temptempBoolIndex);
tempNONNAN_meta_CF_trialLevel = meta_trialLevel(tempNONNANBoolIndex123&temptempBoolIndex);

temptempBoolIndex = trialBoolIndex_memoryPrecisionLow_metaHigh_choiceMemoryError;
tempNONNAN_meta_CME_trialLevel_baseline = meta_trialLevel_baseline(tempNONNANBoolIndex123&temptempBoolIndex);
tempNONNAN_memoryPrecision_CME_trialLevel = memoryPrecision_trialLevel(tempNONNANBoolIndex123&temptempBoolIndex);
tempNONNAN_meta_CME_trialLevel = meta_trialLevel(tempNONNANBoolIndex123&temptempBoolIndex);

% temptempBoolIndex = meta_trialLevel_baseline > lowThreshold_meta;
temptempBoolIndex = meta_trialLevel_baseline <= lowThreshold_meta;
temp1 = prctile(meta_trialLevel_baseline,25);%25,15
temp2 = prctile(meta_trialLevel_baseline,75);%75,85
% temptempBoolIndex = (meta_trialLevel_baseline < temp1) | (meta_trialLevel_baseline > temp2);

temp1 = prctile(meta_trialLevel_baseline,25);%25,35
temp2 = prctile(meta_trialLevel_baseline,75);%75,65
% temptempBoolIndex = (meta_trialLevel_baseline > temp1) & (meta_trialLevel_baseline < temp2);

temp1 = prctile(memoryPrecision_trialLevel,25);%25,15
temp2 = prctile(memoryPrecision_trialLevel,75);%75,85
% temptempBoolIndex = (memoryPrecision_trialLevel < temp1) | (memoryPrecision_trialLevel > temp2);
% temptempBoolIndex = (memoryPrecision_trialLevel > temp1) & (memoryPrecision_trialLevel < temp2);

temp1 = prctile(memoryPrecision_trialLevel,33);%25,15
temp2 = prctile(memoryPrecision_trialLevel,66);%75,85
% temptempBoolIndex = memoryPrecision_trialLevel < temp1;
% temptempBoolIndex = (memoryPrecision_trialLevel > temp1) & (memoryPrecision_trialLevel < temp2);
% temptempBoolIndex = memoryPrecision_trialLevel > temp2;

tempNONNAN_meta_HBSL_trialLevel_baseline = meta_trialLevel_baseline(tempNONNANBoolIndex123&temptempBoolIndex);
tempNONNAN_memoryPrecision_HBSL_trialLevel = memoryPrecision_trialLevel(tempNONNANBoolIndex123&temptempBoolIndex);
tempNONNAN_meta_HBSL_trialLevel = meta_trialLevel(tempNONNANBoolIndex123&temptempBoolIndex);

if if_compute == 1
    if if_resample == 0
        
        temp_x1 = tempNONNAN_memoryPrecision_trialLevel;
        temp_x2 = tempNONNAN_meta_trialLevel_baseline;
        temp_y = tempNONNAN_meta_trialLevel;
        
        temp_x1_zscored = (temp_x1-mean(temp_x1))./std(temp_x1);
        temp_x2_zscored = (temp_x2-mean(temp_x2))./std(temp_x2);
        temp_y_zscored = (temp_y-mean(temp_y))./std(temp_y);
        
        % Case A
        x = temp_x1_zscored;
        y = temp_y_zscored;
        temp_mdl_caseA = fitglm(x,y,'linear');
        r2_caseA = temp_mdl_caseA.Rsquared.Adjusted;
        beta0_caseA = temp_mdl_caseA.Coefficients.Estimate(1);
        beta1_caseA = temp_mdl_caseA.Coefficients.Estimate(2);
        
        r2_caseA_resampled = r2_caseA+0.01*rand(1,100);
        
        % Case B
        x = temp_x2_zscored;
        y = temp_y_zscored;
        temp_mdl_caseB = fitglm(x,y,'linear');
        r2_caseB = temp_mdl_caseB.Rsquared.Adjusted;
        beta0_caseB = temp_mdl_caseB.Coefficients.Estimate(1);
        beta1_caseB = temp_mdl_caseB.Coefficients.Estimate(2);
        
        r2_caseB_resampled = r2_caseB+0.01*rand(1,100);
        
        % Case C
        x = [temp_x1_zscored,temp_x2_zscored];
        y = temp_y_zscored;
        %temp_mdl_caseC = fitglm(x,y,'linear');
        temp_mdl_caseC = fitglm(x,y,'interactions');
        r2_caseC = temp_mdl_caseC.Rsquared.Adjusted;
        beta0_caseC = temp_mdl_caseC.Coefficients.Estimate(1);
        beta1_caseC= temp_mdl_caseC.Coefficients.Estimate(2);
        beta2_caseC= temp_mdl_caseC.Coefficients.Estimate(3);
        
        r2_caseC_resampled = r2_caseC+0.01*rand(1,100);
        
        if false
            % Case C
            x = tempNONNAN_memoryPrecision_trialLevel.*tempNONNAN_meta_trialLevel_baseline; %#ok<*UNRCH>
            y = tempNONNAN_meta_trialLevel;
            temp_mdl_caseC = fitglm(x,y,'linear');
            r2_caseC = temp_mdl_caseC.Rsquared.Adjusted;
            
        end
        
        
    elseif if_resample == 1
        tempNONNAN_meta_trialLevel;
        resampleIterCount_GLM;
        
        temp_trialNum = length(tempNONNAN_meta_trialLevel);
        
        r2_caseA_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseA_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseA_resampled = nan(1,resampleIterCount_GLM);
        BIC_caseA_resampled = nan(1,resampleIterCount_GLM);
        
        r2_caseB_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseB_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseB_resampled = nan(1,resampleIterCount_GLM);
        BIC_caseB_resampled = nan(1,resampleIterCount_GLM);
        
        r2_caseC_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseC_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseC_resampled = nan(1,resampleIterCount_GLM);
        beta2_caseC_resampled = nan(1,resampleIterCount_GLM);
        beta3_caseC_resampled = nan(1,resampleIterCount_GLM);
        BIC_caseC_resampled = nan(1,resampleIterCount_GLM);
        
        
        temp_trialNum_CME_CF = length(tempNONNAN_meta_CME_CF_trialLevel_baseline);
        r2_caseD_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseD_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseD_resampled = nan(1,resampleIterCount_GLM);
        beta2_caseD_resampled = nan(1,resampleIterCount_GLM);
        beta3_caseD_resampled = nan(1,resampleIterCount_GLM);
        
        r2_caseD2_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseD2_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseD2_resampled = nan(1,resampleIterCount_GLM);
        beta2_caseD2_resampled = nan(1,resampleIterCount_GLM);
        beta3_caseD2_resampled = nan(1,resampleIterCount_GLM);
        
        r2_caseD2A_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseD2A_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseD2A_resampled = nan(1,resampleIterCount_GLM);
        
        r2_caseD2B_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseD2B_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseD2B_resampled = nan(1,resampleIterCount_GLM);
        
        
        temp_trialNum_CMC = length(tempNONNAN_meta_CMC_trialLevel_baseline);
        r2_caseE_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseE_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseE_resampled = nan(1,resampleIterCount_GLM);
        beta2_caseE_resampled = nan(1,resampleIterCount_GLM);
        beta3_caseE_resampled = nan(1,resampleIterCount_GLM);
        
        
        temp_trialNum_CM = length(tempNONNAN_meta_CM_trialLevel_baseline);
        r2_caseE2_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseE2_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseE2_resampled = nan(1,resampleIterCount_GLM);
        beta2_caseE2_resampled = nan(1,resampleIterCount_GLM);
        beta3_caseE2_resampled = nan(1,resampleIterCount_GLM);
        
        
        temp_trialNum_CF = length(tempNONNAN_meta_CF_trialLevel_baseline);
        r2_caseF_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseF_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseF_resampled = nan(1,resampleIterCount_GLM);
        beta2_caseF_resampled = nan(1,resampleIterCount_GLM);
        beta3_caseF_resampled = nan(1,resampleIterCount_GLM);
        
        
        temp_trialNum_CME = length(tempNONNAN_meta_CME_trialLevel_baseline);
        r2_caseG_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseG_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseG_resampled = nan(1,resampleIterCount_GLM);
        beta2_caseG_resampled = nan(1,resampleIterCount_GLM);
        beta3_caseG_resampled = nan(1,resampleIterCount_GLM);
        
        
        temp_trialNum_HBSL = length(tempNONNAN_meta_HBSL_trialLevel_baseline);
        r2_caseH_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseH_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseH_resampled = nan(1,resampleIterCount_GLM);
        beta2_caseH_resampled = nan(1,resampleIterCount_GLM);
        beta3_caseH_resampled = nan(1,resampleIterCount_GLM);
        
        
        %for iter=1:resampleIterCount_GLM
        parfor iter=1:resampleIterCount_GLM
            temp_trialIndex_resample = sort(randi(temp_trialNum,1,temp_trialNum));
            
            tempNONNAN_meta_trialLevel_baseline_resampled = tempNONNAN_meta_trialLevel_baseline(temp_trialIndex_resample); %#ok<*PFBNS>
            tempNONNAN_memoryPrecision_trialLevel_resampled = tempNONNAN_memoryPrecision_trialLevel(temp_trialIndex_resample);
            tempNONNAN_meta_trialLevel_resampled = tempNONNAN_meta_trialLevel(temp_trialIndex_resample);
            
            temp_x1 = tempNONNAN_memoryPrecision_trialLevel_resampled;
            temp_x2 = tempNONNAN_meta_trialLevel_baseline_resampled;
            temp_y = tempNONNAN_meta_trialLevel_resampled;
            
            temp_x1_zscored = (temp_x1-mean(temp_x1))./std(temp_x1);
            temp_x2_zscored = (temp_x2-mean(temp_x2))./std(temp_x2);
            temp_y_zscored = (temp_y-mean(temp_y))./std(temp_y);
            
            
            % Case A
            %x = temp_x1;
            %y = temp_y;
            
            x = temp_x1_zscored;
            y = temp_y_zscored;
            
            temp_mdl_caseA = fitglm(x,y,'linear');
            r2_caseA_resampled(iter) = temp_mdl_caseA.Rsquared.Adjusted;
            beta0_caseA_resampled(iter) = temp_mdl_caseA.Coefficients.Estimate(1);
            beta1_caseA_resampled(iter) = temp_mdl_caseA.Coefficients.Estimate(2);
            if if_BIC0_AIC1 == 0
                BIC_caseA_resampled(iter) = temp_mdl_caseA.ModelCriterion.BIC;
            elseif if_BIC0_AIC1 == 1
                BIC_caseA_resampled(iter) = temp_mdl_caseA.ModelCriterion.AIC;
            end
            
            % Case B
            %x = temp_x2;
            %y = temp_y;
            
            x = temp_x2_zscored;
            y = temp_y_zscored;
            
            temp_mdl_caseB = fitglm(x,y,'linear');
            r2_caseB_resampled(iter) = temp_mdl_caseB.Rsquared.Adjusted;
            beta0_caseB_resampled(iter) = temp_mdl_caseB.Coefficients.Estimate(1);
            beta1_caseB_resampled(iter) = temp_mdl_caseB.Coefficients.Estimate(2);
            if if_BIC0_AIC1 == 0
                BIC_caseB_resampled(iter) = temp_mdl_caseB.ModelCriterion.BIC;
            elseif if_BIC0_AIC1 == 1
                BIC_caseB_resampled(iter) = temp_mdl_caseB.ModelCriterion.AIC;
            end
            
            % Case C
            %x = [temp_x1,temp_x2];
            %y = temp_y;
            
            x = [temp_x1_zscored,temp_x2_zscored];
            y = temp_y_zscored;
            
            %temp_mdl_caseC = fitglm(x,y,'linear');
            temp_mdl_caseC = fitglm(x,y,'interactions');
            r2_caseC_resampled(iter) = temp_mdl_caseC.Rsquared.Adjusted;
            beta0_caseC_resampled(iter) = temp_mdl_caseC.Coefficients.Estimate(1);
            beta1_caseC_resampled(iter) = temp_mdl_caseC.Coefficients.Estimate(2);
            beta2_caseC_resampled(iter) = temp_mdl_caseC.Coefficients.Estimate(3);
            beta3_caseC_resampled(iter) = temp_mdl_caseC.Coefficients.Estimate(4);
            if if_BIC0_AIC1 == 0
                BIC_caseC_resampled(iter) = temp_mdl_caseC.ModelCriterion.BIC;
            elseif if_BIC0_AIC1 == 1
                BIC_caseC_resampled(iter) = temp_mdl_caseC.ModelCriterion.AIC;
            end
            
            
            % CME_CF
            temp_trialIndex_CME_CF_resample = sort(randi(temp_trialNum_CME_CF,1,temp_trialNum_CME_CF));
            
            tempNONNAN_meta_CME_CF_trialLevel_baseline_resampled = tempNONNAN_meta_CME_CF_trialLevel_baseline(temp_trialIndex_CME_CF_resample); %#ok<*PFBNS>
            tempNONNAN_memoryPrecision_CME_CF_trialLevel_resampled = tempNONNAN_memoryPrecision_CME_CF_trialLevel(temp_trialIndex_CME_CF_resample);
            tempNONNAN_meta_CME_CF_trialLevel_resampled = tempNONNAN_meta_CME_CF_trialLevel(temp_trialIndex_CME_CF_resample);
            
            tempNONNAN_boolIndex_CME_resampled = tempNONNAN_boolIndex_CME(temp_trialIndex_CME_CF_resample);
            
            temp_x1_CME_CF = tempNONNAN_memoryPrecision_CME_CF_trialLevel_resampled;
            temp_x2_CME_CF = tempNONNAN_meta_CME_CF_trialLevel_baseline_resampled;
            temp_y_CME_CF = tempNONNAN_meta_CME_CF_trialLevel_resampled;
            
            temp_x1_CME_CF_zscored = (temp_x1_CME_CF-mean(temp_x1_CME_CF))./std(temp_x1_CME_CF);
            temp_x2_CME_CF_zscored = (temp_x2_CME_CF-mean(temp_x2_CME_CF))./std(temp_x2_CME_CF);
            temp_y_CME_CF_zscored = (temp_y_CME_CF-mean(temp_y_CME_CF))./std(temp_y_CME_CF);
            
            temp_y2_CME_CF = tempNONNAN_boolIndex_CME_resampled;
            temp_y2_CME_CF_zscored = (temp_y2_CME_CF-mean(temp_y2_CME_CF))./std(temp_y2_CME_CF);
            
            
            % Case D
            x = [temp_x1_CME_CF_zscored,temp_x2_CME_CF_zscored];
            y = temp_y_CME_CF_zscored;
            
            %temp_mdl_caseD = fitglm(x,y,'linear');
            temp_mdl_caseD = fitglm(x,y,'interactions');
            r2_caseD_resampled(iter) = temp_mdl_caseD.Rsquared.Adjusted;
            beta0_caseD_resampled(iter) = temp_mdl_caseD.Coefficients.Estimate(1);
            beta1_caseD_resampled(iter) = temp_mdl_caseD.Coefficients.Estimate(2);
            beta2_caseD_resampled(iter) = temp_mdl_caseD.Coefficients.Estimate(3);
            beta3_caseD_resampled(iter) = temp_mdl_caseD.Coefficients.Estimate(4);
            
            
            % Case D2
            x = [temp_x1_CME_CF_zscored,temp_x2_CME_CF_zscored];
            y = temp_y2_CME_CF_zscored;
            
            temp_mdl_caseD2 = fitglm(x,y,'interactions');
            r2_caseD2_resampled(iter) = temp_mdl_caseD2.Rsquared.Adjusted;
            beta0_caseD2_resampled(iter) = temp_mdl_caseD2.Coefficients.Estimate(1);
            beta1_caseD2_resampled(iter) = temp_mdl_caseD2.Coefficients.Estimate(2);
            beta2_caseD2_resampled(iter) = temp_mdl_caseD2.Coefficients.Estimate(3);
            beta3_caseD2_resampled(iter) = temp_mdl_caseD2.Coefficients.Estimate(4);
            
            % Case D2A
            x = temp_x1_CME_CF_zscored;
            y = temp_y2_CME_CF_zscored;
            
            temp_mdl_caseD2A = fitglm(x,y,'linear');
            r2_caseD2A_resampled(iter) = temp_mdl_caseD2A.Rsquared.Adjusted;
            beta0_caseD2A_resampled(iter) = temp_mdl_caseD2A.Coefficients.Estimate(1);
            beta1_caseD2A_resampled(iter) = temp_mdl_caseD2A.Coefficients.Estimate(2);
            
            % Case D2B
            x = temp_x2_CME_CF_zscored;
            y = temp_y2_CME_CF_zscored;
            
            temp_mdl_caseD2B = fitglm(x,y,'linear');
            r2_caseD2B_resampled(iter) = temp_mdl_caseD2B.Rsquared.Adjusted;
            beta0_caseD2B_resampled(iter) = temp_mdl_caseD2B.Coefficients.Estimate(1);
            beta1_caseD2B_resampled(iter) = temp_mdl_caseD2B.Coefficients.Estimate(2);
            
            
            % CMC
            temp_trialIndex_CMC_resample = sort(randi(temp_trialNum_CMC,1,temp_trialNum_CMC));
            
            tempNONNAN_meta_CMC_trialLevel_baseline_resampled = tempNONNAN_meta_CMC_trialLevel_baseline(temp_trialIndex_CMC_resample); %#ok<*PFBNS>
            tempNONNAN_memoryPrecision_CMC_trialLevel_resampled = tempNONNAN_memoryPrecision_CMC_trialLevel(temp_trialIndex_CMC_resample);
            tempNONNAN_meta_CMC_trialLevel_resampled = tempNONNAN_meta_CMC_trialLevel(temp_trialIndex_CMC_resample);
            
            temp_x1_CMC = tempNONNAN_memoryPrecision_CMC_trialLevel_resampled;
            temp_x2_CMC = tempNONNAN_meta_CMC_trialLevel_baseline_resampled;
            temp_y_CMC = tempNONNAN_meta_CMC_trialLevel_resampled;
            
            temp_x1_CMC_zscored = (temp_x1_CMC-mean(temp_x1_CMC))./std(temp_x1_CMC);
            temp_x2_CMC_zscored = (temp_x2_CMC-mean(temp_x2_CMC))./std(temp_x2_CMC);
            temp_y_CMC_zscored = (temp_y_CMC-mean(temp_y_CMC))./std(temp_y_CMC);
            
            
            % Case E
            x = [temp_x1_CMC_zscored,temp_x2_CMC_zscored];
            y = temp_y_CMC_zscored;
            
            %temp_mdl_caseE = fitglm(x,y,'linear');
            temp_mdl_caseE = fitglm(x,y,'interactions');
            r2_caseE_resampled(iter) = temp_mdl_caseE.Rsquared.Adjusted;
            beta0_caseE_resampled(iter) = temp_mdl_caseE.Coefficients.Estimate(1);
            beta1_caseE_resampled(iter) = temp_mdl_caseE.Coefficients.Estimate(2);
            beta2_caseE_resampled(iter) = temp_mdl_caseE.Coefficients.Estimate(3);
            beta3_caseE_resampled(iter) = temp_mdl_caseE.Coefficients.Estimate(4);
            
            
            % CM
            temp_trialIndex_CM_resample = sort(randi(temp_trialNum_CM,1,temp_trialNum_CM));
            
            tempNONNAN_meta_CM_trialLevel_baseline_resampled = tempNONNAN_meta_CM_trialLevel_baseline(temp_trialIndex_CM_resample); %#ok<*PFBNS>
            tempNONNAN_memoryPrecision_CM_trialLevel_resampled = tempNONNAN_memoryPrecision_CM_trialLevel(temp_trialIndex_CM_resample);
            tempNONNAN_meta_CM_trialLevel_resampled = tempNONNAN_meta_CM_trialLevel(temp_trialIndex_CM_resample);
            
            temp_x1_CM = tempNONNAN_memoryPrecision_CM_trialLevel_resampled;
            temp_x2_CM = tempNONNAN_meta_CM_trialLevel_baseline_resampled;
            temp_y_CM = tempNONNAN_meta_CM_trialLevel_resampled;
            
            temp_x1_CM_zscored = (temp_x1_CM-mean(temp_x1_CM))./std(temp_x1_CM);
            temp_x2_CM_zscored = (temp_x2_CM-mean(temp_x2_CM))./std(temp_x2_CM);
            temp_y_CM_zscored = (temp_y_CM-mean(temp_y_CM))./std(temp_y_CM);
            
            
            % Case E2
            x = [temp_x1_CM_zscored,temp_x2_CM_zscored];
            y = temp_y_CM_zscored;
            
            %temp_mdl_caseE2 = fitglm(x,y,'linear');
            temp_mdl_caseE2 = fitglm(x,y,'interactions');
            r2_caseE2_resampled(iter) = temp_mdl_caseE2.Rsquared.Adjusted;
            beta0_caseE2_resampled(iter) = temp_mdl_caseE2.Coefficients.Estimate(1);
            beta1_caseE2_resampled(iter) = temp_mdl_caseE2.Coefficients.Estimate(2);
            beta2_caseE2_resampled(iter) = temp_mdl_caseE2.Coefficients.Estimate(3);
            beta3_caseE2_resampled(iter) = temp_mdl_caseE2.Coefficients.Estimate(4);
            
            
            % CF
            temp_trialIndex_CF_resample = sort(randi(temp_trialNum_CF,1,temp_trialNum_CF));
            
            tempNONNAN_meta_CF_trialLevel_baseline_resampled = tempNONNAN_meta_CF_trialLevel_baseline(temp_trialIndex_CF_resample); %#ok<*PFBNS>
            tempNONNAN_memoryPrecision_CF_trialLevel_resampled = tempNONNAN_memoryPrecision_CF_trialLevel(temp_trialIndex_CF_resample);
            tempNONNAN_meta_CF_trialLevel_resampled = tempNONNAN_meta_CF_trialLevel(temp_trialIndex_CF_resample);
            
            temp_x1_CF = tempNONNAN_memoryPrecision_CF_trialLevel_resampled;
            temp_x2_CF = tempNONNAN_meta_CF_trialLevel_baseline_resampled;
            temp_y_CF = tempNONNAN_meta_CF_trialLevel_resampled;
            
            temp_x1_CF_zscored = (temp_x1_CF-mean(temp_x1_CF))./std(temp_x1_CF);
            temp_x2_CF_zscored = (temp_x2_CF-mean(temp_x2_CF))./std(temp_x2_CF);
            temp_y_CF_zscored = (temp_y_CF-mean(temp_y_CF))./std(temp_y_CF);
            
            
            % Case F
            x = [temp_x1_CF_zscored,temp_x2_CF_zscored];
            y = temp_y_CF_zscored;
            
            %temp_mdl_caseF = fitglm(x,y,'linear');
            temp_mdl_caseF = fitglm(x,y,'interactions');
            r2_caseF_resampled(iter) = temp_mdl_caseF.Rsquared.Adjusted;
            beta0_caseF_resampled(iter) = temp_mdl_caseF.Coefficients.Estimate(1);
            beta1_caseF_resampled(iter) = temp_mdl_caseF.Coefficients.Estimate(2);
            beta2_caseF_resampled(iter) = temp_mdl_caseF.Coefficients.Estimate(3);
            beta3_caseF_resampled(iter) = temp_mdl_caseF.Coefficients.Estimate(4);
            
            
            
            % CME
            temp_trialIndex_CME_resample = sort(randi(temp_trialNum_CME,1,temp_trialNum_CME));
            
            tempNONNAN_meta_CME_trialLevel_baseline_resampled = tempNONNAN_meta_CME_trialLevel_baseline(temp_trialIndex_CME_resample); %#ok<*PFBNS>
            tempNONNAN_memoryPrecision_CME_trialLevel_resampled = tempNONNAN_memoryPrecision_CME_trialLevel(temp_trialIndex_CME_resample);
            tempNONNAN_meta_CME_trialLevel_resampled = tempNONNAN_meta_CME_trialLevel(temp_trialIndex_CME_resample);
            
            temp_x1_CME = tempNONNAN_memoryPrecision_CME_trialLevel_resampled;
            temp_x2_CME = tempNONNAN_meta_CME_trialLevel_baseline_resampled;
            temp_y_CME = tempNONNAN_meta_CME_trialLevel_resampled;
            
            temp_x1_CME_zscored = (temp_x1_CME-mean(temp_x1_CME))./std(temp_x1_CME);
            temp_x2_CME_zscored = (temp_x2_CME-mean(temp_x2_CME))./std(temp_x2_CME);
            temp_y_CME_zscored = (temp_y_CME-mean(temp_y_CME))./std(temp_y_CME);
            
            
            % Case G
            x = [temp_x1_CME_zscored,temp_x2_CME_zscored];
            y = temp_y_CME_zscored;
            
            %temp_mdl_caseG = fitglm(x,y,'linear');
            temp_mdl_caseG = fitglm(x,y,'interactions');
            r2_caseG_resampled(iter) = temp_mdl_caseG.Rsquared.Adjusted;
            beta0_caseG_resampled(iter) = temp_mdl_caseG.Coefficients.Estimate(1);
            beta1_caseG_resampled(iter) = temp_mdl_caseG.Coefficients.Estimate(2);
            beta2_caseG_resampled(iter) = temp_mdl_caseG.Coefficients.Estimate(3);
            beta3_caseG_resampled(iter) = temp_mdl_caseG.Coefficients.Estimate(4);
            
            
            % HBSL
            temp_trialIndex_HBSL_resample = sort(randi(temp_trialNum_HBSL,1,temp_trialNum_HBSL));
            
            tempNONNAN_meta_HBSL_trialLevel_baseline_resampled = tempNONNAN_meta_HBSL_trialLevel_baseline(temp_trialIndex_HBSL_resample); %#ok<*PFBNS>
            tempNONNAN_memoryPrecision_HBSL_trialLevel_resampled = tempNONNAN_memoryPrecision_HBSL_trialLevel(temp_trialIndex_HBSL_resample);
            tempNONNAN_meta_HBSL_trialLevel_resampled = tempNONNAN_meta_HBSL_trialLevel(temp_trialIndex_HBSL_resample);
            
            temp_x1_HBSL = tempNONNAN_memoryPrecision_HBSL_trialLevel_resampled;
            temp_x2_HBSL = tempNONNAN_meta_HBSL_trialLevel_baseline_resampled;
            temp_y_HBSL = tempNONNAN_meta_HBSL_trialLevel_resampled;
            
            temp_x1_HBSL_zscored = (temp_x1_HBSL-mean(temp_x1_HBSL))./std(temp_x1_HBSL);
            temp_x2_HBSL_zscored = (temp_x2_HBSL-mean(temp_x2_HBSL))./std(temp_x2_HBSL);
            temp_y_HBSL_zscored = (temp_y_HBSL-mean(temp_y_HBSL))./std(temp_y_HBSL);
            
            
            % Case H
            x = [temp_x1_HBSL_zscored,temp_x2_HBSL_zscored];
            y = temp_y_HBSL_zscored;
            
            %temp_mdl_caseH = fitglm(x,y,'linear');
            temp_mdl_caseH = fitglm(x,y,'interactions');
            r2_caseH_resampled(iter) = temp_mdl_caseH.Rsquared.Adjusted;
            beta0_caseH_resampled(iter) = temp_mdl_caseH.Coefficients.Estimate(1);
            beta1_caseH_resampled(iter) = temp_mdl_caseH.Coefficients.Estimate(2);
            beta2_caseH_resampled(iter) = temp_mdl_caseH.Coefficients.Estimate(3);
            beta3_caseH_resampled(iter) = temp_mdl_caseH.Coefficients.Estimate(4);
            
        end
    end
    
end


BIC_caseA_resampled;
BIC_caseB_resampled;
BIC_caseC_resampled;

r2_caseA_resampled;
r2_caseB_resampled;
r2_caseC_resampled;

% fprintf("R2 = [%.3f, %.3f, %.3f].\n",mean(r2_caseA_resampled),mean(r2_caseB_resampled),mean(r2_caseC_resampled));
% if if_BIC0_AIC1 == 0
%     fprintf("BIC = [%.0f, %.0f, %.0f].\n",mean(BIC_caseA_resampled),mean(BIC_caseB_resampled),mean(BIC_caseC_resampled));
% elseif if_BIC0_AIC1 == 1
%     fprintf("AIC = [%.0f, %.0f, %.0f].\n",mean(BIC_caseA_resampled),mean(BIC_caseB_resampled),mean(BIC_caseC_resampled));
% end

%% Cross-validation test
if false
    for iter=1:resampleIterCount_GLM
        %parfor iter=1:resampleIterCount_GLM
        temp_trialIndex_resample = sort(randi(temp_trialNum,1,temp_trialNum));
        
        tempNONNAN_meta_trialLevel_baseline_resampled = tempNONNAN_meta_trialLevel_baseline(temp_trialIndex_resample); %#ok<*PFBNS>
        tempNONNAN_memoryPrecision_trialLevel_resampled = tempNONNAN_memoryPrecision_trialLevel(temp_trialIndex_resample);
        tempNONNAN_meta_trialLevel_resampled = tempNONNAN_meta_trialLevel(temp_trialIndex_resample);
        
        temp_x1 = tempNONNAN_memoryPrecision_trialLevel_resampled;
        temp_x2 = tempNONNAN_meta_trialLevel_baseline_resampled;
        temp_y = tempNONNAN_meta_trialLevel_resampled;
        
        temp_x1_zscored = (temp_x1-mean(temp_x1))./std(temp_x1);
        temp_x2_zscored = (temp_x2-mean(temp_x2))./std(temp_x2);
        temp_y_zscored = (temp_y-mean(temp_y))./std(temp_y);
        
        
        % Case C
        %x = [temp_x1,temp_x2];
        %y = temp_y;
        
        x = [temp_x1_zscored,temp_x2_zscored];
        y = temp_y_zscored;
        
        %temp_mdl_caseC = fitglm(x,y,'linear');
        temp_mdl_caseC = fitglm(x,y,'interactions');
        r2_caseC_resampled(iter) = temp_mdl_caseC.Rsquared.Adjusted;
        beta0_caseC_resampled(iter) = temp_mdl_caseC.Coefficients.Estimate(1);
        beta1_caseC_resampled(iter) = temp_mdl_caseC.Coefficients.Estimate(2);
        beta2_caseC_resampled(iter) = temp_mdl_caseC.Coefficients.Estimate(3);
        beta3_caseC_resampled(iter) = temp_mdl_caseC.Coefficients.Estimate(4);
        
    end
    median(r2_caseC_resampled)
end



%% test
if false
    
    x = tempNONNAN_meta_trialLevel_baseline;
    y = tempNONNAN_memoryPrecision_trialLevel;
    temp_mdl_caseX = fitglm(x,y,'linear');
    r2_caseX = temp_mdl_caseX.Rsquared.Adjusted;
    beta0_caseX = temp_mdl_caseX.Coefficients.Estimate(1);
    beta1_caseX = temp_mdl_caseX.Coefficients.Estimate(2);
    
    
end

%%

[~,pAB,~,~] = ttest(r2_caseA_resampled,r2_caseB_resampled);
[~,pAC,~,~] = ttest(r2_caseA_resampled,r2_caseC_resampled);
[~,pBC,~,~] = ttest(r2_caseB_resampled,r2_caseC_resampled);

[~,pA,~,~] = ttest(r2_caseA_resampled,0);
[~,pB,~,~] = ttest(r2_caseB_resampled,0);
[~,pC,~,~] = ttest(r2_caseC_resampled,0);


[~,pAB_BIC,~,~] = ttest(BIC_caseA_resampled,BIC_caseB_resampled);
[~,pAC_BIC,~,~] = ttest(BIC_caseA_resampled,BIC_caseC_resampled);
[~,pBC_BIC,~,~] = ttest(BIC_caseB_resampled,BIC_caseC_resampled);


[~,p12_bothBeta,~,~] = ttest(beta1_caseC_resampled,beta2_caseC_resampled);

[~,p1_bothBeta,~,~] = ttest(beta1_caseC_resampled,0,'tail','right');
[~,p2_bothBeta,~,~] = ttest(beta2_caseC_resampled,0,'tail','right');
[~,p3_bothBeta,~,~] = ttest(beta3_caseC_resampled,0,'tail','right');


[~,p12_caseD,~,~] = ttest(beta1_caseD_resampled,beta2_caseD_resampled);

[~,p1_caseD,~,~] = ttest(beta1_caseD_resampled,0,'tail','right');
[~,p2_caseD,~,~] = ttest(beta2_caseD_resampled,0,'tail','right');
[~,p3_caseD,~,~] = ttest(beta3_caseD_resampled,0,'tail','right');

[~,p12_caseD2,~,~] = ttest(beta1_caseD2_resampled,beta2_caseD2_resampled);

[~,p1_caseD2,~,~] = ttest(beta1_caseD2_resampled,0,'tail','right');
[~,p2_caseD2,~,~] = ttest(beta2_caseD2_resampled,0,'tail','right');
[~,p3_caseD2,~,~] = ttest(beta3_caseD2_resampled,0,'tail','right');

[~,p12_caseD2AB,~,~] = ttest(r2_caseD2A_resampled,r2_caseD2B_resampled);


[~,p12_caseE,~,~] = ttest(beta1_caseE_resampled,beta2_caseE_resampled);

[~,p1_caseE,~,~] = ttest(beta1_caseE_resampled,0,'tail','right');
[~,p2_caseE,~,~] = ttest(beta2_caseE_resampled,0,'tail','right');
[~,p3_caseE,~,~] = ttest(beta3_caseE_resampled,0,'tail','right');


[~,p12_caseE2,~,~] = ttest(beta1_caseE2_resampled,beta2_caseE2_resampled);

[~,p1_caseE2,~,~] = ttest(beta1_caseE2_resampled,0,'tail','right');
[~,p2_caseE2,~,~] = ttest(beta2_caseE2_resampled,0,'tail','right');
[~,p3_caseE2,~,~] = ttest(beta3_caseE2_resampled,0,'tail','right');


[~,p12_caseF,~,~] = ttest(beta1_caseF_resampled,beta2_caseF_resampled);

[~,p1_caseF,~,~] = ttest(beta1_caseF_resampled,0,'tail','right');
[~,p2_caseF,~,~] = ttest(beta2_caseF_resampled,0,'tail','right');
[~,p3_caseF,~,~] = ttest(beta3_caseF_resampled,0,'tail','right');


[~,p12_caseG,~,~] = ttest(beta1_caseG_resampled,beta2_caseG_resampled);

[~,p1_caseG,~,~] = ttest(beta1_caseG_resampled,0,'tail','right');
[~,p2_caseG,~,~] = ttest(beta2_caseG_resampled,0,'tail','right');
[~,p3_caseG,~,~] = ttest(beta3_caseG_resampled,0,'tail','right');


[~,p12_caseH,~,~] = ttest(beta1_caseH_resampled,beta2_caseH_resampled);

[~,p1_caseH,~,~] = ttest(beta1_caseH_resampled,0,'tail','right');
[~,p2_caseH,~,~] = ttest(beta2_caseH_resampled,0,'tail','right');
[~,p3_caseH,~,~] = ttest(beta3_caseH_resampled,0,'tail','right');


%% Plot
if if_plot == 1
    close all
    
    
    %% 3 regression A
    if true
        fig = figure('Name','locDistri','NumberTitle','off');
        %set(gcf,'Position',[50+0 50+0 240 252*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 264*1.04 240*1.04]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 264*1.06 480*1.01]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 264*1.06 480*1.01*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 264*1.06*1.22*0.99 480*1.01*0.9*0.86*0.97]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 264*1.06*1.22*0.99*0.73 480*1.01*0.9*0.86*0.97]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 246.7*0.78 364.0*0.93]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 246.7*0.78*0.9 364.0*0.93]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 246.7*0.78*0.9 364.0*0.93*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[50+0 50+0 246.7*0.78*0.9*0.9*0.95 364.0*0.93*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_p12 = pAB;
        temp_p13 = pBC;
        temp_p23 = pAC;
        
        temp_p1 = pB;
        temp_p2 = pA;
        temp_p3 = pC;
        
        temp_1 = r2_caseB_resampled;
        temp_2 = r2_caseA_resampled;
        temp_3 = r2_caseC_resampled;
        
        
        %temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        %temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        %temp3_SEM = std(temp_3)/sqrt(length(temp_3));
        
        temp1_SEM = std(temp_1);
        temp2_SEM = std(temp_2);
        temp3_SEM = std(temp_3);
        
        temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM,mean(temp_3)-temp3_SEM]);
        temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM]);
        
        temp_y_max12 = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
        temp_y_max13 = max([mean(temp_1)+temp1_SEM,mean(temp_3)+temp3_SEM]);
        temp_y_max23 = max([mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM]);
        
        
        %bar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)], ...
        %    'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
        bar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)], ...
            'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0],'LineWidth',1);
        hold on
        errorbar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)],[temp1_SEM temp2_SEM temp3_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 10);%15
        
        % temp_p12
        tempTxt = sprintf('');
        if temp_p12 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p12 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p12 < 0.05
            tempTxt = sprintf('*');
        end
        %text(0.5,temp_y_min+(temp_y_max12-temp_y_min)*1.06,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        %    'HorizontalAlignment','center');
        %plot([0 1],temp_y_min+(temp_y_max12-temp_y_min)*1.05*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        %hold on
        
        % temp_p13
        tempTxt = sprintf('');
        if temp_p13 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p13 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p13 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y_min+(temp_y_max13-temp_y_min)*1.16,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([0 2],temp_y_min+(temp_y_max13-temp_y_min)*1.15*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        % temp_p23
        tempTxt = sprintf('');
        if temp_p23 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p23 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p23 < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_min+(temp_y_max23-temp_y_min)*1.06,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1 2],temp_y_min+(temp_y_max23-temp_y_min)*1.05*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        % temp_p1
        tempTxt = sprintf('');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        %text(0,temp_y_min-(temp_y_max-temp_y_min)*0.10,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        %    'HorizontalAlignment','center');
        
        % temp_p2
        tempTxt = sprintf('');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        %text(1,temp_y_min-(temp_y_max-temp_y_min)*0.10,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        %    'HorizontalAlignment','center');
        
        % temp_p3
        tempTxt = sprintf('');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        %text(2,temp_y_min-(temp_y_max-temp_y_min)*0.10,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        %    'HorizontalAlignment','center');
        
        
        set(gca,'linewidth',1.5)
        
        %ylim([max(0,temp_y_min-(temp_y_max-temp_y_min)*0.25) temp_y_max+(temp_y_max-temp_y_min)*0.4]);
        ylim([max(0,temp_y_min-(temp_y_max-temp_y_min)*0.25) temp_y_max+(temp_y_max-temp_y_min)*0.25]);
        %ylim([0 1]);
        %ylim([0 1.1]);
        xlim([-0.65 2.65]);
        set(gca, 'FontSize', 8) %14
        %set(gca,'XTickLabel', ["Memory"; "Baseline"; "Both"],'FontSize', 12);%给坐标加标签
        %set(gca,'YTick',0:0.2:1,'FontSize', 12);%给坐标加标签
        
        
        
        %xtl = ["Memory", "Baseline", "Both"];
        %xtl = ["Precision", "Baseline", "Both"];
        xtl = ["Baseline","WM strength","Both"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %         if if_monkey_D0_Z1 == 0
        %             %xtext_yp=(yt(1))*ones(1,length(xt))-0.03;%0.02
        %             %xtext_yp=(yt(1))*ones(1,length(xt))-0.04;%0.02
        %             xtext_yp=(yt(1))*ones(1,length(xt))-0.01;%0.02
        %         elseif if_monkey_D0_Z1 == 1
        %             %xtext_yp=(yt(1))*ones(1,length(xt))-0.03;%0.02
        %             %xtext_yp=(yt(1))*ones(1,length(xt))-0.006;%0.02
        %             xtext_yp=(yt(1))*ones(1,length(xt))-0.02;%0.02
        %         end
        
        if if_monkey_D0_Z1 == 0
            xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.205;%0.155
        elseif if_monkey_D0_Z1 == 1
            xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.175;%0.155
        end
        
        %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%9
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9
        
        set(gca,'xticklabel','');
        
        
        temp_stepNum = 6;
        
        %temptemp_range = linspace(temp_y_min,temp_y_max,4);
        temptemp_range = linspace(0,temp_y_max,temp_stepNum);
        temptemp_range2 = ceil((temptemp_range*100))/100;
        temptemp_step = temptemp_range2(2)-temptemp_range2(1);
        %     temptemp_range3 = [temptemp_range2(1)-temptemp_step,temptemp_range2,temptemp_range2(end)+temptemp_step];
        %     set(gca,'ytick',temptemp_range3);
        %set(gca,'ytick',temptemp_range2);
        set(gca,'ytick',0:temptemp_step:temptemp_step*temp_stepNum);
        
        if if_monkey_D0_Z1 == 0 && currentSessionIndex_AB == 8
            yticks([0:0.1:0.3]); %#ok<*NBRAK>
        end
        
        set(gca,'box','off');% 取消右、上边框
        %ylabel('r2', 'FontSize', 10, 'FontWeight', 'normal');
        ylabel('Explained variance', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf('Meta-WM regression'),'fontsize',8.5);%9
        %temp_title.Interpreter = 'none';
        
    end
    
    %% 3 regression A, BIC
    if true
        fig = figure('Name','locDistri','NumberTitle','off');
        %set(gcf,'Position',[50+0 50+0 246.7*0.78*0.9*0.9*0.95 364.0*0.93*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 246.7*0.78*0.9*0.9*0.95*1.2 364.0*0.93*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[50+0 50+0 246.7*0.78*0.9*0.9*0.95*1.2*0.8 364.0*0.93*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_p12 = pAB_BIC;
        temp_p13 = pBC_BIC;
        temp_p23 = pAC_BIC;
        
        temp_1 = BIC_caseB_resampled;
        temp_2 = BIC_caseA_resampled;
        temp_3 = BIC_caseC_resampled;
        
        
        %temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        %temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        %temp3_SEM = std(temp_3)/sqrt(length(temp_3));
        
        temp1_SEM = std(temp_1);
        temp2_SEM = std(temp_2);
        temp3_SEM = std(temp_3);
        
        temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM,mean(temp_3)-temp3_SEM]);
        temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM]);
        
        temp_y_min = 3060-5-1-0.45-0.3*2;
        temp_y_max = 3360+2+0.5+0.45+0.3*3;
        
        temp_y_max12 = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
        temp_y_max13 = max([mean(temp_1)+temp1_SEM,mean(temp_3)+temp3_SEM]);
        temp_y_max23 = max([mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM]);
        
        
        %bar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)],'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0],'LineWidth',1);
        bar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)],'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0],'LineWidth',0.5);
        hold on
        %errorbar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)],[temp1_SEM temp2_SEM temp3_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 10);%15,10
        errorbar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)],[temp1_SEM temp2_SEM temp3_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 1.5, 'CapSize', 8);%15,10
        
        % temp_p12
        tempTxt = sprintf('');
        if temp_p12 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p12 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p12 < 0.05
            tempTxt = sprintf('*');
        end
        %text(0.5,temp_y_min+(temp_y_max12-temp_y_min)*1.06,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        %    'HorizontalAlignment','center');
        %plot([0 1],temp_y_min+(temp_y_max12-temp_y_min)*1.05*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        %hold on
        
        % temp_p13
        tempTxt = sprintf('');
        if temp_p13 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p13 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p13 < 0.05
            tempTxt = sprintf('*');
        end
        %text(1,temp_y_min+(temp_y_max13-temp_y_min)*1.16,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        %    'HorizontalAlignment','center');
        %plot([0 2],temp_y_min+(temp_y_max13-temp_y_min)*1.13*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        text(1,temp_y_min+(temp_y_max-temp_y_min)*1.16,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
           'HorizontalAlignment','center');
        plot([0 2],temp_y_min+(temp_y_max-temp_y_min)*1.13*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);        
        hold on
        
        % temp_p23
        tempTxt = sprintf('');
        if temp_p23 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p23 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p23 < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_min+(temp_y_max23-temp_y_min)*1.31,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1 2],temp_y_min+(temp_y_max23-temp_y_min)*1.23*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        
        set(gca,'linewidth',1.5)
        
        %ylim([max(0,temp_y_min-(temp_y_max-temp_y_min)*0.25) temp_y_max+(temp_y_max-temp_y_min)*0.4]);
        ylim([max(0,temp_y_min-(temp_y_max-temp_y_min)*0.25) temp_y_max+(temp_y_max-temp_y_min)*0.25]);
        %ylim([0 1]);
        %ylim([0 1.1]);
        xlim([-0.65 2.65]);
        set(gca, 'FontSize', 8) %14
        %set(gca,'XTickLabel', ["Memory"; "Baseline"; "Both"],'FontSize', 12);%给坐标加标签
        %set(gca,'YTick',0:0.2:1,'FontSize', 12);%给坐标加标签
        
        
        
        xtl = ["Baseline","WM strength","Both"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        
        
        if if_monkey_D0_Z1 == 0
            xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.305;%0.155
        elseif if_monkey_D0_Z1 == 1
            xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.305;%0.175
        end
        
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9,8,7.5
        
        set(gca,'xticklabel','');
        
        
        %         temp_stepNum = 6;
        %
        %         temptemp_range = linspace(0,temp_y_max,temp_stepNum);
        %         temptemp_range2 = ceil((temptemp_range*100))/100;
        %         temptemp_step = temptemp_range2(2)-temptemp_range2(1);
        %         set(gca,'ytick',0:temptemp_step:temptemp_step*temp_stepNum);
        %
        %         if if_monkey_D0_Z1 == 0 && currentSessionIndex_AB == 8
        %             yticks([0:0.1:0.3]); %#ok<*NBRAK>
        %         end
        
        set(gca,'box','off');% 取消右、上边框
        if if_BIC0_AIC1 == 0
            ylabel('BIC', 'FontSize', 9, 'FontWeight', 'normal');
        elseif if_BIC0_AIC1 == 1
            ylabel('AIC', 'FontSize', 9, 'FontWeight', 'normal');
        end
        %temp_title = title(sprintf('Meta-WM regression'),'fontsize',8.5);%9
        %temp_title.Interpreter = 'none';
        
    end
    
    
    
    %% 3 regression B
    if false
        fig = figure('Name','locDistri','NumberTitle','off');
        %set(gcf,'Position',[50+0 50+0 240 252*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 264*1.04 240*1.04]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 264*1.06 480*1.01]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 264*1.06 480*1.01*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 264*1.06*1.22*0.99 480*1.01*0.9*0.86*0.97]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 264*1.06*1.22*0.99*0.73 480*1.01*0.9*0.86*0.97]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 246.7*0.78 364.0*0.93]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 246.7*0.78*0.9 364.0*0.93]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[50+0 50+0 246.7*0.78*0.9 364.0*0.93*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[50+0 50+0 246.7*0.78*0.9*0.9*0.95*0.95 364.0*0.93*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_p12 = pAB;
        temp_p13 = pBC;
        temp_p23 = pAC;
        
        temp_p1 = pB;
        temp_p2 = pA;
        temp_p3 = pC;
        
        temp_1 = r2_caseB_resampled;
        temp_2 = r2_caseA_resampled;
        temp_3 = r2_caseC_resampled;
        
        
        %temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        %temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        %temp3_SEM = std(temp_3)/sqrt(length(temp_3));
        
        temp1_SEM = std(temp_1);
        temp2_SEM = std(temp_2);
        temp3_SEM = std(temp_3);
        
        temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM,mean(temp_3)-temp3_SEM]);
        temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM]);
        
        temp_y_max12 = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
        temp_y_max13 = max([mean(temp_1)+temp1_SEM,mean(temp_3)+temp3_SEM]);
        temp_y_max23 = max([mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM]);
        
        
        bar([0 1 2], [mean(temp_3) mean(temp_1) mean(temp_2)], 0.75, ...
            'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0],'LineWidth',1);
        hold on
        errorbar([0 1 2], [mean(temp_3) mean(temp_1) mean(temp_2)],[temp3_SEM temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 1.5, 'CapSize', 8);%15
        
        
        % temp_p13
        tempTxt = sprintf('');
        if temp_p13 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p13 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p13 < 0.05
            tempTxt = sprintf('*');
        end
        text(0.5,temp_y_min+(temp_y_max13-temp_y_min)*1.06,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([0 1],temp_y_min+(temp_y_max13-temp_y_min)*1.05*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        % temp_p23
        tempTxt = sprintf('');
        if temp_p23 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p23 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p23 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y_min+(temp_y_max23-temp_y_min)*1.16,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([0 2],temp_y_min+(temp_y_max23-temp_y_min)*1.15*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        
        set(gca,'linewidth',1.5)
        
        %ylim([max(0,temp_y_min-(temp_y_max-temp_y_min)*0.25) temp_y_max+(temp_y_max-temp_y_min)*0.4]);
        ylim([max(0,temp_y_min-(temp_y_max-temp_y_min)*0.25) temp_y_max+(temp_y_max-temp_y_min)*0.25]);
        %ylim([0 1]);
        %ylim([0 1.1]);
        xlim([-0.65 2.65]);
        set(gca, 'FontSize', 8) %14
        %set(gca,'XTickLabel', ["Memory"; "Baseline"; "Both"],'FontSize', 12);%给坐标加标签
        %set(gca,'YTick',0:0.2:1,'FontSize', 12);%给坐标加标签
        
        
        
        %xtl = ["Memory", "Baseline", "Both"];
        %xtl = ["Precision", "Baseline", "Both"];
        xtl = ["Both","Baseline","WM strength"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %         if if_monkey_D0_Z1 == 0
        %             %xtext_yp=(yt(1))*ones(1,length(xt))-0.03;%0.02
        %             %xtext_yp=(yt(1))*ones(1,length(xt))-0.04;%0.02
        %             xtext_yp=(yt(1))*ones(1,length(xt))-0.01;%0.02
        %         elseif if_monkey_D0_Z1 == 1
        %             %xtext_yp=(yt(1))*ones(1,length(xt))-0.03;%0.02
        %             %xtext_yp=(yt(1))*ones(1,length(xt))-0.006;%0.02
        %             xtext_yp=(yt(1))*ones(1,length(xt))-0.02;%0.02
        %         end
        
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.155;
        
        %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%9
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9
        
        set(gca,'xticklabel','');
        
        
        temp_stepNum = 6;
        
        %temptemp_range = linspace(temp_y_min,temp_y_max,4);
        temptemp_range = linspace(0,temp_y_max,temp_stepNum);
        temptemp_range2 = ceil((temptemp_range*100))/100;
        temptemp_step = temptemp_range2(2)-temptemp_range2(1);
        %     temptemp_range3 = [temptemp_range2(1)-temptemp_step,temptemp_range2,temptemp_range2(end)+temptemp_step];
        %     set(gca,'ytick',temptemp_range3);
        %set(gca,'ytick',temptemp_range2);
        set(gca,'ytick',0:temptemp_step:temptemp_step*temp_stepNum);
        
        if if_monkey_D0_Z1 == 0 && currentSessionIndex_AB == 8
            %yticks([0:0.1:0.3]); %#ok<*NBRAK>
            yticks([0.1:0.1:0.3]); %#ok<*NBRAK>
        end
        
        set(gca,'box','off');% 取消右、上边框
        %ylabel('r2', 'FontSize', 10, 'FontWeight', 'normal');
        ylabel('Explained variance', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf('Meta-WM regression'),'fontsize',8.5);%9
        %temp_title.Interpreter = 'none';
        
    end
    
    
    
    %% One regression (both)
    if true
        
        temp_chanceLevel = 0;
        
        fig = figure('Name','locDistri','NumberTitle','off');
        %set(gcf,'Position',[150+0 50+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[150+0 50+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05*1.05*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[150+0 50+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05*1.05*1.05*0.9*0.95 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[150+0 50+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05*1.05*1.05*0.9*0.95*0.95 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[150+0 50+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05*1.05*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[150+0 50+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05*1.05*1.05*0.80 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_p12 = p12_bothBeta;
        
        temp_p1 = p2_bothBeta;
        temp_p2 = p1_bothBeta;
        temp_p3 = p3_bothBeta;
        
        temp_1 = beta2_caseC_resampled;
        temp_2 = beta1_caseC_resampled;
        temp_3 = beta3_caseC_resampled;
        
        temp_r2 = mean(r2_caseC_resampled);
        
        
        temp_y_min = min([temp_1 temp_2 temp_3]); %#ok<*NASGU>
        temp_y_max = max([temp_1 temp_2 temp_3]);
        
        temp_y_min = -0.1;
        temp_y_max = 0.4790;
        
        [~,temp_y1_max] = bounds(temp_1);
        [~,temp_y2_max] = bounds(temp_2);
        [~,temp_y3_max] = bounds(temp_3);
        
        if if_plot_violinplot0_pairline1 == 0
            temp_data = [temp_1';temp_2';temp_3'];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            g3 = repmat({'C'},length(temp_3),1);
            
            temp_label = [g1;g2;g3];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 3, 1);
            
            violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'};{'C'}]);
            
        elseif if_plot_violinplot0_pairline1 == 1
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            
            scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
        end
        
        
        plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.25,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1.0 2.0],temp_y_max+(temp_y_max-temp_y_min)*0.22*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        
        
        tempTxt = sprintf('');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        xticks([1 2 3]);
        
        xtl = ["Baseline", "WM strength", "Interaction"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.155;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9,12,8,7.5
        
        set(gca,'xticklabel','');
        
        
        set(gca,'linewidth',1.5)
        
        xlim([0.35 3.65]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.2
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        
        ylabel('Beta', 'FontSize', 9);
        temp_title = title(sprintf('Meta-WM regression'),'fontsize',8);
        temp_title.Interpreter = 'none';
        %subtitle(sprintf('r2 = %.3f',temp_r2),'fontsize',8);%9
        subtitle(sprintf('all trials, r2 = %.3f',temp_r2),'fontsize',8);%9
        
    end
    
    
    %% One regression (caseD, CME_CF)
    if true
        
        temp_chanceLevel = 0;
        
        fig = figure('Name','locDistri','NumberTitle','off');
        set(gcf,'Position',[350+0 50+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_p12 = p12_caseD;
        
        temp_p1 = p2_caseD;
        temp_p2 = p1_caseD;
        temp_p3 = p3_caseD;
        
        temp_1 = beta2_caseD_resampled;
        temp_2 = beta1_caseD_resampled;
        temp_3 = beta3_caseD_resampled;
        
        temp_r2 = mean(r2_caseD_resampled);
        
        
        temp_y_min = min([temp_1 temp_2 temp_3]); %#ok<*NASGU>
        temp_y_max = max([temp_1 temp_2 temp_3]);
        
        temp_y_min = -0.1;
        temp_y_max = 0.4790;
        
        [~,temp_y1_max] = bounds(temp_1);
        [~,temp_y2_max] = bounds(temp_2);
        [~,temp_y3_max] = bounds(temp_3);
        
        
        if if_plot_violinplot0_pairline1 == 0
            temp_data = [temp_1';temp_2';temp_3'];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            g3 = repmat({'C'},length(temp_3),1);
            
            temp_label = [g1;g2;g3];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 3, 1);
            
            violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'};{'C'}]);
            
        elseif if_plot_violinplot0_pairline1 == 1
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            
            scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
        end
        
        
        plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.25,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1.0 2.0],temp_y_max+(temp_y_max-temp_y_min)*0.22*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        
        
        tempTxt = sprintf('');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        xticks([1 2 3]);
        
        xtl = ["Baseline", "WM strength", "Interaction"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.155;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
        
        set(gca,'xticklabel','');
        
        
        set(gca,'linewidth',1.5)
        
        xlim([0.35 3.65]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.3
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        
        ylabel('Beta', 'FontSize', 9);
        temp_title = title(sprintf('Meta-WM regression'),'fontsize',8);
        temp_title.Interpreter = 'none';
        subtitle(sprintf('low-strength mismatch\n & offload, r2 = %.3f',temp_r2),'fontsize',8);%9
        
    end
    
    
    %% One regression (caseE, CMC)
    if true
        
        temp_chanceLevel = 0;
        
        fig = figure('Name','locDistri','NumberTitle','off');
        set(gcf,'Position',[550+0 50+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_p12 = p12_caseE;
        
        temp_p1 = p2_caseE;
        temp_p2 = p1_caseE;
        temp_p3 = p3_caseE;
        
        temp_1 = beta2_caseE_resampled;
        temp_2 = beta1_caseE_resampled;
        temp_3 = beta3_caseE_resampled;
        
        temp_r2 = mean(r2_caseE_resampled);
        
        
        temp_y_min = min([temp_1 temp_2 temp_3]); %#ok<*NASGU>
        temp_y_max = max([temp_1 temp_2 temp_3]);
        
        temp_y_min = -0.1;
        temp_y_max = 0.4790;
        
        [~,temp_y1_max] = bounds(temp_1);
        [~,temp_y2_max] = bounds(temp_2);
        [~,temp_y3_max] = bounds(temp_3);
        
        
        if if_plot_violinplot0_pairline1 == 0
            temp_data = [temp_1';temp_2';temp_3'];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            g3 = repmat({'C'},length(temp_3),1);
            
            temp_label = [g1;g2;g3];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 3, 1);
            
            violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'};{'C'}]);
            
        elseif if_plot_violinplot0_pairline1 == 1
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            
            scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
        end
        
        
        plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.25,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1.0 2.0],temp_y_max+(temp_y_max-temp_y_min)*0.22*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        
        
        tempTxt = sprintf('');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        xticks([1 2 3]);
        
        xtl = ["Baseline", "WM strength", "Interaction"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.155;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
        
        set(gca,'xticklabel','');
        
        
        set(gca,'linewidth',1.5)
        
        xlim([0.35 3.65]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.3
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        
        ylabel('Beta', 'FontSize', 9);
        temp_title = title(sprintf('Meta-WM regression'),'fontsize',8);
        temp_title.Interpreter = 'none';
        subtitle(sprintf('memory correct\nr2 = %.3f',temp_r2),'fontsize',8);%9
        
    end
    
    
    %% One regression (caseE2, CM)
    if true
        
        temp_chanceLevel = 0;
        
        fig = figure('Name','locDistri','NumberTitle','off');
        set(gcf,'Position',[550+0 50+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_p12 = p12_caseE2;
        
        temp_p1 = p2_caseE2;
        temp_p2 = p1_caseE2;
        temp_p3 = p3_caseE2;
        
        temp_1 = beta2_caseE2_resampled;
        temp_2 = beta1_caseE2_resampled;
        temp_3 = beta3_caseE2_resampled;
        
        temp_r2 = mean(r2_caseE2_resampled);
        
        
        temp_y_min = min([temp_1 temp_2 temp_3]); %#ok<*NASGU>
        temp_y_max = max([temp_1 temp_2 temp_3]);
        
        temp_y_min = -0.1;
        temp_y_max = 0.4790;
        
        [~,temp_y1_max] = bounds(temp_1);
        [~,temp_y2_max] = bounds(temp_2);
        [~,temp_y3_max] = bounds(temp_3);
        
        
        if if_plot_violinplot0_pairline1 == 0
            temp_data = [temp_1';temp_2';temp_3'];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            g3 = repmat({'C'},length(temp_3),1);
            
            temp_label = [g1;g2;g3];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 3, 1);
            
            violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'};{'C'}]);
            
        elseif if_plot_violinplot0_pairline1 == 1
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            
            scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
        end
        
        
        plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.25,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1.0 2.0],temp_y_max+(temp_y_max-temp_y_min)*0.22*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        
        
        tempTxt = sprintf('');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        xticks([1 2 3]);
        
        xtl = ["Baseline", "WM strength", "Interaction"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.155;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
        
        set(gca,'xticklabel','');
        
        
        set(gca,'linewidth',1.5)
        
        xlim([0.35 3.65]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.3
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        
        ylabel('Beta', 'FontSize', 9);
        temp_title = title(sprintf('Meta-WM regression'),'fontsize',8);
        temp_title.Interpreter = 'none';
        %subtitle(sprintf('memory\nr2 = %.3f',temp_r2),'fontsize',8);%9
        subtitle(sprintf('memory, r2 = %.3f',temp_r2),'fontsize',8);%9
        
    end
    
    
    
    %% One regression (caseF, CF)
    if true
        
        temp_chanceLevel = 0;
        
        fig = figure('Name','locDistri','NumberTitle','off');
        set(gcf,'Position',[750+0 50+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_p12 = p12_caseF;
        
        temp_p1 = p2_caseF;
        temp_p2 = p1_caseF;
        temp_p3 = p3_caseF;
        
        temp_1 = beta2_caseF_resampled;
        temp_2 = beta1_caseF_resampled;
        temp_3 = beta3_caseF_resampled;
        
        temp_r2 = mean(r2_caseF_resampled);
        
        
        temp_y_min = min([temp_1 temp_2 temp_3]); %#ok<*NASGU>
        temp_y_max = max([temp_1 temp_2 temp_3]);
        
        temp_y_min = -0.1;
        temp_y_max = 0.4790;
        
        [~,temp_y1_max] = bounds(temp_1);
        [~,temp_y2_max] = bounds(temp_2);
        [~,temp_y3_max] = bounds(temp_3);
        
        
        if if_plot_violinplot0_pairline1 == 0
            temp_data = [temp_1';temp_2';temp_3'];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            g3 = repmat({'C'},length(temp_3),1);
            
            temp_label = [g1;g2;g3];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 3, 1);
            
            violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'};{'C'}]);
            
        elseif if_plot_violinplot0_pairline1 == 1
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            
            scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
        end
        
        
        plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.25,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1.0 2.0],temp_y_max+(temp_y_max-temp_y_min)*0.22*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        
        
        tempTxt = sprintf('');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        xticks([1 2 3]);
        
        xtl = ["Baseline", "WM strength", "Interaction"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.155;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
        
        set(gca,'xticklabel','');
        
        
        set(gca,'linewidth',1.5)
        
        xlim([0.35 3.65]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.3
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        
        ylabel('Beta', 'FontSize', 9);
        temp_title = title(sprintf('Meta-WM regression'),'fontsize',8);
        temp_title.Interpreter = 'none';
        %subtitle(sprintf('offload\nr2 = %.3f',temp_r2),'fontsize',8);%9
        subtitle(sprintf('offload, r2 = %.3f',temp_r2),'fontsize',8);%9
        
    end
    
    
    %% One regression (caseG, CME)
    if true
        
        temp_chanceLevel = 0;
        
        fig = figure('Name','locDistri','NumberTitle','off');
        set(gcf,'Position',[950+0 50+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_p12 = p12_caseG;
        
        temp_p1 = p2_caseG;
        temp_p2 = p1_caseG;
        temp_p3 = p3_caseG;
        
        temp_1 = beta2_caseG_resampled;
        temp_2 = beta1_caseG_resampled;
        temp_3 = beta3_caseG_resampled;
        
        temp_r2 = mean(r2_caseG_resampled);
        
        
        temp_y_min = min([temp_1 temp_2 temp_3]); %#ok<*NASGU>
        temp_y_max = max([temp_1 temp_2 temp_3]);
        
        temp_y_min = -0.1;
        temp_y_max = 0.4790;
        
        [~,temp_y1_max] = bounds(temp_1);
        [~,temp_y2_max] = bounds(temp_2);
        [~,temp_y3_max] = bounds(temp_3);
        
        
        if if_plot_violinplot0_pairline1 == 0
            temp_data = [temp_1';temp_2';temp_3'];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            g3 = repmat({'C'},length(temp_3),1);
            
            temp_label = [g1;g2;g3];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 3, 1);
            
            violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'};{'C'}]);
            
        elseif if_plot_violinplot0_pairline1 == 1
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            
            scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
        end
        
        
        plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.25,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1.0 2.0],temp_y_max+(temp_y_max-temp_y_min)*0.22*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        
        
        tempTxt = sprintf('');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        xticks([1 2 3]);
        
        xtl = ["Baseline", "WM strength", "Interaction"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.155;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
        
        set(gca,'xticklabel','');
        
        
        set(gca,'linewidth',1.5)
        
        xlim([0.35 3.65]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.3
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        
        ylabel('Beta', 'FontSize', 9);
        temp_title = title(sprintf('Meta-WM regression'),'fontsize',8);
        temp_title.Interpreter = 'none';
        subtitle(sprintf('low-strength mismatch\nr2 = %.3f',temp_r2),'fontsize',8);%9
        
    end
    
    
    
    %% One regression (caseH, HBSL)
    if true
        
        temp_chanceLevel = 0;
        
        fig = figure('Name','locDistri','NumberTitle','off');
        set(gcf,'Position',[1150+0 50+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        %temp_p12 = p12_caseH;
        temp_p = p12_caseH;
        
        temp_p1 = p2_caseH;
        temp_p2 = p1_caseH;
        temp_p3 = p3_caseH;
        
        temp_1 = beta2_caseH_resampled;
        temp_2 = beta1_caseH_resampled;
        temp_3 = beta3_caseH_resampled;
        
        temp_r2 = mean(r2_caseH_resampled);
        
        
        temp_y_min = min([temp_1 temp_2 temp_3]); %#ok<*NASGU>
        temp_y_max = max([temp_1 temp_2 temp_3]);
        
        temp_y_min = -0.1;
        temp_y_max = 0.4790;
        
        [~,temp_y1_max] = bounds(temp_1);
        [~,temp_y2_max] = bounds(temp_2);
        [~,temp_y3_max] = bounds(temp_3);
        
        
        if if_plot_violinplot0_pairline1 == 0
            temp_data = [temp_1';temp_2';temp_3'];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            g3 = repmat({'C'},length(temp_3),1);
            
            temp_label = [g1;g2;g3];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 3, 1);
            
            violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'};{'C'}]);
            
        elseif if_plot_violinplot0_pairline1 == 1
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            
            scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
        end
        
        
        plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.25,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1.0 2.0],temp_y_max+(temp_y_max-temp_y_min)*0.22*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        
        
        tempTxt = sprintf('');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        xticks([1 2 3]);
        
        xtl = ["Baseline", "WM strength", "Interaction"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.155;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
        
        set(gca,'xticklabel','');
        
        
        set(gca,'linewidth',1.5)
        
        xlim([0.35 3.65]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.3
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        
        ylabel('Beta', 'FontSize', 9);
        temp_title = title(sprintf('Meta-WM regression'),'fontsize',8);
        temp_title.Interpreter = 'none';
        %subtitle(sprintf('High baseline\nr2 = %.3f',temp_r2),'fontsize',8);%9
        subtitle(sprintf('Low baseline\nr2 = %.3f',temp_r2),'fontsize',8);%9
        %subtitle(sprintf('Extreme baseline\nr2 = %.3f',temp_r2),'fontsize',8);%9
        %subtitle(sprintf('Medium baseline\nr2 = %.3f',temp_r2),'fontsize',8);%9
        %subtitle(sprintf('Extreme WM strength\nr2 = %.3f',temp_r2),'fontsize',8);%9
        %subtitle(sprintf('Medium WM strength\nr2 = %.3f',temp_r2),'fontsize',8);%9
        %subtitle(sprintf('Low 1/3 WM strength\nr2 = %.3f',temp_r2),'fontsize',8);%9
        %subtitle(sprintf('Middle 1/3 WM strength\nr2 = %.3f',temp_r2),'fontsize',8);%9
        %subtitle(sprintf('High 1/3 WM strength\nr2 = %.3f',temp_r2),'fontsize',8);%9
        
    end
    
    
    
end


a = 1;
%%
if if_onlyGLM == 0
    
    %% Baseline & meta, Baseline & mismatch
    % if true
    if exist('if_compute_summary','var') == 0
        if if_monkey_D0_Z1 == 0 && currentSessionIndex_AB == 8
            
            
            memoryPrecision_trialLevel;
            
            
            temptempBoolIndex = choiceBoolIndex' & (~isnan(memoryPrecision_trialLevel));
            
            temp_memoryPrecision = memoryPrecision_trialLevel(temptempBoolIndex);
            temp_meta = meta_trialLevel(temptempBoolIndex);
            temp_baselineMeta = meta_trialLevel_baseline(temptempBoolIndex);
            
            temp_binNum = 100;
            
            temp1 = linspace(0,1,temp_binNum+1);
            temp_step = temp1(2)-temp1(1);
            temp2 = temp1 + temp_step;
            temp_range = [temp1;temp2]';
            temp_range(end,:) = [];
            
            temp_count = nan(1,temp_binNum);
            
            temp12_index_cell = cell(1,temp_binNum);
            
            temp_r = nan(1,temp_binNum);
            
            for tempi=1:temp_binNum
                temptemp_range = temp_range(tempi,:);
                
                temp1 = temp_memoryPrecision>=temptemp_range(1);
                temp2 = temp_memoryPrecision<temptemp_range(2);
                temp12 = temp1 & temp2;
                
                temp_count(tempi) = sum(temp12);
                
                if sum(temp12) <= 1
                    continue
                end
                a = 1;
                
                temp12_index = find(temp12==true);
                
                temp12_index_cell{tempi} = temp12_index;
                
                temptemp_memoryPrecision = temp_memoryPrecision(temp12_index);
                temptemp_meta = temp_meta(temp12_index);
                temptemp_baselineMeta = temp_baselineMeta(temp12_index);
                
                [temp_r(tempi),~] = corr(temptemp_baselineMeta,temptemp_meta);
                
            end
            temp_r_valid = temp_r(~isnan(temp_r));
            
            temp1 = 1:temp_binNum;
            temp_bin_valid = temp1(~isnan(temp_r));
            temp2 = mean(temp_range,2);
            temp_xmesh = temp2(~isnan(temp_r));
            
            [~,temp_p_sameWM] = ttest(temp_r_valid,0); %#ok<*ASGLU>
            
            [M,I] = max(temp_count);
            
            temp_I = 69;
            
            temp12_index = temp12_index_cell{temp_I};
            
            temptemp_baselineMeta = temp_baselineMeta(temp12_index);
            temptemp_meta = temp_meta(temp12_index);
            
            temptemp_range = temp_range(temp_I,:);
            
            if if_plot
                fig = figure('Name','Demo, same strength trials correlation','NumberTitle','off');
                set(gcf,'Position',[510 390 196.4*0.905*0.915*0.93*0.95*0.94 227.5*0.905*0.915*0.93*0.95]);
                
                t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
                
                nexttile
                
                x = temptemp_baselineMeta;
                y = temptemp_meta;
                
                
                [r,p] = corr(x,y);
                
                mdl = fitglm(x,y);
                
                
                temp_size = 10;%10,3
                temp_h = scatter(x, y, ...
                    temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
                hold on
                
                [temp_xmin,temp_xmax] = bounds(x);
                [temp_ymin,temp_ymax] = bounds(y); %#ok<*ASGLU>
                
                x_fit = temp_xmin:0.001:temp_xmax;
                y_fit = predict(mdl,x_fit')';
                
                plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
                hold on
                
                xlim([-0.05 1.05]);
                ylim([-0.05 1.05]);
                
                xticks([0 0.5 1]);
                yticks([0 1]);
                
                
                tempTxt = 'p > 0.05';
                if p < 0.05
                    tempTxt = 'p < 0.05';
                end
                if p < 0.01
                    tempTxt = 'p < 0.01';
                end
                if p < 0.001
                    tempTxt = 'p < 0.001';
                end
                
                set(gca,'linewidth',1.5)
                set(gca, 'FontSize', 8)
                set(gca,'box','off');% 取消右、上边框
                xlabel(sprintf('Baseline meta'), 'FontSize', 9, 'FontWeight', 'normal');
                ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
                
                title(sprintf('WM strength [0.68,0.69]'), 'FontSize', 8, 'FontWeight', 'normal');
                subtitle(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 7, 'FontWeight', 'normal');
                
            end
            
            
            if if_plot
                fig = figure('Name','Same strength trials correlation','NumberTitle','off');
                set(gcf,'Position',[810 390 196.4*0.905*0.915*2*0.75*0.8*1.08*1.04 227.5*0.905*0.915*0.93*0.95]);
                
                t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
                
                nexttile
                
                x = temp_xmesh;
                y = temp_r_valid;
                
                p = temp_p_sameWM;
                r = median(temp_r_valid);
                
                
                temp_size = 10;%10,3
                temp_h = scatter(x, y, ...
                    temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
                hold on
                
                [temp_xmin,temp_xmax] = bounds(x);
                [temp_ymin,temp_ymax] = bounds(y); %#ok<*ASGLU>
                
                
                plot([0 1],0*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
                hold on
                
                
                %         xlim([-0.05 1.05]);
                ylim([temp_ymin-(temp_ymax-temp_ymin)*0.05 temp_ymax+(temp_ymax-temp_ymin)*0.05]);
                %
                %         xticks([0 1]);
                %         yticks([0 1]);
                
                
                tempTxt = 'p > 0.05';
                if p < 0.05
                    tempTxt = 'p < 0.05';
                end
                if p < 0.01
                    tempTxt = 'p < 0.01';
                end
                if p < 0.001
                    tempTxt = 'p < 0.001';
                end
                
                set(gca,'linewidth',1.5)
                set(gca, 'FontSize', 8)
                set(gca,'box','off');% 取消右、上边框
                xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
                ylabel(sprintf('Correlation (r)'), 'FontSize', 9, 'FontWeight', 'normal');
                
                title(sprintf('Within each WM strength bin'), 'FontSize', 8, 'FontWeight', 'normal');
                subtitle(sprintf('r(median) = %.3f, %s',r,tempTxt), 'FontSize', 7, 'FontWeight', 'normal');
                
            end
            
            
            
            a = 1;
        end
    end
    
    
    
    %% Regression (one), using prior & WM strength to predict trial labels (CME or CF)
    if if_plot
        
        temp_chanceLevel = 0;
        
        fig = figure('Name','locDistri','NumberTitle','off');
        set(gcf,'Position',[350+0 50+350 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_p12 = p12_caseD2;
        
        temp_p1 = p2_caseD2;
        temp_p2 = p1_caseD2;
        temp_p3 = p3_caseD2;
        
        temp_1 = beta2_caseD2_resampled;
        temp_2 = beta1_caseD2_resampled;
        temp_3 = beta3_caseD2_resampled;
        
        temp_r2 = mean(r2_caseD2_resampled);
        
        
        temp_y_min = min([temp_1 temp_2 temp_3]); %#ok<*NASGU>
        temp_y_max = max([temp_1 temp_2 temp_3]);
        
        temp_y_min = -0.1;
        temp_y_max = 0.4790;
        
        [~,temp_y1_max] = bounds(temp_1);
        [~,temp_y2_max] = bounds(temp_2);
        [~,temp_y3_max] = bounds(temp_3);
        
        
        if if_plot_violinplot0_pairline1 == 0
            temp_data = [temp_1';temp_2';temp_3'];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            g3 = repmat({'C'},length(temp_3),1);
            
            temp_label = [g1;g2;g3];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 3, 1);
            
            violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'};{'C'}]);
            
        elseif if_plot_violinplot0_pairline1 == 1
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            
            scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
        end
        
        
        plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.25,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1.0 2.0],temp_y_max+(temp_y_max-temp_y_min)*0.22*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        
        
        tempTxt = sprintf('');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        xticks([1 2 3]);
        
        xtl = ["Baseline", "WM strength", "Interaction"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.155;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
        
        set(gca,'xticklabel','');
        
        
        set(gca,'linewidth',1.5)
        
        xlim([0.35 3.65]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.3
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        
        ylabel('Beta', 'FontSize', 9);
        temp_title = title(sprintf('Mismatch regression'),'fontsize',8);
        temp_title.Interpreter = 'none';
        subtitle(sprintf('low-strength mismatch\n & offload, r2 = %.3f',temp_r2),'fontsize',8);%9
        
    end
    
    %% Regression (two), using prior & WM strength to predict trial labels (CME or CF)
    if if_plot
        
        temp_chanceLevel = 0;
        
        fig = figure('Name','locDistri','NumberTitle','off');
        set(gcf,'Position',[550+0 50+350 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_p12 = p12_caseD2AB;
        
        temp_1 = r2_caseD2B_resampled;
        temp_2 = r2_caseD2A_resampled;
        
        
        temp_y_min = min([temp_1 temp_2]); %#ok<*NASGU>
        temp_y_max = max([temp_1 temp_2]);
        
        %temp_y_min = -0.1;
        %temp_y_max = 0.4790;
        
        [~,temp_y1_max] = bounds(temp_1);
        [~,temp_y2_max] = bounds(temp_2);
        
        
        if if_plot_violinplot0_pairline1 == 0
            temp_data = [temp_1';temp_2'];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            
            temp_label = [g1;g2];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 2, 1);
            
            violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'}]);
            
        elseif if_plot_violinplot0_pairline1 == 1
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            
            scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
        end
        
        
        plot([0 3],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.25,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        plot([1.0 2.0],temp_y_max+(temp_y_max-temp_y_min)*0.22*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        hold on
        
        
        xticks([1 2]);
        
        xtl = ["Baseline", "WM strength"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.155;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
        
        set(gca,'xticklabel','');
        
        
        set(gca,'linewidth',1.5)
        
        xlim([0.35 2.65]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.3
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        
        ylabel('Explained variance', 'FontSize', 9);
        temp_title = title(sprintf('Mismatch regression'),'fontsize',8);
        temp_title.Interpreter = 'none';
        %subtitle(sprintf('low-strength mismatch\n & offload'),'fontsize',8);%9
        
    end
    
    
    %% if_usePupilSize
    if if_usePupilSize == 1
        close all
        
        if_predict_prior0_decision1 = 0;%0
        
        
        total_ptbData;
        
        pupilSize_forAnalysis;
        pupilSizeBaseline;% filtered
        
        currentSessionIndex_AB;
        %temptemp_multiFOV_matrix_summary = nan(size(temp_multiFOV_matrix_summary));
        
        if if_monkey_D0_Z1 == 0
            temptemp_multiFOV_matrix_summary = ...
                [1,2;
                3,4;
                5,7;
                8,10;
                11,12;
                13,14;
                15,16;
                17,18;
                19,20;
                21,22;
                23,24;
                25,26];
            
        elseif if_monkey_D0_Z1 == 1
            temptemp_multiFOV_matrix_summary = reshape(1:32,[2,16])';
        end
        
        temptemp_FOV = temptemp_multiFOV_matrix_summary(currentSessionIndex_AB,:);
        
        temp_trial_num_multi = pupilSize_forAnalysis.trial_num_multi;
        
        temp_trial_num_multi_cumulative = nan(size(temp_trial_num_multi));
        for tempi=1:length(temp_trial_num_multi_cumulative)
            temp_trial_num_multi_cumulative(tempi) = sum(temp_trial_num_multi(1:tempi));
        end
        
        if currentSessionIndex_AB == 1
            temp_trialRange = [1,temp_trial_num_multi_cumulative(temptemp_FOV(2))];
        else
            temp_trialRange = ...
                [temp_trial_num_multi_cumulative(temptemp_FOV(1)-1)+1,...
                temp_trial_num_multi_cumulative(temptemp_FOV(2))];
        end
        temptemp_FOV;
        temp_trialRange;
        temp_choice_g_index_collapsed = (pupilSize_forAnalysis.choice_g_index_collapsed(temp_trialRange(1):temp_trialRange(2)))==1;
        temp_choice_r_index_collapsed = (pupilSize_forAnalysis.choice_r_index_collapsed(temp_trialRange(1):temp_trialRange(2)))==1;
        
        %total_ptbData
        temp_if_outlierTrial_collapsed = [];
        for file_index=temptemp_FOV(1):temptemp_FOV(2)
            ptbData = eval(['total_ptbData.file',sprintf('%d',file_index)]);
            temp_if_outlierTrial = ptbData.basic_para.if_outlierTrial;
            temp_if_outlierTrial_collapsed = [temp_if_outlierTrial_collapsed temp_if_outlierTrial]; %#ok<*AGROW>
        end
        temp_if_outlierTrial_collapsed = temp_if_outlierTrial_collapsed == 1;
        
        temp1 = choiceMemoryBoolIndex(~temp_if_outlierTrial_collapsed);
        temp2 = choiceOffloadBoolIndex(~temp_if_outlierTrial_collapsed);
        
        temp1_B = temp1 == temp_choice_g_index_collapsed;
        temp2_B = temp2 == temp_choice_r_index_collapsed;
        
        %% Extract pupil size
        pupilSize_forAnalysis = pupilSize_forAnalysis; %#ok<*ASGSL>
        %temp_eyeBaseline_pupilSize_cell_raw = pupilSize_forAnalysis.eyeBaseline_pupilSize_collapsed(temp_trialRange(1):temp_trialRange(2));
        %temp_eyeBaseline_pupilSize_cell_raw = pupilSize_forAnalysis.eyeDelay1_pupilSize_collapsed(temp_trialRange(1):temp_trialRange(2));
        temp_eyeBaseline_pupilSize_cell_raw = pupilSizeBaseline(temp_trialRange(1):temp_trialRange(2));
        
        %temp_size = 1100;
        
        %temp_eyeBaseline_pupilSize_raw = nan(length(temp_eyeBaseline_pupilSize_cell_raw),length(temp_eyeBaseline_pupilSize_cell_raw{1}));
        temp_eyeBaseline_pupilSize_raw = nan(length(temp_eyeBaseline_pupilSize_cell_raw),500);
        for tempi=1:size(temp_eyeBaseline_pupilSize_raw,1)
            %temp_eyeBaseline_pupilSize_raw(tempi,:) = temp_eyeBaseline_pupilSize_cell_raw{tempi};
            temp_eyeBaseline_pupilSize_raw(tempi,:) = temp_eyeBaseline_pupilSize_cell_raw{tempi}(601:1100);
        end
        temp_eyeBaselineMean_pupilSize_raw = mean(temp_eyeBaseline_pupilSize_raw,2);
        
        temp_eyeBaselineMean_pupilSize_C = temp_eyeBaselineMean_pupilSize_raw .* pupilSize_factor;%arbitray unit --> mm
        
        temp_eyeBaselineMean_pupilSize_C = ...
            (temp_eyeBaselineMean_pupilSize_C-mean(temp_eyeBaselineMean_pupilSize_C))./std(temp_eyeBaselineMean_pupilSize_C);
        
        temp_eyeBaselineMean_pupilSize = nan(length(choiceMemoryBoolIndex),1);
        temp_eyeBaselineMean_pupilSize(~temp_if_outlierTrial_collapsed) = temp_eyeBaselineMean_pupilSize_C;
        
        temp1 = mean(temp_eyeBaselineMean_pupilSize(choiceMemoryBoolIndex),'omitnan');
        temp2 = mean(temp_eyeBaselineMean_pupilSize(choiceOffloadBoolIndex),'omitnan');
        
        temptempIndex_A = find(choiceMemoryBoolIndex==true);
        temptempIndex_B = find(choiceOffloadBoolIndex==true);
        
        %temptempIndex_A(temptempIndex_A > 1035) = [];
        %temptempIndex_B(temptempIndex_B > 1035) = [];
        
        %temptempIndex_A(temptempIndex_A <= 1035) = [];
        %temptempIndex_B(temptempIndex_B <= 1035) = [];
        
        temp1 = temp_eyeBaselineMean_pupilSize(temptempIndex_A);
        temp2 = temp_eyeBaselineMean_pupilSize(temptempIndex_B);
        
        temp1_mean = mean(temp1,'omitnan');
        temp2_mean = mean(temp2,'omitnan');
        fprintf('pupilSize_baseline(memory, offload) = [%.3f,%.3f].\n',temp1_mean,temp2_mean);
        
        temp_eyeBaselineMean_pupilSize;
        
        
        %%% Further computation
        %% Regression: Prior ~ Pupil size + trial history
        meta_trialLevel_baseline;
        temp_eyeBaselineMean_pupilSize;
        historyWeightedReward;
        
        %temp1 = mean(historyWeightedReward(choiceMemoryBoolIndex),'omitnan');
        %temp2 = mean(historyWeightedReward(choiceOffloadBoolIndex),'omitnan');
        
        tempNONNANBoolIndex123;
        
        tempNONNANBoolIndex_new = tempNONNANBoolIndex123 & (~isnan(temp_eyeBaselineMean_pupilSize)) & (~isnan(historyWeightedReward));
        
        if if_predict_prior0_decision1 == 1
            temptempBoolIndex_y = nan(length(choiceBoolIndex),1);
            temptempBoolIndex_y(choiceMemoryBoolIndex) = true;
            temptempBoolIndex_y(choiceOffloadBoolIndex) = false;
            temptempBoolIndex_y(~tempNONNANBoolIndex_new) = nan;
            
            tempNONNANBoolIndex_new = tempNONNANBoolIndex_new & (~isnan(temptempBoolIndex_y));
        end
        
        tempNONNAN_new_eyeBaselineMean_pupilSize = temp_eyeBaselineMean_pupilSize(tempNONNANBoolIndex_new);
        %tempNONNAN_new_eyeBaselineMean_pupilSize = temp_eyeBaselineMean_pupilSize(tempNONNANBoolIndex_new) .* -1;
        tempNONNAN_new_historyWeightedReward = historyWeightedReward(tempNONNANBoolIndex_new);
        
        if if_predict_prior0_decision1 == 0
            tempNONNAN_new_meta_trialLevel_baseline = meta_trialLevel_baseline(tempNONNANBoolIndex_new);
        elseif if_predict_prior0_decision1 == 1
            tempNONNAN_new_meta_trialLevel_baseline = temptempBoolIndex_y(tempNONNANBoolIndex_new);
        end
        
        
        temp_trialNum = length(tempNONNAN_new_meta_trialLevel_baseline);
        
        
        r2_caseZ_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseZ_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseZ_resampled = nan(1,resampleIterCount_GLM);
        beta2_caseZ_resampled = nan(1,resampleIterCount_GLM);
        beta3_caseZ_resampled = nan(1,resampleIterCount_GLM);
        
        
        %for iter=1:resampleIterCount_GLM
        parfor iter=1:resampleIterCount_GLM
            temp_trialIndex_resample = sort(randi(temp_trialNum,1,temp_trialNum));
            
            tempNONNAN_eyeBaselineMean_pupilSize_resampled = tempNONNAN_new_eyeBaselineMean_pupilSize(temp_trialIndex_resample); %#ok<*PFBNS>
            tempNONNAN_historyWeightedReward_resampled = tempNONNAN_new_historyWeightedReward(temp_trialIndex_resample);
            tempNONNAN_meta_trialLevel_baseline_resampled = tempNONNAN_new_meta_trialLevel_baseline(temp_trialIndex_resample);
            
            temp_x1 = tempNONNAN_eyeBaselineMean_pupilSize_resampled;
            temp_x2 = tempNONNAN_historyWeightedReward_resampled;
            temp_y = tempNONNAN_meta_trialLevel_baseline_resampled;
            
            temp_x1_zscored = (temp_x1-mean(temp_x1))./std(temp_x1);
            temp_x2_zscored = (temp_x2-mean(temp_x2))./std(temp_x2);
            temp_y_zscored = (temp_y-mean(temp_y))./std(temp_y);
            
            
            x = [temp_x1_zscored,temp_x2_zscored];
            y = temp_y_zscored;
            
            temp_mdl = fitglm(x,y,'interactions');
            r2_caseZ_resampled(iter) = temp_mdl.Rsquared.Adjusted;
            beta0_caseZ_resampled(iter) = temp_mdl.Coefficients.Estimate(1);
            beta1_caseZ_resampled(iter) = temp_mdl.Coefficients.Estimate(2);
            beta2_caseZ_resampled(iter) = temp_mdl.Coefficients.Estimate(3);
            beta3_caseZ_resampled(iter) = temp_mdl.Coefficients.Estimate(4);
        end
        
        [~,p1_caseZ,~,~] = ttest(beta1_caseZ_resampled,0);
        [~,p2_caseZ,~,~] = ttest(beta2_caseZ_resampled,0);
        [~,p3_caseZ,~,~] = ttest(beta3_caseZ_resampled,0);
        
        
        %% Plot. Regression: Prior ~ Pupil size + trial history
        if if_plot == 1
            
            temp_chanceLevel = 0;
            
            fig = figure('Name','locDistri','NumberTitle','off');
            %set(gcf,'Position',[1150+0 400+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05*1.05*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            set(gcf,'Position',[1150+0 400+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05*1.05*1.05*0.8 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
            nexttile
            
            temp_p1 = p1_caseZ;
            temp_p2 = p2_caseZ;
            temp_p3 = p3_caseZ;
            
            temp_1 = beta1_caseZ_resampled;
            temp_2 = beta2_caseZ_resampled;
            temp_3 = beta3_caseZ_resampled;
            
            temp_r2 = mean(r2_caseZ_resampled);
            fprintf('r2 of prior regression = %.3f.\n',temp_r2);
            
            temp_y_min = min([temp_1 temp_2 temp_3]); %#ok<*NASGU>
            temp_y_max = max([temp_1 temp_2 temp_3]);
            
            %temp_y_min = -0.1;
            %temp_y_max = 0.4790;
            
            %temp_y_min = -0.2+0.02*2+0.01-0.00115;
            %temp_y_max = 0.25-0.02*4+(0.005/3);
            
            temp_y_min = -0.18;
            temp_y_max = 0.21+0.007;
            
            
            [~,temp_y1_max] = bounds(temp_1);
            [~,temp_y2_max] = bounds(temp_2);
            [~,temp_y3_max] = bounds(temp_3);
            
            %temp_if_plot_violinplot0_pairline1 = if_plot_violinplot0_pairline1;
            temp_if_plot_violinplot0_pairline1 = 0;
            
            if temp_if_plot_violinplot0_pairline1 == 0
                temp_data = [temp_1';temp_2';temp_3'];
                
                g1 = repmat({'A'},length(temp_1),1);
                g2 = repmat({'B'},length(temp_2),1);
                g3 = repmat({'C'},length(temp_3),1);
                
                temp_label = [g1;g2;g3];
                
                temptemp_color1 = [1 1 1]*0.5;
                temptemp_color2 = repmat(temptemp_color1, 3, 1);
                
                violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                    'GroupOrder',[{'A'};{'B'};{'C'}]);
                
            elseif temp_if_plot_violinplot0_pairline1 == 1
                %for tempi=1:length(temp_1)
                %    plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                %    hold on
                %end
                
                plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                
                
                scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
                scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
                scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
            end
            
            
            plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
            hold on
            
            
            tempTxt = sprintf('');
            if temp_p1 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p1 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p1 < 0.05
                tempTxt = sprintf('*');
            end
            text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p2 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p2 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p2 < 0.05
                tempTxt = sprintf('*');
            end
            text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p3 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p3 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p3 < 0.05
                tempTxt = sprintf('*');
            end
            text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            xticks([1 2 3]);
            
            xtl = ["Pupil size", "History reward", "Interaction"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
            xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.255;
            text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9,12,8,7.5
            
            set(gca,'xticklabel','');
            
            
            set(gca,'linewidth',1.5)
            
            xlim([0.35 3.65]);
            %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.2
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.2 temp_y_max+(temp_y_max-temp_y_min)*0.2]);%0.2
            set(gca, 'FontSize', 8) %14
            set(gca,'box','off');% 取消右、上边框
            
            
            
            ylabel('Beta', 'FontSize', 9);
            if if_predict_prior0_decision1 == 0
                temp_title = title(sprintf('Prior regression'),'fontsize',8);
            elseif if_predict_prior0_decision1 == 1
                temp_title = title(sprintf('Decision regression'),'fontsize',8);
            end
            temp_title.Interpreter = 'none';
            %subtitle(sprintf('r2 = %.3f',temp_r2),'fontsize',8);%9
            subtitle(sprintf('all trials, r2 = %.3f',temp_r2),'fontsize',8);%9
            
        end
        
        
        
        %% Regression: Pupil size ~ WM strength + meta-WM
        temp_eyeBaselineMean_pupilSize;
        meta_trialLevel;
        memoryPrecision_trialLevel;
        
        
        tempNONNANBoolIndex_new = tempNONNANBoolIndex123 & (~isnan(temp_eyeBaselineMean_pupilSize));
        
        
        tempNONNAN_new_eyeBaselineMean_pupilSize = temp_eyeBaselineMean_pupilSize(tempNONNANBoolIndex_new);
        tempNONNAN_new_memoryPrecision_trialLevel = memoryPrecision_trialLevel(tempNONNANBoolIndex_new);
        tempNONNAN_new_meta_trialLevel = meta_trialLevel(tempNONNANBoolIndex_new);
        
        
        temp_trialNum = length(tempNONNAN_new_meta_trialLevel);
        
        
        r2_caseZ2_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseZ2_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseZ2_resampled = nan(1,resampleIterCount_GLM);
        beta2_caseZ2_resampled = nan(1,resampleIterCount_GLM);
        beta3_caseZ2_resampled = nan(1,resampleIterCount_GLM);
        
        
        %for iter=1:resampleIterCount_GLM
        parfor iter=1:resampleIterCount_GLM
            temp_trialIndex_resample = sort(randi(temp_trialNum,1,temp_trialNum));
            
            tempNONNAN_eyeBaselineMean_pupilSize_resampled = tempNONNAN_new_eyeBaselineMean_pupilSize(temp_trialIndex_resample); %#ok<*PFBNS>
            tempNONNAN_memoryPrecision_trialLevel_baseline = tempNONNAN_new_memoryPrecision_trialLevel(temp_trialIndex_resample);
            tempNONNAN_meta_trialLevel_resampled = tempNONNAN_new_meta_trialLevel(temp_trialIndex_resample);
            
            temp_x1 = tempNONNAN_memoryPrecision_trialLevel_baseline;
            temp_x2 = tempNONNAN_meta_trialLevel_resampled;
            temp_y = tempNONNAN_eyeBaselineMean_pupilSize_resampled;
            
            temp_x1_zscored = (temp_x1-mean(temp_x1))./std(temp_x1);
            temp_x2_zscored = (temp_x2-mean(temp_x2))./std(temp_x2);
            temp_y_zscored = (temp_y-mean(temp_y))./std(temp_y);
            
            
            x = [temp_x1_zscored,temp_x2_zscored];
            y = temp_y_zscored;
            
            temp_mdl = fitglm(x,y,'interactions');
            r2_caseZ2_resampled(iter) = temp_mdl.Rsquared.Adjusted;
            beta0_caseZ2_resampled(iter) = temp_mdl.Coefficients.Estimate(1);
            beta1_caseZ2_resampled(iter) = temp_mdl.Coefficients.Estimate(2);
            beta2_caseZ2_resampled(iter) = temp_mdl.Coefficients.Estimate(3);
            beta3_caseZ2_resampled(iter) = temp_mdl.Coefficients.Estimate(4);
        end
        
        [~,p1_caseZ2,~,~] = ttest(beta1_caseZ2_resampled,0);
        [~,p2_caseZ2,~,~] = ttest(beta2_caseZ2_resampled,0);
        [~,p3_caseZ2,~,~] = ttest(beta3_caseZ2_resampled,0);
        
        [~,p12_caseZ2,~,~] = ttest(beta1_caseZ2_resampled,beta2_caseZ2_resampled);
        
        
        %% Plot. Regression: Pupil size ~ WM strength + meta-WM
        if if_plot == 1
            
            temp_chanceLevel = 0;
            
            fig = figure('Name','locDistri','NumberTitle','off');
            set(gcf,'Position',[1150+0 400+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05*1.05*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
            nexttile
            
            temp_p1 = p1_caseZ2;
            temp_p2 = p2_caseZ2;
            temp_p3 = p3_caseZ2;
            temp_p12 = p12_caseZ2;
            
            
            temp_1 = beta1_caseZ2_resampled;
            temp_2 = beta2_caseZ2_resampled;
            temp_3 = beta3_caseZ2_resampled;
            
            temp_r2 = mean(r2_caseZ2_resampled);
            
            
            temp_y_min = min([temp_1 temp_2 temp_3]); %#ok<*NASGU>
            temp_y_max = max([temp_1 temp_2 temp_3]);
            
            %temp_y_min = -0.1;
            %temp_y_max = 0.4790;
            
            [~,temp_y1_max] = bounds(temp_1);
            [~,temp_y2_max] = bounds(temp_2);
            [~,temp_y3_max] = bounds(temp_3);
            
            %temp_if_plot_violinplot0_pairline1 = if_plot_violinplot0_pairline1;
            temp_if_plot_violinplot0_pairline1 = 0;
            
            if temp_if_plot_violinplot0_pairline1 == 0
                temp_data = [temp_1';temp_2';temp_3'];
                
                g1 = repmat({'A'},length(temp_1),1);
                g2 = repmat({'B'},length(temp_2),1);
                g3 = repmat({'C'},length(temp_3),1);
                
                temp_label = [g1;g2;g3];
                
                temptemp_color1 = [1 1 1]*0.5;
                temptemp_color2 = repmat(temptemp_color1, 3, 1);
                
                violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                    'GroupOrder',[{'A'};{'B'};{'C'}]);
                
            elseif temp_if_plot_violinplot0_pairline1 == 1
                %for tempi=1:length(temp_1)
                %    plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                %    hold on
                %end
                
                plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                
                
                scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
                scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
                scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
            end
            
            
            plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
            hold on
            
            
            tempTxt = sprintf('');
            if temp_p1 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p1 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p1 < 0.05
                tempTxt = sprintf('*');
            end
            text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p2 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p2 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p2 < 0.05
                tempTxt = sprintf('*');
            end
            text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p3 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p3 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p3 < 0.05
                tempTxt = sprintf('*');
            end
            text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p12 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p12 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p12 < 0.05
                tempTxt = sprintf('*');
            end
            text(1.5,temp_y_min-(temp_y_max-temp_y_min)*0.21,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            plot([1 2],temp_y_min-(temp_y_max-temp_y_min)*0.1*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
            hold on
            
            
            xticks([1 2 3]);
            
            xtl = ["WM strength", "Meta-WM", "Interaction"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
            xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.37;
            text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
            
            set(gca,'xticklabel','');
            
            
            set(gca,'linewidth',1.5)
            
            xlim([0.35 3.65]);
            %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.2
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.3 temp_y_max+(temp_y_max-temp_y_min)*0.3]);%0.2
            set(gca, 'FontSize', 8) %14
            set(gca,'box','off');% 取消右、上边框
            
            
            
            ylabel('Beta', 'FontSize', 9);
            temp_title = title(sprintf('Pupil size regression'),'fontsize',8);
            temp_title.Interpreter = 'none';
            %subtitle(sprintf('r2 = %.3f',temp_r2),'fontsize',8);%9
            %subtitle(sprintf('all trials, r2 = %.3f',temp_r2),'fontsize',8);%9
            subtitle(sprintf('r2 = %.3f',temp_r2),'fontsize',8);%9
            
        end
        
        
        
        %% Regression: Meta-WM ~ Baseline meta + pupil size + trial history
        meta_trialLevel;
        meta_trialLevel_baseline;
        temp_eyeBaselineMean_pupilSize;
        historyWeightedReward;
        
        
        
        tempNONNANBoolIndex_new = tempNONNANBoolIndex123 & (~isnan(temp_eyeBaselineMean_pupilSize)) & (~isnan(historyWeightedReward));
        
        
        tempNONNAN_new_meta_trialLevel = meta_trialLevel(tempNONNANBoolIndex_new);
        tempNONNAN_new_meta_trialLevel_baseline = meta_trialLevel_baseline(tempNONNANBoolIndex_new);
        tempNONNAN_new_eyeBaselineMean_pupilSize = temp_eyeBaselineMean_pupilSize(tempNONNANBoolIndex_new);
        tempNONNAN_new_historyWeightedReward = historyWeightedReward(tempNONNANBoolIndex_new);
        
        
        temp_trialNum = length(tempNONNAN_new_meta_trialLevel);
        
        
        r2_caseZ3_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseZ3_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseZ3_resampled = nan(1,resampleIterCount_GLM);
        beta2_caseZ3_resampled = nan(1,resampleIterCount_GLM);
        beta3_caseZ3_resampled = nan(1,resampleIterCount_GLM);
        
        
        %for iter=1:resampleIterCount_GLM
        parfor iter=1:resampleIterCount_GLM
            temp_trialIndex_resample = sort(randi(temp_trialNum,1,temp_trialNum));
            
            tempNONNAN_meta_trialLevel_resampled = tempNONNAN_new_meta_trialLevel(temp_trialIndex_resample);
            tempNONNAN_meta_trialLevel_baseline_resampled = tempNONNAN_new_meta_trialLevel_baseline(temp_trialIndex_resample);
            tempNONNAN_eyeBaselineMean_pupilSize_resampled = tempNONNAN_new_eyeBaselineMean_pupilSize(temp_trialIndex_resample); %#ok<*PFBNS>
            tempNONNAN_historyWeightedReward_resampled = tempNONNAN_new_historyWeightedReward(temp_trialIndex_resample); %#ok<*PFBNS>
            
            temp_x1 = tempNONNAN_meta_trialLevel_baseline_resampled;
            temp_x2 = tempNONNAN_eyeBaselineMean_pupilSize_resampled;
            temp_x3 = tempNONNAN_historyWeightedReward_resampled;
            temp_y = tempNONNAN_meta_trialLevel_resampled;
            
            temp_x1_zscored = (temp_x1-mean(temp_x1))./std(temp_x1);
            temp_x2_zscored = (temp_x2-mean(temp_x2))./std(temp_x2);
            temp_x3_zscored = (temp_x3-mean(temp_x3))./std(temp_x3);
            temp_y_zscored = (temp_y-mean(temp_y))./std(temp_y);
            
            
            x = [temp_x1_zscored,temp_x2_zscored,temp_x3_zscored];
            y = temp_y_zscored;
            
            %temp_mdl = fitglm(x,y,'interactions');
            temp_mdl = fitglm(x,y,'linear');
            r2_caseZ3_resampled(iter) = temp_mdl.Rsquared.Adjusted;
            beta0_caseZ3_resampled(iter) = temp_mdl.Coefficients.Estimate(1);
            beta1_caseZ3_resampled(iter) = temp_mdl.Coefficients.Estimate(2);
            beta2_caseZ3_resampled(iter) = temp_mdl.Coefficients.Estimate(3);
            beta3_caseZ3_resampled(iter) = temp_mdl.Coefficients.Estimate(4);
        end
        
        [~,p1_caseZ3,~,~] = ttest(beta1_caseZ3_resampled,0);
        [~,p2_caseZ3,~,~] = ttest(beta2_caseZ3_resampled,0);
        [~,p3_caseZ3,~,~] = ttest(beta3_caseZ3_resampled,0);
        
        %% Plot. Regression: Meta-WM ~ Baseline meta + pupil size + trial history
        if if_plot == 1
            
            temp_chanceLevel = 0;
            
            fig = figure('Name','locDistri','NumberTitle','off');
            set(gcf,'Position',[1150+0 400+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05*1.05*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
            nexttile
            
            temp_p1 = p1_caseZ3;
            temp_p2 = p2_caseZ3;
            temp_p3 = p3_caseZ3;
            
            
            temp_1 = beta1_caseZ3_resampled;
            temp_2 = beta2_caseZ3_resampled;
            temp_3 = beta3_caseZ3_resampled;
            
            temp_r2 = mean(r2_caseZ3_resampled);
            
            
            temp_y_min = min([temp_1 temp_2 temp_3]); %#ok<*NASGU>
            temp_y_max = max([temp_1 temp_2 temp_3]);
            
            %temp_y_min = -0.1;
            %temp_y_max = 0.4790;
            
            [~,temp_y1_max] = bounds(temp_1);
            [~,temp_y2_max] = bounds(temp_2);
            [~,temp_y3_max] = bounds(temp_3);
            
            %temp_if_plot_violinplot0_pairline1 = if_plot_violinplot0_pairline1;
            temp_if_plot_violinplot0_pairline1 = 0;
            
            if temp_if_plot_violinplot0_pairline1 == 0
                temp_data = [temp_1';temp_2';temp_3'];
                
                g1 = repmat({'A'},length(temp_1),1);
                g2 = repmat({'B'},length(temp_2),1);
                g3 = repmat({'C'},length(temp_3),1);
                
                temp_label = [g1;g2;g3];
                
                temptemp_color1 = [1 1 1]*0.5;
                temptemp_color2 = repmat(temptemp_color1, 3, 1);
                
                violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                    'GroupOrder',[{'A'};{'B'};{'C'}]);
                
            elseif temp_if_plot_violinplot0_pairline1 == 1
                %for tempi=1:length(temp_1)
                %    plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                %    hold on
                %end
                
                plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                
                
                scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
                scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
                scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
            end
            
            
            plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
            hold on
            
            
            tempTxt = sprintf('');
            if temp_p1 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p1 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p1 < 0.05
                tempTxt = sprintf('*');
            end
            text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p2 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p2 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p2 < 0.05
                tempTxt = sprintf('*');
            end
            text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p3 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p3 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p3 < 0.05
                tempTxt = sprintf('*');
            end
            text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            
            xticks([1 2 3]);
            
            xtl = ["Baseline", "Pupil size", "History reward"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
            xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.255;
            text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
            
            set(gca,'xticklabel','');
            
            
            set(gca,'linewidth',1.5)
            
            xlim([0.35 3.65]);
            %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.2
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.2 temp_y_max+(temp_y_max-temp_y_min)*0.2]);%0.2
            set(gca, 'FontSize', 8) %14
            set(gca,'box','off');% 取消右、上边框
            
            
            
            ylabel('Beta', 'FontSize', 9);
            temp_title = title(sprintf('Meta-WM regression'),'fontsize',8);
            temp_title.Interpreter = 'none';
            %subtitle(sprintf('r2 = %.3f',temp_r2),'fontsize',8);%9
            %subtitle(sprintf('all trials, r2 = %.3f',temp_r2),'fontsize',8);%9
            subtitle(sprintf('r2 = %.3f',temp_r2),'fontsize',8);%9
            
        end
        
        
        %% Regression: Pupil size ~ Baseline meta + Meta-WM
        meta_trialLevel;
        meta_trialLevel_baseline;
        temp_eyeBaselineMean_pupilSize;
        
        
        
        tempNONNANBoolIndex_new = tempNONNANBoolIndex123 & (~isnan(temp_eyeBaselineMean_pupilSize));
        
        
        tempNONNAN_new_meta_trialLevel = meta_trialLevel(tempNONNANBoolIndex_new);
        tempNONNAN_new_meta_trialLevel_baseline = meta_trialLevel_baseline(tempNONNANBoolIndex_new);
        tempNONNAN_new_eyeBaselineMean_pupilSize = temp_eyeBaselineMean_pupilSize(tempNONNANBoolIndex_new);
        
        
        temp_trialNum = length(tempNONNAN_new_meta_trialLevel);
        
        
        r2_caseZ4_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseZ4_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseZ4_resampled = nan(1,resampleIterCount_GLM);
        beta2_caseZ4_resampled = nan(1,resampleIterCount_GLM);
        beta3_caseZ4_resampled = nan(1,resampleIterCount_GLM);
        
        
        %for iter=1:resampleIterCount_GLM
        parfor iter=1:resampleIterCount_GLM
            temp_trialIndex_resample = sort(randi(temp_trialNum,1,temp_trialNum));
            
            tempNONNAN_meta_trialLevel_resampled = tempNONNAN_new_meta_trialLevel(temp_trialIndex_resample);
            tempNONNAN_meta_trialLevel_baseline_resampled = tempNONNAN_new_meta_trialLevel_baseline(temp_trialIndex_resample);
            tempNONNAN_eyeBaselineMean_pupilSize_resampled = tempNONNAN_new_eyeBaselineMean_pupilSize(temp_trialIndex_resample); %#ok<*PFBNS>
            
            temp_x1 = tempNONNAN_meta_trialLevel_baseline_resampled;
            temp_x2 = tempNONNAN_meta_trialLevel_resampled;
            temp_y = tempNONNAN_eyeBaselineMean_pupilSize_resampled;
            
            temp_x1_zscored = (temp_x1-mean(temp_x1))./std(temp_x1);
            temp_x2_zscored = (temp_x2-mean(temp_x2))./std(temp_x2);
            temp_y_zscored = (temp_y-mean(temp_y))./std(temp_y);
            
            
            x = [temp_x1_zscored,temp_x2_zscored];
            y = temp_y_zscored;
            
            temp_mdl = fitglm(x,y,'interactions');
            %temp_mdl = fitglm(x,y,'linear');
            r2_caseZ4_resampled(iter) = temp_mdl.Rsquared.Adjusted;
            beta0_caseZ4_resampled(iter) = temp_mdl.Coefficients.Estimate(1);
            beta1_caseZ4_resampled(iter) = temp_mdl.Coefficients.Estimate(2);
            beta2_caseZ4_resampled(iter) = temp_mdl.Coefficients.Estimate(3);
            beta3_caseZ4_resampled(iter) = temp_mdl.Coefficients.Estimate(4);
        end
        
        [~,p1_caseZ4,~,~] = ttest(beta1_caseZ4_resampled,0);
        [~,p2_caseZ4,~,~] = ttest(beta2_caseZ4_resampled,0);
        [~,p3_caseZ4,~,~] = ttest(beta3_caseZ4_resampled,0);
        
        [~,p12_caseZ4,~,~] = ttest(beta1_caseZ4_resampled,beta2_caseZ4_resampled);
        
        %% Plot. Regression: Pupil size ~ Baseline meta + Meta-WM
        if if_plot == 1
            
            temp_chanceLevel = 0;
            
            fig = figure('Name','locDistri','NumberTitle','off');
            set(gcf,'Position',[1150+0 400+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05*1.05*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
            nexttile
            
            temp_p1 = p1_caseZ4;
            temp_p2 = p2_caseZ4;
            temp_p3 = p3_caseZ4;
            
            
            temp_1 = beta1_caseZ4_resampled;
            temp_2 = beta2_caseZ4_resampled;
            temp_3 = beta3_caseZ4_resampled;
            
            temp_r2 = mean(r2_caseZ4_resampled);
            
            
            temp_y_min = min([temp_1 temp_2 temp_3]); %#ok<*NASGU>
            temp_y_max = max([temp_1 temp_2 temp_3]);
            
            %temp_y_min = -0.1;
            %temp_y_max = 0.4790;
            
            [~,temp_y1_max] = bounds(temp_1);
            [~,temp_y2_max] = bounds(temp_2);
            [~,temp_y3_max] = bounds(temp_3);
            
            %temp_if_plot_violinplot0_pairline1 = if_plot_violinplot0_pairline1;
            temp_if_plot_violinplot0_pairline1 = 0;
            
            if temp_if_plot_violinplot0_pairline1 == 0
                temp_data = [temp_1';temp_2';temp_3'];
                
                g1 = repmat({'A'},length(temp_1),1);
                g2 = repmat({'B'},length(temp_2),1);
                g3 = repmat({'C'},length(temp_3),1);
                
                temp_label = [g1;g2;g3];
                
                temptemp_color1 = [1 1 1]*0.5;
                temptemp_color2 = repmat(temptemp_color1, 3, 1);
                
                violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                    'GroupOrder',[{'A'};{'B'};{'C'}]);
                
            elseif temp_if_plot_violinplot0_pairline1 == 1
                %for tempi=1:length(temp_1)
                %    plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                %    hold on
                %end
                
                plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                
                
                scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
                scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
                scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
            end
            
            
            plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
            hold on
            
            
            tempTxt = sprintf('');
            if temp_p1 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p1 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p1 < 0.05
                tempTxt = sprintf('*');
            end
            text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p2 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p2 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p2 < 0.05
                tempTxt = sprintf('*');
            end
            text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p3 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p3 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p3 < 0.05
                tempTxt = sprintf('*');
            end
            text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p12 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p12 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p12 < 0.05
                tempTxt = sprintf('*');
            end
            text(1.5,temp_y_min-(temp_y_max-temp_y_min)*0.21,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            plot([1 2],temp_y_min-(temp_y_max-temp_y_min)*0.1*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
            hold on
            
            
            xticks([1 2 3]);
            
            xtl = ["Baseline", "Meta", "Interaction"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
            xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.37;
            text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
            
            set(gca,'xticklabel','');
            
            
            set(gca,'linewidth',1.5)
            
            xlim([0.35 3.65]);
            %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.2
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.3 temp_y_max+(temp_y_max-temp_y_min)*0.3]);%0.2
            set(gca, 'FontSize', 8) %14
            set(gca,'box','off');% 取消右、上边框
            
            
            
            ylabel('Beta', 'FontSize', 9);
            temp_title = title(sprintf('Pupil size regression'),'fontsize',8);
            temp_title.Interpreter = 'none';
            %subtitle(sprintf('r2 = %.3f',temp_r2),'fontsize',8);%9
            %subtitle(sprintf('all trials, r2 = %.3f',temp_r2),'fontsize',8);%9
            subtitle(sprintf('r2 = %.3f',temp_r2),'fontsize',8);%9
            
        end
        
        
        %% Regression: Trial history ~ Baseline meta + Meta-WM
        meta_trialLevel;
        meta_trialLevel_baseline;
        historyWeightedReward;
        
        
        
        tempNONNANBoolIndex_new = tempNONNANBoolIndex123 & (~isnan(historyWeightedReward));
        
        
        tempNONNAN_new_meta_trialLevel = meta_trialLevel(tempNONNANBoolIndex_new);
        tempNONNAN_new_meta_trialLevel_baseline = meta_trialLevel_baseline(tempNONNANBoolIndex_new);
        tempNONNAN_new_historyWeightedReward = historyWeightedReward(tempNONNANBoolIndex_new);
        
        
        temp_trialNum = length(tempNONNAN_new_meta_trialLevel);
        
        
        r2_caseZ5_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseZ5_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseZ5_resampled = nan(1,resampleIterCount_GLM);
        beta2_caseZ5_resampled = nan(1,resampleIterCount_GLM);
        beta3_caseZ5_resampled = nan(1,resampleIterCount_GLM);
        
        
        %for iter=1:resampleIterCount_GLM
        parfor iter=1:resampleIterCount_GLM
            temp_trialIndex_resample = sort(randi(temp_trialNum,1,temp_trialNum));
            
            tempNONNAN_meta_trialLevel_resampled = tempNONNAN_new_meta_trialLevel(temp_trialIndex_resample);
            tempNONNAN_meta_trialLevel_baseline_resampled = tempNONNAN_new_meta_trialLevel_baseline(temp_trialIndex_resample);
            tempNONNAN_historyWeightedReward_resampled = tempNONNAN_new_historyWeightedReward(temp_trialIndex_resample); %#ok<*PFBNS>
            
            temp_x1 = tempNONNAN_meta_trialLevel_baseline_resampled;
            temp_x2 = tempNONNAN_meta_trialLevel_resampled;
            temp_y = tempNONNAN_historyWeightedReward_resampled;
            
            temp_x1_zscored = (temp_x1-mean(temp_x1))./std(temp_x1);
            temp_x2_zscored = (temp_x2-mean(temp_x2))./std(temp_x2);
            temp_y_zscored = (temp_y-mean(temp_y))./std(temp_y);
            
            
            x = [temp_x1_zscored,temp_x2_zscored];
            y = temp_y_zscored;
            
            temp_mdl = fitglm(x,y,'interactions');
            %temp_mdl = fitglm(x,y,'linear');
            r2_caseZ5_resampled(iter) = temp_mdl.Rsquared.Adjusted;
            beta0_caseZ5_resampled(iter) = temp_mdl.Coefficients.Estimate(1);
            beta1_caseZ5_resampled(iter) = temp_mdl.Coefficients.Estimate(2);
            beta2_caseZ5_resampled(iter) = temp_mdl.Coefficients.Estimate(3);
            beta3_caseZ5_resampled(iter) = temp_mdl.Coefficients.Estimate(4);
        end
        
        [~,p1_caseZ5,~,~] = ttest(beta1_caseZ5_resampled,0);
        [~,p2_caseZ5,~,~] = ttest(beta2_caseZ5_resampled,0);
        [~,p3_caseZ5,~,~] = ttest(beta3_caseZ5_resampled,0);
        
        [~,p12_caseZ5,~,~] = ttest(beta1_caseZ5_resampled,beta2_caseZ5_resampled);
        
        
        %% Plot. Regression: Pupil size ~ Baseline meta + Meta-WM
        if if_plot == 1
            
            temp_chanceLevel = 0;
            
            fig = figure('Name','locDistri','NumberTitle','off');
            set(gcf,'Position',[1150+0 400+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05*1.05*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
            nexttile
            
            temp_p1 = p1_caseZ5;
            temp_p2 = p2_caseZ5;
            temp_p3 = p3_caseZ5;
            
            
            temp_1 = beta1_caseZ5_resampled;
            temp_2 = beta2_caseZ5_resampled;
            temp_3 = beta3_caseZ5_resampled;
            
            temp_r2 = mean(r2_caseZ5_resampled);
            
            
            temp_y_min = min([temp_1 temp_2 temp_3]); %#ok<*NASGU>
            temp_y_max = max([temp_1 temp_2 temp_3]);
            
            %temp_y_min = -0.1;
            %temp_y_max = 0.4790;
            
            [~,temp_y1_max] = bounds(temp_1);
            [~,temp_y2_max] = bounds(temp_2);
            [~,temp_y3_max] = bounds(temp_3);
            
            %temp_if_plot_violinplot0_pairline1 = if_plot_violinplot0_pairline1;
            temp_if_plot_violinplot0_pairline1 = 0;
            
            if temp_if_plot_violinplot0_pairline1 == 0
                temp_data = [temp_1';temp_2';temp_3'];
                
                g1 = repmat({'A'},length(temp_1),1);
                g2 = repmat({'B'},length(temp_2),1);
                g3 = repmat({'C'},length(temp_3),1);
                
                temp_label = [g1;g2;g3];
                
                temptemp_color1 = [1 1 1]*0.5;
                temptemp_color2 = repmat(temptemp_color1, 3, 1);
                
                violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                    'GroupOrder',[{'A'};{'B'};{'C'}]);
                
            elseif temp_if_plot_violinplot0_pairline1 == 1
                %for tempi=1:length(temp_1)
                %    plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                %    hold on
                %end
                
                plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                
                
                scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
                scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
                scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
            end
            
            
            plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
            hold on
            
            
            tempTxt = sprintf('');
            if temp_p1 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p1 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p1 < 0.05
                tempTxt = sprintf('*');
            end
            text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p2 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p2 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p2 < 0.05
                tempTxt = sprintf('*');
            end
            text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p3 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p3 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p3 < 0.05
                tempTxt = sprintf('*');
            end
            text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p12 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p12 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p12 < 0.05
                tempTxt = sprintf('*');
            end
            text(1.5,temp_y_min-(temp_y_max-temp_y_min)*0.21,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            plot([1 2],temp_y_min-(temp_y_max-temp_y_min)*0.1*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
            hold on
            
            
            xticks([1 2 3]);
            
            xtl = ["Baseline", "Meta", "Interaction"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
            xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.37;
            text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
            
            set(gca,'xticklabel','');
            
            
            set(gca,'linewidth',1.5)
            
            xlim([0.35 3.65]);
            %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.2
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.3 temp_y_max+(temp_y_max-temp_y_min)*0.3]);%0.2
            set(gca, 'FontSize', 8) %14
            set(gca,'box','off');% 取消右、上边框
            
            
            
            ylabel('Beta', 'FontSize', 9);
            temp_title = title(sprintf('History reward regression'),'fontsize',8);
            temp_title.Interpreter = 'none';
            %subtitle(sprintf('r2 = %.3f',temp_r2),'fontsize',8);%9
            %subtitle(sprintf('all trials, r2 = %.3f',temp_r2),'fontsize',8);%9
            subtitle(sprintf('r2 = %.3f',temp_r2),'fontsize',8);%9
            
        end
        
        
        %% Regression: Baseline meta ~ WM strength + Meta-WM
        meta_trialLevel;
        memoryPrecision_trialLevel;
        meta_trialLevel_baseline;
        
        
        
        
        tempNONNANBoolIndex_new = tempNONNANBoolIndex123;
        
        
        tempNONNAN_new_meta_trialLevel = meta_trialLevel(tempNONNANBoolIndex_new);
        tempNONNAN_new_memoryPrecision_trialLevel = memoryPrecision_trialLevel(tempNONNANBoolIndex_new);
        tempNONNAN_new_meta_trialLevel_baseline = meta_trialLevel_baseline(tempNONNANBoolIndex_new);
        
        
        temp_trialNum = length(tempNONNAN_new_meta_trialLevel);
        
        
        r2_caseZ6_resampled = nan(1,resampleIterCount_GLM);
        beta0_caseZ6_resampled = nan(1,resampleIterCount_GLM);
        beta1_caseZ6_resampled = nan(1,resampleIterCount_GLM);
        beta2_caseZ6_resampled = nan(1,resampleIterCount_GLM);
        beta3_caseZ6_resampled = nan(1,resampleIterCount_GLM);
        
        
        %for iter=1:resampleIterCount_GLM
        parfor iter=1:resampleIterCount_GLM
            temp_trialIndex_resample = sort(randi(temp_trialNum,1,temp_trialNum));
            
            tempNONNAN_meta_trialLevel_resampled = tempNONNAN_new_meta_trialLevel(temp_trialIndex_resample);
            tempNONNAN_memoryPrecision_trialLevel_resampled = tempNONNAN_new_memoryPrecision_trialLevel(temp_trialIndex_resample);
            tempNONNAN_meta_trialLevel_baseline_resampled = tempNONNAN_new_meta_trialLevel_baseline(temp_trialIndex_resample);
            
            temp_x1 = tempNONNAN_memoryPrecision_trialLevel_resampled;
            temp_x2 = tempNONNAN_meta_trialLevel_resampled;
            temp_y = tempNONNAN_meta_trialLevel_baseline_resampled;
            
            temp_x1_zscored = (temp_x1-mean(temp_x1))./std(temp_x1);
            temp_x2_zscored = (temp_x2-mean(temp_x2))./std(temp_x2);
            temp_y_zscored = (temp_y-mean(temp_y))./std(temp_y);
            
            
            x = [temp_x1_zscored,temp_x2_zscored];
            y = temp_y_zscored;
            
            temp_mdl = fitglm(x,y,'interactions');
            %temp_mdl = fitglm(x,y,'linear');
            r2_caseZ6_resampled(iter) = temp_mdl.Rsquared.Adjusted;
            beta0_caseZ6_resampled(iter) = temp_mdl.Coefficients.Estimate(1);
            beta1_caseZ6_resampled(iter) = temp_mdl.Coefficients.Estimate(2);
            beta2_caseZ6_resampled(iter) = temp_mdl.Coefficients.Estimate(3);
            beta3_caseZ6_resampled(iter) = temp_mdl.Coefficients.Estimate(4);
        end
        
        [~,p1_caseZ6,~,~] = ttest(beta1_caseZ6_resampled,0);
        [~,p2_caseZ6,~,~] = ttest(beta2_caseZ6_resampled,0);
        [~,p3_caseZ6,~,~] = ttest(beta3_caseZ6_resampled,0);
        
        [~,p12_caseZ6,~,~] = ttest(beta1_caseZ6_resampled,beta2_caseZ6_resampled);
        
        
        %% Plot. Regression: Baseline meta ~ WM strength + Meta-WM
        if if_plot == 1
            
            temp_chanceLevel = 0;
            
            fig = figure('Name','locDistri','NumberTitle','off');
            set(gcf,'Position',[1150+0 400+0 240*0.80*1.2*0.9*0.9*1.02*0.8*1.05*1.05*1.05 336*1.11*0.9*0.91*0.95*0.95*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
            nexttile
            
            temp_p1 = p1_caseZ6;
            temp_p2 = p2_caseZ6;
            temp_p3 = p3_caseZ6;
            
            
            temp_1 = beta1_caseZ6_resampled;
            temp_2 = beta2_caseZ6_resampled;
            temp_3 = beta3_caseZ6_resampled;
            
            temp_r2 = mean(r2_caseZ6_resampled);
            
            
            temp_y_min = min([temp_1 temp_2 temp_3]); %#ok<*NASGU>
            temp_y_max = max([temp_1 temp_2 temp_3]);
            
            %temp_y_min = -0.1;
            %temp_y_max = 0.4790;
            
            [~,temp_y1_max] = bounds(temp_1);
            [~,temp_y2_max] = bounds(temp_2);
            [~,temp_y3_max] = bounds(temp_3);
            
            %temp_if_plot_violinplot0_pairline1 = if_plot_violinplot0_pairline1;
            temp_if_plot_violinplot0_pairline1 = 0;
            
            if temp_if_plot_violinplot0_pairline1 == 0
                temp_data = [temp_1';temp_2';temp_3'];
                
                g1 = repmat({'A'},length(temp_1),1);
                g2 = repmat({'B'},length(temp_2),1);
                g3 = repmat({'C'},length(temp_3),1);
                
                temp_label = [g1;g2;g3];
                
                temptemp_color1 = [1 1 1]*0.5;
                temptemp_color2 = repmat(temptemp_color1, 3, 1);
                
                violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                    'GroupOrder',[{'A'};{'B'};{'C'}]);
                
            elseif temp_if_plot_violinplot0_pairline1 == 1
                %for tempi=1:length(temp_1)
                %    plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                %    hold on
                %end
                
                plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                
                
                scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
                scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
                scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
                hold on
            end
            
            
            plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
            hold on
            
            
            tempTxt = sprintf('');
            if temp_p1 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p1 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p1 < 0.05
                tempTxt = sprintf('*');
            end
            text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p2 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p2 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p2 < 0.05
                tempTxt = sprintf('*');
            end
            text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p3 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p3 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p3 < 0.05
                tempTxt = sprintf('*');
            end
            text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            tempTxt = sprintf('');
            if temp_p12 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p12 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p12 < 0.05
                tempTxt = sprintf('*');
            end
            text(1.5,temp_y_min-(temp_y_max-temp_y_min)*0.21,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            plot([1 2],temp_y_min-(temp_y_max-temp_y_min)*0.1*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
            hold on
            
            
            xticks([1 2 3]);
            
            xtl = ["WM strength", "Meta", "Interaction"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
            xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.37;
            text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
            
            set(gca,'xticklabel','');
            
            
            set(gca,'linewidth',1.5)
            
            xlim([0.35 3.65]);
            %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.2
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.3 temp_y_max+(temp_y_max-temp_y_min)*0.3]);%0.2
            set(gca, 'FontSize', 8) %14
            set(gca,'box','off');% 取消右、上边框
            
            
            
            ylabel('Beta', 'FontSize', 9);
            temp_title = title(sprintf('Baseline meta regression'),'fontsize',8);
            temp_title.Interpreter = 'none';
            %subtitle(sprintf('r2 = %.3f',temp_r2),'fontsize',8);%9
            %subtitle(sprintf('all trials, r2 = %.3f',temp_r2),'fontsize',8);%9
            subtitle(sprintf('r2 = %.3f',temp_r2),'fontsize',8);%9
            
        end
        
        
    end
    
    a = 1;
    
    %% Neuronal linear regression, using neuron activity to estimate task variables linearly.
    if true
        t_linearAxis = tic;
        
        if_meta_decoder0_behavior1 = 0;%0,1
        
        if if_meta_decoder0_behavior1 == 0
            temp_meta_trialLevel_baseline = meta_trialLevel_baseline;
            temp_meta_trialLevel = meta_trialLevel;
        elseif if_meta_decoder0_behavior1 == 1
            temp_meta_trialLevel_baseline = double(choiceMemoryBoolIndex');
            temp_meta_trialLevel = double(choiceMemoryBoolIndex');
        end
        
        memoryPrecision_trialLevel = memoryPrecision_trialLevel;
        
        
        if_fitglm0_lasso1 = 1;%0
        
        if_resample0_twoSession1 = 0;%1,0
        
        temp_resampleNum = 16*7;%16*7(63 secs in lasso),16*2(20 secs in lasso)
        if_resample_replace = 0;%0
        
        temp_if_allTrials0_partTrials1 = 0;%0
        
        if if_meta_decoder0_behavior1 == 1
            temp_if_allTrials0_partTrials1 = 1;
        end
        
        temp_if_zscore = 1;%1
        
        temp_if_resample_balance = 0;%0
        temp_ifBalance_lowHigh0_seq1 = 1;
        
        if if_resample0_twoSession1 == 1
            temp_resampleNum = 1;
            if_resample_replace = 0;
        end
        
        %fprintf('if_fitglm0_lasso1=%d,if_resample0_twoSession1=%d,if_resample_replace=%d, temp_if_allTrials0_partTrials1=%d, temp_if_zscore=%d, temp_if_resample_balance=%d, temp_ifBalance_lowHigh0_seq1=%d.\n',...
        %    if_fitglm0_lasso1,if_resample0_twoSession1,if_resample_replace,temp_if_allTrials0_partTrials1,temp_if_zscore,temp_if_resample_balance,temp_ifBalance_lowHigh0_seq1);
        
        fprintf('if_meta_decoder0_behavior1=%d,if_fitglm0_lasso1=%d,if_resample0_twoSession1=%d,temp_if_allTrials0_partTrials1=%d,temp_if_zscore=%d,',...
            if_meta_decoder0_behavior1,if_fitglm0_lasso1,if_resample0_twoSession1,temp_if_allTrials0_partTrials1,temp_if_zscore);
        
        if if_resample0_twoSession1 == 0
            fprintf('temp_resampleNum=%d,if_resample_replace=%d,temp_if_resample_balance=%d,temp_ifBalance_lowHigh0_seq1=%d.\n',...
                temp_resampleNum,if_resample_replace,temp_if_resample_balance,temp_ifBalance_lowHigh0_seq1);
        end
        
        
        temp_roiNum = size(F_dff_decisionPeriodA,1);
        
        lasso_linearAxis_jjb_Name_v = autoGetFunName_myScripts('lasso_linearAxis_jjb', [targetPATH '\functions']);
        fun_lasso_linearAxis_jjb = str2func(lasso_linearAxis_jjb_Name_v);
        
        
        temp1 = memoryPrecision_trialLevel(~isnan(memoryPrecision_trialLevel));
        temp2 = temp_meta_trialLevel(~isnan(memoryPrecision_trialLevel));
        temp3 = corr(temp1,temp2);
        
        
        F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
        
        temptemp_range = (baselinePeriod_interval(2)-8):baselinePeriod_interval(3);
        F_dff_baselineBin = mean(F_dff_baselinePeriod(:,:,temptemp_range),3);
        
        if temp_if_zscore == 1
            temp_F_dff_decisionBin1 = F_dff_decisionBin1_zScore;
            temp_F_dff_baselineBin = F_dff_baselineBin_zScore;
        elseif temp_if_zscore == 0
            temp_F_dff_decisionBin1 = F_dff_decisionBin1;
            temp_F_dff_baselineBin = F_dff_baselineBin;
        end
        
        %% Baseline meta-WM, neuronal linear regression, using neuron activity to estimate task variables linearly
        if temp_if_allTrials0_partTrials1 == 0
            x = temp_F_dff_baselineBin';
            y = temp_meta_trialLevel_baseline;
        elseif temp_if_allTrials0_partTrials1 == 1
            x = temp_F_dff_baselineBin(:,choiceBoolIndex)';
            y = temp_meta_trialLevel_baseline(choiceBoolIndex);
        end
        if temp_if_zscore == 1
            y_zscore = (y-mean(y,'omitnan'))./std(y,'omitnan');%zscore
            if temp_if_resample_balance==1 && temp_ifBalance_lowHigh0_seq1==0
                a = 1;
            else
                y = y_zscore;
            end
        end
        
        if if_fitglm0_lasso1 == 0
            temp_mdl = fitglm(x,y);
            r2_prior_linearAxis = temp_mdl.Rsquared.Adjusted;
            beta_prior_linearAxis = temp_mdl.Coefficients.Estimate(2:end);
            beta0_prior_linearAxis = temp_mdl.Coefficients.Estimate(1);
        elseif if_fitglm0_lasso1 == 1
            [beta_prior_linearAxis,r2_prior_linearAxis,beta0_prior_linearAxis] = fun_lasso_linearAxis_jjb(x,y);
        end
        
        %% WM strength, neuronal linear regression, using neuron activity to estimate task variables linearly
        if temp_if_allTrials0_partTrials1 == 0
            x = temp_F_dff_decisionBin1';
            y = memoryPrecision_trialLevel;
        elseif temp_if_allTrials0_partTrials1 == 1
            %x = temp_F_dff_decisionBin1(:,trialIndex_bool_memoryCorrect)';
            %y = memoryPrecision_trialLevel(trialIndex_bool_memoryCorrect);
            x = temp_F_dff_decisionBin1(:,choiceBoolIndex)';
            y = memoryPrecision_trialLevel(choiceBoolIndex);
        end
        if temp_if_zscore == 1
            y_zscore = (y-mean(y,'omitnan'))./std(y,'omitnan');%zscore
            if temp_if_resample_balance==1 && temp_ifBalance_lowHigh0_seq1==0
            else
                y = y_zscore;
            end
        end
        
        if if_fitglm0_lasso1 == 0
            temp_mdl = fitglm(x,y);
            r2_memory_linearAxis = temp_mdl.Rsquared.Adjusted;
            beta_memory_linearAxis = temp_mdl.Coefficients.Estimate(2:end);
            beta0_memory_linearAxis = temp_mdl.Coefficients.Estimate(1);
        elseif if_fitglm0_lasso1 == 1
            [beta_memory_linearAxis,r2_memory_linearAxis,beta0_memory_linearAxis] = fun_lasso_linearAxis_jjb(x,y);
        end
        
        
        %% Meta-WM, neuronal linear regression, using neuron activity to estimate task variables linearly
        if temp_if_allTrials0_partTrials1 == 0
            x = temp_F_dff_decisionBin1';
            y = temp_meta_trialLevel;
        elseif temp_if_allTrials0_partTrials1 == 1
            x = temp_F_dff_decisionBin1(:,choiceBoolIndex)';
            y = temp_meta_trialLevel(choiceBoolIndex);
        end
        if temp_if_zscore == 1
            y_zscore = (y-mean(y,'omitnan'))./std(y,'omitnan');%zscore
            if temp_if_resample_balance==1 && temp_ifBalance_lowHigh0_seq1==0
            else
                y = y_zscore;
            end
        end
        
        if if_fitglm0_lasso1 == 0
            temp_mdl = fitglm(x,y);
            r2_meta_linearAxis = temp_mdl.Rsquared.Adjusted;
            beta_meta_linearAxis = temp_mdl.Coefficients.Estimate(2:end);
            beta0_meta_linearAxis = temp_mdl.Coefficients.Estimate(1);
        elseif if_fitglm0_lasso1 == 1
            [beta_meta_linearAxis,r2_meta_linearAxis,beta0_meta_linearAxis] = fun_lasso_linearAxis_jjb(x,y);
        end
        
        %%
        beta_prior_linearAxis;
        beta_memory_linearAxis;
        beta_meta_linearAxis;
        
        temp_vectors = [beta_prior_linearAxis,beta_memory_linearAxis,beta_meta_linearAxis];
        
        a = temp_vectors(:,1);
        b = temp_vectors(:,2);
        c = temp_vectors(:,3);
        
        temp_degree_ab = subspace(a,b)*180/pi;
        temp_degree_ac = subspace(a,c)*180/pi;
        temp_degree_bc = subspace(b,c)*180/pi;
        
        
        
        %% Resample beta_prior_linearAxis in all trials (from half trials)
        if temp_if_allTrials0_partTrials1 == 0
            x = temp_F_dff_baselineBin';
            y = temp_meta_trialLevel_baseline;
        elseif temp_if_allTrials0_partTrials1 == 1
            x = temp_F_dff_baselineBin(:,choiceBoolIndex)';
            y = temp_meta_trialLevel_baseline(choiceBoolIndex);
        end
        if temp_if_zscore == 1
            y_zscore = (y-mean(y,'omitnan'))./std(y,'omitnan');%zscore
            if temp_if_resample_balance==1 && temp_ifBalance_lowHigh0_seq1==0
            else
                y = y_zscore;
            end
        end
        
        beta_prior_linearAxis_halfA_resampled = cell(temp_resampleNum,1);
        beta_prior_linearAxis_halfB_resampled = cell(temp_resampleNum,1);
        
        parfor tempIter=1:temp_resampleNum
            warning('off');
            temptempBoolIndex_B = [];
            beta_A = [];
            beta_B = [];
            
            temptemp_trialNum = size(x,1);
            %temptempIndex = sort(randperm(temptemp_trialNum,round(temptemp_trialNum/2)));
            temptempIndex = sort(randperm(temptemp_trialNum,round(temptemp_trialNum*0.5)));
            
            temptempBoolIndex = false(temptemp_trialNum,1);
            temptempBoolIndex(temptempIndex) = true;
            temptempBoolIndex_A = temptempBoolIndex;
            
            
            if if_resample0_twoSession1 == 1
                temptempBoolIndex = false(temptemp_trialNum,1);
                temptempBoolIndex(1:round(temptemp_trialNum/2)) = true;
                temptempBoolIndex_A = temptempBoolIndex;
            end
            
            
            if if_resample_replace == 0
                temptempBoolIndex_B = ~temptempBoolIndex;
            elseif if_resample_replace == 1
                %temptempIndex = sort(randperm(temptemp_trialNum,round(temptemp_trialNum/2)));
                temptempIndex = sort(randperm(temptemp_trialNum,round(temptemp_trialNum*0.5)));
                
                temptempBoolIndex = false(temptemp_trialNum,1);
                temptempBoolIndex(temptempIndex) = true;
                temptempBoolIndex_B = temptempBoolIndex;
            end
            
            
            
            if if_fitglm0_lasso1 == 0
                temp_mdl_A = fitglm(x(temptempBoolIndex_A,:),y(temptempBoolIndex_A));
                beta_A = temp_mdl_A.Coefficients.Estimate(2:end);
                temp_mdl_B = fitglm(x(temptempBoolIndex_B,:),y(temptempBoolIndex_B));
                beta_B = temp_mdl_B.Coefficients.Estimate(2:end);
            elseif if_fitglm0_lasso1 == 1
                [beta_A,~] = fun_lasso_linearAxis_jjb(x(temptempBoolIndex_A,:),y(temptempBoolIndex_A));
                [beta_B,~] = fun_lasso_linearAxis_jjb(x(temptempBoolIndex_B,:),y(temptempBoolIndex_B));
            end
            
            
            beta_prior_linearAxis_halfA_resampled{tempIter} = beta_A;
            beta_prior_linearAxis_halfB_resampled{tempIter} = beta_B;
            
        end
        
        betaAngle_prior_linearAxis_halfAB_resampled = nan(temp_resampleNum,1);
        for tempi=1:temp_resampleNum
            a = beta_prior_linearAxis_halfA_resampled{tempi};
            b = beta_prior_linearAxis_halfB_resampled{tempi};
            betaAngle_prior_linearAxis_halfAB_resampled(tempi) = subspace(a,b)*180/pi;
        end
        
        %% Resample beta_memory_linearAxis in all trials (from half trials)
        if temp_if_allTrials0_partTrials1 == 0
            x = temp_F_dff_decisionBin1';
            y = memoryPrecision_trialLevel;
            tempBoolIndex_valid = true(size(y));
        elseif temp_if_allTrials0_partTrials1 == 1
            %tempBoolIndex_valid = trialIndex_bool_memoryCorrect';
            tempBoolIndex_valid = choiceBoolIndex';
            
            x = temp_F_dff_decisionBin1(:,tempBoolIndex_valid)';
            y = memoryPrecision_trialLevel(tempBoolIndex_valid);
        end
        if temp_if_zscore == 1
            y_zscore = (y-mean(y,'omitnan'))./std(y,'omitnan');%zscore
            if temp_if_resample_balance==1 && temp_ifBalance_lowHigh0_seq1==0
            else
                y = y_zscore;
            end
        end
        
        
        if true
            %lowThreshold_meta;
            y_threshold = lowThreshold_memoryPrecision;
            
            tempBoolIndex_y_low = y < y_threshold;
            tempBoolIndex_y_high = y >= y_threshold;
            
            %         temp_seqNum = 41;%seq-based
            %         tempBoolIndex_y_category = false(length(y),temp_categoryNum);
            %
            %         seqIndex;
            %         for tempi=1:temp_seqNum
            %
            %         end
        end
        
        beta_memory_linearAxis_halfA_resampled = cell(temp_resampleNum,1);
        beta_memory_linearAxis_halfB_resampled = cell(temp_resampleNum,1);
        
        parfor tempIter=1:temp_resampleNum
            %for tempIter=1:temp_resampleNum
            warning('off');
            temptempBoolIndex_B = [];
            temptempIndex = [];
            beta_A = [];
            beta_B = [];
            
            temptemp_trialNum = size(x,1);
            temptempIndex = sort(randperm(temptemp_trialNum,round(temptemp_trialNum/2)));
            
            
            temptempBoolIndex = false(temptemp_trialNum,1);
            temptempBoolIndex(temptempIndex) = true;
            temptempBoolIndex_A = temptempBoolIndex;
            %temptempBoolIndex_B = ~temptempBoolIndex;
            
            if if_resample0_twoSession1 == 1
                temptempBoolIndex = false(temptemp_trialNum,1);
                temptempBoolIndex(1:round(temptemp_trialNum/2)) = true;
                %temptempBoolIndex(1:500) = true;
                temptempBoolIndex_A = temptempBoolIndex;
            end
            
            if if_resample_replace == 0
                
                
                if temp_if_resample_balance == 0
                    temptempBoolIndex_B = ~temptempBoolIndex;
                    
                elseif temp_if_resample_balance == 1
                    %% Balance seq
                    if temp_ifBalance_lowHigh0_seq1 == 1
                        %seqIndex = seqIndex;
                        temp_seqNum = 41;%seq-based
                        
                        temptempBoolIndex_A_allSeq = false(length(y),temp_seqNum);
                        temptempBoolIndex_B_allSeq = temptempBoolIndex_A_allSeq;
                        for temptempi=1:temp_seqNum
                            temptempBoolIndex_seq = seqIndex(tempBoolIndex_valid)' == temptempi;
                            
                            temptemp_trialNum_seq = sum(temptempBoolIndex_seq);
                            
                            temptempIndex_seq_A = sort(randperm(temptemp_trialNum_seq,round(temptemp_trialNum_seq/2)));
                            
                            tempIndex_y_seq = find(temptempBoolIndex_seq==true);
                            temp1 = tempIndex_y_seq(temptempIndex_seq_A);
                            temptempBoolIndex_seq_A = false(temptemp_trialNum,1);
                            temptempBoolIndex_seq_A(temp1) = true;
                            
                            temptempBoolIndex_seq_B = (~temptempBoolIndex_seq_A)&temptempBoolIndex_seq;
                            %sum(temptempBoolIndex_seq_B&temptempBoolIndex_seq)
                            %sum(temptempBoolIndex_seq_A&temptempBoolIndex_seq_B)
                            
                            
                            temptempBoolIndex_A_allSeq(:,temptempi) = temptempBoolIndex_seq_A;
                            temptempBoolIndex_B_allSeq(:,temptempi) = temptempBoolIndex_seq_B;
                        end
                        
                        temptempBoolIndex_A = sum(temptempBoolIndex_A_allSeq,2) > 0;
                        temptempBoolIndex_B = sum(temptempBoolIndex_B_allSeq,2) > 0;
                    end
                    
                    
                    %% Balance low WM strength and high WM strength
                    if temp_ifBalance_lowHigh0_seq1 == 0
                        tempBoolIndex_y_low;
                        tempBoolIndex_y_high;
                        
                        temptemp_trialNum_low = sum(tempBoolIndex_y_low);
                        temptemp_trialNum_high = sum(tempBoolIndex_y_high);
                        
                        temptempIndex_low_A = sort(randperm(temptemp_trialNum_low,round(temptemp_trialNum_low/2)));
                        
                        tempIndex_y_low = find(tempBoolIndex_y_low==true);
                        temp1 = tempIndex_y_low(temptempIndex_low_A);
                        temptempBoolIndex_low_A = false(temptemp_trialNum,1);
                        temptempBoolIndex_low_A(temp1) = true;
                        
                        temptempBoolIndex_low_B = (~temptempBoolIndex_low_A)&tempBoolIndex_y_low;
                        %sum(temptempBoolIndex_low_B&tempBoolIndex_y_low)
                        
                        
                        
                        temptempIndex_high_A = sort(randperm(temptemp_trialNum_high,round(temptemp_trialNum_high/2)));
                        
                        tempIndex_y_high = find(tempBoolIndex_y_high==true);
                        temp1 = tempIndex_y_high(temptempIndex_high_A);
                        temptempBoolIndex_high_A = false(temptemp_trialNum,1);
                        temptempBoolIndex_high_A(temp1) = true;
                        
                        temptempBoolIndex_high_B = (~temptempBoolIndex_high_A)&tempBoolIndex_y_high;
                        %sum(temptempBoolIndex_high_B&tempBoolIndex_y_high)
                        
                        
                        
                        temptempBoolIndex_A = temptempBoolIndex_low_A | temptempBoolIndex_high_A;
                        temptempBoolIndex_B = temptempBoolIndex_low_B | temptempBoolIndex_high_B;
                        
                    end
                    
                end
                
                
            elseif if_resample_replace == 1
                %temptempIndex = sort(randperm(temptemp_trialNum,round(temptemp_trialNum/2)));
                temptempIndex = sort(randperm(temptemp_trialNum,round(temptemp_trialNum*0.5)));
                
                temptempBoolIndex = false(temptemp_trialNum,1);
                temptempBoolIndex(temptempIndex) = true;
                temptempBoolIndex_B = temptempBoolIndex;
            end
            
            temp_x_A = x(temptempBoolIndex_A,:);
            temp_y_A = y(temptempBoolIndex_A);
            temp_x_B = x(temptempBoolIndex_B,:);
            temp_y_B = y(temptempBoolIndex_B);
            
            
            %         if true % test for AB z-score
            %             temp_y_A = (temp_y_A-mean(temp_y_A,'omitnan'))./std(temp_y_A,'omitnan');%zscore
            %             temp_y_B = (temp_y_B-mean(temp_y_B,'omitnan'))./std(temp_y_B,'omitnan');%zscore
            %
            %             temp_x_A_zScore = nan(size(temp_x_A));
            %             for temptempi=1:size(temp_x_A_zScore,2)
            %                 temp1 = temp_x_A(:,temptempi);
            %                 temp2 =(temp1-mean(temp1))/std(temp1); % z-score
            %                 temp_x_A_zScore(:,temptempi) = temp2;
            %             end
            %             temp_x_A = temp_x_A_zScore;
            %
            %             temp_x_B_zScore = nan(size(temp_x_B));
            %             for temptempi=1:size(temp_x_B_zScore,2)
            %                 temp1 = temp_x_B(:,temptempi);
            %                 temp2 =(temp1-mean(temp1))/std(temp1); % z-score
            %                 temp_x_B_zScore(:,temptempi) = temp2;
            %             end
            %             temp_x_B = temp_x_B_zScore;
            %         end
            
            
            if if_fitglm0_lasso1 == 0
                temp_mdl_A = fitglm(temp_x_A,temp_y_A);
                beta_A = temp_mdl_A.Coefficients.Estimate(2:end);
                temp_mdl_B = fitglm(temp_x_B,temp_y_B);
                beta_B = temp_mdl_B.Coefficients.Estimate(2:end);
            elseif if_fitglm0_lasso1 == 1
                [beta_A,~] = fun_lasso_linearAxis_jjb(temp_x_A,temp_y_A);
                [beta_B,~] = fun_lasso_linearAxis_jjb(temp_x_B,temp_y_B);
            end
            
            beta_memory_linearAxis_halfA_resampled{tempIter} = beta_A;
            beta_memory_linearAxis_halfB_resampled{tempIter} = beta_B;
            
            %corr(beta_A,beta_B)
            %subspace(beta_A,beta_B)*180/pi
            %tempi = tempIter;
        end
        
        betaAngle_memory_linearAxis_halfAB_resampled = nan(temp_resampleNum,1);
        for tempi=1:temp_resampleNum
            a = beta_memory_linearAxis_halfA_resampled{tempi};
            b = beta_memory_linearAxis_halfB_resampled{tempi};
            betaAngle_memory_linearAxis_halfAB_resampled(tempi) = subspace(a,b)*180/pi;
        end
        
        
        %% Resample beta_meta_linearAxis in all trials (from half trials)
        if temp_if_allTrials0_partTrials1 == 0
            x = temp_F_dff_decisionBin1';
            y = temp_meta_trialLevel;
        elseif temp_if_allTrials0_partTrials1 == 1
            x = temp_F_dff_decisionBin1(:,choiceBoolIndex)';
            y = temp_meta_trialLevel(choiceBoolIndex);
        end
        if temp_if_zscore == 1
            y_zscore = (y-mean(y,'omitnan'))./std(y,'omitnan');%zscore
            if temp_if_resample_balance==1 && temp_ifBalance_lowHigh0_seq1==0
            else
                y = y_zscore;
            end
        end
        
        beta_meta_linearAxis_halfA_resampled = cell(temp_resampleNum,1);
        beta_meta_linearAxis_halfB_resampled = cell(temp_resampleNum,1);
        
        parfor tempIter=1:temp_resampleNum
            warning('off');
            temptempBoolIndex_B = [];
            beta_A = [];
            beta_B = [];
            
            temptemp_trialNum = size(x,1);
            temptempIndex = sort(randperm(temptemp_trialNum,round(temptemp_trialNum/2)));
            
            temptempBoolIndex = false(temptemp_trialNum,1);
            temptempBoolIndex(temptempIndex) = true;
            temptempBoolIndex_A = temptempBoolIndex;
            %temptempBoolIndex_B = ~temptempBoolIndex;
            
            
            if if_resample0_twoSession1 == 1
                temptempBoolIndex = false(temptemp_trialNum,1);
                temptempBoolIndex(1:round(temptemp_trialNum/2)) = true;
                temptempBoolIndex_A = temptempBoolIndex;
            end
            
            
            if if_resample_replace == 0
                temptempBoolIndex_B = ~temptempBoolIndex;
            elseif if_resample_replace == 1
                %temptempIndex = sort(randperm(temptemp_trialNum,round(temptemp_trialNum/2)));
                temptempIndex = sort(randperm(temptemp_trialNum,round(temptemp_trialNum*0.5)));
                
                temptempBoolIndex = false(temptemp_trialNum,1);
                temptempBoolIndex(temptempIndex) = true;
                temptempBoolIndex_B = temptempBoolIndex;
            end
            
            
            if if_fitglm0_lasso1 == 0
                temp_mdl_A = fitglm(x(temptempBoolIndex_A,:),y(temptempBoolIndex_A));
                beta_A = temp_mdl_A.Coefficients.Estimate(2:end);
                temp_mdl_B = fitglm(x(temptempBoolIndex_B,:),y(temptempBoolIndex_B));
                beta_B = temp_mdl_B.Coefficients.Estimate(2:end);
            elseif if_fitglm0_lasso1 == 1
                [beta_A,~] = fun_lasso_linearAxis_jjb(x(temptempBoolIndex_A,:),y(temptempBoolIndex_A));
                [beta_B,~] = fun_lasso_linearAxis_jjb(x(temptempBoolIndex_B,:),y(temptempBoolIndex_B));
            end
            
            beta_meta_linearAxis_halfA_resampled{tempIter} = beta_A;
            beta_meta_linearAxis_halfB_resampled{tempIter} = beta_B;
            
        end
        
        betaAngle_meta_linearAxis_halfAB_resampled = nan(temp_resampleNum,1);
        for tempi=1:temp_resampleNum
            a = beta_meta_linearAxis_halfA_resampled{tempi};
            b = beta_meta_linearAxis_halfB_resampled{tempi};
            betaAngle_meta_linearAxis_halfAB_resampled(tempi) = subspace(a,b)*180/pi;
        end
        
        %% Summary of halfAB
        temp_degree_ab;
        temp_degree_ac;
        temp_degree_bc;
        
        beta_prior_linearAxis;
        beta_memory_linearAxis;
        beta_meta_linearAxis;
        
        temp_degree_aa = betaAngle_prior_linearAxis_halfAB_resampled;
        temp_degree_aa_mean = mean(temp_degree_aa);
        temp_degree_aa_std = std(temp_degree_aa);
        
        temp_degree_bb = betaAngle_memory_linearAxis_halfAB_resampled;
        temp_degree_bb_mean = mean(temp_degree_bb);
        temp_degree_bb_std = std(temp_degree_bb);
        
        temp_degree_cc = betaAngle_meta_linearAxis_halfAB_resampled;
        temp_degree_cc_mean = mean(temp_degree_cc);
        temp_degree_cc_std = std(temp_degree_cc);
        
        
        [~,temp_p_aa_ab] = ttest(temp_degree_aa,temp_degree_ab);
        [~,temp_p_aa_ac] = ttest(temp_degree_aa,temp_degree_ac);
        
        
        [~,temp_p_bb_ab] = ttest(temp_degree_bb,temp_degree_ab);
        [~,temp_p_bb_bc] = ttest(temp_degree_bb,temp_degree_bc);
        
        [~,temp_p_cc_ac] = ttest(temp_degree_cc,temp_degree_ac);
        [~,temp_p_cc_bc] = ttest(temp_degree_cc,temp_degree_bc);
        
        %fprintf('Degree between beta vectors: \nbetween [%.1f, %.1f, %.1f], within [%.1f, %.1f, %.1f].\n',...
        %    temp_degree_ab,temp_degree_ac,temp_degree_bc,temp_degree_aa_mean,temp_degree_bb_mean,temp_degree_cc_mean);
        
        
        fprintf('Degree between beta vectors (subspace view): within [%.1f, %.1f, %.1f], between [%.1f, %.1f, %.1f].\n',...
            temp_degree_aa_mean,temp_degree_bb_mean,temp_degree_cc_mean,temp_degree_ab,temp_degree_ac,temp_degree_bc);
        
        temp_degree_aa_prior = temp_degree_aa;
        temp_degree_bb_memory = temp_degree_bb;
        temp_degree_cc_meta = temp_degree_cc;
        
        temp_degree_ab_priorMemory = temp_degree_ab;
        temp_degree_ac_priorMeta = temp_degree_ac;
        temp_degree_bc_memoryMeta = temp_degree_bc;
        
        
        beta_prior_linearAxis;
        beta_memory_linearAxis;
        beta_meta_linearAxis;
        
        temp_degree_aa_prior;
        temp_degree_bb_memory;
        temp_degree_cc_meta;
        
        temp_degree_ab_priorMemory;
        temp_degree_ac_priorMeta;
        temp_degree_bc_memoryMeta;
        
        if false
            %temp_beta_A = beta_prior_linearAxis;
            %temp_beta_A = beta_prior_linearAxis([7 8]);
            
            temp_beta_A = [1 2]';
            
            temp_beta_B = [-1 -2]';
            %temp_beta_B = [2 -1]';
            
            temp_beta_B_unit = temp_beta_B / norm(temp_beta_B);
            
            temp1_B = temp_beta_B_unit * temp_beta_B_unit';
            
            temp2 = var(temp1_B*temp_beta_A);
            temp3 = var(temp_beta_A);
            
            temp_ratio = temp2/temp3;
        end
        
        
        
        if true
            %temp_VAF_ratio_priorMemory
            temp_beta_A = beta_prior_linearAxis;
            temp_beta_B = beta_memory_linearAxis;
            temp_beta_B_unit = temp_beta_B / norm(temp_beta_B);
            temp1_B = temp_beta_B_unit * temp_beta_B_unit';
            temp_VAF_ratio_AB = var(temp1_B*temp_beta_A)/var(temp_beta_A);
            temp_VAF_ratio_priorMemory = temp_VAF_ratio_AB;
            
            %temp_VAF_ratio_priorMeta
            temp_beta_A = beta_prior_linearAxis;
            temp_beta_B = beta_meta_linearAxis;
            temp_beta_B_unit = temp_beta_B / norm(temp_beta_B);
            temp1_B = temp_beta_B_unit * temp_beta_B_unit';
            temp_VAF_ratio_AB = var(temp1_B*temp_beta_A)/var(temp_beta_A);
            temp_VAF_ratio_priorMeta = temp_VAF_ratio_AB;
            
            %temp_VAF_ratio_memoryMeta
            temp_beta_A = beta_memory_linearAxis;
            temp_beta_B = beta_meta_linearAxis;
            temp_beta_B_unit = temp_beta_B / norm(temp_beta_B);
            temp1_B = temp_beta_B_unit * temp_beta_B_unit';
            temp_VAF_ratio_AB = var(temp1_B*temp_beta_A)/var(temp_beta_A);
            temp_VAF_ratio_memoryMeta = temp_VAF_ratio_AB;
            
            %temp_VAF_ratio_priorPrior
            temp_VAF_ratio_priorPrior = nan(temp_resampleNum,1);
            for tempi=1:temp_resampleNum
                temp_beta_A = beta_prior_linearAxis_halfA_resampled{tempi};
                temp_beta_B = beta_prior_linearAxis_halfB_resampled{tempi};
                temp_beta_B_unit = temp_beta_B / norm(temp_beta_B);
                temp1_B = temp_beta_B_unit * temp_beta_B_unit';
                temp_VAF_ratio_AB = var(temp1_B*temp_beta_A)/var(temp_beta_A);
                temp_VAF_ratio_priorPrior(tempi) = temp_VAF_ratio_AB;
            end
            
            %temp_VAF_ratio_memoryMemory
            temp_VAF_ratio_memoryMemory = nan(temp_resampleNum,1);
            for tempi=1:temp_resampleNum
                temp_beta_A = beta_memory_linearAxis_halfA_resampled{tempi};
                temp_beta_B = beta_memory_linearAxis_halfB_resampled{tempi};
                temp_beta_B_unit = temp_beta_B / norm(temp_beta_B);
                temp1_B = temp_beta_B_unit * temp_beta_B_unit';
                temp_VAF_ratio_AB = var(temp1_B*temp_beta_A)/var(temp_beta_A);
                temp_VAF_ratio_memoryMemory(tempi) = temp_VAF_ratio_AB;
            end
            
            
            %temp_VAF_ratio_metaMeta
            temp_VAF_ratio_metaMeta = nan(temp_resampleNum,1);
            for tempi=1:temp_resampleNum
                temp_beta_A = beta_meta_linearAxis_halfA_resampled{tempi};
                temp_beta_B = beta_meta_linearAxis_halfB_resampled{tempi};
                temp_beta_B_unit = temp_beta_B / norm(temp_beta_B);
                temp1_B = temp_beta_B_unit * temp_beta_B_unit';
                temp_VAF_ratio_AB = var(temp1_B*temp_beta_A)/var(temp_beta_A);
                temp_VAF_ratio_metaMeta(tempi) = temp_VAF_ratio_AB;
            end
            
            fprintf('VAF ratio between beta vectors (subspace view): within [%.2f, %.2f, %.2f], between [%.2f, %.2f, %.2f].\n',...
                mean(temp_VAF_ratio_priorPrior),mean(temp_VAF_ratio_memoryMemory),mean(temp_VAF_ratio_metaMeta),...
                temp_VAF_ratio_priorMemory,temp_VAF_ratio_priorMeta,temp_VAF_ratio_memoryMeta);
            
            temp_VAF_ratio_priorPrior;
            temp_VAF_ratio_memoryMemory;
            temp_VAF_ratio_metaMeta;
            temp_VAF_ratio_priorMemory;
            temp_VAF_ratio_priorMeta;
            temp_VAF_ratio_memoryMeta;
            
            
        end
        
        fprintf('t_linearAxis = %.1f secs.\n\n',toc(t_linearAxis));
    end
    
    
    % %% VAF ratio of neuron activity in three axes
    % if true
    %     beta_prior_linearAxis;
    %     beta_memory_linearAxis;
    %     beta_meta_linearAxis;
    %
    %
    %     %temp_VAF_ratio_neuronPrior
    %     temp_x = temp_F_dff_baselineBin;
    %     temp_beta = beta_prior_linearAxis;
    %     temp_beta_unit = temp_beta / norm(temp_beta);
    %     temp1 = temp_beta_unit * temp_beta_unit';
    %     temp_VAF_ratio_AB = var(temp1*temp_x,0,'all')/var(temp_x,0,'all');
    %     temp_VAF_ratio_neuronPrior = temp_VAF_ratio_AB;
    %
    %     %temp_VAF_ratio_neuronMemory
    %     temp_x = temp_F_dff_decisionBin1;
    %     temp_beta = beta_memory_linearAxis;
    %     temp_beta_unit = temp_beta / norm(temp_beta);
    %     temp1 = temp_beta_unit * temp_beta_unit';
    %     temp_VAF_ratio_AB = var(temp1*temp_x,0,'all')/var(temp_x,0,'all');
    %     temp_VAF_ratio_neuronMemory = temp_VAF_ratio_AB;
    %
    %     %temp_VAF_ratio_neuronMeta
    %     temp_x = temp_F_dff_decisionBin1;
    %     %temp_x = beta_meta_linearAxis;
    %     %temp_x = [beta_meta_linearAxis-0.01 beta_meta_linearAxis+0.01];
    %     temp_beta = beta_meta_linearAxis;
    %     temp_beta_unit = temp_beta / norm(temp_beta);
    %     temp1 = temp_beta_unit * temp_beta_unit';
    %     temp_VAF_ratio_AB = var(temp1*temp_x,0,'all')/var(temp_x,0,'all');
    %     temp_VAF_ratio_neuronMeta = temp_VAF_ratio_AB;
    %
    %
    %     temp_VAF_ratio_neuronPrior;
    %     temp_VAF_ratio_neuronMemory;
    %     temp_VAF_ratio_neuronMeta;
    %
    %     fprintf('VAF ratio of neuron activity by three linear axes (subspace view): [%.4f, %.4f, %.4f].\n',...
    %         temp_VAF_ratio_neuronPrior,temp_VAF_ratio_neuronMemory,temp_VAF_ratio_neuronMeta);
    %
    % end
    
    
    %% Explained variance of neuron activity in three axes
    if true
        beta_prior_linearAxis;
        beta_memory_linearAxis;
        beta_meta_linearAxis;
        
        
        %temp_EV_neuronPrior
        temp_x = temp_F_dff_baselineBin;
        temp_beta = beta_prior_linearAxis;
        temp_beta_unit = temp_beta / norm(temp_beta);
        temp1 = temp_beta_unit * temp_beta_unit';
        %temp_EV_AB = var(temp1*temp_x,0,'all')/var(temp_x,0,'all');
        temp_EV_AB = 1 - var(temp_x-temp1*temp_x,0,'all')/var(temp_x,0,'all');
        temp_EV_neuronPrior = temp_EV_AB;
        
        %temp_EV_neuronMemory
        temp_x = temp_F_dff_decisionBin1;
        temp_beta = beta_memory_linearAxis;
        temp_beta_unit = temp_beta / norm(temp_beta);
        temp1 = temp_beta_unit * temp_beta_unit';
        %temp_EV_AB = var(temp1*temp_x,0,'all')/var(temp_x,0,'all');
        temp_EV_AB = 1 - var(temp_x-temp1*temp_x,0,'all')/var(temp_x,0,'all');
        temp_EV_neuronMemory = temp_EV_AB;
        
        %temp_EV_neuronMeta
        temp_x = temp_F_dff_decisionBin1;
        %temp_x = beta_meta_linearAxis;
        %temp_x = [beta_meta_linearAxis-0.01 beta_meta_linearAxis+0.01];
        temp_beta = beta_meta_linearAxis;
        temp_beta_unit = temp_beta / norm(temp_beta);
        temp1 = temp_beta_unit * temp_beta_unit';
        %temp_EV_AB = var(temp1*temp_x,0,'all')/var(temp_x,0,'all');
        temp_EV_AB = 1 - var(temp_x-temp1*temp_x,0,'all')/var(temp_x,0,'all');
        temp_EV_neuronMeta = temp_EV_AB;
        
        
        temp_EV_neuronPrior;
        temp_EV_neuronMemory;
        temp_EV_neuronMeta;
        
        fprintf('Explained variance of neuron activity by three linear axes (subspace view): [%.4f, %.4f, %.4f].\n',...
            temp_EV_neuronPrior,temp_EV_neuronMemory,temp_EV_neuronMeta);
        
    end
    
    
    %% PCA results of neuron activity
    if true
        temp_F_dff_baselineBin;
        temp_F_dff_decisionBin1;
        
        if_halfTrial = 0;
        
        temptemp_trialNum = length(memoryPrecision_trialLevel);
        
        
        tempBoolIndex_valid = ~isnan(memoryPrecision_trialLevel);
        if if_halfTrial == 1
            tempBoolIndex_valid(round(temptemp_trialNum/2)+1:end) = false;%only use first session trials
        end
        
        [coeff_baseline,score_baseline,latent_baseline,tsquared_baseline,explained_baseline,mu_baseline] = pca(temp_F_dff_baselineBin(:,tempBoolIndex_valid)');
        [coeff_delay,score_delay,latent_delay,tsquared_delay,explained_delay,mu_delay] = pca(temp_F_dff_decisionBin1(:,tempBoolIndex_valid)');
        
        %temp_dff_A = temp_F_dff_baselineBin(:,tempBoolIndex_valid)';
        %temp_dff_B = temp_F_dff_decisionBin1(:,tempBoolIndex_valid)';
        %temp_dff_AB = cat(3,temp_dff_A,temp_dff_B);
        %[coeff_baselineDelay,score_baselineDelay,latent_baselineDelay,tsquared_baselineDelay,explained_baselineDelay,mu_baselineDelay] = pca(temp_dff_AB);
        
        
        temp_prior = meta_trialLevel_baseline(tempBoolIndex_valid);
        temp_memoryPrecision = memoryPrecision_trialLevel(tempBoolIndex_valid);
        temp_meta = meta_trialLevel(tempBoolIndex_valid);
        
        temp_trialIndex = (1:sum(tempBoolIndex_valid))';
        
        
        temp_eventVectors = [temp_prior,temp_memoryPrecision,temp_meta,temp_trialIndex];
        
        temp_PC_range = [1,length(explained_delay)];
        
        temp_r_eachPC_VS_eventVector_baseline = nan(temp_PC_range(2),size(temp_eventVectors,2));
        temp_r_eachPC_VS_eventVector_delay = nan(temp_PC_range(2),size(temp_eventVectors,2));
        for tempi=temp_PC_range(1):temp_PC_range(2)
            score_baseline_PCx = score_baseline(:,tempi);
            [temp_r,temp_p] = corr(score_baseline_PCx,temp_eventVectors);
            temp_r_eachPC_VS_eventVector_baseline(tempi,:) = temp_r.*(temp_p < 0.001);
            
            score_delay_PCx = score_delay(:,tempi);
            [temp_r,temp_p] = corr(score_delay_PCx,temp_eventVectors);
            temp_r_eachPC_VS_eventVector_delay(tempi,:) = temp_r.*(temp_p < 0.001);
        end
        
        explained_baseline;
        explained_delay;
        temp_r_eachPC_VS_eventVector_baseline;
        temp_r_eachPC_VS_eventVector_delay;
        
        
        locDistri_trialLevel;
        
        temmp_boolIndex_location_trial = false(numFrames,length(seqIndex));
        for tempi=1:length(seqIndex_valid)
            currentSequence = target_seqSet_inOne{seqIndex(tempi)};
            temmp_boolIndex_location_trial(currentSequence,tempi) = true;
        end
        temmp_boolIndex_location_trial_T = temmp_boolIndex_location_trial';
        
        locDistri_trialLevel;
        temmp_boolIndex_location_trial_T;
        temp_r = nan(size(locDistri_trialLevel,1),1);
        temp_p = nan(size(locDistri_trialLevel,1),1);
        for tempi=1:size(locDistri_trialLevel,1)
            [temp_r(tempi),temp_p(tempi)] = corr(locDistri_trialLevel(tempi,:)',temmp_boolIndex_location_trial_T(tempi,:)');
        end
        temp_r_median = median(temp_r,'omitnan');
        
        
        temp_locDistri_trialLevel = locDistri_trialLevel(tempBoolIndex_valid,:);
        
        temp_r_eachPC_VS_locDistri_baseline = nan(temp_PC_range(2),size(temp_locDistri_trialLevel,2));
        temp_r_eachPC_VS_locDistri_delay = nan(temp_PC_range(2),size(temp_locDistri_trialLevel,2));
        for tempi=temp_PC_range(1):temp_PC_range(2)
            score_baseline_PCx = score_baseline(:,tempi);
            [temp_r,temp_p] = corr(score_baseline_PCx,temp_locDistri_trialLevel);
            temp_r_eachPC_VS_locDistri_baseline(tempi,:) = temp_r.*(temp_p < 0.001);
            
            score_delay_PCx = score_delay(:,tempi);
            [temp_r,temp_p] = corr(score_delay_PCx,temp_locDistri_trialLevel);
            temp_r_eachPC_VS_locDistri_delay(tempi,:) = temp_r.*(temp_p < 0.001);
        end
        
        temp_r_eachPC_VS_locDistri_baseline;
        temp_r_eachPC_VS_locDistri_delay;
        
        %%
        explained_baseline;
        explained_delay;
        temp_r_eachPC_VS_eventVector_baseline;
        temp_r_eachPC_VS_eventVector_delay;
        
        temp_r_eachPC_VS_locDistri_baseline;
        temp_r_eachPC_VS_locDistri_delay;
        
        temp_pca_explained_baseline_example = explained_baseline;
        temp_pca_r_eachPC_VS_eventVector_baseline_example = temp_r_eachPC_VS_eventVector_baseline;
        temp_pca_r_eachPC_VS_locDistri_baseline_example = temp_r_eachPC_VS_locDistri_baseline;
        
        temp_pca_explained_delay_example = explained_delay;
        temp_pca_r_eachPC_VS_eventVector_delay_example = temp_r_eachPC_VS_eventVector_delay;
        temp_pca_r_eachPC_VS_locDistri_delay_example = temp_r_eachPC_VS_locDistri_delay;
        
        
        
        temp_PCLimit = 25;
        
        %% Correlations between features and PCs of example FOV (baseline)
        if if_plot == 1
            fig = figure('Name','Correlations between features and PCs of example FOV (baseline)','NumberTitle','off');
            %set(gcf,'Position',[10 50 390 103*1.9]);
            set(gcf,'Position',[10 50 390 103*1.9*0.85]);
            
            t = tiledlayout(4,31,'TileSpacing','tight','Padding','compact');
            
            nexttile([2 3])
            set(gca,'Visible','off');
            
            nexttile([2 27])
            
            temp_EV = temp_pca_explained_baseline_example(1:temp_PCLimit);
            
            C = temp_pca_r_eachPC_VS_eventVector_baseline_example(1:temp_PCLimit,:)';
            %C_max = 1;
            C = abs(C);
            %C_max = max(C,[],'all');
            %C_min = min(C,[],'all');
            C_max = 1;
            C_min = 0;
            
            
            imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
            hold on
            
            set(gca,'linewidth',1.5)
            set(gca,'box','off');% 取消右、上边框
            set(gca, 'FontSize', 8) %12
            
            
            %set(gca,'XTick',1:size(C,2));
            set(gca,'xticklabel','');
            
            set(gca,'YTick',1:size(C,1));
            ytl = ["Prior","WM strength","Meta","Time"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
            ytext_xp=(xt(1))*ones(1,length(yt))-4.8;%-1.6
            ytext_yp=yt;
            text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',7.5);%8,7
            set(gca,'yticklabel','');
            
            set(gca,'TickLength',[0 0]);
            
            x_lim = xlim;
            y_lim = ylim;
            xlabel(sprintf('PC (1~25, EV = %.0f%%)',sum(temp_EV)),'Position',[mean(x_lim) max(y_lim)+0.2],'FontSize', 9);
            temp_title = title(sprintf('Correlations with PCs of one FOV (baseline)'),'FontSize',9);%11
            
            c = colorbar('westoutside','FontSize',6.5);
            c.Position = c.Position+[0.832 0.056 -0.02 -0.065];
            c.Ticks = [0 1];
            
            colormap(coolwarm());
            
            
            
            nexttile([2 3])
            set(gca,'Visible','off');
            
            nexttile([2 27])
            
            temp_1 = temp_EV;
            
            temp_x = 1:length(temp_1);
            
            [temp_x_min,temp_x_max] = bounds(temp_x);
            [temp_y_min,temp_y_max] = bounds(temp_1);
            
            temp_width = 0.8;
            temp_x_min = temp_x_min - temp_width/2;
            temp_x_max = temp_x_max + temp_width/2;
            
            bar(temp_x, temp_1,temp_width, ...
                'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0],'LineWidth',1);
            hold on
            
            
            xlim([temp_x_min temp_x_max]);
            ylim([max(0,temp_y_min-(temp_y_max-temp_y_min)*0.25) temp_y_max+(temp_y_max-temp_y_min)*0.25]);
            
            set(gca,'xticklabel','');
            h=gca;
            h.XAxis.TickLength = [0 0];
            h.YAxis.TickLength = [0 0];
            
            set(gca,'linewidth',1.5)
            set(gca,'box','off');% 取消右、上边框
            set(gca, 'FontSize', 8) %14
            ylabel('EV (%)', 'FontSize', 9, 'FontWeight', 'normal');
            
        end
        
        
        
        %% Correlations between features and PCs of example FOV (baseline)
        if if_plot == 1
            fig = figure('Name','Correlations between features and PCs of example FOV (baseline)','NumberTitle','off');
            %set(gcf,'Position',[10 50 390 103*1.9]);
            set(gcf,'Position',[10 450 390 103*1.9*0.85*1.5*0.85]);
            
            t = tiledlayout(15,31,'TileSpacing','none','Padding','compact');
            
            nexttile([3 4])
            set(gca,'Visible','off');
            
            nexttile([3 26])
            
            temp_EV = temp_pca_explained_baseline_example(1:temp_PCLimit);
            
            C = temp_pca_r_eachPC_VS_eventVector_baseline_example(1:temp_PCLimit,1:3)';
            %C_max = 1;
            C = abs(C);
            %C_max = max(C,[],'all');
            %C_min = min(C,[],'all');
            C_max = 1;
            C_min = 0;
            
            
            imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
            hold on
            
            set(gca,'linewidth',1.5)
            set(gca,'box','off');% 取消右、上边框
            set(gca, 'FontSize', 8) %12
            
            
            %set(gca,'XTick',1:size(C,2));
            set(gca,'xticklabel','');
            
            set(gca,'YTick',1:size(C,1));
            ytl = ["Prior","WM strength","Meta"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
            ytext_xp=(xt(1))*ones(1,length(yt))-4.82;%-1.6
            ytext_yp=yt;
            text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',8);%8,7
            set(gca,'yticklabel','');
            
            set(gca,'TickLength',[0 0]);
            
            x_lim = xlim;
            y_lim = ylim;
            %xlabel(sprintf('PC (1~25, EV = %.0f%%)',sum(temp_EV)),'Position',[mean(x_lim) max(y_lim)+0.2],'FontSize', 9);
            temp_title = title(sprintf('Correlations with PCs of one FOV (baseline)'),'FontSize',9);%11
            
            c = colorbar('westoutside','FontSize',6.5);
            c.Position = c.Position+[0.822 0.016 -0.02 -0.025];
            c.Ticks = [0 1];
            
            colormap(coolwarm());
            
            
            nexttile([1 31])
            set(gca,'Visible','off');
            
            nexttile([6 4])
            set(gca,'Visible','off');
            
            nexttile([6 26])
            
            temp_EV = temp_pca_explained_baseline_example(1:temp_PCLimit);
            
            C = temp_pca_r_eachPC_VS_locDistri_baseline_example(1:temp_PCLimit,:)';
            %C_max = 1;
            C = abs(C);
            %C_max = max(C,[],'all');
            %C_min = min(C,[],'all');
            C_max = 1;
            C_min = 0;
            
            
            imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
            hold on
            
            set(gca,'linewidth',1.5)
            set(gca,'box','off');% 取消右、上边框
            set(gca, 'FontSize', 8) %12
            
            
            %set(gca,'XTick',1:size(C,2));
            set(gca,'xticklabel','');
            
            set(gca,'YTick',1:size(C,1));
            ytl = ["1","2","3","4","5","6"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
            ytext_xp=(xt(1))*ones(1,length(yt))-4.8;%-1.6
            ytext_yp=yt;
            text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',7.5);%8,7
            set(gca,'yticklabel','');
            
            set(gca,'TickLength',[0 0]);
            
            x_lim = xlim;
            y_lim = ylim;
            %xlabel(sprintf('PC (1~25, EV = %.0f%%)',sum(temp_EV)),'Position',[mean(x_lim) max(y_lim)+0.2],'FontSize', 9);
            ylabel('Location','FontSize',9,'Position',[min(x_lim)-1 mean(y_lim)]);
            %temp_title = title(sprintf('Correlations with PCs of one FOV (baseline)'),'FontSize',9);%11
            
            %c = colorbar('westoutside','FontSize',6.5);
            %c.Position = c.Position+[0.832 0.056 -0.02 -0.065];
            %c.Ticks = [0 1];
            
            colormap(coolwarm());
            
            
            nexttile([1 31])
            set(gca,'Visible','off');
            
            
            nexttile([4 4])
            set(gca,'Visible','off');
            
            nexttile([4 26])
            
            temp_1 = temp_EV;
            
            temp_x = 1:length(temp_1);
            
            [temp_x_min,temp_x_max] = bounds(temp_x);
            [temp_y_min,temp_y_max] = bounds(temp_1);
            
            temp_width = 0.8;
            temp_x_min = temp_x_min - temp_width/2;
            temp_x_max = temp_x_max + temp_width/2;
            
            bar(temp_x, temp_1,temp_width, ...
                'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0],'LineWidth',1);
            hold on
            
            
            xlim([temp_x_min temp_x_max]);
            ylim([max(0,temp_y_min-(temp_y_max-temp_y_min)*0.25) temp_y_max+(temp_y_max-temp_y_min)*0.25]);
            
            set(gca,'xticklabel','');
            h=gca;
            h.XAxis.TickLength = [0 0];
            h.YAxis.TickLength = [0 0];
            
            set(gca,'linewidth',1.5)
            set(gca,'box','off');% 取消右、上边框
            set(gca, 'FontSize', 8) %14
            ylabel('EV (%)', 'FontSize', 9, 'FontWeight', 'normal');
            
            x_lim = xlim;
            y_lim = ylim;
            
            xlabel(sprintf('PC (1~25, EV = %.0f%%)',sum(temp_EV)),'Position',[mean(x_lim) min(y_lim)-0.5],'FontSize', 9);
            
        end
        
        
        
        %% Correlations between features and PCs of example FOV (delay)
        if if_plot == 1
            fig = figure('Name','Correlations between features and PCs of example FOV (delay)','NumberTitle','off');
            set(gcf,'Position',[450 50 390 103*1.9*0.85]);
            
            t = tiledlayout(4,31,'TileSpacing','tight','Padding','compact');
            
            nexttile([2 3])
            set(gca,'Visible','off');
            
            nexttile([2 27])
            
            temp_EV = temp_pca_explained_delay_example(1:temp_PCLimit);
            
            C = temp_pca_r_eachPC_VS_eventVector_delay_example(1:temp_PCLimit,:)';
            %C_max = 1;
            C = abs(C);
            %C_max = max(C,[],'all');
            %C_min = min(C,[],'all');
            
            C_max = 1;
            C_min = 0;
            
            imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
            hold on
            
            set(gca,'linewidth',1.5)
            set(gca,'box','off');% 取消右、上边框
            set(gca, 'FontSize', 8) %12
            
            
            %set(gca,'XTick',1:size(C,2));
            set(gca,'xticklabel','');
            
            set(gca,'YTick',1:size(C,1));
            ytl = ["Prior","WM strength","Meta","Time"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
            ytext_xp=(xt(1))*ones(1,length(yt))-4.8;%-1.6
            ytext_yp=yt;
            text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',7.5);%8,7
            set(gca,'yticklabel','');
            
            set(gca,'TickLength',[0 0]);
            
            x_lim = xlim;
            y_lim = ylim;
            xlabel(sprintf('PC (1~25, EV = %.0f%%)',sum(temp_EV)),'Position',[mean(x_lim) max(y_lim)+0.2],'FontSize', 9);
            temp_title = title(sprintf('Correlations with PCs of one FOV (delay)'),'FontSize',9);%11
            
            c = colorbar('westoutside','FontSize',6.5);
            c.Position = c.Position+[0.832 0.056 -0.02 -0.065];
            c.Ticks = [0 1];
            
            colormap(coolwarm());
            
            
            
            nexttile([2 3])
            set(gca,'Visible','off');
            
            nexttile([2 27])
            
            temp_1 = temp_EV;
            
            temp_x = 1:length(temp_1);
            
            [temp_x_min,temp_x_max] = bounds(temp_x);
            [temp_y_min,temp_y_max] = bounds(temp_1);
            
            temp_width = 0.8;
            temp_x_min = temp_x_min - temp_width/2;
            temp_x_max = temp_x_max + temp_width/2;
            
            bar(temp_x, temp_1,temp_width, ...
                'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0],'LineWidth',1);
            hold on
            
            
            xlim([temp_x_min temp_x_max]);
            ylim([max(0,temp_y_min-(temp_y_max-temp_y_min)*0.25) temp_y_max+(temp_y_max-temp_y_min)*0.25]);
            
            set(gca,'xticklabel','');
            h=gca;
            h.XAxis.TickLength = [0 0];
            h.YAxis.TickLength = [0 0];
            
            set(gca,'linewidth',1.5)
            set(gca,'box','off');% 取消右、上边框
            set(gca, 'FontSize', 8) %14
            ylabel('EV (%)', 'FontSize', 9, 'FontWeight', 'normal');
            
        end
        
        
        %% Correlations between features and PCs of example FOV (delay)
        if if_plot == 1
            fig = figure('Name','Correlations between features and PCs of example FOV (delay)','NumberTitle','off');
            %set(gcf,'Position',[10 50 390 103*1.9]);
            set(gcf,'Position',[450 450 390 103*1.9*0.85*1.5*0.85]);
            
            t = tiledlayout(15,31,'TileSpacing','none','Padding','compact');
            
            nexttile([3 4])
            set(gca,'Visible','off');
            
            nexttile([3 26])
            
            temp_EV = temp_pca_explained_delay_example(1:temp_PCLimit);
            
            C = temp_pca_r_eachPC_VS_eventVector_delay_example(1:temp_PCLimit,1:3)';
            %C_max = 1;
            C = abs(C);
            %C_max = max(C,[],'all');
            %C_min = min(C,[],'all');
            C_max = 1;
            C_min = 0;
            
            
            imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
            hold on
            
            set(gca,'linewidth',1.5)
            set(gca,'box','off');% 取消右、上边框
            set(gca, 'FontSize', 8) %12
            
            
            %set(gca,'XTick',1:size(C,2));
            set(gca,'xticklabel','');
            
            set(gca,'YTick',1:size(C,1));
            ytl = ["Prior","WM strength","Meta"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
            ytext_xp=(xt(1))*ones(1,length(yt))-4.82;%-1.6
            ytext_yp=yt;
            text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',8);%8,7
            set(gca,'yticklabel','');
            
            set(gca,'TickLength',[0 0]);
            
            x_lim = xlim;
            y_lim = ylim;
            %xlabel(sprintf('PC (1~25, EV = %.0f%%)',sum(temp_EV)),'Position',[mean(x_lim) max(y_lim)+0.2],'FontSize', 9);
            temp_title = title(sprintf('Correlations with PCs of one FOV (delay)'),'FontSize',9);%11
            
            c = colorbar('westoutside','FontSize',6.5);
            c.Position = c.Position+[0.822 0.016 -0.02 -0.025];
            c.Ticks = [0 1];
            
            colormap(coolwarm());
            
            
            nexttile([1 31])
            set(gca,'Visible','off');
            
            nexttile([6 4])
            set(gca,'Visible','off');
            
            nexttile([6 26])
            
            temp_EV = temp_pca_explained_delay_example(1:temp_PCLimit);
            
            C = temp_pca_r_eachPC_VS_locDistri_delay_example(1:temp_PCLimit,:)';
            %C_max = 1;
            C = abs(C);
            %C_max = max(C,[],'all');
            %C_min = min(C,[],'all');
            C_max = 1;
            C_min = 0;
            
            
            imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
            hold on
            
            set(gca,'linewidth',1.5)
            set(gca,'box','off');% 取消右、上边框
            set(gca, 'FontSize', 8) %12
            
            
            %set(gca,'XTick',1:size(C,2));
            set(gca,'xticklabel','');
            
            set(gca,'YTick',1:size(C,1));
            ytl = ["1","2","3","4","5","6"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
            ytext_xp=(xt(1))*ones(1,length(yt))-4.8;%-1.6
            ytext_yp=yt;
            text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',7.5);%8,7
            set(gca,'yticklabel','');
            
            set(gca,'TickLength',[0 0]);
            
            x_lim = xlim;
            y_lim = ylim;
            %xlabel(sprintf('PC (1~25, EV = %.0f%%)',sum(temp_EV)),'Position',[mean(x_lim) max(y_lim)+0.2],'FontSize', 9);
            ylabel('Location','FontSize',9,'Position',[min(x_lim)-1 mean(y_lim)]);
            %temp_title = title(sprintf('Correlations with PCs of one FOV (baseline)'),'FontSize',9);%11
            
            %c = colorbar('westoutside','FontSize',6.5);
            %c.Position = c.Position+[0.832 0.056 -0.02 -0.065];
            %c.Ticks = [0 1];
            
            colormap(coolwarm());
            
            
            nexttile([1 31])
            set(gca,'Visible','off');
            
            
            nexttile([4 4])
            set(gca,'Visible','off');
            
            nexttile([4 26])
            
            temp_1 = temp_EV;
            
            temp_x = 1:length(temp_1);
            
            [temp_x_min,temp_x_max] = bounds(temp_x);
            [temp_y_min,temp_y_max] = bounds(temp_1);
            
            temp_width = 0.8;
            temp_x_min = temp_x_min - temp_width/2;
            temp_x_max = temp_x_max + temp_width/2;
            
            bar(temp_x, temp_1,temp_width, ...
                'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0],'LineWidth',1);
            hold on
            
            
            xlim([temp_x_min temp_x_max]);
            ylim([max(0,temp_y_min-(temp_y_max-temp_y_min)*0.25) temp_y_max+(temp_y_max-temp_y_min)*0.25]);
            
            set(gca,'xticklabel','');
            h=gca;
            h.XAxis.TickLength = [0 0];
            h.YAxis.TickLength = [0 0];
            
            set(gca,'linewidth',1.5)
            set(gca,'box','off');% 取消右、上边框
            set(gca, 'FontSize', 8) %14
            ylabel('EV (%)', 'FontSize', 9, 'FontWeight', 'normal');
            
            x_lim = xlim;
            y_lim = ylim;
            
            xlabel(sprintf('PC (1~25, EV = %.0f%%)',sum(temp_EV)),'Position',[mean(x_lim) min(y_lim)-0.5],'FontSize', 9);
            
        end
        
        
        %%
        
        if if_monkey_D0_Z1 == 0 && currentSessionIndex_AB == 8
            seqSet_inOne_inOne;
            
            
            temp_PCIndex_location = 5;%5
            temp_PCIndex_WMStrength = 4;
            temp_PCIndex_meta = 6;
            
            exampleSeq = 7;%16,7,15
            
            %[location 5,strength 4,seq 7]
            %[location 5,meta 6,seq 7]
            
            if true
                temp_dff_baseline = temp_F_dff_baselineBin(:,tempBoolIndex_valid)';
                temp_dff_delay = temp_F_dff_decisionBin1(:,tempBoolIndex_valid)';
                
                score_delay_hat = (temp_dff_delay-mu_delay)*coeff_delay;
                temp_error = score_delay_hat-score_delay;
                temp_errorSum = sum(temp_error,'all');
                
                temp_dff_delay_mean = mean(temp_dff_delay,1);
                
                score_baselineIntoDelay_proj = (temp_dff_baseline-mu_baseline)*coeff_delay;
            end
            
            temptempBoolIndex = seqIndex' == exampleSeq;
            temptempBoolIndex_exampleSeq = temptempBoolIndex(tempBoolIndex_valid);
            tempExampleSeq_memoryPrecision = temp_memoryPrecision(temptempBoolIndex_exampleSeq);
            
            temp_lowThreshold_memoryPrecision = lowThreshold_memoryPrecision;
            %temp_lowThreshold_memoryPrecision = median(tempExampleSeq_memoryPrecision);
            
            tempBoolIndex_lowMemoryPrecision = temp_memoryPrecision < temp_lowThreshold_memoryPrecision;
            
            
            temptempBoolIndex_lowMemoryPrecision_exampleSeq = temptempBoolIndex_exampleSeq & tempBoolIndex_lowMemoryPrecision;
            temptempBoolIndex_highMemoryPrecision_exampleSeq = temptempBoolIndex_exampleSeq & (~tempBoolIndex_lowMemoryPrecision);
            
            trialA_locationPC = score_delay(temptempBoolIndex_lowMemoryPrecision_exampleSeq,temp_PCIndex_location);
            trialB_locationPC = score_delay(temptempBoolIndex_highMemoryPrecision_exampleSeq,temp_PCIndex_location);
            trialA_WMStrengthPC = score_delay(temptempBoolIndex_lowMemoryPrecision_exampleSeq,temp_PCIndex_WMStrength);
            trialB_WMStrengthPC = score_delay(temptempBoolIndex_highMemoryPrecision_exampleSeq,temp_PCIndex_WMStrength);
            
            trialA_locationPC_baselineProj = score_baselineIntoDelay_proj(temptempBoolIndex_lowMemoryPrecision_exampleSeq,temp_PCIndex_location);
            trialB_locationPC_baselineProj = score_baselineIntoDelay_proj(temptempBoolIndex_highMemoryPrecision_exampleSeq,temp_PCIndex_location);
            trialA_WMStrengthPC_baselineProj = score_baselineIntoDelay_proj(temptempBoolIndex_lowMemoryPrecision_exampleSeq,temp_PCIndex_WMStrength);
            trialB_WMStrengthPC_baselineProj = score_baselineIntoDelay_proj(temptempBoolIndex_highMemoryPrecision_exampleSeq,temp_PCIndex_WMStrength);
            
            [~,temp_p_AB_locationPC] = ttest2(trialA_locationPC,trialB_locationPC);
            [~,temp_p_AB_WMStrengthPC] = ttest2(trialA_WMStrengthPC,trialB_WMStrengthPC);
            
            [~,temp_p_AprojBproj_locationPC] = ttest2(trialA_locationPC_baselineProj,trialB_locationPC_baselineProj);
            [~,temp_p_AprojBproj_WMStrengthPC] = ttest2(trialA_WMStrengthPC_baselineProj,trialB_WMStrengthPC_baselineProj);
            
            fprintf('Seq [%d],temp_p_AB_locationPC=%.3f,temp_p_AB_WMStrengthPC=%.3f.\n',...
                seqSet_inOne_inOne(exampleSeq),temp_p_AB_locationPC,temp_p_AB_WMStrengthPC);
            
            
            
            
            temp_lowThreshold_meta = lowThreshold_meta;
            
            tempBoolIndex_lowMeta = temp_meta < temp_lowThreshold_meta;
            
            
            temptempBoolIndex_lowMeta_exampleSeq = temptempBoolIndex_exampleSeq & tempBoolIndex_lowMeta;
            temptempBoolIndex_highMeta_exampleSeq = temptempBoolIndex_exampleSeq & (~tempBoolIndex_lowMeta);
            
            trialA2_locationPC = score_delay(temptempBoolIndex_lowMeta_exampleSeq,temp_PCIndex_location);
            trialB2_locationPC = score_delay(temptempBoolIndex_highMeta_exampleSeq,temp_PCIndex_location);
            trialA2_metaPC = score_delay(temptempBoolIndex_lowMeta_exampleSeq,temp_PCIndex_meta);
            trialB2_metaPC = score_delay(temptempBoolIndex_highMeta_exampleSeq,temp_PCIndex_meta);
            
            trialA2_locationPC_baselineProj = score_baselineIntoDelay_proj(temptempBoolIndex_lowMeta_exampleSeq,temp_PCIndex_location);
            trialB2_locationPC_baselineProj = score_baselineIntoDelay_proj(temptempBoolIndex_highMeta_exampleSeq,temp_PCIndex_location);
            trialA2_metaPC_baselineProj = score_baselineIntoDelay_proj(temptempBoolIndex_lowMeta_exampleSeq,temp_PCIndex_meta);
            trialB2_metaPC_baselineProj = score_baselineIntoDelay_proj(temptempBoolIndex_highMeta_exampleSeq,temp_PCIndex_meta);
            
            
            [~,temp_p_AB2_locationPC] = ttest2(trialA2_locationPC,trialB2_locationPC);
            [~,temp_p_AB2_metaPC] = ttest2(trialA2_metaPC,trialB2_metaPC);
            
            [~,temp_p_AprojBproj2_locationPC] = ttest2(trialA2_locationPC_baselineProj,trialB2_locationPC_baselineProj);
            [~,temp_p_AprojBproj2_metaPC] = ttest2(trialA2_metaPC_baselineProj,trialB2_metaPC_baselineProj);
            
            fprintf('Seq [%d],temp_p_AB2_locationPC=%.3f,temp_p_AB2_metaPC=%.3f.\n',...
                seqSet_inOne_inOne(exampleSeq),temp_p_AB2_locationPC,temp_p_AB2_metaPC);
            
            
            
            trialA_locationPC;
            trialB_locationPC;
            trialA_WMStrengthPC;
            trialB_WMStrengthPC;
            
            trialA_locationPC_baselineProj;
            trialB_locationPC_baselineProj;
            trialA_WMStrengthPC_baselineProj;
            trialB_WMStrengthPC_baselineProj;
            
            trialA_locationPC_mean = mean(trialA_locationPC);
            trialB_locationPC_mean = mean(trialB_locationPC);
            trialA_WMStrengthPC_mean = mean(trialA_WMStrengthPC);
            trialB_WMStrengthPC_mean = mean(trialB_WMStrengthPC);
            trialA_locationPC_baselineProj_mean = mean(trialA_locationPC_baselineProj);
            trialB_locationPC_baselineProj_mean = mean(trialB_locationPC_baselineProj);
            trialA_WMStrengthPC_baselineProj_mean = mean(trialA_WMStrengthPC_baselineProj);
            trialB_WMStrengthPC_baselineProj_mean = mean(trialB_WMStrengthPC_baselineProj);
            
            trialA_locationPC_SEM = std(trialA_locationPC)./sqrt(length(trialA_locationPC));
            trialB_locationPC_SEM = std(trialB_locationPC)./sqrt(length(trialB_locationPC));
            trialA_WMStrengthPC_SEM = std(trialA_WMStrengthPC)./sqrt(length(trialA_WMStrengthPC));
            trialB_WMStrengthPC_SEM = std(trialB_WMStrengthPC)./sqrt(length(trialB_WMStrengthPC));
            trialA_locationPC_baselineProj_SEM = std(trialA_locationPC_baselineProj)./sqrt(length(trialA_locationPC_baselineProj));
            trialB_locationPC_baselineProj_SEM = std(trialB_locationPC_baselineProj)./sqrt(length(trialB_locationPC_baselineProj));
            trialA_WMStrengthPC_baselineProj_SEM = std(trialA_WMStrengthPC_baselineProj)./sqrt(length(trialA_WMStrengthPC_baselineProj));
            trialB_WMStrengthPC_baselineProj_SEM = std(trialB_WMStrengthPC_baselineProj)./sqrt(length(trialB_WMStrengthPC_baselineProj));
            
            
            trialA_2d_mean = [trialA_locationPC_mean,trialA_WMStrengthPC_mean];
            trialB_2d_mean = [trialB_locationPC_mean,trialB_WMStrengthPC_mean];
            trialA_2d_baselineProj_mean = [trialA_locationPC_baselineProj_mean,trialA_WMStrengthPC_baselineProj_mean];
            trialB_2d_baselineProj_mean = [trialB_locationPC_baselineProj_mean,trialB_WMStrengthPC_baselineProj_mean];
            
            trialA_2d_SEM = [trialA_locationPC_SEM,trialA_WMStrengthPC_SEM];
            trialB_2d_SEM = [trialB_locationPC_SEM,trialB_WMStrengthPC_SEM];
            trialA_2d_baselineProj_SEM = [trialA_locationPC_baselineProj_SEM,trialA_WMStrengthPC_baselineProj_SEM];
            trialB_2d_baselineProj_SEM = [trialB_locationPC_baselineProj_SEM,trialB_WMStrengthPC_baselineProj_SEM];
            
            
            %% Plot
            if true
                fig = figure('Name','asd','NumberTitle','off'); %#ok<*NASGU>
                set(gcf,'Position',[450+0 50+0 200 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                
                
                temp1_mean = trialA_2d_mean;
                temp2_mean = trialB_2d_mean;
                temp3_mean = trialA_2d_baselineProj_mean;
                temp4_mean = trialB_2d_baselineProj_mean;
                
                temp1_SEM = trialA_2d_SEM;
                temp2_SEM = trialB_2d_SEM;
                temp3_SEM = trialA_2d_baselineProj_SEM;
                temp4_SEM = trialB_2d_baselineProj_SEM;
                
                x = [temp1_mean(1),temp2_mean(1),temp3_mean(1),temp4_mean(1)];
                y = [temp1_mean(2),temp2_mean(2),temp3_mean(2),temp4_mean(2)];
                
                x_SEM = [temp1_SEM(1),temp2_SEM(1),temp3_SEM(1),temp4_SEM(1)];
                y_SEM = [temp1_SEM(2),temp2_SEM(2),temp3_SEM(2),temp4_SEM(2)];
                
                %plot(x([1,3]),y([1,3]),'color',[0.25 0.25 0.25],'linewidth',1.5);
                errorbar(x([1,3]),y([1,3]),x_SEM([1,3]),'horizontal','color',[0.25 0.25 0.25],'linewidth',1.5);
                hold on
                errorbar(x([1,3]),y([1,3]),y_SEM([1,3]),'color',[0.25 0.25 0.25],'linewidth',1.5);
                hold on
                
                %plot(x([2,4]),y([2,4]),'color',[0.25 0.25 0.25],'linewidth',1.5);
                errorbar(x([2,4]),y([2,4]),x_SEM([2,4]),'horizontal','color',[0.25 0.25 0.25],'linewidth',1.5);
                hold on
                errorbar(x([2,4]),y([2,4]),y_SEM([2,4]),'color',[0.25 0.25 0.25],'linewidth',1.5);
                hold on
                
                [temp_xmin,temp_xmax] = bounds(x);
                [temp_ymin,temp_ymax] = bounds(y);
                
                
                xlim([temp_xmin-(temp_xmax-temp_xmin)*0.4 temp_xmax+(temp_xmax-temp_xmin)*0.4]);
                ylim([temp_ymin-(temp_ymax-temp_ymin)*0.4 temp_ymax+(temp_ymax-temp_ymin)*0.4]);
                
                set(gca,'linewidth',1.5)
                set(gca, 'FontSize', 8)
                set(gca,'box','off');% 取消右、上边框
                
                
                xtickangle(0);
                
                xlabel(sprintf('PC%d (Location)',temp_PCIndex_location), 'FontSize', 9);
                ylabel(sprintf('PC%d (WM strength)',temp_PCIndex_WMStrength), 'FontSize', 9);
                
                
                %temp_title = title(sprintf('Across seqs'), 'FontSize', 9);
                
            end
            
            
            trialA2_locationPC;
            trialB2_locationPC;
            trialA2_metaPC;
            trialB2_metaPC;
            
            trialA2_locationPC_baselineProj;
            trialB2_locationPC_baselineProj;
            trialA2_metaPC_baselineProj;
            trialB2_metaPC_baselineProj;
            
            trialA2_locationPC_mean = mean(trialA2_locationPC);
            trialB2_locationPC_mean = mean(trialB2_locationPC);
            trialA2_metaPC_mean = mean(trialA2_metaPC);
            trialB2_metaPC_mean = mean(trialB2_metaPC);
            trialA2_locationPC_baselineProj_mean = mean(trialA2_locationPC_baselineProj);
            trialB2_locationPC_baselineProj_mean = mean(trialB2_locationPC_baselineProj);
            trialA2_metaPC_baselineProj_mean = mean(trialA2_metaPC_baselineProj);
            trialB2_metaPC_baselineProj_mean = mean(trialB2_metaPC_baselineProj);
            
            trialA2_locationPC_SEM = std(trialA2_locationPC)./sqrt(length(trialA2_locationPC));
            trialB2_locationPC_SEM = std(trialB2_locationPC)./sqrt(length(trialB2_locationPC));
            trialA2_metaPC_SEM = std(trialA2_metaPC)./sqrt(length(trialA2_metaPC));
            trialB2_metaPC_SEM = std(trialB2_metaPC)./sqrt(length(trialB2_metaPC));
            trialA2_locationPC_baselineProj_SEM = std(trialA2_locationPC_baselineProj)./sqrt(length(trialA2_locationPC_baselineProj));
            trialB2_locationPC_baselineProj_SEM = std(trialB2_locationPC_baselineProj)./sqrt(length(trialB2_locationPC_baselineProj));
            trialA2_metaPC_baselineProj_SEM = std(trialA2_metaPC_baselineProj)./sqrt(length(trialA2_metaPC_baselineProj));
            trialB2_metaPC_baselineProj_SEM = std(trialB2_metaPC_baselineProj)./sqrt(length(trialB2_metaPC_baselineProj));
            
            
            trialA2_2d_mean = [trialA2_locationPC_mean,trialA2_metaPC_mean];
            trialB2_2d_mean = [trialB2_locationPC_mean,trialB2_metaPC_mean];
            trialA2_2d_baselineProj_mean = [trialA2_locationPC_baselineProj_mean,trialA2_metaPC_baselineProj_mean];
            trialB2_2d_baselineProj_mean = [trialB2_locationPC_baselineProj_mean,trialB2_metaPC_baselineProj_mean];
            
            trialA2_2d_SEM = [trialA2_locationPC_SEM,trialA2_metaPC_SEM];
            trialB2_2d_SEM = [trialB2_locationPC_SEM,trialB2_metaPC_SEM];
            trialA2_2d_baselineProj_SEM = [trialA2_locationPC_baselineProj_SEM,trialA2_metaPC_baselineProj_SEM];
            trialB2_2d_baselineProj_SEM = [trialB2_locationPC_baselineProj_SEM,trialB2_metaPC_baselineProj_SEM];
            
            
            
            %% Plot
            if true
                fig = figure('Name','asd','NumberTitle','off'); %#ok<*NASGU>
                set(gcf,'Position',[700+0 50+0 200 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                
                
                temp1_mean = trialA2_2d_mean;
                temp2_mean = trialB2_2d_mean;
                temp3_mean = trialA2_2d_baselineProj_mean;
                temp4_mean = trialB2_2d_baselineProj_mean;
                
                temp1_SEM = trialA2_2d_SEM;
                temp2_SEM = trialB2_2d_SEM;
                temp3_SEM = trialA2_2d_baselineProj_SEM;
                temp4_SEM = trialB2_2d_baselineProj_SEM;
                
                x = [temp1_mean(1),temp2_mean(1),temp3_mean(1),temp4_mean(1)];
                y = [temp1_mean(2),temp2_mean(2),temp3_mean(2),temp4_mean(2)];
                
                x_SEM = [temp1_SEM(1),temp2_SEM(1),temp3_SEM(1),temp4_SEM(1)];
                y_SEM = [temp1_SEM(2),temp2_SEM(2),temp3_SEM(2),temp4_SEM(2)];
                
                %plot(x([1,3]),y([1,3]),'color',[0.25 0.25 0.25],'linewidth',1.5);
                errorbar(x([1,3]),y([1,3]),x_SEM([1,3]),'horizontal','color',[0.25 0.25 0.25],'linewidth',1.5);
                hold on
                errorbar(x([1,3]),y([1,3]),y_SEM([1,3]),'color',[0.25 0.25 0.25],'linewidth',1.5);
                hold on
                
                %plot(x([2,4]),y([2,4]),'color',[0.25 0.25 0.25],'linewidth',1.5);
                errorbar(x([2,4]),y([2,4]),x_SEM([2,4]),'horizontal','color',[0.25 0.25 0.25],'linewidth',1.5);
                hold on
                errorbar(x([2,4]),y([2,4]),y_SEM([2,4]),'color',[0.25 0.25 0.25],'linewidth',1.5);
                hold on
                
                [temp_xmin,temp_xmax] = bounds(x);
                [temp_ymin,temp_ymax] = bounds(y);
                
                
                xlim([temp_xmin-(temp_xmax-temp_xmin)*0.4 temp_xmax+(temp_xmax-temp_xmin)*0.4]);
                ylim([temp_ymin-(temp_ymax-temp_ymin)*0.4 temp_ymax+(temp_ymax-temp_ymin)*0.4]);
                
                set(gca,'linewidth',1.5)
                set(gca, 'FontSize', 8)
                set(gca,'box','off');% 取消右、上边框
                
                
                xtickangle(0);
                
                xlabel(sprintf('PC%d (Location)',temp_PCIndex_location), 'FontSize', 9);
                ylabel(sprintf('PC%d (Meta)',temp_PCIndex_meta), 'FontSize', 9);
                
                
                %temp_title = title(sprintf('Across seqs'), 'FontSize', 9);
                
            end
            
        end
    end
    
    
    %% Get selective neuron from subspace view
    if true
        
        beta_prior_linearAxis_halfA_resampled;
        
        temp_minusFold = 5;%5
        temp_resampleRatio = 0.99;%0.95,0.99
        
        temp_resampleNum = length(beta_prior_linearAxis_halfA_resampled);
        temp_roiNum = length(beta_prior_linearAxis_halfA_resampled{1});
        
        
        %% tempSelectiveBoolIndex_priorSubspace
        betaA = beta_prior_linearAxis_halfA_resampled;
        betaB = beta_prior_linearAxis_halfB_resampled;
        
        betaABMean_acrossReample_2d = nan(temp_roiNum,temp_resampleNum);
        betaABMinus_acrossReample_2d = nan(temp_roiNum,temp_resampleNum);
        betaABStd_acrossReample_2d = nan(temp_roiNum,temp_resampleNum);
        for tempi=1:temp_roiNum
            temp_betaA_acrossReample = nan(temp_resampleNum,1);
            temp_betaB_acrossReample = nan(temp_resampleNum,1);
            for tempj=1:temp_resampleNum
                temp_betaA_acrossReample(tempj) = betaA{tempj}(tempi);
                temp_betaB_acrossReample(tempj) = betaB{tempj}(tempi);
            end
            betaABMean_acrossReample_2d(tempi,:) = (temp_betaA_acrossReample+temp_betaB_acrossReample)/2;
            betaABMinus_acrossReample_2d(tempi,:) = (temp_betaA_acrossReample-temp_betaB_acrossReample);
            betaABStd_acrossReample_2d(tempi,:) = std([temp_betaA_acrossReample,temp_betaB_acrossReample],0,2);
        end
        betaABMinus_acrossReampleMean_2d = mean(betaABMinus_acrossReample_2d,2);
        betaABStd_acrossReampleMean_2d = mean(betaABStd_acrossReample_2d,2);
        temp1 = abs(betaABMean_acrossReample_2d) > temp_minusFold*abs(betaABMinus_acrossReampleMean_2d);
        %temp1 = abs(betaABMean_acrossReample_2d) > temp_minusFold*abs(betaABStd_acrossReampleMean_2d);
        temp2 = sum(temp1,2);
        temp3 = temp2 > (temp_resampleNum*temp_resampleRatio);
        
        tempSelectiveBoolIndex_priorSubspace = temp3;
        
        
        %% tempSelectiveBoolIndex_metaSubspace
        betaA = beta_memory_linearAxis_halfA_resampled;
        betaB = beta_memory_linearAxis_halfB_resampled;
        
        betaABMean_acrossReample_2d = nan(temp_roiNum,temp_resampleNum);
        betaABMinus_acrossReample_2d = nan(temp_roiNum,temp_resampleNum);
        betaABStd_acrossReample_2d = nan(temp_roiNum,temp_resampleNum);
        for tempi=1:temp_roiNum
            temp_betaA_acrossReample = nan(temp_resampleNum,1);
            temp_betaB_acrossReample = nan(temp_resampleNum,1);
            for tempj=1:temp_resampleNum
                temp_betaA_acrossReample(tempj) = betaA{tempj}(tempi);
                temp_betaB_acrossReample(tempj) = betaB{tempj}(tempi);
            end
            betaABMean_acrossReample_2d(tempi,:) = (temp_betaA_acrossReample+temp_betaB_acrossReample)/2;
            betaABMinus_acrossReample_2d(tempi,:) = (temp_betaA_acrossReample-temp_betaB_acrossReample);
            betaABStd_acrossReample_2d(tempi,:) = std([temp_betaA_acrossReample,temp_betaB_acrossReample],0,2);
        end
        betaABMinus_acrossReampleMean_2d = mean(betaABMinus_acrossReample_2d,2);
        betaABStd_acrossReampleMean_2d = mean(betaABStd_acrossReample_2d,2);
        temp1 = abs(betaABMean_acrossReample_2d) > temp_minusFold*abs(betaABMinus_acrossReampleMean_2d);
        %temp1 = abs(betaABMean_acrossReample_2d) > temp_minusFold*abs(betaABStd_acrossReampleMean_2d);
        temp2 = sum(temp1,2);
        temp3 = temp2 > (temp_resampleNum*temp_resampleRatio);
        
        tempSelectiveBoolIndex_memorySubspace = temp3;
        
        
        %% tempSelectiveBoolIndex_metaSubspace
        betaA = beta_meta_linearAxis_halfA_resampled;
        betaB = beta_meta_linearAxis_halfB_resampled;
        
        betaABMean_acrossReample_2d = nan(temp_roiNum,temp_resampleNum);
        betaABMinus_acrossReample_2d = nan(temp_roiNum,temp_resampleNum);
        betaABStd_acrossReample_2d = nan(temp_roiNum,temp_resampleNum);
        for tempi=1:temp_roiNum
            temp_betaA_acrossReample = nan(temp_resampleNum,1);
            temp_betaB_acrossReample = nan(temp_resampleNum,1);
            for tempj=1:temp_resampleNum
                temp_betaA_acrossReample(tempj) = betaA{tempj}(tempi);
                temp_betaB_acrossReample(tempj) = betaB{tempj}(tempi);
            end
            betaABMean_acrossReample_2d(tempi,:) = (temp_betaA_acrossReample+temp_betaB_acrossReample)/2;
            betaABMinus_acrossReample_2d(tempi,:) = (temp_betaA_acrossReample-temp_betaB_acrossReample);
            betaABStd_acrossReample_2d(tempi,:) = std([temp_betaA_acrossReample,temp_betaB_acrossReample],0,2);
        end
        betaABMinus_acrossReampleMean_2d = mean(betaABMinus_acrossReample_2d,2);
        betaABStd_acrossReampleMean_2d = mean(betaABStd_acrossReample_2d,2);
        temp1 = abs(betaABMean_acrossReample_2d) > temp_minusFold*abs(betaABMinus_acrossReampleMean_2d);
        %temp1 = abs(betaABMean_acrossReample_2d) > temp_minusFold*abs(betaABStd_acrossReampleMean_2d);
        temp2 = sum(temp1,2);
        temp3 = temp2 > (temp_resampleNum*temp_resampleRatio);
        
        tempSelectiveBoolIndex_metaSubspace = temp3;
        
        %% Other
        tempSelectiveBoolIndex_priorSubspace;
        tempSelectiveBoolIndex_memorySubspace;
        tempSelectiveBoolIndex_metaSubspace;
        
        fprintf('Subspace-based selective number [Prior,WM,Meta] = [%d,%d,%d]. (minusFold=%d,resampleRatio=%.3f)\n',...
            sum(tempSelectiveBoolIndex_priorSubspace),sum(tempSelectiveBoolIndex_memorySubspace),sum(tempSelectiveBoolIndex_metaSubspace),...
            temp_minusFold,temp_resampleRatio);
        
        %     betaABMinus_acrossReampleMean_2d;
        %     betaABMinus_acrossReample_2d;
        %     temp_betaA_acrossReample;
        %     temp_betaB_acrossReample;
        %     b1 = 1;
        %     b2 = 2;
        %     %b12Mean = (b1+b2)/2;
        %     %b12Std = sqrt(((b1-b12Mean).^2 + (b2-b12Mean).^2)/(2-1));
        %     b12Mean = mean([b1,b2]);
        %     b12Std = std([b1,b2]);
        
        
    end
    
    
    %% Dynamics from subspace view
    if true
        
        %% Prepare F_dff
        if_delay1_forward0_backward1 = 1;
        if exist('if_compute_summary','var') == 1
            if_delay1_forward0_backward1 = 1;
        end
        
        if if_monkey_D0_Z1 == 0 && currentSessionIndex_AB == 8
            if_delay1_forward0_backward1 = 0;
        end
        
        if if_delay1_forward0_backward1 == 1
            F_dff_delay = F_dff_decisionPeriodA;
        elseif if_delay1_forward0_backward1 == 0
            F_dff_delay = F_dff_decisionPeriodA;
            F_dff_delay(:,:,1:decisionPeriodB_interval(2)) = F_dff_decisionPeriodB;
        end
        F_dff_delay = F_dff_delay(:,:,1:34);
        
        
        F_dff_lastT = nan(size(F_dff_delay,1),size(F_dff_delay,2),18);
        temp_seqLength = nan(size(F_dff_delay,2),1);
        for tempi=1:size(F_dff_delay,2)
            temp_seqLength(tempi) = sum(boolIndex_location_seq_T(seqIndex(tempi),:));
        end
        F_dff_lastT(:,temp_seqLength==1,:) = F_dff_length1_sample(:,temp_seqLength==1,end-17:end);
        F_dff_lastT(:,temp_seqLength==2,:) = F_dff_length2_sample(:,temp_seqLength==2,end-17:end);
        F_dff_lastT(:,temp_seqLength==3,:) = F_dff_length3_sample(:,temp_seqLength==3,end-17:end);
        
        F_dff_lastTDelay = cat(3,F_dff_lastT,F_dff_delay);
        temp_F_dff_lastTDelay = F_dff_lastTDelay;
        
        %temptemp_range = (baselinePeriod_interval(2)-8):baselinePeriod_interval(3);
        temptemp_range = baselinePeriod_interval(1):baselinePeriod_interval(3);
        temp_F_dff_baseline = F_dff_baselinePeriod(:,:,temptemp_range);
        
        temp_F_dff_baseline_raw = temp_F_dff_baseline;
        temp_F_dff_lastTDelay_raw = temp_F_dff_lastTDelay;
        if temp_if_zscore == 1
            dff_mean_baselineBin = nan(temp_roiNum,1);
            dff_std_baselineBin = nan(temp_roiNum,1);
            dff_mean_delayBin = nan(temp_roiNum,1);
            dff_std_delayBin = nan(temp_roiNum,1);
            for tempi=1:temp_roiNum
                temp1 = F_dff_baselineBin(tempi,:);
                dff_mean_baselineBin(tempi) = mean(temp1);
                dff_std_baselineBin(tempi) = std(temp1);
                
                temp2 = F_dff_decisionBin1(tempi,:);
                dff_mean_delayBin(tempi) = mean(temp2);
                dff_std_delayBin(tempi) = std(temp2);
            end
            temp_F_dff_baseline = (temp_F_dff_baseline_raw-dff_mean_baselineBin)./dff_std_baselineBin;% z-scored
            temp_F_dff_baseline_Zdelay = (temp_F_dff_baseline_raw-dff_mean_delayBin)./dff_std_delayBin;% z-scored
            
            temp_F_dff_lastTDelay_Zbaseline = (temp_F_dff_lastTDelay_raw-dff_mean_baselineBin)./dff_std_baselineBin;% z-scored
            temp_F_dff_lastTDelay = (temp_F_dff_lastTDelay_raw-dff_mean_delayBin)./dff_std_delayBin;% z-scored
            
        elseif temp_if_zscore == 0
            temp_F_dff_baseline_Zdelay = temp_F_dff_baseline_raw;
            temp_F_dff_lastTDelay_Zbaseline = temp_F_dff_lastTDelay_raw;
        end
        temp_F_dff_baseline;
        temp_F_dff_baseline_Zdelay;
        temp_F_dff_lastTDelay_Zbaseline;
        temp_F_dff_lastTDelay;
        %temp1 = mean(temp_F_dff_baseline,3);
        
        
        %% Proj to beta_prior_linearAxis
        temp_x1 = temp_F_dff_baseline;
        temp_x2 = temp_F_dff_lastTDelay_Zbaseline;
        temp_y = temp_meta_trialLevel_baseline;
        temp_beta = beta_prior_linearAxis;
        temp_beta0 = beta0_prior_linearAxis;
        
        temp_y1_hat = squeeze(pagemtimes(temp_beta',temp_x1) + temp_beta0);
        temp_y1_hat_raw = temp_y1_hat;
        temp_y2_hat = squeeze(pagemtimes(temp_beta',temp_x2) + temp_beta0);
        temp_y2_hat_raw = temp_y2_hat;
        if temp_if_zscore == 1
            temp_y_mean = mean(temp_y,'omitnan');
            temp_y_std = std(temp_y,'omitnan');
            temp_y1_hat = (temp_y1_hat_raw.*temp_y_std)+temp_y_mean;
            temp_y2_hat = (temp_y2_hat_raw.*temp_y_std)+temp_y_mean;
        end
        temp_y1_hat_mean = mean(temp_y1_hat,2);
        temp_y2_hat_mean = mean(temp_y2_hat,2);
        
        temp_r1 = corr(temp_y(~isnan(temp_y)),temp_y1_hat_mean(~isnan(temp_y)));
        temp_r2 = corr(temp_y(~isnan(temp_y)),temp_y2_hat_mean(~isnan(temp_y)));
        
        prior1_dynamicProj = temp_y1_hat;
        prior2_dynamicProj = temp_y2_hat;
        
        %% Proj to beta_memory_linearAxis
        temp_x1 = temp_F_dff_baseline_Zdelay;
        temp_x2 = temp_F_dff_lastTDelay;
        temp_y = memoryPrecision_trialLevel;
        temp_beta = beta_memory_linearAxis;
        temp_beta0 = beta0_memory_linearAxis;
        
        temp_y1_hat = squeeze(pagemtimes(temp_beta',temp_x1) + temp_beta0);
        temp_y1_hat_raw = temp_y1_hat;
        temp_y2_hat = squeeze(pagemtimes(temp_beta',temp_x2) + temp_beta0);
        temp_y2_hat_raw = temp_y2_hat;
        if temp_if_zscore == 1
            temp_y_mean = mean(temp_y,'omitnan');
            temp_y_std = std(temp_y,'omitnan');
            temp_y1_hat = (temp_y1_hat_raw.*temp_y_std)+temp_y_mean;
            temp_y2_hat = (temp_y2_hat_raw.*temp_y_std)+temp_y_mean;
        end
        temp_y1_hat_mean = mean(temp_y1_hat,2);
        temp_y2_hat_mean = mean(temp_y2_hat,2);
        
        temp_r1 = corr(temp_y(~isnan(temp_y)),temp_y1_hat_mean(~isnan(temp_y)));
        temp_r2 = corr(temp_y(~isnan(temp_y)),temp_y2_hat_mean(~isnan(temp_y)));
        
        memoryPrecision1_dynamicProj = temp_y1_hat;
        memoryPrecision2_dynamicProj = temp_y2_hat;
        
        
        
        %% Proj to beta_meta_linearAxis
        temp_x1 = temp_F_dff_baseline_Zdelay;
        temp_x2 = temp_F_dff_lastTDelay;
        temp_y = temp_meta_trialLevel;
        temp_beta = beta_meta_linearAxis;
        temp_beta0 = beta0_meta_linearAxis;
        
        temp_y1_hat = squeeze(pagemtimes(temp_beta',temp_x1) + temp_beta0);
        temp_y1_hat_raw = temp_y1_hat;
        temp_y2_hat = squeeze(pagemtimes(temp_beta',temp_x2) + temp_beta0);
        temp_y2_hat_raw = temp_y2_hat;
        if temp_if_zscore == 1
            temp_y_mean = mean(temp_y,'omitnan');
            temp_y_std = std(temp_y,'omitnan');
            temp_y1_hat = (temp_y1_hat_raw.*temp_y_std)+temp_y_mean;
            temp_y2_hat = (temp_y2_hat_raw.*temp_y_std)+temp_y_mean;
        end
        temp_y1_hat_mean = mean(temp_y1_hat,2);
        temp_y2_hat_mean = mean(temp_y2_hat,2);
        
        temp_r1 = corr(temp_y(~isnan(temp_y)),temp_y1_hat_mean(~isnan(temp_y)));
        temp_r2 = corr(temp_y(~isnan(temp_y)),temp_y2_hat_mean(~isnan(temp_y)));
        
        meta1_dynamicProj = temp_y1_hat;
        meta2_dynamicProj = temp_y2_hat;
        
        %%
        prior_dynamicProj = [prior1_dynamicProj,prior2_dynamicProj];
        memoryPrecision_dynamicProj = [memoryPrecision1_dynamicProj,memoryPrecision2_dynamicProj];
        meta_dynamicProj = [meta1_dynamicProj,meta2_dynamicProj];
        
        %% Plot
        if true
            temp_if_setYBound_all = 0;
            temp_if_setYBound_all_diff = 0;
            
            y_min_all = 0.1;%0
            y_max_all = 0.9;%1
            
            y_min_all_diff = -0.1;%0
            y_max_all_diff = 0.4;%1
            
            temp1A = prior_dynamicProj(choiceMemoryBoolIndex,:);
            temp1B = prior_dynamicProj(choiceOffloadBoolIndex,:);
            temp2A = memoryPrecision_dynamicProj(choiceMemoryBoolIndex,:);
            temp2B = memoryPrecision_dynamicProj(choiceOffloadBoolIndex,:);
            temp3A = meta_dynamicProj(choiceMemoryBoolIndex,:);
            temp3B = meta_dynamicProj(choiceOffloadBoolIndex,:);
            
            temp1A_mean = mean(temp1A,1,'omitnan');
            temp2A_mean = mean(temp2A,1,'omitnan');
            temp3A_mean = mean(temp3A,1,'omitnan');
            temp1B_mean = mean(temp1B,1,'omitnan');
            temp2B_mean = mean(temp2B,1,'omitnan');
            temp3B_mean = mean(temp3B,1,'omitnan');
            
            temp1_diff_mean = temp1A_mean - temp1B_mean;
            temp2_diff_mean = temp2A_mean - temp2B_mean;
            temp3_diff_mean = temp3A_mean - temp3B_mean;
            
            temp_timeNum = size(prior_dynamicProj,2);
            
            %         temp_p1AB = nan(1,temp_timeNum);
            %         temp_p2AB = temp_p1AB;
            %         temp_p3AB = temp_p1AB;
            %         for tempi=1:temp_timeNum
            %             [~,temp_p1AB(tempi)] = ttest2(temp1A(:,tempi),temp1B(:,tempi));
            %             [~,temp_p2AB(tempi)] = ttest2(temp2A(:,tempi),temp2B(:,tempi));
            %             [~,temp_p3AB(tempi)] = ttest2(temp3A(:,tempi),temp3B(:,tempi));
            %         end
            %         temp_hugeDiff1 = temp_p1AB < 0.01;
            %         temp_hugeDiff2 = temp_p2AB < 0.01;
            %         temp_hugeDiff3 = temp_p3AB < 0.01;
            
            %         temp_timeNum1 = size(prior1_dynamicProj,2);
            %         temp_timeNum2 = size(prior2_dynamicProj,2);
            %         temp_p1AB = nan(1,temp_timeNum2);
            %         temp_p2AB = temp_p1AB;
            %         temp_p3AB = temp_p1AB;
            %         for tempi=1:temp_timeNum2
            %             temp_t = tempi + temp_timeNum1;
            %             [~,temp_p1AB(tempi)] = ttest2(temp1_diff_mean(1:temp_timeNum1),temp1_diff_mean(temp_t),'tail','left');
            %             [~,temp_p2AB(tempi)] = ttest2(temp2_diff_mean(1:temp_timeNum1),temp2_diff_mean(temp_t),'tail','left');
            %             [~,temp_p3AB(tempi)] = ttest2(temp3_diff_mean(1:temp_timeNum1),temp3_diff_mean(temp_t),'tail','left');
            %         end
            %         temp_hugeDiff1 = temp_p1AB < 0.01;
            %         temp_hugeDiff2 = temp_p2AB < 0.01;
            %         temp_hugeDiff3 = temp_p3AB < 0.01;
            
            
            if if_plot == 1
                backgounrdColor = [1 1 1];
                
                temp_blankSize = 10;
                
                fig = figure('Name','asd','NumberTitle','off');
                temptemp1 = 360*0.975*0.74*1.37*0.97*1.02*1.015*0.84;
                temptemp2 = 200*0.8*0.8*0.87;
                
                set(gcf,'Position',[400 50+350 temptemp1 temptemp2*2.5]);
                
                t = tiledlayout(3,1,'TileSpacing','loose','Padding','Compact');
                
                
                %% Baseline meta
                nexttile
                
                x = 1:baselinePeriod_interval(end);
                y = temp1A_mean(1:baselinePeriod_interval(end));
                plot(x,y,'color',color_choiceMemory,'linewidth',1.5);
                hold on
                
                x = ((baselinePeriod_interval(end)+1):size(temp1A_mean,2))+temp_blankSize;
                y = temp1A_mean((baselinePeriod_interval(end)+1):size(temp1A_mean,2));
                plot(x,y,'color',color_choiceMemory,'linewidth',1.5);
                hold on
                
                x = 1:baselinePeriod_interval(end);
                y = temp1B_mean(1:baselinePeriod_interval(end));
                plot(x,y,'color',color_choiceOffload,'linewidth',1.5);
                hold on
                
                x = ((baselinePeriod_interval(end)+1):size(temp1A_mean,2))+temp_blankSize;
                y = temp1B_mean((baselinePeriod_interval(end)+1):size(temp1A_mean,2));
                plot(x,y,'color',color_choiceOffload,'linewidth',1.5);
                hold on
                
                
                
                [y_min,y_max] = bounds([temp1A_mean,temp1B_mean]);
                if temp_if_setYBound_all == 1
                    y_min = y_min_all;
                    y_max = y_max_all;
                end
                
                temp_interval_all = [baselinePeriod_interval(1),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)) + temp_blankSize,...
                    (baselinePeriod_interval(end) + temp_blankSize + length1_sample_interval(end) + decisionPeriodA_interval)];
                
                text(mean(temp_interval_all(2:3)),y_min-(y_max-y_min)*0.1,...
                    sprintf('%s','//'),'HorizontalAlignment','center','fontsize',8);
                
                xticks(temp_interval_all([1 2 3 4]));
                %xticklabels({'Fixation','T1','LastT','Delay-on'});
                xticklabels([]);
                
                %yticks([-1 0 0.5]);
                
                xlim([temp_interval_all(1) temp_interval_all(end-1)]);
                
                set(gca,'linewidth',1.5)
                set(gca, 'FontSize', 7.5)%12,8
                set(gca,'color',backgounrdColor);
                set(gca,'box','off');% 取消右、上边框
                
                ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);
                
                xtickangle(0);
                
                ylabel('Baseline meta', 'FontSize', 8, 'FontWeight', 'normal');%
                %title('Meta-memory time course','FontSize',9);
                
                
                %% WM strength
                nexttile
                
                x = 1:baselinePeriod_interval(end);
                y = temp2A_mean(1:baselinePeriod_interval(end));
                plot(x,y,'color',color_choiceMemory,'linewidth',1.5);
                hold on
                
                x = ((baselinePeriod_interval(end)+1):size(temp1A_mean,2))+temp_blankSize;
                y = temp2A_mean((baselinePeriod_interval(end)+1):size(temp1A_mean,2));
                plot(x,y,'color',color_choiceMemory,'linewidth',1.5);
                hold on
                
                x = 1:baselinePeriod_interval(end);
                y = temp2B_mean(1:baselinePeriod_interval(end));
                plot(x,y,'color',color_choiceOffload,'linewidth',1.5);
                hold on
                
                x = ((baselinePeriod_interval(end)+1):size(temp1A_mean,2))+temp_blankSize;
                y = temp2B_mean((baselinePeriod_interval(end)+1):size(temp1A_mean,2));
                plot(x,y,'color',color_choiceOffload,'linewidth',1.5);
                hold on
                
                
                [y_min,y_max] = bounds([temp2A_mean,temp2B_mean]);
                if temp_if_setYBound_all == 1
                    y_min = y_min_all;
                    y_max = y_max_all;
                end
                
                temp_interval_all = [baselinePeriod_interval(1),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)) + temp_blankSize,...
                    (baselinePeriod_interval(end) + temp_blankSize + length1_sample_interval(end) + decisionPeriodA_interval)];
                
                text(mean(temp_interval_all(2:3)),y_min-(y_max-y_min)*0.1,...
                    sprintf('%s','//'),'HorizontalAlignment','center','fontsize',8);
                
                xticks(temp_interval_all([1 2 3 4]));
                %xticklabels({'Fixation','T1','LastT','Delay-on'});
                xticklabels([]);
                
                
                %yticks([-1 0 0.5]);
                
                xlim([temp_interval_all(1) temp_interval_all(end-1)]);
                
                set(gca,'linewidth',1.5)
                set(gca, 'FontSize', 7.5)%12,8
                set(gca,'color',backgounrdColor);
                set(gca,'box','off');% 取消右、上边框
                
                ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);
                
                xtickangle(0);
                
                ylabel('WM strength', 'FontSize', 8, 'FontWeight', 'normal');%
                %title('Meta-memory time course','FontSize',9);
                
                
                %% Meta
                nexttile
                
                x = 1:baselinePeriod_interval(end);
                y = temp3A_mean(1:baselinePeriod_interval(end));
                plot(x,y,'color',color_choiceMemory,'linewidth',1.5);
                hold on
                
                x = ((baselinePeriod_interval(end)+1):size(temp1A_mean,2))+temp_blankSize;
                y = temp3A_mean((baselinePeriod_interval(end)+1):size(temp1A_mean,2));
                plot(x,y,'color',color_choiceMemory,'linewidth',1.5);
                hold on
                
                
                x = 1:baselinePeriod_interval(end);
                y = temp3B_mean(1:baselinePeriod_interval(end));
                plot(x,y,'color',color_choiceOffload,'linewidth',1.5);
                hold on
                
                x = ((baselinePeriod_interval(end)+1):size(temp1A_mean,2))+temp_blankSize;
                y = temp3B_mean((baselinePeriod_interval(end)+1):size(temp1A_mean,2));
                plot(x,y,'color',color_choiceOffload,'linewidth',1.5);
                hold on
                
                
                [y_min,y_max] = bounds([temp3A_mean,temp3B_mean]);
                if temp_if_setYBound_all == 1
                    y_min = y_min_all;
                    y_max = y_max_all;
                end
                
                
                temp_interval_all = [baselinePeriod_interval(1),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)) + temp_blankSize,...
                    (baselinePeriod_interval(end) + temp_blankSize + length1_sample_interval(end) + decisionPeriodA_interval)];
                
                text(mean(temp_interval_all(2:3)),y_min-(y_max-y_min)*0.1,...
                    sprintf('%s','//'),'HorizontalAlignment','center','fontsize',8);
                
                xticks(temp_interval_all([1 2 3 4]));
                xticklabels({'Fixation','T1','LastT','Delay-on'});
                
                %yticks([-1 0 0.5]);
                
                xlim([temp_interval_all(1) temp_interval_all(end-1)]);
                
                set(gca,'linewidth',1.5)
                set(gca, 'FontSize', 7.5)%12,8
                set(gca,'color',backgounrdColor);
                set(gca,'box','off');% 取消右、上边框
                
                ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);
                
                xtickangle(0);
                
                ylabel('Meta', 'FontSize', 8, 'FontWeight', 'normal');%
                %title('Meta-memory time course','FontSize',9);
                
                
                
            end
            
            
            if if_plot == 1
                backgounrdColor = [1 1 1];
                
                temp_blankSize = 10;
                
                fig = figure('Name','asd','NumberTitle','off');
                temptemp1 = 360*0.975*0.74*1.37*0.97*1.02*1.015*0.84;
                temptemp2 = 200*0.8*0.8*0.87;
                
                set(gcf,'Position',[700 50+350 temptemp1 temptemp2*2.5]);
                
                t = tiledlayout(3,1,'TileSpacing','loose','Padding','Compact');
                
                
                %% Baseline meta
                nexttile
                
                x = 1:baselinePeriod_interval(end);
                y = temp1_diff_mean(1:baselinePeriod_interval(end));
                plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1.5);
                hold on
                
                x = ((baselinePeriod_interval(end)+1):size(temp1A_mean,2))+temp_blankSize;
                y = temp1_diff_mean((baselinePeriod_interval(end)+1):size(temp1A_mean,2));
                plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1.5);
                hold on
                
                
                
                [y_min,y_max] = bounds(temp1_diff_mean);
                if temp_if_setYBound_all_diff == 1
                    y_min = y_min_all_diff;
                    y_max = y_max_all_diff;
                end
                
                temp_interval_all = [baselinePeriod_interval(1),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)) + temp_blankSize,...
                    (baselinePeriod_interval(end) + temp_blankSize + length1_sample_interval(end) + decisionPeriodA_interval)];
                
                text(mean(temp_interval_all(2:3)),y_min-(y_max-y_min)*0.1,...
                    sprintf('%s','//'),'HorizontalAlignment','center','fontsize',8);
                
                xticks(temp_interval_all([1 2 3 4]));
                %xticklabels({'Fixation','T1','LastT','Delay-on'});
                xticklabels([]);
                
                %yticks([-1 0 0.5]);
                
                xlim([temp_interval_all(1) temp_interval_all(end-1)]);
                
                set(gca,'linewidth',1.5)
                set(gca, 'FontSize', 7.5)%12,8
                set(gca,'color',backgounrdColor);
                set(gca,'box','off');% 取消右、上边框
                
                ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);
                
                xtickangle(0);
                
                ylabel('Baseline meta', 'FontSize', 8, 'FontWeight', 'normal');%
                %title('Meta-memory time course','FontSize',9);
                
                
                %% WM strength
                nexttile
                
                x = 1:baselinePeriod_interval(end);
                y = temp2_diff_mean(1:baselinePeriod_interval(end));
                plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1.5);
                hold on
                
                x = ((baselinePeriod_interval(end)+1):size(temp1A_mean,2))+temp_blankSize;
                y = temp2_diff_mean((baselinePeriod_interval(end)+1):size(temp1A_mean,2));
                plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1.5);
                hold on
                
                
                
                [y_min,y_max] = bounds(temp2_diff_mean);
                if temp_if_setYBound_all_diff == 1
                    y_min = y_min_all_diff;
                    y_max = y_max_all_diff;
                end
                
                
                temp_interval_all = [baselinePeriod_interval(1),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)) + temp_blankSize,...
                    (baselinePeriod_interval(end) + temp_blankSize + length1_sample_interval(end) + decisionPeriodA_interval)];
                
                text(mean(temp_interval_all(2:3)),y_min-(y_max-y_min)*0.1,...
                    sprintf('%s','//'),'HorizontalAlignment','center','fontsize',8);
                
                xticks(temp_interval_all([1 2 3 4]));
                %xticklabels({'Fixation','T1','LastT','Delay-on'});
                xticklabels([]);
                
                
                %yticks([-1 0 0.5]);
                
                xlim([temp_interval_all(1) temp_interval_all(end-1)]);
                
                set(gca,'linewidth',1.5)
                set(gca, 'FontSize', 7.5)%12,8
                set(gca,'color',backgounrdColor);
                set(gca,'box','off');% 取消右、上边框
                
                ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);
                
                xtickangle(0);
                
                ylabel('WM strength', 'FontSize', 8, 'FontWeight', 'normal');%
                %title('Meta-memory time course','FontSize',9);
                
                
                %% Meta
                nexttile
                
                x = 1:baselinePeriod_interval(end);
                y = temp3_diff_mean(1:baselinePeriod_interval(end));
                plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1.5);
                hold on
                
                x = ((baselinePeriod_interval(end)+1):size(temp1A_mean,2))+temp_blankSize;
                y = temp3_diff_mean((baselinePeriod_interval(end)+1):size(temp1A_mean,2));
                plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1.5);
                hold on
                
                
                [y_min,y_max] = bounds(temp3_diff_mean);
                if temp_if_setYBound_all_diff == 1
                    y_min = y_min_all_diff;
                    y_max = y_max_all_diff;
                end
                
                
                temp_interval_all = [baselinePeriod_interval(1),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)) + temp_blankSize,...
                    (baselinePeriod_interval(end) + temp_blankSize + length1_sample_interval(end) + decisionPeriodA_interval)];
                
                text(mean(temp_interval_all(2:3)),y_min-(y_max-y_min)*0.1,...
                    sprintf('%s','//'),'HorizontalAlignment','center','fontsize',8);
                
                xticks(temp_interval_all([1 2 3 4]));
                xticklabels({'Fixation','T1','LastT','Delay-on'});
                
                %yticks([-1 0 0.5]);
                
                xlim([temp_interval_all(1) temp_interval_all(end-1)]);
                
                set(gca,'linewidth',1.5)
                set(gca, 'FontSize', 7.5)%12,8
                set(gca,'color',backgounrdColor);
                set(gca,'box','off');% 取消右、上边框
                
                ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);
                
                xtickangle(0);
                
                ylabel('Meta', 'FontSize', 8, 'FontWeight', 'normal');%
                %title('Meta-memory time course','FontSize',9);
                
                
                
            end
            
            
            %% AUROC time course
            
            temp_choiceBoolIndex_valid = (~isnan(memoryPrecision_trialLevel)) & choiceBoolIndex';
            %temp_choiceBoolIndex_valid = temp_choiceBoolIndex_valid & (trial_para_isCorrect==1)';
            
            lowThreshold_memoryPrecision;
            memoryPrecision_trialLevel;
            highMemoryPrecisionBoolIndex = memoryPrecision_trialLevel > lowThreshold_memoryPrecision;
            highMetaBoolIndex = meta_trialLevel > lowThreshold_meta;
            
            temp_thresholdRange = 0.01:0.001:0.99;
            trueLabel1 = choiceMemoryBoolIndex(temp_choiceBoolIndex_valid)';
            trueLabel3 = choiceMemoryBoolIndex(temp_choiceBoolIndex_valid)';
            %trueLabel1 = highMetaBoolIndex(temp_choiceBoolIndex_valid)';% huge time diff
            %trueLabel3 = highMetaBoolIndex(temp_choiceBoolIndex_valid)';% huge time diff
            
            %trueLabel2 = choiceMemoryBoolIndex(temp_choiceBoolIndex_valid)';
            trueLabel2 = highMemoryPrecisionBoolIndex(temp_choiceBoolIndex_valid)';% more resonable
            
            AUROC_prior_dynamicProj = nan(1,temp_timeNum);
            AUROC_memoryPrecision_dynamicProj = nan(1,temp_timeNum);
            AUROC_meta_dynamicProj = nan(1,temp_timeNum);
            for tempi=1:temp_timeNum
                temp1_svm_choiceMemory_hit = zeros(1,length(temp_thresholdRange));
                temp1_svm_choiceMemory_falseAlarm = temp1_svm_choiceMemory_hit;
                temp2_svm_choiceMemory_hit = temp1_svm_choiceMemory_hit;
                temp2_svm_choiceMemory_falseAlarm = temp1_svm_choiceMemory_hit;
                temp3_svm_choiceMemory_hit = temp1_svm_choiceMemory_hit;
                temp3_svm_choiceMemory_falseAlarm = temp1_svm_choiceMemory_hit;
                for temptempi=1:length(temp_thresholdRange)
                    temp_threshold = temp_thresholdRange(temptempi);
                    
                    predictLabel1 = prior_dynamicProj(temp_choiceBoolIndex_valid,tempi) > temp_threshold;
                    predictLabel2 = memoryPrecision_dynamicProj(temp_choiceBoolIndex_valid,tempi) > temp_threshold;
                    predictLabel3 = meta_dynamicProj(temp_choiceBoolIndex_valid,tempi) > temp_threshold;
                    
                    temp1_svm_choiceMemory_hit(temptempi) = sum(predictLabel1(trueLabel1==true))/sum(trueLabel1==true);
                    temp1_svm_choiceMemory_falseAlarm(temptempi) = sum(predictLabel1(trueLabel1==false))/sum(trueLabel1==false);
                    temp2_svm_choiceMemory_hit(temptempi) = sum(predictLabel2(trueLabel2==true))/sum(trueLabel2==true);
                    temp2_svm_choiceMemory_falseAlarm(temptempi) = sum(predictLabel2(trueLabel2==false))/sum(trueLabel2==false);
                    temp3_svm_choiceMemory_hit(temptempi) = sum(predictLabel3(trueLabel3==true))/sum(trueLabel3==true);
                    temp3_svm_choiceMemory_falseAlarm(temptempi) = sum(predictLabel3(trueLabel3==false))/sum(trueLabel3==false);
                    
                    
                end
                x = temp1_svm_choiceMemory_falseAlarm(end:-1:1);
                y = temp1_svm_choiceMemory_hit(end:-1:1);
                temp_AUROC1 = trapz(x,y);
                x = temp2_svm_choiceMemory_falseAlarm(end:-1:1);
                y = temp2_svm_choiceMemory_hit(end:-1:1);
                temp_AUROC2 = trapz(x,y);
                x = temp3_svm_choiceMemory_falseAlarm(end:-1:1);
                y = temp3_svm_choiceMemory_hit(end:-1:1);
                temp_AUROC3 = trapz(x,y);
                
                AUROC_prior_dynamicProj(tempi) = temp_AUROC1;
                AUROC_memoryPrecision_dynamicProj(tempi) = temp_AUROC2;
                AUROC_meta_dynamicProj(tempi) = temp_AUROC3;
                
            end
            temp_timeNum1 = size(prior1_dynamicProj,2);
            temp_timeNum2 = size(prior2_dynamicProj,2);
            temp_p1AB = nan(1,temp_timeNum2);
            temp_p2AB = temp_p1AB;
            temp_p3AB = temp_p1AB;
            for tempi=1:temp_timeNum2
                temp_t = tempi + temp_timeNum1;
                [~,temp_p1AB(tempi)] = ttest2(AUROC_prior_dynamicProj(1:temp_timeNum1),AUROC_prior_dynamicProj(temp_t),'tail','left');
                [~,temp_p2AB(tempi)] = ttest2(AUROC_memoryPrecision_dynamicProj(1:temp_timeNum1),AUROC_memoryPrecision_dynamicProj(temp_t),'tail','left');
                [~,temp_p3AB(tempi)] = ttest2(AUROC_meta_dynamicProj(1:temp_timeNum1),AUROC_meta_dynamicProj(temp_t),'tail','left');
            end
            temp_hugeDiff1 = temp_p1AB < 0.01;
            temp_hugeDiff2 = temp_p2AB < 0.01;
            temp_hugeDiff3 = temp_p3AB < 0.01;
            
            
            %% Plot
            temp1_AUROC = AUROC_prior_dynamicProj;
            temp2_AUROC = AUROC_memoryPrecision_dynamicProj;
            temp3_AUROC = AUROC_meta_dynamicProj;
            
            temp_if_setYBound_all_AUROC = 1;
            y_min_all_AUROC = 0.45;%0.47
            y_max_all_AUROC = 0.85;%0.78
            
            if if_plot == 1
                backgounrdColor = [1 1 1];
                
                temp_blankSize = 10;
                
                fig = figure('Name','asd','NumberTitle','off');
                temptemp1 = 360*0.975*0.74*1.37*0.97*1.02*1.015*0.84;
                temptemp2 = 200*0.8*0.8*0.87;
                
                set(gcf,'Position',[1000 50+350 temptemp1 temptemp2*2.5*0.9]);
                
                t = tiledlayout(3,1,'TileSpacing','loose','Padding','Compact');
                
                
                %% Baseline meta
                nexttile
                
                x = 1:baselinePeriod_interval(end);
                y = temp1_AUROC(1:baselinePeriod_interval(end));
                plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1.5);
                hold on
                
                x = ((baselinePeriod_interval(end)+1):size(temp1A_mean,2))+temp_blankSize;
                y = temp1_AUROC((baselinePeriod_interval(end)+1):size(temp1A_mean,2));
                plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1.5);
                hold on
                
                
                
                [y_min,y_max] = bounds(temp1_AUROC);
                if temp_if_setYBound_all_AUROC == 1
                    y_min = y_min_all_AUROC;
                    y_max = y_max_all_AUROC;
                end
                
                temp_interval_all = [baselinePeriod_interval(1),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)) + temp_blankSize,...
                    (baselinePeriod_interval(end) + temp_blankSize + length1_sample_interval(end) + decisionPeriodA_interval)];
                
                text(mean(temp_interval_all(2:3)),y_min-(y_max-y_min)*0.1,...
                    sprintf('%s','//'),'HorizontalAlignment','center','fontsize',8);
                
                xticks(temp_interval_all([1 2 3 4]));
                %xticklabels({'Fixation','T1','LastT','Delay-on'});
                xticklabels([]);
                
                %yticks([-1 0 0.5]);
                
                xlim([temp_interval_all(1) temp_interval_all(end-1)]);
                
                set(gca,'linewidth',1.5)
                set(gca, 'FontSize', 7.5)%12,8
                set(gca,'color',backgounrdColor);
                set(gca,'box','off');% 取消右、上边框
                
                ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);
                
                xtickangle(0);
                
                ylabel('AUROC', 'FontSize', 8, 'FontWeight', 'normal');%
                title('Baseline meta subspace','FontSize',9);
                
                
                %% WM strength
                nexttile
                
                x = 1:baselinePeriod_interval(end);
                y = temp2_AUROC(1:baselinePeriod_interval(end));
                plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1.5);
                hold on
                
                x = ((baselinePeriod_interval(end)+1):size(temp1A_mean,2))+temp_blankSize;
                y = temp2_AUROC((baselinePeriod_interval(end)+1):size(temp1A_mean,2));
                plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1.5);
                hold on
                
                
                
                [y_min,y_max] = bounds(temp2_AUROC);
                if temp_if_setYBound_all_AUROC == 1
                    y_min = y_min_all_AUROC;
                    y_max = y_max_all_AUROC;
                end
                
                
                temp_interval_all = [baselinePeriod_interval(1),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)) + temp_blankSize,...
                    (baselinePeriod_interval(end) + temp_blankSize + length1_sample_interval(end) + decisionPeriodA_interval)];
                
                text(mean(temp_interval_all(2:3)),y_min-(y_max-y_min)*0.1,...
                    sprintf('%s','//'),'HorizontalAlignment','center','fontsize',8);
                
                xticks(temp_interval_all([1 2 3 4]));
                %xticklabels({'Fixation','T1','LastT','Delay-on'});
                xticklabels([]);
                
                
                %yticks([-1 0 0.5]);
                
                xlim([temp_interval_all(1) temp_interval_all(end-1)]);
                
                set(gca,'linewidth',1.5)
                set(gca, 'FontSize', 7.5)%12,8
                set(gca,'color',backgounrdColor);
                set(gca,'box','off');% 取消右、上边框
                
                ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);
                
                xtickangle(0);
                
                ylabel('AUROC', 'FontSize', 8, 'FontWeight', 'normal');%
                title('WM strength subspace','FontSize',9);
                
                
                %% Meta
                nexttile
                
                x = 1:baselinePeriod_interval(end);
                y = temp3_AUROC(1:baselinePeriod_interval(end));
                plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1.5);
                hold on
                
                x = ((baselinePeriod_interval(end)+1):size(temp1A_mean,2))+temp_blankSize;
                y = temp3_AUROC((baselinePeriod_interval(end)+1):size(temp1A_mean,2));
                plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1.5);
                hold on
                
                
                [y_min,y_max] = bounds(temp3_AUROC);
                if temp_if_setYBound_all_AUROC == 1
                    y_min = y_min_all_AUROC;
                    y_max = y_max_all_AUROC;
                end
                
                
                temp_interval_all = [baselinePeriod_interval(1),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)),...
                    (baselinePeriod_interval(end) + length1_sample_interval(1:end-1)) + temp_blankSize,...
                    (baselinePeriod_interval(end) + temp_blankSize + length1_sample_interval(end) + decisionPeriodA_interval)];
                
                text(mean(temp_interval_all(2:3)),y_min-(y_max-y_min)*0.1,...
                    sprintf('%s','//'),'HorizontalAlignment','center','fontsize',8);
                
                xticks(temp_interval_all([1 2 3 4]));
                xticklabels({'Fixation','T1','LastT','Delay-on'});
                
                %yticks([-1 0 0.5]);
                
                xlim([temp_interval_all(1) temp_interval_all(end-1)]);
                
                set(gca,'linewidth',1.5)
                set(gca, 'FontSize', 7.5)%12,8
                set(gca,'color',backgounrdColor);
                set(gca,'box','off');% 取消右、上边框
                
                ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);
                
                xtickangle(0);
                
                ylabel('AUROC', 'FontSize', 8, 'FontWeight', 'normal');%
                title('Meta-WM subspace','FontSize',9);
                
                
                
            end
            
            
        end
        
    end
    
    
    %% Rotation dynamics of choice axis
    if true
        %t_ratationDynamics = tic;
        
        if_meta_decoder0_behavior1 = 1;%0,1
        
        temp_meta_trialLevel_timebin = double(choiceMemoryBoolIndex');
        
        
        if_fitglm0_lasso1 = 1;%0
        
        if_resample0_twoSession1 = 0;%1,0
        
        temp_resampleNum = 16*7;%16*7(63 secs in lasso),16*2(20 secs in lasso)
        if_resample_replace = 0;%0
        
        temp_if_allTrials0_partTrials1 = 1;%1
        
        if if_meta_decoder0_behavior1 == 1
            temp_if_allTrials0_partTrials1 = 1;
        end
        
        temp_if_zscore = 1;%1
        
        temp_if_resample_balance = 0;%0
        temp_ifBalance_lowHigh0_seq1 = 1;
        
        if if_resample0_twoSession1 == 1
            temp_resampleNum = 1;
            if_resample_replace = 0;
        end
        
        fprintf('Rotation dynamics of choice axis: if_meta_decoder0_behavior1=%d,if_fitglm0_lasso1=%d,if_resample0_twoSession1=%d,temp_if_allTrials0_partTrials1=%d,temp_if_zscore=%d,',...
            if_meta_decoder0_behavior1,if_fitglm0_lasso1,if_resample0_twoSession1,temp_if_allTrials0_partTrials1,temp_if_zscore);
        
        if if_resample0_twoSession1 == 0
            fprintf('temp_resampleNum=%d,if_resample_replace=%d,temp_if_resample_balance=%d,temp_ifBalance_lowHigh0_seq1=%d.\n',...
                temp_resampleNum,if_resample_replace,temp_if_resample_balance,temp_ifBalance_lowHigh0_seq1);
        end
        
        
        temp_roiNum = size(F_dff_decisionPeriodA,1);
        
        temp_F_dff_baseline;
        temp_F_dff_lastTDelay_Zbaseline;
        temp_F_dff_baselineLastTDelay_Zbaseline = cat(3,temp_F_dff_baseline,temp_F_dff_lastTDelay_Zbaseline);
        temp_F_dff_baselineLastTDelay_Zbaseline;
        
        temp_timeBinNum = size(temp_F_dff_baselineLastTDelay_Zbaseline,3);
        %temp_timeBinNum = 16;
        
        
        %% Rotation dynamics (WM strength): all choice trials
        if true
            t_ratationDynamics = tic;
            
            r2_memoryDynamics_linearAxis = nan(1,temp_timeBinNum);
            beta_memoryDynamics_linearAxis = nan(temp_roiNum,temp_timeBinNum);
            beta0_memoryDynamicsr_linearAxis = nan(1,temp_timeBinNum);
            
            %for tempi=1:temp_timeBinNum
            parfor tempi=1:temp_timeBinNum
                
                x = temp_F_dff_baselineLastTDelay_Zbaseline(:,choiceBoolIndex,tempi)';
                y = memoryPrecision_trialLevel(choiceBoolIndex);
                if temp_if_zscore == 1
                    y_zscore = (y-mean(y,'omitnan'))./std(y,'omitnan');%zscore
                    if temp_if_resample_balance==1 && temp_ifBalance_lowHigh0_seq1==0
                        a = 1;
                    else
                        y = y_zscore;
                    end
                end
                
                if if_fitglm0_lasso1 == 0
                    temp_mdl = fitglm(x,y);
                    r2_memoryDynamics_linearAxis(tempi) = temp_mdl.Rsquared.Adjusted;
                    beta_memoryDynamics_linearAxis(:,tempi) = temp_mdl.Coefficients.Estimate(2:end);
                    beta0_memoryDynamicsr_linearAxis(tempi) = temp_mdl.Coefficients.Estimate(1);
                elseif if_fitglm0_lasso1 == 1
                    [beta_memoryDynamics_linearAxis(:,tempi),r2_memoryDynamics_linearAxis(tempi),beta0_memoryDynamicsr_linearAxis(tempi)] = fun_lasso_linearAxis_jjb(x,y);
                end
                
            end
            
            r2_memoryDynamics_linearAxis;
            beta_memoryDynamics_linearAxis;
                        
            temp_degree_rotationDynamics = nan(1,temp_timeBinNum);
            for tempi=1:temp_timeBinNum
                temptemp_vectors = [mean(beta_memoryDynamics_linearAxis(:,53:end),2),beta_memoryDynamics_linearAxis(:,tempi)];
                
                temp1 = temptemp_vectors(:,1);
                temp2 = temptemp_vectors(:,2);
                
                temp_degree_rotationDynamics(tempi) = subspace(temp1,temp2)*180/pi;
            end
            
            temp1 = median(temp_degree_rotationDynamics(1:34));
            temp2 = median(temp_degree_rotationDynamics(35:52));
            temp3 = median(temp_degree_rotationDynamics(53:end));
            
            temp123_min = min([temp1,temp2,temp3]);
            temp123_max = max([temp1,temp2,temp3]);
            
            temp1_n11n = (temp1-temp123_min) ./ (temp123_max-temp123_min);
            temp2_n11n = (temp2-temp123_min) ./ (temp123_max-temp123_min);
            temp3_n11n = (temp3-temp123_min) ./ (temp123_max-temp123_min);
            
            fprintf('Rotation dynamics (WM strength): [baseline, lastT, delay] = [%.3f, %.3f, %.3f] (%.3f, %.3f, %.3f).\n',temp1,temp2,temp3,temp1_n11n,temp2_n11n,temp3_n11n);
            
            
            fprintf('t_ratationDynamics = %.1f secs.\n\n',toc(t_ratationDynamics));
            
        end        
        
        %% Rotation dynamics: all choice trials
        if true
            t_ratationDynamics = tic;
            
            r2_choiceDynamics_linearAxis = nan(1,temp_timeBinNum);
            beta_choiceDynamics_linearAxis = nan(temp_roiNum,temp_timeBinNum);
            beta0_choiceDynamicsr_linearAxis = nan(1,temp_timeBinNum);
            
            %for tempi=1:temp_timeBinNum
            parfor tempi=1:temp_timeBinNum
                
                x = temp_F_dff_baselineLastTDelay_Zbaseline(:,choiceBoolIndex,tempi)';
                y = temp_meta_trialLevel_timebin(choiceBoolIndex);
                if temp_if_zscore == 1
                    y_zscore = (y-mean(y,'omitnan'))./std(y,'omitnan');%zscore
                    if temp_if_resample_balance==1 && temp_ifBalance_lowHigh0_seq1==0
                        a = 1;
                    else
                        y = y_zscore;
                    end
                end
                
                if if_fitglm0_lasso1 == 0
                    temp_mdl = fitglm(x,y);
                    r2_choiceDynamics_linearAxis(tempi) = temp_mdl.Rsquared.Adjusted;
                    beta_choiceDynamics_linearAxis(:,tempi) = temp_mdl.Coefficients.Estimate(2:end);
                    beta0_choiceDynamicsr_linearAxis(tempi) = temp_mdl.Coefficients.Estimate(1);
                elseif if_fitglm0_lasso1 == 1
                    [beta_choiceDynamics_linearAxis(:,tempi),r2_choiceDynamics_linearAxis(tempi),beta0_choiceDynamicsr_linearAxis(tempi)] = fun_lasso_linearAxis_jjb(x,y);
                    %[a,b,c] = fun_lasso_linearAxis_jjb(x,y);
                end
                
            end
            
            r2_choiceDynamics_linearAxis;
            beta_choiceDynamics_linearAxis;
            
            beta_prior_linearAxis;
            
            temp_degree = nan(temp_timeBinNum,temp_timeBinNum);
            
            for tempi=1:temp_timeBinNum
                for tempj=1:temp_timeBinNum
                    temptemp_vectors = [beta_choiceDynamics_linearAxis(:,tempj),beta_choiceDynamics_linearAxis(:,tempi)];
                    
                    temp1 = temptemp_vectors(:,1);
                    temp2 = temptemp_vectors(:,2);
                    temp_degree(tempi,tempj) = subspace(temp1,temp2)*180/pi;
                end
            end
            %temp_degree(35:end,:) = 0;
            %temp_degree(:,35:end) = 0;
            temp_degree(1:34,:) = 0;
            temp_degree(:,1:34) = 0;
            
            temp_degree_sum = sum(temp_degree);
            %[M,I] = min(temp_degree_sum(1:34));
            temp_degree_sum(1:34) = nan;
            [M,I] = min(temp_degree_sum);
            
            temp_degree_rotationDynamics = nan(1,temp_timeBinNum);
            for tempi=1:temp_timeBinNum
                %temptemp_vectors = [beta_prior_linearAxis,beta_choiceDynamics_linearAxis(:,tempi)];
                %temptemp_vectors = [beta_choiceDynamics_linearAxis(:,28),beta_choiceDynamics_linearAxis(:,tempi)];
                %temptemp_vectors = [mean(beta_choiceDynamics_linearAxis(:,1:34),2),beta_choiceDynamics_linearAxis(:,tempi)];
                %temptemp_vectors = [beta_choiceDynamics_linearAxis(:,I),beta_choiceDynamics_linearAxis(:,tempi)];
                %temptemp_vectors = [beta_meta_linearAxis,beta_choiceDynamics_linearAxis(:,tempi)];
                temptemp_vectors = [mean(beta_choiceDynamics_linearAxis(:,53:end),2),beta_choiceDynamics_linearAxis(:,tempi)];
                
                temp1 = temptemp_vectors(:,1);
                temp2 = temptemp_vectors(:,2);
                
                temp_degree_rotationDynamics(tempi) = subspace(temp1,temp2)*180/pi;
            end
            
            temp1 = median(temp_degree_rotationDynamics(1:34));
            temp2 = median(temp_degree_rotationDynamics(35:52));
            temp3 = median(temp_degree_rotationDynamics(53:end));
            
            temp123_min = min([temp1,temp2,temp3]);
            temp123_max = max([temp1,temp2,temp3]);
            
            temp1_n11n = (temp1-temp123_min) ./ (temp123_max-temp123_min);
            temp2_n11n = (temp2-temp123_min) ./ (temp123_max-temp123_min);
            temp3_n11n = (temp3-temp123_min) ./ (temp123_max-temp123_min);
            
            %fprintf('Rotation dynamics: [baseline, lastT, delay] = [%.3f, %.3f, %.3f].\n',temp1,temp2,temp3);
            fprintf('Rotation dynamics: [baseline, lastT, delay] = [%.3f, %.3f, %.3f] (%.3f, %.3f, %.3f).\n',temp1,temp2,temp3,temp1_n11n,temp2_n11n,temp3_n11n);
            
            
            fprintf('t_ratationDynamics = %.1f secs.\n\n',toc(t_ratationDynamics));
            
        end
        
        
        
        %% Rotation dynamics: choice & high WM strength trials
        if true
            t_ratationDynamics_2 = tic;
            
            r2_choiceDynamics_highStrength_linearAxis = nan(1,temp_timeBinNum);
            beta_choiceDynamics_highStrength_linearAxis = nan(temp_roiNum,temp_timeBinNum);
            beta0_choiceDynamicsr_highStrength_linearAxis = nan(1,temp_timeBinNum);
            
            choiceBoolIndex = choiceBoolIndex;
            
            meta_trialLevel_baseline = meta_trialLevel_baseline;
            
            %temp_highStrengthBoolIndex = memoryPrecision_trialLevel > -1;            
            %temp_highStrengthBoolIndex = memoryPrecision_trialLevel > lowThreshold_memoryPrecision;
            %temp_highStrengthBoolIndex = memoryPrecision_trialLevel <= lowThreshold_memoryPrecision;
            %temp_highStrengthBoolIndex = meta_trialLevel_baseline > lowThreshold_meta;
            temp_highStrengthBoolIndex = meta_trialLevel_baseline <= lowThreshold_meta;
            
            temp1 = prctile(meta_trialLevel_baseline,25);%25,15
            temp2 = prctile(meta_trialLevel_baseline,75);%75,85
            %temp_highStrengthBoolIndex = (meta_trialLevel_baseline < temp1) | (meta_trialLevel_baseline > temp2);%Extreme baseline
            %temp_highStrengthBoolIndex = (meta_trialLevel_baseline > temp1) & (meta_trialLevel_baseline < temp2);%Medium baseline
            
            temp1 = prctile(memoryPrecision_trialLevel,25);%25,15
            temp2 = prctile(memoryPrecision_trialLevel,75);%75,85
            %temp_highStrengthBoolIndex = (memoryPrecision_trialLevel < temp1) | (memoryPrecision_trialLevel > temp2);%Extreme strength
            %temp_highStrengthBoolIndex = (memoryPrecision_trialLevel > temp1) & (memoryPrecision_trialLevel < temp2);%Medium strength
            
            temp1 = prctile(memoryPrecision_trialLevel,16.5);%25,15,33
            temp2 = prctile(memoryPrecision_trialLevel,83.5);%75,85,66
            %temp_highStrengthBoolIndex = memoryPrecision_trialLevel < temp1;%Low 1/3 strength
            %temp_highStrengthBoolIndex = (memoryPrecision_trialLevel > temp1) & (memoryPrecision_trialLevel < temp2);%Middle 1/3 strength
            %temp_highStrengthBoolIndex = memoryPrecision_trialLevel > temp2;%High 1/3 strength
           
            %for tempi=1:temp_timeBinNum
            parfor tempi=1:temp_timeBinNum
                
                x = temp_F_dff_baselineLastTDelay_Zbaseline(:,choiceBoolIndex&temp_highStrengthBoolIndex',tempi)';
                y = temp_meta_trialLevel_timebin(choiceBoolIndex&temp_highStrengthBoolIndex');
                if temp_if_zscore == 1
                    y_zscore = (y-mean(y,'omitnan'))./std(y,'omitnan');%zscore
                    if temp_if_resample_balance==1 && temp_ifBalance_lowHigh0_seq1==0
                        a = 1;
                    else
                        y = y_zscore;
                    end
                end
                
                if if_fitglm0_lasso1 == 0
                    temp_mdl = fitglm(x,y);
                    r2_choiceDynamics_highStrength_linearAxis(tempi) = temp_mdl.Rsquared.Adjusted;
                    beta_choiceDynamics_highStrength_linearAxis(:,tempi) = temp_mdl.Coefficients.Estimate(2:end);
                    beta0_choiceDynamicsr_highStrength_linearAxis(tempi) = temp_mdl.Coefficients.Estimate(1);
                elseif if_fitglm0_lasso1 == 1
                    [beta_choiceDynamics_highStrength_linearAxis(:,tempi),r2_choiceDynamics_highStrength_linearAxis(tempi),beta0_choiceDynamicsr_highStrength_linearAxis(tempi)] = fun_lasso_linearAxis_jjb(x,y);
                end
                
            end
            
            r2_choiceDynamics_highStrength_linearAxis;
            beta_choiceDynamics_highStrength_linearAxis;
                        
            temp_degree_rotationDynamics_highStrength = nan(1,temp_timeBinNum);
            for tempi=1:temp_timeBinNum
                % align to self
                %temptemp_vectors = [mean(beta_choiceDynamics_highStrength_linearAxis(:,53:end),2),beta_choiceDynamics_highStrength_linearAxis(:,tempi)];
                
                % align to all
                %temptemp_vectors = [mean(beta_choiceDynamics_linearAxis(:,53:end),2),beta_choiceDynamics_highStrength_linearAxis(:,tempi)];
                
                % align to WM strength axis (delay)
                %temptemp_vectors = [beta_memory_linearAxis,beta_choiceDynamics_highStrength_linearAxis(:,tempi)];
                %temptemp_vectors = [mean(beta_memoryDynamics_linearAxis(:,53:end),2),beta_choiceDynamics_highStrength_linearAxis(:,tempi)];
                
                % align to WM strength axis (timebin)
                %temptemp_vectors = [beta_memoryDynamics_linearAxis(:,tempi),beta_choiceDynamics_highStrength_linearAxis(:,tempi)];
                
                % align to WM strength axis (sample)
                temptemp_vectors = [mean(beta_memoryDynamics_linearAxis(:,35:52),2),beta_choiceDynamics_highStrength_linearAxis(:,tempi)];
                
                temp1 = temptemp_vectors(:,1);
                temp2 = temptemp_vectors(:,2);
                
                temp_degree_rotationDynamics_highStrength(tempi) = subspace(temp1,temp2)*180/pi;
            end
            
            temp1 = median(temp_degree_rotationDynamics_highStrength(1:34));
            temp2 = median(temp_degree_rotationDynamics_highStrength(35:52));
            temp3 = median(temp_degree_rotationDynamics_highStrength(53:end));
            
            temp123_min = min([temp1,temp2,temp3]);
            temp123_max = max([temp1,temp2,temp3]);
            
            temp1_n11n = (temp1-temp123_min) ./ (temp123_max-temp123_min);
            temp2_n11n = (temp2-temp123_min) ./ (temp123_max-temp123_min);
            temp3_n11n = (temp3-temp123_min) ./ (temp123_max-temp123_min);
           
            %fprintf('Rotation dynamics (highStrength): [baseline, lastT, delay] = [%.3f, %.3f, %.3f].\n',temp1,temp2,temp3);
            %fprintf('Rotation dynamics (lowStrength): [baseline, lastT, delay] = [%.3f, %.3f, %.3f].\n',temp1,temp2,temp3);
            fprintf('Rotation dynamics (specific condition): [baseline, lastT, delay] = [%.3f, %.3f, %.3f] (%.3f, %.3f, %.3f).\n',temp1,temp2,temp3,temp1_n11n,temp2_n11n,temp3_n11n);
            
            
            fprintf('t_ratationDynamics_2 = %.1f secs.\n\n',toc(t_ratationDynamics_2));
            
        end
        
        
    end
    
    
end



% %% Behavioral correlation between offloading rate & recall performance in high & low baseline meta
% meta_trialLevel_baseline;
% lowThreshold_meta;
% 
% 


%% End