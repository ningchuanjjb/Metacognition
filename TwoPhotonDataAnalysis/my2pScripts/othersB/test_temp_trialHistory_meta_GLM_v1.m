cumulativeScore_trialHistory;

a1 = cumulativeScore_trialHistory(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh,:);
a2 = cumulativeScore_trialHistory(trialBoolIndex_memoryPrecisionLow_choiceOffloadLow,:);

a3 = cumulativeScore_trialHistory(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh,:);
a4 = cumulativeScore_trialHistory(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow,:);

a5 = cumulativeScore_trialHistory(trialBoolIndex_metaHigh_choice_baseline,:);
a6 = cumulativeScore_trialHistory(trialBoolIndex_metaLow_choice_baseline,:);

r2_caseB_rangeHistory = zeros(1,length(range_trialHistory));

for tempi=1:size(cumulativeScore_trialHistory,2)
    tempNANBoolIndex1 = isnan(cumulativeScore_trialHistory(:,tempi));
    tempNANBoolIndex2 = isnan(memoryPrecision_trialLevel);
    tempNANBoolIndex3 = isnan(meta_trialLevel);
    
    tempNANBoolIndex123 = tempNANBoolIndex1 | tempNANBoolIndex2 | tempNANBoolIndex3;
    
    tempNONNANBoolIndex123 = ~tempNANBoolIndex123;
    
    tempNONNAN_cumulativeScore_trialHistory = cumulativeScore_trialHistory(tempNONNANBoolIndex123,tempi);
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
    x = tempNONNAN_cumulativeScore_trialHistory;
    y = tempNONNAN_meta_trialLevel;
    temp_mdl_caseB = fitglm(x,y,'linear');
    r2_caseB = temp_mdl_caseB.Rsquared.Adjusted;
    beta0_caseB = temp_mdl_caseB.Coefficients.Estimate(1);
    beta1_caseB = temp_mdl_caseB.Coefficients.Estimate(2);
    
    r2_caseB_rangeHistory(tempi) = r2_caseB;
    
    % Case C
    x = [tempNONNAN_memoryPrecision_trialLevel,tempNONNAN_cumulativeScore_trialHistory];
    y = tempNONNAN_meta_trialLevel;
    temp_mdl_caseC = fitglm(x,y,'linear');
    r2_caseC = temp_mdl_caseC.Rsquared.Adjusted;
    beta0_caseC = temp_mdl_caseC.Coefficients.Estimate(1);
    beta1_caseC= temp_mdl_caseC.Coefficients.Estimate(2);
    beta2_caseC= temp_mdl_caseC.Coefficients.Estimate(3);    
end
max(r2_caseB_rangeHistory)


