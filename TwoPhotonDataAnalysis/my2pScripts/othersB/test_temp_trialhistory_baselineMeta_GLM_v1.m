r2_rangeHistory = zeros(1,length(range_trialHistory));

for tempi=1:size(cumulativeScore_trialHistory,2)
    tempNANBoolIndex1 = isnan(cumulativeScore_trialHistory(:,tempi));
    tempNANBoolIndex3 = isnan(meta_trialLevel_baseline);
    
    tempNANBoolIndex12 = tempNANBoolIndex1 | tempNANBoolIndex2;
    
    tempNONNANBoolIndex12_raw = ~tempNANBoolIndex12;
    
    tempNONNANBoolIndex12 = tempNONNANBoolIndex12_raw;    
%     tempNONNANBoolIndex12 = tempNONNANBoolIndex12_raw & choiceBoolIndex';
%     tempNONNANBoolIndex12 = tempNONNANBoolIndex12_raw & ~choiceBoolIndex';    
    
    
    
    tempNONNAN_cumulativeScore_trialHistory = cumulativeScore_trialHistory(tempNONNANBoolIndex12,tempi);
    tempNONNAN_meta_trialLevel_baseline = meta_trialLevel_baseline(tempNONNANBoolIndex12);
    
    x = tempNONNAN_cumulativeScore_trialHistory;
    y = tempNONNAN_meta_trialLevel_baseline;
    temp_mdl_trialHistory = fitglm(x,y,'linear');
    r2_trialHistory = temp_mdl_trialHistory.Rsquared.Adjusted;
    beta0_trialHistory = temp_mdl_trialHistory.Coefficients.Estimate(1);
    beta1_trialHistory = temp_mdl_trialHistory.Coefficients.Estimate(2);
    
    r2_rangeHistory(tempi) = r2_trialHistory;
end
max(r2_rangeHistory)

