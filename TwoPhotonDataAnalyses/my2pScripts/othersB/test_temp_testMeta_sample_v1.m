F_dff_length1_sample;
length1_sample_interval;
F_dff_length2_sample;
length2_sample_interval;
F_dff_length3_sample;
length3_sample_interval;





    temptemp_range = (baselinePeriod_interval(2)-8):baselinePeriod_interval(3);
    F_dff_baselineBin = mean(F_dff_baselinePeriod(:,:,temptemp_range),3);
    F_dff_baselineBin = double(F_dff_baselineBin);
    F_dff_baselineBin = F_dff_baselineBin + eps;
    
    F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
    F_dff_decisionBin1 = double(F_dff_decisionBin1);
    F_dff_decisionBin1 = F_dff_decisionBin1 + eps;
    
    fun_testMeta_currentTime_Name_v = autoGetFunName_myScripts('fun_testMeta_currentTime', [targetPATH '\functions']);
    fun_testMeta_currentTime = str2func(fun_testMeta_currentTime_Name_v);
    
    options_testMeta = struct;
    options_testMeta.KFold_num = KFold_num;
    options_testMeta.svm_Meta = svm_Meta;
    options_testMeta.choiceBoolIndex = choiceBoolIndex;
    options_testMeta.seqIndex = seqIndex;
    
    
    F_dff_baselinePeriod;
    F_dff_length1_sample;
    F_dff_length2_sample;
    F_dff_length3_sample;
    F_dff_decisionPeriodA;
        
    %F_dff_current = F_dff_baselinePeriod;
    F_dff_current = F_dff_length1_sample;    
    %F_dff_current = F_dff_decisionPeriodA;
    
    
    temptemp_range = 1:size(F_dff_current,3);
    
    meta_trialLevel_crossTime = zeros(length(seqIndex),length(temptemp_range));
    %for tempi=temptemp_range(1):temptemp_range(end)
    parfor tempi=temptemp_range(1):temptemp_range(end)                
        F_dff_currentTime = F_dff_current(:,:,temptemp_range(tempi)); %#ok<*PFBNS>
        testMeta_output = fun_testMeta_currentTime(F_dff_currentTime,options_testMeta);        
        meta_trialLevel_crossTime(:,tempi) = testMeta_output.meta_trialLevel_currentTime;
    end
