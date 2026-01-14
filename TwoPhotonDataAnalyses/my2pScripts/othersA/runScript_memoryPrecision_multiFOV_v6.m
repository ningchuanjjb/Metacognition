%% Initialization
% clear
close all

home;

if_runScript_trainDecoder = 1;
if_runScript_testCorrectFew = 1;
if_runScript_testError = 1;
if_runScript_computeLocDistriResample = 1;
if_runScript_save = 0;



currentSession_multi = string;


% currentSession_multi = [currentSession_multi; '113Recording_20230426A_Ding_Site16'];
% currentSession_multi = [currentSession_multi; '113Recording_20230427A_Ding_Site16_sameFOV0426'];
% currentSession_multi = [currentSession_multi; '113Recording_20230502A_Ding_Site13'];
% currentSession_multi = [currentSession_multi; '113Recording_20230503A_Ding_Site13_sameFOV0502'];
currentSession_multi = [currentSession_multi; '113Recording_20230504A_Ding_Site02'];
% currentSession_multi = [currentSession_multi; '113Recording_20230508A_Ding_Site02_sameFOV0509'];
% currentSession_multi = [currentSession_multi; '113Recording_20230509A_Ding_Site02'];% 660000 frames, easy to crash
% 
% currentSession_multi = [currentSession_multi; '113Recording_20230510A_Ding_Site05_sameFOV0511'];
% currentSession_multi = [currentSession_multi; '113Recording_20230510B_Ding_Site05_sameFOV0511'];
% currentSession_multi = [currentSession_multi; '113Recording_20230511A_Ding_Site05'];
% currentSession_multi = [currentSession_multi; '113Recording_20230512A_Ding_Site09'];
% currentSession_multi = [currentSession_multi; '113Recording_20230513A_Ding_Site09_sameFOV0512'];
% 
% currentSession_multi = [currentSession_multi; '113Recording_20230515A_Ding_Site24_sameFOV0516'];
% currentSession_multi = [currentSession_multi; '113Recording_20230516A_Ding_Site24'];
% currentSession_multi = [currentSession_multi; '113Recording_20230517A_Ding_Site16B'];
% currentSession_multi = [currentSession_multi; '113Recording_20230522A_Ding_Site05B'];
% currentSession_multi = [currentSession_multi; '113Recording_20230525A_Ding_Site09B'];
% currentSession_multi = [currentSession_multi; '113Recording_20230526A_Ding_Site09B_sameFOV0525'];
% 
% currentSession_multi = [currentSession_multi; '113Recording_20230527A_Ding_Site02B'];
% % currentSession_multi = [currentSession_multi; '113Recording_20230528A_Ding_Site02B_sameFOV0527'];
% currentSession_multi = [currentSession_multi; '113Recording_20230530A_Ding_Site05C'];
% currentSession_multi = [currentSession_multi; '113Recording_20230531A_Ding_Site05C_sameFOV0530'];
% currentSession_multi = [currentSession_multi; '113Recording_20230601A_Ding_Site13B'];
% currentSession_multi = [currentSession_multi; '113Recording_20230602A_Ding_Site13B_sameFOV0601'];
% 
% currentSession_multi = [currentSession_multi; '113Recording_20230604A_Ding_Site07'];
% currentSession_multi = [currentSession_multi; '113Recording_20230605A_Ding_Site07_sameFOV0604'];
% currentSession_multi = [currentSession_multi; '113Recording_20230612A_Ding_Site14'];
% currentSession_multi = [currentSession_multi; '113Recording_20230614A_Ding_Site05D'];
% currentSession_multi = [currentSession_multi; '113Recording_20230615A_Ding_Site05D_sameFOV0614'];
% currentSession_multi = [currentSession_multi; '113Recording_20230619A_Ding_Site02C'];
% currentSession_multi = [currentSession_multi; '113Recording_20230620A_Ding_Site05E'];

currentSession_multi(1) = [];
num_FOV = length(currentSession_multi);



targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)


stepD3_train123_test123_Name_v = autoGetFunName_myScripts('stepD3_train123_test123', targetPATH);
script_stepD3_train123_test123 = str2func(stepD3_train123_test123_Name_v);
test_testingSet_correct_fewTrialCount_Name_v = autoGetFunName_myScripts('test_testingSet_correct_fewTrialCount', targetPATH);
script_test_testingSet_correct_fewTrialCount = str2func(test_testingSet_correct_fewTrialCount_Name_v);
test_testingSet_error_Name_v = autoGetFunName_myScripts('test_testingSet_error', targetPATH);
script_test_testingSet_error = str2func(test_testingSet_error_Name_v);
test_testingSet_error_eachResampleIter_Name_v = autoGetFunName_myScripts('test_testingSet_error_eachResampleIter', targetPATH);
script_test_testingSet_error_eachResampleIter = str2func(test_testingSet_error_eachResampleIter_Name_v);

a = 1; %#ok<*NASGU>



time0 = tic;
for tempSessionIndex=1:num_FOV
    currentSession = currentSession_multi{tempSessionIndex};
    fprintf('currentSession = %s.\n',currentSession);
    
    
    %% stepD3_train123_test123
    if if_runScript_trainDecoder == 1
        fprintf('\nNow runing %s.  ------> \n', stepD3_train123_test123_Name_v);
        
        
        %clear decodingData
        
        script_stepD3_train123_test123();
        
        behaviorAcc = [gAcc_length1';gAcc_length2';gAcc_length3'];
        decodingAcc_correct = [svm_train_length1_outputs.svm_posterior_lengthx';svm_train_length2_outputs.svm_posterior_lengthx';svm_train_length3_outputs.svm_posterior_lengthx'];
        pProd_correct = [svm_train_length1_outputs.p_seq_prod';svm_train_length2_outputs.p_seq_prod';svm_train_length3_outputs.p_seq_prod'];
        
        tempBoolIndex_nonNAN = ~isnan(decodingAcc_correct);
        [r_decodingAcc_correct,p_decodingAcc_correct] = corr(decodingAcc_correct(tempBoolIndex_nonNAN),behaviorAcc(tempBoolIndex_nonNAN));
        [r_pProd_correct,p_pProd_correct] = corr(pProd_correct(tempBoolIndex_nonNAN),behaviorAcc(tempBoolIndex_nonNAN));
        
        %home;
        fprintf('r_decodingAcc_correct=%.4f, p_decodingAcc_correct=%.4f.\n',r_decodingAcc_correct,p_decodingAcc_correct);
        fprintf('r_pProd_correct=%.4f, p_pProd_correct=%.4f.\n',r_pProd_correct,p_pProd_correct);
    end
    
    a = 1;
    %% test_testingSet_correct_fewTrialCount
    if if_runScript_testCorrectFew == 1
        fprintf('\nNow runing %s.  ------> \n', test_testingSet_correct_fewTrialCount_Name_v);
        script_test_testingSet_correct_fewTrialCount();
        decodingAcc_correctFew = [svm_posterior_length1_decodingAcc';svm_posterior_length2_decodingAcc';svm_posterior_length3_decodingAcc'];
        pProd_correctFew = [svm_posterior_length1_pProd';svm_posterior_length2_pProd';svm_posterior_length3_pProd'];
        
        tempBoolIndex_nonNAN = ~isnan(decodingAcc_correctFew);
        [r_decodingAcc_correctFew,p_decodingAcc_correctFew] = corr(decodingAcc_correctFew(tempBoolIndex_nonNAN),behaviorAcc(tempBoolIndex_nonNAN));
        [r_pProd_correctFew,p_pProd_correctFew] = corr(pProd_correctFew(tempBoolIndex_nonNAN),behaviorAcc(tempBoolIndex_nonNAN));
        
        %home;
        fprintf('r_decodingAcc_correctFew=%.4f, p_decodingAcc_correctFew=%.4f.\n',r_decodingAcc_correctFew,p_decodingAcc_correctFew);
        fprintf('r_pProd_correctFew=%.4f, p_pProd_correctFew=%.4f.\n',r_pProd_correctFew,p_pProd_correctFew);
    end
    
    %% test_testingSet_error
    if if_runScript_testError == 1
        fprintf('\nNow runing %s.  ------> \n', test_testingSet_error_Name_v);
        script_test_testingSet_error();
        
        decodingAcc_stimuliError = svm_posterior_stimuliErrorTrial';
        pProd_stimuliError = p_seq_prod_stimuli';
        
        tempBoolIndex_nonNAN = ~isnan(decodingAcc_stimuliError);
        [r_decodingAcc_stimuliError,p_decodingAcc_stimuliError] = corr(decodingAcc_stimuliError(tempBoolIndex_nonNAN),behaviorAcc(tempBoolIndex_nonNAN));
        [r_pProd_stimuliError,p_pProd_stimuliError] = corr(pProd_stimuliError(tempBoolIndex_nonNAN),behaviorAcc(tempBoolIndex_nonNAN));
        
        %home;
        fprintf('r_decodingAcc_stimuliError=%.4f, p_decodingAcc_stimuliError=%.4f.\n',r_decodingAcc_stimuliError,p_decodingAcc_stimuliError);
        fprintf('r_pProd_stimuliError=%.4f, p_pProd_stimuliError=%.4f.\n',r_pProd_stimuliError,p_pProd_stimuliError);
    end
    a = 1;
    %% test_testingSet_error_eachResampleIter
    if if_runScript_computeLocDistriResample == 1
        fprintf('\nNow runing %s.  ------> \n', test_testingSet_error_eachResampleIter_Name_v);
        script_test_testingSet_error_eachResampleIter();
        signiSeqProportion_resample_correctFew = signiSeqProportion_resample_correct;
        signiSeqProportion_resample_stimuliError;
        signiSeqProportion_resample_responseError;
        
        signiSeqNum_resample_correctFew = signiSeqNum_resample_correct;
        signiSeqNum_resample_stimuliError;
        signiSeqNum_resample_responseError;
        
        seqNum_resample_correctFew = sum(~isnan(p_n11n_resampleMean_correct));
        seqNum_resample_stimuliError = sum(~isnan(p_n11n_resampleMean_stimuliError));
        seqNum_resample_responseError = sum(~isnan(p_n11n_resampleMean_responseError));
        
        
        fewTrialCountSeqBoolIndex_raw = fewTrialCountSeqBoolIndex;
        clear fewTrialCountSeqBoolIndex
        script_test_testingSet_error_eachResampleIter();
        signiSeqProportion_resample_correct;
        fewTrialCountSeqBoolIndex = fewTrialCountSeqBoolIndex_raw;
        
        signiSeqNum_resample_correct;
        seqNum_resample_correct = sum(~isnan(p_n11n_resampleMean_correct));
        
        signiSeqProportion_resample_correct;
        signiSeqProportion_resample_correctFew;
        signiSeqProportion_resample_stimuliError;
        signiSeqProportion_resample_responseError;
        
        %home;
        fprintf('prop(p_posterior_seq_n11n<0.05, correct)=%.1f/%d=%.3f.\n',mean(signiSeqNum_resample_correct),seqNum_resample_correct,mean(signiSeqProportion_resample_correct));
        fprintf('prop(p_posterior_seq_n11n<0.05, correctFew)=%.1f/%d=%.3f.\n',mean(signiSeqNum_resample_correctFew),seqNum_resample_correctFew,mean(signiSeqProportion_resample_correctFew));
        fprintf('prop(p_posterior_seq_n11n<0.05, stimuliError)=%.1f/%d=%.3f.\n',mean(signiSeqNum_resample_stimuliError),seqNum_resample_stimuliError,mean(signiSeqProportion_resample_stimuliError));
        fprintf('prop(p_posterior_seq_n11n<0.05, responseError)=%.1f/%d=%.3f.\n',mean(signiSeqNum_resample_responseError),seqNum_resample_responseError,mean(signiSeqProportion_resample_responseError));
    end
    
    
    %% Save
    if if_runScript_save == 1
        behaviorAcc;
        decodingAcc_correct;
        decodingAcc_correctFew;
        decodingAcc_stimuliError;
        pProd_correct;
        pProd_correctFew;
        pProd_stimuliError;
        
        signiSeqProportion_resample_correct;
        signiSeqProportion_resample_correctFew;
        signiSeqProportion_resample_stimuliError;
        signiSeqProportion_resample_responseError;
        
        a = 1;
        
        memoryPrecision_singleFOV = struct;        
        memoryPrecision_singleFOV.behaviorAcc = behaviorAcc;
        memoryPrecision_singleFOV.decodingAcc_correct = decodingAcc_correct;
        memoryPrecision_singleFOV.decodingAcc_correctFew = decodingAcc_correctFew;
        memoryPrecision_singleFOV.decodingAcc_stimuliError = decodingAcc_stimuliError;
        memoryPrecision_singleFOV.pProd_correct = pProd_correct;
        memoryPrecision_singleFOV.pProd_correctFew = pProd_correctFew;
        memoryPrecision_singleFOV.pProd_stimuliError = pProd_stimuliError;
        
        memoryPrecision_singleFOV.signiSeqProportion_resample_correct = signiSeqProportion_resample_correct;
        memoryPrecision_singleFOV.signiSeqProportion_resample_correctFew = signiSeqProportion_resample_correctFew;
        memoryPrecision_singleFOV.signiSeqProportion_resample_stimuliError = signiSeqProportion_resample_stimuliError;
        memoryPrecision_singleFOV.signiSeqProportion_resample_responseError = signiSeqProportion_resample_responseError;
        
        memoryPrecision_singleFOV.resampleTrialCount = resampleTrialCount;
        
        a = 1;
        memoryPrecision_singleFOV;
        
        temp_str = ['memoryPrecision',currentSession(14:22)];
        
        save(['C:\ASDROOT\STUDY\temp\memoryPrecision\',temp_str],'memoryPrecision_singleFOV');
    end
    
    
    
    
    
    fprintf('Time is %.1f secs.\n',toc(t0));
end

fprintf('Time of the stepB_analysis is %.1f secs.\n',toc(time0));
fprintf('num_FOV=%d.\n',num_FOV);



%% End