% Chuan's 17th script (20251214)
% This script: pipeline to conduct analyses across each FOV, related to many figures in the paper.
%% Initialization
% clear
close all

home;

if_monkey_D0_Z1 = 1;%1


if_runScript_stepF1 = 1; 
if_runScript_stepF2 = 1;

if_compute_summary = 1;

if_part2_summary = 1;

if_cellIndexMapping = 1;

if_part3_summary = 1;
if_part4_summary = 1;%1
if_part5_summary = 1;%1
if_part6_summary = 1;%


if_singleFOV0_multiFOV1 = 1;%1
currentSessionIndex_AB_summary = 8;%Ding: 8,3,4,1
% currentSessionIndex_AB_summary = 10;%Zelku: 9,1,8,11,12




if if_monkey_D0_Z1 == 0
    multiFOV_matrix_summary = ...
        [1,2;
        3,4;
        5,7;
        8,10;
        11,12;
        13,14;
        17,18;
        19,20;
        21,22;
        23,24;
        25,26;
        28,29];
    
elseif if_monkey_D0_Z1 == 1
    multiFOV_matrix_summary = ...
        [7,8;
        9,10;
        11,12;
        13,14;
        15,16;
        17,18;
        19,20;
        21,22;
        23,24;
        25,26;
        27,28;
        29,30;
        31,32;
        33,34;
        35,36;
        37,38];
end

numFOV_AB_summary = size(multiFOV_matrix_summary,1);





targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)

%% Preparation
stepF1_test_glm_locProd_eachLength_multiFOV_Name_v = autoGetFunName_myScripts('stepF1_test_glm_locProd_eachLength_multiFOV', targetPATH);
script_stepF1_test_glm_locProd_eachLength_multiFOV = str2func(stepF1_test_glm_locProd_eachLength_multiFOV_Name_v);

stepF2_test_locTuning_Name_v = autoGetFunName_myScripts('stepF2_test_locTuning', targetPATH);
script_stepF2_test_locTuning = str2func(stepF2_test_locTuning_Name_v);

stepD3_train123_test123_Name_v = autoGetFunName_myScripts('stepD3_train123_test123', targetPATH);
script_stepD3_train123_test123 = str2func(stepD3_train123_test123_Name_v);

testingSet_error_compute_Name_v = autoGetFunName_myScripts('test_testingSet_error_compute', targetPATH);
script_testingSet_error_compute = str2func(testingSet_error_compute_Name_v);

testingSet_error_eachResampleIter_Name_v = autoGetFunName_myScripts('test_testingSet_error_eachResampleIter', targetPATH);
script_testingSet_error_eachResampleIter = str2func(testingSet_error_eachResampleIter_Name_v);

stepF4_memoryMetaMismatch_singleNeuron_Name_v = autoGetFunName_myScripts('stepF4A_memoryMetaMismatch_singleNeuron', targetPATH);
script_stepF4_memoryMetaMismatch_singleNeuron = str2func(stepF4_memoryMetaMismatch_singleNeuron_Name_v);

stepF4B_test_plot_FOV_wholeImage_complexTuning_Name_v = autoGetFunName_myScripts('stepF4B_test_plot_FOV_wholeImage_complexTuning', targetPATH);
script_stepF4B_test_plot_FOV_wholeImage_complexTuning = str2func(stepF4B_test_plot_FOV_wholeImage_complexTuning_Name_v);

stepDM1_test_decodingMeta_Name_v = autoGetFunName_myScripts('stepDM1_test_decodingMeta', targetPATH);
script_stepDM1_test_decodingMeta = str2func(stepDM1_test_decodingMeta_Name_v);

stepF5A_memoryMetaMismatch_twoDecoder_Name_v = autoGetFunName_myScripts('stepF5A_memoryMetaMismatch_twoDecoder', targetPATH);
script_stepF5A_memoryMetaMismatch_twoDecoder = str2func(stepF5A_memoryMetaMismatch_twoDecoder_Name_v);

test_temp_baselineMeta_meta_GLM_Name_v = autoGetFunName_myScripts('test_temp_baselineMeta_meta_GLM', targetPATH);
script_test_temp_baselineMeta_meta_GLM = str2func(test_temp_baselineMeta_meta_GLM_Name_v);

stepF5B_test_trialHistory_Name_v = autoGetFunName_myScripts('stepF5B_test_trialHistory', targetPATH);
script_stepF5B_test_trialHistory = str2func(stepF5B_test_trialHistory_Name_v);


stepF5C_test_mismatchMechnism_bin_Name_v = autoGetFunName_myScripts('stepF5C_test_mismatchMechnism_bin', targetPATH);
script_stepF5C_test_mismatchMechnism_bin = str2func(stepF5C_test_mismatchMechnism_bin_Name_v);

stepF5D_test_mismatchMechnism_meta_crossTime_Name_v = autoGetFunName_myScripts('stepF5D_test_mismatchMechnism_meta_crossTime', targetPATH);
script_stepF5D_test_mismatchMechnism_meta_crossTime = str2func(stepF5D_test_mismatchMechnism_meta_crossTime_Name_v);

stepF5E_test_mismatchMechnism_memoryPrecision_crossTime_Name_v = autoGetFunName_myScripts('stepF5E_test_mismatchMechnism_memoryPrecision_crossTime', targetPATH);
script_stepF5E_test_mismatchMechnism_memoryPrecision_crossTime = str2func(stepF5E_test_mismatchMechnism_memoryPrecision_crossTime_Name_v);


stepD4_onlyLoad_Name_v = autoGetFunName_myScripts('stepD4_onlyLoad', targetPATH);
script_stepD4_onlyLoad = str2func(stepD4_onlyLoad_Name_v);



t0_runScript = tic;



%% Compute: Step F1
if if_runScript_stepF1 == 1
    fprintf('\nNow runing %s.  ------> \n', stepF1_test_glm_locProd_eachLength_multiFOV_Name_v);
    script_stepF1_test_glm_locProd_eachLength_multiFOV();
end

%% Compute: Step F2
if if_runScript_stepF2 == 1
    fprintf('\nNow runing %s.  ------> \n', stepF2_test_locTuning_Name_v);
    script_stepF2_test_locTuning();
end


part2_summary_mulFOV = struct;
part3_summary_mulFOV = struct;
part4_summary_mulFOV = struct;
part5_summary_mulFOV = struct;
part6_summary_mulFOV = struct;

%part2_summary_mulFOV
part2_summary_mulFOV.r_precision_summary = [];
part2_summary_mulFOV.p_precision_summary = [];

part2_summary_mulFOV.r_precision_allMemory_summary = [];
part2_summary_mulFOV.p_precision_allMemory_summary = [];
part2_summary_mulFOV.r_precision_allNoChoice_summary = [];
part2_summary_mulFOV.p_precision_allNoChoice_summary = [];

part2_summary_mulFOV.p12_locationCorr_summary = [];
part2_summary_mulFOV.p13_locationCorr_summary = [];
part2_summary_mulFOV.p23_locationCorr_summary = [];
part2_summary_mulFOV.r_mean_correct_summary = [];
part2_summary_mulFOV.r_mean_stimuliError_summary = [];
part2_summary_mulFOV.r_mean_responseError_summary = [];
part2_summary_mulFOV.r_mean_chanceLevel_summary = [];

part2_summary_mulFOV.r_mean_correct_theoretical_summary = [];

part2_summary_mulFOV.seqAccuracy_memoryPrecisionLow_allMemoryStimuli_mean_summary = [];
part2_summary_mulFOV.seqAccuracy_memoryPrecisionHigh_allMemoryStimuli_mean_summary = [];
part2_summary_mulFOV.p_precisionLowHigh_stimuli_summary = [];
part2_summary_mulFOV.seqAccuracy_memoryPrecisionLow_allMemoryResponse_mean_summary = [];
part2_summary_mulFOV.seqAccuracy_memoryPrecisionHigh_allMemoryResponse_mean_summary = [];
part2_summary_mulFOV.p_precisionLowHigh_response_summary = [];

part2_summary_mulFOV.memoryPrecision_allSeq_allMemoryCorrectMean_mean_summary = [];
part2_summary_mulFOV.memoryPrecision_allSeq_allMemoryStimuliErrorMean_mean_summary = [];
part2_summary_mulFOV.memoryPrecision_allSeq_allMemoryResponseErrorMean_mean_summary = [];
part2_summary_mulFOV.p_correctErrorPrecision_stimuli_summary = [];
part2_summary_mulFOV.p_correctErrorPrecision_response_summary = [];


%part3_summary_mulFOV
part3_summary_mulFOV.AUROC_meta_delay1_summary = [];
part3_summary_mulFOV.r_meta_seqLevel_choice_delay1_summary = [];
part3_summary_mulFOV.p_meta_seqLevel_choice_delay1_summary = [];
part3_summary_mulFOV.r_meta_seqLevel_noChoice_delay1_summary = [];
part3_summary_mulFOV.p_meta_seqLevel_noChoice_delay1_summary = [];

part3_summary_mulFOV.r_meta_seqLevel_choiceMemory_delay1_summary = [];
part3_summary_mulFOV.r_meta_seqLevel_choiceOffload_delay1_summary = [];

part3_summary_mulFOV.seqRProb_metaLow_choice_mean_summary = [];
part3_summary_mulFOV.seqRProb_metaHigh_choice_mean_summary = [];
part3_summary_mulFOV.p_metaLowHigh_summary = [];

part3_summary_mulFOV.meta_allSeq_choiceOffloadMean_mean_summary = [];
part3_summary_mulFOV.meta_allSeq_choiceMemoryStimuliMean_mean_summary = [];
part3_summary_mulFOV.p_choiceMeta_stimuli_summary = [];

part3_summary_mulFOV.r2_seqPrecision_summary = [];
part3_summary_mulFOV.r2_rProb_summary = [];

part3_summary_mulFOV.r2_seqPrecision_cell_summary = [];
part3_summary_mulFOV.r2_rProb_cell_summary = [];

part3_summary_mulFOV.r_precisionMeta_seqLevel_summary = [];
part3_summary_mulFOV.p_precisionMeta_seqLevel_summary = [];
part3_summary_mulFOV.seqRProb_memoryPrecisionLow_choiceResponse_mean_summary = [];
part3_summary_mulFOV.seqRProb_memoryPrecisionHigh_choiceResponse_mean_summary = [];
part3_summary_mulFOV.p_seqRProb_precisionLowHigh_response_summary = [];

part3_summary_mulFOV.memoryPrecision_allSeq_choiceOffloadMean_mean_summary = [];
part3_summary_mulFOV.memoryPrecision_allSeq_choiceMemoryResponseMean_mean_summary = [];
part3_summary_mulFOV.p_choicePrecision_response_summary = [];

part3_summary_mulFOV.memoryPrecision_allSeq_choiceMemoryStimuliMean_mean_summary = [];
part3_summary_mulFOV.p_choicePrecision_stimuli_summary = [];

part3_summary_mulFOV.seqMeta_memoryPrecisionLow_choiceStimuli_mean_summary = [];
part3_summary_mulFOV.seqMeta_memoryPrecisionHigh_choiceStimuli_mean_summary = [];
part3_summary_mulFOV.p_seqMeta_precisionLowHigh_stimuli_summary = [];
part3_summary_mulFOV.seqMeta_memoryPrecisionLow_choiceResponse_mean_summary = [];
part3_summary_mulFOV.seqMeta_memoryPrecisionHigh_choiceResponse_mean_summary = [];
part3_summary_mulFOV.p_seqMeta_precisionLowHigh_response_summary = [];

part3_summary_mulFOV.seqCorr_precisionMeta_choiceStimuli_mean_summary = [];
part3_summary_mulFOV.seqCorr_precisionMeta_choiceResponse_mean_summary = [];

part3_summary_mulFOV.trialCorr_precisionMeta_choice_summary = [];
part3_summary_mulFOV.trialCorr_precisionMeta_choiceMemory_summary = [];
part3_summary_mulFOV.trialCorr_precisionMeta_choiceOffload_summary = [];

%part4_summary
part4_summary_mulFOV.trialNum_memoryPrecisionLow_metaLow_choice_summary = [];
part4_summary_mulFOV.trialNum_memoryPrecisionHigh_metaHigh_choice_summary = [];
part4_summary_mulFOV.trialNum_memoryPrecisionLow_metaHigh_choice_summary = [];
part4_summary_mulFOV.trialNum_memoryPrecisionHigh_metaLow_choice_summary = [];

part4_summary_mulFOV.signifTimeStampCount_meta_mean_summary = [];
part4_summary_mulFOV.signifTimeStampCount_memory_mean_summary = [];
part4_summary_mulFOV.overTrialProp_resampled_mean_summary = [];
part4_summary_mulFOV.overTrialProp_resampled_median_summary = [];
part4_summary_mulFOV.overTrialProp_summary = [];


part4_summary_mulFOV.temptempBoolIndex_CMC_summary = [];
part4_summary_mulFOV.temptempBoolIndex_CF_summary = [];
part4_summary_mulFOV.temptempBoolIndex_CME_summary = [];

part4_summary_mulFOV.proj1_summary = [];
part4_summary_mulFOV.proj2_summary = [];
part4_summary_mulFOV.proj12_summary = [];
part4_summary_mulFOV.proj3_summary = [];
        
%part5_summary
part5_summary_mulFOV.AUROC_meta_baseline_summary = [];
part5_summary_mulFOV.p_baseline_choiceMemory_choiceOffload_summary = [];
part5_summary_mulFOV.r_meta_seqLevel_choice_baseline_summary = [];
part5_summary_mulFOV.p_meta_seqLevel_choice_baseline_summary = [];
part5_summary_mulFOV.r_meta_seqLevel_noChoice_baseline_summary = [];
part5_summary_mulFOV.p_meta_seqLevel_noChoice_baseline_summary = [];


part5_summary_mulFOV.baselineMeta_highMatch_summary = [];
part5_summary_mulFOV.baselineMeta_lowMatch_summary = [];
part5_summary_mulFOV.baselineMeta_overMismatch_summary = [];
part5_summary_mulFOV.baselineMeta_underMismatch_summary = [];
part5_summary_mulFOV.p_baselineMeta_highMatch_underMismatch_summary = [];
part5_summary_mulFOV.p_baselineMeta_lowMatch_overMismatch_summary = [];

part5_summary_mulFOV.baselineMeta_highMeta_summary = [];
part5_summary_mulFOV.baselineMeta_lowMeta_summary = [];
part5_summary_mulFOV.p_baselineMeta_highMeta_lowMeta_summary = [];


part5_summary_mulFOV.meta_trialLevel_baseline_summary = [];
part5_summary_mulFOV.memoryPrecision_trialLevel_summary = [];
part5_summary_mulFOV.meta_trialLevel_summary = [];

part5_summary_mulFOV.meta_trialLevel_baseline_cell_summary = [];
part5_summary_mulFOV.memoryPrecision_trialLevel_cell_summary = [];
part5_summary_mulFOV.meta_trialLevel_cell_summary = [];

part5_summary_mulFOV.r2_metaRegress_caseA_resampled_mean_summary = [];
part5_summary_mulFOV.r2_metaRegress_caseB_resampled_mean_summary = [];
part5_summary_mulFOV.r2_metaRegress_caseC_resampled_mean_summary = [];
part5_summary_mulFOV.pAB_metaRegress_summary = [];
part5_summary_mulFOV.pAC_metaRegress_summary = [];
part5_summary_mulFOV.pBC_metaRegress_summary = [];

part5_summary_mulFOV.beta_precision_metaRegress_caseC_resampled_mean_summary = [];
part5_summary_mulFOV.beta_baseline_metaRegress_caseC_resampled_mean_summary = [];
part5_summary_mulFOV.beta_interaction_metaRegress_caseC_resampled_mean_summary = [];


part5_summary_mulFOV.r2_metaRegress_caseE2_CM_resampled_mean_summary = [];
part5_summary_mulFOV.beta_precision_metaRegress_caseE2_CM_resampled_mean_summary = [];
part5_summary_mulFOV.beta_baseline_metaRegress_caseE2_CM_resampled_mean_summary = [];
part5_summary_mulFOV.beta_interaction_metaRegress_caseE2_CM_resampled_mean_summary = [];

part5_summary_mulFOV.r2_metaRegress_caseF_CF_resampled_mean_summary = [];
part5_summary_mulFOV.beta_precision_metaRegress_caseF_CF_resampled_mean_summary = [];
part5_summary_mulFOV.beta_baseline_metaRegress_caseF_CF_resampled_mean_summary = [];
part5_summary_mulFOV.beta_interaction_metaRegress_caseF_CF_CM_resampled_mean_summary = [];

part5_summary_mulFOV.r2_mismatchRegress_caseD2B_baseline_resampled_mean_summary = [];
part5_summary_mulFOV.r2_mismatchRegress_caseD2A_precision_resampled_mean_summary = [];


part5_summary_mulFOV.memoryPrecision_CME_baseline_mean_summary = [];
part5_summary_mulFOV.memoryPrecision_CF_baseline_mean_summary = [];
part5_summary_mulFOV.meta_CME_baseline_mean_summary = [];
part5_summary_mulFOV.meta_CF_baseline_mean_summary = [];

part5_summary_mulFOV.memoryPrecision_CME_baseline_median_summary = [];
part5_summary_mulFOV.memoryPrecision_CF_baseline_median_summary = [];
part5_summary_mulFOV.meta_CME_baseline_median_summary = [];
part5_summary_mulFOV.meta_CF_baseline_median_summary = [];


part5_summary_mulFOV.AUROC_optimal_trialHistory_summary = [];
part5_summary_mulFOV.normStd_optimal_trialHistory_summary = [];

part5_summary_mulFOV.weight_trialHistory_summary = [];

part5_summary_mulFOV.p_historyReward_baselineLowHigh_summary = [];
part5_summary_mulFOV.historyRewardMean_baselineHigh_summary = [];
part5_summary_mulFOV.historyRewardMean_baselineLow_summary = [];

part5_summary_mulFOV.p_historyReward_metaLowHigh_summary = [];
part5_summary_mulFOV.historyRewardMean_metaHigh_summary = [];
part5_summary_mulFOV.historyRewardMean_metaLow_summary = [];

part5_summary_mulFOV.p_reward_highMatch_underMismatch_summary = [];
part5_summary_mulFOV.p_reward_lowMatch_overMismatch_summary = [];

part5_summary_mulFOV.p_linear_historyReward_4Types_summary = [];
part5_summary_mulFOV.historyRewardMean_HighMatch_summary = [];
part5_summary_mulFOV.historyRewardMean_LowMatch_summary = [];
part5_summary_mulFOV.historyRewardMean_OverMismatch_summary = [];
part5_summary_mulFOV.historyRewardMean_UnderMismatch_summary = [];
        
%part6_summary_mulFOV
part6_summary_mulFOV.beta_memoryPrecision_summary = [];
part6_summary_mulFOV.beta_choiceMemory_summary = [];
part6_summary_mulFOV.beta_choiceMemory_baseline_summary = [];

part6_summary_mulFOV.beta_memoryPrecision_cell_summary = [];
part6_summary_mulFOV.beta_choiceMemory_cell_summary = [];
part6_summary_mulFOV.beta_choiceMemory_baseline_cell_summary = [];

part6_summary_mulFOV.beta_seqPrecision_summary = [];
part6_summary_mulFOV.beta_gProb_summary = [];
part6_summary_mulFOV.beta_gProb_baseline_summary = [];

part6_summary_mulFOV.beta_seqPrecision_cell_summary = [];
part6_summary_mulFOV.beta_gProb_cell_summary = [];
part6_summary_mulFOV.beta_gProb_baseline_cell_summary = [];

part6_summary_mulFOV.r2_memoryPrecision_summary = [];
part6_summary_mulFOV.r2_choiceMemory_summary = [];
part6_summary_mulFOV.r2_choiceMemory_baseline_summary = [];
part6_summary_mulFOV.r2_seqPrecision_summary = [];
part6_summary_mulFOV.r2_gProb_summary = [];
part6_summary_mulFOV.r2_gProb_baseline_summary = [];
part6_summary_mulFOV.beta_6loc_summary = [];
part6_summary_mulFOV.r2_6loc_summary = [];
     
part6_summary_mulFOV.r2_memoryPrecision_cell_summary = [];
part6_summary_mulFOV.r2_choiceMemory_cell_summary = [];
part6_summary_mulFOV.r2_choiceMemory_baseline_cell_summary = [];
part6_summary_mulFOV.r2_seqPrecision_cell_summary = [];
part6_summary_mulFOV.r2_gProb_cell_summary = [];
part6_summary_mulFOV.r2_gProb_baseline_cell_summary = [];
part6_summary_mulFOV.beta_6loc_cell_summary = [];
part6_summary_mulFOV.r2_6loc_cell_summary = [];


part6_summary_mulFOV.betaAngle_6loc_halfAB_norm_unitView_summary = [];
part6_summary_mulFOV.betaAngle_gProb_halfAB_norm_unitView_summary = [];
part6_summary_mulFOV.betaAngle_choiceMemory_baseline_halfAB_norm_unitView_summary = [];


part6_summary_mulFOV.beta_prior_linearAxis_summary = [];
part6_summary_mulFOV.beta_memory_linearAxis_summary = [];
part6_summary_mulFOV.beta_meta_linearAxis_summary = [];

part6_summary_mulFOV.r2_prior_linearAxis_summary = [];
part6_summary_mulFOV.r2_memory_linearAxis_summary = [];
part6_summary_mulFOV.r2_meta_linearAxis_summary = [];

        
part6_summary_mulFOV.temp_degree_aa_prior_summary = [];
part6_summary_mulFOV.temp_degree_bb_memory_summary = [];
part6_summary_mulFOV.temp_degree_cc_meta_summary = [];

part6_summary_mulFOV.temp_degree_ab_priorMemory_summary = [];
part6_summary_mulFOV.temp_degree_ac_priorMeta_summary = [];
part6_summary_mulFOV.temp_degree_bc_memoryMeta_summary = [];

part6_summary_mulFOV.temp_VAF_ratio_priorPrior_summary = [];
part6_summary_mulFOV.temp_VAF_ratio_memoryMemory_summary = [];
part6_summary_mulFOV.temp_VAF_ratio_metaMeta_summary = [];
part6_summary_mulFOV.temp_VAF_ratio_priorMemory_summary = [];
part6_summary_mulFOV.temp_VAF_ratio_priorMeta_summary = [];
part6_summary_mulFOV.temp_VAF_ratio_memoryMeta_summary = [];


part6_summary_mulFOV.temp_EV_neuronPrior_summary = [];
part6_summary_mulFOV.temp_EV_neuronMemory_summary = [];
part6_summary_mulFOV.temp_EV_neuronMeta_summary = [];
part6_summary_mulFOV.pca_explained_baseline_cell_summary = [];
part6_summary_mulFOV.pca_explained_delay_cell_summary = [];
part6_summary_mulFOV.pca_r_eachPC_VS_eventVector_baseline_cell_summary = [];
part6_summary_mulFOV.pca_r_eachPC_VS_eventVector_delay_cell_summary = [];


part6_summary_mulFOV.clusterIndex_memoryPrecision_disBin_mean_summary = [];
part6_summary_mulFOV.clusterIndex_choiceMemory_disBin_mean_summary = [];
part6_summary_mulFOV.clusterIndex_choiceMemory_baseline_disBin_mean_summary = [];
part6_summary_mulFOV.prctile_high_memoryPrecision_shuffled_summary = [];
part6_summary_mulFOV.prctile_high_choiceMemory_shuffled_summary = [];
part6_summary_mulFOV.prctile_high_choiceMemory_baseline_shuffled_summary = [];

part6_summary_mulFOV.clusterIndex_memoryPrecision_shuffled_meanB_summary = [];
part6_summary_mulFOV.clusterIndex_choiceMemory_shuffled_meanB_summary = [];
part6_summary_mulFOV.clusterIndex_choiceMemory_baseline_shuffled_meanB_summary = [];

part6_summary_mulFOV.centriodDis_all_shuffled_prctileA_summary = [];
part6_summary_mulFOV.centriodDis_all_shuffled_prctileB_summary = [];

part6_summary_mulFOV.centriodDis_all_shuffled_mean_summary = [];

part6_summary_mulFOV.centriodDis_all_summary = [];
part6_summary_mulFOV.centriodDis_all_AC_shuffled_prctileA_summary = [];
part6_summary_mulFOV.centriodDis_all_AC_shuffled_prctileB_summary = [];

part6_summary_mulFOV.centriodDis_all_AC_shuffled_mean_summary = [];

part6_summary_mulFOV.centriodDis_AC_all_summary = [];
part6_summary_mulFOV.centriodDis_all_BC_shuffled_prctileA_summary = [];
part6_summary_mulFOV.centriodDis_all_BC_shuffled_prctileB_summary = [];

part6_summary_mulFOV.centriodDis_all_BC_shuffled_mean_summary = [];

part6_summary_mulFOV.centriodDis_BC_all_summary = [];


part6_summary_mulFOV.spatialDisNew_1to1_mean_summary = [];
part6_summary_mulFOV.spatialDisNew_2to2_mean_summary = [];
part6_summary_mulFOV.spatialDisNew_3to3_mean_summary = [];
part6_summary_mulFOV.spatialDisNew_1toOthers_mean_summary = [];
part6_summary_mulFOV.spatialDisNew_2toOthers_mean_summary = [];
part6_summary_mulFOV.spatialDisNew_3toOthers_mean_summary = [];
part6_summary_mulFOV.spatialCentriodDis_1toOthers_summary = [];
part6_summary_mulFOV.spatialCentriodDis_2toOthers_summary = [];
part6_summary_mulFOV.spatialCentriodDis_3toOthers_summary = [];


%% 
a = 1;
b = 2;

temp_cellIndexMapping_AB_mulFOV = [];
p_6loc_mulFOV = [];
p_memoryPrecision_mulFOV = [];
p_choiceMemory_mulFOV = [];
p_choiceMemory_baseline_mulFOV = [];

p_seqPrecision_mulFOV = [];
p_gProb_mulFOV = [];
p_gProb_baseline_mulFOV = [];
        
%%        
for tempIndexFOV_AB=1:numFOV_AB_summary
% for tempIndexFOV_AB=2:numFOV_AB_summary
  
    if if_singleFOV0_multiFOV1 == 0
        tempIndexFOV_AB = currentSessionIndex_AB_summary; %#ok<*FXSET>        
        fprintf('Single FOV mode (AB), id = %d.\n',tempIndexFOV_AB);
    elseif if_singleFOV0_multiFOV1 == 1
        currentSessionIndex_AB_summary = tempIndexFOV_AB;      
    end
    
    %% Compute: Part 2
    if if_part2_summary == 1
        %% Compute: Step D3
        if if_compute_summary == 1
            fprintf('\nNow runing %s.  ------> \n', stepD3_train123_test123_Name_v);
            script_stepD3_train123_test123();
        end
        
        
        %% Compute: testingSet_error
        if if_compute_summary == 1
            fprintf('\nNow runing %s.  ------> \n', testingSet_error_compute_Name_v);
            script_testingSet_error_compute();
        end
        
        %% Compute: testingSet_error_eachResampleIter
        if if_compute_summary == 1
            fprintf('\nNow runing %s.  ------> \n', testingSet_error_eachResampleIter_Name_v);
            script_testingSet_error_eachResampleIter();
        end
        
        %% Compute: Step F4
        if if_compute_summary == 1
            fprintf('\nNow runing %s.  ------> \n', stepF4_memoryMetaMismatch_singleNeuron_Name_v);
            
            if if_memoryPrecision_accuracy0_sigma1 == 0
                if_memeoryPrecision_stimuli0_response1_summary = 0; %#ok<*NASGU>
                script_stepF4_memoryMetaMismatch_singleNeuron();
                
                if_memeoryPrecision_stimuli0_response1_summary = 1;
                script_stepF4_memoryMetaMismatch_singleNeuron();
                
            elseif if_memoryPrecision_accuracy0_sigma1 == 1
                if_memeoryPrecision_stimuli0_response1_summary = 1;
                script_stepF4_memoryMetaMismatch_singleNeuron();
                
                if_memeoryPrecision_stimuli0_response1_summary = 0;
                script_stepF4_memoryMetaMismatch_singleNeuron();
            end

            %if_memeoryPrecision_stimuli0_response1_summary = 0;
            %script_stepF4_memoryMetaMismatch_singleNeuron();
            
        end
        
        part2_summary = struct;
    end
    
    
    %% Part 2H: Sequence-level memory precision evidence
    if if_part2_summary == 1
        % memory correct trials
        temp_x = seqPrecision_neuron;
        temp_y = seqPrecision_behavior;
        tempBoolIndex = ~isnan(temp_x);
        
        x = temp_x(~isnan(temp_x));
        y = temp_y(1:sum(numSeq(1:valid_length)));
        y = y(~isnan(temp_x));
        
        [r_precision_summary,p_precision_summary] = corr(x,y);
        
        % all (correct + error) memory trials
        temp_x = temp_seqPrecision_neuron;
        temp_y = seqPrecision_behavior;
        tempBoolIndex = ~isnan(temp_x);
        
        x = temp_x(~isnan(temp_x));
        y = temp_y(1:sum(numSeq(1:valid_length)));
        y = y(~isnan(temp_x));
        
        [r_precision_allMemory_summary,p_precision_allMemory_summary] = corr(x,y);
        
        
        % all (correct + error) noChoice trials
        temp_x = temp2_seqPrecision_neuron;
        temp_y = seqPrecision_behavior;
        tempBoolIndex = ~isnan(temp_x);
        
        x = temp_x(~isnan(temp_x));
        y = temp_y(1:sum(numSeq(1:valid_length)));
        y = y(~isnan(temp_x));
        
        [r_precision_allNoChoice_summary,p_precision_allNoChoice_summary] = corr(x,y);
        
        part2_summary.r_precision_summary = r_precision_summary;
        part2_summary.p_precision_summary = p_precision_summary;

        part2_summary.r_precision_allMemory_summary = r_precision_allMemory_summary;
        part2_summary.p_precision_allMemory_summary = p_precision_allMemory_summary;
        
        part2_summary.r_precision_allNoChoice_summary = r_precision_allNoChoice_summary;
        part2_summary.p_precision_allNoChoice_summary = p_precision_allNoChoice_summary;        
    end
    
    
    %% Part 2G: Location correlation of correct, stimuliError and responseError trials
    if if_part2_summary == 1
        [~,p12_locationCorr_summary,~,~] = ttest(r_n11n_resampleMean_correct,r_n11n_resampleMean_stimuliError);
        [~,p13_locationCorr_summary,~,~] = ttest(r_n11n_resampleMean_correct,r_n11n_resampleMean_responseError);
        [~,p23_locationCorr_summary,~,~] = ttest(r_n11n_resampleMean_stimuliError,r_n11n_resampleMean_responseError);
        
        
        r_mean_correct_summary = mean(r_n11n_resampleMean_correct,'omitnan');
        r_mean_stimuliError_summary = mean(r_n11n_resampleMean_stimuliError,'omitnan');
        r_mean_responseError_summary = mean(r_n11n_resampleMean_responseError,'omitnan');
        r_mean_chanceLevel_summary = mean(r_n11n_resampleMean_correct_chanceLevel,'omitnan');
        r_mean_correct_theoretical_summary = mean(r_n11n_resampleMean_correct_theoretical,'omitnan');
        
        part2_summary.p12_locationCorr_summary = p12_locationCorr_summary;
        part2_summary.p13_locationCorr_summary = p13_locationCorr_summary;
        part2_summary.p23_locationCorr_summary = p23_locationCorr_summary;
        part2_summary.r_mean_correct_summary = r_mean_correct_summary;
        part2_summary.r_mean_stimuliError_summary = r_mean_stimuliError_summary;
        part2_summary.r_mean_responseError_summary = r_mean_responseError_summary;
        part2_summary.r_mean_chanceLevel_summary = r_mean_chanceLevel_summary;
        part2_summary.r_mean_correct_theoretical_summary = r_mean_correct_theoretical_summary;
        
    end
    
    
    
    %% Part 2I, Right: Trial-level memory precision evidence
    if if_part2_summary == 1
        b1 = seqAccuracy_memoryPrecisionLow_allMemoryStimuli;
        b2 = seqAccuracy_memoryPrecisionHigh_allMemoryStimuli;
        
        b3 = seqAccuracy_memoryPrecisionLow_allMemoryResponse;
        b4 = seqAccuracy_memoryPrecisionHigh_allMemoryResponse;
        
        seqAccuracy_memoryPrecisionLow_allMemoryStimuli_mean_summary = mean(b1,'omitnan');
        seqAccuracy_memoryPrecisionHigh_allMemoryStimuli_mean_summary = mean(b2,'omitnan');
        seqAccuracy_memoryPrecisionLow_allMemoryResponse_mean_summary = mean(b3,'omitnan');
        seqAccuracy_memoryPrecisionHigh_allMemoryResponse_mean_summary = mean(b4,'omitnan');
        
        [~,p_precisionLowHigh_stimuli_summary] = ttest(b1,b2);
        [~,p_precisionLowHigh_response_summary] = ttest(b3,b4);
        
        part2_summary.seqAccuracy_memoryPrecisionLow_allMemoryStimuli_mean_summary = seqAccuracy_memoryPrecisionLow_allMemoryStimuli_mean_summary;
        part2_summary.seqAccuracy_memoryPrecisionHigh_allMemoryStimuli_mean_summary = seqAccuracy_memoryPrecisionHigh_allMemoryStimuli_mean_summary;
        part2_summary.p_precisionLowHigh_stimuli_summary = p_precisionLowHigh_stimuli_summary;
        part2_summary.seqAccuracy_memoryPrecisionLow_allMemoryResponse_mean_summary = seqAccuracy_memoryPrecisionLow_allMemoryResponse_mean_summary;
        part2_summary.seqAccuracy_memoryPrecisionHigh_allMemoryResponse_mean_summary = seqAccuracy_memoryPrecisionHigh_allMemoryResponse_mean_summary;
        part2_summary.p_precisionLowHigh_response_summary = p_precisionLowHigh_response_summary;
        
    end
    
    
    %% (New) Part 2I, Right: Trial-level memory precision evidence
    if if_part2_summary == 1
        
        d1 = memoryPrecision_trialLevel_allSeq_allMemoryCorrectMean;
        d2 = memoryPrecision_trialLevel_allSeq_allMemoryStimuliErrorMean;
        d3 = memoryPrecision_trialLevel_allSeq_allMemoryResponseErrorMean;
                
        memoryPrecision_allSeq_allMemoryCorrectMean_mean_summary = mean(d1,'omitnan');
        memoryPrecision_allSeq_allMemoryStimuliErrorMean_mean_summary = mean(d2,'omitnan');
        memoryPrecision_allSeq_allMemoryResponseErrorMean_mean_summary = mean(d3,'omitnan');
        
        [~,p_correctErrorPrecision_stimuli_summary] = ttest(d1,d2);
        [~,p_correctErrorPrecision_response_summary] = ttest(d1,d3);
        
        
        part2_summary.memoryPrecision_allSeq_allMemoryCorrectMean_mean_summary = memoryPrecision_allSeq_allMemoryCorrectMean_mean_summary;
        part2_summary.memoryPrecision_allSeq_allMemoryStimuliErrorMean_mean_summary = memoryPrecision_allSeq_allMemoryStimuliErrorMean_mean_summary;
        part2_summary.memoryPrecision_allSeq_allMemoryResponseErrorMean_mean_summary = memoryPrecision_allSeq_allMemoryResponseErrorMean_mean_summary;
        part2_summary.p_correctErrorPrecision_stimuli_summary = p_correctErrorPrecision_stimuli_summary;
        part2_summary.p_correctErrorPrecision_response_summary = p_correctErrorPrecision_response_summary;
                
    end    
    
    %% Summary: Part 2
    if if_part2_summary == 1
        part2_summary_mulFOV.r_precision_summary = [part2_summary_mulFOV.r_precision_summary r_precision_summary];
        part2_summary_mulFOV.p_precision_summary = [part2_summary_mulFOV.p_precision_summary p_precision_summary];
  
        part2_summary_mulFOV.r_precision_allMemory_summary = [part2_summary_mulFOV.r_precision_allMemory_summary r_precision_allMemory_summary];
        part2_summary_mulFOV.p_precision_allMemory_summary = [part2_summary_mulFOV.p_precision_allMemory_summary p_precision_allMemory_summary];
        part2_summary_mulFOV.r_precision_allNoChoice_summary = [part2_summary_mulFOV.r_precision_allNoChoice_summary r_precision_allNoChoice_summary];
        part2_summary_mulFOV.p_precision_allNoChoice_summary = [part2_summary_mulFOV.p_precision_allNoChoice_summary p_precision_allNoChoice_summary];
        
        part2_summary_mulFOV.p12_locationCorr_summary = [part2_summary_mulFOV.p12_locationCorr_summary p12_locationCorr_summary];
        part2_summary_mulFOV.p13_locationCorr_summary = [part2_summary_mulFOV.p13_locationCorr_summary p13_locationCorr_summary];
        part2_summary_mulFOV.p23_locationCorr_summary = [part2_summary_mulFOV.p23_locationCorr_summary p23_locationCorr_summary];
        part2_summary_mulFOV.r_mean_correct_summary = [part2_summary_mulFOV.r_mean_correct_summary r_mean_correct_summary];
        part2_summary_mulFOV.r_mean_stimuliError_summary = [part2_summary_mulFOV.r_mean_stimuliError_summary r_mean_stimuliError_summary];
        part2_summary_mulFOV.r_mean_responseError_summary = [part2_summary_mulFOV.r_mean_responseError_summary r_mean_responseError_summary];
        part2_summary_mulFOV.r_mean_chanceLevel_summary = [part2_summary_mulFOV.r_mean_chanceLevel_summary r_mean_chanceLevel_summary];
        
        part2_summary_mulFOV.r_mean_correct_theoretical_summary = [part2_summary_mulFOV.r_mean_correct_theoretical_summary r_mean_correct_theoretical_summary];        
        
        part2_summary_mulFOV.seqAccuracy_memoryPrecisionLow_allMemoryStimuli_mean_summary = [part2_summary_mulFOV.seqAccuracy_memoryPrecisionLow_allMemoryStimuli_mean_summary seqAccuracy_memoryPrecisionLow_allMemoryStimuli_mean_summary];
        part2_summary_mulFOV.seqAccuracy_memoryPrecisionHigh_allMemoryStimuli_mean_summary = [part2_summary_mulFOV.seqAccuracy_memoryPrecisionHigh_allMemoryStimuli_mean_summary seqAccuracy_memoryPrecisionHigh_allMemoryStimuli_mean_summary];
        part2_summary_mulFOV.p_precisionLowHigh_stimuli_summary = [part2_summary_mulFOV.p_precisionLowHigh_stimuli_summary p_precisionLowHigh_stimuli_summary];
        part2_summary_mulFOV.seqAccuracy_memoryPrecisionLow_allMemoryResponse_mean_summary = [part2_summary_mulFOV.seqAccuracy_memoryPrecisionLow_allMemoryResponse_mean_summary seqAccuracy_memoryPrecisionLow_allMemoryResponse_mean_summary];
        part2_summary_mulFOV.seqAccuracy_memoryPrecisionHigh_allMemoryResponse_mean_summary = [part2_summary_mulFOV.seqAccuracy_memoryPrecisionHigh_allMemoryResponse_mean_summary seqAccuracy_memoryPrecisionHigh_allMemoryResponse_mean_summary];
        part2_summary_mulFOV.p_precisionLowHigh_response_summary = [part2_summary_mulFOV.p_precisionLowHigh_response_summary p_precisionLowHigh_response_summary];        
        
        part2_summary_mulFOV.memoryPrecision_allSeq_allMemoryCorrectMean_mean_summary = [part2_summary_mulFOV.memoryPrecision_allSeq_allMemoryCorrectMean_mean_summary memoryPrecision_allSeq_allMemoryCorrectMean_mean_summary];
        part2_summary_mulFOV.memoryPrecision_allSeq_allMemoryStimuliErrorMean_mean_summary = [part2_summary_mulFOV.memoryPrecision_allSeq_allMemoryStimuliErrorMean_mean_summary memoryPrecision_allSeq_allMemoryStimuliErrorMean_mean_summary];
        part2_summary_mulFOV.memoryPrecision_allSeq_allMemoryResponseErrorMean_mean_summary = [part2_summary_mulFOV.memoryPrecision_allSeq_allMemoryResponseErrorMean_mean_summary memoryPrecision_allSeq_allMemoryResponseErrorMean_mean_summary];
        part2_summary_mulFOV.p_correctErrorPrecision_stimuli_summary = [part2_summary_mulFOV.p_correctErrorPrecision_stimuli_summary p_correctErrorPrecision_stimuli_summary];
        part2_summary_mulFOV.p_correctErrorPrecision_response_summary = [part2_summary_mulFOV.p_correctErrorPrecision_response_summary p_correctErrorPrecision_response_summary];
        
    end

    %% CellIndexMapping
    if if_cellIndexMapping == 1
        
        if if_part2_summary == 0
            %% Compute: Step D3
            if if_compute_summary == 1
                fprintf('\nNow runing %s.  ------> \n', stepD4_onlyLoad_Name_v);
                script_stepD4_onlyLoad();
            end
        end
        
        decodingDataSimplified;
        currentSessionIndex_B;
        
        
        cellIndex_suite2p = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,1);
        
        
        temp_range_B = FOVAllCellRange_multiFOV(currentSessionIndex_B,1):FOVAllCellRange_multiFOV(currentSessionIndex_B,2);
        cellIndex_suite2p_B_raw = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(temp_range_B);
        cellIndex_suite2p_B = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,end);
        temp_range_AB = nan(1,length(cellIndex_suite2p_B));
        for tempi=1:length(cellIndex_suite2p_B)
            temp_range_AB(tempi) = find(cellIndex_suite2p_B_raw==cellIndex_suite2p_B(tempi));
        end
        temp_range_AB = temp_range_AB + temp_range_B(1) - 1;
        
        temp_cellIndexMapping_AB = temp_range_AB';
        
        temp_cellIndexMapping_AB_mulFOV = [temp_cellIndexMapping_AB_mulFOV;{temp_cellIndexMapping_AB}]; %#ok<*AGROW>
        
        
        
        p_memoryPrecision;
        p_choiceMemory;
        p_choiceMemory_baseline;
        
        p_seqPrecision;
        p_gProb;
        p_gProb_baseline;
        
        
        p_6loc_mulFOV = [p_6loc_mulFOV;{p_6loc}];        
        p_memoryPrecision_mulFOV = [p_memoryPrecision_mulFOV;{p_memoryPrecision}];
        p_choiceMemory_mulFOV = [p_choiceMemory_mulFOV;{p_choiceMemory}];
        p_choiceMemory_baseline_mulFOV = [p_choiceMemory_baseline_mulFOV;{p_choiceMemory_baseline}];

        p_seqPrecision_mulFOV = [p_seqPrecision_mulFOV;{p_seqPrecision}];
        p_gProb_mulFOV = [p_gProb_mulFOV;{p_gProb}];
        p_gProb_baseline_mulFOV = [p_gProb_baseline_mulFOV;{p_gProb_baseline}];
        

    end
    
    %% Compute: Part 3
    if if_part3_summary == 1
        %% Compute: Step DM1
        if if_compute_summary == 1
            fprintf('\nNow runing %s.  ------> \n', stepDM1_test_decodingMeta_Name_v);
            
            if_trainMeta_0baseline_1delay1 = 0; %#ok<*NASGU>
            script_stepDM1_test_decodingMeta();
            
            if_trainMeta_0baseline_1delay1 = 1;
            script_stepDM1_test_decodingMeta();
        end
        
        %% Compute: Step F4
        % no need
        
        %% Compute: Step F5A
        if if_compute_summary == 1
            fprintf('\nNow runing %s.  ------> \n', stepF5A_memoryMetaMismatch_twoDecoder_Name_v);
            script_stepF5A_memoryMetaMismatch_twoDecoder();
        end
        
        part3_summary = struct;
    end
    
    
    %% Part 3C: Sequence-level meta-memory evidence
    if if_part3_summary == 1
        AUROC_meta_delay1_summary = AUROC_meta_delay1;
        
        x = meta_seqLevel_choice_delay1;
        y = offloadingProb_inOne(1:sum(numSeq(1:3)))';
        [r_meta_seqLevel_choice_delay1_summary,p_meta_seqLevel_choice_delay1_summary] = corr(x(~isnan(x)),y(~isnan(x)));
        
        x = meta_seqLevel_noChoice_delay1;
        y = offloadingProb_inOne(1:sum(numSeq(1:3)))';
        [r_meta_seqLevel_noChoice_delay1_summary,p_meta_seqLevel_noChoice_delay1_summary] = corr(x(~isnan(x)),y(~isnan(x)));
        
        x = meta_seqLevel_choiceMemory_delay1;
        y = offloadingProb_inOne(1:sum(numSeq(1:3)))';
        [r_meta_seqLevel_choiceMemory_delay1_summary,~] = corr(x(~isnan(x)),y(~isnan(x)));
        
        x = meta_seqLevel_choiceOffload_delay1;
        y = offloadingProb_inOne(1:sum(numSeq(1:3)))';
        [r_meta_seqLevel_choiceOffload_delay1_summary,~] = corr(x(~isnan(x)),y(~isnan(x)));
        
        
        part3_summary.AUROC_meta_delay1_summary = AUROC_meta_delay1_summary;
        part3_summary.r_meta_seqLevel_choice_delay1_summary = r_meta_seqLevel_choice_delay1_summary;
        part3_summary.p_meta_seqLevel_choice_delay1_summary = p_meta_seqLevel_choice_delay1_summary;
        part3_summary.r_meta_seqLevel_noChoice_delay1_summary = r_meta_seqLevel_noChoice_delay1_summary;
        part3_summary.p_meta_seqLevel_noChoice_delay1_summary = p_meta_seqLevel_noChoice_delay1_summary;
        
        part3_summary.r_meta_seqLevel_choiceMemory_delay1_summary = r_meta_seqLevel_choiceMemory_delay1_summary;
        part3_summary.r_meta_seqLevel_choiceOffload_delay1_summary = r_meta_seqLevel_choiceOffload_delay1_summary;
        
    end
    
    
    
    %% Part 3D, Right: Trial-level meta-memory evidence
    if if_part3_summary == 1
        b1 = seqRProb_metaLow_choice;
        b2 = seqRProb_metaHigh_choice;
        
        seqRProb_metaLow_choice_mean_summary = mean(b1,'omitnan');
        seqRProb_metaHigh_choice_mean_summary = mean(b2,'omitnan');
        
        [~,p_metaLowHigh_summary] = ttest(b1,b2);
        
        part3_summary.seqRProb_metaLow_choice_mean_summary = seqRProb_metaLow_choice_mean_summary;
        part3_summary.seqRProb_metaHigh_choice_mean_summary = seqRProb_metaHigh_choice_mean_summary;
        part3_summary.p_metaLowHigh_summary = p_metaLowHigh_summary;
        
    end
    
    %% (New) Part 3D, Right: Trial-level meta-memory evidence
    if if_part3_summary == 1        
        d1 = meta_trialLevel_allSeq_choiceOffloadMean;
        d2 = meta_trialLevel_allSeq_choiceMemoryStimuliMean;
                
        meta_allSeq_choiceOffloadMean_mean_summary = mean(d1,'omitnan');
        meta_allSeq_choiceMemoryStimuliMean_mean_summary = mean(d2,'omitnan');
        
        [~,p_choiceMeta_stimuli_summary] = ttest(d1,d2);
                
        part3_summary.meta_allSeq_choiceOffloadMean_mean_summary = meta_allSeq_choiceOffloadMean_mean_summary;
        part3_summary.meta_allSeq_choiceMemoryStimuliMean_mean_summary = meta_allSeq_choiceMemoryStimuliMean_mean_summary;
        part3_summary.p_choiceMeta_stimuli_summary = p_choiceMeta_stimuli_summary;
    end      
    

    
    %% Part 3F: Single neuron tuning property of memory precision and meta-memory
    if if_part3_summary == 1
        r2_seqPrecision_summary = r2_seqPrecision;
        r2_rProb_summary = r2_rProb;
        
        part3_summary.r2_seqPrecision_summary = r2_seqPrecision_summary;
        part3_summary.r2_rProb_summary = r2_rProb_summary;
        
    end
    
    %% Part 3G: Correlation between memory precision and meta-memory
    if if_part3_summary == 1
        %x = seqPrecision_neuron;
        x = seqPrecision_neuron_choice;
        y = meta_seqLevel_choice_delay1;
        [r_precisionMeta_seqLevel_summary,p_precisionMeta_seqLevel_summary] = corr(x(~isnan(x)),y(~isnan(x)));
        
        part3_summary.r_precisionMeta_seqLevel_summary = r_precisionMeta_seqLevel_summary;
        part3_summary.p_precisionMeta_seqLevel_summary = p_precisionMeta_seqLevel_summary;
        
    end
    
    
    %% Part 3H, Right: Trial-level memory precision VS. offloading rate evidence
    if if_part3_summary == 1
        c3 = seqRProb_memoryPrecisionLow_choiceResponse;
        c4 = seqRProb_memoryPrecisionHigh_choiceResponse;
        
        seqRProb_memoryPrecisionLow_choiceResponse_mean_summary = mean(c3,'omitnan');
        seqRProb_memoryPrecisionHigh_choiceResponse_mean_summary = mean(c4,'omitnan');
        
        [~,p_seqRProb_precisionLowHigh_response_summary] = ttest(c3,c4);
        
        part3_summary.seqRProb_memoryPrecisionLow_choiceResponse_mean_summary = seqRProb_memoryPrecisionLow_choiceResponse_mean_summary;
        part3_summary.seqRProb_memoryPrecisionHigh_choiceResponse_mean_summary = seqRProb_memoryPrecisionHigh_choiceResponse_mean_summary;
        part3_summary.p_seqRProb_precisionLowHigh_response_summary = p_seqRProb_precisionLowHigh_response_summary;
        
    end
    
    %% (New) Part 3H, Right: Trial-level memory precision VS. offloading rate evidence
    if if_part3_summary == 1        
        d4 = memoryPrecision_trialLevel_allSeq_choiceOffloadMean;
        d5 = memoryPrecision_trialLevel_allSeq_choiceMemoryStimuliMean;
        d6 = memoryPrecision_trialLevel_allSeq_choiceMemoryResponseMean;
                        
        memoryPrecision_allSeq_choiceOffloadMean_mean_summary = mean(d4,'omitnan');
        memoryPrecision_allSeq_choiceMemoryStimuliMean_mean_summary = mean(d5,'omitnan');
        memoryPrecision_allSeq_choiceMemoryResponseMean_mean_summary = mean(d6,'omitnan');
        
        [~,p_choicePrecision_stimuli_summary] = ttest(d4,d5);
        [~,p_choicePrecision_response_summary] = ttest(d4,d6);
                
        part3_summary.memoryPrecision_allSeq_choiceOffloadMean_mean_summary = memoryPrecision_allSeq_choiceOffloadMean_mean_summary;
        part3_summary.memoryPrecision_allSeq_choiceMemoryResponseMean_mean_summary = memoryPrecision_allSeq_choiceMemoryResponseMean_mean_summary;
        part3_summary.p_choicePrecision_response_summary = p_choicePrecision_response_summary;
        
        part3_summary.memoryPrecision_allSeq_choiceMemoryStimuliMean_mean_summary = memoryPrecision_allSeq_choiceMemoryStimuliMean_mean_summary;
        part3_summary.p_choicePrecision_stimuli_summary = p_choicePrecision_stimuli_summary;
        
    end        
    
    %% (New) Part 3: Trial-level memory precision VS. meta-memory evidence
    if if_part3_summary == 1
        e1 = seqMeta_memoryPrecisionLow_choiceStimuli;
        e2 = seqMeta_memoryPrecisionHigh_choiceStimuli;
        e3 = seqMeta_memoryPrecisionLow_choiceResponse;
        e4 = seqMeta_memoryPrecisionHigh_choiceResponse;
        
        seqMeta_memoryPrecisionLow_choiceStimuli_mean_summary = mean(e1,'omitnan');
        seqMeta_memoryPrecisionHigh_choiceStimuli_mean_summary = mean(e2,'omitnan');
        seqMeta_memoryPrecisionLow_choiceResponse_mean_summary = mean(e3,'omitnan');
        seqMeta_memoryPrecisionHigh_choiceResponse_mean_summary = mean(e4,'omitnan');
        
        [~,p_seqMeta_precisionLowHigh_stimuli_summary] = ttest(e1,e2);
        [~,p_seqMeta_precisionLowHigh_response_summary] = ttest(e3,e4);
        
        part3_summary.seqMeta_memoryPrecisionLow_choiceStimuli_mean_summary = seqMeta_memoryPrecisionLow_choiceStimuli_mean_summary;
        part3_summary.seqMeta_memoryPrecisionHigh_choiceStimuli_mean_summary = seqMeta_memoryPrecisionHigh_choiceStimuli_mean_summary;
        part3_summary.p_seqMeta_precisionLowHigh_stimuli_summary = p_seqMeta_precisionLowHigh_stimuli_summary;
        part3_summary.seqMeta_memoryPrecisionLow_choiceResponse_mean_summary = seqMeta_memoryPrecisionLow_choiceResponse_mean_summary;
        part3_summary.seqMeta_memoryPrecisionHigh_choiceResponse_mean_summary = seqMeta_memoryPrecisionHigh_choiceResponse_mean_summary;
        part3_summary.p_seqMeta_precisionLowHigh_response_summary = p_seqMeta_precisionLowHigh_response_summary;
        
        
        seqCorr_precisionMeta_choiceStimuli_mean_summary = mean(seqCorr_precisionMeta_choiceStimuli,'omitnan');
        seqCorr_precisionMeta_choiceResponse_mean_summary = mean(seqCorr_precisionMeta_choiceResponse,'omitnan');
  
        part3_summary.seqCorr_precisionMeta_choiceStimuli_mean_summary = seqCorr_precisionMeta_choiceStimuli_mean_summary;
        part3_summary.seqCorr_precisionMeta_choiceResponse_mean_summary = seqCorr_precisionMeta_choiceResponse_mean_summary;
                        
    end 
    
    
    %% (New) Part 3: Trial-level correlation between memory precision and meta-memory
    if if_part3_summary == 1
        % choice
        temptempBoolIndex_A = choiceBoolIndex';
        temptempBoolIndex_B = (~isnan(memoryPrecision_trialLevel)) & (~isnan(meta_trialLevel));
        
        temptempBoolIndex = temptempBoolIndex_A & temptempBoolIndex_B;
        
        [trialCorr_precisionMeta_choice_summary,~] = corr(memoryPrecision_trialLevel(temptempBoolIndex),meta_trialLevel(temptempBoolIndex));

        % choice-memory
        temptempBoolIndex_A = choiceMemoryBoolIndex';
        temptempBoolIndex_B = (~isnan(memoryPrecision_trialLevel)) & (~isnan(meta_trialLevel));
        
        temptempBoolIndex = temptempBoolIndex_A & temptempBoolIndex_B;
        
        [trialCorr_precisionMeta_choiceMemory_summary,~] = corr(memoryPrecision_trialLevel(temptempBoolIndex),meta_trialLevel(temptempBoolIndex));
        
        
        % choice-offload
        temptempBoolIndex_A = choiceOffloadBoolIndex';
        temptempBoolIndex_B = (~isnan(memoryPrecision_trialLevel)) & (~isnan(meta_trialLevel));
        
        temptempBoolIndex = temptempBoolIndex_A & temptempBoolIndex_B;
        
        [trialCorr_precisionMeta_choiceOffload_summary,~] = corr(memoryPrecision_trialLevel(temptempBoolIndex),meta_trialLevel(temptempBoolIndex));
        
        %
        part3_summary.trialCorr_precisionMeta_choice_summary = trialCorr_precisionMeta_choice_summary;
        part3_summary.trialCorr_precisionMeta_choiceMemory_summary = trialCorr_precisionMeta_choiceMemory_summary;
        part3_summary.trialCorr_precisionMeta_choiceOffload_summary = trialCorr_precisionMeta_choiceOffload_summary;
        
    end
    
    
    % Summary: Part 3
    if if_part3_summary == 1
        part3_summary_mulFOV.AUROC_meta_delay1_summary = [part3_summary_mulFOV.AUROC_meta_delay1_summary AUROC_meta_delay1_summary];
        part3_summary_mulFOV.r_meta_seqLevel_choice_delay1_summary = [part3_summary_mulFOV.r_meta_seqLevel_choice_delay1_summary r_meta_seqLevel_choice_delay1_summary];
        part3_summary_mulFOV.p_meta_seqLevel_choice_delay1_summary = [part3_summary_mulFOV.p_meta_seqLevel_choice_delay1_summary p_meta_seqLevel_choice_delay1_summary];
        part3_summary_mulFOV.r_meta_seqLevel_noChoice_delay1_summary = [part3_summary_mulFOV.r_meta_seqLevel_noChoice_delay1_summary r_meta_seqLevel_noChoice_delay1_summary];
        part3_summary_mulFOV.p_meta_seqLevel_noChoice_delay1_summary = [part3_summary_mulFOV.p_meta_seqLevel_noChoice_delay1_summary p_meta_seqLevel_noChoice_delay1_summary];
        
        part3_summary_mulFOV.r_meta_seqLevel_choiceMemory_delay1_summary = [part3_summary_mulFOV.r_meta_seqLevel_choiceMemory_delay1_summary r_meta_seqLevel_choiceMemory_delay1_summary];
        part3_summary_mulFOV.r_meta_seqLevel_choiceOffload_delay1_summary = [part3_summary_mulFOV.r_meta_seqLevel_choiceOffload_delay1_summary r_meta_seqLevel_choiceOffload_delay1_summary];
        
        part3_summary_mulFOV.seqRProb_metaLow_choice_mean_summary = [part3_summary_mulFOV.seqRProb_metaLow_choice_mean_summary seqRProb_metaLow_choice_mean_summary];
        part3_summary_mulFOV.seqRProb_metaHigh_choice_mean_summary = [part3_summary_mulFOV.seqRProb_metaHigh_choice_mean_summary seqRProb_metaHigh_choice_mean_summary];
        part3_summary_mulFOV.p_metaLowHigh_summary = [part3_summary_mulFOV.p_metaLowHigh_summary p_metaLowHigh_summary];
        
        part3_summary_mulFOV.meta_allSeq_choiceOffloadMean_mean_summary = [part3_summary_mulFOV.meta_allSeq_choiceOffloadMean_mean_summary meta_allSeq_choiceOffloadMean_mean_summary];
        part3_summary_mulFOV.meta_allSeq_choiceMemoryStimuliMean_mean_summary = [part3_summary_mulFOV.meta_allSeq_choiceMemoryStimuliMean_mean_summary meta_allSeq_choiceMemoryStimuliMean_mean_summary];
        part3_summary_mulFOV.p_choiceMeta_stimuli_summary = [part3_summary_mulFOV.p_choiceMeta_stimuli_summary p_choiceMeta_stimuli_summary];        
        
        %part3_summary_mulFOV.r2_seqPrecision_summary = [part3_summary_mulFOV.r2_seqPrecision_summary r2_seqPrecision_summary];
        %part3_summary_mulFOV.r2_rProb_summary = [part3_summary_mulFOV.r2_rProb_summary r2_rProb_summary];
        part3_summary_mulFOV.r2_seqPrecision_summary = [part3_summary_mulFOV.r2_seqPrecision_summary;r2_seqPrecision_summary];
        part3_summary_mulFOV.r2_rProb_summary = [part3_summary_mulFOV.r2_rProb_summary;r2_rProb_summary];        
        
        part3_summary_mulFOV.r2_seqPrecision_cell_summary = [part3_summary_mulFOV.r2_seqPrecision_cell_summary;{r2_seqPrecision_summary}];
        part3_summary_mulFOV.r2_rProb_cell_summary = [part3_summary_mulFOV.r2_rProb_cell_summary;{r2_rProb_summary}];        
        
        part3_summary_mulFOV.r_precisionMeta_seqLevel_summary = [part3_summary_mulFOV.r_precisionMeta_seqLevel_summary r_precisionMeta_seqLevel_summary];
        part3_summary_mulFOV.p_precisionMeta_seqLevel_summary = [part3_summary_mulFOV.p_precisionMeta_seqLevel_summary p_precisionMeta_seqLevel_summary];
        part3_summary_mulFOV.seqRProb_memoryPrecisionLow_choiceResponse_mean_summary = [part3_summary_mulFOV.seqRProb_memoryPrecisionLow_choiceResponse_mean_summary seqRProb_memoryPrecisionLow_choiceResponse_mean_summary];
        part3_summary_mulFOV.seqRProb_memoryPrecisionHigh_choiceResponse_mean_summary = [part3_summary_mulFOV.seqRProb_memoryPrecisionHigh_choiceResponse_mean_summary seqRProb_memoryPrecisionHigh_choiceResponse_mean_summary];
        part3_summary_mulFOV.p_seqRProb_precisionLowHigh_response_summary = [part3_summary_mulFOV.p_seqRProb_precisionLowHigh_response_summary p_seqRProb_precisionLowHigh_response_summary];        
        
        part3_summary_mulFOV.memoryPrecision_allSeq_choiceOffloadMean_mean_summary = [part3_summary_mulFOV.memoryPrecision_allSeq_choiceOffloadMean_mean_summary memoryPrecision_allSeq_choiceOffloadMean_mean_summary];
        part3_summary_mulFOV.memoryPrecision_allSeq_choiceMemoryResponseMean_mean_summary = [part3_summary_mulFOV.memoryPrecision_allSeq_choiceMemoryResponseMean_mean_summary memoryPrecision_allSeq_choiceMemoryResponseMean_mean_summary];
        part3_summary_mulFOV.p_choicePrecision_response_summary = [part3_summary_mulFOV.p_choicePrecision_response_summary p_choicePrecision_response_summary];

        part3_summary_mulFOV.memoryPrecision_allSeq_choiceMemoryStimuliMean_mean_summary = [part3_summary_mulFOV.memoryPrecision_allSeq_choiceMemoryStimuliMean_mean_summary memoryPrecision_allSeq_choiceMemoryStimuliMean_mean_summary];
        part3_summary_mulFOV.p_choicePrecision_stimuli_summary = [part3_summary_mulFOV.p_choicePrecision_stimuli_summary p_choicePrecision_stimuli_summary];
        
        part3_summary_mulFOV.seqMeta_memoryPrecisionLow_choiceStimuli_mean_summary = [part3_summary_mulFOV.seqMeta_memoryPrecisionLow_choiceStimuli_mean_summary seqMeta_memoryPrecisionLow_choiceStimuli_mean_summary];
        part3_summary_mulFOV.seqMeta_memoryPrecisionHigh_choiceStimuli_mean_summary = [part3_summary_mulFOV.seqMeta_memoryPrecisionHigh_choiceStimuli_mean_summary seqMeta_memoryPrecisionHigh_choiceStimuli_mean_summary];
        part3_summary_mulFOV.p_seqMeta_precisionLowHigh_stimuli_summary = [part3_summary_mulFOV.p_seqMeta_precisionLowHigh_stimuli_summary p_seqMeta_precisionLowHigh_stimuli_summary];
        part3_summary_mulFOV.seqMeta_memoryPrecisionLow_choiceResponse_mean_summary = [part3_summary_mulFOV.seqMeta_memoryPrecisionLow_choiceResponse_mean_summary seqMeta_memoryPrecisionLow_choiceResponse_mean_summary];
        part3_summary_mulFOV.seqMeta_memoryPrecisionHigh_choiceResponse_mean_summary = [part3_summary_mulFOV.seqMeta_memoryPrecisionHigh_choiceResponse_mean_summary seqMeta_memoryPrecisionHigh_choiceResponse_mean_summary];
        part3_summary_mulFOV.p_seqMeta_precisionLowHigh_response_summary = [part3_summary_mulFOV.p_seqMeta_precisionLowHigh_response_summary p_seqMeta_precisionLowHigh_response_summary];
        
        part3_summary_mulFOV.seqCorr_precisionMeta_choiceStimuli_mean_summary = [part3_summary_mulFOV.seqCorr_precisionMeta_choiceStimuli_mean_summary seqCorr_precisionMeta_choiceStimuli_mean_summary];
        part3_summary_mulFOV.seqCorr_precisionMeta_choiceResponse_mean_summary = [part3_summary_mulFOV.seqCorr_precisionMeta_choiceResponse_mean_summary seqCorr_precisionMeta_choiceResponse_mean_summary];
        
        part3_summary_mulFOV.trialCorr_precisionMeta_choice_summary = [part3_summary_mulFOV.trialCorr_precisionMeta_choice_summary trialCorr_precisionMeta_choice_summary];
        part3_summary_mulFOV.trialCorr_precisionMeta_choiceMemory_summary = [part3_summary_mulFOV.trialCorr_precisionMeta_choiceMemory_summary trialCorr_precisionMeta_choiceMemory_summary];
        part3_summary_mulFOV.trialCorr_precisionMeta_choiceOffload_summary = [part3_summary_mulFOV.trialCorr_precisionMeta_choiceOffload_summary trialCorr_precisionMeta_choiceOffload_summary];
    end
    
    %% Compute: Part 4
    if if_part4_summary == 1
        %% Compute: Step F5A
        % no need
        
        %% Compute: Step F5D
        if if_compute_summary == 1
            fprintf('\nNow runing %s.  ------> \n', stepF5D_test_mismatchMechnism_meta_crossTime_Name_v);
            script_stepF5D_test_mismatchMechnism_meta_crossTime();
        end
        
        %% Compute: Step F5E
        if if_compute_summary == 1
            fprintf('\nNow runing %s.  ------> \n', stepF5E_test_mismatchMechnism_memoryPrecision_crossTime_Name_v);
            script_stepF5E_test_mismatchMechnism_memoryPrecision_crossTime();
        end
        
        part4_summary = struct;
    end
    
    
    %% Part 4B, Left: Trial number of OverMismatch, HighMatch, LowMatch, UnderMismatch
    if if_part4_summary == 1
        trialBoolIndex_memoryPrecisionLow_metaLow_choice = trialBoolIndex_memoryPrecisionLow_choice & trialBoolIndex_metaLow_choice;
        trialBoolIndex_memoryPrecisionHigh_metaHigh_choice = trialBoolIndex_memoryPrecisionHigh_choice & trialBoolIndex_metaHigh_choice;
        trialBoolIndex_memoryPrecisionLow_metaHigh_choice = trialBoolIndex_memoryPrecisionLow_choice & trialBoolIndex_metaHigh_choice;
        trialBoolIndex_memoryPrecisionHigh_metaLow_choice = trialBoolIndex_memoryPrecisionHigh_choice & trialBoolIndex_metaLow_choice;
        
        trialNum_memoryPrecisionLow_metaLow_choice_summary = sum(trialBoolIndex_memoryPrecisionLow_metaLow_choice);
        trialNum_memoryPrecisionHigh_metaHigh_choice_summary = sum(trialBoolIndex_memoryPrecisionHigh_metaHigh_choice);
        trialNum_memoryPrecisionLow_metaHigh_choice_summary = sum(trialBoolIndex_memoryPrecisionLow_metaHigh_choice);
        trialNum_memoryPrecisionHigh_metaLow_choice_summary = sum(trialBoolIndex_memoryPrecisionHigh_metaLow_choice);
        
        part4_summary.trialNum_memoryPrecisionLow_metaLow_choice_summary = trialNum_memoryPrecisionLow_metaLow_choice_summary;
        part4_summary.trialNum_memoryPrecisionHigh_metaHigh_choice_summary = trialNum_memoryPrecisionHigh_metaHigh_choice_summary;
        part4_summary.trialNum_memoryPrecisionLow_metaHigh_choice_summary = trialNum_memoryPrecisionLow_metaHigh_choice_summary;
        part4_summary.trialNum_memoryPrecisionHigh_metaLow_choice_summary = trialNum_memoryPrecisionHigh_metaLow_choice_summary;
        
    end
    
    %% New Part 4
    if if_part4_summary == 1
        signifTimeStampCount_meta_mean_summary = temp_signifTimeStampCount_meta_mean;
        signifTimeStampCount_memory_mean_summary = temp_signifTimeStampCount_memory_mean;
        
        overTrialProp_resampled_mean_summary = mean(temp_overTrialProp_resampled);
        overTrialProp_resampled_median_summary = median(temp_overTrialProp_resampled);
        overTrialProp_summary = temp_overTrialProp; 
        
        
        temptempBoolIndex_A = (~isnan(memoryPrecision_trialLevel)) & (~isnan(meta_trialLevel));
        temptempBoolIndex_B = choiceMemoryCorrectBoolIndex';
        temptempBoolIndex_C = choiceMemoryErrorBoolIndex';
        temptempBoolIndex_D = choiceOffloadBoolIndex';
        
        temptempBoolIndex_CMC = temptempBoolIndex_A & temptempBoolIndex_B;
        temptempBoolIndex_CF = temptempBoolIndex_A & temptempBoolIndex_D;
        temptempBoolIndex_CME = temptempBoolIndex_A & temptempBoolIndex_C;
        
        temptempBoolIndex_CMC_summary = temptempBoolIndex_CMC;
        temptempBoolIndex_CF_summary = temptempBoolIndex_CF;
        temptempBoolIndex_CME_summary = temptempBoolIndex_CME;
        
        proj1_summary = proj1;
        proj2_summary = proj2;
        proj12_summary = proj12;
        proj3_summary = proj3;
        
        
        part4_summary.signifTimeStampCount_meta_mean_summary = signifTimeStampCount_meta_mean_summary;
        part4_summary.signifTimeStampCount_memory_mean_summary = signifTimeStampCount_memory_mean_summary;
        part4_summary.overTrialProp_resampled_mean_summary = overTrialProp_resampled_mean_summary;
        part4_summary.overTrialProp_resampled_median_summary = overTrialProp_resampled_median_summary;
        part4_summary.overTrialProp_summary = overTrialProp_summary;
        
        part4_summary.temptempBoolIndex_CMC_summary = temptempBoolIndex_CMC_summary;
        part4_summary.temptempBoolIndex_CF_summary = temptempBoolIndex_CF_summary;
        part4_summary.temptempBoolIndex_CME_summary = temptempBoolIndex_CME_summary;
        
        part4_summary.proj1_summary = proj1_summary;
        part4_summary.proj2_summary = proj2_summary;
        part4_summary.proj12_summary = proj12_summary;
        part4_summary.proj3_summary = proj3_summary;
        
    end
    
    %% Summary: Part 4
    if if_part4_summary == 1
        part4_summary_mulFOV.trialNum_memoryPrecisionLow_metaLow_choice_summary = [part4_summary_mulFOV.trialNum_memoryPrecisionLow_metaLow_choice_summary trialNum_memoryPrecisionLow_metaLow_choice_summary];
        part4_summary_mulFOV.trialNum_memoryPrecisionHigh_metaHigh_choice_summary = [part4_summary_mulFOV.trialNum_memoryPrecisionHigh_metaHigh_choice_summary trialNum_memoryPrecisionHigh_metaHigh_choice_summary];
        part4_summary_mulFOV.trialNum_memoryPrecisionLow_metaHigh_choice_summary = [part4_summary_mulFOV.trialNum_memoryPrecisionLow_metaHigh_choice_summary trialNum_memoryPrecisionLow_metaHigh_choice_summary];
        part4_summary_mulFOV.trialNum_memoryPrecisionHigh_metaLow_choice_summary = [part4_summary_mulFOV.trialNum_memoryPrecisionHigh_metaLow_choice_summary trialNum_memoryPrecisionHigh_metaLow_choice_summary];        
        
        
        part4_summary_mulFOV.signifTimeStampCount_meta_mean_summary = [part4_summary_mulFOV.signifTimeStampCount_meta_mean_summary signifTimeStampCount_meta_mean_summary];
        part4_summary_mulFOV.signifTimeStampCount_memory_mean_summary = [part4_summary_mulFOV.signifTimeStampCount_memory_mean_summary signifTimeStampCount_memory_mean_summary];
        part4_summary_mulFOV.overTrialProp_resampled_mean_summary = [part4_summary_mulFOV.overTrialProp_resampled_mean_summary overTrialProp_resampled_mean_summary];
        part4_summary_mulFOV.overTrialProp_resampled_median_summary = [part4_summary_mulFOV.overTrialProp_resampled_median_summary overTrialProp_resampled_median_summary];
        part4_summary_mulFOV.overTrialProp_summary = [part4_summary_mulFOV.overTrialProp_summary overTrialProp_summary];
                
        part4_summary_mulFOV.temptempBoolIndex_CMC_summary = [part4_summary_mulFOV.temptempBoolIndex_CMC_summary {temptempBoolIndex_CMC_summary}];
        part4_summary_mulFOV.temptempBoolIndex_CF_summary = [part4_summary_mulFOV.temptempBoolIndex_CF_summary {temptempBoolIndex_CF_summary}];
        part4_summary_mulFOV.temptempBoolIndex_CME_summary = [part4_summary_mulFOV.temptempBoolIndex_CME_summary {temptempBoolIndex_CME_summary}];
        
        part4_summary_mulFOV.proj1_summary = [part4_summary_mulFOV.proj1_summary {proj1_summary}];
        part4_summary_mulFOV.proj2_summary = [part4_summary_mulFOV.proj2_summary {proj2_summary}];
        part4_summary_mulFOV.proj12_summary = [part4_summary_mulFOV.proj12_summary {proj12_summary}];
        part4_summary_mulFOV.proj3_summary = [part4_summary_mulFOV.proj3_summary {proj3_summary}];
        
    end
    
    
    %% Compute: Part 5
    if if_part5_summary == 1
        %% Compute: Step DM1
        % no need
        
        %% Compute: baselineMeta_meta_GLM
        if if_compute_summary == 1
            fprintf('\nNow runing %s.  ------> \n', test_temp_baselineMeta_meta_GLM_Name_v);
            script_test_temp_baselineMeta_meta_GLM();
        end
        
        %% Compute: Step F5B
        if if_compute_summary == 1
            fprintf('\nNow runing %s.  ------> \n', stepF5B_test_trialHistory_Name_v);
            script_stepF5B_test_trialHistory();
        end
        
        part5_summary = struct;
    end
    
    a = 1;
    
    %% Part 5B: Baseline meta-memory evidence
    if if_part5_summary == 1
        AUROC_meta_baseline_summary = AUROC_meta_baseline;
        
        meta_trialLevel_baseline_choiceMemory = meta_trialLevel_baseline(choiceMemoryBoolIndex_validLength);
        meta_trialLevel_baseline_choiceOffload = meta_trialLevel_baseline(choiceOffloadBoolIndex_validLength);
        [~,p_baseline_choiceMemory_choiceOffload_summary] = ttest2(meta_trialLevel_baseline_choiceMemory,meta_trialLevel_baseline_choiceOffload);
        
        
        x = meta_seqLevel_choice_baseline;
        y = offloadingProb_inOne(1:sum(numSeq(1:3)))';
        [r_meta_seqLevel_choice_baseline_summary,p_meta_seqLevel_choice_baseline_summary] = corr(x(~isnan(x)),y(~isnan(x)));
        
        x = meta_seqLevel_noChoice_baseline;
        y = offloadingProb_inOne(1:sum(numSeq(1:3)))';
        [r_meta_seqLevel_noChoice_baseline_summary,p_meta_seqLevel_noChoice_baseline_summary] = corr(x(~isnan(x)),y(~isnan(x)));
        
        
        
        
        z1 = meta_trialLevel_baseline(temp_trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh,:);
        z2 = meta_trialLevel_baseline(temp_trialBoolIndex_memoryPrecisionLow_choiceOffloadLow,:);
        z3 = meta_trialLevel_baseline(temp_trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh,:);
        z4 = meta_trialLevel_baseline(temp_trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow,:);
                
        [~,p_z1_z4] = ttest2(z1,z4);         
        [~,p_z2_z3] = ttest2(z2,z3);        
        
        z5 = meta_trialLevel_baseline(trialBoolIndex_metaHigh_choice,:);
        z6 = meta_trialLevel_baseline(trialBoolIndex_metaLow_choice,:);
        
        [~,p_z5_z6] = ttest2(z5,z6);
        
        z1_mean = mean(z1,1,'omitnan');
        z2_mean = mean(z2,1,'omitnan');
        z3_mean = mean(z3,1,'omitnan');
        z4_mean = mean(z4,1,'omitnan');        
        z5_mean = mean(z5,1,'omitnan');
        z6_mean = mean(z6,1,'omitnan');
        
        baselineMeta_highMatch_summary = z1_mean;
        baselineMeta_lowMatch_summary = z2_mean;
        baselineMeta_overMismatch_summary = z3_mean;
        baselineMeta_underMismatch_summary = z4_mean;
        
        p_baselineMeta_highMatch_underMismatch_summary = p_z1_z4;
        p_baselineMeta_lowMatch_overMismatch_summary = p_z2_z3;
        
        baselineMeta_highMeta_summary = z5_mean;
        baselineMeta_lowMeta_summary = z6_mean;
        
        p_baselineMeta_highMeta_lowMeta_summary = p_z5_z6;
        
        
        part5_summary.AUROC_meta_baseline_summary = AUROC_meta_baseline_summary;
        part5_summary.p_baseline_choiceMemory_choiceOffload_summary = p_baseline_choiceMemory_choiceOffload_summary;
        part5_summary.r_meta_seqLevel_choice_baseline_summary = r_meta_seqLevel_choice_baseline_summary;
        part5_summary.p_meta_seqLevel_choice_baseline_summary = p_meta_seqLevel_choice_baseline_summary;
        part5_summary.r_meta_seqLevel_noChoice_baseline_summary = r_meta_seqLevel_noChoice_baseline_summary;
        part5_summary.p_meta_seqLevel_noChoice_baseline_summary = p_meta_seqLevel_noChoice_baseline_summary;
        
        part5_summary.baselineMeta_highMatch_summary = baselineMeta_highMatch_summary;
        part5_summary.baselineMeta_lowMatch_summary = baselineMeta_lowMatch_summary;
        part5_summary.baselineMeta_overMismatch_summary = baselineMeta_overMismatch_summary;
        part5_summary.baselineMeta_underMismatch_summary = baselineMeta_underMismatch_summary;        
        part5_summary.p_baselineMeta_highMatch_underMismatch_summary = p_baselineMeta_highMatch_underMismatch_summary;
        part5_summary.p_baselineMeta_lowMatch_overMismatch_summary = p_baselineMeta_lowMatch_overMismatch_summary;
        
        part5_summary.baselineMeta_highMeta_summary = baselineMeta_highMeta_summary;
        part5_summary.baselineMeta_lowMeta_summary = baselineMeta_lowMeta_summary;
        part5_summary.p_baselineMeta_highMeta_lowMeta_summary = p_baselineMeta_highMeta_lowMeta_summary;
        
    end
    
    %% Part 5C: Linear regression: baseline meta-memory, memory precision, meta-memory
    if if_part5_summary == 1
        r2_caseA_resampled;
        r2_caseB_resampled;
        r2_caseC_resampled;
        
        r2_metaRegress_caseA_resampled_mean_summary = mean(r2_caseA_resampled);
        r2_metaRegress_caseB_resampled_mean_summary = mean(r2_caseB_resampled);
        r2_metaRegress_caseC_resampled_mean_summary = mean(r2_caseC_resampled);
        
        beta_precision_metaRegress_caseC_resampled_mean_summary = mean(beta1_caseC_resampled);
        beta_baseline_metaRegress_caseC_resampled_mean_summary = mean(beta2_caseC_resampled);
        beta_interaction_metaRegress_caseC_resampled_mean_summary = mean(beta3_caseC_resampled);
        
        [~,pAB_metaRegress_summary,~,~] = ttest(r2_caseA_resampled,r2_caseB_resampled);
        [~,pAC_metaRegress_summary,~,~] = ttest(r2_caseA_resampled,r2_caseC_resampled);
        [~,pBC_metaRegress_summary,~,~] = ttest(r2_caseB_resampled,r2_caseC_resampled);
        
        meta_trialLevel_baseline;
        memoryPrecision_trialLevel;
        meta_trialLevel;
        
        part5_summary.meta_trialLevel_baseline = meta_trialLevel_baseline;
        part5_summary.memoryPrecision_trialLevel = memoryPrecision_trialLevel;
        part5_summary.meta_trialLevel = meta_trialLevel;
        part5_summary.r2_metaRegress_caseA_resampled_mean_summary = r2_metaRegress_caseA_resampled_mean_summary;
        part5_summary.r2_metaRegress_caseB_resampled_mean_summary = r2_metaRegress_caseB_resampled_mean_summary;
        part5_summary.r2_metaRegress_caseC_resampled_mean_summary = r2_metaRegress_caseC_resampled_mean_summary;
        part5_summary.pAB_metaRegress_summary = pAB_metaRegress_summary;
        part5_summary.pAC_metaRegress_summary = pAC_metaRegress_summary;
        part5_summary.pBC_metaRegress_summary = pBC_metaRegress_summary;

        part5_summary.beta_precision_metaRegress_caseC_resampled_mean_summary = beta_precision_metaRegress_caseC_resampled_mean_summary;
        part5_summary.beta_baseline_metaRegress_caseC_resampled_mean_summary = beta_baseline_metaRegress_caseC_resampled_mean_summary;
        part5_summary.beta_interaction_metaRegress_caseC_resampled_mean_summary = beta_interaction_metaRegress_caseC_resampled_mean_summary;
                
    end
    
    %% Part 5: Supp: Linear regression of baseline meta-memory, memory precision, meta-memory
    if if_part5_summary == 1

        r2_metaRegress_caseE2_CM_resampled_mean_summary = mean(r2_caseE2_resampled);
        beta_precision_metaRegress_caseE2_CM_resampled_mean_summary = mean(beta1_caseE2_resampled);
        beta_baseline_metaRegress_caseE2_CM_resampled_mean_summary = mean(beta2_caseE2_resampled);
        beta_interaction_metaRegress_caseE2_CM_resampled_mean_summary = mean(beta3_caseE2_resampled);
        
        r2_metaRegress_caseF_CF_resampled_mean_summary = mean(r2_caseF_resampled);
        beta_precision_metaRegress_caseF_CF_resampled_mean_summary = mean(beta1_caseF_resampled);
        beta_baseline_metaRegress_caseF_CF_resampled_mean_summary = mean(beta2_caseF_resampled);
        beta_interaction_metaRegress_caseF_CF_CM_resampled_mean_summary = mean(beta3_caseF_resampled);
        
        part5_summary.r2_metaRegress_caseE2_CM_resampled_mean_summary = r2_metaRegress_caseE2_CM_resampled_mean_summary;
        part5_summary.beta_precision_metaRegress_caseE2_CM_resampled_mean_summary = beta_precision_metaRegress_caseE2_CM_resampled_mean_summary;
        part5_summary.beta_baseline_metaRegress_caseE2_CM_resampled_mean_summary = beta_baseline_metaRegress_caseE2_CM_resampled_mean_summary;
        part5_summary.beta_interaction_metaRegress_caseE2_CM_resampled_mean_summary = beta_interaction_metaRegress_caseE2_CM_resampled_mean_summary;

        part5_summary.r2_metaRegress_caseF_CF_resampled_mean_summary = r2_metaRegress_caseF_CF_resampled_mean_summary;
        part5_summary.beta_precision_metaRegress_caseF_CF_resampled_mean_summary = beta_precision_metaRegress_caseF_CF_resampled_mean_summary;
        part5_summary.beta_baseline_metaRegress_caseF_CF_resampled_mean_summary = beta_baseline_metaRegress_caseF_CF_resampled_mean_summary;
        part5_summary.beta_interaction_metaRegress_caseF_CF_CM_resampled_mean_summary = beta_interaction_metaRegress_caseF_CF_CM_resampled_mean_summary;
        
        
        
        r2_caseD2B_resampled;
        r2_caseD2A_resampled;
        
        r2_mismatchRegress_caseD2B_baseline_resampled_mean_summary = mean(r2_caseD2B_resampled);
        r2_mismatchRegress_caseD2A_precision_resampled_mean_summary = mean(r2_caseD2A_resampled);
        
        part5_summary.r2_mismatchRegress_caseD2B_baseline_resampled_mean_summary = r2_mismatchRegress_caseD2B_baseline_resampled_mean_summary;
        part5_summary.r2_mismatchRegress_caseD2A_precision_resampled_mean_summary = r2_mismatchRegress_caseD2A_precision_resampled_mean_summary;
        
        
        % Offload trial VS. Low-strength mismatch trials
        meta_CME_baseline_summary = meta_trialLevel_baseline(trialBoolIndex_memoryPrecisionLow_metaHigh_choiceMemoryError,:);
        meta_CF_baseline_summary = meta_trialLevel_baseline(choiceOffloadBoolIndex,:);
                
        memoryPrecision_CME_baseline_summary = memoryPrecision_trialLevel(trialBoolIndex_memoryPrecisionLow_metaHigh_choiceMemoryError,:);
        memoryPrecision_CF_baseline_summary = memoryPrecision_trialLevel(choiceOffloadBoolIndex,:);
        
        memoryPrecision_CME_baseline_mean_summary = mean(memoryPrecision_CME_baseline_summary,'omitnan');
        memoryPrecision_CF_baseline_mean_summary = mean(memoryPrecision_CF_baseline_summary,'omitnan');
        meta_CME_baseline_mean_summary = mean(meta_CME_baseline_summary,'omitnan');
        meta_CF_baseline_mean_summary = mean(meta_CF_baseline_summary,'omitnan');
        
        memoryPrecision_CME_baseline_median_summary = median(memoryPrecision_CME_baseline_summary,'omitnan');
        memoryPrecision_CF_baseline_median_summary = median(memoryPrecision_CF_baseline_summary,'omitnan');
        meta_CME_baseline_median_summary = median(meta_CME_baseline_summary,'omitnan');
        meta_CF_baseline_median_summary = median(meta_CF_baseline_summary,'omitnan');
       
        part5_summary.memoryPrecision_CME_baseline_mean_summary = memoryPrecision_CME_baseline_mean_summary;
        part5_summary.memoryPrecision_CF_baseline_mean_summary = memoryPrecision_CF_baseline_mean_summary;
        part5_summary.meta_CME_baseline_mean_summary = meta_CME_baseline_mean_summary;
        part5_summary.meta_CF_baseline_mean_summary = meta_CF_baseline_mean_summary;
       
        part5_summary.memoryPrecision_CME_baseline_median_summary = memoryPrecision_CME_baseline_median_summary;
        part5_summary.memoryPrecision_CF_baseline_median_summary = memoryPrecision_CF_baseline_median_summary;
        part5_summary.meta_CME_baseline_median_summary = meta_CME_baseline_median_summary;
        part5_summary.meta_CF_baseline_median_summary = meta_CF_baseline_median_summary;
        
    end
    
    
    %% Part 5E: Trial history weight
    if if_part5_summary == 1
        AUROC_optimal_trialHistory_summary = AUROC_optimal_trialHistory;
        if exist('normStd_optimal','var') == 1
            normStd_optimal_trialHistory_summary = normStd_optimal;
        else
            normStd_optimal_trialHistory_summary = 0;
        end
        
        weight_trialHistory_summary = weight_trialHistory;
        
        part5_summary.AUROC_optimal_trialHistory_summary = AUROC_optimal_trialHistory_summary;
        part5_summary.normStd_optimal_trialHistory_summary = normStd_optimal_trialHistory_summary;
        part5_summary.weight_trialHistory_summary = weight_trialHistory_summary;
        
    end
    
    
    %% Part 5F, Left: Trial history of baseline meta
    if if_part5_summary == 1
        a5 = score_trialHistory(trialBoolIndex_metaHigh_choice_baseline,:);
        a6 = score_trialHistory(trialBoolIndex_metaLow_choice_baseline,:);
        
        [~,p_a5_a6] = ttest2(a5,a6);
        p_historyReward_baselineLowHigh_summary = p_a5_a6(end);
        
        a5_mean = mean(a5,1,'omitnan');
        a6_mean = mean(a6,1,'omitnan');
        
        historyRewardMean_baselineHigh_summary = a5_mean(end);
        historyRewardMean_baselineLow_summary = a6_mean(end);
        
        part5_summary.p_historyReward_baselineLowHigh_summary = p_historyReward_baselineLowHigh_summary;
        part5_summary.historyRewardMean_baselineHigh_summary = historyRewardMean_baselineHigh_summary;
        part5_summary.historyRewardMean_baselineLow_summary = historyRewardMean_baselineLow_summary;
        
    end
    
    %% Part 5F, Right: Trial history of meta
    if if_part5_summary == 1
        a9 = score_trialHistory(trialBoolIndex_metaHigh_choice,:);
        a10 = score_trialHistory(trialBoolIndex_metaLow_choice,:);
        
        [~,p_a9_a10] = ttest2(a9,a10);
        p_historyReward_metaLowHigh_summary = p_a9_a10(end);
        
        a9_mean = mean(a9,1,'omitnan');
        a10_mean = mean(a10,1,'omitnan');
        
        historyRewardMean_metaHigh_summary = a9_mean(end);
        historyRewardMean_metaLow_summary = a10_mean(end);
        
        part5_summary.p_historyReward_metaLowHigh_summary = p_historyReward_metaLowHigh_summary;
        part5_summary.historyRewardMean_metaHigh_summary = historyRewardMean_metaHigh_summary;
        part5_summary.historyRewardMean_metaLow_summary = historyRewardMean_metaLow_summary;
        
    end
    
    %% Part 5F, Right: Trial history of mismatch
    if if_part5_summary == 1
        a1 = score_trialHistory(temp_trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh,:);
        a2 = score_trialHistory(temp_trialBoolIndex_memoryPrecisionLow_choiceOffloadLow,:);
        
        a3 = score_trialHistory(temp_trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh,:);
        a4 = score_trialHistory(temp_trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow,:);
        
        % overMismatch-->highMatch-->lowMatch-->underMismatch
        x = [1*ones(size(a3,1),1);2*ones(size(a1,1),1);3*ones(size(a2,1),1);4*ones(size(a4,1),1)];
        y = [a3;a1;a2;a4];
        
        
        p_linear = nan(1,length(range_trialHistory));
        for tempi=1:length(range_trialHistory)
            temp_mdl = fitglm(x,y(:,tempi),'linear');
            p_linear(tempi) = temp_mdl.Coefficients.pValue(2);
        end
        
        [M,I]=min(p_linear);
                
        p_linear_historyReward_4Types_summary = p_linear(end);
        
        
        [~,p_a1_a4] = ttest2(a1,a4,'Tail','right');
        [~,p_a2_a3] = ttest2(a2,a3,'Tail','left');
        
        p_reward_highMatch_underMismatch_summary = p_a1_a4(end);
        p_reward_lowMatch_overMismatch_summary = p_a2_a3(end);
        
        
        temp_1_raw = a3(:,end);
        temp_2_raw = a1(:,end);
        temp_3_raw = a2(:,end);
        temp_4_raw = a4(:,end);
        
        temp_1 = temp_1_raw(~isnan(temp_1_raw));
        temp_2 = temp_2_raw(~isnan(temp_2_raw));
        temp_3 = temp_3_raw(~isnan(temp_3_raw));
        temp_4 = temp_4_raw(~isnan(temp_4_raw));
        
        historyRewardMean_OverMismatch_summary = mean(temp_1);        
        historyRewardMean_HighMatch_summary = mean(temp_2);
        historyRewardMean_LowMatch_summary = mean(temp_3);
        historyRewardMean_UnderMismatch_summary = mean(temp_4);
        
        
        part5_summary.p_reward_highMatch_underMismatch_summary = p_reward_highMatch_underMismatch_summary;
        part5_summary.p_reward_lowMatch_overMismatch_summary = p_reward_lowMatch_overMismatch_summary;        
        part5_summary.p_linear_historyReward_4Types_summary = p_linear_historyReward_4Types_summary;
        part5_summary.historyRewardMean_OverMismatch_summary = historyRewardMean_OverMismatch_summary;        
        part5_summary.historyRewardMean_HighMatch_summary = historyRewardMean_HighMatch_summary;
        part5_summary.historyRewardMean_LowMatch_summary = historyRewardMean_LowMatch_summary;
        part5_summary.historyRewardMean_UnderMismatch_summary = historyRewardMean_UnderMismatch_summary;
        
    end
    
    
    %% Summary: Part 5
    if if_part5_summary == 1
        part5_summary_mulFOV.AUROC_meta_baseline_summary = [part5_summary_mulFOV.AUROC_meta_baseline_summary AUROC_meta_baseline_summary];
        part5_summary_mulFOV.p_baseline_choiceMemory_choiceOffload_summary = [part5_summary_mulFOV.p_baseline_choiceMemory_choiceOffload_summary p_baseline_choiceMemory_choiceOffload_summary];
        part5_summary_mulFOV.r_meta_seqLevel_choice_baseline_summary = [part5_summary_mulFOV.r_meta_seqLevel_choice_baseline_summary r_meta_seqLevel_choice_baseline_summary];
        part5_summary_mulFOV.p_meta_seqLevel_choice_baseline_summary = [part5_summary_mulFOV.p_meta_seqLevel_choice_baseline_summary p_meta_seqLevel_choice_baseline_summary];
        part5_summary_mulFOV.r_meta_seqLevel_noChoice_baseline_summary = [part5_summary_mulFOV.r_meta_seqLevel_noChoice_baseline_summary r_meta_seqLevel_noChoice_baseline_summary];
        part5_summary_mulFOV.p_meta_seqLevel_noChoice_baseline_summary = [part5_summary_mulFOV.p_meta_seqLevel_noChoice_baseline_summary p_meta_seqLevel_noChoice_baseline_summary];
        
        
        part5_summary_mulFOV.baselineMeta_highMatch_summary = [part5_summary_mulFOV.baselineMeta_highMatch_summary baselineMeta_highMatch_summary];
        part5_summary_mulFOV.baselineMeta_lowMatch_summary = [part5_summary_mulFOV.baselineMeta_lowMatch_summary baselineMeta_lowMatch_summary];
        part5_summary_mulFOV.baselineMeta_overMismatch_summary = [part5_summary_mulFOV.baselineMeta_overMismatch_summary baselineMeta_overMismatch_summary];
        part5_summary_mulFOV.baselineMeta_underMismatch_summary = [part5_summary_mulFOV.baselineMeta_underMismatch_summary baselineMeta_underMismatch_summary];
        part5_summary_mulFOV.p_baselineMeta_highMatch_underMismatch_summary = [part5_summary_mulFOV.p_baselineMeta_highMatch_underMismatch_summary p_baselineMeta_highMatch_underMismatch_summary];
        part5_summary_mulFOV.p_baselineMeta_lowMatch_overMismatch_summary = [part5_summary_mulFOV.p_baselineMeta_lowMatch_overMismatch_summary p_baselineMeta_lowMatch_overMismatch_summary];
        
        part5_summary_mulFOV.baselineMeta_highMeta_summary = [part5_summary_mulFOV.baselineMeta_highMeta_summary baselineMeta_highMeta_summary];
        part5_summary_mulFOV.baselineMeta_lowMeta_summary = [part5_summary_mulFOV.baselineMeta_lowMeta_summary baselineMeta_lowMeta_summary];
        part5_summary_mulFOV.p_baselineMeta_highMeta_lowMeta_summary = [part5_summary_mulFOV.p_baselineMeta_highMeta_lowMeta_summary p_baselineMeta_highMeta_lowMeta_summary];
        
        
        part5_summary_mulFOV.meta_trialLevel_baseline_summary = [part5_summary_mulFOV.meta_trialLevel_baseline_summary; meta_trialLevel_baseline];
        part5_summary_mulFOV.memoryPrecision_trialLevel_summary = [part5_summary_mulFOV.memoryPrecision_trialLevel_summary; memoryPrecision_trialLevel];
        part5_summary_mulFOV.meta_trialLevel_summary = [part5_summary_mulFOV.meta_trialLevel_summary; meta_trialLevel];
        
        part5_summary_mulFOV.meta_trialLevel_baseline_cell_summary = [part5_summary_mulFOV.meta_trialLevel_baseline_cell_summary; {meta_trialLevel_baseline}];
        part5_summary_mulFOV.memoryPrecision_trialLevel_cell_summary = [part5_summary_mulFOV.memoryPrecision_trialLevel_cell_summary; {memoryPrecision_trialLevel}];
        part5_summary_mulFOV.meta_trialLevel_cell_summary = [part5_summary_mulFOV.meta_trialLevel_cell_summary; {meta_trialLevel}];
        
        part5_summary_mulFOV.r2_metaRegress_caseA_resampled_mean_summary = [part5_summary_mulFOV.r2_metaRegress_caseA_resampled_mean_summary r2_metaRegress_caseA_resampled_mean_summary];
        part5_summary_mulFOV.r2_metaRegress_caseB_resampled_mean_summary = [part5_summary_mulFOV.r2_metaRegress_caseB_resampled_mean_summary r2_metaRegress_caseB_resampled_mean_summary];
        part5_summary_mulFOV.r2_metaRegress_caseC_resampled_mean_summary = [part5_summary_mulFOV.r2_metaRegress_caseC_resampled_mean_summary r2_metaRegress_caseC_resampled_mean_summary];                
        part5_summary_mulFOV.pAB_metaRegress_summary = [part5_summary_mulFOV.pAB_metaRegress_summary pAB_metaRegress_summary];
        part5_summary_mulFOV.pAC_metaRegress_summary = [part5_summary_mulFOV.pAC_metaRegress_summary pAC_metaRegress_summary];
        part5_summary_mulFOV.pBC_metaRegress_summary = [part5_summary_mulFOV.pBC_metaRegress_summary pBC_metaRegress_summary];        
        
        part5_summary_mulFOV.beta_precision_metaRegress_caseC_resampled_mean_summary = [part5_summary_mulFOV.beta_precision_metaRegress_caseC_resampled_mean_summary beta_precision_metaRegress_caseC_resampled_mean_summary];
        part5_summary_mulFOV.beta_baseline_metaRegress_caseC_resampled_mean_summary = [part5_summary_mulFOV.beta_baseline_metaRegress_caseC_resampled_mean_summary beta_baseline_metaRegress_caseC_resampled_mean_summary];
        part5_summary_mulFOV.beta_interaction_metaRegress_caseC_resampled_mean_summary = [part5_summary_mulFOV.beta_interaction_metaRegress_caseC_resampled_mean_summary beta_interaction_metaRegress_caseC_resampled_mean_summary];
        
        
        part5_summary_mulFOV.r2_metaRegress_caseE2_CM_resampled_mean_summary = [part5_summary_mulFOV.r2_metaRegress_caseE2_CM_resampled_mean_summary r2_metaRegress_caseE2_CM_resampled_mean_summary];
        part5_summary_mulFOV.beta_precision_metaRegress_caseE2_CM_resampled_mean_summary = [part5_summary_mulFOV.beta_precision_metaRegress_caseE2_CM_resampled_mean_summary beta_precision_metaRegress_caseE2_CM_resampled_mean_summary];
        part5_summary_mulFOV.beta_baseline_metaRegress_caseE2_CM_resampled_mean_summary = [part5_summary_mulFOV.beta_baseline_metaRegress_caseE2_CM_resampled_mean_summary beta_baseline_metaRegress_caseE2_CM_resampled_mean_summary];
        part5_summary_mulFOV.beta_interaction_metaRegress_caseE2_CM_resampled_mean_summary = [part5_summary_mulFOV.beta_interaction_metaRegress_caseE2_CM_resampled_mean_summary beta_interaction_metaRegress_caseE2_CM_resampled_mean_summary];
        
        part5_summary_mulFOV.r2_metaRegress_caseF_CF_resampled_mean_summary = [part5_summary_mulFOV.r2_metaRegress_caseF_CF_resampled_mean_summary r2_metaRegress_caseF_CF_resampled_mean_summary];
        part5_summary_mulFOV.beta_precision_metaRegress_caseF_CF_resampled_mean_summary = [part5_summary_mulFOV.beta_precision_metaRegress_caseF_CF_resampled_mean_summary beta_precision_metaRegress_caseF_CF_resampled_mean_summary];
        part5_summary_mulFOV.beta_baseline_metaRegress_caseF_CF_resampled_mean_summary = [part5_summary_mulFOV.beta_baseline_metaRegress_caseF_CF_resampled_mean_summary beta_baseline_metaRegress_caseF_CF_resampled_mean_summary];
        part5_summary_mulFOV.beta_interaction_metaRegress_caseF_CF_CM_resampled_mean_summary = [part5_summary_mulFOV.beta_interaction_metaRegress_caseF_CF_CM_resampled_mean_summary beta_interaction_metaRegress_caseF_CF_CM_resampled_mean_summary];
        
        part5_summary_mulFOV.r2_mismatchRegress_caseD2B_baseline_resampled_mean_summary = [part5_summary_mulFOV.r2_mismatchRegress_caseD2B_baseline_resampled_mean_summary r2_mismatchRegress_caseD2B_baseline_resampled_mean_summary];
        part5_summary_mulFOV.r2_mismatchRegress_caseD2A_precision_resampled_mean_summary = [part5_summary_mulFOV.r2_mismatchRegress_caseD2A_precision_resampled_mean_summary r2_mismatchRegress_caseD2A_precision_resampled_mean_summary];
        
        
        part5_summary_mulFOV.memoryPrecision_CME_baseline_mean_summary = [part5_summary_mulFOV.memoryPrecision_CME_baseline_mean_summary memoryPrecision_CME_baseline_mean_summary];
        part5_summary_mulFOV.memoryPrecision_CF_baseline_mean_summary = [part5_summary_mulFOV.memoryPrecision_CF_baseline_mean_summary memoryPrecision_CF_baseline_mean_summary];
        part5_summary_mulFOV.meta_CME_baseline_mean_summary = [part5_summary_mulFOV.meta_CME_baseline_mean_summary meta_CME_baseline_mean_summary];
        part5_summary_mulFOV.meta_CF_baseline_mean_summary = [part5_summary_mulFOV.meta_CF_baseline_mean_summary meta_CF_baseline_mean_summary];
        
        part5_summary_mulFOV.memoryPrecision_CME_baseline_median_summary = [part5_summary_mulFOV.memoryPrecision_CME_baseline_median_summary memoryPrecision_CME_baseline_median_summary];
        part5_summary_mulFOV.memoryPrecision_CF_baseline_median_summary = [part5_summary_mulFOV.memoryPrecision_CF_baseline_median_summary memoryPrecision_CF_baseline_median_summary];
        part5_summary_mulFOV.meta_CME_baseline_median_summary = [part5_summary_mulFOV.meta_CME_baseline_median_summary meta_CME_baseline_median_summary];
        part5_summary_mulFOV.meta_CF_baseline_median_summary = [part5_summary_mulFOV.meta_CF_baseline_median_summary meta_CF_baseline_median_summary];
        
        
        part5_summary_mulFOV.AUROC_optimal_trialHistory_summary = [part5_summary_mulFOV.AUROC_optimal_trialHistory_summary AUROC_optimal_trialHistory_summary];
        part5_summary_mulFOV.normStd_optimal_trialHistory_summary = [part5_summary_mulFOV.normStd_optimal_trialHistory_summary normStd_optimal_trialHistory_summary];
        
        part5_summary_mulFOV.weight_trialHistory_summary = [part5_summary_mulFOV.weight_trialHistory_summary; weight_trialHistory_summary];
        
        part5_summary_mulFOV.p_historyReward_baselineLowHigh_summary = [part5_summary_mulFOV.p_historyReward_baselineLowHigh_summary p_historyReward_baselineLowHigh_summary];
        part5_summary_mulFOV.historyRewardMean_baselineHigh_summary = [part5_summary_mulFOV.historyRewardMean_baselineHigh_summary historyRewardMean_baselineHigh_summary];
        part5_summary_mulFOV.historyRewardMean_baselineLow_summary = [part5_summary_mulFOV.historyRewardMean_baselineLow_summary historyRewardMean_baselineLow_summary];
        
        part5_summary_mulFOV.p_historyReward_metaLowHigh_summary = [part5_summary_mulFOV.p_historyReward_metaLowHigh_summary p_historyReward_metaLowHigh_summary];
        part5_summary_mulFOV.historyRewardMean_metaHigh_summary = [part5_summary_mulFOV.historyRewardMean_metaHigh_summary historyRewardMean_metaHigh_summary];
        part5_summary_mulFOV.historyRewardMean_metaLow_summary = [part5_summary_mulFOV.historyRewardMean_metaLow_summary historyRewardMean_metaLow_summary];
        
        part5_summary_mulFOV.p_reward_highMatch_underMismatch_summary = [part5_summary_mulFOV.p_reward_highMatch_underMismatch_summary p_reward_highMatch_underMismatch_summary];
        part5_summary_mulFOV.p_reward_lowMatch_overMismatch_summary = [part5_summary_mulFOV.p_reward_lowMatch_overMismatch_summary p_reward_lowMatch_overMismatch_summary];
        
        part5_summary_mulFOV.p_linear_historyReward_4Types_summary = [part5_summary_mulFOV.p_linear_historyReward_4Types_summary p_linear_historyReward_4Types_summary];
        part5_summary_mulFOV.historyRewardMean_OverMismatch_summary = [part5_summary_mulFOV.historyRewardMean_OverMismatch_summary historyRewardMean_OverMismatch_summary];        
        part5_summary_mulFOV.historyRewardMean_HighMatch_summary = [part5_summary_mulFOV.historyRewardMean_HighMatch_summary historyRewardMean_HighMatch_summary];
        part5_summary_mulFOV.historyRewardMean_LowMatch_summary = [part5_summary_mulFOV.historyRewardMean_LowMatch_summary historyRewardMean_LowMatch_summary];
        part5_summary_mulFOV.historyRewardMean_UnderMismatch_summary = [part5_summary_mulFOV.historyRewardMean_UnderMismatch_summary historyRewardMean_UnderMismatch_summary];
    end
    
    
    
    %% Compute: Part 6
    if if_part6_summary == 1
        %% Compute: Step F4A
        % no need
        
        %% Compute: Step F4B
        if if_compute_summary == 1
            fprintf('\nNow runing %s.  ------> \n', stepF4B_test_plot_FOV_wholeImage_complexTuning_Name_v);
            script_stepF4B_test_plot_FOV_wholeImage_complexTuning();
        end
        
        %% Compute: baselineMeta_meta_GLM
        %no need
        
        part6_summary = struct;
    end
    
    
    %% Part 6A-6D: Neuron tuning structure
    if if_part6_summary == 1
        beta_memoryPrecision_summary = beta_memoryPrecision;
        beta_choiceMemory_summary = beta_choiceMemory;
        beta_choiceMemory_baseline_summary = beta_choiceMemory_baseline;
        
        beta_seqPrecision_summary = beta_seqPrecision;
        beta_gProb_summary = beta_gProb;
        beta_gProb_baseline_summary = beta_gProb_baseline;     
        
        r2_memoryPrecision_summary = r2_memoryPrecision;
        r2_choiceMemory_summary = r2_choiceMemory;
        r2_choiceMemory_baseline_summary = r2_choiceMemory_baseline;
        
        r2_seqPrecision_summary = r2_seqPrecision;
        r2_gProb_summary = r2_gProb;
        r2_gProb_baseline_summary = r2_gProb_baseline;         
                
        beta_6loc_summary = beta_6loc;
        r2_6loc_summary = r2_6loc;      
        
        
        betaAngle_6loc_halfAB_norm_unitView_summary = betaAngle_6loc_halfAB_norm_resampled_unitView;
        betaAngle_gProb_halfAB_norm_unitView_summary = betaAngle_gProb_halfAB_norm_resampled_unitView;
        betaAngle_choiceMemory_baseline_halfAB_norm_unitView_summary = betaAngle_choiceMemory_baseline_halfAB_norm_resampled_unitView;
        
        
        part6_summary.beta_memoryPrecision_summary = beta_memoryPrecision_summary;
        part6_summary.beta_choiceMemory_summary = beta_choiceMemory_summary;
        part6_summary.beta_choiceMemory_baseline_summary = beta_choiceMemory_baseline_summary;
        part6_summary.beta_seqPrecision_summary = beta_seqPrecision_summary;
        part6_summary.beta_gProb_summary = beta_gProb_summary;
        part6_summary.beta_gProb_baseline_summary = beta_gProb_baseline_summary;
        
        part6_summary.r2_memoryPrecision_summary = r2_memoryPrecision_summary;
        part6_summary.r2_choiceMemory_summary = r2_choiceMemory_summary;
        part6_summary.r2_choiceMemory_baseline_summary = r2_choiceMemory_baseline_summary;
        part6_summary.r2_seqPrecision_summary = r2_seqPrecision_summary;
        part6_summary.r2_gProb_summary = r2_gProb_summary;
        part6_summary.r2_gProb_baseline_summary = r2_gProb_baseline_summary;

        part6_summary.beta_6loc_summary = beta_6loc_summary;
        part6_summary.r2_6loc_summary = r2_6loc_summary;   
        
        
        part6_summary.betaAngle_6loc_halfAB_norm_unitView_summary = betaAngle_6loc_halfAB_norm_unitView_summary;
        part6_summary.betaAngle_gProb_halfAB_norm_unitView_summary = betaAngle_gProb_halfAB_norm_unitView_summary;
        part6_summary.betaAngle_choiceMemory_baseline_halfAB_norm_unitView_summary = betaAngle_choiceMemory_baseline_halfAB_norm_unitView_summary;
        
        
    end
    
    
    %% Part 6 new: Subspace centered view of neuron tuning
    if if_part6_summary == 1                
        beta_prior_linearAxis_summary = beta_prior_linearAxis;
        beta_memory_linearAxis_summary = beta_memory_linearAxis;
        beta_meta_linearAxis_summary = beta_meta_linearAxis;
        
        r2_prior_linearAxis_summary = r2_prior_linearAxis;
        r2_memory_linearAxis_summary = r2_memory_linearAxis;
        r2_meta_linearAxis_summary = r2_meta_linearAxis;
        
        temp_degree_aa_prior_summary = temp_degree_aa_prior;
        temp_degree_bb_memory_summary = temp_degree_bb_memory;
        temp_degree_cc_meta_summary = temp_degree_cc_meta;
        
        temp_degree_ab_priorMemory_summary = temp_degree_ab_priorMemory;
        temp_degree_ac_priorMeta_summary = temp_degree_ac_priorMeta;
        temp_degree_bc_memoryMeta_summary = temp_degree_bc_memoryMeta;
        
        
        temp_VAF_ratio_priorPrior_summary = temp_VAF_ratio_priorPrior;
        temp_VAF_ratio_memoryMemory_summary = temp_VAF_ratio_memoryMemory;
        temp_VAF_ratio_metaMeta_summary = temp_VAF_ratio_metaMeta;
        temp_VAF_ratio_priorMemory_summary = temp_VAF_ratio_priorMemory;
        temp_VAF_ratio_priorMeta_summary = temp_VAF_ratio_priorMeta;
        temp_VAF_ratio_memoryMeta_summary = temp_VAF_ratio_memoryMeta;
        
        
        temp_EV_neuronPrior_summary = temp_EV_neuronPrior;
        temp_EV_neuronMemory_summary = temp_EV_neuronMemory;
        temp_EV_neuronMeta_summary = temp_EV_neuronMeta;
        pca_explained_baseline_summary = explained_baseline;
        pca_explained_delay_summary = explained_delay;
        pca_r_eachPC_VS_eventVector_baseline_summary = temp_r_eachPC_VS_eventVector_baseline;
        pca_r_eachPC_VS_eventVector_delay_summary = temp_r_eachPC_VS_eventVector_delay;
        
    
        part6_summary.beta_prior_linearAxis_summary = beta_prior_linearAxis_summary;
        part6_summary.beta_memory_linearAxis_summary = beta_memory_linearAxis_summary;
        part6_summary.beta_meta_linearAxis_summary = beta_meta_linearAxis_summary;
        
        part6_summary.r2_prior_linearAxis_summary = r2_prior_linearAxis_summary;
        part6_summary.r2_memory_linearAxis_summary = r2_memory_linearAxis_summary;
        part6_summary.r2_meta_linearAxis_summary = r2_meta_linearAxis_summary;
        
        
        part6_summary.temp_degree_aa_prior_summary = temp_degree_aa_prior_summary;
        part6_summary.temp_degree_bb_memory_summary = temp_degree_bb_memory_summary;
        part6_summary.temp_degree_cc_meta_summary = temp_degree_cc_meta_summary;
        
        part6_summary.temp_degree_ab_priorMemory_summary = temp_degree_ab_priorMemory_summary;
        part6_summary.temp_degree_ac_priorMeta_summary = temp_degree_ac_priorMeta_summary;
        part6_summary.temp_degree_bc_memoryMeta_summary = temp_degree_bc_memoryMeta_summary;

        
        part6_summary.temp_VAF_ratio_priorPrior_summary = temp_VAF_ratio_priorPrior_summary;
        part6_summary.temp_VAF_ratio_memoryMemory_summary = temp_VAF_ratio_memoryMemory_summary;
        part6_summary.temp_VAF_ratio_metaMeta_summary = temp_VAF_ratio_metaMeta_summary;
        part6_summary.temp_VAF_ratio_priorMemory_summary = temp_VAF_ratio_priorMemory_summary;
        part6_summary.temp_VAF_ratio_priorMeta_summary = temp_VAF_ratio_priorMeta_summary;
        part6_summary.temp_VAF_ratio_memoryMeta_summary = temp_VAF_ratio_memoryMeta_summary;
        
        part6_summary.temp_EV_neuronPrior_summary = temp_EV_neuronPrior_summary;
        part6_summary.temp_EV_neuronMemory_summary = temp_EV_neuronMemory_summary;
        part6_summary.temp_EV_neuronMeta_summary = temp_EV_neuronMeta_summary;
        part6_summary.pca_explained_baseline_summary = pca_explained_baseline_summary;
        part6_summary.pca_explained_delay_summary = pca_explained_delay_summary;
        part6_summary.pca_r_eachPC_VS_eventVector_baseline_summary = pca_r_eachPC_VS_eventVector_baseline_summary;
        part6_summary.pca_r_eachPC_VS_eventVector_delay_summary = pca_r_eachPC_VS_eventVector_delay_summary;
        
    end
            
    %% Part 6E: Neuron spatial organization, clustering index
    if if_part6_summary == 1
        %% Summary
        clusterIndex_memoryPrecision_disBin_mean_summary = clusterIndex_memoryPrecision_disBin_mean;
        clusterIndex_choiceMemory_disBin_mean_summary = clusterIndex_choiceMemory_disBin_mean;
        clusterIndex_choiceMemory_baseline_disBin_mean_summary = clusterIndex_choiceMemory_baseline_disBin_mean;
                
        prctile_high_memoryPrecision_shuffled_summary = prctile(clusterIndex_memoryPrecision_disBin_shuffled_mean,threshold_prctile_clusterB_high,2);
        prctile_high_choiceMemory_shuffled_summary = prctile(clusterIndex_choiceMemory_disBin_shuffled_mean,threshold_prctile_clusterB_high,2);
        prctile_high_choiceMemory_baseline_shuffled_summary = prctile(clusterIndex_choiceMemory_baseline_disBin_shuffled_mean,threshold_prctile_clusterB_high,2);
        
        clusterIndex_memoryPrecision_shuffled_meanB_summary = clusterIndex_memoryPrecision_disBin_shuffled_meanB;
        clusterIndex_choiceMemory_shuffled_meanB_summary = clusterIndex_choiceMemory_disBin_shuffled_meanB;
        clusterIndex_choiceMemory_baseline_shuffled_meanB_summary =clusterIndex_choiceMemory_baseline_disBin_shuffled_meanB;        
        
        part6_summary.clusterIndex_memoryPrecision_disBin_mean_summary = clusterIndex_memoryPrecision_disBin_mean_summary;
        part6_summary.clusterIndex_choiceMemory_disBin_mean_summary = clusterIndex_choiceMemory_disBin_mean_summary;
        part6_summary.clusterIndex_choiceMemory_baseline_disBin_mean_summary = clusterIndex_choiceMemory_baseline_disBin_mean_summary;
        part6_summary.prctile_high_memoryPrecision_shuffled_summary = prctile_high_memoryPrecision_shuffled_summary;
        part6_summary.prctile_high_choiceMemory_shuffled_summary = prctile_high_choiceMemory_shuffled_summary;
        part6_summary.prctile_high_choiceMemory_baseline_shuffled_summary = prctile_high_choiceMemory_baseline_shuffled_summary;
        part6_summary.clusterIndex_memoryPrecision_shuffled_meanB_summary = clusterIndex_memoryPrecision_shuffled_meanB_summary;
        part6_summary.clusterIndex_choiceMemory_shuffled_meanB_summary = clusterIndex_choiceMemory_shuffled_meanB_summary;
        part6_summary.clusterIndex_choiceMemory_baseline_shuffled_meanB_summary = clusterIndex_choiceMemory_baseline_shuffled_meanB_summary;
        
    end
    
    
    %% Part 6F: Neuron spatial organization, centroid distance
    if if_part6_summary == 1
        %% Summary
        centriodDis_all_shuffled_prctileA_summary = centriodDis_all_shuffled_prctileA;
        centriodDis_all_shuffled_prctileB_summary = centriodDis_all_shuffled_prctileB;
        centriodDis_all_shuffled_mean_summary = mean(centriodDis_all_shuffled);
        centriodDis_all_summary = centriodDis_all;
        
        centriodDis_all_AC_shuffled_prctileA_summary = centriodDis_all_AC_shuffled_prctileA;
        centriodDis_all_AC_shuffled_prctileB_summary = centriodDis_all_AC_shuffled_prctileB;
        centriodDis_all_AC_shuffled_mean_summary = mean(centriodDis_all_AC_shuffled);        
        centriodDis_AC_all_summary = centriodDis_AC_all;
        
        centriodDis_all_BC_shuffled_prctileA_summary = centriodDis_all_BC_shuffled_prctileA;
        centriodDis_all_BC_shuffled_prctileB_summary = centriodDis_all_BC_shuffled_prctileB;
        centriodDis_all_BC_shuffled_mean_summary = mean(centriodDis_all_BC_shuffled);                
        centriodDis_BC_all_summary = centriodDis_BC_all;
        
        %% New spatial analysis
        spatialDisNew_1to1_mean_summary = temp_disNew_1to1_mean;
        spatialDisNew_2to2_mean_summary = temp_disNew_2to2_mean;
        spatialDisNew_3to3_mean_summary = temp_disNew_3to3_mean;
        
        spatialDisNew_1toOthers_mean_summary = temp_disNew_1toOthers_mean;
        spatialDisNew_2toOthers_mean_summary = temp_disNew_2toOthers_mean;
        spatialDisNew_3toOthers_mean_summary = temp_disNew_3toOthers_mean;
        
        spatialCentriodDis_1toOthers_summary = centriodDis_1toOthers;
        spatialCentriodDis_2toOthers_summary = centriodDis_2toOthers;
        spatialCentriodDis_3toOthers_summary = centriodDis_3toOthers;

        %%
        part6_summary.centriodDis_all_shuffled_prctileA_summary = centriodDis_all_shuffled_prctileA_summary;
        part6_summary.centriodDis_all_shuffled_prctileB_summary = centriodDis_all_shuffled_prctileB_summary;
        part6_summary.centriodDis_all_shuffled_mean_summary = centriodDis_all_shuffled_mean_summary;        
        part6_summary.centriodDis_all_summary = centriodDis_all_summary;
        part6_summary.centriodDis_all_AC_shuffled_prctileA_summary = centriodDis_all_AC_shuffled_prctileA_summary;
        part6_summary.centriodDis_all_AC_shuffled_prctileB_summary = centriodDis_all_AC_shuffled_prctileB_summary;
        part6_summary.centriodDis_all_AC_shuffled_mean_summary = centriodDis_all_AC_shuffled_mean_summary;                
        part6_summary.centriodDis_AC_all_summary = centriodDis_AC_all_summary;
        part6_summary.centriodDis_all_BC_shuffled_prctileA_summary = centriodDis_all_BC_shuffled_prctileA_summary;
        part6_summary.centriodDis_all_BC_shuffled_prctileB_summary = centriodDis_all_BC_shuffled_prctileB_summary;
        part6_summary.centriodDis_all_BC_shuffled_mean_summary = centriodDis_all_BC_shuffled_mean_summary;                        
        part6_summary.centriodDis_BC_all_summary = centriodDis_BC_all_summary;
        
        
        part6_summary.spatialDisNew_1to1_mean_summary = spatialDisNew_1to1_mean_summary;
        part6_summary.spatialDisNew_2to2_mean_summary = spatialDisNew_2to2_mean_summary;
        part6_summary.spatialDisNew_3to3_mean_summary = spatialDisNew_3to3_mean_summary;
        part6_summary.spatialDisNew_1toOthers_mean_summary = spatialDisNew_1toOthers_mean_summary;
        part6_summary.spatialDisNew_2toOthers_mean_summary = spatialDisNew_2toOthers_mean_summary;
        part6_summary.spatialDisNew_3toOthers_mean_summary = spatialDisNew_3toOthers_mean_summary;
        part6_summary.spatialCentriodDis_1toOthers_summary = spatialCentriodDis_1toOthers_summary;
        part6_summary.spatialCentriodDis_2toOthers_summary = spatialCentriodDis_2toOthers_summary;
        part6_summary.spatialCentriodDis_3toOthers_summary = spatialCentriodDis_3toOthers_summary;
        
    end
    
    %% Summary: Part 6
    if if_part6_summary == 1    
        %part6_summary_mulFOV.beta_memoryPrecision_summary = [part6_summary_mulFOV.beta_memoryPrecision_summary beta_memoryPrecision_summary];
        %part6_summary_mulFOV.beta_choiceMemory_summary = [part6_summary_mulFOV.beta_choiceMemory_summary beta_choiceMemory_summary];
        %part6_summary_mulFOV.beta_choiceMemory_baseline_summary = [part6_summary_mulFOV.beta_choiceMemory_baseline_summary beta_choiceMemory_baseline_summary];
        part6_summary_mulFOV.beta_memoryPrecision_summary = [part6_summary_mulFOV.beta_memoryPrecision_summary; beta_memoryPrecision_summary];
        part6_summary_mulFOV.beta_choiceMemory_summary = [part6_summary_mulFOV.beta_choiceMemory_summary; beta_choiceMemory_summary];
        part6_summary_mulFOV.beta_choiceMemory_baseline_summary = [part6_summary_mulFOV.beta_choiceMemory_baseline_summary; beta_choiceMemory_baseline_summary];        
        
        part6_summary_mulFOV.beta_memoryPrecision_cell_summary = [part6_summary_mulFOV.beta_memoryPrecision_cell_summary; {beta_memoryPrecision_summary}];
        part6_summary_mulFOV.beta_choiceMemory_cell_summary = [part6_summary_mulFOV.beta_choiceMemory_cell_summary; {beta_choiceMemory_summary}];
        part6_summary_mulFOV.beta_choiceMemory_baseline_cell_summary = [part6_summary_mulFOV.beta_choiceMemory_baseline_cell_summary; {beta_choiceMemory_baseline_summary}];        
        
        part6_summary_mulFOV.beta_seqPrecision_summary = [part6_summary_mulFOV.beta_seqPrecision_summary; beta_seqPrecision_summary];
        part6_summary_mulFOV.beta_gProb_summary = [part6_summary_mulFOV.beta_gProb_summary; beta_gProb_summary];
        part6_summary_mulFOV.beta_gProb_baseline_summary = [part6_summary_mulFOV.beta_gProb_baseline_summary; beta_gProb_baseline_summary];        
        
        part6_summary_mulFOV.beta_seqPrecision_cell_summary = [part6_summary_mulFOV.beta_seqPrecision_cell_summary; {beta_seqPrecision_summary}];
        part6_summary_mulFOV.beta_gProb_cell_summary = [part6_summary_mulFOV.beta_gProb_cell_summary; {beta_gProb_summary}];
        part6_summary_mulFOV.beta_gProb_baseline_cell_summary = [part6_summary_mulFOV.beta_gProb_baseline_cell_summary; {beta_gProb_baseline_summary}];        
         
        
        part6_summary_mulFOV.r2_memoryPrecision_summary = [part6_summary_mulFOV.r2_memoryPrecision_summary; r2_memoryPrecision_summary];
        part6_summary_mulFOV.r2_choiceMemory_summary = [part6_summary_mulFOV.r2_choiceMemory_summary; r2_choiceMemory_summary];
        part6_summary_mulFOV.r2_choiceMemory_baseline_summary = [part6_summary_mulFOV.r2_choiceMemory_baseline_summary; r2_choiceMemory_baseline_summary];
        part6_summary_mulFOV.r2_seqPrecision_summary = [part6_summary_mulFOV.r2_seqPrecision_summary; r2_seqPrecision_summary];
        part6_summary_mulFOV.r2_gProb_summary = [part6_summary_mulFOV.r2_gProb_summary; r2_gProb_summary];
        part6_summary_mulFOV.r2_gProb_baseline_summary = [part6_summary_mulFOV.r2_gProb_baseline_summary; r2_gProb_baseline_summary];
        part6_summary_mulFOV.beta_6loc_summary = [part6_summary_mulFOV.beta_6loc_summary; beta_6loc_summary];
        part6_summary_mulFOV.r2_6loc_summary = [part6_summary_mulFOV.r2_6loc_summary; r2_6loc_summary];
        
        part6_summary_mulFOV.r2_memoryPrecision_cell_summary = [part6_summary_mulFOV.r2_memoryPrecision_cell_summary; {r2_memoryPrecision_summary}];
        part6_summary_mulFOV.r2_choiceMemory_cell_summary = [part6_summary_mulFOV.r2_choiceMemory_cell_summary; {r2_choiceMemory_summary}];
        part6_summary_mulFOV.r2_choiceMemory_baseline_cell_summary = [part6_summary_mulFOV.r2_choiceMemory_baseline_cell_summary; {r2_choiceMemory_baseline_summary}];
        part6_summary_mulFOV.r2_seqPrecision_cell_summary = [part6_summary_mulFOV.r2_seqPrecision_cell_summary; {r2_seqPrecision_summary}];
        part6_summary_mulFOV.r2_gProb_cell_summary = [part6_summary_mulFOV.r2_gProb_cell_summary; {r2_gProb_summary}];
        part6_summary_mulFOV.r2_gProb_baseline_cell_summary = [part6_summary_mulFOV.r2_gProb_baseline_cell_summary; {r2_gProb_baseline_summary}];
        part6_summary_mulFOV.beta_6loc_cell_summary = [part6_summary_mulFOV.beta_6loc_cell_summary; {beta_6loc_summary}];
        part6_summary_mulFOV.r2_6loc_cell_summary = [part6_summary_mulFOV.r2_6loc_cell_summary; {r2_6loc_summary}];
        
                
        part6_summary_mulFOV.betaAngle_6loc_halfAB_norm_unitView_summary = [part6_summary_mulFOV.betaAngle_6loc_halfAB_norm_unitView_summary; betaAngle_6loc_halfAB_norm_unitView_summary];
        part6_summary_mulFOV.betaAngle_gProb_halfAB_norm_unitView_summary = [part6_summary_mulFOV.betaAngle_gProb_halfAB_norm_unitView_summary; betaAngle_gProb_halfAB_norm_unitView_summary];
        part6_summary_mulFOV.betaAngle_choiceMemory_baseline_halfAB_norm_unitView_summary = [part6_summary_mulFOV.betaAngle_choiceMemory_baseline_halfAB_norm_unitView_summary; betaAngle_choiceMemory_baseline_halfAB_norm_unitView_summary];
                
        
        part6_summary_mulFOV.r2_prior_linearAxis_summary = [part6_summary_mulFOV.r2_prior_linearAxis_summary; r2_prior_linearAxis_summary];
        part6_summary_mulFOV.r2_memory_linearAxis_summary = [part6_summary_mulFOV.r2_memory_linearAxis_summary; r2_memory_linearAxis_summary];
        part6_summary_mulFOV.r2_meta_linearAxis_summary = [part6_summary_mulFOV.r2_meta_linearAxis_summary; r2_meta_linearAxis_summary];
        
        part6_summary_mulFOV.beta_prior_linearAxis_summary = [part6_summary_mulFOV.beta_prior_linearAxis_summary; beta_prior_linearAxis_summary];
        part6_summary_mulFOV.beta_memory_linearAxis_summary = [part6_summary_mulFOV.beta_memory_linearAxis_summary; beta_memory_linearAxis_summary];
        part6_summary_mulFOV.beta_meta_linearAxis_summary = [part6_summary_mulFOV.beta_meta_linearAxis_summary; beta_meta_linearAxis_summary];

        part6_summary_mulFOV.temp_degree_aa_prior_summary = [part6_summary_mulFOV.temp_degree_aa_prior_summary; temp_degree_aa_prior_summary];
        part6_summary_mulFOV.temp_degree_bb_memory_summary = [part6_summary_mulFOV.temp_degree_bb_memory_summary; temp_degree_bb_memory_summary];
        part6_summary_mulFOV.temp_degree_cc_meta_summary = [part6_summary_mulFOV.temp_degree_cc_meta_summary; temp_degree_cc_meta_summary];

        part6_summary_mulFOV.temp_degree_ab_priorMemory_summary = [part6_summary_mulFOV.temp_degree_ab_priorMemory_summary; temp_degree_ab_priorMemory_summary];
        part6_summary_mulFOV.temp_degree_ac_priorMeta_summary = [part6_summary_mulFOV.temp_degree_ac_priorMeta_summary; temp_degree_ac_priorMeta_summary];
        part6_summary_mulFOV.temp_degree_bc_memoryMeta_summary = [part6_summary_mulFOV.temp_degree_bc_memoryMeta_summary; temp_degree_bc_memoryMeta_summary];

        
        part6_summary_mulFOV.temp_VAF_ratio_priorPrior_summary = [part6_summary_mulFOV.temp_VAF_ratio_priorPrior_summary; temp_VAF_ratio_priorPrior_summary];
        part6_summary_mulFOV.temp_VAF_ratio_memoryMemory_summary = [part6_summary_mulFOV.temp_VAF_ratio_memoryMemory_summary; temp_VAF_ratio_memoryMemory_summary];
        part6_summary_mulFOV.temp_VAF_ratio_metaMeta_summary = [part6_summary_mulFOV.temp_VAF_ratio_metaMeta_summary; temp_VAF_ratio_metaMeta_summary];
        part6_summary_mulFOV.temp_VAF_ratio_priorMemory_summary = [part6_summary_mulFOV.temp_VAF_ratio_priorMemory_summary; temp_VAF_ratio_priorMemory_summary];
        part6_summary_mulFOV.temp_VAF_ratio_priorMeta_summary = [part6_summary_mulFOV.temp_VAF_ratio_priorMeta_summary; temp_VAF_ratio_priorMeta_summary];
        part6_summary_mulFOV.temp_VAF_ratio_memoryMeta_summary = [part6_summary_mulFOV.temp_VAF_ratio_memoryMeta_summary; temp_VAF_ratio_memoryMeta_summary];
        

        part6_summary_mulFOV.temp_EV_neuronPrior_summary = [part6_summary_mulFOV.temp_EV_neuronPrior_summary; temp_EV_neuronPrior_summary];
        part6_summary_mulFOV.temp_EV_neuronMemory_summary = [part6_summary_mulFOV.temp_EV_neuronMemory_summary; temp_EV_neuronMemory_summary];
        part6_summary_mulFOV.temp_EV_neuronMeta_summary = [part6_summary_mulFOV.temp_EV_neuronMeta_summary; temp_EV_neuronMeta_summary];        
        part6_summary_mulFOV.pca_explained_baseline_cell_summary = [part6_summary_mulFOV.pca_explained_baseline_cell_summary; {pca_explained_baseline_summary}];
        part6_summary_mulFOV.pca_explained_delay_cell_summary = [part6_summary_mulFOV.pca_explained_delay_cell_summary; {pca_explained_delay_summary}];
        part6_summary_mulFOV.pca_r_eachPC_VS_eventVector_baseline_cell_summary = [part6_summary_mulFOV.pca_r_eachPC_VS_eventVector_baseline_cell_summary; {pca_r_eachPC_VS_eventVector_baseline_summary}];
        part6_summary_mulFOV.pca_r_eachPC_VS_eventVector_delay_cell_summary = [part6_summary_mulFOV.pca_r_eachPC_VS_eventVector_delay_cell_summary; {pca_r_eachPC_VS_eventVector_delay_summary}];
        
        
        part6_summary_mulFOV.clusterIndex_memoryPrecision_disBin_mean_summary = [part6_summary_mulFOV.clusterIndex_memoryPrecision_disBin_mean_summary clusterIndex_memoryPrecision_disBin_mean_summary];
        part6_summary_mulFOV.clusterIndex_choiceMemory_disBin_mean_summary = [part6_summary_mulFOV.clusterIndex_choiceMemory_disBin_mean_summary clusterIndex_choiceMemory_disBin_mean_summary];
        part6_summary_mulFOV.clusterIndex_choiceMemory_baseline_disBin_mean_summary = [part6_summary_mulFOV.clusterIndex_choiceMemory_baseline_disBin_mean_summary clusterIndex_choiceMemory_baseline_disBin_mean_summary];
        part6_summary_mulFOV.prctile_high_memoryPrecision_shuffled_summary = [part6_summary_mulFOV.prctile_high_memoryPrecision_shuffled_summary prctile_high_memoryPrecision_shuffled_summary];
        part6_summary_mulFOV.prctile_high_choiceMemory_shuffled_summary = [part6_summary_mulFOV.prctile_high_choiceMemory_shuffled_summary prctile_high_choiceMemory_shuffled_summary];
        part6_summary_mulFOV.prctile_high_choiceMemory_baseline_shuffled_summary = [part6_summary_mulFOV.prctile_high_choiceMemory_baseline_shuffled_summary prctile_high_choiceMemory_baseline_shuffled_summary];
        
        part6_summary_mulFOV.clusterIndex_memoryPrecision_shuffled_meanB_summary = [part6_summary_mulFOV.clusterIndex_memoryPrecision_shuffled_meanB_summary clusterIndex_memoryPrecision_shuffled_meanB_summary];
        part6_summary_mulFOV.clusterIndex_choiceMemory_shuffled_meanB_summary = [part6_summary_mulFOV.clusterIndex_choiceMemory_shuffled_meanB_summary clusterIndex_choiceMemory_shuffled_meanB_summary];
        part6_summary_mulFOV.clusterIndex_choiceMemory_baseline_shuffled_meanB_summary = [part6_summary_mulFOV.clusterIndex_choiceMemory_baseline_shuffled_meanB_summary clusterIndex_choiceMemory_baseline_shuffled_meanB_summary];
                
        part6_summary_mulFOV.centriodDis_all_shuffled_prctileA_summary = [part6_summary_mulFOV.centriodDis_all_shuffled_prctileA_summary centriodDis_all_shuffled_prctileA_summary];
        part6_summary_mulFOV.centriodDis_all_shuffled_prctileB_summary = [part6_summary_mulFOV.centriodDis_all_shuffled_prctileB_summary centriodDis_all_shuffled_prctileB_summary];
        
        part6_summary_mulFOV.centriodDis_all_shuffled_mean_summary = [part6_summary_mulFOV.centriodDis_all_shuffled_mean_summary centriodDis_all_shuffled_mean_summary];
        
        part6_summary_mulFOV.centriodDis_all_summary = [part6_summary_mulFOV.centriodDis_all_summary centriodDis_all_summary];
        part6_summary_mulFOV.centriodDis_all_AC_shuffled_prctileA_summary = [part6_summary_mulFOV.centriodDis_all_AC_shuffled_prctileA_summary centriodDis_all_AC_shuffled_prctileA_summary];
        part6_summary_mulFOV.centriodDis_all_AC_shuffled_prctileB_summary = [part6_summary_mulFOV.centriodDis_all_AC_shuffled_prctileB_summary centriodDis_all_AC_shuffled_prctileB_summary];
        
        part6_summary_mulFOV.centriodDis_all_AC_shuffled_mean_summary = [part6_summary_mulFOV.centriodDis_all_AC_shuffled_mean_summary centriodDis_all_AC_shuffled_mean_summary];
        
        part6_summary_mulFOV.centriodDis_AC_all_summary = [part6_summary_mulFOV.centriodDis_AC_all_summary centriodDis_AC_all_summary];
        part6_summary_mulFOV.centriodDis_all_BC_shuffled_prctileA_summary = [part6_summary_mulFOV.centriodDis_all_BC_shuffled_prctileA_summary centriodDis_all_BC_shuffled_prctileA_summary];
        part6_summary_mulFOV.centriodDis_all_BC_shuffled_prctileB_summary = [part6_summary_mulFOV.centriodDis_all_BC_shuffled_prctileB_summary centriodDis_all_BC_shuffled_prctileB_summary];
        
        part6_summary_mulFOV.centriodDis_all_BC_shuffled_mean_summary = [part6_summary_mulFOV.centriodDis_all_BC_shuffled_mean_summary centriodDis_all_BC_shuffled_mean_summary];
        
        part6_summary_mulFOV.centriodDis_BC_all_summary = [part6_summary_mulFOV.centriodDis_BC_all_summary centriodDis_BC_all_summary];
        
        
        part6_summary_mulFOV.spatialDisNew_1to1_mean_summary = [part6_summary_mulFOV.spatialDisNew_1to1_mean_summary {spatialDisNew_1to1_mean_summary}];
        part6_summary_mulFOV.spatialDisNew_2to2_mean_summary = [part6_summary_mulFOV.spatialDisNew_2to2_mean_summary {spatialDisNew_2to2_mean_summary}];
        part6_summary_mulFOV.spatialDisNew_3to3_mean_summary = [part6_summary_mulFOV.spatialDisNew_3to3_mean_summary {spatialDisNew_3to3_mean_summary}];
        
        part6_summary_mulFOV.spatialDisNew_1toOthers_mean_summary = [part6_summary_mulFOV.spatialDisNew_1toOthers_mean_summary {spatialDisNew_1toOthers_mean_summary}];
        part6_summary_mulFOV.spatialDisNew_2toOthers_mean_summary = [part6_summary_mulFOV.spatialDisNew_2toOthers_mean_summary {spatialDisNew_2toOthers_mean_summary}];
        part6_summary_mulFOV.spatialDisNew_3toOthers_mean_summary = [part6_summary_mulFOV.spatialDisNew_3toOthers_mean_summary {spatialDisNew_3toOthers_mean_summary}];
        
        part6_summary_mulFOV.spatialCentriodDis_1toOthers_summary = [part6_summary_mulFOV.spatialCentriodDis_1toOthers_summary spatialCentriodDis_1toOthers_summary];
        part6_summary_mulFOV.spatialCentriodDis_2toOthers_summary = [part6_summary_mulFOV.spatialCentriodDis_2toOthers_summary spatialCentriodDis_2toOthers_summary];
        part6_summary_mulFOV.spatialCentriodDis_3toOthers_summary = [part6_summary_mulFOV.spatialCentriodDis_3toOthers_summary spatialCentriodDis_3toOthers_summary];
        
    end
    
    
    
    %% Tail of a loop
    fprintf('Progress %d/%d, t=%.1f secs.\n',tempIndexFOV_AB,numFOV_AB_summary,toc(t0_runScript));
    
    %part2_summary;
    %part3_summary;
    %part4_summary;
    %part5_summary;
    %part6_summary;
    
    
    if if_singleFOV0_multiFOV1 == 0
        break
    end
    
end

fprintf('Time is %.1f secs.\n',toc(t0_runScript));


part2_summary_mulFOV;
part3_summary_mulFOV;
part4_summary_mulFOV;
part5_summary_mulFOV;
part6_summary_mulFOV;

temp_cellIndexMapping_AB_mulFOV;

allPart_summary_mulFOV = struct;
allPart_summary_mulFOV.part2_summary_mulFOV = part2_summary_mulFOV;
allPart_summary_mulFOV.part3_summary_mulFOV = part3_summary_mulFOV;
allPart_summary_mulFOV.part4_summary_mulFOV = part4_summary_mulFOV;
allPart_summary_mulFOV.part5_summary_mulFOV = part5_summary_mulFOV;
allPart_summary_mulFOV.part6_summary_mulFOV = part6_summary_mulFOV;


allPart_summary_mulFOV.temp_cellIndexMapping_AB_mulFOV = temp_cellIndexMapping_AB_mulFOV;

allPart_summary_mulFOV.p_tuning = struct;

allPart_summary_mulFOV.p_tuning.p_6loc_mulFOV = p_6loc_mulFOV;
allPart_summary_mulFOV.p_tuning.p_memoryPrecision_mulFOV = p_memoryPrecision_mulFOV;
allPart_summary_mulFOV.p_tuning.p_choiceMemory_mulFOV = p_choiceMemory_mulFOV;
allPart_summary_mulFOV.p_tuning.p_choiceMemory_baseline_mulFOV = p_choiceMemory_baseline_mulFOV;

allPart_summary_mulFOV.p_tuning.p_seqPrecision_mulFOV = p_seqPrecision_mulFOV;
allPart_summary_mulFOV.p_tuning.p_gProb_mulFOV = p_gProb_mulFOV;
allPart_summary_mulFOV.p_tuning.p_gProb_baseline_mulFOV = p_gProb_baseline_mulFOV;

% temp_cellIndexMapping_AB_mulFOV_collapsed = [];
% for tempi=1:length(temp_cellIndexMapping_AB_mulFOV)
%     temp_cellIndexMapping_AB_mulFOV_collapsed = [temp_cellIndexMapping_AB_mulFOV_collapsed; temp_cellIndexMapping_AB_mulFOV{tempi}];
% end

% a1 = tempBoolIndex123_or(temp_cellIndexMapping_AB_mulFOV_collapsed);

% allPart_summary_mulFOV.p_tuning.p_loc_mulFOV = tempBoolIndex123_or(temp_cellIndexMapping_AB_mulFOV_collapsed);

boolIndex_loc_mulFOV = [];
for tempi=1:length(temp_cellIndexMapping_AB_mulFOV)    
    boolIndex_loc = tempBoolIndex123_or(temp_cellIndexMapping_AB_mulFOV{tempi});
    boolIndex_loc_mulFOV = [boolIndex_loc_mulFOV;{boolIndex_loc}];
end

allPart_summary_mulFOV.p_tuning.boolIndex_loc_mulFOV = boolIndex_loc_mulFOV;



if false
   allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_correct_theoretical_summary = ...
       part2_summary_mulFOV.r_mean_correct_theoretical_summary; %#ok<*UNRCH>  
end

if false
    allPart_summary_mulFOV.part3_summary_mulFOV.seqMeta_memoryPrecisionLow_choiceStimuli_mean_summary = part3_summary_mulFOV.seqMeta_memoryPrecisionLow_choiceStimuli_mean_summary;
    allPart_summary_mulFOV.part3_summary_mulFOV.seqMeta_memoryPrecisionHigh_choiceStimuli_mean_summary = part3_summary_mulFOV.seqMeta_memoryPrecisionHigh_choiceStimuli_mean_summary;
    allPart_summary_mulFOV.part3_summary_mulFOV.p_seqMeta_precisionLowHigh_stimuli_summary = part3_summary_mulFOV.p_seqMeta_precisionLowHigh_stimuli_summary;
    allPart_summary_mulFOV.part3_summary_mulFOV.seqMeta_memoryPrecisionLow_choiceResponse_mean_summary = part3_summary_mulFOV.seqMeta_memoryPrecisionLow_choiceResponse_mean_summary;
    allPart_summary_mulFOV.part3_summary_mulFOV.seqMeta_memoryPrecisionHigh_choiceResponse_mean_summary = part3_summary_mulFOV.seqMeta_memoryPrecisionHigh_choiceResponse_mean_summary;
    allPart_summary_mulFOV.part3_summary_mulFOV.p_seqMeta_precisionLowHigh_response_summary = part3_summary_mulFOV.p_seqMeta_precisionLowHigh_response_summary;
    
    allPart_summary_mulFOV.part3_summary_mulFOV.seqCorr_precisionMeta_choiceStimuli_mean_summary = part3_summary_mulFOV.seqCorr_precisionMeta_choiceStimuli_mean_summary;
    allPart_summary_mulFOV.part3_summary_mulFOV.seqCorr_precisionMeta_choiceResponse_mean_summary = part3_summary_mulFOV.seqCorr_precisionMeta_choiceResponse_mean_summary;

end

%% End