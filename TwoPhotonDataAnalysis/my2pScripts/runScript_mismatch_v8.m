% Chuan's 18th script (20251214)
%% Initialization
% clear
close all

home;

if_monkey_D0_Z1 = 0;%1

if_runScript_stepF1 = 1;
if_runScript_stepF2 = 1;
if_runScript_stepD3 = 1;
if_runScript_stepF4 = 1;
if_runScript_stepDM1 = 1;
if_runScript_stepF5A = 1;
if_runScript_stepF5B = 1;
if_runScript_stepF5C = 0;
if_runScript_stepF5D = 0;
if_runScript_stepF5E = 0;% Time-consuming


targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)

%% Preparation
stepF1_test_glm_locProd_eachLength_multiFOV_Name_v = autoGetFunName_myScripts('stepF1_test_glm_locProd_eachLength_multiFOV', targetPATH);
script_stepF1_test_glm_locProd_eachLength_multiFOV = str2func(stepF1_test_glm_locProd_eachLength_multiFOV_Name_v);

stepF2_test_locTuning_Name_v = autoGetFunName_myScripts('stepF2_test_locTuning', targetPATH);
script_stepF2_test_locTuning = str2func(stepF2_test_locTuning_Name_v);

stepD3_train123_test123_Name_v = autoGetFunName_myScripts('stepD3_train123_test123', targetPATH);
script_stepD3_train123_test123 = str2func(stepD3_train123_test123_Name_v);

stepF4_memoryMetaMismatch_singleNeuron_Name_v = autoGetFunName_myScripts('stepF4A_memoryMetaMismatch_singleNeuron', targetPATH);
script_stepF4_memoryMetaMismatch_singleNeuron = str2func(stepF4_memoryMetaMismatch_singleNeuron_Name_v);

stepDM1_test_decodingMeta_Name_v = autoGetFunName_myScripts('stepDM1_test_decodingMeta', targetPATH);
script_stepDM1_test_decodingMeta = str2func(stepDM1_test_decodingMeta_Name_v);

stepF5A_memoryMetaMismatch_twoDecoder_Name_v = autoGetFunName_myScripts('stepF5A_memoryMetaMismatch_twoDecoder', targetPATH);
script_stepF5A_memoryMetaMismatch_twoDecoder = str2func(stepF5A_memoryMetaMismatch_twoDecoder_Name_v);

stepF5B_test_trialHistory_Name_v = autoGetFunName_myScripts('stepF5B_test_trialHistory', targetPATH);
script_stepF5B_test_trialHistory = str2func(stepF5B_test_trialHistory_Name_v);

stepF5C_test_mismatchMechnism_bin_Name_v = autoGetFunName_myScripts('stepF5C_test_mismatchMechnism_bin', targetPATH);
script_stepF5C_test_mismatchMechnism_bin = str2func(stepF5C_test_mismatchMechnism_bin_Name_v);

stepF5D_test_mismatchMechnism_meta_crossTime_Name_v = autoGetFunName_myScripts('stepF5D_test_mismatchMechnism_meta_crossTime', targetPATH);
script_stepF5D_test_mismatchMechnism_meta_crossTime = str2func(stepF5D_test_mismatchMechnism_meta_crossTime_Name_v);

stepF5E_test_mismatchMechnism_memoryPrecision_crossTime_Name_v = autoGetFunName_myScripts('stepF5E_test_mismatchMechnism_memoryPrecision_crossTime', targetPATH);
script_stepF5E_test_mismatchMechnism_memoryPrecision_crossTime = str2func(stepF5E_test_mismatchMechnism_memoryPrecision_crossTime_Name_v);


%% Run script
t0_runScript = tic;
if if_runScript_stepF1 == 1
    fprintf('\nNow runing %s.  ------> \n', stepF1_test_glm_locProd_eachLength_multiFOV_Name_v);
    script_stepF1_test_glm_locProd_eachLength_multiFOV();
end
if if_runScript_stepF2 == 1
    fprintf('\nNow runing %s.  ------> \n', stepF2_test_locTuning_Name_v);
    script_stepF2_test_locTuning();
end
if if_runScript_stepD3 == 1
    fprintf('\nNow runing %s.  ------> \n', stepD3_train123_test123_Name_v);
    script_stepD3_train123_test123();
end

if if_runScript_stepF4 == 1
    fprintf('\nNow runing %s.  ------> \n', stepF4_memoryMetaMismatch_singleNeuron_Name_v);    
    script_stepF4_memoryMetaMismatch_singleNeuron();        
end
if if_runScript_stepDM1 == 1
    fprintf('\nNow runing %s.  ------> \n', stepDM1_test_decodingMeta_Name_v);
    
    if_trainMeta_0baseline_1delay1 = 0; %#ok<*NASGU>
    script_stepDM1_test_decodingMeta();
    
    if_trainMeta_0baseline_1delay1 = 1;
    script_stepDM1_test_decodingMeta();    
end
if if_runScript_stepF5A == 1
    fprintf('\nNow runing %s.  ------> \n', stepF5A_memoryMetaMismatch_twoDecoder_Name_v);
    script_stepF5A_memoryMetaMismatch_twoDecoder();
end
if if_runScript_stepF5B == 1
    fprintf('\nNow runing %s.  ------> \n', stepF5B_test_trialHistory_Name_v);
    script_stepF5B_test_trialHistory();
end
if if_runScript_stepF5C == 1
    fprintf('\nNow runing %s.  ------> \n', stepF5C_test_mismatchMechnism_bin_Name_v);
    script_stepF5C_test_mismatchMechnism_bin();
end
if if_runScript_stepF5D == 1
    fprintf('\nNow runing %s.  ------> \n', stepF5D_test_mismatchMechnism_meta_crossTime_Name_v);
    script_stepF5D_test_mismatchMechnism_meta_crossTime();
end
if if_runScript_stepF5E == 1
    fprintf('\nNow runing %s.  ------> \n', stepF5E_test_mismatchMechnism_memoryPrecision_crossTime_Name_v);
    script_stepF5E_test_mismatchMechnism_memoryPrecision_crossTime();
end


fprintf('Time is %.1f secs.\n',toc(t0_runScript));




%% End