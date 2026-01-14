fun_rProb_glm_frame_Name_v = autoGetFunName_myScripts('fun_rProb_glm_frame', [targetPATH '\functions']);
fun_rProb_glm_frame = str2func(fun_rProb_glm_frame_Name_v);

rProb_glm_options = struct;
rProb_glm_options.offloadingProb_inOne = offloadingProb_inOne;
rProb_glm_options.numSeq = numSeq;
rProb_glm_options.target_seqSet = target_seqSet;
rProb_glm_options.trial_para = trial_para;


%temp_F_dff = F_dff_baseline;
temp_F_dff = F_dff_delay1;
% temp_F_dff = F_dff_delay2;

rProb_glm_options.temp_F_dff = temp_F_dff;
rProb_glm_output = fun_rProb_glm_frame(rProb_glm_options);

selectiveCellBoolIndex_rProb_glm = rProb_glm_output.selectiveCellBoolIndex_rProb_glm;
temp_sum = sum(rProb_glm_output.selectiveCellBoolIndex_rProb_glm);
a = 1;


