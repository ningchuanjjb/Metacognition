threshold_locRatio = 0.055;%0.06-->0.05-->0.045
threshold_locStd = 0.026;%0.0278-->0.025-->0.026

fun_locBasedResample_seqCount_Name_v = autoGetFunName_myScripts('fun_locBasedResample_seqCount', [targetPATH '\functions']);
fun_locBasedResample_seqCount = str2func(fun_locBasedResample_seqCount_Name_v);

options_locBased.threshold_locRatio = threshold_locRatio;
options_locBased.threshold_locStd = threshold_locStd;
options_locBased.numSeq = numSeq;
options_locBased.numFrames = numFrames;
options_locBased.seqIndex_valid = seqIndex_valid;
options_locBased.boolIndex_location_seq_T = boolIndex_location_seq_T;

resampleTrialCount_seq = fun_locBasedResample_seqCount(options_locBased);

fun_resampleXY_seqCountBased_Name_v = autoGetFunName_myScripts('fun_resampleXY_seqCountBased', [targetPATH '\functions']);
fun_resampleXY_seqCountBased = str2func(fun_resampleXY_seqCountBased_Name_v);

options_resampleXY.numSeq = numSeq;
options_resampleXY.resampleIterCount = resampleIterCount;
options_resampleXY.resampleTrialCount_seq = resampleTrialCount_seq;
options_resampleXY.seqIndex_valid = seqIndex_valid;
options_resampleXY.temp_F_dff_decisionBin = temp_F_dff_decisionBin;

output_resampleXY = fun_resampleXY_seqCountBased(options_resampleXY);

temp_F_dff_decisionBin_resample = output_resampleXY.temp_F_dff_decisionBin_resample;
seqIndex_valid_resample = output_resampleXY.seqIndex_valid_resample;




