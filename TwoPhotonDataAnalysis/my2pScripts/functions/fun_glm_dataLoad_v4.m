function glm_dataLoad_output = fun_glm_dataLoad_v4(glm_dataLoad_options)

tempData = glm_dataLoad_options.tempData;
if_load_allTrial0_memoryCorrect1 = glm_dataLoad_options.if_load_allTrial0_memoryCorrect1;
plot_lengthFlag = glm_dataLoad_options.plot_lengthFlag;

if if_load_allTrial0_memoryCorrect1 == 1
    sequence_length1_onehot = tempData.sequence_length1_memoryCorrect_onehot;
    sequence_length2_onehot = tempData.sequence_length2_memoryCorrect_onehot;
    sequence_length3_onehot = tempData.sequence_length3_memoryCorrect_onehot;
    sequence_length4_onehot = tempData.sequence_length4_memoryCorrect_onehot;
    responseSeq_length1_onehot = tempData.sequence_length1_memoryCorrect_onehot;
    responseSeq_length2_onehot = tempData.sequence_length2_memoryCorrect_onehot;
    responseSeq_length3_onehot = tempData.sequence_length3_memoryCorrect_onehot;
    responseSeq_length4_onehot = tempData.sequence_length4_memoryCorrect_onehot;
elseif if_load_allTrial0_memoryCorrect1 == 0
    sequence_length1_onehot = tempData.sequence_length1_onehot;
    sequence_length2_onehot = tempData.sequence_length2_onehot;
    sequence_length3_onehot = tempData.sequence_length3_onehot;
    sequence_length4_onehot = tempData.sequence_length4_onehot;
    responseSeq_length1_onehot = tempData.responseSeq_length1_onehot;
    responseSeq_length2_onehot = tempData.responseSeq_length2_onehot;
    responseSeq_length3_onehot = tempData.responseSeq_length3_onehot;
    responseSeq_length4_onehot = tempData.responseSeq_length4_onehot;
end


if plot_lengthFlag == 1
    F_dff_lengthx_delay1Bin = tempData.F_dff_length1_delay1Bin;
    F_dff_lengthx_delay2Bin = tempData.F_dff_length1_delay2Bin;
    F_dff_lengthx_baselineBin = tempData.F_dff_length1_baselineBin;
    sequence_lengthx_onehot = sequence_length1_onehot;
    responseSeq_lengthx_onehot = responseSeq_length1_onehot;
elseif plot_lengthFlag == 2
    F_dff_lengthx_delay1Bin = tempData.F_dff_length2_delay1Bin;
    F_dff_lengthx_delay2Bin = tempData.F_dff_length2_delay2Bin;
    F_dff_lengthx_baselineBin = tempData.F_dff_length2_baselineBin;
    sequence_lengthx_onehot = sequence_length2_onehot;
    responseSeq_lengthx_onehot = responseSeq_length2_onehot;
elseif plot_lengthFlag == 3
    F_dff_lengthx_delay1Bin = tempData.F_dff_length3_delay1Bin;
    F_dff_lengthx_delay2Bin = tempData.F_dff_length3_delay2Bin;
    F_dff_lengthx_baselineBin = tempData.F_dff_length3_baselineBin;
    sequence_lengthx_onehot = sequence_length3_onehot;
    responseSeq_lengthx_onehot = responseSeq_length3_onehot;
elseif plot_lengthFlag == 4
    F_dff_lengthx_delay1Bin = tempData.F_dff_length4_delay1Bin;
    F_dff_lengthx_delay2Bin = tempData.F_dff_length4_delay2Bin;
    F_dff_lengthx_baselineBin = tempData.F_dff_length4_baselineBin;
    sequence_lengthx_onehot = sequence_length4_onehot;
    responseSeq_lengthx_onehot = responseSeq_length4_onehot;    
elseif plot_lengthFlag == 0
    F_dff_lengthx_delay1Bin = [tempData.F_dff_length1_delay1Bin tempData.F_dff_length2_delay1Bin tempData.F_dff_length3_delay1Bin tempData.F_dff_length4_delay1Bin];
    F_dff_lengthx_delay2Bin = [tempData.F_dff_length1_delay2Bin tempData.F_dff_length2_delay2Bin tempData.F_dff_length3_delay2Bin tempData.F_dff_length4_delay2Bin];
    F_dff_lengthx_baselineBin = [tempData.F_dff_length1_baselineBin tempData.F_dff_length2_baselineBin tempData.F_dff_length3_baselineBin tempData.F_dff_length4_baselineBin];
    sequence_lengthx_onehot = [tempData.sequence_length1_onehot; tempData.sequence_length2_onehot; tempData.sequence_length3_onehot; tempData.sequence_length4_onehot];
    responseSeq_lengthx_onehot = [tempData.responseSeq_length1_onehot; tempData.responseSeq_length2_onehot; tempData.responseSeq_length3_onehot; tempData.responseSeq_length4_onehot];    
end

glm_dataLoad_output = struct;
glm_dataLoad_output.F_dff_lengthx_delay1Bin = F_dff_lengthx_delay1Bin;
glm_dataLoad_output.F_dff_lengthx_delay2Bin = F_dff_lengthx_delay2Bin;
glm_dataLoad_output.F_dff_lengthx_baselineBin = F_dff_lengthx_baselineBin;
glm_dataLoad_output.sequence_lengthx_onehot = sequence_lengthx_onehot;
glm_dataLoad_output.responseSeq_lengthx_onehot = responseSeq_lengthx_onehot;

