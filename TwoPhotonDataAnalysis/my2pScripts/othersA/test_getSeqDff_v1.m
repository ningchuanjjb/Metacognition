

glm_dataLoad_options = struct;
glm_dataLoad_options.tempData = tempData;
glm_dataLoad_options.if_load_allTrial0_memoryCorrect1 = if_load_allTrial0_memoryCorrect1;
glm_dataLoad_options.plot_lengthFlag = plot_lengthFlag;

glm_dataLoad_output = fun_glm_dataLoad(glm_dataLoad_options);

F_dff_lengthx_delay1Bin = glm_dataLoad_output.F_dff_lengthx_delay1Bin;
sequence_lengthx_onehot = glm_dataLoad_output.sequence_lengthx_onehot;
responseSeq_lengthx_onehot = glm_dataLoad_output.responseSeq_lengthx_onehot;

[temp_roiNum,temp_trialNum] = size(F_dff_lengthx_delay1Bin);

boolIndex_location_seq;
boolIndex_location_seq_T = boolIndex_location_seq';

seqNum = size(boolIndex_location_seq_T,1);
seqIndex = zeros(temp_trialNum,1);
for tempi=1:temp_trialNum
    temp_seq_onehot = sequence_lengthx_onehot(tempi,:);    
    for tempj=1:seqNum
        if sum(temp_seq_onehot == boolIndex_location_seq_T(tempj,:)) == numFrames
            break        
        end
    end
    seqIndex(tempi) = tempj;    
end
isCorrect = sum(sequence_lengthx_onehot == responseSeq_lengthx_onehot,2) == numFrames;

F_dff_lengthx_delay1Bin_seq = zeros(temp_roiNum,seqNum,'single');
for tempi=1:seqNum
    %temptempBoolIndex = seqIndex==tempi;
    temptempBoolIndex = (seqIndex==tempi) & isCorrect;    
    if sum(temptempBoolIndex) > 0
        F_dff_lengthx_delay1Bin_seq(:,tempi) = mean(F_dff_lengthx_delay1Bin(:,temptempBoolIndex),2);
    end
end

