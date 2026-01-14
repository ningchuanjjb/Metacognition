
% trialIndex_bool_allMemory = ~(trial_para.ifSelectOffloading==1);
% % trialIndex_bool_offload = trial_para.ifSelectOffloading==1;
% F_dff_decisionBin1;
% 
% gAcc_allLength = cell(1,3);
% seqCount_allLength_memoryCorrect = cell(1,3);
% seqCount_allLength_allMemory = cell(1,3);
% seqCount_allLength_gAcc = cell(1,3);
% for target_length=1:3
%     
%     temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
%     gAcc = gAcc_target_collapsed_inOne(temp_range);
%     
%     svm_options.gAcc = gAcc;
%     svm_options.seq_range = temp_range;
%     svm_options.target_length = target_length;
%     
%     tempBoolIndex = ismember(seqIndex_valid,temp_range);
%     x = temp_F_dff_decisionBin(:,tempBoolIndex);
%     y = boolIndex_location_trial(:,tempBoolIndex);
%     
%     
%     temp_uniqueSeq = unique(seqIndex_valid(tempBoolIndex));
%     validSeqBoolIndex = ismember(temp_range,temp_uniqueSeq);
%     
%     temp_seqCount_memoryCorrect = zeros(1,numSeq(target_length));
%     temp_seqCount_allMemory = zeros(1,numSeq(target_length));    
%     for tempi=1:numSeq(target_length)
%         tempj = temp_range(tempi);
%         temp_seqCount_memoryCorrect(tempi) = sum(seqIndex(temp_trialIndex_bool)==tempj,'all');
%         temp_seqCount_allMemory(tempi) = sum(seqIndex(trialIndex_bool_allMemory)==tempj,'all');        
%     end
%     
%     seqCount_allLength_memoryCorrect{target_length} = temp_seqCount_memoryCorrect;
%     seqCount_allLength_allMemory{target_length} = temp_seqCount_allMemory;    
%     seqCount_allLength_gAcc{target_length} = temp_seqCount_memoryCorrect./temp_seqCount_allMemory;    
%     gAcc_allLength{target_length} = gAcc;
% end
a = 1;
%%
for target_length=1:3
    
    if target_length == 1
        svm_posterior_lengthx = svm_posterior_length1;        
    elseif target_length == 2
        svm_posterior_lengthx = svm_posterior_length2;        
    elseif target_length == 3        
        svm_posterior_lengthx = svm_posterior_length3;
    end
    
    fprintf('seqCount_memoryCorrect{%d}, seqCount_allMemory{%d}, seqCount_allLength_gAcc{%d}, svm_posterior_length%d, gAcc_allLength{%d}: \n',target_length,target_length,target_length,target_length,target_length);
    temp_seqNum = length(seqCount_allLength_memoryCorrect{target_length});
    for tempi=1:temp_seqNum
        fprintf('  %2d, ',seqCount_allLength_memoryCorrect{target_length}(tempi));
    end
    fprintf('\n');
    for tempi=1:temp_seqNum
        fprintf('  %2d, ',seqCount_allLength_allMemory{target_length}(tempi));
    end
    fprintf('\n');
    for tempi=1:temp_seqNum
        fprintf('%.2f, ',seqCount_allLength_gAcc{target_length}(tempi));
    end
    fprintf('\n');
    for tempi=1:temp_seqNum
        if isnan(svm_posterior_lengthx(tempi)) == false
            fprintf('%.2f, ',svm_posterior_lengthx(tempi));
        else
            fprintf(' NaN, ');
        end
    end
    fprintf('\n');
    for tempi=1:temp_seqNum
        fprintf('%.2f, ',gAcc_allLength{target_length}(tempi));
    end
    fprintf('\n');    
    
end

%%

