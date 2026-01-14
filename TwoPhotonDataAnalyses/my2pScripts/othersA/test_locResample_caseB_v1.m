seqCount1;
temp_svm_Y_valid;

seq_length = sum(temp_svm_Y_valid,1);
trial_num = size(temp_svm_Y_valid,2);
locCount_length = zeros(4,numFrames);
for tempi=1:trial_num
    locCount_length(seq_length(tempi),temp_svm_Y_valid(:,tempi)) = ...
        locCount_length(seq_length(tempi),temp_svm_Y_valid(:,tempi)) + 1;
end

numSeq = [6,15,20,15];
resampleCount_seq = zeros(1,sum(numSeq));
for tempi=1:4
    temp_range = (sum(numSeq(1:(tempi-1)))+1):sum(numSeq(1:tempi));
        
    if tempi == 1
        resampleCount_seq(temp_range) = min(locCount_length(tempi,:));
    end
    
    if tempi == 2
        temp_trial_num = sum(seq_length == 2);
        temp_locCount_length = locCount_length(tempi,:);
        
        
    end
    
end