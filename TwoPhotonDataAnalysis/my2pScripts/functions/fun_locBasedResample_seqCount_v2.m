function resampleTrialCount_seq = fun_locBasedResample_seqCount_v1(options)
% threshold_locRatio = 0.055;%0.06-->0.05-->0.045
% threshold_locStd = 0.026;%0.0278-->0.025-->0.026

threshold_locRatio = options.threshold_locRatio;
threshold_locStd = options.threshold_locStd;

numSeq = options.numSeq;
numFrames = options.numFrames;
seqIndex_valid = options.seqIndex_valid;
boolIndex_location_seq_T = options.boolIndex_location_seq_T;




seqCount_length_locResampled = zeros(1,sum(numSeq(1:3)));
for target_length=1:3
    temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));

    temp_seqCount = zeros(1,numSeq(target_length));
    for tempi=1:numSeq(target_length)
        tempj = temp_range(tempi);
        temp_seqCount(tempi) = sum(seqIndex_valid==tempj,'all');
    end
        
    temp_boolIndex_location_seq_T = boolIndex_location_seq_T(temp_range,:);
        
    temp_seqCount2 = temp_seqCount;
    
    temp1 = temp_boolIndex_location_seq_T .* temp_seqCount';
    temp_locCount = sum(temp1,1);
    temp_locCount2 = temp_locCount;
    
    for tempLoopCount=1:sum(temp_seqCount)
        [MA,IA] = sort(temp_seqCount2','descend'); %#ok<*ASGLU,*UDIM>
        [MB,IB] = sort(temp_locCount2,'descend');
        IB2 = find(temp_boolIndex_location_seq_T(:,IB(1))==1 & temp_boolIndex_location_seq_T(:,IB(numFrames))==0);

        for temptempi=1:length(IA)
            if ismember(IA(temptempi),IB2) == true
                break
            end
        end
        IA2 = IA(temptempi);
        temp_seqCount2(IA2) = temp_seqCount2(IA2) - 1;
        temp1 = temp_boolIndex_location_seq_T .* temp_seqCount2';
        temp_locCount2 = sum(temp1,1);
        a = 1;
        
        [temp_min,temp_max] = bounds(temp_locCount2);
        temp_ratio = (temp_max-temp_min)/sum(temp_locCount2);
        temp_std = std(temp_locCount2/sum(temp_locCount2));
        
        if temp_ratio <= threshold_locRatio
            if temp_std < threshold_locStd
                break
            end
        end
                        
        %test in 20240125
        if max(temp_locCount2) <= numSeq(target_length)*2
            break
        end
        
    end
    
    fprintf('length%d,trialNum(valid/raw)=%3d/%3d. ',target_length,sum(temp_seqCount2),sum(temp_seqCount));

    fprintf('locCount =[');
    fprintf('%d,',temp_locCount);
    fprintf('] --> ');
    fprintf('locCount2=[');
    fprintf('%d,',temp_locCount2);
    fprintf('].\n');
    
    seqCount_length_locResampled(temp_range) = temp_seqCount2;
end

seqCount_length_locResampled_cell = cell(3,1);
fprintf('seqCount_length_locResampled=\n');
for tempi=1:3
    fprintf('[');
    temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
    seqCount_length_locResampled_cell{tempi} = seqCount_length_locResampled(temp_range);
    fprintf('%2d,',seqCount_length_locResampled_cell{tempi});
    fprintf(']\n');
end

seqCount_length_locResampled; %#ok<*VUNUS>
resampleTrialCount_seq = seqCount_length_locResampled;