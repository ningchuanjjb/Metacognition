function target_seqSet = get_target_seqSet_v1(numFrames,pointKindsNum)

% numFrames = 6;
% pointKindsNum = 4;
target_seqSet = cell(1, pointKindsNum);

for target_seqLength=1:4
    
    target_numSeq = nchoosek(numFrames, target_seqLength);
    
    target_seqSet{target_seqLength}= cell(target_numSeq, 1);
    temp_index = 0;
    if target_seqLength == 1
        for tempi=1:numFrames+1-target_seqLength
            temp_index = temp_index + 1;
            target_seqSet{target_seqLength}{temp_index} = [tempi]; %#ok<*NBRAK>
        end
    elseif target_seqLength == 2
        for tempi=1:numFrames+1-target_seqLength
            for tempj=tempi+1:numFrames+1-target_seqLength+1
                temp_index = temp_index + 1;
                target_seqSet{target_seqLength}{temp_index} = [tempi tempj];
            end
        end
    elseif target_seqLength == 3
        for tempi=1:numFrames+1-target_seqLength
            for tempj=tempi+1:numFrames+1-target_seqLength+1
                for tempk=tempj+1:numFrames+1-target_seqLength+2
                    temp_index = temp_index + 1;
                    target_seqSet{target_seqLength}{temp_index} = [tempi tempj tempk];
                end
            end
        end
        
    elseif target_seqLength == 4
        for tempi=1:numFrames+1-target_seqLength
            for tempj=tempi+1:numFrames+1-target_seqLength+1
                for tempk=tempj+1:numFrames+1-target_seqLength+2
                    for tempz=tempk+1:numFrames+1-target_seqLength+3
                        temp_index = temp_index + 1;
                        target_seqSet{target_seqLength}{temp_index} = [tempi tempj tempk tempz];
                    end
                end
            end
        end
    end
end