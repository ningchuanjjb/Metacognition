function features = getMappingFeaturesUnified(sequence, errorRate, numFrames_rangeTail)
%features = [seqOnehot seq_length_trial numChunk_trial rightOffset_trial ifSymmetry_trial seq_range_trial errorRate];    

% numFrames_rangeTail = 6; or 7

%sequence onehot
seqOnehot = zeros(1, numFrames_rangeTail);
seqOnehot(sequence) = 1;

%seq_length_trial
seq_length_trial = length(sequence);

%numChunk_trial
numChunk_trial = 1;
if seq_length_trial > 1
    for target_seqLength=1:seq_length_trial-1
        if sequence(target_seqLength) ~= sequence(target_seqLength+1)-1
            numChunk_trial = numChunk_trial + 1;
        end
    end
end

%rightOffset_trial
rightOffset_trial = sum(sequence)/seq_length_trial - (1+numFrames_rangeTail)/2;

%ifSymmetry_trial
ifSymmetry_trial = 0;
if sum(ismember([1 numFrames_rangeTail], sequence)) == 2 ||...
        sum(ismember([2 numFrames_rangeTail-1], sequence)) == 2 ||...
        sum(ismember([3 numFrames_rangeTail-2], sequence)) == 2 ||...
        sum(ismember([4 numFrames_rangeTail-3], sequence)) == 2
    ifSymmetry_trial = 1;
    
end

%seq_range_trial
seq_range_trial = sequence(end) - sequence(1);

features = [seqOnehot seq_length_trial numChunk_trial rightOffset_trial ifSymmetry_trial seq_range_trial errorRate];    

end