function output = fun_markerDownSampling_v4(marker_X_sampleIndex,valid_frames_sampleIndex)

validRecordingBlock_range = zeros(1, 2);
validRecordingBlock_range(1) = valid_frames_sampleIndex(1);
validRecordingBlock_range(2) = valid_frames_sampleIndex(end);

valid_marker_X_trialIndex = ...
    marker_X_sampleIndex >= validRecordingBlock_range(1)...
    & marker_X_sampleIndex <= validRecordingBlock_range(2);

valid_marker_X_sampleIndex = marker_X_sampleIndex(valid_marker_X_trialIndex);

valid_marker_X_frameIndex = zeros(1,sum(valid_marker_X_trialIndex));
for tempi=1:length(valid_marker_X_frameIndex)  
    temp_frameIndex = find(...
        valid_marker_X_sampleIndex(tempi) >= valid_frames_sampleIndex,...
        1,'last');    
    valid_marker_X_frameIndex(tempi) = temp_frameIndex;
    
    % find the nearest frame for downsampling, 20230609
%     temp1 = abs(valid_frames_sampleIndex(temp_frameIndex)-valid_marker_X_sampleIndex(tempi));
%     temp2 = abs(valid_frames_sampleIndex(temp_frameIndex+1)-valid_marker_X_sampleIndex(tempi));
%     if temp1 >= temp2
%         valid_marker_X_frameIndex(tempi) = temp_frameIndex+1;
%     elseif temp1 < temp2
%         valid_marker_X_frameIndex(tempi) = temp_frameIndex;
%     end

end
valid_marker_X_frameIndex = valid_marker_X_frameIndex + 1; % compensate 1 frame for downsampling, 20230609
valid_marker_X_frameIndex = valid_marker_X_frameIndex + 0.5; % compensate 0.5 frame for 2p marker delay, 20230609
valid_marker_X_frameIndex = valid_marker_X_frameIndex + 0.5; % compensate 0.5 frame for ptb display delay, 20230609
valid_marker_X_frameIndex = valid_marker_X_frameIndex + 1; % compensate 1 frame for F rising to peak delay, 20230609

output = struct;
output.valid_marker_X_trialIndex = valid_marker_X_trialIndex;
output.valid_marker_X_sampleIndex = valid_marker_X_sampleIndex;
output.valid_marker_X_frameIndex = valid_marker_X_frameIndex;
