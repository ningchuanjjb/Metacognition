
if fasle
    PTB_marker_content_trial_para = trial_para.marker.content;
    
    % PTB_marker_content_edf = edf.Messages.info(11:end);
    PTB_marker_sampleIndex_edf = edf.Messages.time(11:end);
    
    PTB_marker_sampleIndex_edf4K_raw = PTB_marker_sampleIndex_edf*4;
    PTB_marker_sampleIndex_edf4K = PTB_marker_sampleIndex_edf4K_raw-PTB_marker_sampleIndex_edf4K_raw(1)+PTB_marker_sampleIndex(1);
    
    temp_diff1 = diff(PTB_marker_sampleIndex);
    temp_diff2 = diff(PTB_marker_sampleIndex_edf4K);
    temp_diff12 = [temp_diff1 temp_diff1(end)] - temp_diff2;
    
    PTB_marker_sampleIndex;
    PTB_marker_content_trial_para;
    
    compensate_sampleIndex = 137;
    compensate_interval = 44;
    
    PTB_marker_sampleIndex_compensate = nan(1,length(PTB_marker_content_trial_para));
    PTB_marker_sampleIndex_compensate(1:(compensate_sampleIndex-1)) = PTB_marker_sampleIndex(1:(compensate_sampleIndex-1));
    
    PTB_marker_sampleIndex_compensate(compensate_sampleIndex) = PTB_marker_sampleIndex(compensate_sampleIndex-1) + compensate_interval;
    
    PTB_marker_sampleIndex_compensate((compensate_sampleIndex+1):end) = PTB_marker_sampleIndex(compensate_sampleIndex:end);
    
    temp_diff1_compensate = diff(PTB_marker_sampleIndex_compensate);
    
    
    PTB_marker_sampleIndex;
    PTB_marker_sampleIndex_compensate;
    
    PTB_marker_sampleIndex_raw = PTB_marker_sampleIndex;
    PTB_marker_sampleIndex = PTB_marker_sampleIndex_compensate;    
end


