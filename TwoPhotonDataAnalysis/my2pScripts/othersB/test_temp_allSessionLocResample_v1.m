if_compute = 0;

%% Get resampleTrialCount_seq
if if_resample_seqBased0_locBased1 == 0
    resampleTrialCount_seq = zeros(1,sum(numSeq(1:3)));
    for tempi=1:3
        temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
        resampleTrialCount_seq(temp_range) = resampleTrialCount_length(tempi);
    end
    
elseif if_resample_seqBased0_locBased1 == 1
    if if_compute == 1
        seqCount_valid_allSessions = cell(1,num_FOV_AB);
        
        for currentSessionIndex=1:num_FOV_AB
            % for currentSessionIndex=1:2
            currentABSession = currentABSession_multi{currentSessionIndex};
            fprintf('currentABSession = %s.\n',currentABSession);
            
            temp_load = load([output_path '\' currentABSession]);
            temp_decodingDataSimplified = temp_load.decodingDataSimplified_AB;
            
            seqIndex = temp_decodingDataSimplified.seqIndex;
            trialIndex_bool_memoryCorrect = temp_decodingDataSimplified.trialIndex_bool_memoryCorrect;
            target_seqSet_inOne = temp_decodingDataSimplified.target_seqSet_inOne;
            
            temp_trialIndex_bool = trialIndex_bool_memoryCorrect; % use memory correct trial
            seqIndex_valid = seqIndex(temp_trialIndex_bool);
            
            temp_seqCount = zeros(1,sum(numSeq(1:3)));
            for tempi=1:sum(numSeq(1:3))
                temp_seqCount(tempi) = sum(seqIndex_valid==tempi,'all');
            end
            seqCount_valid_allSessions{currentSessionIndex} = temp_seqCount;
        end
        
    end
    
    
    seqCount_valid_allSessions_collapsed = zeros(length(seqCount_valid_allSessions{1}),num_FOV_AB);
    for tempi=1:num_FOV_AB
        seqCount_valid_allSessions_collapsed(:,tempi) = seqCount_valid_allSessions{tempi};
    end
    
    seqCount_valid_allSessions_collapsed_min = min(seqCount_valid_allSessions_collapsed,[],2);
    
    a = 1;
    
    seqIndex_valid_fromAllSessionsMin = zeros(sum(seqCount_valid_allSessions_collapsed_min),1);
    for tempi=1:length(seqCount_valid_allSessions_collapsed_min)
        temp_range = (sum(seqCount_valid_allSessions_collapsed_min(1:tempi-1))+1):sum(seqCount_valid_allSessions_collapsed_min(1:tempi));
        
        seqIndex_valid_fromAllSessionsMin(temp_range) = tempi;
    end
    
    boolIndex_location_seq = false(numFrames,length(target_seqSet_inOne));
    for tempi=1:length(target_seqSet_inOne)
        currentSequence = target_seqSet_inOne{tempi};
        boolIndex_location_seq(currentSequence,tempi) = true;
    end
    boolIndex_location_seq_T = boolIndex_location_seq';
    
    
    threshold_locRatio = 0.055;%0.06-->0.05-->0.045
    threshold_locStd = 0.026;%0.0278-->0.025-->0.026
    
    fun_locBasedResample_seqCount_Name_v = autoGetFunName_myScripts('fun_locBasedResample_seqCount', [targetPATH '\functions']);
    fun_locBasedResample_seqCount = str2func(fun_locBasedResample_seqCount_Name_v);
    
    options_locBased.threshold_locRatio = threshold_locRatio;
    options_locBased.threshold_locStd = threshold_locStd;
    options_locBased.numSeq = numSeq;
    options_locBased.numFrames = numFrames;
    options_locBased.seqIndex_valid = seqIndex_valid_fromAllSessionsMin;
    options_locBased.boolIndex_location_seq_T = boolIndex_location_seq_T;
    
    resampleTrialCount_seq = fun_locBasedResample_seqCount(options_locBased);    
end
resampleTrialCount_seq;