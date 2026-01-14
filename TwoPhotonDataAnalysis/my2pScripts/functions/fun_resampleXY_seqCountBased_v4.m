function output = fun_resampleXY_seqCountBased_v4(options)

numSeq = options.numSeq;
resampleIterCount = options.resampleIterCount;
resampleTrialCount_seq = options.resampleTrialCount_seq;
seqIndex_valid = options.seqIndex_valid;
temp_F_dff_decisionBin = options.temp_F_dff_decisionBin;


%temp_F_dff_decisionBin_resample = zeros(resampleIterCount,size(temp_F_dff_decisionBin,1),sum(numSeq(1:3))*resampleTrialCount);
%seqIndex_valid_resample = zeros(resampleIterCount,sum(numSeq(1:3))*resampleTrialCount);
for tempIter=1:resampleIterCount
    
    tempIndex_resample_cell = cell(1,sum(numSeq(1:3)));
    for tempi=1:sum(numSeq(1:3))
        temp_resampleTrialCount = resampleTrialCount_seq(tempi);
        
        tempIndex_resample_cell{tempi} = zeros(1,temp_resampleTrialCount);
        tempIndex = find(seqIndex_valid==tempi);
        temp_length = length(tempIndex);
        if temp_length == 0
            continue
        end
        if temp_length >= temp_resampleTrialCount
            tempI = sort(randperm(temp_length,temp_resampleTrialCount));
        elseif temp_length < temp_resampleTrialCount
            temp_resampleTrialCount2 = temp_resampleTrialCount;
            temp_intLoopCount = floor(temp_resampleTrialCount2/temp_length);
            temp_resampleTrialCount2 = temp_resampleTrialCount2-temp_length*temp_intLoopCount;
            tempI1 = [];
            for temptempi=1:temp_intLoopCount
                tempI1 = [tempI1 (1:temp_length)]; %#ok<*AGROW>
            end
            tempI2 = sort(randperm(temp_length,temp_resampleTrialCount2));
            
            tempI = [tempI1 tempI2]; %#ok<*NASGU>
        end
        tempI = tempI(randperm(length(tempI))); % shuffle the resampled trial index
        
        tempIndex_resample_cell{tempi} = tempIndex(tempI);
    end
    
    tempIndex_resample_cell; %#ok<*VUNUS>
    tempIndex_resample_1d = [];
    for tempi=1:length(tempIndex_resample_cell)
        tempIndex_resample_1d = [tempIndex_resample_1d,tempIndex_resample_cell{tempi}];
    end
    a = 1;
    
    temptemp_isvalid = tempIndex_resample_1d > 0;
    
    % -1 means no trial in this sequence
    temptemp_F_dff_decisionBin_resample = -1*ones(size(temp_F_dff_decisionBin,1),length(tempIndex_resample_1d)); %#ok<*USENS>
    tempseqIndex_valid_resample = -1*ones(1,length(tempIndex_resample_1d));
    
    temptemp_F_dff_decisionBin_resample(:,temptemp_isvalid) = temp_F_dff_decisionBin(:,tempIndex_resample_1d(temptemp_isvalid));
    tempseqIndex_valid_resample(temptemp_isvalid) = seqIndex_valid(tempIndex_resample_1d(temptemp_isvalid));
    
    if tempIter == 1
        temp_F_dff_decisionBin_resample = zeros(resampleIterCount,size(temp_F_dff_decisionBin,1),length(tempIndex_resample_1d));
        seqIndex_valid_resample = zeros(resampleIterCount,length(tempIndex_resample_1d));
        temp_trialIndex_valid_resample = zeros(resampleIterCount,length(tempIndex_resample_1d));        
    end
    temp_F_dff_decisionBin_resample(tempIter,:,:) = temptemp_F_dff_decisionBin_resample;
    seqIndex_valid_resample(tempIter,:) = tempseqIndex_valid_resample;
    temp_trialIndex_valid_resample(tempIter,:) = tempIndex_resample_1d(temptemp_isvalid);
end
seqIndex_valid_resample_raw = seqIndex_valid_resample;
seqIndex_valid_resample = seqIndex_valid_resample(1,:);

output = struct;
output.temp_F_dff_decisionBin_resample = temp_F_dff_decisionBin_resample;
output.seqIndex_valid_resample = seqIndex_valid_resample;
output.temp_trialIndex_valid_resample = temp_trialIndex_valid_resample;