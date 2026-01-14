tempBoolIndex123_and;

temp_glm_peakLoc_length123;

temp_glm_peakBeta_length123;

temp_glm_std_length123;


a = 1;

roiNum_length123 = sum(tempBoolIndex123_and);

a = 1;

temp_dff_seq_length123 = F_dff_lengthx_delay1Bin_seq_multiFOV(tempBoolIndex123_and,:);
boolIndex_location_seq_T;
a = 1;

temp_dff_peakLocSeq_eachLength_length123 = cell(roiNum_length123,seq_length_max);
temp_dffMinIndex_peakLocSeq_eachLength_length123 = zeros(roiNum_length123,seq_length_max);
for tempi=1:roiNum_length123
    temptemp_peakLoc = temp_glm_peakLoc_length123(tempi,1);    
    for tempj=1:seq_length_max
        temp_range = (sum(numSeq(1:(tempj-1)))+1):sum(numSeq(1:tempj));
        temp_range_boolIndex = false(size(boolIndex_location_seq_T,1),1);
        temp_range_boolIndex(temp_range) = true;
        
        temptempBoolIndex = boolIndex_location_seq_T(:,temptemp_peakLoc) & temp_range_boolIndex;
        
        temptemp_dff = temp_dff_seq_length123(tempi,temptempBoolIndex);        
        temp_dff_peakLocSeq_eachLength_length123{tempi,tempj} = temptemp_dff;
        
        temptempIndex = find(temptempBoolIndex==true);
        temptempIndex2 = temptempIndex(temptemp_dff~=0);
        [M,I] = min(temptemp_dff(temptemp_dff~=0));
        tempMinSeqIndex = temptempIndex2(I);
        temp_dffMinIndex_peakLocSeq_eachLength_length123(tempi,tempj) = tempMinSeqIndex;
    end
end
temp_dff_peakLocSeq_eachLength_length123;
temp_dffMinIndex_peakLocSeq_eachLength_length123;
a = 1;

temp_dffSEM_peakLocSeq_eachLength_length123 = zeros(roiNum_length123,seq_length_max);
for tempi=1:roiNum_length123
    for tempj=1:seq_length_max
        temptemp_dff = temp_dff_peakLocSeq_eachLength_length123{tempi,tempj};
        temptemp_dff2 = temptemp_dff(temptemp_dff~=0);
        %temptemp_dffSEM = std(temptemp_dff2)/sqrt(length(temptemp_dff2));
        temptemp_dffSEM = std(temptemp_dff2);
        
        temp_dffSEM_peakLocSeq_eachLength_length123(tempi,tempj) = temptemp_dffSEM;
    end
end

% [~,p23_peakLocSeq_dffSEM,~,~] = ttest2(temp_dffSEM_peakLocSeq_eachLength_length123(:,2),temp_dffSEM_peakLocSeq_eachLength_length123(:,3),'Tail','right');
% [~,p24_peakLocSeq_dffSEM,~,~] = ttest2(temp_dffSEM_peakLocSeq_eachLength_length123(:,2),temp_dffSEM_peakLocSeq_eachLength_length123(:,4),'Tail','right');
% [~,p34_peakLocSeq_dffSEM,~,~] = ttest2(temp_dffSEM_peakLocSeq_eachLength_length123(:,3),temp_dffSEM_peakLocSeq_eachLength_length123(:,4),'Tail','right');
% 
% a2 = mean(temp_dffSEM_peakLocSeq_eachLength_length123(:,2));
% a3 = mean(temp_dffSEM_peakLocSeq_eachLength_length123(:,3));
% a4 = mean(temp_dffSEM_peakLocSeq_eachLength_length123(:,4));

a = 1;
% F_dff_lengthx_delay1Bin_seq_trial_multiFOV
temp_dff_seq_trial_length123 = F_dff_lengthx_delay1Bin_seq_trial_multiFOV(tempBoolIndex123_and,:);
temp_dff_trialAnovaP_peakLocSeq_eachLength_length123 = zeros(roiNum_length123,seq_length_max);
for tempi=1:roiNum_length123
    temptemp_peakLoc = temp_glm_peakLoc_length123(tempi,1);
    temptemp_dff_cell = temp_dff_seq_trial_length123(tempi,:);
    for tempj=1:seq_length_max
        temp_range = (sum(numSeq(1:(tempj-1)))+1):sum(numSeq(1:tempj));
        temp_range_boolIndex = false(size(boolIndex_location_seq_T,1),1);
        temp_range_boolIndex(temp_range) = true;
        
        temptempBoolIndex = boolIndex_location_seq_T(:,temptemp_peakLoc) & temp_range_boolIndex;        
        temptemp_dff_cell2 = temptemp_dff_cell(temptempBoolIndex);
        
        x = [];
        y = [];
        for tempk=1:length(temptemp_dff_cell2)
            if length(temptemp_dff_cell2{tempk}) <= 1
                continue
            end
            x = [x temptemp_dff_cell2{tempk}]; %#ok<*AGROW>
            y = [y tempk*ones(1,length(temptemp_dff_cell2{tempk}))];
        end
                
        [temptemp_p,~,~] = anova1(x, y,'off');
        %[temptemp_p,~,~] = anova1(x, y);
        temp_dff_trialAnovaP_peakLocSeq_eachLength_length123(tempi,tempj) = temptemp_p;
    end    
end
a = 1;
temp_dff_trialAnovaBool_peakLocSeq_eachLength_length123 = temp_dff_trialAnovaP_peakLocSeq_eachLength_length123 < 0.05;
a2 = sum(temp_dff_trialAnovaBool_peakLocSeq_eachLength_length123(:,2)); %#ok<*NASGU>
a3 = sum(temp_dff_trialAnovaBool_peakLocSeq_eachLength_length123(:,3));
a4 = sum(temp_dff_trialAnovaBool_peakLocSeq_eachLength_length123(:,4));
fprintf('Seq of max tuning location anova roi num = [%d, %d, %d].\n',a2,a3,a4);

% tempAnovaBoolIndex = temp_dff_trialAnovaBool_peakLocSeq_eachLength_length123(:,2) & temp_dff_trialAnovaBool_peakLocSeq_eachLength_length123(:,3);
tempAnovaBoolIndex = true(roiNum_length123,1);

% [~,p23_peakLocSeq_dffSEM,~,~] = ttest2(temp_dffSEM_peakLocSeq_eachLength_length123(:,2),temp_dffSEM_peakLocSeq_eachLength_length123(:,3),'Tail','right');
% [~,p24_peakLocSeq_dffSEM,~,~] = ttest2(temp_dffSEM_peakLocSeq_eachLength_length123(:,2),temp_dffSEM_peakLocSeq_eachLength_length123(:,4),'Tail','right');
% [~,p34_peakLocSeq_dffSEM,~,~] = ttest2(temp_dffSEM_peakLocSeq_eachLength_length123(:,3),temp_dffSEM_peakLocSeq_eachLength_length123(:,4),'Tail','right');
[~,p23_peakLocSeq_dffSEM,~,~] = ttest2(temp_dffSEM_peakLocSeq_eachLength_length123(tempAnovaBoolIndex,2),temp_dffSEM_peakLocSeq_eachLength_length123(tempAnovaBoolIndex,3),'Tail','right');
[~,p24_peakLocSeq_dffSEM,~,~] = ttest2(temp_dffSEM_peakLocSeq_eachLength_length123(tempAnovaBoolIndex,2),temp_dffSEM_peakLocSeq_eachLength_length123(tempAnovaBoolIndex,4),'Tail','right');
[~,p34_peakLocSeq_dffSEM,~,~] = ttest2(temp_dffSEM_peakLocSeq_eachLength_length123(tempAnovaBoolIndex,3),temp_dffSEM_peakLocSeq_eachLength_length123(tempAnovaBoolIndex,4),'Tail','right');

x = [];
y = [];
for tempi=2:4
    x = [x;tempi*ones(size(temp_dffSEM_peakLocSeq_eachLength_length123(tempAnovaBoolIndex,tempi),1),1)];    %#ok<*AGROW>
    y = [y;temp_dffSEM_peakLocSeq_eachLength_length123(tempAnovaBoolIndex,tempi)];
end
%cftool(x,y);
temp_mdl = fitglm(x,y,'linear');
p234_peakLocSeq_dffSEM = temp_mdl.Coefficients.pValue(2);

% a2 = mean(temp_dffSEM_peakLocSeq_eachLength_length123(:,2));
% a3 = mean(temp_dffSEM_peakLocSeq_eachLength_length123(:,3));
% a4 = mean(temp_dffSEM_peakLocSeq_eachLength_length123(:,4));
a2 = mean(temp_dffSEM_peakLocSeq_eachLength_length123(tempAnovaBoolIndex,2));
a3 = mean(temp_dffSEM_peakLocSeq_eachLength_length123(tempAnovaBoolIndex,3));
a4 = mean(temp_dffSEM_peakLocSeq_eachLength_length123(tempAnovaBoolIndex,4));

a = 1;

temp_dffMinIndex_peakLocSeq_eachLength_length123;

a = 1;

minPairLoc_peakLocSeq_eachLength_length123 = zeros(numFrames,numFrames,seq_length_max);
boolIndex_location_seq_T;
for tempi=1:roiNum_length123
    temptemp_dffMinIndex = temp_dffMinIndex_peakLocSeq_eachLength_length123(tempi,:);
    temptemp_peakLoc = temptemp_dffMinIndex(1);    
    for tempj=2:seq_length_max
        temp_seqBoolIndex = boolIndex_location_seq_T(temptemp_dffMinIndex(tempj),:);
        temp_seqBoolIndex(temptemp_peakLoc) = false;
        minPairLoc_peakLocSeq_eachLength_length123(temptemp_peakLoc,temp_seqBoolIndex,tempj) = ...
            minPairLoc_peakLocSeq_eachLength_length123(temptemp_peakLoc,temp_seqBoolIndex,tempj) + 1;                        
    end
end


