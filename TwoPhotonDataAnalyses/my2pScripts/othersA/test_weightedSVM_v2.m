Posterior_2d_n11n_mean;

boolIndex_location_seq;
boolIndex_location_seq_T = boolIndex_location_seq';

locDistri_length = zeros(3,numFrames);
for tempi=1:3
    temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));    
    temp_seqBool_range = boolIndex_location_seq_T(temp_range,:);    
    temp_Posterior_2d_n11n_mean = Posterior_2d_n11n_mean(temp_range,:);        
    for tempj=1:numFrames
        temptempBoolIndex = temp_seqBool_range(:,tempj);        
        locDistri_length(tempi,tempj) = mean(temp_Posterior_2d_n11n_mean(temptempBoolIndex,tempj),'omitnan');
    end
end


locCount_length = zeros(3,numFrames);
seqCount_length = cell(3,1);
for tempi=1:3
    temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));    
    seqCount_length{tempi} = zeros(1,length(temp_range));
    for tempj=1:length(temp_range)
        temptempIndex = temp_range(tempj);
        temptempSeqBoolIndex = boolIndex_location_seq_T(temptempIndex,:);
        
        seqCount_length{tempi}(tempj) = sum(seqIndex_valid==temp_range(tempj),'all');        
        locCount_length(tempi,temptempSeqBoolIndex) = locCount_length(tempi,temptempSeqBoolIndex)+seqCount_length{tempi}(tempj);        
    end
end

a = 1;
seqCount_length;
locCount_length;
locDistri_length;

temp_r = zeros(1,3);
temp_p = zeros(1,3);
for tempi=1:3
    [temp_r(tempi),temp_p(tempi)] = corr(locCount_length(tempi,:)',locDistri_length(tempi,:)');
end

% temp_r
temp_p
% [temp_r_length123,temp_p_length123] = corr(reshape(locCount_length',[],1),reshape(locDistri_length',[],1));
a = 1;

% seqCount_length_inOne = zeros(1,sum(numSeq(1:3)));
seqCount_length_inOne = [];
for tempi=1:3
    seqCount_length_inOne = [seqCount_length_inOne; seqCount_length{tempi}']; %#ok<*AGROW>
end


temp2_r = zeros(1,numFrames);
temp2_p = zeros(1,numFrames);
for tempi=1:numFrames
    
    temptempBoolIndex1 = boolIndex_location_seq_T(1:sum(numSeq(1:3)),tempi);
       
    temptempBoolIndex2 = seqCount_length_inOne > 0;
    
    temptempBoolIndex = temptempBoolIndex1 & temptempBoolIndex2;
    
    [temp2_r(tempi),temp2_p(tempi)] = corr(seqCount_length_inOne(temptempBoolIndex),Posterior_2d_n11n_mean(temptempBoolIndex,tempi));
end
a = 1;






%% Correlation of r_n11n and trial count
r_n11n = zeros(1,sum(numSeq(1:valid_length)));
p_n11n = zeros(1,sum(numSeq(1:valid_length)));
for tempi=1:sum(numSeq(1:valid_length))
    [r_n11n(tempi),p_n11n(tempi)] = corr(Posterior_2d_n11n_mean(tempi,:)',Posterior_2d_model(tempi,:)');
end

r_n11n;
seqTrialCount';

[temptemp_r,temptemp_p] = corr(r_n11n(~isnan(r_n11n))',seqTrialCount(~isnan(r_n11n)));


temp_range = 1:6;
temp_range = 7:21;
temp_range = 22:41;

temp_r_n11n = r_n11n(temp_range);
temp_seqTrialCount = seqTrialCount(temp_range);

[temptemp_r,temptemp_p] = corr(temp_r_n11n(~isnan(temp_r_n11n))',temp_seqTrialCount(~isnan(temp_r_n11n)));