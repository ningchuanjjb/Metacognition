function [svm_accuracy_seqLevel,validBoolIndex_seqLevel] = ...
    seqProbSVM_56seq_trialLevel_v10(F_dff_raw,temp_svm_Y_raw,boolIndex_location_seq,t_SVM,numFrames,F_dff_littleBin_raw) %#ok<*INUSD>

%% Initialization
F_dff = F_dff_raw; %#ok<*NASGU>
F_dff_littleBin = F_dff_littleBin_raw; %#ok<*NASGU>

if_resample = 1;
if_resample_normal0_littleBin1 = 0;%1

temp_svm_Y_raw2 = zeros(1,size(temp_svm_Y_raw,2));
for tempi=1:size(boolIndex_location_seq,2)
    tempSeq_bool = boolIndex_location_seq(:,tempi)';
    temp1_boolIndex = sum(temp_svm_Y_raw' == tempSeq_bool,2) == numFrames;
    temp_svm_Y_raw2(temp1_boolIndex) = tempi;
end

a = 1;

if if_resample == 1
    uniqueSeq = unique(temp_svm_Y_raw2);
    seqCount = zeros(1,length(uniqueSeq));
    for tempi=1:length(uniqueSeq)
        seqCount(tempi) = sum(temp_svm_Y_raw2==uniqueSeq(tempi),'all');
    end
    max_seqCount = max(seqCount);
    
    min_seqCount = min(seqCount);
    
    
    % resampleCount = round(max_seqCount*0.5);%2-->1
    if min_seqCount < 6
        resampleCount = 6;
    else
        resampleCount = min_seqCount;
    end
    %resampleCount = 10;
    
    [temp_svm_Y_sort,I] = sort(temp_svm_Y_raw2);
    F_dff_sort = F_dff(:,I);
    F_dff_littleBin_sort = F_dff_littleBin(:,I,:);
    if if_resample_normal0_littleBin1 == 0
        %% normal resample
        temp_svm_YIndex_resample_2d = zeros(length(seqCount),resampleCount);
        for tempi=1:length(seqCount)
            currentSeq = uniqueSeq(tempi);
            currentSeq_Y = find(temp_svm_Y_sort==currentSeq);
            currentSeq_trialCount = seqCount(tempi);
            
            tempIndex_resample = randsample(currentSeq_trialCount,resampleCount,true)';
            temp_svm_YIndex_resample_2d(tempi,:) = currentSeq_Y(tempIndex_resample);
        end
        temp_svm_YIndex_resample = reshape(temp_svm_YIndex_resample_2d',1,[]);
        
        temp_svm_Y_resample = temp_svm_Y_sort(temp_svm_YIndex_resample);
        F_dff_resample = F_dff_sort(:,temp_svm_YIndex_resample);        
        
    elseif if_resample_normal0_littleBin1 == 1
        %% resample of little time bin
        littleBinNum = size(F_dff_littleBin_sort,3);
        repeatNum = 1;
        
        temp_svm_Y_resample_2d = zeros(length(seqCount),resampleCount);
        for tempi=1:length(seqCount)
            temp_svm_Y_resample_2d(tempi,:) = uniqueSeq(tempi);
        end
        temp_svm_Y_resample = reshape(temp_svm_Y_resample_2d',1,[]);
                
        F_dff_resample = zeros(size(F_dff,1),length(seqCount)*resampleCount);
        for tempi=1:length(seqCount)
            currentSeq = uniqueSeq(tempi);
            currentSeq_Y = find(temp_svm_Y_sort==currentSeq);
            currentSeq_trialCount = seqCount(tempi);
            
            tempIndex_resample = zeros(littleBinNum,repeatNum,resampleCount);
            for tempj=1:littleBinNum
                for tempk=1:repeatNum
                    tempIndex_resample(tempj,tempk,:) = randsample(currentSeq_trialCount,resampleCount,true)';
                end
            end
            temp_YIndex_resample = currentSeq_Y(tempIndex_resample);
            
            temp_dff_littleBin = zeros(size(F_dff_littleBin_sort,1),littleBinNum,resampleCount);
            for tempj=1:littleBinNum
                temptemp_YIndex = squeeze(temp_YIndex_resample(tempj,:,:));
                if repeatNum == 1
                    temptemp_YIndex = temptemp_YIndex';
                end
                for tempk=1:resampleCount
                    temptemp_dff = F_dff_littleBin_sort(:,temptemp_YIndex(:,tempk),tempj);
                    temp_dff_littleBin(:,tempj,tempk) = mean(temptemp_dff,2);
                end
            end
            temp_dff_littleBinMean = squeeze(mean(temp_dff_littleBin,2));
            
            tempIndex = find(temp_svm_Y_resample == uniqueSeq(tempi));
            F_dff_resample(:,tempIndex) = temp_dff_littleBinMean;    %#ok<*FNDSB>
            
            a = 1;
        end
                
    end
            
    temp_svm_Y = temp_svm_Y_resample;
    F_dff = F_dff_resample;
    
elseif if_resample == 0    
    temp_svm_Y = temp_svm_Y_raw2;
    F_dff = F_dff_raw;
    
end

temp_svm_Y_T = temp_svm_Y';
F_dff_T = F_dff';


a = 1;
%% KFold_num
% t_KFold_num = tic;
currentLabel = temp_svm_Y;
uniqueSeq = unique(currentLabel);
seqCount = zeros(1,length(uniqueSeq));
for tempi=1:length(uniqueSeq)
    seqCount(tempi) = sum(currentLabel==uniqueSeq(tempi),'all');
end
KFold_num = min(10,min(seqCount));%10-->8-->5-->15-->10
KFold_num = max(KFold_num,3);%2-->3
% fprintf('t_KFold_num = %.1f secs.\n',toc(t_KFold_num));


a = 1;

%% temp_Mdl_CV
a = 1;
currentLabel = temp_svm_Y;
c = cvpartition(currentLabel,'KFold',KFold_num);

temp_svm_X = F_dff_T; %#ok<*PFIIN>
temp_Mdl = fitcecoc(temp_svm_X,currentLabel,'Learners',t_SVM); %#ok<*PFOUS>
temp_Mdl_CV = crossval(temp_Mdl,'CVPartition',c); % Very time-consuming!!!

svm_predictLabel = kfoldPredict(temp_Mdl_CV);

a = 1;

%% svm_accuracy_seqLevel
% t_svm_accuracy_seqLevel = tic;
svm_accuracy_seqLevel = zeros(1,size(boolIndex_location_seq,2));
validBoolIndex_seqLevel = true(1,size(boolIndex_location_seq,2));
a = 1;
for tempi=1:size(boolIndex_location_seq,2)
    temp1_boolIndex = temp_svm_Y_T == tempi;
    
    if sum(temp1_boolIndex) == 0
        validBoolIndex_seqLevel(tempi) = false;
        continue
    end
    
    temp_svm_accuracy_seqLevel = sum(svm_predictLabel(temp1_boolIndex)==tempi)/sum(temp1_boolIndex);
    svm_accuracy_seqLevel(tempi) = temp_svm_accuracy_seqLevel;
end
% fprintf('t_svm_accuracy_seqLevel = %.1f secs.\n',toc(t_svm_accuracy_seqLevel));
a = 1;

% prob_mean = mean(svm_accuracy_seqLevel)
% prob_median = median(svm_accuracy_seqLevel)
% a = 1;
% figure(2)
% plot(svm_accuracy_seqLevel)

%% End
