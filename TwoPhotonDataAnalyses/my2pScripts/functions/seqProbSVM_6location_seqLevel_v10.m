function svm_posterior_seqLevel = seqProbSVM_6location_seqLevel_v10(F_dff_raw,temp_svm_Y,t_SVM,numFrames) %#ok<*INUSD>

%% Initialization
F_dff = F_dff_raw; %#ok<*NASGU>

sum_F_dff = sum(F_dff,1);
validSeqBoolIndex = ~isnan(sum_F_dff);

temp_svm_Y_valid = temp_svm_Y(:,validSeqBoolIndex);

F_dff = F_dff(:,validSeqBoolIndex);

F_dff_T = F_dff';


a = 1;
%% KFold_num
KFold_num = zeros(1,numFrames);
for temploc=1:numFrames
    a = 1;
    currentLabel = temp_svm_Y_valid(temploc,:);
    
    uniqueSeq = unique(currentLabel);
    seqCount = zeros(1,length(uniqueSeq));
    for tempi=1:length(uniqueSeq)
        seqCount(tempi) = sum(currentLabel==uniqueSeq(tempi),'all');
    end
    temp_KFold_num = min(15,min(seqCount));%10-->8-->5
    temp_KFold_num = max(temp_KFold_num,3);%2
    
    KFold_num(temploc) = temp_KFold_num;
end
KFold_num_raw = KFold_num;
KFold_num = min(KFold_num_raw);

% KFold_num = 15;

temp_Mdl_CV_binary = cell(numFrames,1);
multi_Posterior_cell = cell(1, KFold_num);
Posterior_2d = zeros(sum(validSeqBoolIndex),numFrames);

%% temp_Mdl_CV_binary
c_all = cell(1,numFrames);
for temploc=1:numFrames
    a = 1;
    %currentLabel_raw = temp_svm_Y(temploc,:);
    %currentLabel = currentLabel_raw(validSeqBoolIndex);
    
    currentLabel = temp_svm_Y_valid(temploc,:);
    
    c = cvpartition(currentLabel,'KFold',KFold_num);
    c_all{temploc} = c;
    
    temp_svm_X = F_dff_T; %#ok<*PFIIN>
    temp_Mdl_binary = fitcecoc(temp_svm_X,currentLabel,'Learners',t_SVM,'FitPosterior',true); %#ok<*PFOUS>
    temp_Mdl_CV_binary{temploc} = crossval(temp_Mdl_binary,'CVPartition',c);
end


%% multi_Posterior_cell
for temploc=1:numFrames
    for tempk=1:KFold_num
        tempTrialBoolIndex_fold = temp_Mdl_CV_binary{temploc}.ModelParameters.Generator.UseObsForIter(:,tempk);
        temp_F_dff_T = F_dff_T(~tempTrialBoolIndex_fold,:); %#ok<*PFBNS>
        temp_F_dff_T_2d = temp_F_dff_T;
        [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary{temploc}.Trained{tempk},temp_F_dff_T_2d);
        tempPosterior_2 = tempPosterior(:,2);
        multi_Posterior_cell{tempk}(:,temploc) = tempPosterior_2;
    end
end


%% Posterior_2d
for tempk=1:KFold_num
    temp_Posterior = multi_Posterior_cell{tempk};
    for temploc=1:numFrames
        tempTrialBoolIndex_fold = temp_Mdl_CV_binary{temploc}.ModelParameters.Generator.UseObsForIter(:,tempk);
        Posterior_2d(~tempTrialBoolIndex_fold,temploc) = temp_Posterior(:,temploc);
    end
end


%% svm_posterior_seqLevel
svm_posterior_seqLevel = zeros(1,sum(validSeqBoolIndex));
a = 1;
for tempi=1:sum(validSeqBoolIndex)
    tempSeq_bool = temp_svm_Y_valid(:,tempi)';
    tempPosterior_locationLevel = Posterior_2d(tempi,:);
    
    temp1 = tempPosterior_locationLevel.*tempSeq_bool;
    temp2 = (1-tempPosterior_locationLevel).*~tempSeq_bool;
    %temp2 = 1.*~tempSeq_bool;
    temp3 = temp1 + temp2;
    temp_svm_posterior_seqLevel = prod(temp3);
        
    svm_posterior_seqLevel(tempi) = temp_svm_posterior_seqLevel;
    
    
%     temp_seqLength = sum(tempSeq_bool);    
%     temp_seqIndex = find(sum(temp_svm_Y_valid,1) == temp_seqLength);    
%     tempj_svm_posterior_seqLevel = zeros(1,length(temp_seqIndex));
%     temptempj = 0;
%     for tempj=temp_seqIndex(1):temp_seqIndex(end)
%         temptempj = temptempj + 1;
%         tempjSeq_bool = temp_svm_Y_valid(:,tempj)';
%         
%         tempj1 = tempPosterior_locationLevel.*tempjSeq_bool;
%         tempj2 = (1-tempPosterior_locationLevel).*~tempjSeq_bool;
%         tempj3 = tempj1 + tempj2;
%         tempj_svm_posterior_seqLevel(temptempj) = prod(tempj3);
%     end
%     temp_svm_posterior_seqLevel_n11n = temp_svm_posterior_seqLevel/sum(tempj_svm_posterior_seqLevel);        
%     svm_posterior_seqLevel(tempi) = temp_svm_posterior_seqLevel_n11n;
end

a = 1;

% prob_mean = mean(svm_posterior_seqLevel)
% prob_median = median(svm_posterior_seqLevel)
% a = 1;
% figure(2)
% plot(svm_posterior_seqLevel)


%% End
