function [svm_accuracy_CV, crossTime_svm_accuracy] = ...
    svm_CCGP_v4(F_dff_lengthx_period1,temp_svm_Y,t_SVM,T_lengthx_period1,numFrames)


uniqueSeq = unique(temp_svm_Y);
seqCount = zeros(1,numFrames);
for tempi=1:length(uniqueSeq)
    seqCount(tempi) = sum(temp_svm_Y==uniqueSeq(tempi),'all');    
end

KFold_num = min(10,min(seqCount));
KFold_num = max(KFold_num,3);%2
c = cvpartition(temp_svm_Y,'KFold',KFold_num);


svm_Mdl = cell(1,T_lengthx_period1);
svm_label = zeros(length(temp_svm_Y),T_lengthx_period1);
svm_accuracy = zeros(1,T_lengthx_period1);
svm_Mdl_CV = cell(1,T_lengthx_period1);
svm_label_CV = zeros(length(temp_svm_Y),T_lengthx_period1);
svm_accuracy_CV = zeros(1,T_lengthx_period1);
parfor tempi=1:T_lengthx_period1
% for tempi=1:T_lengthx_period1
    temp_svm_X = F_dff_lengthx_period1(:,:,tempi); %#ok<*PFIIN>    
    svm_Mdl{tempi} = fitcecoc(temp_svm_X,temp_svm_Y,'Learners',t_SVM); %#ok<*PFOUS> 
    svm_label(:,tempi) = predict(svm_Mdl{tempi},temp_svm_X);
    svm_accuracy(tempi) = sum(svm_label(:,tempi) == temp_svm_Y')/length(temp_svm_Y);    
    svm_Mdl_CV{tempi} = crossval(svm_Mdl{tempi},'CVPartition',c);
    svm_label_CV(:,tempi) = kfoldPredict(svm_Mdl_CV{tempi});
    svm_accuracy_CV(tempi) = sum(svm_label_CV(:,tempi) == temp_svm_Y')/length(temp_svm_Y);
    a = 1;
end

crossTime_svm_accuracy = zeros(T_lengthx_period1,T_lengthx_period1);
F_dff_lengthx_period1_T = permute(F_dff_lengthx_period1,[1 3 2]);
F_dff_lengthx_period1_T_2d = reshape(F_dff_lengthx_period1_T,[],size(F_dff_lengthx_period1_T,3)); %#ok<*NASGU>


parfor tempj=1:T_lengthx_period1
% for tempj=1:T_lengthx_period1
    multi_NegLoss_cell = cell(1, KFold_num);
    for tempk=1:KFold_num
        tempTrialBoolIndex_fold = svm_Mdl_CV{tempj}.Ensemble.ModelParams.Generator.UseObsForIter(:,tempk);
        temp_F_dff_lengthx_period1_T = F_dff_lengthx_period1_T(~tempTrialBoolIndex_fold,:,:); %#ok<*PFBNS>
        temp_F_dff_lengthx_period1_T_2d = reshape(temp_F_dff_lengthx_period1_T,[],size(temp_F_dff_lengthx_period1_T,3));
        [~,multi_NegLoss_cell{tempk}] = predict(svm_Mdl_CV{tempj}.Trained{tempk},temp_F_dff_lengthx_period1_T_2d);
        a = 1;
    end
    
    NegLoss_3d = zeros(length(temp_svm_Y),T_lengthx_period1,numFrames);
    for tempk=1:KFold_num
        tempTrialBoolIndex_fold = svm_Mdl_CV{tempj}.Ensemble.ModelParams.Generator.UseObsForIter(:,tempk);
        temp_NegLoss = multi_NegLoss_cell{tempk};
        temp_NegLoss_3d = reshape(temp_NegLoss,[],T_lengthx_period1,numFrames);
        NegLoss_3d(~tempTrialBoolIndex_fold,:,:) = temp_NegLoss_3d;
    end
    [~,temp_svm_label_2d] = max(NegLoss_3d,[],3);
    
    crossTime_svm_accuracy(:,tempj) = sum(temp_svm_label_2d == temp_svm_Y',1)/length(temp_svm_Y);
end
