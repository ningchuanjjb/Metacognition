function [svm_accuracy_CV, crossTime_svm_accuracy] = ...
    svm_CCGP_eachLoc_v15(F_dff_lengthx_period1,temp_svm_Y,t_SVM,T_lengthx_period1,numFrames,responseBool_lengthx,sequenceBool_lengthx) %#ok<*INUSD>

seqLength = size(temp_svm_Y,1);

if seqLength > 1
    temp_svm_Y_merge = zeros(1,size(temp_svm_Y,2));
    for tempi=1:length(temp_svm_Y_merge)
        temp_str = [];
        for tempj=1:seqLength
            temp_str = [temp_str, num2str(temp_svm_Y(tempj,tempi))]; %#ok<*AGROW>
        end
        temp_svm_Y_merge(tempi) = str2double(temp_str);
    end
else
    temp_svm_Y_merge = temp_svm_Y;
end

uniqueSeq = unique(temp_svm_Y_merge);
seqCount = zeros(1,numFrames);
for tempi=1:length(uniqueSeq)
    seqCount(tempi) = sum(temp_svm_Y_merge==uniqueSeq(tempi),'all');    
end


KFold_num = min(10,min(seqCount));
KFold_num = max(KFold_num,3);%2
c = cvpartition(temp_svm_Y_merge,'KFold',KFold_num);


crossTime_svm_accuracy = zeros(T_lengthx_period1,T_lengthx_period1);
F_dff_lengthx_period1_T = permute(F_dff_lengthx_period1,[1 3 2]);

pause(2);% wait for 2 seconds, to avoid cpu overheat

parfor tempj=1:T_lengthx_period1
% for tempj=1:T_lengthx_period1
    
    %tempj = 45;%96

    temp_Mdl_CV_binary = cell(numFrames,1);
    for temploc=1:numFrames
        %temp_svm_Y_binary = (temp_svm_Y == temploc);
        temp_binary = ismember(temp_svm_Y,temploc);
        temp_svm_Y_binary = (sum(temp_binary,1)==1);
        
        temp_svm_X = F_dff_lengthx_period1(:,:,tempj); %#ok<*PFIIN>
        temp_Mdl_binary = fitcecoc(temp_svm_X,temp_svm_Y_binary,'Learners',t_SVM,'FitPosterior',true); %#ok<*PFOUS>
        temp_Mdl_CV_binary{temploc} = crossval(temp_Mdl_binary,'CVPartition',c);
    end

    multi_Posterior_cell = cell(1, KFold_num);
    for tempk=1:KFold_num
        tempTrialBoolIndex_fold = temp_Mdl_CV_binary{1}.ModelParameters.Generator.UseObsForIter(:,tempk);
        temp_F_dff_lengthx_period1_T = F_dff_lengthx_period1_T(~tempTrialBoolIndex_fold,:,:); %#ok<*PFBNS>
        temp_F_dff_lengthx_period1_T_2d = reshape(temp_F_dff_lengthx_period1_T,[],size(temp_F_dff_lengthx_period1_T,3));
        multi_Posterior_cell{tempk} = zeros(size(temp_F_dff_lengthx_period1_T_2d,1),numFrames);
        for temploc=1:numFrames
            [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary{temploc}.Trained{tempk},temp_F_dff_lengthx_period1_T_2d);
            tempPosterior_2 = tempPosterior(:,2);
            multi_Posterior_cell{tempk}(:,temploc) = tempPosterior_2;
        end
    end
    
    Posterior_3d = zeros(length(temp_svm_Y),T_lengthx_period1,numFrames);
    for tempk=1:KFold_num
        tempTrialBoolIndex_fold = temp_Mdl_CV_binary{1}.ModelParameters.Generator.UseObsForIter(:,tempk);
        temp_Posterior = multi_Posterior_cell{tempk};
        temp_Posterior_3d = reshape(temp_Posterior,[],T_lengthx_period1,numFrames);
        Posterior_3d(~tempTrialBoolIndex_fold,:,:) = temp_Posterior_3d;
    end
    %[~,temp_svm_label_2d] = max(Posterior_3d,[],3);   
    
    [~,I] = sort(Posterior_3d,3,'descend');
    I2 = sort(I(:,:,1:seqLength),3,'ascend');    
    I3 = reshape(I2,[],seqLength);
    I4 = cell(1,seqLength);
    for tempk=1:seqLength
        I4{tempk} = num2str(I3(:,tempk));        
    end
    I5 = [];
    for tempk=1:seqLength
        I5 = [I5,I4{tempk}];        
    end    
    I6 = reshape(str2num(I5),size(Posterior_3d,1),size(Posterior_3d,2)); %#ok<*NASGU,*ST2NM>
    temp_svm_label_2d = I6;
    
    
    crossTime_svm_accuracy(:,tempj) = sum(temp_svm_label_2d == temp_svm_Y_merge',1)/length(temp_svm_Y_merge);
    %crossTime_svm_accuracy(:,tempj) = sum(temp_svm_label_2d == temp_svm_Y',1)/length(temp_svm_Y);
end


temp_Bool = false(T_lengthx_period1,T_lengthx_period1);
for tempi=1:T_lengthx_period1
    temp_Bool(tempi,tempi) = true;
end
temp_acc = crossTime_svm_accuracy(temp_Bool)';
svm_accuracy_CV = temp_acc;

%% End
