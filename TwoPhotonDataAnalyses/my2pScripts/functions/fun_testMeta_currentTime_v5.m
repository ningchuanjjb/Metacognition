function testMeta_output = fun_testMeta_currentTime_v5(F_dff_currentTime,options_testMeta)


% F_dff_currentTime = options_testMeta.F_dff_currentTime;
KFold_num = options_testMeta.KFold_num;
svm_Meta = options_testMeta.svm_Meta;
choiceBoolIndex = options_testMeta.choiceBoolIndex;
% seqIndex = options_testMeta.seqIndex;
noChoiceBoolIndex = options_testMeta.noChoiceBoolIndex;

% F_dff_currentTime = F_dff_decisionBin1;
%%  Test meta-memory on training set trials (choice trials)
% multi_Posterior_cell
multi_Posterior_cell = cell(1, KFold_num);
temp_Mdl_CV_binary = svm_Meta.temp_Mdl_CV_binary;
validSeqBoolIndex = svm_Meta.validSeqBoolIndex;
x_raw = F_dff_currentTime';
x = x_raw(choiceBoolIndex,:);
% F_dff_T = x;
F_dff_T = x(validSeqBoolIndex,:);


for tempk=1:KFold_num
    tempTrialBoolIndex_fold = temp_Mdl_CV_binary.ModelParameters.Generator.UseObsForIter(:,tempk);
    temp_F_dff_T = F_dff_T(~tempTrialBoolIndex_fold,:); %#ok<*PFBNS>
    temp_F_dff_T_2d = temp_F_dff_T;
    [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary.Trained{tempk},temp_F_dff_T_2d);
    if size(tempPosterior,2) == 1
        tempPosterior_2 = tempPosterior(:,1);
    else
        tempPosterior_2 = tempPosterior(:,2);
    end
    multi_Posterior_cell{tempk} = tempPosterior_2;
end

% Posterior_2d
% Posterior_2d = zeros(size(F_dff_T,1),1);
% Posterior_2d = nan(size(F_dff_T,1),1);
Posterior_2d = zeros(sum(validSeqBoolIndex),1);
for tempk=1:KFold_num
    temp_Posterior = multi_Posterior_cell{tempk};
    tempTrialBoolIndex_fold = temp_Mdl_CV_binary.ModelParameters.Generator.UseObsForIter(:,tempk);
    Posterior_2d(~tempTrialBoolIndex_fold) = temp_Posterior;
end
Posterior_2d_n11n = Posterior_2d;

% meta_trialLevel_choice_currentTime = Posterior_2d_n11n;
meta_trialLevel_choice_currentTime = nan(size(validSeqBoolIndex,2),1);
meta_trialLevel_choice_currentTime(validSeqBoolIndex) = Posterior_2d_n11n;

%%  Test meta-memory on other trials (noChoice trials)
x_raw = F_dff_currentTime';
x = x_raw(noChoiceBoolIndex,:);

temp_Mdl_CV_binary = svm_Meta.temp_Mdl_CV_binary;
temptemp_Posterior_2d_kfold = zeros(KFold_num,size(x,1));
for tempk=1:KFold_num
    [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary.Trained{tempk},x);% Very time-consuming!!!
    if size(tempPosterior,2) == 1
        tempPosterior_2 = tempPosterior(:,1);
    else
        tempPosterior_2 = tempPosterior(:,2);
    end
    temptemp_Posterior_2d_kfold(tempk,:) = tempPosterior_2;
end
temptemp_Posterior_2d_kfoldMean = squeeze(mean(temptemp_Posterior_2d_kfold,1));
temptemp_Posterior_2d_n11n = temptemp_Posterior_2d_kfoldMean';

meta_trialLevel_noChoice_currentTime = temptemp_Posterior_2d_n11n;

%% Test meta-memory on all trials
% meta_trialLevel_currentTime = nan(length(seqIndex),1);
meta_trialLevel_currentTime = nan(length(choiceBoolIndex),1);
meta_trialLevel_currentTime(choiceBoolIndex) = meta_trialLevel_choice_currentTime;
meta_trialLevel_currentTime(noChoiceBoolIndex) = meta_trialLevel_noChoice_currentTime;


%%
testMeta_output = struct;
testMeta_output.meta_trialLevel_choice_currentTime = meta_trialLevel_choice_currentTime;
testMeta_output.meta_trialLevel_noChoice_currentTime = meta_trialLevel_noChoice_currentTime;
testMeta_output.meta_trialLevel_currentTime = meta_trialLevel_currentTime;


