function svm_outs = metaDecoder_v10(x,y,svm_options) %#ok<*INUSD>

%% Initialization
F_dff = x'; %#ok<*NASGU>
temp_svm_Y = y';
numSeq = svm_options.numSeq;
seqIndex_choice = svm_options.seqIndex_choice;
KFold_num = svm_options.KFold_num;
t_SVM = svm_options.t_decoder;
% metaDecoderPrctileThreshold = svm_options.metaDecoderPrctileThreshold;
if_resample = svm_options.if_resample;

sum_F_dff = sum(F_dff,1);
validSeqBoolIndex = ~isnan(sum_F_dff);

temp_svm_Y_valid = temp_svm_Y(:,validSeqBoolIndex);
F_dff = F_dff(:,validSeqBoolIndex);

temp_svm_Y_valid_T = temp_svm_Y_valid';
F_dff_T = F_dff';

num_roi = size(F_dff,1);

W = ones(1,size(temp_svm_Y_valid_T,1));
temp1 = sum(temp_svm_Y_valid_T)./length(temp_svm_Y_valid_T);
temp2 = 1-temp1;
% W(temp_svm_Y_valid_T) = temp2;
% W(~temp_svm_Y_valid_T) = temp1;

% temp_n = 7;%4(ROC771),5(ROC805),6(ROC809),7(ROC815),8(ROC816)
% temp_n = 7;%4(ROC),5(ROC787),6(ROC803),7(ROC818),8(ROC810)
temp_n = 7;%4(ROC),5(ROC799),6(ROC811),7(ROC819),8(ROC811)
W(temp_svm_Y_valid_T) = temp2^temp_n;
W(~temp_svm_Y_valid_T) = temp1^temp_n;
% W(temp_svm_Y_valid_T) = 16;
% W(~temp_svm_Y_valid_T) = 1;


%% To setup kfold params
c = cvpartition(size(temp_svm_Y_valid,2),'KFold',KFold_num);

%% temp_Mdl_CV_binary
    temp_svm_X = F_dff_T;
    currentLabel = temp_svm_Y_valid;        
    %temp_Mdl_binary= fitcecoc(temp_svm_X,currentLabel,'Learners',t_SVM,'FitPosterior',true); %#ok<*PFOUS>
    %temp_Mdl_binary= fitcecoc(temp_svm_X,currentLabel,'Learners',t_SVM,'FitPosterior',true,'Cost',[0 1;1 0]); %#ok<*PFOUS>    
    %temp_Mdl_binary= fitcecoc(temp_svm_X,currentLabel,'Learners',t_SVM,'FitPosterior',true,'Weights',W,'Cost',[0 2;1 0]); %#ok<*PFOUS>    
    
%     temp_Mdl_binary= fitcecoc(temp_svm_X,currentLabel,'Learners',t_SVM,'FitPosterior',true); %#ok<*PFOUS>
%     temp_Mdl_binary= fitcecoc(temp_svm_X,currentLabel,'Learners',t_SVM,'FitPosterior',true,'Weights',W); %#ok<*PFOUS>    

    %temp_cost = 10;%10(ROC815),15(ROC816)
    %temp_cost = 15;%10(ROC818),15(ROC822),20(ROC812)
    temp_cost = 14;%10(ROC816),13(ROC818),14(ROC821),15(ROC819),16(ROC816),20(ROC815)
    
    if if_resample == 0
        temp_Mdl_binary= fitcecoc(temp_svm_X,currentLabel,'Learners',t_SVM,'FitPosterior',true,'Weights',W,'Cost',[0 temp_cost;1 0]); %#ok<*PFOUS>
    elseif if_resample == 1
        temp_Mdl_binary= fitcecoc(temp_svm_X,currentLabel,'Learners',t_SVM,'FitPosterior',true);
    end
    temp_Mdl_CV_binary = crossval(temp_Mdl_binary,'CVPartition',c); % Very time-consuming!!!
    
%% multi_Posterior_cell
multi_Posterior_cell = cell(1, KFold_num);
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

%% Posterior_2d
Posterior_2d = zeros(size(temp_svm_Y_valid,2),1);
for tempk=1:KFold_num
    temp_Posterior = multi_Posterior_cell{tempk};
    tempTrialBoolIndex_fold = temp_Mdl_CV_binary.ModelParameters.Generator.UseObsForIter(:,tempk);
    Posterior_2d(~tempTrialBoolIndex_fold) = temp_Posterior;
end
Posterior_2d_n11n = Posterior_2d;


temp_Posterior_2d = Posterior_2d_n11n;

a = 1;

predictLabel_boolIndex = temp_Posterior_2d > 0.5;%0.5
% predictLabel_boolIndex = temp_Posterior_2d > 0.65;%0.5


% temp_threshold = prctile(temp_Posterior_2d,metaDecoderPrctileThreshold*100);
% a1 = sum(temp_Posterior_2d > temp_threshold);
% 
% predictLabel_boolIndex = temp_Posterior_2d > temp_threshold;%0.5


predict_boolIndex = predictLabel_boolIndex == temp_svm_Y_valid_T;

%% svm_posterior_seqLevel
svm_posterior_seqLevel = zeros(1,sum(numSeq(1:3)));
for tempi=1:size(svm_posterior_seqLevel,2)
    temptempBoolIndex_all = seqIndex_choice(validSeqBoolIndex) == tempi;
    temptempBoolIndex_hit = predict_boolIndex(temptempBoolIndex_all);
    svm_posterior_seqLevel(tempi) = sum(temptempBoolIndex_hit)./length(temptempBoolIndex_hit);
end

svm_posterior_allTrial = sum(predict_boolIndex)./length(predict_boolIndex);

temptempBoolIndex = seqIndex_choice(validSeqBoolIndex)<=sum(numSeq(1:3));
predict_boolIndex_length123 = predict_boolIndex(temptempBoolIndex);
svm_posterior_allTrial_length123 = sum(predict_boolIndex_length123)./length(predict_boolIndex_length123);


svm_outs = struct;
svm_outs.svm_posterior_seqLevel = svm_posterior_seqLevel;
svm_outs.svm_posterior_allTrial = svm_posterior_allTrial;
svm_outs.svm_posterior_allTrial_length123 = svm_posterior_allTrial_length123;
svm_outs.Posterior_2d_n11n = Posterior_2d_n11n;
svm_outs.temp_svm_Y_valid_T = temp_svm_Y_valid_T;
svm_outs.seqIndex_choice = seqIndex_choice;
svm_outs.temp_Mdl_CV_binary = temp_Mdl_CV_binary;
svm_outs.temp_Mdl_binary = temp_Mdl_binary;
svm_outs.validSeqBoolIndex = validSeqBoolIndex;

