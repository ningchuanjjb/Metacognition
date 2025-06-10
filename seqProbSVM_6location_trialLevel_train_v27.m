function svm_outs = ...
    seqProbSVM_6location_trialLevel_train_v27(F_dff_raw,temp_svm_Y,svm_options,svm_options2) %#ok<*INUSD>

%% Initialization

boolIndex_location_seq = svm_options.boolIndex_location_seq;
t_SVM = svm_options.t_decoder;
numFrames = svm_options.numFrames;
targetPATH = svm_options.targetPATH;
gAcc = svm_options.gAcc;
seq_range = svm_options.seq_range;
target_length = svm_options.target_length;
if_normal0_resample1_leaveOneSeqOut2 = svm_options.if_normal0_resample1_leaveOneSeqOut2;
if_decodingAcc_threshold0_sort1 = svm_options.if_decodingAcc_threshold0_sort1;

F_dff = F_dff_raw; %#ok<*NASGU>

sum_F_dff = sum(F_dff,1);
validSeqBoolIndex = ~isnan(sum_F_dff);

temp_svm_Y_valid = temp_svm_Y(:,validSeqBoolIndex);
F_dff = F_dff(:,validSeqBoolIndex);

temp_svm_Y_valid_T = temp_svm_Y_valid';
F_dff_T = F_dff';

num_roi = size(F_dff,1);


%% Some preparation
boolIndex_location_seq_T = boolIndex_location_seq';

seqIndex_valid = zeros(1,size(temp_svm_Y_valid_T,1));
for tempi=1:size(temp_svm_Y_valid_T,1)
    temp_seq = temp_svm_Y_valid_T(tempi,:);
    for tempj=1:size(boolIndex_location_seq_T,1)
        if boolIndex_location_seq_T(tempj,:) == temp_seq
            break
        end
    end
    seqIndex_valid(tempi) = tempj;
end

temp_uniqueSeq = seq_range;
temp_seqCount = zeros(1,length(seq_range));
for tempi=1:length(seq_range)
    tempj = seq_range(tempi);
    temp_seqCount(tempi) = sum(seqIndex_valid==tempj,'all');
end


temp_locCount = zeros(1,numFrames);
for tempi=1:length(seq_range)
    temptempIndex = seq_range(tempi);
    temptempSeqBoolIndex = boolIndex_location_seq_T(temptempIndex,:);
    temp_locCount(temptempSeqBoolIndex) = temp_locCount(temptempSeqBoolIndex)+temp_seqCount(tempi);
end

a = 1;

%% To setup kfold params
% if exist('svm_options.c','var') == 1
if if_normal0_resample1_leaveOneSeqOut2 == 1
    % if it's resample mode, to control rng seed stable in parfor, cvpartition has been set already, before parfor loop
    c = svm_options2.c;
    KFold_num = svm_options2.KFold_num;
elseif if_normal0_resample1_leaveOneSeqOut2 == 0
    KFold_num = 10;%10-->5-->4
    % this would lead to same CV partition in each location
    c = cvpartition(size(temp_svm_Y_valid,2),'KFold',KFold_num);    
elseif if_normal0_resample1_leaveOneSeqOut2 == 2
    a = 1;
end

%% Weight
if_weight_seq0_loc1 = 1;
if if_weight_seq0_loc1 == 0
    if isfield(svm_options,'if_weighted') == 1
        if svm_options.if_weighted == 1
            
            temp_seqCount_weight = 1./temp_seqCount;
            
            W = ones(1,size(temp_svm_Y_valid_T,1));
            for tempi=1:length(seq_range)
                tempj = seq_range(tempi);
                temptempBoolIndex = seqIndex_valid==tempj;
                
                W(temptempBoolIndex) = temp_seqCount_weight(tempi);
            end
            %W = W.^2;
            W = W./sum(W);
        elseif svm_options.if_weighted == 0
            W = ones(1,size(temp_svm_Y_valid_T,1));
        end
    else
        W = ones(1,size(temp_svm_Y_valid_T,1));
    end
    
elseif if_weight_seq0_loc1 == 1    
    if isfield(svm_options,'if_weighted') == 1
        if svm_options.if_weighted == 1
            
            temp_locCount_weight = 1./temp_locCount;
            
            W = ones(numFrames,size(temp_svm_Y_valid_T,1));
            for tempi=1:numFrames
                temptempBoolIndex = temp_svm_Y_valid_T(:,tempi) == true;
                
                tempw1 = 1./sum(temptempBoolIndex==true);
                tempw2 = 1./sum(temptempBoolIndex==false);
                
                W(tempi,temptempBoolIndex) = tempw1;
                W(tempi,~temptempBoolIndex) = tempw2;
            end
            W = W./sum(W,2);
            a1 = sum(W,2);

        elseif svm_options.if_weighted == 0
            W = ones(numFrames,size(temp_svm_Y_valid_T,1));
        end
    else
        W = ones(numFrames,size(temp_svm_Y_valid_T,1));
    end
end

%% temp_Mdl_CV_binary
temp_Mdl_binary = cell(numFrames,1);
temp_Mdl_CV_binary = cell(numFrames,1);
% warning('off');
for tempi=1:numFrames
    temp_svm_X = F_dff_T;
    currentLabel = temp_svm_Y_valid(tempi,:);
    
    % this would lead to different CV partition in each location
    %c = cvpartition(size(temp_svm_Y_valid,2),'KFold',KFold_num);
    %c = cvpartition(currentLabel,'KFold',KFold_num);
    
    %temp_Mdl_binary{tempi} = fitcecoc(temp_svm_X,currentLabel,'Learners',t_SVM,'FitPosterior',true); %#ok<*PFOUS>
    if if_weight_seq0_loc1 == 0
        tempW = W;
    elseif if_weight_seq0_loc1 == 1
        tempW = W(tempi,:);
    end
    %if tempi == 4
    if tempi == 3
        a = 1;
    end
    
    temp_Mdl_binary{tempi} = fitcecoc(temp_svm_X,currentLabel,'Learners',t_SVM,'FitPosterior',true,'Weights',tempW); %#ok<*PFOUS>
    temp_Mdl_CV_binary{tempi} = crossval(temp_Mdl_binary{tempi},'CVPartition',c); % Very time-consuming!!!
end
% warning('on');

%% multi_Posterior_cell
multi_Posterior_cell = cell(1, KFold_num);
for temploc=1:numFrames
    for tempk=1:KFold_num
        tempTrialBoolIndex_fold = temp_Mdl_CV_binary{temploc}.ModelParameters.Generator.UseObsForIter(:,tempk);
        temp_F_dff_T = F_dff_T(~tempTrialBoolIndex_fold,:); %#ok<*PFBNS>
        temp_F_dff_T_2d = temp_F_dff_T;
        [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary{temploc}.Trained{tempk},temp_F_dff_T_2d);
        %tempPosterior_2 = tempPosterior(:,2);
        if size(tempPosterior,2) == 1
            tempPosterior_2 = tempPosterior(:,1);
        else
            tempPosterior_2 = tempPosterior(:,2);
        end
        multi_Posterior_cell{tempk}(:,temploc) = tempPosterior_2;
    end
end
% find(~tempTrialBoolIndex_fold==1)'
a = 1;
%% Posterior_2d
Posterior_2d = zeros(size(temp_svm_Y_valid,2),numFrames);

for tempk=1:KFold_num
    temp_Posterior = multi_Posterior_cell{tempk};
    for temploc=1:numFrames
        tempTrialBoolIndex_fold = temp_Mdl_CV_binary{temploc}.ModelParameters.Generator.UseObsForIter(:,tempk);
        Posterior_2d(~tempTrialBoolIndex_fold,temploc) = temp_Posterior(:,temploc);
    end
end

if_n11n = 0;%0

if if_n11n == 1
    Posterior_2d_n11n = (Posterior_2d./sum(Posterior_2d,2))*target_length;
    a = 0;
    while true
        tempindex = Posterior_2d_n11n>1;
        if sum(tempindex,'all') == 0
            break
        end
        Posterior_2d_n11n(Posterior_2d_n11n>1) = 1;
        Posterior_2d_n11n = (Posterior_2d_n11n./sum(Posterior_2d_n11n,2))*target_length;
        a = a + 1;
    end
    a;
    
    p0_full = [0 0 0 0];
    %p0_full = [0.0639    0.0649    0.0878    0.2773];%kwp0_lp0_v5_loopN11N_20230904A
    %p0_full = [0.0447    0.0494    0.0740    0.9998];%kwp0_lp0_v5_loopN11N_onlyDiagonal_20230904B
    p0 = p0_full(target_length);
    
    Posterior_2d_n11n = Posterior_2d_n11n + p0;
    Posterior_2d_n11n = (Posterior_2d_n11n./sum(Posterior_2d_n11n,2))*target_length;
    a = 0;
    while true
        tempindex = Posterior_2d_n11n>1;
        if sum(tempindex,'all') == 0
            break
        end
        Posterior_2d_n11n(Posterior_2d_n11n>1) = 1;
        Posterior_2d_n11n = (Posterior_2d_n11n./sum(Posterior_2d_n11n,2))*target_length;
        a = a + 1;
    end
    a;
elseif if_n11n == 0
    Posterior_2d_n11n = Posterior_2d;
end


% gAcc = [0.9,0.9,0.9,0.9,0.75,0.8];


fun_tempModel_Name_v = autoGetFunName_myScripts('fun_tempModel', [targetPATH '\functions']);
fun_tempModel = str2func(fun_tempModel_Name_v);

tempModel_options.Posterior_2d_n11n = Posterior_2d_n11n;
tempModel_options.temp_svm_Y_valid_T = temp_svm_Y_valid_T;
tempModel_options.numFrames = numFrames;
tempModel_options.boolIndex_location_seq = boolIndex_location_seq;
tempModel_options.gAcc = gAcc;
tempModel_options.seq_range = seq_range;
tempModel_options.target_length = target_length;
tempModel_options.if_n11n = if_n11n;
tempModel_options.if_decodingAcc_threshold0_sort1 = if_decodingAcc_threshold0_sort1;

w_initial = ones(1,numFrames);
coeff_power_initial = 1;
coeff_w_powered_initial = w_initial.^coeff_power_initial;
[Loss_initial,q_initial,Posterior_2d_w_initial] = fun_tempModel(coeff_w_powered_initial,tempModel_options);
a = 1;
% w = gAcc./q_initial;
w = w_initial;

if_w = 0;
loopNum = 2;%1

if if_w == 1
    temp_boolIndex_location_seq = boolIndex_location_seq(:,seq_range)';
    for tempi=1:numFrames
        tempBoolIndex = temp_boolIndex_location_seq(:,tempi) == 1;
        isnanBoolIndex = (isnan(gAcc) | isnan(q_initial))';
        tempBoolIndex2 = tempBoolIndex & ~isnanBoolIndex;
        w(tempi) = sum(gAcc(tempBoolIndex2)) / sum(q_initial(tempBoolIndex2));
    end
elseif if_w == 0
    loopNum = 0;
end

if loopNum > 0
    for tempLoopCount = 1:loopNum
        if tempLoopCount > 1
            a = 1;
            %w = w.*(gAcc./svm_posterior_lengthx);
            
            temp_w = ones(1,numFrames);
            for tempi=1:numFrames
                tempBoolIndex = temp_boolIndex_location_seq(:,tempi) == 1;
                isnanBoolIndex = (isnan(gAcc) | isnan(svm_posterior_lengthx))';
                tempBoolIndex2 = tempBoolIndex & ~isnanBoolIndex;
                temp_w(tempi) = sum(gAcc(tempBoolIndex2)) / sum(svm_posterior_lengthx(tempBoolIndex2));
            end
            w = w.*temp_w;
            
        end
        temp_range_coeff_power = 0.1:0.001:10;
        Loss = zeros(length(temp_range_coeff_power),1);
        q = zeros(length(temp_range_coeff_power),length(seq_range));
        for tempi=1:length(temp_range_coeff_power)
            temp_coeff_power = temp_range_coeff_power(tempi);
            temp_coeff_w_powered = w.^temp_coeff_power;
            [Loss(tempi),q(tempi,:)] = fun_tempModel(temp_coeff_w_powered,tempModel_options);
        end
        [M,I] = min(Loss); %#ok<*ASGLU>
        q_optimal = q(I,:);
        coeff_power = temp_range_coeff_power(I);
        coeff_w_powered = w.^coeff_power;
        
        [Loss_optimal,svm_posterior_lengthx,Posterior_2d_w] = fun_tempModel(coeff_w_powered,tempModel_options);
    end
    gAcc; %#ok<*VUNUS>
    Loss_optimal;
    coeff_w_powered;
elseif loopNum == 0
    svm_posterior_lengthx = q_initial;
    coeff_w_powered = coeff_w_powered_initial;
    Posterior_2d_w = Posterior_2d_w_initial;
end



% Loss_optimal

svm_outs = struct;
svm_outs.svm_posterior_lengthx = svm_posterior_lengthx;
svm_outs.Posterior_2d_n11n = Posterior_2d_n11n;
svm_outs.Posterior_2d_w = Posterior_2d_w;
svm_outs.temp_svm_Y_valid_T = temp_svm_Y_valid_T;
svm_outs.coeff_w_power = coeff_w_powered;
svm_outs.boolIndex_location_seq = boolIndex_location_seq;
svm_outs.seq_range = seq_range;
svm_outs.temp_Mdl_CV_binary = temp_Mdl_CV_binary;
svm_outs.temp_Mdl_binary = temp_Mdl_binary;

%% End
