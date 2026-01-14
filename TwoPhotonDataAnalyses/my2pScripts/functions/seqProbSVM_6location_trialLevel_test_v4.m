function svm_outs = ...
    seqProbSVM_6location_trialLevel_test_v4(F_dff_raw,temp_svm_Y,svm_options) %#ok<*INUSD>

%% Initialization

if_resample = 0;


boolIndex_location_seq = svm_options.boolIndex_location_seq;
t_SVM = svm_options.t_decoder;
numFrames = svm_options.numFrames;
targetPATH = svm_options.targetPATH;
gAcc = svm_options.gAcc;
seq_range = svm_options.seq_range;
target_length = svm_options.target_length;

coeff_w_power = svm_options.coeff_w_power;
temp_Mdl_CV_binary = svm_options.temp_Mdl_CV_binary;
temp_Mdl_binary = svm_options.temp_Mdl_binary;


F_dff = F_dff_raw; %#ok<*NASGU>

sum_F_dff = sum(F_dff,1);
validSeqBoolIndex = ~isnan(sum_F_dff);

temp_svm_Y_valid = temp_svm_Y(:,validSeqBoolIndex);
F_dff = F_dff(:,validSeqBoolIndex);

temp_svm_Y_valid_T = temp_svm_Y_valid';
F_dff_T = F_dff';

num_roi = size(F_dff,1);


coeff_w_power;
temp_Mdl_CV_binary; %#ok<*VUNUS>
temp_Mdl_binary;
a = 1;

%% Posterior_2d
Posterior_2d = zeros(size(F_dff_T,1),numFrames);
for temploc=1:numFrames
    temp_svm_X = F_dff_T;
    currentLabel = temp_svm_Y_valid(temploc,:);
    
    [~,~,~,tempPosterior] = predict(temp_Mdl_binary{temploc},temp_svm_X);
    tempPosterior_2 = tempPosterior(:,2);
    Posterior_2d(:,temploc) = tempPosterior_2;         %#ok<*AGROW>
end   

Posterior_2d_raw = Posterior_2d;

a = 1;
% Posterior_2d_n11n = Posterior_2d_raw./sum(Posterior_2d_raw,2);
Posterior_2d_n11n = (Posterior_2d_raw./sum(Posterior_2d_raw,2))*target_length;

a = 0;
while true
    tempindex = Posterior_2d_n11n>1;
    if sum(tempindex,'all') == 0
        break
    end
    Posterior_2d_n11n(Posterior_2d_n11n>1) = 1;
    %break
    Posterior_2d_n11n = (Posterior_2d_n11n./sum(Posterior_2d_n11n,2))*target_length; %#ok<*UNRCH>
    a = a + 1;
    %break
end
a;



a = 1;
coeff_w_power;

fun_tempModel_Name_v = autoGetFunName_myScripts('fun_tempModel', [targetPATH '\functions']);
fun_tempModel = str2func(fun_tempModel_Name_v);

tempModel_options.Posterior_2d_n11n = Posterior_2d_n11n;
tempModel_options.temp_svm_Y_valid_T = temp_svm_Y_valid_T;
tempModel_options.numFrames = numFrames;
tempModel_options.boolIndex_location_seq = boolIndex_location_seq;
tempModel_options.gAcc = gAcc;
tempModel_options.seq_range = seq_range;
tempModel_options.target_length = target_length;

[Loss,svm_posterior_lengthx,Posterior_2d_w] = fun_tempModel(coeff_w_power,tempModel_options); %#ok<*ASGLU>

a = 1;

svm_outs = struct;
svm_outs.svm_posterior_lengthx = svm_posterior_lengthx;
svm_outs.Posterior_2d_n11n = Posterior_2d_n11n;
svm_outs.Posterior_2d_w = Posterior_2d_w;
svm_outs.temp_svm_Y_valid_T = temp_svm_Y_valid_T;
svm_outs.coeff_w_power = coeff_w_power; 
svm_outs.boolIndex_location_seq = boolIndex_location_seq;
svm_outs.seq_range = seq_range;





%% End
