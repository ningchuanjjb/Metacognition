% Chuan's 11th script (20251214)
% This script: To conduct analysis of relationship between memory strength and meta-memory, related to figure 4.
%% Initialization
memoryPrecision_trialLevel;
meta_trialLevel;

x_raw = memoryPrecision_trialLevel(choiceBoolIndex);
y_raw = meta_trialLevel(choiceBoolIndex);
x = x_raw(~isnan(x_raw));
y = y_raw(~isnan(x_raw));
% median(x)
% median(y)
lowThreshold_memoryPrecision;
% lowThreshold_meta = median(y);%0.5
% lowThreshold_meta = 0.5;

if if_twoThreshold_median0_optimal1 == 0
    lowThreshold_meta = median(y);
elseif if_twoThreshold_median0_optimal1 == 1
    lowThreshold_meta = metaDecoderThreshold_delay1;
end

highThreshold_meta = lowThreshold_meta;

if_smoothHistogram = 0;
histogram_numBins = 30;%10-->30
smoothCoeff = 128;%32

if_plot = 1;
if_plot_trialDistriA = 1;
if_plot_trialDistriB = 0;
if_plot_trialDistriC = 1;
if_plot_trialDistriD = 0;
if_plot_trialDistriE = 1;
if_plot_trialDistriE2 = 0;
if_plot_trialDistriF = 1;

if_plot_seqDistri = 0;
if_plot_exampleSeqDistri = 1;

if_colormap_loadEnhanced = 0;

if_trialDistriC_neuronLabel0_doubleLabel1 = 0;

if_plot_fineTuningMismatch = 0;%1

if_plot_memoryPrecision_trialLevelEvidence = if_plot_memoryPrecision_trialLevelEvidence; %#ok<*ASGSL>


% color_choiceMemoryHigh = [127,201,127]/255;%[0.5,1,0.5]
% color_choiceMemoryLow = [0,0.5,0];%[0,0.5,0]
% color_choiceOffloadHigh = [0.5,0,0];%[1,0.5,0.5]
% color_choiceOffloadLow = [1,0.5,0.5];%[0.5,0,0]

color_choiceMemoryHigh = [128,205,193]/255; %[146,197,222]/255
color_choiceMemoryLow = [1,133,113]/255; %[5,113,176]/255
color_choiceOffloadHigh = [166,97,26]/255; %[202,0,32]/255
color_choiceOffloadLow = [223,194,125]/255; %[244,165,130]/255


color_choiceMemory = [1,133+20,113]/255;
color_choiceOffload = [166+40,97,26]/255;
color_choiceMemoryError = [0.3 0.3 0.3];


%% Compute fine-tuning mismatch trial
choiceOffloadBoolIndex;
choiceMemoryBoolIndex;
choiceMemoryCorrectBoolIndex;
choiceMemoryErrorBoolIndex;

% trialBoolIndex_metaLow_choice_baseline
trialBoolIndex_metaLow_choice_baseline = false(trial_num,1);
trialBoolIndex_metaHigh_choice_baseline = false(trial_num,1);
for tempi=1:trial_num
    if choiceBoolIndex(tempi) == 0
        continue
    end
    %if meta_trialLevel_baseline(tempi) <= median(meta_trialLevel_baseline,'omitnan')
    if meta_trialLevel_baseline(tempi) <= metaDecoderThreshold_baseline
        trialBoolIndex_metaLow_choice_baseline(tempi) = true;
    end
    %if meta_trialLevel_baseline(tempi) > median(meta_trialLevel_baseline,'omitnan')
    if meta_trialLevel_baseline(tempi) > metaDecoderThreshold_baseline
        trialBoolIndex_metaHigh_choice_baseline(tempi) = true;
    end
end

a1_baseline = sum(trialBoolIndex_metaLow_choice_baseline);
a2_baseline = sum(trialBoolIndex_metaHigh_choice_baseline);

% trialBoolIndex_metaLow_choice
trialBoolIndex_metaLow_choice = false(trial_num,1);
trialBoolIndex_metaHigh_choice = false(trial_num,1);
for tempi=1:trial_num
    if choiceBoolIndex(tempi) == 0
        continue
    end
    if isnan(memoryPrecision_trialLevel(tempi))
        continue
    end
    if meta_trialLevel(tempi) <= lowThreshold_meta
        trialBoolIndex_metaLow_choice(tempi) = true;
    end
    if meta_trialLevel(tempi) > highThreshold_meta
        trialBoolIndex_metaHigh_choice(tempi) = true;
    end
end

a1 = sum(trialBoolIndex_metaLow_choice);
a2 = sum(trialBoolIndex_metaHigh_choice);


% match trials
trialBoolIndex_memoryPrecisionLow_choiceOffload;
trialBoolIndex_memoryPrecisionHigh_choiceMemory;
trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory;

trialBoolIndex_memoryPrecisionLow_choiceOffloadLow = trialBoolIndex_memoryPrecisionLow_choiceOffload & trialBoolIndex_metaLow_choice;
trialBoolIndex_memoryPrecisionHigh_choiceMemoryHigh = trialBoolIndex_memoryPrecisionHigh_choiceMemory & trialBoolIndex_metaHigh_choice;
trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh = trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory & trialBoolIndex_metaHigh_choice;

% mismatch trials
trialBoolIndex_memoryPrecisionLow_choiceMemory;
trialBoolIndex_memoryPrecisionLowError_choiceMemory;
trialBoolIndex_memoryPrecisionHigh_choiceOffload;

trialBoolIndex_memoryPrecisionLow_choiceMemoryHigh = trialBoolIndex_memoryPrecisionLow_choiceMemory & trialBoolIndex_metaHigh_choice;
trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh = trialBoolIndex_memoryPrecisionLowError_choiceMemory & trialBoolIndex_metaHigh_choice;
trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow = trialBoolIndex_memoryPrecisionHigh_choiceOffload & trialBoolIndex_metaLow_choice;


a1 = sum(trialBoolIndex_memoryPrecisionLow_choiceOffloadLow);
a2 = sum(trialBoolIndex_memoryPrecisionHigh_choiceMemoryHigh);
a3 = sum(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh);
a4 = sum(trialBoolIndex_memoryPrecisionLow_choiceMemoryHigh);
a5 = sum(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh);
a6 = sum(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow);


memoryMetaMismatch.trialBoolIndex_memoryPrecisionLow_choiceOffloadLow = trialBoolIndex_memoryPrecisionLow_choiceOffloadLow;
memoryMetaMismatch.trialBoolIndex_memoryPrecisionHigh_choiceMemoryHigh = trialBoolIndex_memoryPrecisionHigh_choiceMemoryHigh;
memoryMetaMismatch.trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh = trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh;
memoryMetaMismatch.trialBoolIndex_memoryPrecisionLow_choiceMemoryHigh = trialBoolIndex_memoryPrecisionLow_choiceMemoryHigh;
memoryMetaMismatch.trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh = trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh;
memoryMetaMismatch.trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow = trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow;


trialBoolIndex_allMatchMisMatch = trialBoolIndex_memoryPrecisionLow_choiceOffloadLow | ...
    trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh | ...
    trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh | ...
    trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow;


% if false
if true
    trialBoolIndex_memoryPrecisionLow_choice; %#ok<*UNRCH>
    trialBoolIndex_memoryPrecisionHigh_choice;
    trialBoolIndex_metaLow_choice;
    trialBoolIndex_metaHigh_choice;
    
    trialBoolIndex_memoryPrecisionLow_metaLow_choice = trialBoolIndex_memoryPrecisionLow_choice & trialBoolIndex_metaLow_choice;
    trialBoolIndex_memoryPrecisionHigh_metaHigh_choice = trialBoolIndex_memoryPrecisionHigh_choice & trialBoolIndex_metaHigh_choice;
    trialBoolIndex_memoryPrecisionLow_metaHigh_choice = trialBoolIndex_memoryPrecisionLow_choice & trialBoolIndex_metaHigh_choice;
    trialBoolIndex_memoryPrecisionHigh_metaLow_choice = trialBoolIndex_memoryPrecisionHigh_choice & trialBoolIndex_metaLow_choice;
    
    b1 = sum(trialBoolIndex_memoryPrecisionLow_metaLow_choice);
    b2 = sum(trialBoolIndex_memoryPrecisionHigh_metaHigh_choice);
    b3 = sum(trialBoolIndex_memoryPrecisionLow_metaHigh_choice);
    b4 = sum(trialBoolIndex_memoryPrecisionHigh_metaLow_choice);
    
end

if true
    
    trialBoolIndex_memoryPrecisionLow_metaLow_choice;
    trialBoolIndex_memoryPrecisionHigh_metaHigh_choice;
    trialBoolIndex_memoryPrecisionLow_metaHigh_choice;
    trialBoolIndex_memoryPrecisionHigh_metaLow_choice;
    
    % memoryPrecisionLow_metaLow
    trialBoolIndex_memoryPrecisionLow_metaLow_choiceMemoryCorrect = trialBoolIndex_memoryPrecisionLow_metaLow_choice & choiceMemoryCorrectBoolIndex';
    trialBoolIndex_memoryPrecisionLow_metaLow_choiceMemoryError = trialBoolIndex_memoryPrecisionLow_metaLow_choice & choiceMemoryErrorBoolIndex';
    trialBoolIndex_memoryPrecisionLow_metaLow_choiceOffload = trialBoolIndex_memoryPrecisionLow_metaLow_choice & choiceOffloadBoolIndex';
    a1 = sum(trialBoolIndex_memoryPrecisionLow_metaLow_choiceMemoryCorrect);
    a2 = sum(trialBoolIndex_memoryPrecisionLow_metaLow_choiceMemoryError);
    a3 = sum(trialBoolIndex_memoryPrecisionLow_metaLow_choiceOffload);
    
    % memoryPrecisionHigh_metaHigh
    trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryCorrect = trialBoolIndex_memoryPrecisionHigh_metaHigh_choice & choiceMemoryCorrectBoolIndex';
    trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError = trialBoolIndex_memoryPrecisionHigh_metaHigh_choice & choiceMemoryErrorBoolIndex';
    trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceOffload = trialBoolIndex_memoryPrecisionHigh_metaHigh_choice & choiceOffloadBoolIndex';
    b1 = sum(trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryCorrect);
    b2 = sum(trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError);
    b3 = sum(trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceOffload);
    
    % memoryPrecisionLow_metaHigh
    trialBoolIndex_memoryPrecisionLow_metaHigh_choiceMemoryCorrect = trialBoolIndex_memoryPrecisionLow_metaHigh_choice & choiceMemoryCorrectBoolIndex';
    trialBoolIndex_memoryPrecisionLow_metaHigh_choiceMemoryError = trialBoolIndex_memoryPrecisionLow_metaHigh_choice & choiceMemoryErrorBoolIndex';
    trialBoolIndex_memoryPrecisionLow_metaHigh_choiceOffload = trialBoolIndex_memoryPrecisionLow_metaHigh_choice & choiceOffloadBoolIndex';
    c1 = sum(trialBoolIndex_memoryPrecisionLow_metaHigh_choiceMemoryCorrect);
    c2 = sum(trialBoolIndex_memoryPrecisionLow_metaHigh_choiceMemoryError);
    c3 = sum(trialBoolIndex_memoryPrecisionLow_metaHigh_choiceOffload);
    
    
    % memoryPrecisionHigh_metaLow
    trialBoolIndex_memoryPrecisionHigh_metaLow_choiceMemoryCorrect = trialBoolIndex_memoryPrecisionHigh_metaLow_choice & choiceMemoryCorrectBoolIndex';
    trialBoolIndex_memoryPrecisionHigh_metaLow_choiceMemoryError = trialBoolIndex_memoryPrecisionHigh_metaLow_choice & choiceMemoryErrorBoolIndex';
    trialBoolIndex_memoryPrecisionHigh_metaLow_choiceOffload = trialBoolIndex_memoryPrecisionHigh_metaLow_choice & choiceOffloadBoolIndex';
    d1 = sum(trialBoolIndex_memoryPrecisionHigh_metaLow_choiceMemoryCorrect);
    d2 = sum(trialBoolIndex_memoryPrecisionHigh_metaLow_choiceMemoryError);
    d3 = sum(trialBoolIndex_memoryPrecisionHigh_metaLow_choiceOffload);
    
    
    % trialBoolIndex_precisionLow_metaHigh_choiceMemoryError_response
    if if_memeoryPrecision_stimuli0_response1 == 1
        trialBoolIndex_precisionLow_metaHigh_choiceMemoryError_response = ...
            trialBoolIndex_memoryPrecisionLow_metaHigh_choiceMemoryError;
    end
    
    % to confirm that over-mismatch trials of response-labeled are stable when in stimuli-labeled
    % to confirm that over-mismatch trials don't distribute only in few seqs but in a broad range of seqs.
    if false
        trialBoolIndex_precisionLow_metaHigh_choiceMemoryError_response;
        
        temp_a2 = sum(trialBoolIndex_memoryPrecisionLow_metaLow_choice & ...
            trialBoolIndex_precisionLow_metaHigh_choiceMemoryError_response);
        
        temp_b2 = sum(trialBoolIndex_memoryPrecisionHigh_metaHigh_choice & ...
            trialBoolIndex_precisionLow_metaHigh_choiceMemoryError_response);
        
        temp_c2 = sum(trialBoolIndex_memoryPrecisionLow_metaHigh_choice & ...
            trialBoolIndex_precisionLow_metaHigh_choiceMemoryError_response);
        
        temp_d2 = sum(trialBoolIndex_memoryPrecisionHigh_metaLow_choice & ...
            trialBoolIndex_precisionLow_metaHigh_choiceMemoryError_response);
        
        temp_bdac2 = [temp_a2+temp_b2+temp_c2+temp_d2 temp_b2 temp_d2 temp_a2 temp_c2];
        temp_bdac2_n11n = [temp_a2+temp_b2+temp_c2+temp_d2 temp_b2 temp_d2 temp_a2 temp_c2]./temp_bdac2(1);
        
        
        temp_trialBoolIndex_c2 = trialBoolIndex_memoryPrecisionLow_metaHigh_choice & ...
            trialBoolIndex_precisionLow_metaHigh_choiceMemoryError_response;
        
        temp_seqIndex_c2 = seqIndex(temp_trialBoolIndex_c2);
        
        temp_seqIndex_c2_sorted = sort(temp_seqIndex_c2);
        temp_boolIndex_location_seq_c2_sorted = boolIndex_location_seq_T(temp_seqIndex_c2_sorted,:);
        
        %sum(temp_boolIndex_location_seq_c2_sorted,1)
        
        temp_seqIndex_c2_unique = unique(temp_seqIndex_c2_sorted);
        
        temp_seqIndex_c2_unique_seqCount = length(temp_seqIndex_c2_unique);
        temp_seqIndex_c2_trialCount = length(temp_seqIndex_c2);
        
    end
    
    
    %     abcd1 = [a1+b1+c1+d1 a1 b1 c1 d1];
    %     abcd2 = [a2+b2+c2+d2 a2 b2 c2 d2];
    %     abcd3 = [a3+b3+c3+d3 a3 b3 c3 d3];
    
    bdac1 = [a1+b1+c1+d1 b1 d1 a1 c1];
    bdac2 = [a2+b2+c2+d2 b2 d2 a2 c2];
    bdac3 = [a3+b3+c3+d3 b3 d3 a3 c3];
    
    bdac1_n11n = [a1+b1+c1+d1 b1 d1 a1 c1]./bdac1(1);
    bdac2_n11n = [a2+b2+c2+d2 b2 d2 a2 c2]./bdac2(1);
    bdac3_n11n = [a3+b3+c3+d3 b3 d3 a3 c3]./bdac3(1);
    
    
    %if false
    if true
        bdac1_n11n_B = bdac1_n11n(2:end);
        bdac2_n11n_B = bdac2_n11n(2:end);
        bdac3_n11n_B = bdac3_n11n(2:end);
        
        bdac1_template = [1 0 0 0];
        bdac2_template = [0.5 0 0 0.5];
        %bdac2_template = [0 0 0 1];
        %bdac2_template = [1 0 0 0];
        
        bdac2_template_A = [1 0 0 0];
        bdac2_template_B = [0.5 0 0 0.5];
        
        
        bdac3_template = [0 0 1 0];
        
        [r_bdac1,p_bdac1] = corr(bdac1_n11n_B',bdac1_template');
        [r_bdac2,p_bdac2] = corr(bdac2_n11n_B',bdac2_template');
        [r_bdac3,p_bdac3] = corr(bdac3_n11n_B',bdac3_template');
        
        [r_bdac2_A,p_bdac2_A] = corr(bdac2_n11n_B',bdac2_template_A');
        [r_bdac2_B,p_bdac2_B] = corr(bdac2_n11n_B',bdac2_template_B');
        
        
        [r_bdac123,p_bdac123] = corr([bdac1_n11n_B,bdac2_n11n_B,bdac3_n11n_B]',...
            [bdac1_template,bdac2_template,bdac3_template]');
        
        b_132 = [bdac1_n11n_B(1),bdac3_n11n_B(1),bdac2_n11n_B(1)];
        a_132 = [bdac1_n11n_B(3),bdac3_n11n_B(3),bdac2_n11n_B(3)];
        c_132 = [bdac1_n11n_B(4),bdac3_n11n_B(4),bdac2_n11n_B(4)];
        
        b_132_template = [bdac1_template(1),bdac3_template(1),bdac2_template(1)];
        a_132_template = [bdac1_template(3),bdac3_template(3),bdac2_template(3)];
        c_132_template = [bdac1_template(4),bdac3_template(4),bdac2_template(4)];
        
        [r_b,p_b] = corr(b_132',b_132_template');
        [r_a,p_a] = corr(a_132',a_132_template');
        [r_c,p_c] = corr(c_132',c_132_template');
        
        
        
        
    end
    
end



% %% Compare accuracy of OverMismatch (correct+error) trials and HighMatch (correct+error) trials
%% Compare ChoiceMemory trial accuracy in LowPrecision and HighPrecision
% choiceMemoryCorrectBoolIndex;
% choiceMemoryErrorBoolIndex;
trialBoolIndex_memoryPrecisionHigh_choiceMemoryHigh;
trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh;

trialBoolIndex_memoryPrecisionLow_choiceMemoryHigh;
trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh;

trialBoolIndex_memoryPrecisionLowCorrect_choiceMemoryHigh = trialBoolIndex_memoryPrecisionLow_choiceMemoryHigh & choiceMemoryCorrectBoolIndex';
trialBoolIndex_memoryPrecisionLowCorrect_choiceMemory = trialBoolIndex_memoryPrecisionLow_choiceMemory & choiceMemoryCorrectBoolIndex';


% temp1 = sum(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh)./sum(trialBoolIndex_memoryPrecisionHigh_choiceMemoryHigh);
% temp2 = sum(trialBoolIndex_memoryPrecisionLowCorrect_choiceMemoryHigh)./sum(trialBoolIndex_memoryPrecisionLow_choiceMemoryHigh);
temp1 = sum(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory)./sum(trialBoolIndex_memoryPrecisionHigh_choiceMemory);
temp2 = sum(trialBoolIndex_memoryPrecisionLowCorrect_choiceMemory)./sum(trialBoolIndex_memoryPrecisionLow_choiceMemory);
fprintf('Accuracy of HighPrecision is %.3f.\n',temp1);
fprintf('Accuracy of LowPrecision is %.3f.\n',temp2);

% temptempIndexA1 = find(trialBoolIndex_memoryPrecisionHigh_choiceMemoryHigh == true);
% temptempIndexA2 = find(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh == true);
temptempIndexA1 = find(trialBoolIndex_memoryPrecisionHigh_choiceMemory == true);
temptempIndexA2 = find(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory == true);
temptempBoolIndexA = ismember(temptempIndexA1,temptempIndexA2);

% temptempIndexB1 = find(trialBoolIndex_memoryPrecisionLow_choiceMemoryHigh == true);
% temptempIndexB2 = find(trialBoolIndex_memoryPrecisionLowCorrect_choiceMemoryHigh == true);
temptempIndexB1 = find(trialBoolIndex_memoryPrecisionLow_choiceMemory == true);
temptempIndexB2 = find(trialBoolIndex_memoryPrecisionLowCorrect_choiceMemory == true);
temptempBoolIndexB = ismember(temptempIndexB1,temptempIndexB2);

correctBoolIndex_HighPrecision = temptempBoolIndexA;
correctBoolIndex_LowPrecision = temptempBoolIndexB;
[~,temp_p_LowPrecision_HighPrecision] = ttest2(correctBoolIndex_HighPrecision,correctBoolIndex_LowPrecision);


%% Compare FreeChoice trial rProb in LowMeta and HighMeta
trialBoolIndex_metaLow_choice;
trialBoolIndex_metaHigh_choice;

trialBoolIndex_metaLow_choiceOffload = trialBoolIndex_metaLow_choice & choiceOffloadBoolIndex';
trialBoolIndex_metaHigh_choiceOffload = trialBoolIndex_metaHigh_choice & choiceOffloadBoolIndex';

temp1 = sum(trialBoolIndex_metaLow_choiceOffload)./sum(trialBoolIndex_metaLow_choice);
temp2 = sum(trialBoolIndex_metaHigh_choiceOffload)./sum(trialBoolIndex_metaHigh_choice);
fprintf('Offloading rate of LowMeta is %.3f.\n',temp1);
fprintf('Offloading rate of HighMeta is %.3f.\n',temp2);

temptempIndexA1 = find(trialBoolIndex_metaLow_choice == true);
temptempIndexA2 = find(trialBoolIndex_metaLow_choiceOffload == true);
temptempBoolIndexA = ismember(temptempIndexA1,temptempIndexA2);

temptempIndexB1 = find(trialBoolIndex_metaHigh_choice == true);
temptempIndexB2 = find(trialBoolIndex_metaHigh_choiceOffload == true);
temptempBoolIndexB = ismember(temptempIndexB1,temptempIndexB2);

choiceOffloadBoolIndex_metaLow = temptempBoolIndexA;
choiceOffloadBoolIndex_metaHigh = temptempBoolIndexB;
[~,temp_p_metaLow_metaHigh] = ttest2(choiceOffloadBoolIndex_metaLow,choiceOffloadBoolIndex_metaHigh);



%% Seq count distribution of mismatch trials
trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh;
trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow;

if if_memeoryPrecision_stimuli0_response1 == 0
    seqIndex_current = seqIndex;
elseif if_memeoryPrecision_stimuli0_response1 == 1
    seqIndex_current = seqIndex_response;
end

% tempSeqIndex1 = seqIndex_current(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh)';
% tempSeqIndex2 = seqIndex_current(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow)';

tempSeqIndex1 = sort(seqIndex_current(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh)'); %#ok<*TRSRT>
tempSeqIndex2 = sort(seqIndex_current(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow)');
tempSeq1 = target_seqSet_inOne(tempSeqIndex1)';
tempSeq2 = target_seqSet_inOne(tempSeqIndex2)';

tempSeqIndexB1 = sort(seqIndex_current(trialBoolIndex_memoryPrecisionLowError_choiceMemory)'); %#ok<*TRSRT>
tempSeqIndexB2 = sort(seqIndex_current(trialBoolIndex_memoryPrecisionHigh_choiceOffload)');
tempSeqB1 = target_seqSet_inOne(tempSeqIndexB1)';
tempSeqB2 = target_seqSet_inOne(tempSeqIndexB2)';

tempSeqCount1 = zeros(sum(numSeq(1:3)),1);
tempSeqCount2 = zeros(sum(numSeq(1:3)),1);
tempSeqCountB1 = zeros(sum(numSeq(1:3)),1);
tempSeqCountB2 = zeros(sum(numSeq(1:3)),1);
for tempi=1:length(tempSeqCount1)
    tempSeqCount1(tempi) = sum(tempSeqIndex1==tempi);
    tempSeqCount2(tempi) = sum(tempSeqIndex2==tempi);
    tempSeqCountB1(tempi) = sum(tempSeqIndexB1==tempi);
    tempSeqCountB2(tempi) = sum(tempSeqIndexB2==tempi);
end

seqCount_memoryPrecisionLowError_choiceMemoryHigh = tempSeqCount1;
seqCount_memoryPrecisionHigh_choiceOffloadLow = tempSeqCount2;
seqCount_memoryPrecisionLowError_choiceMemory = tempSeqCountB1;
seqCount_memoryPrecisionHigh_choiceOffload = tempSeqCountB2;

a1 = sum(seqCount_memoryPrecisionLowError_choiceMemoryHigh);
a2 = sum(seqCount_memoryPrecisionHigh_choiceOffloadLow);
a3 = sum(seqCount_memoryPrecisionLowError_choiceMemory);
a4 = sum(seqCount_memoryPrecisionHigh_choiceOffload);



%% seqSet_inOne_inOne
seqSet_inOne = cell(1,4);
for target_seqLength = 1:4
    seqSet_inOne{target_seqLength} = zeros(1, numSeq(target_seqLength));
    for tempi=1:numSeq(target_seqLength)
        for tempj=1:target_seqLength
            seqSet_inOne{target_seqLength}(tempi) = 10^(target_seqLength-tempj) * ...
                target_seqSet{target_seqLength}{tempi}(tempj) + seqSet_inOne{target_seqLength}(tempi);
        end
    end
end
seqSet_inOne_inOne = [];
for tempi=1:4
    seqSet_inOne_inOne = [seqSet_inOne_inOne; seqSet_inOne{tempi}']; %#ok<*AGROW>
end


%% Plot
if if_plot == 1
    close all
    
    %     %% All trial
    %     if if_plot_trialDistriA == 1
    %         fig = figure('Name','Trial distributionA','NumberTitle','off'); %#ok<*NASGU>
    %         % set(gcf,'Position',[10 50 1680 830]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %         %set(gcf,'Position',[10 50 1250 830]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %         set(gcf,'Position',[10 50 720 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %
    %         t = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
    %         t.Title.String = sprintf(sprintf('Trial distribution, %s',FOVName_currentFOV2));
    %         t.Title.FontSize = 12;
    %         t.Title.Interpreter = 'none';
    %
    %         for loopCount=1:5
    %             nexttile
    %
    %             if loopCount==1
    %                 tempStr = 'free choice';
    %                 x_raw = memoryPrecision_trialLevel(choiceBoolIndex);
    %                 y_raw = meta_trialLevel(choiceBoolIndex);
    %             elseif loopCount==2
    %                 tempStr = 'choiceMemory';
    %                 x_raw = memoryPrecision_trialLevel(choiceMemoryBoolIndex);
    %                 y_raw = meta_trialLevel(choiceMemoryBoolIndex);
    %             elseif loopCount==3
    %                 tempStr = 'choiceOffload';
    %                 x_raw = memoryPrecision_trialLevel(choiceOffloadBoolIndex);
    %                 y_raw = meta_trialLevel(choiceOffloadBoolIndex);
    %             elseif loopCount==4
    %                 tempStr = 'choiceMemoryCorrect';
    %                 x_raw = memoryPrecision_trialLevel(choiceMemoryCorrectBoolIndex);
    %                 y_raw = meta_trialLevel(choiceMemoryCorrectBoolIndex);
    %             elseif loopCount==5
    %                 tempStr = 'choiceMemoryError';
    %                 x_raw = memoryPrecision_trialLevel(choiceMemoryErrorBoolIndex);
    %                 y_raw = meta_trialLevel(choiceMemoryErrorBoolIndex);
    %             end
    %
    %             x = x_raw(~isnan(x_raw));
    %             y = y_raw(~isnan(x_raw));
    %
    %             [x_min,x_max] = bounds(x);
    %             x_max = x_max + 0.01;
    %
    %             if if_smoothHistogram == 0
    %                 XYBinLimits = [0 1];
    %
    %                 h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
    %                 h.NumBins = [10 10];
    %                 h.XBinLimits = XYBinLimits;
    %                 h.YBinLimits = XYBinLimits;
    %                 hold on
    %
    %                 plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],XYBinLimits,...
    %                     'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
    %                 hold on
    %                 plot(XYBinLimits,[lowThreshold_meta lowThreshold_meta],...
    %                     'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
    %                 hold on
    %
    %                 xticks(0:0.2:1);
    %                 %xticks(x_min:0.2:x_max);
    %                 yticks(0:0.2:1);
    %
    %             elseif if_smoothHistogram == 1
    %                 n=smoothCoeff;%50
    %                 n=2^ceil(log2(n)); % round up n to the next power of 2;
    %                 gridx1 = linspace(0,1,n);
    %                 %gridx1 = linspace(0,x_max,n);
    %                 gridx2 = linspace(0,1,n);
    %
    %                 XYBinLimits = [1 n];
    %
    %                 [tempx1,tempx2] = meshgrid(gridx1, gridx2);
    %                 tempx1 = tempx1(:);
    %                 tempx2 = tempx2(:);
    %                 xi = [tempx1 tempx2];
    %                 [pdf1,xmesh1,bandwidth1] = ksdensity([x,y],xi,'Function','pdf','Support',[0 0;1 1],'BoundaryCorrection','Reflection');
    %                 %[pdf1,xmesh1,bandwidth1] = ksdensity([x,y],xi,'Function','pdf','Support',[0 0;x_max 1],'BoundaryCorrection','Reflection');
    %
    %                 pdf1_2d = reshape(pdf1,n,n);
    %
    %                 if loopCount == 1
    %                     %pdf_max = max(pdf1_2d,[],'all') * 1.1;
    %                     %pdf_max = max(pdf1_2d,[],'all') * 1;
    %
    %                     % edit in 20240126
    %                     pdf_max = max(pdf1_2d,[],'all') * 0.5;%0.5
    %                 end
    %
    %                 C = pdf1_2d;
    %                 %imagesc([1 size(C,2)],[size(C,1) 1],C);
    %                 %imagesc([1 size(C,2)],[size(C,1) 1],C,[0 5]);
    %                 %imagesc([1 size(C,2)],[size(C,1) 1],C,[0 pdf_max]);
    %                 imagesc([1 size(C,2)],[1 size(C,1)],C,[0 pdf_max]);
    %                 hold on
    %                 plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision].*(n+1),XYBinLimits,...
    %                     'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
    %                 hold on
    %                 %plot(XYBinLimits,[lowThreshold_meta lowThreshold_meta].*n,...
    %                 %    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
    %                 %plot(XYBinLimits,[1-lowThreshold_meta 1-lowThreshold_meta].*(n+1),...
    %                 %    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
    %                 plot(XYBinLimits,[lowThreshold_meta lowThreshold_meta].*(n+1),...
    %                    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
    %                 hold on
    %
    %                 xticks(linspace(1,n,6));
    %                 yticks(linspace(1,n,6));
    %                 xticklabels({'0','0.2','0.4','0.6','0.8','1'});
    %                 yticklabels({'0','0.2','0.4','0.6','0.8','1'});
    %                 %yticklabels({'1','0.8','0.6','0.4','0.2','0'});
    %             end
    %
    %             axis equal
    %             set(gca,'Ydir','normal') %reverse
    %             xlim(XYBinLimits);
    %             ylim(XYBinLimits);
    %             xtickangle(0);
    %             set(gca, 'FontSize', 11) %12
    %             if loopCount == 5
    %                 c = colorbar('FontSize',10,'location','layout');
    %                 c.Position(1) = c.Position(1) + 0.15;
    %             end
    %             %c.Position(1) = c.Position(1) -0.005;
    %             %if loopCount == 4 || loopCount == 5
    %                 xlabel('Memory precision', 'FontSize', 12, 'FontWeight', 'bold');
    %             %end
    %             if loopCount == 1 || loopCount == 4
    %                 ylabel('Meta-memory', 'FontSize', 12, 'FontWeight', 'bold');
    %             end
    %
    %             temp_title = title(sprintf('%s, %d trials',tempStr,length(x)),'FontSize',10,'FontWeight','bold');
    %             temp_title.Interpreter = 'none';
    %         end
    %
    %     end
    
    
    
    %% All trial
    if if_plot_trialDistriA == 1
        fig = figure('Name','Trial distributionA','NumberTitle','off'); %#ok<*NASGU>
        %set(gcf,'Position',[10 50 720 500*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 50 720*0.92*0.95 500*0.85*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        %         %t = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
        %         t = tiledlayout(2,3,'TileSpacing','tight','Padding','tight');
        %         t.Title.String = sprintf(sprintf('Trial distribution, %s',FOVName_currentFOV2));
        %         t.Title.FontSize = 12;
        %         t.Title.Interpreter = 'none';
        
        temp_n = 6;%6
        
        t = tiledlayout(2,3*temp_n,'TileSpacing','loose','Padding','compact');
        
        
        for loopCount=1:5
            nexttile([1 temp_n-1])
            %nexttile
            
            if loopCount==1
                tempStr = 'Choice';
                x_raw = memoryPrecision_trialLevel(choiceBoolIndex);
                y_raw = meta_trialLevel(choiceBoolIndex);
            elseif loopCount==2
                tempStr = 'Memory';
                x_raw = memoryPrecision_trialLevel(choiceMemoryBoolIndex);
                y_raw = meta_trialLevel(choiceMemoryBoolIndex);
            elseif loopCount==3
                tempStr = 'Offload';
                x_raw = memoryPrecision_trialLevel(choiceOffloadBoolIndex);
                y_raw = meta_trialLevel(choiceOffloadBoolIndex);
            elseif loopCount==4
                tempStr = 'Memory-correct';
                x_raw = memoryPrecision_trialLevel(choiceMemoryCorrectBoolIndex);
                y_raw = meta_trialLevel(choiceMemoryCorrectBoolIndex);
            elseif loopCount==5
                tempStr = 'Memory-error';
                x_raw = memoryPrecision_trialLevel(choiceMemoryErrorBoolIndex);
                y_raw = meta_trialLevel(choiceMemoryErrorBoolIndex);
            end
            
            x = x_raw(~isnan(x_raw));
            y = y_raw(~isnan(x_raw));
            
            [x_min,x_max] = bounds(x);
            x_max = x_max + 0.01;
            
            if if_smoothHistogram == 0
                XYBinLimits = [0 1];
                
                %h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
                h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','count');
                %h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','off','Normalization','count');
                %h.NumBins = [10 10];
                h.XBinLimits = XYBinLimits;
                h.YBinLimits = XYBinLimits;
                h.NumBins = [1 1].*histogram_numBins;
                hold on
                
                plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],XYBinLimits,...
                    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
                hold on
                plot(XYBinLimits,[lowThreshold_meta lowThreshold_meta],...
                    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
                hold on
                
                xticks(0:0.2:1);
                yticks(0:0.2:1);
                
            elseif if_smoothHistogram == 1
                n=smoothCoeff;%50
                n=2^ceil(log2(n)); % round up n to the next power of 2;
                gridx1 = linspace(0,1,n);
                gridx2 = linspace(0,1,n);
                
                XYBinLimits = [1 n];
                
                [tempx1,tempx2] = meshgrid(gridx1, gridx2);
                tempx1 = tempx1(:);
                tempx2 = tempx2(:);
                xi = [tempx1 tempx2];
                %[pdf1,xmesh1,bandwidth1] = ksdensity([x,y],xi,'Function','pdf','Support',[0 0;1 1],'BoundaryCorrection','Reflection');
                [pdf1,xmesh1,bandwidth1] = ksdensity([x,y],xi,'Function','pdf','Support',[-0.01 -0.01;1.01 1.01],'BoundaryCorrection','Reflection');
                
                pdf1_2d = reshape(pdf1,n,n);
                
                %if loopCount == 1
                %pdf_max = max(pdf1_2d,[],'all') * 1.1;
                pdf_max = max(pdf1_2d,[],'all') * 1;
                
                % edit in 20240126
                %pdf_max = max(pdf1_2d,[],'all') * 0.5;%0.5
                %end
                
                C = pdf1_2d;
                imagesc([1 size(C,2)],[1 size(C,1)],C,[0 pdf_max]);
                hold on
                plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision].*(n+1),XYBinLimits,...
                    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
                hold on
                plot(XYBinLimits,[lowThreshold_meta lowThreshold_meta].*(n+1),...
                    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
                hold on
                
                xticks(linspace(1,n,6));
                yticks(linspace(1,n,6));
                %if loopCount == 1 || loopCount == 2 || loopCount == 3
                %    xticklabels('');
                %end
                %if loopCount == 4 || loopCount == 5
                xticklabels({'0','0.2','0.4','0.6','0.8','1'});
                %end
                yticklabels({'0','0.2','0.4','0.6','0.8','1'});
            end
            
            %axis equal
            set(gca,'Ydir','normal') %reverse
            xlim(XYBinLimits);
            ylim(XYBinLimits);
            xtickangle(0);
            set(gca, 'FontSize', 8) %12
            %if loopCount == 5
            %c = colorbar('FontSize',10,'location','layout');%10
            %c.Position(1) = c.Position(1) + 0.15;%0.15
            %end
            
            c = colorbar('FontSize',6.5);
            %             %c.Position(1) = c.Position(1) + 0.01;%0.15
            %             %c.Position = c.Position+[0 0 0 -0.05];
            %
            %             c = colorbar('eastoutside','FontSize',6.5);
            %             %c.Position = c.Position+[+0.05 0 0 0];
            %
            %c.Position = c.Position+[0 0 0.0001 0];
            
            if loopCount == 1
                c.Position = c.Position+[0.015 0.006 -0.005 -0.032];
            elseif loopCount == 2
                c.Position = c.Position+[0.008 0.018 -0.005 -0.0195];
            elseif loopCount == 3
                c.Position = c.Position+[-0.018 0.018 -0.005 -0.0195];
            elseif loopCount == 4
                c.Position = c.Position+[0.035 0.039 -0.005 -0.021];
            elseif loopCount == 5
                c.Position = c.Position+[0.010 0 -0.005 -0.000];
            end
            
            if if_colormap_loadEnhanced == 1
                load('parula_enhanced');
                colormap(parula_enhanced);
            elseif if_colormap_loadEnhanced == 0
                %colormap parula
                colormap(coolwarm());
                
                %temp1 = gray;
                %temp1 = temp1(end:-1:1,:);
                %colormap(temp1);
            end
            
            if loopCount == 4 || loopCount == 5
                xlabel('WM strength', 'FontSize', 9);
            end
            if loopCount == 1 || loopCount == 4
                ylabel('Meta-WM', 'FontSize', 9);
            end
            
            temp_title = title(sprintf('%s, %d trials',tempStr,length(x)),'FontSize',9);
            temp_title.Interpreter = 'none';
            
            %if loopCount ~= 3
            nexttile([1 1])
            set(gca, 'visible', 'off')
            %end
            %if loopCount == 3
            %    a = 1;
            %end
            
        end
        
    end
    
    %% plot_trialDistriB
    if if_plot_trialDistriB == 1
        fig = figure('Name','Trial distribution B','NumberTitle','off'); %#ok<*NASGU>
        % set(gcf,'Position',[10 50 1680 830]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 50 720 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 50 720 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 50 720 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        
        t = tiledlayout(3,4,'TileSpacing','loose','Padding','Compact');
        t.Title.String = sprintf(sprintf('Trial distribution, %s',FOVName_currentFOV2));
        t.Title.FontSize = 12;
        t.Title.Interpreter = 'none';
        
        for loopCount=1:4
            %nexttile
            
            if loopCount==1
                nexttile(1)
                %nexttile(4)
                %tempStr = 'ChoiceMemoryError';
                tempStr = 'MemoryError';
                x_raw = memoryPrecision_trialLevel(choiceMemoryErrorBoolIndex);
                y_raw = meta_trialLevel(choiceMemoryErrorBoolIndex);
            elseif loopCount==2
                nexttile(5)
                %nexttile(8)
                %tempStr = 'ChoiceMemoryCorrect';
                tempStr = 'MemoryCorrect';
                x_raw = memoryPrecision_trialLevel(choiceMemoryCorrectBoolIndex);
                y_raw = meta_trialLevel(choiceMemoryCorrectBoolIndex);
            elseif loopCount==3
                nexttile(9)
                %nexttile(12)
                %tempStr = 'Choice offload';
                tempStr = 'Offload';
                x_raw = memoryPrecision_trialLevel(choiceOffloadBoolIndex);
                y_raw = meta_trialLevel(choiceOffloadBoolIndex);
            elseif loopCount==4
                nexttile([3 3])
                tempStr = 'Merged';
                x_raw = memoryPrecision_trialLevel(choiceOffloadBoolIndex);
                y_raw = meta_trialLevel(choiceOffloadBoolIndex);
            end
            
            x = x_raw(~isnan(x_raw));
            y = y_raw(~isnan(x_raw));
            
            if if_smoothHistogram == 0
                XYBinLimits = [0 1];
                
                h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
                h.NumBins = [10 10];
                h.XBinLimits = XYBinLimits;
                h.YBinLimits = XYBinLimits;
                hold on
                
                plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],XYBinLimits,...
                    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
                hold on
                plot(XYBinLimits,[lowThreshold_meta lowThreshold_meta],...
                    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
                hold on
                
                xticks(0:0.2:1);
                yticks(0:0.2:1);
                
            elseif if_smoothHistogram == 1
                n=smoothCoeff;%50
                n=2^ceil(log2(n)); % round up n to the next power of 2;
                gridx1 = linspace(0,1,n);
                gridx2 = linspace(0,1,n);
                
                XYBinLimits = [1 n];
                
                [tempx1,tempx2] = meshgrid(gridx1, gridx2);
                tempx1 = tempx1(:);
                tempx2 = tempx2(:);
                xi = [tempx1 tempx2];
                %[pdf1,xmesh1,bandwidth1] = ksdensity([x,y],xi,'Function','pdf','Support',[0 0;1 1],'BoundaryCorrection','Reflection');
                [pdf1,xmesh1,bandwidth1] = ksdensity([x,y],xi,'Function','pdf','Support',[-0.01 -0.01;1.01 1.01],'BoundaryCorrection','Reflection');
                
                pdf1_2d = reshape(pdf1,n,n);
                
                if loopCount == 1
                    %pdf_max = max(pdf1_2d,[],'all') * 1.1;
                    
                    % edit in 20240126
                    pdf_max = max(pdf1_2d,[],'all') * 0.85;
                end
                
                if loopCount < 4
                    C = pdf1_2d;
                elseif loopCount == 4
                    C = C_muilti(:,:,1) .* template_boolIndex_multi(:,:,1) + ...
                        C_muilti(:,:,2) .* template_boolIndex_multi(:,:,2) + ...
                        C_muilti(:,:,3) .* template_boolIndex_multi(:,:,3);
                end
                
                
                
                %imagesc([1 size(C,2)],[size(C,1) 1],C);
                %imagesc([1 size(C,2)],[size(C,1) 1],C,[0 5]);
                %imagesc([1 size(C,2)],[size(C,1) 1],C,[0 pdf_max]);
                imagesc([1 size(C,2)],[1 size(C,1)],C,[0 pdf_max]);
                hold on
                
                if loopCount == 1
                    C_muilti = nan(size(C,1),size(C,1),3);
                end
                if loopCount < 4
                    C_muilti(:,:,loopCount) = C;
                end
                
                
                % template
                
                %temp_threshold_meta = round((1-lowThreshold_meta)*(n+1));
                temp_threshold_meta = ceil(lowThreshold_meta*(n+1));
                temp_threshold_memoryPrecision = ceil((lowThreshold_memoryPrecision)*(n+1));
                
                if loopCount < 4
                    %C_template = zeros(size(C));
                    C_template = C;
                    
                    template_boolIndex = false(size(C));
                    
                    if loopCount == 1
                        template_boolIndex_multi = false(size(C,1),size(C,1),3);
                        %C_muilti = nan(size(C,1),size(C,1),3);
                        
                        template_boolIndex(temp_threshold_meta:end,1:temp_threshold_memoryPrecision) = true;
                        template_boolIndex_multi(:,:,1) = template_boolIndex;
                        %C_muilti(:,:,1) = C;
                    elseif loopCount == 2
                        template_boolIndex(temp_threshold_meta:end,temp_threshold_memoryPrecision:end) = true;
                        template_boolIndex_multi(:,:,2) = template_boolIndex;
                        %C_muilti(:,:,1) = C;
                    elseif loopCount == 3
                        template_boolIndex(1:temp_threshold_meta,:) = true;
                        template_boolIndex_multi(:,:,3) = template_boolIndex;
                        %C_muilti(:,:,1) = C;
                    end
                    C_template = C_template .* template_boolIndex;
                    
                    
                    %imagesc('XData',[1 size(C_template,2)],'YData',[size(C_template,1) 1],...
                    %    'CData',C_template,'AlphaData',0.5,[0 pdf_max]);
                    imagesc('XData',[1 size(C_template,2)],'YData',[1 size(C_template,1)],...
                        'CData',C_template,'AlphaData',0.5,[0 pdf_max]);
                    hold on
                end
                
                
                temp_offset = 3;%1
                temp_lineWidth = 4;%5
                
                if loopCount == 1 || loopCount == 4
                    plot([1 1].*(temp_threshold_memoryPrecision-temp_offset),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                    plot([1 1].*(XYBinLimits(1)+temp_offset),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                    plot([1+temp_offset temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(temp_threshold_meta+temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                    plot([1+temp_offset temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(XYBinLimits(2)-temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                end
                if loopCount == 2 || loopCount == 4
                    plot([1 1].*(temp_threshold_memoryPrecision+temp_offset),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                    plot([1 1].*(XYBinLimits(2)-temp_offset),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset],...
                        [1 1].*(temp_threshold_meta+temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset],...
                        [1 1].*(XYBinLimits(2)-temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                end
                if loopCount == 3 || loopCount == 4
                    plot([1 1].*(temp_threshold_memoryPrecision-temp_offset),...
                        [XYBinLimits(1)+temp_offset temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    plot([1 1].*(XYBinLimits(1)+temp_offset),...
                        [XYBinLimits(1)+temp_offset temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    plot([1+temp_offset temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(temp_threshold_meta-temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    plot([1+temp_offset temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(XYBinLimits(1)+temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    
                    
                    plot([1 1].*(temp_threshold_memoryPrecision+temp_offset),...
                        [XYBinLimits(1)+temp_offset temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                    plot([1 1].*(XYBinLimits(2)-temp_offset),...
                        [XYBinLimits(1)+temp_offset temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset],...
                        [1 1].*(temp_threshold_meta-temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset],...
                        [1 1].*(XYBinLimits(1)+temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                end
                
                
                %plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision].*(n+1),XYBinLimits,...
                %    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
                %hold on
                %plot(XYBinLimits,[1-lowThreshold_meta 1-lowThreshold_meta].*(n+1),...
                %    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
                %hold on
                
                
                if loopCount < 4
                    xticks(linspace(1,n,2));
                    yticks(linspace(1,n,2));
                    %xticklabels({'0','1'});
                    %yticklabels({'1','0'});
                    set(gca,'xticklabel','');
                    set(gca,'yticklabel','');
                else
                    xticks(linspace(1,n,6));
                    yticks(linspace(1,n,6));
                    xticklabels({'0','0.2','0.4','0.6','0.8','1'});
                    %yticklabels({'1','0.8','0.6','0.4','0.2','0'});
                    yticklabels({'0','0.2','0.4','0.6','0.8','1'});
                end
            end
            
            set(gca,'Ydir','normal') %reverse
            set(gca,'linewidth',1.5)
            axis equal
            xlim(XYBinLimits);
            ylim(XYBinLimits);
            set(gca, 'FontSize', 12) %12
            if loopCount == 4
                c = colorbar('FontSize',10);
                %c = colorbar('FontSize',10,'Location','westoutside');
            end
            
            colormap parula
            %colormap copper
            %colormap pink
            %colormap gray
            %colormap hot
            %colormap cool
            %colormap jet
            %colormap turbo
            
            %c.Position(1) = c.Position(1) -0.005;
            if loopCount == 4
                xlabel('Memory strength', 'FontSize', 12, 'FontWeight', 'bold');
                ylabel('Meta-memory', 'FontSize', 12, 'FontWeight', 'bold');
            end
            
            %temp_title = title(sprintf('%s, %d trials',tempStr,length(x)),'FontSize',12,'FontWeight','bold');
            %temp_title = title(sprintf('%s, %d trials',tempStr,length(x)),'FontSize',12);
            %temp_title = title(sprintf(' %s \n %d trials ',tempStr,length(x)),'FontSize',12);
            temp_title = title(sprintf('%s',tempStr),'FontSize',12);
            temp_title.Interpreter = 'none';
        end
    end
    
    
    %% plot_trialDistriC
    if if_plot_trialDistriC == 1
        fig = figure('Name','Trial distribution C','NumberTitle','off'); %#ok<*NASGU>
        %set(gcf,'Position',[10 50 720 364*0.82]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %t = tiledlayout(1,2,'TileSpacing','tight','Padding','Compact');
        
        set(gcf,'Position',[10 50 720*0.5*0.98 364*0.82*1.05*0.99*1.005]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,1,'TileSpacing','tight','Padding','Compact');
        
        %for loopCount=1:5
        for loopCount=1:5
            %nexttile
            
            temp_if_smoothHistogram = if_smoothHistogram;
            
            if loopCount==1
                %nexttile(1)
                %nexttile(7)
                %tempStr = 'ChoiceMemoryError';
                tempStr = 'MemoryError';
                x_raw = memoryPrecision_trialLevel(choiceMemoryErrorBoolIndex);
                y_raw = meta_trialLevel(choiceMemoryErrorBoolIndex);
            elseif loopCount==2
                %nexttile(5)
                %nexttile(14)
                %tempStr = 'ChoiceMemoryCorrect';
                tempStr = 'MemoryCorrect';
                x_raw = memoryPrecision_trialLevel(choiceMemoryCorrectBoolIndex);
                y_raw = meta_trialLevel(choiceMemoryCorrectBoolIndex);
            elseif loopCount==3
                %nexttile(9)
                %nexttile(21)
                %tempStr = 'Choice offload';
                tempStr = 'Offload';
                x_raw = memoryPrecision_trialLevel(choiceOffloadBoolIndex);
                y_raw = meta_trialLevel(choiceOffloadBoolIndex);
            elseif loopCount==4
                %nexttile([3 3])
                %nexttile
                
                temp_if_smoothHistogram = 1;
                
                tempStr = 'Theoretical';
                %x_raw = memoryPrecision_trialLevel(choiceOffloadBoolIndex);
                %y_raw = meta_trialLevel(choiceOffloadBoolIndex);
                
                % Generate some theretical data point
                %tempN = 10000;
                tempN = sum(trialBoolIndex_allMatchMisMatch);
                
                temp1 = lowThreshold_memoryPrecision*lowThreshold_meta;
                temp2 = (1-lowThreshold_memoryPrecision)*(1-lowThreshold_meta);
                
                tempN1 = round(tempN * temp1 / (temp1+temp2));
                tempN2 = tempN - tempN1;
                
                % N1
                x_raw_N1 = nan(tempN1,1);
                y_raw_N1 = nan(tempN1,1);
                temp1A = round(lowThreshold_memoryPrecision*histogram_numBins);
                temp2A = round(lowThreshold_meta*histogram_numBins);
                tempStep = 1/histogram_numBins;
                tempRepNum = round(tempN1/(temp1A*temp2A));
                for tempi=1:temp1A
                    for tempj=1:temp2A
                        temptemp1 = (((tempi-1)*temp2A + tempj)-1)*tempRepNum+1;
                        temptemp2 = ((tempi-1)*temp2A + tempj)*tempRepNum;
                        temp_range = temptemp1:temptemp2;
                        
                        temp_x_value = tempStep * (tempi-0.5);
                        temp_y_value = tempStep * (tempj-0.5);
                        x_raw_N1(temp_range) = temp_x_value;
                        y_raw_N1(temp_range) = temp_y_value;
                    end
                end
                a = 1;
                
                % N2
                x_raw_N2 = nan(tempN2,1);
                y_raw_N2 = nan(tempN2,1);
                temp1A = histogram_numBins-round(lowThreshold_memoryPrecision*histogram_numBins);
                temp2A = histogram_numBins-round(lowThreshold_meta*histogram_numBins);
                tempStep = 1/histogram_numBins;
                %tempRepNum = round(tempN2/(temp1A*temp2A));
                %tempRepNum = tempRepNum;
                for tempi=1:temp1A
                    for tempj=1:temp2A
                        temptemp1 = (((tempi-1)*temp2A + tempj)-1)*tempRepNum+1;
                        temptemp2 = ((tempi-1)*temp2A + tempj)*tempRepNum;
                        temp_range = temptemp1:temptemp2;
                        
                        temp_x_value = tempStep * (tempi-0.5) + lowThreshold_memoryPrecision;
                        temp_y_value = tempStep * (tempj-0.5) + lowThreshold_meta;
                        x_raw_N2(temp_range) = temp_x_value;
                        y_raw_N2(temp_range) = temp_y_value;
                    end
                end
                a = 1;
                
                %
                x_raw = [x_raw_N1;x_raw_N2];
                y_raw = [y_raw_N1;y_raw_N2];
                
            elseif loopCount==5
                %nexttile([3 3])
                nexttile
                %tempStr = 'Data (merged)';
                %x_raw = memoryPrecision_trialLevel(choiceOffloadBoolIndex);
                %y_raw = meta_trialLevel(choiceOffloadBoolIndex);
                
                if if_trialDistriC_neuronLabel0_doubleLabel1 == 1
                    tempStr = 'Data (merged)';
                    temptempBoolIndex = trialBoolIndex_allMatchMisMatch;
                    x_raw = memoryPrecision_trialLevel(temptempBoolIndex);
                    y_raw = meta_trialLevel(temptempBoolIndex);
                elseif if_trialDistriC_neuronLabel0_doubleLabel1 == 0
                    %tempStr = 'Data (FreeChoice)';
                    tempStr = 'Joint distribution';
                    x_raw = memoryPrecision_trialLevel(choiceBoolIndex);
                    y_raw = meta_trialLevel(choiceBoolIndex);
                end
            end
            
            x = x_raw(~isnan(x_raw));
            y = y_raw(~isnan(x_raw));
            
            if temp_if_smoothHistogram == 0
                
            elseif temp_if_smoothHistogram == 1
                template_boolIndex = false(size(C));
                if loopCount == 1
                    template_boolIndex_multi = false(size(C,1),size(C,1),3);
                    template_boolIndex(temp_threshold_meta:end,1:temp_threshold_memoryPrecision) = true;
                    template_boolIndex_multi(:,:,1) = template_boolIndex;
                elseif loopCount == 2
                    template_boolIndex(temp_threshold_meta:end,temp_threshold_memoryPrecision:end) = true;
                    template_boolIndex_multi(:,:,2) = template_boolIndex;
                elseif loopCount == 3
                    template_boolIndex(1:temp_threshold_meta,:) = true;
                    template_boolIndex_multi(:,:,3) = template_boolIndex;
                end
            end
            
            %if loopCount < 4
            if loopCount < 5
                continue
            end
            
            
            if temp_if_smoothHistogram == 0
                XYBinLimits = [0 1];
                
                %if loopCount == 4
                %    h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','count');
                %    h.NumBins = [10 10];
                %    h.XBinLimits = XYBinLimits;
                %    h.YBinLimits = XYBinLimits;
                %    a = 1;
                %end
                %if loopCount == 5
                
                %if loopCount == 5
                %    temp_coeff = 2;%0.8,0.5
                %    x = x.^temp_coeff;
                %    y = y.^temp_coeff;
                %end
                
                h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','count');
                h.XBinLimits = XYBinLimits;
                h.YBinLimits = XYBinLimits;
                %h.NumBins = [10 10];
                h.NumBins = [1 1].* histogram_numBins;
                %end
                
                XYBinLimits = [-0.03 1.03];
                
                hold on
                
                %temp_max = 25;
                %a1 = h.Values;
                %a1(a1>temp_max) = temp_max;
                %h.Values = a1;
                
                caxis([0 15])
                
                %plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],XYBinLimits,...
                %    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);%[0.9290 0.6940 0.1250]
                %hold on
                %plot(XYBinLimits,[lowThreshold_meta lowThreshold_meta],...
                %    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);%[0.9290 0.6940 0.1250]
                %hold on
                
                plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],XYBinLimits,...
                    'LineWidth',3.5,'color',[1 1 1]);%[0.9290 0.6940 0.1250]
                hold on
                plot(XYBinLimits,[lowThreshold_meta lowThreshold_meta],...
                    'LineWidth',3.5,'color',[1 1 1]);%[0.9290 0.6940 0.1250]
                hold on
                
                
                xticks(0:0.2:1);
                yticks(0:0.2:1);
                
                xticklabels({'0','Low','','','High','1'});
                yticklabels({'0','Low','','','High','1'});
                
                xtickangle(0);
                ytickangle(0);
                
                
                
                
                temp_offset = 0.02;%0.02
                temp_offset2 = 0.015;
                temp_lineWidth = 4;%2.5
                
                temp_threshold_memoryPrecision = lowThreshold_memoryPrecision;
                temp_threshold_meta = lowThreshold_meta;
                
                if loopCount == 1 || loopCount == 4 || loopCount == 5
                    plot([1 1].*(temp_threshold_memoryPrecision-temp_offset),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                    plot([1 1].*(XYBinLimits(1)+temp_offset2),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                    plot([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(temp_threshold_meta+temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                    plot([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(XYBinLimits(2)-temp_offset2),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                end
                if loopCount == 2 || loopCount == 4 || loopCount == 5
                    plot([1 1].*(temp_threshold_memoryPrecision+temp_offset),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                    plot([1 1].*(XYBinLimits(2)-temp_offset2),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
                        [1 1].*(temp_threshold_meta+temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
                        [1 1].*(XYBinLimits(2)-temp_offset2),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                end
                if loopCount == 3 || loopCount == 4 || loopCount == 5
                    plot([1 1].*(temp_threshold_memoryPrecision-temp_offset),...
                        [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    plot([1 1].*(XYBinLimits(1)+temp_offset2),...
                        [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    plot([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(temp_threshold_meta-temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    plot([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(XYBinLimits(1)+temp_offset2),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    
                    
                    plot([1 1].*(temp_threshold_memoryPrecision+temp_offset),...
                        [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                    plot([1 1].*(XYBinLimits(2)-temp_offset2),...
                        [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
                        [1 1].*(temp_threshold_meta-temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
                        [1 1].*(XYBinLimits(1)+temp_offset2),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                end
                
                
                
                
                
            elseif temp_if_smoothHistogram == 1
                n=smoothCoeff;%50
                n=2^ceil(log2(n)); % round up n to the next power of 2;
                gridx1 = linspace(0,1,n);
                gridx2 = linspace(0,1,n);
                
                XYBinLimits = [1 n];
                
                [tempx1,tempx2] = meshgrid(gridx1, gridx2);
                tempx1 = tempx1(:);
                tempx2 = tempx2(:);
                xi = [tempx1 tempx2];
                [pdf1,xmesh1,bandwidth1] = ksdensity([x,y],xi,'Function','pdf','Support',[0 0;1 1],'BoundaryCorrection','Reflection');
                
                
                temp_threshold_meta = ceil(lowThreshold_meta*(n+1));
                temp_threshold_memoryPrecision = ceil((lowThreshold_memoryPrecision)*(n+1));
                
                
                
                pdf1_2d = reshape(pdf1,n,n);
                
                %if loopCount == 1
                %pdf_max = max(pdf1_2d,[],'all') * 1.1;
                
                % edit in 20240126
                %pdf_max = max(pdf1_2d,[],'all') * 0.8;
                %pdf_max = 2;
                %end
                
                if loopCount < 4
                    C = pdf1_2d;
                elseif loopCount == 4
                    C = zeros(size(pdf1_2d));
                    
                    template_boolIndex = false(size(C));
                    template_boolIndex(temp_threshold_meta:end,temp_threshold_memoryPrecision:end) = true;
                    template_boolIndex(1:temp_threshold_meta,1:temp_threshold_memoryPrecision) = true;
                    
                    if if_smoothHistogram == 0
                        pdf_max = 10;
                    end
                    C(template_boolIndex) = pdf_max;
                    
                elseif loopCount == 5
                    C = C_muilti(:,:,1) .* template_boolIndex_multi(:,:,1) + ...
                        C_muilti(:,:,2) .* template_boolIndex_multi(:,:,2) + ...
                        C_muilti(:,:,3) .* template_boolIndex_multi(:,:,3);
                end
                
                
                
                %imagesc([1 size(C,2)],[size(C,1) 1],C);
                %imagesc([1 size(C,2)],[size(C,1) 1],C,[0 5]);
                %imagesc([1 size(C,2)],[size(C,1) 1],C,[0 pdf_max]);
                imagesc([1 size(C,2)],[1 size(C,1)],C,[0 pdf_max]);
                hold on
                
                if loopCount == 1
                    C_muilti = nan(size(C,1),size(C,1),3);
                end
                if loopCount < 4
                    C_muilti(:,:,loopCount) = C;
                end
                
                
                % template
                
                %temp_threshold_meta = round((1-lowThreshold_meta)*(n+1));
                temp_threshold_meta = ceil(lowThreshold_meta*(n+1));
                temp_threshold_memoryPrecision = ceil((lowThreshold_memoryPrecision)*(n+1));
                
                if loopCount < 4
                    %C_template = zeros(size(C));
                    C_template = C;
                    
                    %                     template_boolIndex = false(size(C));
                    %
                    %                     if loopCount == 1
                    %                         template_boolIndex_multi = false(size(C,1),size(C,1),3);
                    %                         %C_muilti = nan(size(C,1),size(C,1),3);
                    %
                    %                         template_boolIndex(temp_threshold_meta:end,1:temp_threshold_memoryPrecision) = true;
                    %                         template_boolIndex_multi(:,:,1) = template_boolIndex;
                    %                         %C_muilti(:,:,1) = C;
                    %                     elseif loopCount == 2
                    %                         template_boolIndex(temp_threshold_meta:end,temp_threshold_memoryPrecision:end) = true;
                    %                         template_boolIndex_multi(:,:,2) = template_boolIndex;
                    %                         %C_muilti(:,:,1) = C;
                    %                     elseif loopCount == 3
                    %                         template_boolIndex(1:temp_threshold_meta,:) = true;
                    %                         template_boolIndex_multi(:,:,3) = template_boolIndex;
                    %                         %C_muilti(:,:,1) = C;
                    %                     end
                    C_template = C_template .* template_boolIndex;
                    
                    
                    %imagesc('XData',[1 size(C_template,2)],'YData',[size(C_template,1) 1],...
                    %    'CData',C_template,'AlphaData',0.5,[0 pdf_max]);
                    imagesc('XData',[1 size(C_template,2)],'YData',[1 size(C_template,1)],...
                        'CData',C_template,'AlphaData',0.5,[0 pdf_max]);
                    hold on
                end
                
                
                temp_offset = 3;%1-->3-->4
                temp_lineWidth = 4;%5-->4-->3
                
                if loopCount == 1 || loopCount == 4 || loopCount == 5
                    plot([1 1].*(temp_threshold_memoryPrecision-temp_offset),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                    plot([1 1].*(XYBinLimits(1)+temp_offset),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                    plot([1+temp_offset temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(temp_threshold_meta+temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                    plot([1+temp_offset temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(XYBinLimits(2)-temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                end
                if loopCount == 2 || loopCount == 4 || loopCount == 5
                    plot([1 1].*(temp_threshold_memoryPrecision+temp_offset),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                    plot([1 1].*(XYBinLimits(2)-temp_offset),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset],...
                        [1 1].*(temp_threshold_meta+temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset],...
                        [1 1].*(XYBinLimits(2)-temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                end
                if loopCount == 3 || loopCount == 4 || loopCount == 5
                    plot([1 1].*(temp_threshold_memoryPrecision-temp_offset),...
                        [XYBinLimits(1)+temp_offset temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    plot([1 1].*(XYBinLimits(1)+temp_offset),...
                        [XYBinLimits(1)+temp_offset temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    plot([1+temp_offset temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(temp_threshold_meta-temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    plot([1+temp_offset temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(XYBinLimits(1)+temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    
                    
                    plot([1 1].*(temp_threshold_memoryPrecision+temp_offset),...
                        [XYBinLimits(1)+temp_offset temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                    plot([1 1].*(XYBinLimits(2)-temp_offset),...
                        [XYBinLimits(1)+temp_offset temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset],...
                        [1 1].*(temp_threshold_meta-temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset],...
                        [1 1].*(XYBinLimits(1)+temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                end
                
                
                %plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision].*(n+1),XYBinLimits,...
                %    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
                %hold on
                %plot(XYBinLimits,[1-lowThreshold_meta 1-lowThreshold_meta].*(n+1),...
                %    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
                %hold on
                
                %                 if loopCount == 4
                %                     plot([1 1].*(temp_threshold_memoryPrecision),...
                %                         [XYBinLimits(1) XYBinLimits(2)],...
                %                         'LineWidth',temp_lineWidth,'color',[1 1 1]);
                %                     hold on
                %                     plot([XYBinLimits(1) XYBinLimits(2)],...
                %                         [temp_threshold_meta temp_threshold_meta],...
                %                         'LineWidth',temp_lineWidth,'color',[1 1 1]);
                %                     hold on
                %                 end
                
                
                if loopCount < 4
                    xticks(linspace(1,n,2));
                    yticks(linspace(1,n,2));
                    %xticklabels({'0','1'});
                    %yticklabels({'1','0'});
                    set(gca,'xticklabel','');
                    set(gca,'yticklabel','');
                else
                    xticks(linspace(1,n,6));
                    yticks(linspace(1,n,6));
                    %xticklabels({'0','0.2','0.4','0.6','0.8','1'});
                    %yticklabels({'0','0.2','0.4','0.6','0.8','1'});
                    xticklabels({'0','Low','','','High','1'});
                    yticklabels({'0','Low','','','High','1'});
                    
                    xtickangle(0);
                    ytickangle(0);
                    
                end
            end
            
            set(gca,'Ydir','normal') %reverse
            set(gca,'linewidth',1.5)
            %axis equal
            if temp_if_smoothHistogram == 1
                xlim(XYBinLimits);
                ylim(XYBinLimits);
            elseif temp_if_smoothHistogram == 0
                xlim([-0.03 1.03]);
                ylim([-0.03 1.03]);
            end
            set(gca, 'FontSize', 12) %12
            if loopCount == 5
                c = colorbar('FontSize',10);
                %c = colorbar('FontSize',10,'Location','southoutside');
                
                %temp_max = max(c.Ticks);
                %temp_max = 25;
                %c.Ticks = (0:5:temp_max);
                %c.Ticks = (0:5:temp_max).^temp_coeff;
                %c.TickLabels = {'0','0.2','0.4','0.6','0.8','1'};
                
                a = 1;
                %c.Limits(2) = 10;
            end
            
            if temp_if_smoothHistogram == 0
                if loopCount == 4
                    %c = colorbar('FontSize',10);
                    caxis([0 77])
                    %caxis([0 16])
                end
            end
            
            if if_colormap_loadEnhanced == 1
                load('parula_enhanced');
                colormap(parula_enhanced);
            elseif if_colormap_loadEnhanced == 0
                colormap parula
            end
            %colormap copper
            %colormap pink
            %colormap gray
            
            %c.Position(1) = c.Position(1) -0.005;
            if loopCount == 4 || loopCount == 5
                xlabel('Memory strength', 'FontSize', 12, 'FontWeight', 'bold');
                ylabel('Meta-memory', 'FontSize', 12, 'FontWeight', 'bold');
            end
            
            %temp_title = title(sprintf('%s, %d trials',tempStr,length(x)),'FontSize',12,'FontWeight','bold');
            %temp_title = title(sprintf('%s, %d trials',tempStr,length(x)),'FontSize',12);
            %temp_title = title(sprintf(' %s \n %d trials ',tempStr,length(x)),'FontSize',12);
            temp_title = title(sprintf('%s',tempStr),'FontSize',12);
            temp_title.Interpreter = 'none';
        end
    end
    
    %% plot_trialDistriD
    if if_plot_trialDistriD == 1
        fig = figure('Name','Trial distribution D','NumberTitle','off'); %#ok<*NASGU>
        set(gcf,'Position',[10 50 300 600]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        
        t = tiledlayout(3,1,'TileSpacing','loose','Padding','Compact');
        %t.Title.String = sprintf(sprintf('Trial distribution, %s',FOVName_currentFOV2));
        %t.Title.FontSize = 12;
        %t.Title.Interpreter = 'none';
        
        for loopCount=1:3
            %nexttile
            
            if loopCount==1
                nexttile
                tempStr = 'MemoryError';
                x_raw = memoryPrecision_trialLevel(choiceMemoryErrorBoolIndex);
                y_raw = meta_trialLevel(choiceMemoryErrorBoolIndex);
            elseif loopCount==2
                nexttile
                tempStr = 'MemoryCorrect';
                x_raw = memoryPrecision_trialLevel(choiceMemoryCorrectBoolIndex);
                y_raw = meta_trialLevel(choiceMemoryCorrectBoolIndex);
            elseif loopCount==3
                nexttile
                tempStr = 'Offload';
                x_raw = memoryPrecision_trialLevel(choiceOffloadBoolIndex);
                y_raw = meta_trialLevel(choiceOffloadBoolIndex);
            elseif loopCount==4
            end
            
            x = x_raw(~isnan(x_raw));
            y = y_raw(~isnan(x_raw));
            
            if if_smoothHistogram == 0
                XYBinLimits = [0 1];
                
                h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
                h.NumBins = [10 10];
                h.XBinLimits = XYBinLimits;
                h.YBinLimits = XYBinLimits;
                hold on
                
                plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],XYBinLimits,...
                    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
                hold on
                plot(XYBinLimits,[lowThreshold_meta lowThreshold_meta],...
                    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
                hold on
                
                xticks(0:0.2:1);
                yticks(0:0.2:1);
                
            elseif if_smoothHistogram == 1
                n=smoothCoeff;%50
                n=2^ceil(log2(n)); % round up n to the next power of 2;
                gridx1 = linspace(0,1,n);
                gridx2 = linspace(0,1,n);
                
                XYBinLimits = [1 n];
                
                [tempx1,tempx2] = meshgrid(gridx1, gridx2);
                tempx1 = tempx1(:);
                tempx2 = tempx2(:);
                xi = [tempx1 tempx2];
                [pdf1,xmesh1,bandwidth1] = ksdensity([x,y],xi,'Function','pdf','Support',[0 0;1 1],'BoundaryCorrection','Reflection');
                
                pdf1_2d = reshape(pdf1,n,n);
                
                if loopCount == 1
                    pdf_max = max(pdf1_2d,[],'all') * 1.1;
                end
                
                if loopCount < 4
                    C = pdf1_2d;
                elseif loopCount == 4
                    C = C_muilti(:,:,1) .* template_boolIndex_multi(:,:,1) + ...
                        C_muilti(:,:,2) .* template_boolIndex_multi(:,:,2) + ...
                        C_muilti(:,:,3) .* template_boolIndex_multi(:,:,3);
                end
                
                
                
                %imagesc([1 size(C,2)],[size(C,1) 1],C);
                %imagesc([1 size(C,2)],[size(C,1) 1],C,[0 5]);
                %imagesc([1 size(C,2)],[size(C,1) 1],C,[0 pdf_max]);
                imagesc([1 size(C,2)],[1 size(C,1)],C,[0 pdf_max]);
                hold on
                
                if loopCount == 1
                    C_muilti = nan(size(C,1),size(C,1),3);
                end
                if loopCount < 4
                    C_muilti(:,:,loopCount) = C;
                end
                
                
                % template
                
                %temp_threshold_meta = round((1-lowThreshold_meta)*(n+1));
                temp_threshold_meta = ceil(lowThreshold_meta*(n+1));
                temp_threshold_memoryPrecision = ceil((lowThreshold_memoryPrecision)*(n+1));
                
                if loopCount < 4
                    %C_template = zeros(size(C));
                    C_template = C;
                    
                    template_boolIndex = false(size(C));
                    
                    if loopCount == 1
                        template_boolIndex_multi = false(size(C,1),size(C,1),3);
                        %C_muilti = nan(size(C,1),size(C,1),3);
                        
                        template_boolIndex(temp_threshold_meta:end,1:temp_threshold_memoryPrecision) = true;
                        template_boolIndex_multi(:,:,1) = template_boolIndex;
                        %C_muilti(:,:,1) = C;
                    elseif loopCount == 2
                        template_boolIndex(temp_threshold_meta:end,temp_threshold_memoryPrecision:end) = true;
                        template_boolIndex_multi(:,:,2) = template_boolIndex;
                        %C_muilti(:,:,1) = C;
                    elseif loopCount == 3
                        template_boolIndex(1:temp_threshold_meta,:) = true;
                        template_boolIndex_multi(:,:,3) = template_boolIndex;
                        %C_muilti(:,:,1) = C;
                    end
                    C_template = C_template .* template_boolIndex;
                    
                    
                    %imagesc('XData',[1 size(C_template,2)],'YData',[size(C_template,1) 1],...
                    %    'CData',C_template,'AlphaData',0.5,[0 pdf_max]);
                    imagesc('XData',[1 size(C_template,2)],'YData',[1 size(C_template,1)],...
                        'CData',C_template,'AlphaData',0.5,[0 pdf_max]);
                    hold on
                end
                
                
                temp_offset = 3;%1
                temp_lineWidth = 4;%5
                
                if loopCount == 1 || loopCount == 4
                    plot([1 1].*(temp_threshold_memoryPrecision-temp_offset),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                    plot([1 1].*(XYBinLimits(1)+temp_offset),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                    plot([1+temp_offset temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(temp_threshold_meta+temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                    plot([1+temp_offset temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(XYBinLimits(2)-temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
                    hold on
                end
                if loopCount == 2 || loopCount == 4
                    plot([1 1].*(temp_threshold_memoryPrecision+temp_offset),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                    plot([1 1].*(XYBinLimits(2)-temp_offset),...
                        [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset],...
                        [1 1].*(temp_threshold_meta+temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset],...
                        [1 1].*(XYBinLimits(2)-temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
                    hold on
                end
                if loopCount == 3 || loopCount == 4
                    plot([1 1].*(temp_threshold_memoryPrecision-temp_offset),...
                        [XYBinLimits(1)+temp_offset temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    plot([1 1].*(XYBinLimits(1)+temp_offset),...
                        [XYBinLimits(1)+temp_offset temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    plot([1+temp_offset temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(temp_threshold_meta-temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    plot([1+temp_offset temp_threshold_memoryPrecision-temp_offset],...
                        [1 1].*(XYBinLimits(1)+temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
                    hold on
                    
                    
                    plot([1 1].*(temp_threshold_memoryPrecision+temp_offset),...
                        [XYBinLimits(1)+temp_offset temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                    plot([1 1].*(XYBinLimits(2)-temp_offset),...
                        [XYBinLimits(1)+temp_offset temp_threshold_meta-temp_offset],...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset],...
                        [1 1].*(temp_threshold_meta-temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                    plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset],...
                        [1 1].*(XYBinLimits(1)+temp_offset),...
                        'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
                    hold on
                end
                
                
                %plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision].*(n+1),XYBinLimits,...
                %    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
                %hold on
                %plot(XYBinLimits,[1-lowThreshold_meta 1-lowThreshold_meta].*(n+1),...
                %    'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
                %hold on
                
                
                if loopCount < 4
                    xticks(linspace(1,n,2));
                    yticks(linspace(1,n,2));
                    %xticklabels({'0','1'});
                    %yticklabels({'1','0'});
                    set(gca,'xticklabel','');
                    set(gca,'yticklabel','');
                else
                    xticks(linspace(1,n,6));
                    yticks(linspace(1,n,6));
                    xticklabels({'0','0.2','0.4','0.6','0.8','1'});
                    %yticklabels({'1','0.8','0.6','0.4','0.2','0'});
                    yticklabels({'0','0.2','0.4','0.6','0.8','1'});
                end
            end
            
            set(gca,'Ydir','normal') %reverse
            set(gca,'linewidth',1.5)
            axis equal
            xlim(XYBinLimits);
            ylim(XYBinLimits);
            set(gca, 'FontSize', 12) %12
            if loopCount == 4
                c = colorbar('FontSize',10);
                %c = colorbar('FontSize',10,'Location','westoutside');
            end
            
            colormap parula
            %colormap copper
            %colormap pink
            %colormap gray
            %colormap hot
            %colormap cool
            %colormap jet
            %colormap turbo
            
            %c.Position(1) = c.Position(1) -0.005;
            if loopCount == 4
                xlabel('Memory strength', 'FontSize', 12, 'FontWeight', 'bold');
                ylabel('Meta-memory', 'FontSize', 12, 'FontWeight', 'bold');
            end
            
            %temp_title = title(sprintf('%s, %d trials',tempStr,length(x)),'FontSize',12,'FontWeight','bold');
            %temp_title = title(sprintf('%s, %d trials',tempStr,length(x)),'FontSize',12);
            %temp_title = title(sprintf(' %s \n %d trials ',tempStr,length(x)),'FontSize',12);
            temp_title = title(sprintf('%s',tempStr),'FontSize',12);
            temp_title.Interpreter = 'none';
        end
    end
    
    %% plot_trialDistriE
    if if_plot_trialDistriE == 1
        fig = figure('Name','Trial distribution E','NumberTitle','off'); %#ok<*NASGU>
        
        set(gcf,'Position',[10 50 353*1.05 200*0.7*1.1*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
        
        for loopCount=1:3
            nexttile
            
            temp_if_smoothHistogram = if_smoothHistogram;
            
            if loopCount==1
                tempStr = 'Memory-correct';
                temptempTrialNum = bdac1(1);
            elseif loopCount==2
                tempStr = 'Offload';
                temptempTrialNum = bdac3(1);
            elseif loopCount==3
                tempStr = 'Memory-error';
                temptempTrialNum = bdac2(1);
            end
            
            
            XYBinLimits = [0 1];
            
            temp_offset = 0.02;%0.02
            temp_offset2 = 0.015;
            temp_lineWidth = 4;%2.5
            
            temp_threshold_memoryPrecision = lowThreshold_memoryPrecision;
            temp_threshold_meta = lowThreshold_meta;
            
            plot([1 1].*(temp_threshold_memoryPrecision-temp_offset),...
                [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
            hold on
            plot([1 1].*(XYBinLimits(1)+temp_offset2),...
                [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
            hold on
            plot([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
                [1 1].*(temp_threshold_meta+temp_offset),...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
            hold on
            plot([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
                [1 1].*(XYBinLimits(2)-temp_offset2),...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
            hold on
            
            
            plot([1 1].*(temp_threshold_memoryPrecision+temp_offset),...
                [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
            hold on
            plot([1 1].*(XYBinLimits(2)-temp_offset2),...
                [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
            hold on
            plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
                [1 1].*(temp_threshold_meta+temp_offset),...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
            hold on
            plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
                [1 1].*(XYBinLimits(2)-temp_offset2),...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
            hold on
            
            
            
            
            plot([1 1].*(temp_threshold_memoryPrecision-temp_offset),...
                [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
            hold on
            plot([1 1].*(XYBinLimits(1)+temp_offset2),...
                [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
            hold on
            plot([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
                [1 1].*(temp_threshold_meta-temp_offset),...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
            hold on
            plot([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
                [1 1].*(XYBinLimits(1)+temp_offset2),...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
            hold on
            
            
            plot([1 1].*(temp_threshold_memoryPrecision+temp_offset),...
                [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
            hold on
            plot([1 1].*(XYBinLimits(2)-temp_offset2),...
                [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
            hold on
            plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
                [1 1].*(temp_threshold_meta-temp_offset),...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
            hold on
            plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
                [1 1].*(XYBinLimits(1)+temp_offset2),...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
            hold on
            
            
            set(gca,'color',[0.8 0.8 0.8]);
            set(gca,'Ydir','normal') %reverse
            set(gca,'linewidth',1.5)
            xticks([]);
            yticks([]);
            
            temp_title = title(sprintf('%s\n%d trials',tempStr,temptempTrialNum),'FontSize',8);
            temp_title.Interpreter = 'none';
        end
    end
    
    
    %% plot_trialDistriE2
    if if_plot_trialDistriE2 == 1
        fig = figure('Name','Trial distribution E2','NumberTitle','off'); %#ok<*NASGU>
        
        %set(gcf,'Position',[10 290 353*1.05 200*0.7*1.1*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 290 125 600*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(4,1,'TileSpacing','Compact','Padding','Compact');
        
        for loopCount=1:4
            nexttile
            
            temp_if_smoothHistogram = if_smoothHistogram;
            
            if loopCount==1
                tempStr = 'Choice-memory correct';
                temptempTrialNum = bdac1(1);
            elseif loopCount==2
                tempStr = 'Choice-offload';
                temptempTrialNum = bdac3(1);
            elseif loopCount==3
                tempStr = 'Choice-memory error';
                temptempTrialNum = bdac2(1);
            elseif loopCount==4
                tempStr = 'Choice-memory error';
                temptempTrialNum = bdac2(1);
            end
            
            
            XYBinLimits = [0 1];
            
            temp_offset = 0.02;%0.02
            temp_offset2 = 0.015;
            temp_lineWidth = 4;%2.5
            
            temp_threshold_memoryPrecision = lowThreshold_memoryPrecision;
            temp_threshold_meta = lowThreshold_meta;
            
            plot([1 1].*(temp_threshold_memoryPrecision-temp_offset),...
                [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
            hold on
            plot([1 1].*(XYBinLimits(1)+temp_offset2),...
                [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
            hold on
            plot([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
                [1 1].*(temp_threshold_meta+temp_offset),...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
            hold on
            plot([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
                [1 1].*(XYBinLimits(2)-temp_offset2),...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
            hold on
            
            
            plot([1 1].*(temp_threshold_memoryPrecision+temp_offset),...
                [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
            hold on
            plot([1 1].*(XYBinLimits(2)-temp_offset2),...
                [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
            hold on
            plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
                [1 1].*(temp_threshold_meta+temp_offset),...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
            hold on
            plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
                [1 1].*(XYBinLimits(2)-temp_offset2),...
                'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
            hold on
            
            
            
            
            plot([1 1].*(temp_threshold_memoryPrecision-temp_offset),...
                [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
            hold on
            plot([1 1].*(XYBinLimits(1)+temp_offset2),...
                [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
            hold on
            plot([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
                [1 1].*(temp_threshold_meta-temp_offset),...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
            hold on
            plot([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
                [1 1].*(XYBinLimits(1)+temp_offset2),...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
            hold on
            
            
            plot([1 1].*(temp_threshold_memoryPrecision+temp_offset),...
                [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
            hold on
            plot([1 1].*(XYBinLimits(2)-temp_offset2),...
                [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
            hold on
            plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
                [1 1].*(temp_threshold_meta-temp_offset),...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
            hold on
            plot([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
                [1 1].*(XYBinLimits(1)+temp_offset2),...
                'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
            hold on
            
            
            set(gca,'color',[0.8 0.8 0.8]);
            set(gca,'Ydir','normal') %reverse
            set(gca,'linewidth',1.5)
            xticks([]);
            yticks([]);
            
            temp_title = title(sprintf('%s\n%d trials',tempStr,temptempTrialNum),'FontSize',8);
            temp_title.Interpreter = 'none';
        end
    end
    
    %% Plot seq count distribution of mismatch trials
    if if_plot_seqDistri == 1
        fig = figure('Name','Seq distribution','NumberTitle','off'); %#ok<*NASGU>
        set(gcf,'Position',[10 50 1200 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        if if_plot_fineTuningMismatch == 0
            y1 = seqCount_memoryPrecisionLowError_choiceMemory;
            y2 = seqCount_memoryPrecisionHigh_choiceOffload;
        elseif if_plot_fineTuningMismatch == 1
            y1 = seqCount_memoryPrecisionLowError_choiceMemoryHigh;
            y2 = seqCount_memoryPrecisionHigh_choiceOffloadLow;
        end
        
        [y_min,y_max] = bounds([y1;y2]);
        
        plot(y1,'color',[0.4660 0.6740 0.1880],'LineWidth',1.5);
        hold on
        plot(y2,'color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
        hold on
        
        if if_plot_fineTuningMismatch == 0
            %le = legend('memoryPrecisionLowError_choiceMemory','memoryPrecisionHigh_choiceOffload',...
            %    'Location','northeast','fontsize',11);
            
            temptempStr1 = sprintf('memoryPrecisionLowError_choiceMemory(%d)',sum(y1));
            temptempStr2 = sprintf('memoryPrecisionHigh_choiceOffload(%d)',sum(y2));
        elseif if_plot_fineTuningMismatch == 1
            %le = legend('memoryPrecisionLowError_choiceMemoryHigh','memoryPrecisionHigh_choiceOffloadLow',...
            %    'Location','northeast','fontsize',11);
            
            temptempStr1 = sprintf('memoryPrecisionLowError_choiceMemoryHigh(%d)',sum(y1));
            temptempStr2 = sprintf('memoryPrecisionHigh_choiceOffloadLow(%d)',sum(y2));
        end
        le = legend(temptempStr1,temptempStr2,'Location','northeast','fontsize',11);
        le.Interpreter = 'none';
        
        
        xlim([0 sum(numSeq(1:3))+1]);
        ylim([y_min y_max+(y_max-y_min)*0.4]);%0.1
        set(gca, 'FontSize', 16)
        set(gca,'box','off');% 取消右、上边框
        
        
        set(gca,'XTick',1:sum(numSeq(1:3)));
        xtl=string(seqSet_inOne_inOne(1:sum(numSeq(1:3))));
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp = zeros(size(xtext_xp)) - (y_max-y_min)*0.04;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',90,'fontsize',15);%9
        set(gca,'xticklabel','');
        
        ylabel('Trial count', 'FontSize', 14, 'FontWeight', 'bold');
        
        if if_memeoryPrecision_stimuli0_response1 == 0
            if if_plot_fineTuningMismatch == 0
                temp_str = 'stimuli, roughTuning';
            elseif if_plot_fineTuningMismatch == 1
                temp_str = 'stimuli, fineTuning';
            end
        elseif if_memeoryPrecision_stimuli0_response1 == 1
            if if_plot_fineTuningMismatch == 0
                temp_str = 'response, roughTuning';
            elseif if_plot_fineTuningMismatch == 1
                temp_str = 'response, fineTuning';
            end
        end
        
        
        temp_title = title(sprintf('Seq distribution, %s, %s',temp_str,FOVName_currentFOV2),...
            'FontSize',14,'FontWeight','bold');
        temp_title.Interpreter = 'none';
        
    end
    
    
    
    %     %% Compare ChoiceMemory trial accuracy in LowPrecision and HighPrecision
    %
    %     fig = figure('Name','','NumberTitle','off'); %#ok<*NASGU>
    %     %set(gcf,'Position',[10 50 400 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %     %set(gcf,'Position',[35+0 40+0 273*0.8 2*233*0.865]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %     set(gcf,'Position',[35+0 40+0 152*0.85 403]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %
    %     t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    %     nexttile
    %
    %     temp_p = temp_p_LowPrecision_HighPrecision;
    %     temp_1 = correctBoolIndex_LowPrecision;
    %     temp_2 = correctBoolIndex_HighPrecision;
    %
    %
    %     temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    %     temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    %
    %
    %     temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM]);
    %     temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
    %
    %         temp_bar = bar([1 2], [mean(temp_1) mean(temp_2)], ...
    %             'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1);
    %         hold on
    %         a = 1;
    %         temp_bar.CData(1,:) = [1 1 1]*0.5;%0.3
    %         temp_bar.CData(2,:) = [1 1 1]*0.5;
    %
    %         errorbar([1 2], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 8);
    %
    %     tempTxt = sprintf('');
    %     if temp_p < 0.001
    %         tempTxt = sprintf('***');
    %     elseif temp_p < 0.01
    %         tempTxt = sprintf('**');
    %     elseif temp_p < 0.05
    %         tempTxt = sprintf('*');
    %     end
    %     text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %         'HorizontalAlignment','center');
    %
    %     set(gca,'linewidth',1.5)
    %     xlim([0.15 2.65])
    %     ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
    %     set(gca, 'FontSize', 10)
    %
    %     %set(gca,'XTickLabel', ["HighMatch"; "OverMismatch"],'FontSize', 12);%给坐标加标签
    %     %xtickangle(20);
    %
    %     yticks([0.7 0.8]);
    %
    %     %xtl = ["HighMatch"; "OverMismatch"];
    %     %xtl = ["OverMismatch"; "HighMatch"];
    %     %xtl = ["LowPrecison"; "HighPrecision"];
    %     xtl = ["Low"; "High"];
    %     xt=get(gca,'XTick');
    %     yt=get(gca,'YTick');
    %     xtext_xp=xt;
    %     xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.57;%0.56,0.4
    %     text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%25
    %     set(gca,'xticklabel','');
    %
    %
    %     set(gca,'box','off');% 取消右、上边框
    %     ylabel('Accuracy', 'FontSize', 11, 'FontWeight', 'bold');
    %     %temp_title = title(sprintf('Precision\nChoiceMemory'),'fontsize',9);
    %     title(sprintf('Precision'),'fontsize',10);
    %     subtitle(sprintf('ChoiceMemory'),'fontsize',9);
    %     a = 1;
    %
    %
    %
    %     %% Compare FreeChoice trial rProb in LowMeta and HighMeta
    %
    %     fig = figure('Name','','NumberTitle','off'); %#ok<*NASGU>
    %     %set(gcf,'Position',[35+0 40+0 273*0.8 2*233*0.865]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %     set(gcf,'Position',[35+0 40+0 152*0.85 403]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %
    %
    %     t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    %     nexttile
    %
    %     temp_p = temp_p_metaLow_metaHigh;
    %     temp_1 = choiceOffloadBoolIndex_metaLow;
    %     temp_2 = choiceOffloadBoolIndex_metaHigh;
    %
    %
    %     temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    %     temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    %
    %
    %     temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM]);
    %     temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
    %
    %     temp_bar = bar([1 2], [mean(temp_1) mean(temp_2)], ...
    %         'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1);
    %     hold on
    %     a = 1;
    %     temp_bar.CData(1,:) = [1 1 1]*0.5;%0.3
    %     temp_bar.CData(2,:) = [1 1 1]*0.5;
    %
    %     errorbar([1 2], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 8);
    %
    %     tempTxt = sprintf('');
    %     if temp_p < 0.001
    %         tempTxt = sprintf('***');
    %     elseif temp_p < 0.01
    %         tempTxt = sprintf('**');
    %     elseif temp_p < 0.05
    %         tempTxt = sprintf('*');
    %     end
    %     text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %         'HorizontalAlignment','center');
    %
    %     set(gca,'linewidth',1.5)
    %     xlim([0.15 2.65])
    %     ylim([max(temp_y_min-(temp_y_max-temp_y_min)*0.25,0) temp_y_max+(temp_y_max-temp_y_min)*0.25]);
    %     set(gca, 'FontSize', 10)
    %
    %
    %     %xtl = ["LowMeta"; "HighMeta"];
    %     xtl = ["Low"; "High"];
    %     xt=get(gca,'XTick');
    %     yt=get(gca,'YTick');
    %     xtext_xp=xt;
    %     %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.555;%0.555
    %     xtext_yp=ones(1,length(xt))*temp_y_min-(temp_y_max-temp_y_min)*0.3;%
    %     text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    %     set(gca,'xticklabel','');
    %
    %
    %     set(gca,'box','off');% 取消右、上边框
    %     ylabel('Offloading rate', 'FontSize', 11, 'FontWeight', 'bold');
    %     %title(sprintf('FreeChoice'),'fontsize',12);
    %     title(sprintf('Meta-memory'),'fontsize',10);
    %     subtitle(sprintf('FreeChoice'),'fontsize',9);
    
    %% plot example seq distribution
    
    
    %exampleSeq = 16;% Accuracy: 17,13,16
    %exampleSeq = 10;% Offloading rate: 17,10,21
    exampleSeq = 18;% 16,18,21,23,26
    
    fprintf('exampleSeq is %d.\n',seqSet_inOne_inOne(exampleSeq));
    
    mildSeqBoolIndex_seqLevel = false(1,length(offloadingProb_inOne));
    mildSeqBoolIndex_seqLevel(exampleSeq) = true;
    
    mildSeqBoolIndex = false(size(choiceBoolIndex));
    for tempi=1:length(mildSeqBoolIndex)
        if mildSeqBoolIndex_seqLevel(seqIndex(tempi)) == true
            mildSeqBoolIndex(tempi) = true;
        end
    end
    
    ChoiceBoolIndex_mildSeq = mildSeqBoolIndex & choiceBoolIndex;
    
    
    if if_plot_exampleSeqDistri == 1
        
        %% 1d distrition
        fig = figure('Name','example seq distribution','NumberTitle','off');
        %set(gcf,'Position',[380 450 300*0.85*1.05*1.03 200*0.7*0.85*1.03]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[380 450 276*1.21 126*1.21]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,2,'TileSpacing','compact','Padding','loose');
        
        % memoryPrecision
        nexttile
        
        x1 = memoryPrecision_trialLevel(ChoiceBoolIndex_mildSeq);
        
        
        h_NumBins = 8;%10
        
        x = x1;
        h1 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
        hold on
        h1.NumBins = h_NumBins;
        h1.EdgeColor = [0 0 0];
        
        y1 = h1.Values;
        
        x_min = 0;
        x_max = 1;
        
        [y_min,y_max] = bounds(y1);
        
        if if_plot_additionalSmooth == 1
            n=100;
            n=2^ceil(log2(n)); % round up n to the next power of 2;
            
            [pdf1,xmesh1,~] = ksdensity(x1,'NumPoints',n,'Function','pdf');
            plot(xmesh1,pdf1*sum(y1)*h1.BinWidth,':','LineWidth',1.5,'color',[0 0 0]);
            hold on
        end
        
        %         plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],[y_min y_max],...
        %             'LineWidth',1,'color',[0.9290 0.6940 0.1250]);
        %         hold on
        
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
        ylim([y_min y_max+(y_max-y_min)*0.1]);%0.1
        %xticks([0 1]);
        set(gca, 'FontSize', 9)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Memory strength', 'FontSize', 10, 'FontWeight', 'normal');
        ylabel('Trial count', 'FontSize', 10, 'FontWeight', 'normal');
        
        
        % meta_trialLevel
        nexttile
        
        x1 = meta_trialLevel(ChoiceBoolIndex_mildSeq);
        
        
        h_NumBins = 8;%10
        
        x = x1;
        h1 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
        hold on
        h1.NumBins = h_NumBins;
        h1.EdgeColor = [0 0 0];
        
        y1 = h1.Values;
        
        x_min = 0;
        x_max = 1;
        
        [y_min,y_max] = bounds(y1);
        
        if if_plot_additionalSmooth == 1
            n=100;
            n=2^ceil(log2(n)); % round up n to the next power of 2;
            
            [pdf1,xmesh1,~] = ksdensity(x1,'NumPoints',n,'Function','pdf');
            plot(xmesh1,pdf1*sum(y1)*h1.BinWidth,':','LineWidth',1.5,'color',[0 0 0]);
            hold on
        end
        
        %         plot([lowThreshold_meta lowThreshold_meta],[y_min y_max],...
        %             'LineWidth',1,'color',[0.9290 0.6940 0.1250]);
        %         hold on
        
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
        ylim([y_min y_max+(y_max-y_min)*0.1]);%0.1
        %xticks([0 1]);
        set(gca, 'FontSize', 9)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta-memory', 'FontSize', 10, 'FontWeight', 'normal');
        ylabel('Trial count', 'FontSize', 10, 'FontWeight', 'normal');
        
        
        %% 2d distrition
        fig = figure('Name','example seq distribution','NumberTitle','off');
        %set(gcf,'Position',[680 450 150 200*0.7*0.85*1.03]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[680 450 150*1.18 200*0.7*0.85*1.03]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[680 450 150*1.18*1.21 200*0.7*0.85*1.03*1.21]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,5,'TileSpacing','compact','Padding','loose');
        
        %temp_histogram_numBins = histogram_numBins;
        temp_histogram_numBins = 10;
        
        nexttile([1 4])
        
        tempStr = 'FreeChoice';
        %x_raw = memoryPrecision_trialLevel(choiceBoolIndex);
        %y_raw = meta_trialLevel(choiceBoolIndex);
        x_raw = memoryPrecision_trialLevel(ChoiceBoolIndex_mildSeq);
        y_raw = meta_trialLevel(ChoiceBoolIndex_mildSeq);
        
        
        x = x_raw(~isnan(x_raw));
        y = y_raw(~isnan(x_raw));
        
        [x_min,x_max] = bounds(x);
        x_max = x_max + 0.01;
        
        XYBinLimits = [0 1];
        
        %h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
        h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','count');
        %h.NumBins = [10 10];
        h.XBinLimits = XYBinLimits;
        h.YBinLimits = XYBinLimits;
        h.NumBins = [1 1].*temp_histogram_numBins;
        hold on
        
        %         plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],XYBinLimits,...
        %             'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
        %         hold on
        %         plot(XYBinLimits,[lowThreshold_meta lowThreshold_meta],...
        %             'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
        %         hold on
        
        %xticks([0 1]);
        %yticks([0 1]);
        xticks('');
        yticks('');
        
        
        %axis equal
        set(gca,'Ydir','normal') %reverse
        xlim(XYBinLimits);
        ylim(XYBinLimits);
        xtickangle(0);
        set(gca, 'FontSize', 9) %12
        
        %c = colorbar('FontSize',7);
        c = colorbar('eastoutside','FontSize',7);
        %c.Position = c.Position+[0.05 0.15 -0.04 -0.15];
        c.Position = c.Position+[0.05 0.10 -0.02 -0.10];
        
        
        xlabel('Memory strength', 'FontSize', 10, 'FontWeight', 'normal');
        ylabel('Meta-memory', 'FontSize', 10, 'FontWeight', 'normal');
        
        
    end
    
    
    %% Plot trial proportion of choice-memory correct, choice-offload, choice-memory error
    if true
        fig = figure('Name','meta_seqLevel','NumberTitle','off');
        %set(gcf,'Position',[610 390 337*0.8*0.99*1.15 176*1.065*1.2*0.85]);
        set(gcf,'Position',[610 390 337*0.8*0.99*1.15*1.1 176*1.065*1.2*0.85]);
        t = tiledlayout(4,1,'TileSpacing','Compact','Padding','Compact');
        
        %[r_bdac123,p_bdac123] = corr([bdac1_n11n_B,bdac2_n11n_B,bdac3_n11n_B]',...
        %    [bdac1_template,bdac2_template,bdac3_template]');
        
        x1 = bdac1_n11n_B;
        x2 = bdac3_n11n_B;
        x3 = bdac2_n11n_B;
        
        y1 = bdac1_template;
        y2 = bdac3_template;
        y3 = bdac2_template;
        
        temp_r = r_bdac123;
        temp_p = p_bdac123;
        
        
        nexttile([3 1])
        
        
        tempStr = 'Correlation';
        
        x = [x1,x2,x3];
        y = [y1,y2,y3];
        
        h = [];
        h = [scatter(-1,-1,60,'x','MarkerEdgeColor',[0 0 0],'LineWidth',2) h];
        hold on
        h = [scatter(-1,-1,60,'s','MarkerEdgeColor',[0 0 0],'LineWidth',2) h];
        hold on
        h = [scatter(-1,-1,60,'o','MarkerEdgeColor',[0 0 0],'LineWidth',2) h];
        hold on
        
        %         h = [plot(-1,-1,'Marker','o','MarkerEdgeColor',[0 0 0],'LineWidth',2) ;h];
        %         hold on
        %         h = [plot(-1,-1,'Marker','s','MarkerEdgeColor',[0 0 0],'LineWidth',2) ;h];
        %         hold on
        %         h = [plot(-1,-1,'Marker','x','MarkerEdgeColor',[0 0 0],'LineWidth',2) ;h];
        %         hold on
        
        for tempi=1:3
            if tempi==1
                temp_x = x1;
                temp_y = y1;
                temp_mkr = 'o';
            elseif tempi==2
                temp_x = x2;
                temp_y = y2;
                temp_mkr = 's';
            elseif tempi==3
                temp_x = x3;
                temp_y = y3;
                temp_mkr = 'x';
            end
            
            for tempj=1:4
                if tempj==1
                    temp_color = color_choiceMemoryHigh;
                elseif tempj==2
                    temp_color = color_choiceOffloadHigh;
                elseif tempj==3
                    temp_color = color_choiceOffloadLow;
                elseif tempj==4
                    temp_color = color_choiceMemoryLow;
                end
                %scatter(temp_x(tempj),temp_y(tempj),...
                %    30,'filled','MarkerFaceColor',temp_color,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.7);
                %scatter(temp_x(tempj),temp_y(tempj),...
                %   60,temp_mkr,'filled','MarkerFaceColor',temp_color,'MarkerEdgeColor',temp_color,'LineWidth',2);
                scatter(temp_x(tempj),temp_y(tempj),...
                    60,temp_mkr,'MarkerEdgeColor',temp_color,'LineWidth',2);
                hold on
            end
        end
        
        [temp_xmin,temp_xmax] = bounds(x);
        [temp_ymin,temp_ymax] = bounds(y);
        
        
        xlim([temp_xmin-(temp_xmax-temp_xmin)*0.08 temp_xmax+(temp_xmax-temp_xmin)*0.08]);
        ylim([temp_ymin-(temp_ymax-temp_ymin)*0.08 temp_ymax+(temp_ymax-temp_ymin)*0.08]);
        
        
        tempTxt = 'p>0.05';
        if temp_p < 0.05
            tempTxt = 'p<0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p<0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p<0.001';
        end
        text(temp_xmin+(temp_xmax-temp_xmin)*0.03,temp_ymin+(temp_ymax-temp_xmin)*0.99,sprintf('r=%.3f',temp_r), 'fontsize',9,'FontWeight','normal');
        text(temp_xmin+(temp_xmax-temp_xmin)*0.03,temp_ymin+(temp_ymax-temp_xmin)*0.87,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
        
        %le = legend(h,'Choice-memory (Correct)','Choice-offload','Choice-memory (Error)',...
        %    'Location','southeast','fontsize',7);
        %le.ItemTokenSize = ones(1,3)*8;
        
        le = legend(h,'Choice-memory(Correct)','Choice-offload','Choice-memory(Error)',...
            'Location','southeast','fontsize',6.5,'NumColumns',3);
        le.ItemTokenSize = ones(1,3)*8;
        
        
        
        %xticks([]);
        %yticks([]);
        set(gca,'linewidth',1.5)
        set(gca, 'FontSize', 10)
        set(gca,'box','off');% 取消右、上边框
        xlabel(sprintf('Proportion (Data)'), 'FontSize', 10, 'FontWeight', 'normal');
        ylabel(sprintf('Proportion (Template)'), 'FontSize', 10, 'FontWeight', 'normal');
        
        %temp_title = title(sprintf('Correlation'), 'FontSize', 10, 'FontWeight', 'normal');
        %temp_title.Interpreter = 'none';
        
    end
    
    
    
    %%
    if if_plot_trialDistriF
        %% 2d distrition
        fig = figure('Name','example seq distribution','NumberTitle','off');
        %set(gcf,'Position',[680 450 150*1.18*1.21*0.9*0.95*1.1 200*0.7*0.85*1.03*1.21*1.1*1.1]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[680 450 201.4*1.08*0.95*0.98*0.97 179.5*1.08*0.97]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,5,'TileSpacing','compact','Padding','loose');
        
        temp_histogram_numBins = histogram_numBins;
        %temp_histogram_numBins = 10;
        
        nexttile([1 4])
        
        tempStr = 'FreeChoice';
        %x_raw = memoryPrecision_trialLevel(choiceBoolIndex);
        %y_raw = meta_trialLevel(choiceBoolIndex);
        %x_raw = memoryPrecision_trialLevel(ChoiceBoolIndex_mildSeq);
        %y_raw = meta_trialLevel(ChoiceBoolIndex_mildSeq);
        x_raw = memoryPrecision_trialLevel(choiceMemoryErrorBoolIndex);
        y_raw = meta_trialLevel(choiceMemoryErrorBoolIndex);
        
        
        x = x_raw(~isnan(x_raw));
        y = y_raw(~isnan(x_raw));
        
        [x_min,x_max] = bounds(x);
        x_max = x_max + 0.01;
        
        XYBinLimits = [0 1];
        
        %h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
        %h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','count');
        h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','off','Normalization','count');
        %h.NumBins = [10 10];
        h.XBinLimits = XYBinLimits;
        h.YBinLimits = XYBinLimits;
        h.NumBins = [1 1].*temp_histogram_numBins;
        hold on
        
        %         plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],XYBinLimits,...
        %             'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
        %         hold on
        %         plot(XYBinLimits,[lowThreshold_meta lowThreshold_meta],...
        %             'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
        %         hold on
        
        xticks([0 1]);
        yticks([0 1]);
        %xticks('');
        %yticks('');
        
        
        %axis equal
        set(gca,'linewidth',1.5)
        set(gca,'Ydir','normal') %reverse
        xlim(XYBinLimits);
        ylim(XYBinLimits);
        xtickangle(0);
        set(gca, 'FontSize', 8) %12
        
        %c = colorbar('FontSize',7);
        c = colorbar('eastoutside','FontSize',8);
        %c.Position = c.Position+[0.05 0.15 -0.04 -0.15];
        c.Position = c.Position+[0.03 0.10 -0.02 -0.13];
        
        if if_colormap_loadEnhanced == 1
            %load('parula_enhanced');
            %colormap(parula_enhanced);
            %load('gray_enhanced');
            %colormap(gray_enhanced);
            
            %load('gray_reversed_enhanced');
            %colormap(gray_reversed_enhanced);
            
            load('coolwarm_enhanced');
            colormap(coolwarm_enhanced);
            
            %load('coolwarm_n11n_enhanced');
            %colormap(coolwarm_n11n_enhanced);
            
        elseif if_colormap_loadEnhanced == 0
            %colormap parula
            %colormap gray
            
            temp1 = gray;
            temp1 = temp1(end:-1:1,:);
            colormap(temp1);
            
            %colormap(coolwarm());
            
            %temp1 = coolwarm(300);
            %temp2 = ((300-256)/2)+1;
            %temp3 = temp1(temp2:temp2+255,:);
            %colormap(temp3);
        end
        
        
        xlabel('Memory strength', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
    end
    
end

%% (Free-choice) Trial-level correlation between precision and meta-memory
memoryPrecision_trialLevel;
meta_trialLevel;

temptempBoolIndex_A = choiceBoolIndex';
temptempBoolIndex_B = (~isnan(memoryPrecision_trialLevel)) & (~isnan(meta_trialLevel));

temptempBoolIndex = temptempBoolIndex_A & temptempBoolIndex_B;

[temp_r,temp_p] = corr(memoryPrecision_trialLevel(temptempBoolIndex),meta_trialLevel(temptempBoolIndex));

if true
    fig = figure('Name','Free-choice trial correlation','NumberTitle','off');
    %set(gcf,'Position',[610 390 337*0.8*0.99*1.15 176*1.065*1.2*0.85]);
    %set(gcf,'Position',[610 390 337*0.8*0.99*1.15*0.64 176*1.065*1.2*0.85*0.85*1.4]);
    %set(gcf,'Position',[610 390 196.4*0.905 227.5*0.905]);
    %set(gcf,'Position',[610 390 196.4*0.905*1.11*0.97 227.5*0.905*1.11*0.97]);
    %set(gcf,'Position',[610 390 196.4*0.905*1.11*0.97 227.5*0.905*1.11*0.97*0.895]);
    set(gcf,'Position',[610 390 196.4*0.905*1.11*0.97*0.85 227.5*0.905*1.11*0.97*0.895]);
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    x = memoryPrecision_trialLevel(temptempBoolIndex);
    y = meta_trialLevel(temptempBoolIndex);
    
    
    [r,p] = corr(x,y);
    
    mdl = fitglm(x,y);
    
    
    temp_size = 3;%10
    temp_h = scatter(x, y, ...
        temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y);
    
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    
    %xlim([temp_xmin-(temp_xmax-temp_xmin)*0.08 temp_xmax+(temp_xmax-temp_xmin)*0.08]);
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    
    xticks([0 1]);
    yticks([0 1]);
    
    if if_monkey_D0_Z1 == 1
        xlim([-0.05 0.50]);
        %xticks([0 0.4]);
        xticks([0 0.45]);        
    end
    
    
    
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    %text(0.76,0.25,sprintf('r=%.3f',r), 'fontsize',9,'FontWeight','normal');
    %text(0.76,0.12,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
    
    %title(sprintf('Free-choice'), 'FontSize', 9, 'FontWeight', 'normal');
    title(sprintf('Across choice trials'), 'FontSize', 9, 'FontWeight', 'normal');
    subtitle(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
    
end


%% (Choice-memory) Trial-level correlation between precision and meta-memory
temptempBoolIndex_A = choiceMemoryBoolIndex';
temptempBoolIndex_B = (~isnan(memoryPrecision_trialLevel)) & (~isnan(meta_trialLevel));

temptempBoolIndex = temptempBoolIndex_A & temptempBoolIndex_B;

[temp_r,temp_p] = corr(memoryPrecision_trialLevel(temptempBoolIndex),meta_trialLevel(temptempBoolIndex));

if true
    fig = figure('Name','Choice-memory trial correlation','NumberTitle','off');
    set(gcf,'Position',[610+200 390 196.4*0.905*1.11*0.97 227.5*0.905*1.11*0.97]);
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    x = memoryPrecision_trialLevel(temptempBoolIndex);
    y = meta_trialLevel(temptempBoolIndex);
    
    
    [r,p] = corr(x,y);
    
    mdl = fitglm(x,y);
    
    
    temp_size = 3;%10
    temp_h = scatter(x, y, ...
        temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y);
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    
    %xlim([temp_xmin-(temp_xmax-temp_xmin)*0.08 temp_xmax+(temp_xmax-temp_xmin)*0.08]);
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    
    xticks([0 1]);
    yticks([0 1]);
    
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    %text(0.76,0.25,sprintf('r=%.3f',r), 'fontsize',9,'FontWeight','normal');
    %text(0.76,0.12,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
    
    title(sprintf('Across memory trials'), 'FontSize', 9, 'FontWeight', 'normal');
    subtitle(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
    
end


%% (Choice-offload) Trial-level correlation between precision and meta-memory
temptempBoolIndex_A = choiceOffloadBoolIndex';
temptempBoolIndex_B = (~isnan(memoryPrecision_trialLevel)) & (~isnan(meta_trialLevel));

temptempBoolIndex = temptempBoolIndex_A & temptempBoolIndex_B;

[temp_r,temp_p] = corr(memoryPrecision_trialLevel(temptempBoolIndex),meta_trialLevel(temptempBoolIndex));

if true
    fig = figure('Name','Choice-offload trial correlation','NumberTitle','off');
    set(gcf,'Position',[610+200*2 390 196.4*0.905*1.11*0.97 227.5*0.905*1.11*0.97]);
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    x = memoryPrecision_trialLevel(temptempBoolIndex);
    y = meta_trialLevel(temptempBoolIndex);
    
    
    [r,p] = corr(x,y);
    
    mdl = fitglm(x,y);
    
    
    temp_size = 3;%10
    temp_h = scatter(x, y, ...
        temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y);
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    
    %xlim([temp_xmin-(temp_xmax-temp_xmin)*0.08 temp_xmax+(temp_xmax-temp_xmin)*0.08]);
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    
    xticks([0 1]);
    yticks([0 1]);
    
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    %text(0.76,0.25,sprintf('r=%.3f',r), 'fontsize',9,'FontWeight','normal');
    %text(0.76,0.12,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
    
    title(sprintf('Across offload trials'), 'FontSize', 9, 'FontWeight', 'normal');
    subtitle(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
    
end




%% 3 types of trials mixed two-dimension distribution
temptempBoolIndex_A = (~isnan(memoryPrecision_trialLevel)) & (~isnan(meta_trialLevel));
temptempBoolIndex_B = choiceMemoryCorrectBoolIndex';
temptempBoolIndex_C = choiceMemoryErrorBoolIndex';
temptempBoolIndex_D = choiceOffloadBoolIndex';

temptempBoolIndex_CMC = temptempBoolIndex_A & temptempBoolIndex_B;
temptempBoolIndex_CME = temptempBoolIndex_A & temptempBoolIndex_C;
temptempBoolIndex_CF = temptempBoolIndex_A & temptempBoolIndex_D;

temptempBoolIndex_1 = temptempBoolIndex_CMC;
temptempBoolIndex_2 = temptempBoolIndex_CF;
temptempBoolIndex_3 = temptempBoolIndex_CME;

tempColor_1 = color_choiceMemory;
tempColor_2 = color_choiceOffload;
tempColor_3 = color_choiceMemoryError;
% tempColor_12 = [0.3010 0.7450 0.9330]+[0.2 -0.1 -0.2];%[0.3010 0.7450 0.9330]
% tempColor_12 = [68 114 196]./255;
tempColor_12 = [0.3010 0.7450 0.9330]+[0 -0 -0];%[0.3010 0.7450 0.9330]

if true
    fig = figure('Name','A Mixed (3) two-dimension distribution','NumberTitle','off');
    set(gcf,'Position',[610 390 196.4*0.905 227.5*0.905]);
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    x1 = memoryPrecision_trialLevel(temptempBoolIndex_1);
    y1 = meta_trialLevel(temptempBoolIndex_1);
    
    x2 = memoryPrecision_trialLevel(temptempBoolIndex_2);
    y2 = meta_trialLevel(temptempBoolIndex_2);
    
    x3 = memoryPrecision_trialLevel(temptempBoolIndex_3);
    y3 = meta_trialLevel(temptempBoolIndex_3);
    
    [temp_xmin,temp_xmax] = bounds([x1;x2;x3]);
    [temp_ymin,temp_ymax] = bounds([y1;y2;y3]);
    
    temp_flag_mean0_median1 = 1;
    temp_flag_std0_sem1 = 0;
    
    if temp_flag_mean0_median1 == 0
        x1_mu = mean(x1);%
        y1_mu = mean(y1);
        x2_mu = mean(x2);
        y2_mu = mean(y2);
        x3_mu = mean(x3);
        y3_mu = mean(y3);
    elseif temp_flag_mean0_median1 == 1
        x1_mu = median(x1);
        y1_mu = median(y1);
        x2_mu = median(x2);
        y2_mu = median(y2);
        x3_mu = median(x3);
        y3_mu = median(y3);
    end
    
    if temp_flag_std0_sem1 == 0
        x1_std = std(x1);
        y1_std = std(y1);
        
        x2_std = std(x2);
        y2_std = std(y2);
        
        x3_std = std(x3);
        y3_std = std(y3);
    elseif temp_flag_std0_sem1 == 1
        x1_std = std(x1)/sqrt(length(x1));
        y1_std = std(y1)/sqrt(length(y1));
        
        x2_std = std(x2)/sqrt(length(x2));
        y2_std = std(y2)/sqrt(length(y2));
        
        x3_std = std(x3)/sqrt(length(x3));
        y3_std = std(y3)/sqrt(length(y3));
    end
    
    
    
    temp_size = 70;%3,5,7,35,20
    temp_MarkerFaceAlpha = 0.1; %0.5,0.3,0.25,0.125,0.175
    %temp_size = 7;%3,5,7,35,20
    %temp_MarkerFaceAlpha = 0.7; %0.5,0.3,0.25,0.125,0.175
    temp_MarkerEdgeAlpha = 0; %0.7,0
    
    temp_h1 = scatter(x1, y1, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_1, ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    hold on
    
    temp_h2 = scatter(x2, y2, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_2, ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    hold on
    
    temp_h3 = scatter(x3, y3, ...
        temp_size, 'filled', 'MarkerFaceColor', [0 0 0], ... %color_choiceMemoryError
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    hold on
    
    
    
    mdl = fitglm([x1_mu,x2_mu],[y1_mu,y2_mu]);
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    plot(x_fit, y_fit, '-', 'LineWidth', 1, 'Color', [0 0 0]); %[0.35 0.35 0.35 0.7]
    hold on
    
    temp_y3_diff = y3_mu-predict(mdl,x3_mu);
    y_fit = predict(mdl,x_fit')' + temp_y3_diff;
    plot(x_fit, y_fit, '-', 'LineWidth', 1, 'Color', [0 0 0]); %[0.35 0.35 0.35 0.7]
    hold on
    
    temp_slope_A = mdl.Coefficients.Estimate(2);
    
    temp_angle_A = atan(temp_slope_A)*180/pi;
    
    temp_angle_B = temp_angle_A + 90;
    
    temp_slope_B = tan(temp_angle_B*pi/180);
    
    
    x1_proj = (y1-temp_slope_A.*x1)./(temp_slope_B-temp_slope_A);
    y1_proj = temp_slope_B.*x1_proj;
    proj1 = x1_proj./(cos((180-temp_angle_B)*pi/180));
    
    x2_proj = (y2-temp_slope_A.*x2)./(temp_slope_B-temp_slope_A);
    y2_proj = temp_slope_B.*x2_proj;
    proj2 = x2_proj./(cos((180-temp_angle_B)*pi/180));
    
    
    x3_proj = (y3-temp_slope_A.*x3)./(temp_slope_B-temp_slope_A);
    y3_proj = temp_slope_B.*x3_proj;
    proj3 = x3_proj./(cos((180-temp_angle_B)*pi/180));
    
    
    
    temp_LineWidth = 1.5; %1
    
    a = x1_std; % horizontal radius
    b = y1_std; % vertical radius
    temp_x0 = x1_mu; % x0,y0 ellipse centre coordinates
    temp_y0 = y1_mu;
    t = -pi:0.01:pi;
    temp_x = temp_x0 + a*cos(t);
    temp_y = temp_y0 + b*sin(t);
    plot(temp_x,temp_y,'LineWidth',temp_LineWidth,'Color',tempColor_1)
    hold on
    
    a = x2_std; % horizontal radius
    b = y2_std; % vertical radius
    temp_x0 = x2_mu; % x0,y0 ellipse centre coordinates
    temp_y0 = y2_mu;
    t = -pi:0.01:pi;
    temp_x = temp_x0 + a*cos(t);
    temp_y = temp_y0 + b*sin(t);
    plot(temp_x,temp_y,'LineWidth',temp_LineWidth,'Color',tempColor_2)
    hold on
    
    a = x3_std; % horizontal radius
    b = y3_std; % vertical radius
    temp_x0 = x3_mu; % x0,y0 ellipse centre coordinates
    temp_y0 = y3_mu;
    t = -pi:0.01:pi;
    temp_x = temp_x0 + a*cos(t);
    temp_y = temp_y0 + b*sin(t);
    plot(temp_x,temp_y,'LineWidth',temp_LineWidth,'Color',tempColor_3)
    hold on
    
    
    temp_size = 70;%70
    temp_MarkerFaceAlpha = 1;
    temp_MarkerEdgeAlpha = 1;
    
    scatter(x1_mu, y1_mu, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_1,'MarkerEdgeColor',[1 1 1], ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    hold on
    scatter(x2_mu, y2_mu, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_2,'MarkerEdgeColor',[1 1 1], ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    hold on
    scatter(x3_mu, y3_mu, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_3,'MarkerEdgeColor',[1 1 1], ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    hold on
    
    
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    
    xticks([0 1]);
    yticks([0 1]);
    
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
    
    title(sprintf('asd'), 'FontSize', 9, 'FontWeight', 'normal');
    subtitle(sprintf('asd'), 'FontSize', 8, 'FontWeight', 'normal');
    
    
    proj1;
    proj2;
    proj3;
    
    proj12 = [proj1;proj2];
    
    [~,temp_p_proj] = ttest2(proj12,proj3);
end


if true
    fig = figure('Name','B Mixed (3) two-dimension distribution','NumberTitle','off');
    %set(gcf,'Position',[610 390 196.4*0.905 227.5*0.905]);
    set(gcf,'Position',[610 390 196.4*0.905*0.9 227.5*0.905*0.9]);
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    x1 = memoryPrecision_trialLevel(temptempBoolIndex_1);
    y1 = meta_trialLevel(temptempBoolIndex_1);
    
    x2 = memoryPrecision_trialLevel(temptempBoolIndex_2);
    y2 = meta_trialLevel(temptempBoolIndex_2);
    
    x3 = memoryPrecision_trialLevel(temptempBoolIndex_3);
    y3 = meta_trialLevel(temptempBoolIndex_3);
    
    [temp_xmin,temp_xmax] = bounds([x1;x2;x3]);
    [temp_ymin,temp_ymax] = bounds([y1;y2;y3]);
    
    temp_flag_mean0_median1 = 1;
    temp_flag_std0_sem1 = 0;
    
    if temp_flag_mean0_median1 == 0
        x1_mu = mean(x1);%
        y1_mu = mean(y1);
        x2_mu = mean(x2);
        y2_mu = mean(y2);
        x3_mu = mean(x3);
        y3_mu = mean(y3);
    elseif temp_flag_mean0_median1 == 1
        x1_mu = median(x1);
        y1_mu = median(y1);
        x2_mu = median(x2);
        y2_mu = median(y2);
        x3_mu = median(x3);
        y3_mu = median(y3);
    end
    
    if temp_flag_std0_sem1 == 0
        x1_std = std(x1);
        y1_std = std(y1);
        
        x2_std = std(x2);
        y2_std = std(y2);
        
        x3_std = std(x3);
        y3_std = std(y3);
    elseif temp_flag_std0_sem1 == 1
        x1_std = std(x1)/sqrt(length(x1));
        y1_std = std(y1)/sqrt(length(y1));
        
        x2_std = std(x2)/sqrt(length(x2));
        y2_std = std(y2)/sqrt(length(y2));
        
        x3_std = std(x3)/sqrt(length(x3));
        y3_std = std(y3)/sqrt(length(y3));
    end
    
    
    
    %temp_size = 70;%3,5,7,35,20
    %temp_MarkerFaceAlpha = 0.1; %0.5,0.3,0.25,0.125,0.175
    temp_size = 15;%3,5,7,35,20,15
    temp_MarkerFaceAlpha = 0.45; %0.5,0.3,0.25,0.125,0.175,0.625
    temp_MarkerEdgeAlpha = 0; %0.7,0
    
    temp_h1 = scatter(x1, y1, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_1, ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    hold on
    
    temp_h2 = scatter(x2, y2, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_2, ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    hold on
    
    temp_h3 = scatter(x3, y3, ...
        temp_size, 'filled', 'MarkerFaceColor', [0 0 0], ... %color_choiceMemoryError
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    hold on
    
    
    
    mdl = fitglm([x1_mu,x2_mu],[y1_mu,y2_mu]);
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', tempColor_12); %[0.35 0.35 0.35 0.7]
    hold on
    
    temp_y3_diff = y3_mu-predict(mdl,x3_mu);
    y_fit = predict(mdl,x_fit')' + temp_y3_diff;
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0 0 0]); %[0.35 0.35 0.35 0.7]
    hold on
    
    temp_slope_A = mdl.Coefficients.Estimate(2);
    
    temp_angle_A = atan(temp_slope_A)*180/pi;
    
    temp_angle_B = temp_angle_A + 90;
    
    temp_slope_B = tan(temp_angle_B*pi/180);
    
    
    x1_proj = (y1-temp_slope_A.*x1)./(temp_slope_B-temp_slope_A);
    y1_proj = temp_slope_B.*x1_proj;
    proj1 = x1_proj./(cos((180-temp_angle_B)*pi/180));
    
    x2_proj = (y2-temp_slope_A.*x2)./(temp_slope_B-temp_slope_A);
    y2_proj = temp_slope_B.*x2_proj;
    proj2 = x2_proj./(cos((180-temp_angle_B)*pi/180));
    
    
    x3_proj = (y3-temp_slope_A.*x3)./(temp_slope_B-temp_slope_A);
    y3_proj = temp_slope_B.*x3_proj;
    proj3 = x3_proj./(cos((180-temp_angle_B)*pi/180));
    
    
    
    x1_mu_proj = (y1_mu-temp_slope_A.*x1_mu)./(temp_slope_B-temp_slope_A);
    y1_mu_proj = temp_slope_B.*x1_mu_proj;
    proj1_mu = x1_mu_proj./(cos((180-temp_angle_B)*pi/180));
    
    x2_mu_proj = (y2_mu-temp_slope_A.*x2_mu)./(temp_slope_B-temp_slope_A);
    y2_mu_proj = temp_slope_B.*x2_mu_proj;
    proj2_mu = x2_mu_proj./(cos((180-temp_angle_B)*pi/180));
    
    x3_mu_proj = (y3_mu-temp_slope_A.*x3_mu)./(temp_slope_B-temp_slope_A);
    y3_mu_proj = temp_slope_B.*x3_mu_proj;
    proj3_mu = x3_mu_proj./(cos((180-temp_angle_B)*pi/180));
    
    proj12_mu = (proj1_mu+proj2_mu)/2;
    
    %     temp_LineWidth = 1.5; %1
    %
    %     a = x1_std; % horizontal radius
    %     b = y1_std; % vertical radius
    %     temp_x0 = x1_mu; % x0,y0 ellipse centre coordinates
    %     temp_y0 = y1_mu;
    %     t = -pi:0.01:pi;
    %     temp_x = temp_x0 + a*cos(t);
    %     temp_y = temp_y0 + b*sin(t);
    %     plot(temp_x,temp_y,'LineWidth',temp_LineWidth,'Color',tempColor_1)
    %     hold on
    %
    %     a = x2_std; % horizontal radius
    %     b = y2_std; % vertical radius
    %     temp_x0 = x2_mu; % x0,y0 ellipse centre coordinates
    %     temp_y0 = y2_mu;
    %     t = -pi:0.01:pi;
    %     temp_x = temp_x0 + a*cos(t);
    %     temp_y = temp_y0 + b*sin(t);
    %     plot(temp_x,temp_y,'LineWidth',temp_LineWidth,'Color',tempColor_2)
    %     hold on
    %
    %     a = x3_std; % horizontal radius
    %     b = y3_std; % vertical radius
    %     temp_x0 = x3_mu; % x0,y0 ellipse centre coordinates
    %     temp_y0 = y3_mu;
    %     t = -pi:0.01:pi;
    %     temp_x = temp_x0 + a*cos(t);
    %     temp_y = temp_y0 + b*sin(t);
    %     plot(temp_x,temp_y,'LineWidth',temp_LineWidth,'Color',tempColor_3)
    %     hold on
    
    
    temp_size = 70;%70
    temp_MarkerFaceAlpha = 1;
    temp_MarkerEdgeAlpha = 1;
    
    scatter(x1_mu, y1_mu, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_1,'MarkerEdgeColor',[1 1 1], ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha,'LineWidth',1);
    hold on
    scatter(x2_mu, y2_mu, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_2,'MarkerEdgeColor',[1 1 1], ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha,'LineWidth',1);
    hold on
    scatter(x3_mu, y3_mu, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_3,'MarkerEdgeColor',[1 1 1], ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha,'LineWidth',1);
    hold on
    
    
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    
    xticks([0 1]);
    yticks([0 1]);
    
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
    
    title(sprintf('asd'), 'FontSize', 9, 'FontWeight', 'normal');
    subtitle(sprintf('asd'), 'FontSize', 8, 'FontWeight', 'normal');
    
    
    proj1;
    proj2;
    proj3;
    
    proj12 = [proj1;proj2];
    
    [~,temp_p_proj] = ttest2(proj12,proj3);
    
    proj123 = [proj1;proj2;proj3];
    [proj123_min,proj123_max] = bounds(proj123);
    
    
    
    if true
        fig = figure('Name','Proj of Mixed (3) two-dimension','NumberTitle','off');
        %set(gcf,'Position',[210 260 100 240*0.9*0.73*0.95*0.95*0.89*1.42*0.8*0.94*1.01*0.88]);
        %set(gcf,'Position',[210 260 100*1.38*1.07*0.84*1.03 120.2*0.46*1.32]);
        set(gcf,'Position',[210 260 100*1.38*1.07*0.84*1.03*0.89 120.2*0.46*1.32]);
        
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        
        x1 = proj1;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf1,xmesh1,bandwidth1] = ksdensity(x1','NumPoints',n,'Function','pdf','Support',[proj123_min-0.01,proj123_max+0.01],'BoundaryCorrection', 'Reflection');
        temp_ratio1 = 1;
        y1 = pdf1*temp_ratio1;
        %h1 = plot(xmesh1,y1,'LineWidth',1.5,'color',tempColor_1);
        %hold on
        
        
        x2 = proj2;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[proj123_min-0.01,proj123_max+0.01],'BoundaryCorrection', 'Reflection');
        temp_ratio2 = 1;
        y2 = pdf2*temp_ratio2;
        %h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',[0 0.4470 0.7410]);
        %hold on
        
        x12 = proj12;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf12,xmesh12,bandwidth12] = ksdensity(x12','NumPoints',n,'Function','pdf','Support',[proj123_min-0.01,proj123_max+0.01],'BoundaryCorrection', 'Reflection');
        temp_ratio12 = 1;
        y12 = pdf12*temp_ratio12;
        h12 = plot(xmesh12,y12,'LineWidth',1.5,'color',tempColor_12);%[0 0.4470 0.7410],[0.3010 0.7450 0.9330]
        hold on
        
        
        x3 = proj3;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[proj123_min-0.01,proj123_max+0.01],'BoundaryCorrection', 'Reflection');
        temp_ratio3 = 1;
        y3 = pdf3*temp_ratio3;
        h3 = plot(xmesh3,y3,'LineWidth',1.5,'color',tempColor_3); %#ok<*NASGU>
        hold on
        
        
        %[x1_min,x1_max] = bounds(xmesh1);
        %[x2_min,x2_max] = bounds(xmesh2);
        [x12_min,x12_max] = bounds(xmesh12);
        [x3_min,x3_max] = bounds(xmesh3);
        %x_min = min([x1_min,x2_min,x3_min]);
        %x_max = max([x1_max,x2_max,x3_max]);
        x_min = min([x12_min,x3_min]);
        x_max = max([x12_max,x3_max]);
        
        %[y1_min,y1_max] = bounds(y1);
        %[y2_min,y2_max] = bounds(y2);
        [y12_min,y12_max] = bounds(y12);
        [y3_min,y3_max] = bounds(y3);
        %y_min = min([y1_min,y2_min,y3_min]);
        %y_max = max([y1_max,y2_max,y3_max]);
        y_min = min([y12_min,y3_min]);
        y_max = max([y12_max,y3_max]);
        
        
        
        plot([proj12_mu proj12_mu],[y_min y_max+(y_max-y_min)*0.1],'LineWidth',1.5,'color',tempColor_12);
        hold on
        plot([proj3_mu proj3_mu],[y_min y_max+(y_max-y_min)*0.1],'LineWidth',1.5,'color',[0 0 0]);
        hold on
        
        
        yticks([0 1 2]);
        
        axis off
        
        set(gca,'linewidth',1.5)
        xlim([x_min+(x_max-x_min)*0.01 x_max-(x_max-x_min)*0.01]);
        ylim([y_min y_max+(y_max-y_min)*0.4]);
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        %xlabel('Proj', 'FontSize', 9);
        ylabel('Pdf', 'FontSize', 9);
        
        
        
    end
    
    if true
        fig = figure('Name','Legend','NumberTitle','off');
        set(gcf,'Position',[210 260 400 200]);
        
        t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        h_handle1 = [];
        
        h_handle1 = [h_handle1 plot([0 1], [0 1], 'LineWidth', 3, 'Color', tempColor_1)];
        hold on
        h_handle1 = [h_handle1 plot([0 1], [0 1], 'LineWidth', 3, 'Color', tempColor_2)];
        hold on
        
        le = legend(h_handle1,'Memory-correct','Offload','Location','northwest','fontsize',6.5,'NumColumns',1);
        le.ItemTokenSize = ones(1,3)*10;
        legend('boxoff');
        
        
        nexttile
        
        h_handle1 = [];
        
        h_handle1 = [h_handle1 plot([0 1], [0 1], 'LineWidth', 3, 'Color', tempColor_12)];
        hold on
        h_handle1 = [h_handle1 plot([0 1], [0 1], 'LineWidth', 3, 'Color', tempColor_3)];
        hold on
        
        le = legend(h_handle1,'Memory-correct + Offload','Memory-error','Location','northwest','fontsize',6.5,'NumColumns',1);
        %le = legend(h_handle1,'Match','Mismatch','Location','northwest','fontsize',6.5,'NumColumns',1);
        le.ItemTokenSize = ones(1,3)*10;
        legend('boxoff');
        
        
    end
    
end



if true
    fig = figure('Name','C Mixed (3) two-dimension distribution','NumberTitle','off');
    %set(gcf,'Position',[610 390 196.4*0.905*0.9 227.5*0.905*0.9]);
    set(gcf,'Position',[610 390 196.4*0.905*0.9*0.89 227.5*0.905*0.9]);
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    x1 = memoryPrecision_trialLevel(temptempBoolIndex_1);
    y1 = meta_trialLevel(temptempBoolIndex_1);
    
    x2 = memoryPrecision_trialLevel(temptempBoolIndex_2);
    y2 = meta_trialLevel(temptempBoolIndex_2);
    
    x3 = memoryPrecision_trialLevel(temptempBoolIndex_3);
    y3 = meta_trialLevel(temptempBoolIndex_3);
    
    [temp_xmin,temp_xmax] = bounds([x1;x2;x3]);
    [temp_ymin,temp_ymax] = bounds([y1;y2;y3]);
    
    temp_flag_mean0_median1 = 1;
    temp_flag_std0_sem1 = 0;
    
    if temp_flag_mean0_median1 == 0
        x1_mu = mean(x1);%
        y1_mu = mean(y1);
        x2_mu = mean(x2);
        y2_mu = mean(y2);
        x3_mu = mean(x3);
        y3_mu = mean(y3);
    elseif temp_flag_mean0_median1 == 1
        x1_mu = median(x1);
        y1_mu = median(y1);
        x2_mu = median(x2);
        y2_mu = median(y2);
        x3_mu = median(x3);
        y3_mu = median(y3);
    end
    
    if temp_flag_std0_sem1 == 0
        x1_std = std(x1);
        y1_std = std(y1);
        
        x2_std = std(x2);
        y2_std = std(y2);
        
        x3_std = std(x3);
        y3_std = std(y3);
    elseif temp_flag_std0_sem1 == 1
        x1_std = std(x1)/sqrt(length(x1));
        y1_std = std(y1)/sqrt(length(y1));
        
        x2_std = std(x2)/sqrt(length(x2));
        y2_std = std(y2)/sqrt(length(y2));
        
        x3_std = std(x3)/sqrt(length(x3));
        y3_std = std(y3)/sqrt(length(y3));
    end
    
    
    
    %temp_size = 70;%3,5,7,35,20
    %temp_MarkerFaceAlpha = 0.1; %0.5,0.3,0.25,0.125,0.175
    temp_size = 15;%3,5,7,35,20,15
    temp_MarkerFaceAlpha = 0.45; %0.5,0.3,0.25,0.125,0.175,0.625
    temp_MarkerEdgeAlpha = 0; %0.7,0
    
    temp_h1 = scatter(x1, y1, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_1, ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    hold on
    
    temp_h2 = scatter(x2, y2, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_2, ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    hold on
    
    %     temp_h3 = scatter(x3, y3, ...
    %         temp_size, 'filled', 'MarkerFaceColor', [0 0 0], ... %color_choiceMemoryError
    %         'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    %     hold on
    
    
    
    mdl = fitglm([x1_mu,x2_mu],[y1_mu,y2_mu]);
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', tempColor_12); %[0.35 0.35 0.35 0.7]
    hold on
    
    %     temp_y3_diff = y3_mu-predict(mdl,x3_mu);
    %     y_fit = predict(mdl,x_fit')' + temp_y3_diff;
    %     plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0 0 0]); %[0.35 0.35 0.35 0.7]
    %     hold on
    
    temp_slope_A = mdl.Coefficients.Estimate(2);
    
    temp_angle_A = atan(temp_slope_A)*180/pi;
    
    temp_angle_B = temp_angle_A + 90;
    
    temp_slope_B = tan(temp_angle_B*pi/180);
    
    
    
    temp_size = 70;%70
    temp_MarkerFaceAlpha = 1;
    temp_MarkerEdgeAlpha = 1;
    
    scatter(x1_mu, y1_mu, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_1,'MarkerEdgeColor',[1 1 1], ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha,'LineWidth',1);
    hold on
    scatter(x2_mu, y2_mu, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_2,'MarkerEdgeColor',[1 1 1], ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha,'LineWidth',1);
    hold on
    %     scatter(x3_mu, y3_mu, ...
    %         temp_size, 'filled', 'MarkerFaceColor', tempColor_3,'MarkerEdgeColor',[1 1 1], ...
    %         'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha,'LineWidth',1);
    %     hold on
    
    
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    
    xticks([0 1]);
    yticks([0 1]);
    
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
    
    title(sprintf('asd'), 'FontSize', 9, 'FontWeight', 'normal');
    subtitle(sprintf('asd'), 'FontSize', 8, 'FontWeight', 'normal');
    
    
end



if true
    fig = figure('Name','D Mixed (3) two-dimension distribution','NumberTitle','off');
    %set(gcf,'Position',[610 390 196.4*0.905*0.9 227.5*0.905*0.9]);
    set(gcf,'Position',[610 390 196.4*0.905*0.9*0.89 227.5*0.905*0.9]);
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    x1 = memoryPrecision_trialLevel(temptempBoolIndex_1);
    y1 = meta_trialLevel(temptempBoolIndex_1);
    
    x2 = memoryPrecision_trialLevel(temptempBoolIndex_2);
    y2 = meta_trialLevel(temptempBoolIndex_2);
    
    x3 = memoryPrecision_trialLevel(temptempBoolIndex_3);
    y3 = meta_trialLevel(temptempBoolIndex_3);
    
    [temp_xmin,temp_xmax] = bounds([x1;x2;x3]);
    [temp_ymin,temp_ymax] = bounds([y1;y2;y3]);
    
    temp_flag_mean0_median1 = 1;
    temp_flag_std0_sem1 = 0;
    
    if temp_flag_mean0_median1 == 0
        x1_mu = mean(x1);%
        y1_mu = mean(y1);
        x2_mu = mean(x2);
        y2_mu = mean(y2);
        x3_mu = mean(x3);
        y3_mu = mean(y3);
    elseif temp_flag_mean0_median1 == 1
        x1_mu = median(x1);
        y1_mu = median(y1);
        x2_mu = median(x2);
        y2_mu = median(y2);
        x3_mu = median(x3);
        y3_mu = median(y3);
    end
    
    if temp_flag_std0_sem1 == 0
        x1_std = std(x1);
        y1_std = std(y1);
        
        x2_std = std(x2);
        y2_std = std(y2);
        
        x3_std = std(x3);
        y3_std = std(y3);
    elseif temp_flag_std0_sem1 == 1
        x1_std = std(x1)/sqrt(length(x1));
        y1_std = std(y1)/sqrt(length(y1));
        
        x2_std = std(x2)/sqrt(length(x2));
        y2_std = std(y2)/sqrt(length(y2));
        
        x3_std = std(x3)/sqrt(length(x3));
        y3_std = std(y3)/sqrt(length(y3));
    end
    
    
    
    %temp_size = 70;%3,5,7,35,20
    %temp_MarkerFaceAlpha = 0.1; %0.5,0.3,0.25,0.125,0.175
    temp_size = 15;%3,5,7,35,20,15
    temp_MarkerFaceAlpha = 0.45; %0.5,0.3,0.25,0.125,0.175,0.625
    temp_MarkerEdgeAlpha = 0; %0.7,0
    
    %     temp_h1 = scatter(x1, y1, ...
    %         temp_size, 'filled', 'MarkerFaceColor', tempColor_1, ...
    %         'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    %     hold on
    %
    %     temp_h2 = scatter(x2, y2, ...
    %         temp_size, 'filled', 'MarkerFaceColor', tempColor_2, ...
    %         'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    %     hold on
    
    temp_h12 = scatter([x1;x2], [y1;y2], ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_12-[0 0 0], ...
        'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    hold on
    
    temp_h3 = scatter(x3, y3, ...
        temp_size, 'filled', 'MarkerFaceColor', [0 0 0], ... %color_choiceMemoryError
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
    hold on
    
    
    
    mdl = fitglm([x1_mu,x2_mu],[y1_mu,y2_mu]);
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', tempColor_12); %[0.35 0.35 0.35 0.7]
    hold on
    
    temp_y3_diff = y3_mu-predict(mdl,x3_mu);
    y_fit = predict(mdl,x_fit')' + temp_y3_diff;
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0 0 0]); %[0.35 0.35 0.35 0.7]
    hold on
    
    temp_slope_A = mdl.Coefficients.Estimate(2);
    
    temp_angle_A = atan(temp_slope_A)*180/pi;
    
    temp_angle_B = temp_angle_A + 90;
    
    temp_slope_B = tan(temp_angle_B*pi/180);
    
    
    
    temp_size = 70;%70
    temp_MarkerFaceAlpha = 1;
    temp_MarkerEdgeAlpha = 1;
    
    scatter(x1_mu, y1_mu, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_1,'MarkerEdgeColor',[1 1 1], ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha,'LineWidth',1);
    hold on
    scatter(x2_mu, y2_mu, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_2,'MarkerEdgeColor',[1 1 1], ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha,'LineWidth',1);
    hold on
    %     scatter((x1_mu+x2_mu)/2, (y1_mu+y2_mu)/2, ...
    %         temp_size, 'filled', 'MarkerFaceColor', tempColor_12,'MarkerEdgeColor',[1 1 1], ...
    %         'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha,'LineWidth',1);
    %     hold on
    scatter(x3_mu, y3_mu, ...
        temp_size, 'filled', 'MarkerFaceColor', tempColor_3,'MarkerEdgeColor',[1 1 1], ...
        'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha,'LineWidth',1);
    hold on
    
    
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    
    xticks([0 1]);
    yticks([0 1]);
    
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
    
    title(sprintf('asd'), 'FontSize', 9, 'FontWeight', 'normal');
    subtitle(sprintf('asd'), 'FontSize', 8, 'FontWeight', 'normal');
    
    
end




%% (Match trials) Trial-level correlation between precision and meta-memory
temptempBoolIndex_A = choiceBoolIndex';
temptempBoolIndex_B = (~isnan(memoryPrecision_trialLevel)) & (~isnan(meta_trialLevel));

temptempBoolIndex_AB = temptempBoolIndex_A & temptempBoolIndex_B;

temptempBoolIndex_C = (memoryPrecision_trialLevel<=lowThreshold_memoryPrecision) & (meta_trialLevel<=lowThreshold_meta);
temptempBoolIndex_D = (memoryPrecision_trialLevel>lowThreshold_memoryPrecision) & (meta_trialLevel>lowThreshold_meta);

temptempBoolIndex_CD = temptempBoolIndex_C | temptempBoolIndex_D;

temptempBoolIndex = temptempBoolIndex_AB & temptempBoolIndex_CD;

if true
    fig = figure('Name','Match trial correlation','NumberTitle','off');
    set(gcf,'Position',[610 390 196.4*0.905*0.915 227.5*0.905*0.915]);
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    x = memoryPrecision_trialLevel(temptempBoolIndex);
    y = meta_trialLevel(temptempBoolIndex);
    
    
    [r,p] = corr(x,y);
    
    mdl = fitglm(x,y);
    
    
    temp_size = 3;%10
    temp_h = scatter(x, y, ...
        temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y); %#ok<*ASGLU>
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    
    xticks([0 1]);
    yticks([0 1]);
    
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
    
    title(sprintf('Match trials'), 'FontSize', 9, 'FontWeight', 'normal');
    subtitle(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
    
end


if true
    fig = figure('Name','Match trial correlation (baseline)','NumberTitle','off');
    %set(gcf,'Position',[610 390 196.4*0.905 227.5*0.905]);
    set(gcf,'Position',[610 390 196.4*0.905*0.915 227.5*0.905*0.915]);
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    x = meta_trialLevel_baseline(temptempBoolIndex);
    y = meta_trialLevel(temptempBoolIndex);
    
    
    [r,p] = corr(x,y);
    
    mdl = fitglm(x,y);
    
    
    temp_size = 3;%10
    temp_h = scatter(x, y, ...
        temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y); %#ok<*ASGLU>
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    
    xticks([0 1]);
    yticks([0 1]);
    
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('Baseline meta'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
    
    title(sprintf('Match trials'), 'FontSize', 9, 'FontWeight', 'normal');
    subtitle(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
    
end



%% (Mismatch trials) Trial-level correlation between precision and meta-memory
temptempBoolIndex_A = choiceBoolIndex';
temptempBoolIndex_B = (~isnan(memoryPrecision_trialLevel)) & (~isnan(meta_trialLevel));

temptempBoolIndex_AB = temptempBoolIndex_A & temptempBoolIndex_B;

temptempBoolIndex_C = (memoryPrecision_trialLevel<=lowThreshold_memoryPrecision) & (meta_trialLevel<=lowThreshold_meta);
temptempBoolIndex_D = (memoryPrecision_trialLevel>lowThreshold_memoryPrecision) & (meta_trialLevel>lowThreshold_meta);

temptempBoolIndex_CD = temptempBoolIndex_C | temptempBoolIndex_D;

temptempBoolIndex = temptempBoolIndex_AB & (~temptempBoolIndex_CD);

if true
    fig = figure('Name','Mismatch trial correlation','NumberTitle','off');
    %set(gcf,'Position',[610 390 196.4*0.905 227.5*0.905]);
    set(gcf,'Position',[610 390 196.4*0.905*0.915 227.5*0.905*0.915]);
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    x = memoryPrecision_trialLevel(temptempBoolIndex);
    y = meta_trialLevel(temptempBoolIndex);
    
    
    [r,p] = corr(x,y);
    
    mdl = fitglm(x,y);
    
    
    temp_size = 3;%10
    temp_h = scatter(x, y, ...
        temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y); %#ok<*ASGLU>
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    
    xticks([0 1]);
    yticks([0 1]);
    
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
    
    title(sprintf('Mismatch trials'), 'FontSize', 9, 'FontWeight', 'normal');
    subtitle(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
    
end


if true
    fig = figure('Name','Mismatch trial correlation (baseline)','NumberTitle','off');
    set(gcf,'Position',[610 390 196.4*0.905*0.915 227.5*0.905*0.915]);
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    x = meta_trialLevel_baseline(temptempBoolIndex);
    y = meta_trialLevel(temptempBoolIndex);
    
    
    [r,p] = corr(x,y);
    
    mdl = fitglm(x,y);
    
    
    temp_size = 3;%10
    temp_h = scatter(x, y, ...
        temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y); %#ok<*ASGLU>
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    
    xticks([0 1]);
    yticks([0 1]);
    
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('Baseline meta'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
    
    title(sprintf('Mismatch trials'), 'FontSize', 9, 'FontWeight', 'normal');
    subtitle(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
    
end




if true
    
    %% 2d distrition
    fig = figure('Name','example seq distribution','NumberTitle','off');
    set(gcf,'Position',[680 450 201.4*1.08*0.95*0.98*0.97 179.5*1.08*0.97]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,5,'TileSpacing','compact','Padding','loose');
    
    temp_histogram_numBins = histogram_numBins;
    
    nexttile([1 4])
    
    x_raw = memoryPrecision_trialLevel(choiceBoolIndex);
    y_raw = meta_trialLevel(choiceBoolIndex);
    
    
    x = x_raw(~isnan(x_raw));
    y = y_raw(~isnan(x_raw));
    
    [x_min,x_max] = bounds(x);
    x_max = x_max + 0.01;
    
    XYBinLimits = [0 1];
    
    h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','count');
    h.XBinLimits = XYBinLimits;
    h.YBinLimits = XYBinLimits;
    h.NumBins = [1 1].*temp_histogram_numBins;
    hold on
    
    xticks([0 1]);
    yticks([0 1]);
    
    set(gca,'Ydir','normal') %reverse
    xlim(XYBinLimits);
    ylim(XYBinLimits);
    xtickangle(0);
    set(gca, 'FontSize', 8) %12
    
    c = colorbar('eastoutside','FontSize',8);
    c.Position = c.Position+[0.03 0.10 -0.02 -0.13];
    
    if if_colormap_loadEnhanced == 1
        load('parula_enhanced');
        colormap(parula_enhanced);
    elseif if_colormap_loadEnhanced == 0
        colormap parula
    end
    
    xlabel('Memory strength', 'FontSize', 9, 'FontWeight', 'normal');
    ylabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
    
end

%% (Free-choice) Trial-level correlation between location correlation and meta-memory
if true
    temptempBoolIndex_A = choiceBoolIndex';
    temptempBoolIndex_B = (~isnan(memoryPrecision_trialLevel)) & (~isnan(meta_trialLevel));
    
    temptempBoolIndex = temptempBoolIndex_A & temptempBoolIndex_B;
    
    fig = figure('Name','meta_seqLevel','NumberTitle','off');
    %set(gcf,'Position',[610 390 337*0.8*0.99*1.15 176*1.065*1.2*0.85]);
    %set(gcf,'Position',[610 390 337*0.8*0.99*1.15*0.64 176*1.065*1.2*0.85*0.85*1.4]);
    set(gcf,'Position',[610 390 196.4*0.905 227.5*0.905]);
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    x = locCorr_trialLevel(temptempBoolIndex);
    y = meta_trialLevel(temptempBoolIndex);
    
    
    [r,p] = corr(x,y);
    
    mdl = fitglm(x,y);
    
    
    temp_size = 3;%10
    temp_h = scatter(x, y, ...
        temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y);
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    
    %xlim([temp_xmin-(temp_xmax-temp_xmin)*0.08 temp_xmax+(temp_xmax-temp_xmin)*0.08]);
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    
    xticks([0 1]);
    yticks([0 1]);
    
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    %text(0.76,0.25,sprintf('r=%.3f',r), 'fontsize',9,'FontWeight','normal');
    %text(0.76,0.12,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('Location Correlation'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-memory'), 'FontSize', 9, 'FontWeight', 'normal');
    
    title(sprintf('Free-choice'), 'FontSize', 9, 'FontWeight', 'normal');
    subtitle(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
    
end

%%
memoryPrecision_trialLevel;
meta_trialLevel;

temptempBoolIndex_A = choiceMemoryCorrectBoolIndex';
temptempBoolIndex_B = (~isnan(memoryPrecision_trialLevel)) & (~isnan(meta_trialLevel));

temptempBoolIndex = temptempBoolIndex_A & temptempBoolIndex_B;

[temp_r,temp_p] = corr(memoryPrecision_trialLevel(temptempBoolIndex),meta_trialLevel(temptempBoolIndex));

%% (Choice-memory correct) Trial-level correlation between precision and meta-memory
if true
    fig = figure('Name','meta_seqLevel','NumberTitle','off');
    %set(gcf,'Position',[610 390 337*0.8*0.99*1.15 176*1.065*1.2*0.85]);
    %set(gcf,'Position',[610 390 337*0.8*0.99*1.15*0.64 176*1.065*1.2*0.85*0.85*1.4]);
    set(gcf,'Position',[610 390 196.4*0.905 227.5*0.905]);
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    x = memoryPrecision_trialLevel(temptempBoolIndex);
    y = meta_trialLevel(temptempBoolIndex);
    
    
    [r,p] = corr(x,y);
    
    mdl = fitglm(x,y);
    
    
    temp_size = 3;%10
    temp_h = scatter(x, y, ...
        temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    [temp_xmin,temp_xmax] = bounds(x);
    [temp_ymin,temp_ymax] = bounds(y);
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    
    %xlim([temp_xmin-(temp_xmax-temp_xmin)*0.08 temp_xmax+(temp_xmax-temp_xmin)*0.08]);
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    
    xticks([0 1]);
    yticks([0 1]);
    
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    %text(0.76,0.25,sprintf('r=%.3f',r), 'fontsize',9,'FontWeight','normal');
    %text(0.76,0.12,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('Memory strength'), 'FontSize', 9, 'FontWeight', 'normal');
    ylabel(sprintf('Meta-memory'), 'FontSize', 9, 'FontWeight', 'normal');
    
    title(sprintf('Choice-memory correct'), 'FontSize', 9, 'FontWeight', 'normal');
    subtitle(sprintf('r = %.3f, %s',r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
    
end


%%
memoryPrecision_trialLevel;
meta_trialLevel;
if_memeoryPrecision_stimuli0_response1;

choiceMemoryErrorBoolIndex;
choiceMemoryCorrectBoolIndex;

temptempBoolIndex_B = (~isnan(memoryPrecision_trialLevel)) & (~isnan(meta_trialLevel));

temptempBoolIndex_A = choiceMemoryErrorBoolIndex';
temptempBoolIndex_error = temptempBoolIndex_A & temptempBoolIndex_B;

temptempBoolIndex_A = choiceMemoryCorrectBoolIndex';
% temptempBoolIndex_A = allMemoryCorrectBoolIndex';
% temptempBoolIndex_A = choiceMemoryBoolIndex';
temptempBoolIndex_correct = temptempBoolIndex_A & temptempBoolIndex_B;

temptempBoolIndex_A = choiceOffloadBoolIndex';
temptempBoolIndex_offload = temptempBoolIndex_A & temptempBoolIndex_B;

temp1 = memoryPrecision_trialLevel(temptempBoolIndex_offload);
temp2 = memoryPrecision_trialLevel(temptempBoolIndex_error);
temp3 = memoryPrecision_trialLevel(temptempBoolIndex_correct);

temp1_mean = mean(temp1);
temp2_mean = mean(temp2);
temp3_mean = mean(temp3);

[~,temp_p_precisionOffloadError] = ttest2(temp1,temp2);
[~,temp_p_precisionErrorCorrect] = ttest2(temp2,temp3);
[~,temp_p_precisionOffloadCorrect] = ttest2(temp1,temp3);

temp4 = meta_trialLevel(temptempBoolIndex_offload);
temp5 = meta_trialLevel(temptempBoolIndex_error);
temp6 = meta_trialLevel(temptempBoolIndex_correct);

temp4_mean = mean(temp4);
temp5_mean = mean(temp5);
temp6_mean = mean(temp6);

[~,temp_p_metaOffloadError] = ttest2(temp4,temp5);
[~,temp_p_metaErrorCorrect] = ttest2(temp4,temp6);
[~,temp_p_metaOffloadCorrect] = ttest2(temp5,temp6);


%% Plot: Compare memory precision of offload, error, correct trials
if false
    
    fig = figure('Name','example seq distribution','NumberTitle','off');
    set(gcf,'Position',[680 450 201.4*1.08*0.95*0.98*0.97 179.5*1.08*0.97*1.5*1.3*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    
    
    nexttile
    
    temp_p12 = temp_p_precisionOffloadError;
    temp_p23 = temp_p_precisionErrorCorrect;
    temp_p13 = temp_p_precisionOffloadCorrect;
    
    temp_1 = temp1;
    temp_2 = temp2;
    temp_3 = temp3;
    
    
    temp_y_min = min([temp_1;temp_2;temp_3]);
    temp_y_max = max([temp_1;temp_2;temp_3]);
    
    
    temp_data = [temp_1;temp_2;temp_3];
    
    g1 = repmat({'A'},length(temp_1),1);
    g2 = repmat({'B'},length(temp_2),1);
    g3 = repmat({'C'},length(temp_3),1);
    
    temp_label = [g1;g2;g3];
    
    temptemp_color1 = [1 1 1]*0.5;
    temptemp_color2 = repmat(temptemp_color1, 3, 1);
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',[{'A'};{'B'};{'C'}]);
    h(1).ViolinPlot.FaceAlpha = 0.1;
    h(2).ViolinPlot.FaceAlpha = 0.1;
    h(3).ViolinPlot.FaceAlpha = 0.1;
    
    
    tempTxt = sprintf('');
    if temp_p12 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p12 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p12 < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.12,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([1.1 1.9],temp_y_max+(temp_y_max-temp_y_min)*0.10*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    tempTxt = sprintf('');
    if temp_p13 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p13 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p13 < 0.05
        tempTxt = sprintf('*');
    end
    text(2,temp_y_max+(temp_y_max-temp_y_min)*0.22,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([1.1 2.9],temp_y_max+(temp_y_max-temp_y_min)*0.20*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    tempTxt = sprintf('');
    if temp_p23 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p23 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p23 < 0.05
        tempTxt = sprintf('*');
    end
    text(2.5,temp_y_max+(temp_y_max-temp_y_min)*0.12,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([2.1 2.9],temp_y_max+(temp_y_max-temp_y_min)*0.10*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    set(gca,'linewidth',1.5)
    xlim([0.15 3.85])
    %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.30]);
    ylim([-0.125 1.30])
    set(gca, 'FontSize', 8);%10
    
    xticks([1 2 3]);
    
    xtl = ["Offload";"Error";"Correct"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.22;%0.56,0.4
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9
    set(gca,'xticklabel','');
    
    
    yticks([0 1]);
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Memory strength', 'FontSize', 9, 'FontWeight', 'normal');%10
    if if_memeoryPrecision_stimuli0_response1 == 0
        title(sprintf('Stimuli-labeled'),'fontsize',9);
    elseif if_memeoryPrecision_stimuli0_response1 == 1
        title(sprintf('Response-labeled'),'fontsize',9);
    end
    
end

%% Plot: Compare memory precision of error, correct trials
if true
    
    fig = figure('Name','example seq distribution','NumberTitle','off');
    %set(gcf,'Position',[480 450 201.4*1.08*0.95*0.98*0.97 179.5*1.08*0.97*1.5*1.3*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[480 450 201.4*1.08*0.95*0.98*0.97*0.65*0.8*1.5 179.5*1.08*0.97*1.5*1.3*1.05*0.43]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[480 450 201.4*1.08*0.95*0.98*0.97*0.65*0.8*1.5 179.5*1.08*0.97*1.5*1.3*1.05*0.43*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[480 450 201.4*1.08*0.95*0.98*0.97*0.65*0.8*1.5*0.84 179.5*1.08*0.97*1.5*1.3*1.05*0.43*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    
    
    nexttile
    
    temp_p12 = temp_p_precisionErrorCorrect;
    
    temp_1 = temp2;
    temp_2 = temp3;
    
    
    temp_y_min = min([temp_1;temp_2]);
    temp_y_max = max([temp_1;temp_2]);
    
    
    temp_data = [temp_1;temp_2];
    
    g1 = repmat({'A'},length(temp_1),1);
    g2 = repmat({'B'},length(temp_2),1);
    
    temp_label = [g1;g2];
    
    temptemp_color1 = [1 1 1]*0.5;
    temptemp_color2 = repmat(temptemp_color1, 2, 1);
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',[{'A'};{'B'}]);
    h(1).ViolinPlot.FaceAlpha = 0.1;
    h(2).ViolinPlot.FaceAlpha = 0.1;
    
    
    tempTxt = sprintf('');
    if temp_p12 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p12 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p12 < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.20,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([1.1 1.9],temp_y_max+(temp_y_max-temp_y_min)*0.16*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    set(gca,'linewidth',1.5)
    xlim([0.45 2.55])
    %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.30]);
    ylim([-0.125 1.30])
    set(gca, 'FontSize', 8);%10
    
    xticks([1 2]);
    
    xtl = ["Error";"Correct"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.22;%0.56,0.4
    %xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.37;
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9
    
    xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.52;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',8);%9
    set(gca,'xticklabel','');
    
    
    yticks([0 1]);
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('WM strength', 'FontSize', 8, 'FontWeight', 'normal');%10
    %     if if_memeoryPrecision_stimuli0_response1 == 0
    %         title(sprintf('Stimuli-labeled'),'fontsize',9);
    %     elseif if_memeoryPrecision_stimuli0_response1 == 1
    %         title(sprintf('Response-labeled'),'fontsize',9);
    %     end
    
end



%% Plot: Compare meta-memory  of offload, error, correct trials
if false
    
    fig = figure('Name','example seq distribution','NumberTitle','off');
    %set(gcf,'Position',[680 450 201.4*1.08*0.95*0.98*0.97 179.5*1.08*0.97*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[680 450 201.4*1.08*0.95*0.98*0.97 179.5*1.08*0.97*1.5*1.3*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp_p12 = temp_p_metaOffloadError;
    temp_p23 = temp_p_metaErrorCorrect;
    temp_p13 = temp_p_metaOffloadCorrect;
    
    temp_1 = temp4;
    temp_2 = temp5;
    temp_3 = temp6;
    
    
    temp_y_min = min([temp_1;temp_2;temp_3]);
    temp_y_max = max([temp_1;temp_2;temp_3]);
    
    
    temp_data = [temp_1;temp_2;temp_3];
    
    g1 = repmat({'A'},length(temp_1),1);
    g2 = repmat({'B'},length(temp_2),1);
    g3 = repmat({'C'},length(temp_3),1);
    
    temp_label = [g1;g2;g3];
    
    temptemp_color1 = [1 1 1]*0.5;
    temptemp_color2 = repmat(temptemp_color1, 3, 1);
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',[{'A'};{'B'};{'C'}]);
    h(1).ViolinPlot.FaceAlpha = 0.1;
    h(2).ViolinPlot.FaceAlpha = 0.1;
    h(3).ViolinPlot.FaceAlpha = 0.1;
    
    
    tempTxt = sprintf('');
    if temp_p12 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p12 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p12 < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.12,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([1.1 1.9],temp_y_max+(temp_y_max-temp_y_min)*0.10*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    tempTxt = sprintf('');
    if temp_p13 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p13 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p13 < 0.05
        tempTxt = sprintf('*');
    end
    text(2,temp_y_max+(temp_y_max-temp_y_min)*0.22,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([1.1 2.9],temp_y_max+(temp_y_max-temp_y_min)*0.20*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    tempTxt = sprintf('');
    if temp_p23 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p23 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p23 < 0.05
        tempTxt = sprintf('*');
    end
    text(2.5,temp_y_max+(temp_y_max-temp_y_min)*0.12,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([2.1 2.9],temp_y_max+(temp_y_max-temp_y_min)*0.10*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    set(gca,'linewidth',1.5)
    xlim([0.15 3.85])
    %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.30]);
    ylim([-0.125 1.30])
    set(gca, 'FontSize', 8);%10
    
    xticks([1 2 3]);
    
    xtl = ["Offload";"Error";"Correct"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.22;%0.56,0.4
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9
    set(gca,'xticklabel','');
    
    
    yticks([0 1]);
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');%10
    if if_memeoryPrecision_stimuli0_response1 == 0
        title(sprintf('Stimuli-labeled'),'fontsize',9);
    elseif if_memeoryPrecision_stimuli0_response1 == 1
        title(sprintf('Response-labeled'),'fontsize',9);
    end
    
end


%% Plot: Compare meta-memory  of offload, error
if true
    
    fig = figure('Name','example seq distribution','NumberTitle','off');
    %set(gcf,'Position',[680 450 201.4*1.08*0.95*0.98*0.97 179.5*1.08*0.97*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[680 450 201.4*1.08*0.95*0.98*0.97 179.5*1.08*0.97*1.5*1.3*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[680 450 201.4*1.08*0.95*0.98*0.97*0.65*0.8*1.5 179.5*1.08*0.97*1.5*1.3*1.05*0.43]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[680 450 201.4*1.08*0.95*0.98*0.97*0.65*0.8*1.5 179.5*1.08*0.97*1.5*1.3*1.05*0.43*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[680 450 201.4*1.08*0.95*0.98*0.97*0.65*0.8*1.5*0.84 179.5*1.08*0.97*1.5*1.3*1.05*0.43*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp_p12 = temp_p_metaOffloadError;
    
    temp_1 = temp4;
    temp_2 = temp5;
    
    
    temp_y_min = min([temp_1;temp_2]);
    temp_y_max = max([temp_1;temp_2]);
    
    
    temp_data = [temp_1;temp_2];
    
    g1 = repmat({'A'},length(temp_1),1);
    g2 = repmat({'B'},length(temp_2),1);
    
    temp_label = [g1;g2];
    
    temptemp_color1 = [1 1 1]*0.5;
    temptemp_color2 = repmat(temptemp_color1, 2, 1);
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',[{'A'};{'B'}]);
    h(1).ViolinPlot.FaceAlpha = 0.1;
    h(2).ViolinPlot.FaceAlpha = 0.1;
    
    
    tempTxt = sprintf('');
    if temp_p12 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p12 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p12 < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.20,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([1.1 1.9],temp_y_max+(temp_y_max-temp_y_min)*0.16*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    
    set(gca,'linewidth',1.5)
    xlim([0.45 2.55])
    %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.30]);
    ylim([-0.125 1.30])
    set(gca, 'FontSize', 8);%10
    
    xticks([1 2]);
    
    xtl = ["Offload";"Error"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.22;%0.56,0.4
    %xtext_yp=(temp_y_min)*ones(1,length(xt))-(1-temp_y_min)*0.20;
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9
    
    xtext_yp=(temp_y_min)*ones(1,length(xt))-(1-temp_y_min)*0.32;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',8);%9
    set(gca,'xticklabel','');
    
    
    yticks([0 1]);
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Meta-WM', 'FontSize', 8, 'FontWeight', 'normal');%10
    %     if if_memeoryPrecision_stimuli0_response1 == 0
    %         title(sprintf('Stimuli-labeled'),'fontsize',9);
    %     elseif if_memeoryPrecision_stimuli0_response1 == 1
    %         title(sprintf('Response-labeled'),'fontsize',9);
    %     end
    
end




if if_plot_memoryPrecision_trialLevelEvidence == 1
    
    %% memoryPrecision & meta trialLevelEvidence
    
    trialBoolIndex_memoryPrecisionLow_choiceStimuli = (memoryPrecision_trialLevel_stimuli<=lowThreshold_memoryPrecision) & choiceBoolIndex';
    trialBoolIndex_memoryPrecisionHigh_choiceStimuli = (memoryPrecision_trialLevel_stimuli>lowThreshold_memoryPrecision) & choiceBoolIndex';
    trialBoolIndex_memoryPrecisionLow_choiceResponse = (memoryPrecision_trialLevel_response<=lowThreshold_memoryPrecision) & choiceBoolIndex';
    trialBoolIndex_memoryPrecisionHigh_choiceResponse = (memoryPrecision_trialLevel_response>lowThreshold_memoryPrecision) & choiceBoolIndex';
    
    seqMeta_memoryPrecisionLow_choiceStimuli = nan(sum(numSeq(1:valid_length)),1);
    seqMeta_memoryPrecisionHigh_choiceStimuli = nan(sum(numSeq(1:valid_length)),1);
    seqMeta_memoryPrecisionLow_choiceResponse = nan(sum(numSeq(1:valid_length)),1);
    seqMeta_memoryPrecisionHigh_choiceResponse = nan(sum(numSeq(1:valid_length)),1);
    
    for tempi=1:sum(numSeq(1:valid_length))
        % Stimuli-labeled trial
        tempTrialBoolIndex_targetSeq = seqIndex'==tempi;
        
        temptempBoolIndex = tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionLow_choiceStimuli;
        tempMeta_low = meta_trialLevel(temptempBoolIndex);
        seqMeta_memoryPrecisionLow_choiceStimuli(tempi) = mean(tempMeta_low);
        
        temptempBoolIndex = tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionHigh_choiceStimuli;
        tempMeta_high = meta_trialLevel(temptempBoolIndex);
        seqMeta_memoryPrecisionHigh_choiceStimuli(tempi) = mean(tempMeta_high);
        
        
        % Response-labeled trial
        tempTrialBoolIndex_targetSeq = seqIndex_response'==tempi;
        
        temptempBoolIndex = tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionLow_choiceResponse;
        tempMeta_low = meta_trialLevel(temptempBoolIndex);
        seqMeta_memoryPrecisionLow_choiceResponse(tempi) = mean(tempMeta_low);
        
        temptempBoolIndex = tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionHigh_choiceResponse;
        tempMeta_high = meta_trialLevel(temptempBoolIndex);
        seqMeta_memoryPrecisionHigh_choiceResponse(tempi) = mean(tempMeta_high);
    end
    
    e1 = seqMeta_memoryPrecisionLow_choiceStimuli;
    e2 = seqMeta_memoryPrecisionHigh_choiceStimuli;
    e3 = seqMeta_memoryPrecisionLow_choiceResponse;
    e4 = seqMeta_memoryPrecisionHigh_choiceResponse;
    
    e1_mean = mean(e1,'omitnan');
    e2_mean = mean(e2,'omitnan');
    e3_mean = mean(e3,'omitnan');
    e4_mean = mean(e4,'omitnan');
    
    [~,temp_p_e_lowHigh_stimuli] = ttest(e1,e2);
    [~,temp_p_e_lowHigh_response] = ttest(e3,e4);
    
    %%
    
    seqCorr_precisionMeta_choiceStimuli = nan(sum(numSeq(1:valid_length)),1);
    seqCorr_precisionMeta_choiceResponse = nan(sum(numSeq(1:valid_length)),1);
    
    for tempi=1:sum(numSeq(1:valid_length))
        % Stimuli-labeled trial
        tempTrialBoolIndex_targetSeq = seqIndex'==tempi;
        
        tempPrecision = memoryPrecision_trialLevel_stimuli(tempTrialBoolIndex_targetSeq);
        tempMeta = meta_trialLevel(tempTrialBoolIndex_targetSeq);
        
        temptempBoolIndex = (~isnan(tempPrecision)) & (~isnan(tempMeta));
        seqCorr_precisionMeta_choiceStimuli(tempi) = corr(tempPrecision(temptempBoolIndex),tempMeta(temptempBoolIndex));
        
        
        % Response-labeled trial
        tempTrialBoolIndex_targetSeq = seqIndex_response'==tempi;
        
        tempPrecision = memoryPrecision_trialLevel_response(tempTrialBoolIndex_targetSeq);
        tempMeta = meta_trialLevel(tempTrialBoolIndex_targetSeq);
        
        temptempBoolIndex = (~isnan(tempPrecision)) & (~isnan(tempMeta));
        seqCorr_precisionMeta_choiceResponse(tempi) = corr(tempPrecision(temptempBoolIndex),tempMeta(temptempBoolIndex));
        
    end
    
    
    %% memoryPrecision & meta trialLevelEvidence in stimuli & response-labeled trials
    if true
        fig = figure('Name','memory precsion pdf (mildSeq)','NumberTitle','off');
        %set(gcf,'Position',[10 450 355*1.5*1.05*0.7*0.9 295]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35*2 295*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,2,'TileSpacing','compact','Padding','loose');
        
        
        % Compare meta-memory of lowPrecision & highPrecision, stimuli-labeled trials
        nexttile
        
        temp_p = temp_p_e_lowHigh_stimuli;
        temptempBoolIndex = (~isnan(e1)) & (~isnan(e2));
        temp_1 = e1(temptempBoolIndex);
        temp_2 = e2(temptempBoolIndex);
        
        
        %temp_y_min = min([temp_1;temp_2]);
        %temp_y_max = max([temp_1;temp_2]);
        temp_y_min = 0;
        temp_y_max = 1;
        
        if if_trialEvidenceAllSeq_violinplot0_pairline1 == 0
            temp_data = [temp_1;temp_2];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            
            temp_label = [g1;g2];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 2, 1);
            
            h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'}]);
            h(1).ViolinPlot.FaceAlpha = 0.1;
            h(2).ViolinPlot.FaceAlpha = 0.1;
            
        elseif if_trialEvidenceAllSeq_violinplot0_pairline1 == 1
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            
        end
        
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.65])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
        set(gca, 'FontSize', 8)
        
        xticks([1 2]);
        xtl = ["Low-precision"; "High-precision"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        % xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.405;%0.56,0.4
        xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.19;%0.56,0.4
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
        set(gca,'xticklabel','');
        
        
        % yticks([0 1]);
        yticks([0:0.2:1]);
        
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
        % title(sprintf('Stimuli-labeled\nAll seqs'),'fontsize',9);
        title(sprintf('Stimuli-labeled'),'fontsize',9);
        
        
        
        % Compare meta-memory of lowPrecision & highPrecision, response-labeled trials
        nexttile
        
        temp_p = temp_p_e_lowHigh_response;
        temptempBoolIndex = (~isnan(e3)) & (~isnan(e4));
        temp_1 = e3(temptempBoolIndex);
        temp_2 = e4(temptempBoolIndex);
        
        
        %temp_y_min = min([temp_1;temp_2]);
        %temp_y_max = max([temp_1;temp_2]);
        temp_y_min = 0;
        temp_y_max = 1;
        
        if if_trialEvidenceAllSeq_violinplot0_pairline1 == 0
            temp_data = [temp_1;temp_2];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            
            temp_label = [g1;g2];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 2, 1);
            
            h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'}]);
            h(1).ViolinPlot.FaceAlpha = 0.1;
            h(2).ViolinPlot.FaceAlpha = 0.1;
            
        elseif if_trialEvidenceAllSeq_violinplot0_pairline1 == 1
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            
        end
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.65])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
        set(gca, 'FontSize', 8)
        
        xticks([1 2]);
        xtl = ["Low-precision"; "High-precision"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        % xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.405;%0.56,0.4
        xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.19;%0.56,0.4
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
        set(gca,'xticklabel','');
        
        
        % yticks([0 1]);
        yticks([0:0.2:1]);
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
        % title(sprintf('Response-labeled\nAll seqs'),'fontsize',9);
        title(sprintf('Response-labeled'),'fontsize',9);
        
    end
    
    
    %% asdasds
    temp_1 = seqCorr_precisionMeta_choiceStimuli';
    temp_1_data = 0;
    
    [~,temp_p1] = ttest(temp_1,temp_1_data);
    
    temp_2 = seqCorr_precisionMeta_choiceResponse';
    temp_2_chanceLevel = 0;
    
    [~,temp_p2] = ttest(temp_2,temp_2_chanceLevel);
    
    if true
        fig = figure('Name','Within seqs correlation','NumberTitle','off');
        %set(gcf,'Position',[10 450 355*1.5*1.05*0.7*0.9 295]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35*2 295*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35*2*1.02 295*0.78*0.94*1.14]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35*2*1.02 295*0.78*0.94*1.14*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35*2*1.02*1.35 295*0.78*0.94*1.14*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35*2*1.02*1.35*0.87 295*0.78*0.94*1.14*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35*2*1.02*1.35*0.87 295*0.78*0.94*1.14*1.5*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35*2*1.02*1.35*0.87*0.7 295*0.78*0.94*1.14*1.5*0.9*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,2,'TileSpacing','compact','Padding','loose');
        
        %stimuli-labeled trials
        nexttile
        
        %temp_y_min = min([temp_1 temp_1_data]);
        %temp_y_max = max(temp_1);
        %temp_y_min = -0.1;%-0.5
        %temp_y_max = 0.65;%1
        temp_y_min = -0.5;%-0.5
        temp_y_max = 1;%1
        
        
        temp_data = temp_1';
        temp_label = repmat({'A'},length(temp_1),1);
        
        temptemp_color1 = [1 1 1]*0.5;
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',{'A'});
        h(1).ViolinPlot.FaceAlpha = 0.1;
        hold on
        
        
        plot([0 2],temp_1_data*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        
        tempTxt = sprintf('');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        
        %xlim([0.5 1.5]);%[0 2]
        %xlim([0.35 1.65]);%[0 2]
        xlim([0 2]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
        set(gca, 'FontSize', 8) %14,12
        %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
        %set(gca,'YTick',0:1,'FontSize', 10);%12
        set(gca,'box','off');% 取消右、上边框
        
        
        xtl = [""]; %#ok<*NBRAK>
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
        
        set(gca,'xticklabel','');
        
        
        ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
        title(sprintf('Quality VS. Meta'),'fontsize',9);
        %subtitle(sprintf('Seq-level, stimuli-labeled'),'fontsize',7);
        subtitle(sprintf('Within seqs, stimuli-labeled'),'fontsize',7);
        
        
        %response-labeled trials
        nexttile
        
        %temp_y_min = min([temp_2 temp_2_chanceLevel]);
        %temp_y_max = max(temp_2);
        temp_y_min = -0.5;
        temp_y_max = 1;
        
        temp_data = temp_2';
        temp_label = repmat({'A'},length(temp_2),1);
        
        temptemp_color1 = [1 1 1]*0.5;
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',{'A'});
        h(1).ViolinPlot.FaceAlpha = 0.1;
        hold on
        
        
        plot([0 2],temp_2_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        
        tempTxt = sprintf('');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        
        xlim([0 2]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
        set(gca, 'FontSize', 8) %14,12
        %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
        %set(gca,'YTick',0:1,'FontSize', 10);%12
        set(gca,'box','off');% 取消右、上边框
        
        
        xtl = [""]; %#ok<*NBRAK>
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
        
        set(gca,'xticklabel','');
        
        
        ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
        title(sprintf('Quality VS. Meta'),'fontsize',9);
        subtitle(sprintf('Seq-level, response-labeled'),'fontsize',7);
        
    end
    
end




%% Try to expand trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError

temp_flag_expandMismatch = 0;

% if true
% if if_plot_memoryPrecision_trialLevelEvidence == 1
if exist('seqDistri_trialLevel_response','var') == 1
    
    temp_flag_expandMismatch = 1;
    
    %% Step1
    trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError;
    
    
    temp_precision = memoryPrecision_trialLevel(trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError);
    temp_precision_response = memoryPrecision_trialLevel_response(trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError);
    
    
    seqDistri_trialLevel;
    %seqDistri_trialLevel_stimuli;
    %seqDistri_trialLevel_response;
    
    %temp_seqDistri = seqDistri_trialLevel(trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError,:);
    temp_seqDistri = seqDistri_trialLevel_response(trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError,:);
    
    seqIndex;
    seqIndex_response;
    
    temp_seqIndex = seqIndex(trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError);
    temp_seqIndex_response = seqIndex_response(trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError);
    
    
    temp_trialNum = sum(trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError);
    temp_seqProb_stimuli = nan(temp_trialNum,1);
    temp_seqProb_response = nan(temp_trialNum,1);
    for tempi=1:temp_trialNum
        temp_seqProb_stimuli(tempi) = temp_seqDistri(tempi,temp_seqIndex(tempi));
        if temp_seqIndex_response(tempi) <= 41
            temp_seqProb_response(tempi) = temp_seqDistri(tempi,temp_seqIndex_response(tempi));
        end
    end
    
    mean(temp_seqProb_stimuli,'omitnan');
    mean(temp_seqProb_response,'omitnan');
    
    temp_seqProb_response;
    
    temp_seqProb_max = max(temp_seqDistri,[],2);
    
    temp_trialBoolIndex_stimuliMax = temp_seqProb_max==temp_seqProb_stimuli;
    temp_trialBoolIndex_responseMax = temp_seqProb_max==temp_seqProb_response;
    temp_trialBoolIndex_noise = ~(temp_trialBoolIndex_stimuliMax|temp_trialBoolIndex_responseMax);
    
    temp_trialNum_stimuliMax = sum(temp_trialBoolIndex_stimuliMax);
    temp_trialNum_responseMax = sum(temp_trialBoolIndex_responseMax);
    temp_trialNum_noise = sum(temp_trialBoolIndex_noise);
    
    fprintf('TrialNum of broad mismatch: [total, stimuliMax, responseMax, noise] = [%d, %d, %d, %d].\n',...
        temp_trialNum,temp_trialNum_stimuliMax,temp_trialNum_responseMax,temp_trialNum_noise);
    
    temp_trialBoolIndex_responseMax;
    
    temp_trialBoolIndex_CME_responseMax = temp_trialBoolIndex_responseMax;
    temp_trialBoolIndex_CME_others = ~temp_trialBoolIndex_responseMax;
    
    
    %% Step2
    temp_trialBoolIndex_CME_responseMax;
    temp_trialBoolIndex_CME_others;
    
    trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError;
    trialIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError = find(trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError == true);
    
    trialIndex_memoryPrecisionHigh_metaHigh_CME_responseMax = trialIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError(temp_trialBoolIndex_CME_responseMax);
    trialIndex_memoryPrecisionHigh_metaHigh_CME_others = trialIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError(temp_trialBoolIndex_CME_others);
    
    trialBoolIndex_memoryPrecisionHigh_metaHigh_CME_responseMax = false(size(trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError));
    trialBoolIndex_memoryPrecisionHigh_metaHigh_CME_others = false(size(trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError));
    
    trialBoolIndex_memoryPrecisionHigh_metaHigh_CME_responseMax(trialIndex_memoryPrecisionHigh_metaHigh_CME_responseMax) = true;
    trialBoolIndex_memoryPrecisionHigh_metaHigh_CME_others(trialIndex_memoryPrecisionHigh_metaHigh_CME_others) = true;
    
    
    trialBoolIndex_memoryPrecisionHigh_metaHigh_CME_responseMax;
    trialBoolIndex_memoryPrecisionHigh_metaHigh_CME_others;
    
    %%
    choiceMemoryErrorBoolIndex;
    
    trialBoolIndex_memoryPrecisionLow_metaHigh_choiceMemoryError;% mismatch

    trialBoolIndex_memoryPrecisionHigh_metaHigh_choiceMemoryError;% broad match
    trialBoolIndex_memoryPrecisionHigh_metaHigh_CME_responseMax;% narrow match
    trialBoolIndex_memoryPrecisionHigh_metaHigh_CME_others;% noise from match
    
    trialBoolIndex_memoryPrecisionHigh_metaLow_choiceMemoryError;% noise
    trialBoolIndex_memoryPrecisionLow_metaLow_choiceMemoryError;% noise
    
    temptempBoolIndex_1 = trialBoolIndex_memoryPrecisionLow_metaHigh_choiceMemoryError;% mismatch
    temptempBoolIndex_2 = trialBoolIndex_memoryPrecisionHigh_metaHigh_CME_responseMax;% match
    temptempBoolIndex_3 = trialBoolIndex_memoryPrecisionHigh_metaHigh_CME_others...% noise
        | trialBoolIndex_memoryPrecisionHigh_metaLow_choiceMemoryError...
        |trialBoolIndex_memoryPrecisionLow_metaLow_choiceMemoryError;
    
    %% Plot expanded mismatch trials
    if true
        fig = figure('Name','Expanded mismatch trials','NumberTitle','off');
        set(gcf,'Position',[1010 390 196.4*0.905*0.9*0.89 227.5*0.905*0.9*0.9*0.92*0.975]);
        
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        x1 = memoryPrecision_trialLevel(temptempBoolIndex_1);
        y1 = meta_trialLevel(temptempBoolIndex_1);
        
        x2 = memoryPrecision_trialLevel(temptempBoolIndex_2);
        y2 = meta_trialLevel(temptempBoolIndex_2);
        
        x3 = memoryPrecision_trialLevel(temptempBoolIndex_3);
        y3 = meta_trialLevel(temptempBoolIndex_3);
        
        [temp_xmin,temp_xmax] = bounds([x1;x2;x3]);
        [temp_ymin,temp_ymax] = bounds([y1;y2;y3]);
        
        temp_xmin = -0.05;
        temp_xmax = 1.05;
        temp_ymin = -0.05;
        temp_ymax = 1.05;
        
        plot([temp_xmin temp_xmax],[1 1].*lowThreshold_meta, '--', 'LineWidth', 1, 'Color', [1 1 1]*0.5);
        hold on
        plot([1 1].*lowThreshold_memoryPrecision,[temp_ymin temp_ymax], '--', 'LineWidth', 1, 'Color',[1 1 1]*0.5);
        hold on
        
        temp_size = 15;%3,5,7,35,20,15
        temp_MarkerFaceAlpha = 0.45; %0.5,0.3,0.25,0.125,0.175,0.625
        temp_MarkerEdgeAlpha = 0; %0.7,0
        
        temp_h1 = scatter(x1, y1, ...
            temp_size, 'filled', 'MarkerFaceColor', [1 1 1]*1, ...%[1 1 1]*0
            'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 1,'MarkerEdgeColor',[0 0 0],'linewidth',1);
        hold on
        
        temp_h2 = scatter(x2, y2, ...
            temp_size, 'filled', 'MarkerFaceColor', color_choiceMemory, ... %color_choiceMemoryError
            'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
        hold on

        
        temp_h3 = scatter(x3, y3, ...
            temp_size, 'filled', 'MarkerFaceColor', [1 1 1].*0.5, ...
            'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
        hold on
        
        %le = legend('Mismatch','Match','Others','Location','northwest','fontsize',6.5,'NumColumns',2);
        %le.ItemTokenSize = ones(1,3)*5;
        %legend('boxoff');
        
        
        
        %xlim([-0.05 1.05]);
        %ylim([-0.05 1.05]);
        xlim([temp_xmin temp_xmax]);
        ylim([temp_ymin temp_ymax]);
        
        xticks([0 1]);
        yticks([0 1]);
        
        
        set(gca,'linewidth',1.5)
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
        ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');
        
        %title(sprintf('asd'), 'FontSize', 9, 'FontWeight', 'normal');
        %subtitle(sprintf('asd'), 'FontSize', 8, 'FontWeight', 'normal');
        
        
    end
    
end

%%
temp_flag_expandMismatch;

% trialBoolIndex_memoryPrecisionHigh_metaHigh_CME_responseMax;% narrow match


%% To test validity of decoders in choice-memory error trials
x_raw = memoryPrecision_trialLevel(choiceBoolIndex);
y_raw = meta_trialLevel(choiceBoolIndex);
temp_memoryPrecision_choice = x_raw(~isnan(x_raw));
temp_meta_choice = y_raw(~isnan(x_raw));

temp_choiceMemoryErrorBoolIndex_valid = choiceMemoryErrorBoolIndex & ~isnan(memoryPrecision_trialLevel)';

temp_memoryPrecision_over = memoryPrecision_trialLevel(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh);
temp_meta_over = meta_trialLevel(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh);

temp_overTrialNum = length(temp_memoryPrecision_over);

temp_resampledLoopNum = 1000;

temp_resampledTrialNum = sum(temp_choiceMemoryErrorBoolIndex_valid);


if temp_flag_expandMismatch == 1
    temp_resampledTrialNum = temp_resampledTrialNum - sum(trialBoolIndex_memoryPrecisionHigh_metaHigh_CME_responseMax);
end

temp_trialNum = length(temp_memoryPrecision_choice);

temp_memoryPrecision_choice_resampled = nan(temp_resampledTrialNum,temp_resampledLoopNum);
temp_meta_choice_resampled = nan(temp_resampledTrialNum,temp_resampledLoopNum);

for tempi=1:temp_resampledLoopNum
    temptempIndex = sort(randperm(temp_trialNum,temp_resampledTrialNum));
    temptempBoolIndex = false(1,temp_trialNum);
    temptempBoolIndex(temptempIndex) = true;
    
    temp_memoryPrecision_choice_resampled(:,tempi) = temp_memoryPrecision_choice(temptempBoolIndex);
    temp_meta_choice_resampled(:,tempi) = temp_meta_choice(temptempBoolIndex);
end


lowThreshold_memoryPrecision;
lowThreshold_meta;

temp_memoryPrecisionLowBool_choice_resampled = temp_memoryPrecision_choice_resampled<lowThreshold_memoryPrecision;
temp_metaHighBool_choice_resampled = temp_meta_choice_resampled>lowThreshold_meta;

temp_memoryPrecisionLowBool_choice_resampled;


temp_overBool_resampled = temp_memoryPrecisionLowBool_choice_resampled & temp_metaHighBool_choice_resampled;
temp_overTrialNum_resampled = sum(temp_overBool_resampled,1);


temp_overTrialProp_resampled = temp_overTrialNum_resampled./temp_resampledTrialNum;
temp_overTrialProp = temp_overTrialNum./temp_resampledTrialNum;

% [~,temp_p] = ttest(temp_overTrialNum_resampled,temp_overTrialNum);

prctile(temp_overTrialNum_resampled,99);
prctile(temp_overTrialNum_resampled,99)./temp_resampledTrialNum;


temp_1 = temp_overTrialNum_resampled./temp_resampledTrialNum;
temp_1_prctile = prctile(temp_overTrialNum_resampled,99)./temp_resampledTrialNum;

temp_1_data = temp_overTrialNum./temp_resampledTrialNum;

[~,temp_p] = ttest(temp_1,temp_1_data);

if true
    fig = figure('Name','Test validity','NumberTitle','off');
    %set(gcf,'Position',[650+0 50+0 158.4*0.94*1.03*1.3*1.13 167.5*0.94*0.95*1.2*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[650+0 50+0 158.4*0.94*1.03*1.3*1.13*0.9 167.5*0.94*0.95*1.2*1.05*1.11]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[650+0 50+0 158.4*0.94*1.03*1.3*1.13*0.9*0.86*0.95 167.5*0.94*0.95*1.2*1.05*1.11]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[650+0 50+0 158.4*0.94*1.03*1.3*1.13*0.9*0.86*0.95*0.86 167.5*0.94*0.95*1.2*1.05*1.11]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[650+0 50+0 158.4*0.94*1.03*1.3*1.13*0.9*0.86*0.95*0.86*1.3*1.2 167.5*0.94*0.95*1.2*1.05*1.11]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[650+0 50+0 158.4*0.94*1.03*1.3*1.13*0.9*0.86*0.95*0.86*1.3*1.2*0.59*1.1 167.5*0.94*0.95*1.2*1.05*1.11*0.86]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    
    x1 = temp_1;
    
    x1_prctile = temp_1_prctile;
    
    x2 = temp_1_data;
    
    h_NumBins = 8;%10
    
    h_handle = [];
    
    x = x1;
    %h1 = histogram(x,'FaceAlpha',1,'Normalization','probability','DisplayStyle','stairs','LineWidth',1.5);
    h1 = histogram(x,'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1.5);
    hold on
    h1.NumBins = h_NumBins;
    h1.EdgeColor = [1 1 1].*0.5;
    
    y1 = h1.Values;
    
    h_handle = h1;
    
    %x_min = 0;
    %x_max = 1;
    
    [x_min,x_max] = bounds([x1,x2]);
    [y_min,y_max] = bounds(y1);
    
    %h_handle = [h_handle plot(x1_prctile*[1 1],[y_min y_max],'--','color',[0.5 0.5 0.5],'linewidth',1)];
    plot(x1_prctile*[1 1],[y_min y_max],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    h_handle = [h_handle plot(x2*[1 1],[y_min y_max],'-','color',[0.25 0.25 0.25],'linewidth',2)];
    hold on
    
    %le = legend('Resampled','99%','Data','Location','northwest','fontsize',6.5,'NumColumns',3);
    %le = legend('Shuffle','99%','Data','Location','northwest','fontsize',6.5,'NumColumns',2);
    %le = legend(h_handle,'Shuffle','Data','Location','northwest','fontsize',6.5,'NumColumns',2);
    le = legend(h_handle,'Shuffle','Data','Location','northwest','fontsize',6.5,'NumColumns',2);    
    le.ItemTokenSize = ones(1,3)*10;%13
    legend('boxoff');
    
    
    set(gca,'linewidth',1.5)
    xlim([x_min-(x_max-x_min)*0.2 x_max+(x_max-x_min)*0.2]);
    %ylim([y_min y_max+(y_max-y_min)*0.2]);%0.1
    ylim([y_min y_max+(y_max-y_min)*0.32]);%0.1
    %xticks([0 1]);
    if if_twoThreshold_median0_optimal1 == 0
        
    elseif if_twoThreshold_median0_optimal1 == 1
        xticks([0.2 0.5]);
    end
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Trial proportion', 'FontSize', 9, 'FontWeight', 'normal');
    %ylabel('Probability', 'FontSize', 9, 'FontWeight', 'normal');
    ylabel('Pdf', 'FontSize', 9, 'FontWeight', 'normal');
    
    %temp_title = title(sprintf('Low strength mismatch trials'),'fontsize',9);
    %temp_title = title(sprintf('Low-strength\nmismatch trials'),'fontsize',9);
    %temp_title = title(sprintf('Low-strength mismatch trials'),'fontsize',9);    
    temp_title = title(sprintf('Mismatch'),'fontsize',9);        
    temp_title.Interpreter = 'none';
    
    
end



%% End