function mappingTwoSessionFall_v6(path_sessionA, path_sessionB,options_dff,options_category)
% med distance < min(20 pix, 2*sqrt(npix) of roi)
% npix is similar 
% aspect_ratio is similar
% contour is overlapped (optional)

% get a reference roi A, fetch top x (up to) 5 rois (roi Bs) which are subject to med distance and similar npix with reference roi A
% then find if any roi Bs' contour is overlapped with roi A's
% clear all %#ok<*CLALL>

% path_sessionA = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230122A_Ding_Site09B\Result20230202T101726\plane0';
% % path_sessionB = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230123A_Ding_Site09B_sameFOV0122\Result20230202T231115\plane0';
% path_sessionB = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230123A_Ding_Site09B_sameFOV0122\Result20230205T170425\plane0';

fileName_Fall = 'Fall.mat';
fileName_iscell = 'iscell.npy';
fullFileName_Fall_sessionA = [path_sessionA,'\',fileName_Fall];
fullFileName_iscell_sessionA = [path_sessionA,'\',fileName_iscell];
fullFileName_Fall_sessionB = [path_sessionB,'\',fileName_Fall];
fullFileName_iscell_sessionB = [path_sessionB,'\',fileName_iscell];

iscell_sessionA = readNPY(fullFileName_iscell_sessionA);
iscell_sessionB = readNPY(fullFileName_iscell_sessionB);


Fall_sessionA = load(fullFileName_Fall_sessionA);
Fall_sessionB = load(fullFileName_Fall_sessionB);

stat_sessionA = Fall_sessionA.stat;
stat_sessionB = Fall_sessionB.stat;

threshold_cellProb = 0.01;

med_sessionA = zeros(2,length(stat_sessionA));
npix_sessionA = zeros(1,length(stat_sessionA));
aspect_ratio_sessionA = zeros(1,length(stat_sessionA));
med_sessionB = zeros(2,length(stat_sessionB));
npix_sessionB = zeros(1,length(stat_sessionB));
aspect_ratio_sessionB = zeros(1,length(stat_sessionB));

for tempi=1:length(stat_sessionA)
    med_sessionA(:,tempi) = stat_sessionA{tempi}.med;
    npix_sessionA(:,tempi) = stat_sessionA{tempi}.npix;
    aspect_ratio_sessionA(:,tempi) = stat_sessionA{tempi}.aspect_ratio;
    
%     if iscell_sessionA(tempi,2) <= threshold_cellProb
%         med_sessionA(:,tempi) = [-500,-500];
%     end      
end
for tempi=1:length(stat_sessionB)
    med_sessionB(:,tempi) = stat_sessionB{tempi}.med;
    npix_sessionB(:,tempi) = stat_sessionB{tempi}.npix;
    aspect_ratio_sessionB(:,tempi) = stat_sessionB{tempi}.aspect_ratio;
    
%     if iscell_sessionB(tempi,2) <= threshold_cellProb
%         med_sessionB(:,tempi) = [-500,-500];
%     end         
end

threshold_med_sqrtRatio = 2;%2
threshold_medDisMax = 20;
weight = struct;
weight.med = 1;%1
weight.npix = 0.5;%0.5
weight.aspect_ratio = 0.5;%0.5
candidateNum = 5;
%% Calculate candidate_sessionA
candidate_sessionA = struct;
candidate_sessionA.index = zeros(candidateNum,length(stat_sessionA));
candidate_sessionA.score = zeros(candidateNum,length(stat_sessionA));
%candidate_sessionA.med_distance_top = zeros(candidateNum,length(stat_sessionA));
for tempi=1:length(stat_sessionA)    
    if iscell_sessionA(tempi,2) <= threshold_cellProb
        continue
    end    
    
    tempA_med = med_sessionA(:,tempi);
    tempA_npix = npix_sessionA(tempi);
    tempA_aspect_ratio = aspect_ratio_sessionA(tempi);
      
    temp_med_distance = sqrt(sum((tempA_med-med_sessionB).^2,1));
    [~,I_med] = sort(temp_med_distance,'ascend');    
    temp_candidateIndex = I_med(1:candidateNum);
    temp_med_distance_top = temp_med_distance(temp_candidateIndex) + 0.1;    
    %temp_candidateScore = (candidateNum:-1:1) * weight.med;
    temp_candidateScore = (max(temp_med_distance_top)./temp_med_distance_top) * weight.med;
    
    temp_npix_distance = abs(tempA_npix-npix_sessionB(temp_candidateIndex));
    [~,I_npix] = sort(temp_npix_distance,'ascend');
    temp_candidateScore(I_npix) = temp_candidateScore(I_npix) + (candidateNum:-1:1) * weight.npix;
    
    temp_aspect_ratio_distance = abs(tempA_aspect_ratio-aspect_ratio_sessionB(temp_candidateIndex));
    [~,I_aspect_ratio] = sort(temp_aspect_ratio_distance,'ascend');
    temp_candidateScore(I_aspect_ratio) = temp_candidateScore(I_aspect_ratio) + (candidateNum:-1:1) * weight.aspect_ratio;
    
    [~,I_scoreMax] = max(temp_candidateScore);
    if temp_med_distance_top(I_scoreMax) > min(threshold_medDisMax,threshold_med_sqrtRatio*sqrt(tempA_npix))
        continue
    end
    
    temp_bad_med_index = temp_med_distance_top > min(threshold_medDisMax,threshold_med_sqrtRatio*sqrt(tempA_npix));
    if sum(temp_bad_med_index) > 0
        a = 1;
        temp_candidateScore(temp_bad_med_index) = 0;
    end
    if iscell_sessionB(temp_candidateIndex(I_scoreMax),2) <= threshold_cellProb
        continue
    end   
    
    candidate_sessionA.index(:,tempi) = temp_candidateIndex;
    candidate_sessionA.score(:,tempi) = temp_candidateScore;
    %candidate_sessionA.med_distance_top(:,tempi) = temp_med_distance_top;
    
    %if temp_med_distance_top(I_scoreMax) > 2*sqrt(tempA_npix)
    %    a = 1;
    %end
    
    if tempi == 614 + 1
        a = 1; %#ok<*NASGU>
    end
    a = 1;
end
validCandidateCount_sessionA = sum(sum(candidate_sessionA.score,1)>0);
fprintf('validCandidateCount_sessionA = %d.\n',validCandidateCount_sessionA);

%% Calculate candidate_sessionB
candidate_sessionB = struct;
candidate_sessionB.index = zeros(candidateNum,length(stat_sessionB));
candidate_sessionB.score = zeros(candidateNum,length(stat_sessionB));
%candidate_sessionB.med_distance_top = zeros(candidateNum,length(stat_sessionB));
for tempi=1:length(stat_sessionB)    
    if iscell_sessionB(tempi,2) <= threshold_cellProb
        continue
    end    
    
    tempB_med = med_sessionB(:,tempi);
    tempB_npix = npix_sessionB(tempi);
    tempB_aspect_ratio = aspect_ratio_sessionB(tempi);
      
    temp_med_distance = sqrt(sum((tempB_med-med_sessionA).^2,1));
    [~,I_med] = sort(temp_med_distance,'ascend');    
    temp_candidateIndex = I_med(1:candidateNum);
    temp_med_distance_top = temp_med_distance(temp_candidateIndex) + 0.1;    
    %temp_candidateScore = (candidateNum:-1:1) * weight.med;
    temp_candidateScore = (max(temp_med_distance_top)./temp_med_distance_top) * weight.med;
    
    temp_npix_distance = abs(tempB_npix-npix_sessionA(temp_candidateIndex));
    [~,I_npix] = sort(temp_npix_distance,'ascend');
    temp_candidateScore(I_npix) = temp_candidateScore(I_npix) + (candidateNum:-1:1) * weight.npix;
    
    temp_aspect_ratio_distance = abs(tempB_aspect_ratio-aspect_ratio_sessionA(temp_candidateIndex));
    [~,I_aspect_ratio] = sort(temp_aspect_ratio_distance,'ascend');
    temp_candidateScore(I_aspect_ratio) = temp_candidateScore(I_aspect_ratio) + (candidateNum:-1:1) * weight.aspect_ratio;
    
    [~,I_scoreMax] = max(temp_candidateScore);
    if temp_med_distance_top(I_scoreMax) > min(threshold_medDisMax,threshold_med_sqrtRatio*sqrt(tempB_npix))
        continue
    end
    
    temp_bad_med_index = temp_med_distance_top > min(threshold_medDisMax,threshold_med_sqrtRatio*sqrt(tempB_npix));
    if sum(temp_bad_med_index) > 0
        a = 1;
        temp_candidateScore(temp_bad_med_index) = 0;
    end
    
    if iscell_sessionA(temp_candidateIndex(I_scoreMax),2) <= threshold_cellProb
        continue
    end   
    
    candidate_sessionB.index(:,tempi) = temp_candidateIndex;
    candidate_sessionB.score(:,tempi) = temp_candidateScore;   
    %candidate_sessionB.med_distance_top(:,tempi) = temp_med_distance_top;
    
    if tempi == 731
        a = 1; %#ok<*NASGU>
    end
end
validCandidateCount_sessionB = sum(sum(candidate_sessionB.score,1)>0);
fprintf('validCandidateCount_sessionB = %d.\n',validCandidateCount_sessionB);

confusionMatrix_AtoB = zeros(length(stat_sessionA),length(stat_sessionB));
for tempi=1:length(stat_sessionA)
    temp_A_index = candidate_sessionA.index(:,tempi);
    if sum(temp_A_index) == 0
        continue
    end
    temp_A_score = candidate_sessionA.score(:,tempi);
    
    %if tempi == 8
    %    a = 1;
    %end    
    [~,I_score] = sort(temp_A_score,'descend');
    confusionMatrix_AtoB(tempi,temp_A_index(I_score)) = (candidateNum:-1:1);        
    confusionMatrix_AtoB(tempi,temp_A_index(temp_A_score<eps)) = 0; 
    
end
confusionMatrix_BtoA = zeros(length(stat_sessionB),length(stat_sessionA));
for tempi=1:length(stat_sessionB)
    temp_B_index = candidate_sessionB.index(:,tempi);
    if sum(temp_B_index) == 0
        continue
    end    
    temp_B_score = candidate_sessionB.score(:,tempi);
    [~,I_score] = sort(temp_B_score,'descend');
    confusionMatrix_BtoA(tempi,temp_B_index(I_score)) = (candidateNum:-1:1);  
    confusionMatrix_BtoA(tempi,temp_B_index(temp_B_score<eps)) = 0;
end
confusionMatrix_BtoA_trans = confusionMatrix_BtoA';

boolMatrix_ABpair = (confusionMatrix_AtoB > 0) & (confusionMatrix_BtoA_trans>0);
a1 = sum(boolMatrix_ABpair,2);
a2 = sum(boolMatrix_ABpair,1);

confusionMatrix_AtoB_shrink = confusionMatrix_AtoB .* boolMatrix_ABpair;
confusionMatrix_BtoA_trans_shrink = confusionMatrix_BtoA_trans .* boolMatrix_ABpair;

% sum(sum(confusionMatrix_AtoB_shrink_maxFilter>0,2)>1)
% sum(sum(confusionMatrix_AtoB_shrink_maxFilter>0,1)>1)


confusionMatrix_AtoB_shrink_raw = confusionMatrix_AtoB_shrink;
confusionMatrix_BtoA_trans_shrink_raw = confusionMatrix_BtoA_trans_shrink;

conflictBoolIndex_AtoB = 1;
while sum(conflictBoolIndex_AtoB) > 0
    temp_max = max(confusionMatrix_AtoB_shrink,[],2);
    confusionMatrix_AtoB_shrink_maxFilter = confusionMatrix_AtoB_shrink;
    confusionMatrix_AtoB_shrink_maxFilter(confusionMatrix_AtoB_shrink < temp_max) = 0;
    conflictBoolIndex_AtoB = sum(confusionMatrix_AtoB_shrink_maxFilter>0,1)>1;
    
    for tempi=1:length(conflictBoolIndex_AtoB)
        if conflictBoolIndex_AtoB(tempi)==0
            continue
        end
        %a = 1;
        temp_conflictIndex_AtoB = find(confusionMatrix_AtoB_shrink_maxFilter(:,tempi));
                
        [~,temp_temp_maxIndex]= max(confusionMatrix_BtoA_trans_shrink(temp_conflictIndex_AtoB,tempi));
        temp_maxIndex = temp_conflictIndex_AtoB(temp_temp_maxIndex);
        temp_otherIndex = temp_conflictIndex_AtoB;
        temp_otherIndex(temp_temp_maxIndex) = [];
        
        confusionMatrix_AtoB_shrink(temp_otherIndex,tempi) = 0;
        break
    end
end
valid_pairCount = sum(confusionMatrix_AtoB_shrink_maxFilter>0,'all');
valid_pairIndex_AtoB = -1*ones(1,length(stat_sessionA));
for tempi=1:length(stat_sessionA)
    temp_findIndex = find(confusionMatrix_AtoB_shrink_maxFilter(tempi,:)>0);
    if isempty(temp_findIndex) == 0
        valid_pairIndex_AtoB(tempi) = temp_findIndex;
    end
end
valid_pairCount_2 = length(unique(valid_pairIndex_AtoB)) - 1;

valid_pairIndex_AandB = zeros(valid_pairCount,2);
valid_pairIndex_AandB(:,1) = find(valid_pairIndex_AtoB>-1);
valid_pairIndex_AandB(:,2) = valid_pairIndex_AtoB(valid_pairIndex_AandB(:,1));

fprintf('valid_pairCount = %d.\n',valid_pairCount);

valid_pairIndex_AandB; %#ok<*VUNUS>
iscell_sessionA(iscell_sessionA == 2) = 1;
iscell_sessionB(iscell_sessionB == 2) = 1;
iscell_sessionA(iscell_sessionA == 1.5) = 1.1;
iscell_sessionB(iscell_sessionB == 1.5) = 1.1;
isdendrite_sessionA = (iscell_sessionA(valid_pairIndex_AandB(:,1),2)==1.1);
isdendrite_sessionB = (iscell_sessionB(valid_pairIndex_AandB(:,2),2)==1.1);

% isdendrite = isdendrite_sessionA | isdendrite_sessionB;
isdendrite = isdendrite_sessionA & isdendrite_sessionB;
iscell_sessionA(valid_pairIndex_AandB(:,1),2) = 2;
iscell_sessionB(valid_pairIndex_AandB(:,2),2) = 2;
iscell_sessionA(valid_pairIndex_AandB(isdendrite,1),2) = 1.5;
iscell_sessionB(valid_pairIndex_AandB(isdendrite,2),2) = 1.5;
% iscell_sessionA(iscell_sessionA == 1.1) = 0;
% iscell_sessionB(iscell_sessionB == 1.1) = 0;
iscell_sessionA(iscell_sessionA == 1.1) = 0.1;
iscell_sessionB(iscell_sessionB == 1.1) = 0.1;

% iscell_sessionA(iscell_sessionA > 1) = 1;
% iscell_sessionB(iscell_sessionB > 1) = 1;
% iscell_sessionA(valid_pairIndex_AandB(:,1),2) = 2;
% iscell_sessionB(valid_pairIndex_AandB(:,2),2) = 2;

iscell_sessionA(:,1) = 0;
iscell_sessionB(:,1) = 0;
iscell_sessionA(valid_pairIndex_AandB(:,1),1) = 1;
iscell_sessionB(valid_pairIndex_AandB(:,2),1) = 1;

writeNPY(iscell_sessionA, fullFileName_iscell_sessionA);
writeNPY(iscell_sessionB, fullFileName_iscell_sessionB);

mappingIndexTwoSessionROI = struct;
mappingIndexTwoSessionROI.sessionIndex = 'A';
mappingIndexTwoSessionROI.valid_pairIndex_AandB = valid_pairIndex_AandB;
save([path_sessionA,'\','mappingIndexTwoSessionROI.mat'],'mappingIndexTwoSessionROI');

mappingIndexTwoSessionROI = struct;
mappingIndexTwoSessionROI.sessionIndex = 'B';
mappingIndexTwoSessionROI.valid_pairIndex_AandB = valid_pairIndex_AandB;
save([path_sessionB,'\','mappingIndexTwoSessionROI.mat'],'mappingIndexTwoSessionROI');

if options_category.if_categorySomaDendrite_sessionA == 1
    fprintf('Now categorySomaDendrite for sessionA.\n');
    options_category.fun(path_sessionA);    
end

%options_dff, sessionA
if options_dff.if_dff_sessionA == 1
    fprintf('Now calculating dff for sessionA.\n');
    options_dff.if_plot_dff = 0;
    options_dff.fun(path_sessionA,options_dff);
end





