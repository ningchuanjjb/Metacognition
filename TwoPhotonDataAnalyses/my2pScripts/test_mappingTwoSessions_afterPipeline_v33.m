% if_compute = 0;
% if_plot = 1;


% plotRoiNum = 100;
% temp_candidateIndex_suite2p;

if_AB_uniquePair = 0;


threshold_med_sqrtRatio = 4;%2-->3-->4
threshold_medDisMax = 10;%20-->30-->15-->10
weight = struct;
weight.med = 1;%1
weight.npix = 0.5;%0.5
weight.xy_ratio = 0.5;%0.5
weight.sseMatrix = 2;%2-->3-->2
candidateNum = 20;%5-->50-->20

%% Compute displacement field
temp_if_max0_min1 = 0;
template_path = autoGetFileName_general('maxProjection*.tif', input_path_A,temp_if_max0_min1);
template_A = double(loadtiff(template_path));
template_path = autoGetFileName_general('maxProjection*.tif', input_path_B,temp_if_max0_min1);
template_B = double(loadtiff(template_path));



temp_targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
singleImage_nonrigid_jjb_Name_v = autoGetFunName('singleImage_nonrigid_jjb', [temp_targetPATH '\functions']);
fun_singleImage_nonrigid_jjb= str2func(singleImage_nonrigid_jjb_Name_v);


% compute
% if if_compute == 1
if true
% if false
    template_T = template_A; %#ok<*UNRCH>
    template = template_B;
    
    [temp_min_template,temp_max_template] = bounds(template,'all');
    if if_compute == 1
        [D_raw,~] = fun_singleImage_nonrigid_jjb(template,template_T); %#ok<*PFBNS>
        D = D_raw;
    end
    
    %     D = zeros(512,512,2);
    %     %D(:,:,1) = 10;% x left move 10 pix
    %     %D(:,:,2) = 2;% y up move 2 pix
    %
    %     D(1:256,:,1) = 15;% x left move 10 pix
    %     D(1:256,:,2) = 3;% y up move 5 pix
    
    
    template_corrected = imwarp(template,gather(D),'cubic');%'nearest' would lose sub-pixel displayment!
    overMax_index = template_corrected > temp_max_template;
    underMin_index = template_corrected < temp_min_template;
    template_corrected = template_corrected.*(~overMax_index) + temp_max_template.*overMax_index;
    template_corrected = template_corrected.*(~underMin_index) + temp_min_template.*underMin_index;
    
    template_B_corrected = template_corrected;
    D_BtoA = D;
end
% end
%
template_B_corrected;
D_BtoA;



%% Load anatomy info
fileName_Fall = 'Fall.mat';
fileName_iscell = 'iscell.npy';

% Load A anatomy info
path_plane_A = [input_path_A,'\plane0'];
fullFileName_Fall_A = [path_plane_A,'\',fileName_Fall];
fullFileName_iscell_A = [path_plane_A,'\',fileName_iscell];

iscell_A = readNPY(fullFileName_iscell_A);

s_A = load(fullFileName_Fall_A,'stat');
roi_stats_raw_A = s_A.stat;
temp_cellIndex_A = find(iscell_A(:,1)==1);
roi_stats_A = roi_stats_raw_A(temp_cellIndex_A); %#ok<*FNDSB>

roi_stats_simplified_A = cell(1,roiNum_A);
for tempi=1:roiNum_A
    % tempi = find(tempExactCellIndex_suite2p_currentFOV_A==40);
    temp_roi_stat = roi_stats_A{tempi};
    
    temp_roi_npix = double(temp_roi_stat.npix);
    temp_roi_xpix = double(temp_roi_stat.xpix);
    temp_roi_ypix = double(temp_roi_stat.ypix);
    temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
    %temp_roi_xy_ratio = double(temp_roi_stat.xy_ratio);
    temp_roi_xy_ratio = (max(temp_roi_xpix)-min(temp_roi_xpix))/(max(temp_roi_ypix)-min(temp_roi_ypix));
    
    temp_roi_stat_simplified = struct;
    temp_roi_stat_simplified.npix = temp_roi_npix;
    temp_roi_stat_simplified.xpix = temp_roi_xpix;
    temp_roi_stat_simplified.ypix = temp_roi_ypix;
    temp_roi_stat_simplified.med = temp_roi_med;
    temp_roi_stat_simplified.xy_ratio = temp_roi_xy_ratio;
    
    roi_stats_simplified_A{tempi} = temp_roi_stat_simplified;
end



% Load B anatomy info
path_plane_B = [input_path_B,'\plane0'];
fullFileName_Fall_B = [path_plane_B,'\',fileName_Fall];
fullFileName_iscell_B = [path_plane_B,'\',fileName_iscell];

iscell_B = readNPY(fullFileName_iscell_B);

s_B = load(fullFileName_Fall_B,'stat');
roi_stats_raw_B = s_B.stat;
temp_cellIndex_B = find(iscell_B(:,1)==1);
roi_stats_B = roi_stats_raw_B(temp_cellIndex_B);

roi_stats_simplified_B = cell(1,roiNum_B);
for tempi=1:roiNum_B
    temp_roi_stat = roi_stats_B{tempi};
    
    temp_roi_npix = double(temp_roi_stat.npix);
    temp_roi_xpix = double(temp_roi_stat.xpix);
    temp_roi_ypix = double(temp_roi_stat.ypix);
    temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
    %temp_roi_xy_ratio = double(temp_roi_stat.xy_ratio);
    temp_roi_xy_ratio = (max(temp_roi_xpix)-min(temp_roi_xpix))/(max(temp_roi_ypix)-min(temp_roi_ypix));
    
    
    temp_roi_stat_simplified = struct;
    temp_roi_stat_simplified.npix = temp_roi_npix;
    temp_roi_stat_simplified.xpix = temp_roi_xpix;
    temp_roi_stat_simplified.ypix = temp_roi_ypix;
    temp_roi_stat_simplified.med = temp_roi_med;
    temp_roi_stat_simplified.xy_ratio = temp_roi_xy_ratio;
    
    roi_stats_simplified_B{tempi} = temp_roi_stat_simplified;
end



roi_stats_simplified_A;
roi_stats_simplified_B;



%% Compute roi_stats_B_corrected
roi_stats_simplified_B_corrected = roi_stats_simplified_B;
for tempi=1:roiNum_B
    temp_roi_stat = roi_stats_simplified_B_corrected{tempi};
    temp_med = temp_roi_stat.med;
    
    temp_D = zeros(1,2);
    temp_D(2) = D_BtoA(temp_med(1),temp_med(2),1); % get x
    temp_D(1) = D_BtoA(temp_med(1),temp_med(2),2); % get y
    
    temp_med_corrected = temp_med - temp_D;
    
    temp_roi_stat.med_raw = temp_med;
    temp_roi_stat.med = temp_med_corrected;
    
    roi_stats_simplified_B_corrected{tempi} = temp_roi_stat;
end

a = 1;



%% Get med, npix, xy_ratio
% med;
% npix;
% xy_ratio;
med_A = zeros(2,roiNum_A);
npix_A = zeros(1,roiNum_A);
xy_ratio_A = zeros(1,roiNum_A);
for tempi=1:roiNum_A
    temp_roi_stat = roi_stats_simplified_A{tempi};
    med_A(:,tempi) = temp_roi_stat.med;
    npix_A(tempi) = temp_roi_stat.npix;
    xy_ratio_A(tempi) = temp_roi_stat.xy_ratio;
end
med_B_corrected = zeros(2,roiNum_B);
npix_B_corrected = zeros(1,roiNum_B);
xy_ratio_B_corrected = zeros(1,roiNum_B);
for tempi=1:roiNum_B
    temp_roi_stat = roi_stats_simplified_B_corrected{tempi};
    med_B_corrected(:,tempi) = temp_roi_stat.med;
    npix_B_corrected(tempi) = temp_roi_stat.npix;
    xy_ratio_B_corrected(tempi) = temp_roi_stat.xy_ratio;
end

%% medDisMatrix_AB
medDisMatrix_AB = zeros(roiNum_A,roiNum_B);
med_A;
med_B_corrected;
for tempi=1:roiNum_A
    for tempj=1:roiNum_B
        temp_error = med_A(:,tempi)-med_B_corrected(:,tempj);
        medDisMatrix_AB(tempi,tempj) = sqrt(sum(temp_error.^2));
    end
end


%% Calculate candidate_A
candidate_A = struct;
candidate_A.index = zeros(candidateNum,roiNum_A);
candidate_A.score = zeros(candidateNum,roiNum_A);
for tempi=1:roiNum_A
    tempA_med = med_A(:,tempi);
    tempA_npix = npix_A(tempi);
    tempA_xy_ratio = xy_ratio_A(tempi);
    
    temp_med_distance = sqrt(sum((tempA_med-med_B_corrected).^2,1));
    [~,I_med] = sort(temp_med_distance,'ascend');
    temp_candidateIndex = I_med(1:candidateNum);
    temp_med_distance_top = temp_med_distance(temp_candidateIndex) + 0.1;
        
    temptempBoolIndex = temp_med_distance_top > min(threshold_medDisMax,threshold_med_sqrtRatio*sqrt(tempA_npix));
    %temptempBoolIndex = temp_med_distance_top > threshold_medDisMax;   
    temp_candidateNum = sum(~temptempBoolIndex);
    
    temp_candidateIndex = I_med(1:temp_candidateNum);
    temp_med_distance_top = temp_med_distance(temp_candidateIndex) + 0.1;
    temp_candidateScore = (max(temp_med_distance_top)./temp_med_distance_top) * weight.med;
    
    temp_npix_distance = abs(tempA_npix-npix_B_corrected(temp_candidateIndex));
    [~,I_npix] = sort(temp_npix_distance,'ascend');
    temp_candidateScore(I_npix) = temp_candidateScore(I_npix) + (temp_candidateNum:-1:1) * weight.npix;
    
    temp_xy_ratio_distance = abs(tempA_xy_ratio-xy_ratio_B_corrected(temp_candidateIndex));
    [~,I_xy_ratio] = sort(temp_xy_ratio_distance,'ascend');
    temp_candidateScore(I_xy_ratio) = temp_candidateScore(I_xy_ratio) + (temp_candidateNum:-1:1) * weight.xy_ratio;
    
    a = 1;
    
    temp_covariance = covMatrix_AB(tempi,temp_candidateIndex);
    %temp_sse = meanSSEMatrix_AB(tempi,temp_candidateIndex);
    temp_sse = stdSSEMatrix_AB(tempi,temp_candidateIndex);        
    [temp_sort,I_sse] = sort(temp_sse,'ascend');    
    temp_addScore = (temp_candidateNum:-1:1) * weight.sseMatrix;
    temp_addScore(temp_covariance(I_sse)<0.1) = 0;%0.3-->0.2-->0.1-->0.2
    temp_candidateScore(I_sse) = temp_candidateScore(I_sse) + temp_addScore;
    
    temp_candidateScore(temp_covariance<0) = 0;
    temp_candidateScore(temp_covariance<eps) = temp_candidateScore(temp_covariance<eps)/2;
    if length(temp_candidateScore) >= 3
        %temp_candidateScore(temp_candidateScore < prctile(temp_candidateScore,50)) = 0;
        temp_candidateScore(temp_candidateScore < prctile(temp_candidateScore,50)) = 0.1;        
    end

    
    if length(temp_candidateScore) > 0 %#ok<*ISMT>
        for tempj=1:length(temp_candidateScore)
            if temp_covariance(tempj) < eps
                
                temptemp_med_threshold = min(7,0.8*sqrt(tempA_npix));%7               
                temp2 = temp_med_distance_top(tempj) > temptemp_med_threshold;                
                %temp3 = temp_npix_distance(tempj) > tempA_npix * 0.56;                            
                %temp3 = false;
                if tempA_npix > 30
                    temp3 = temp_npix_distance(tempj) > tempA_npix * 0.56;    
                    temp4 = temp_xy_ratio_distance(tempj) > tempA_xy_ratio * 0.35;
                else
                    temp3 = false;
                    temp4 = false;
                end
                temp234 = temp2 | temp3 | temp4;                
                if temp234 == 1
                    %temp_candidateScore(tempj) = 0;
                    temp_candidateScore(tempj) = 0.1;
                end
            end
        end
    end
        
    
    [~,I_scoreMax] = max(temp_candidateScore);
    
    %temp_bad_med_index = temp_med_distance_top > min(threshold_medDisMax,threshold_med_sqrtRatio*sqrt(tempA_npix));
    %if sum(temp_bad_med_index) > 0
    %    a = 1;
    %    temp_candidateScore(temp_bad_med_index) = 0;
    %end
        
    candidate_A.index(1:temp_candidateNum,tempi) = temp_candidateIndex;
    candidate_A.score(1:temp_candidateNum,tempi) = temp_candidateScore;
    
    %if tempi == 614 + 1
    if tempi == 421
       a = 1; %#ok<*NASGU>
    end
    
    temptempIndex_suite2p = tempExactCellIndex_suite2p_currentFOV_A(tempi);
    if temptempIndex_suite2p == 698 %775
        a = 1;
        temp_candidateIndex_suite2p = tempExactCellIndex_suite2p_currentFOV_B(temp_candidateIndex);
        a = 1;
    end
    
    
    a = 1;
end
validCandidateCount_A = sum(sum(candidate_A.score,1)>0);
fprintf('validCandidateCount_A = %d.\n',validCandidateCount_A);

a = 1;

%% Calculate candidate_B_corrected
candidate_B_corrected = struct;
candidate_B_corrected.index = zeros(candidateNum,roiNum_B);
candidate_B_corrected.score = zeros(candidateNum,roiNum_B);
for tempi=1:roiNum_B
    tempB_corrected_med = med_B_corrected(:,tempi);
    tempB_corrected_npix = npix_B_corrected(tempi);
    tempB_corrected_xy_ratio = xy_ratio_B_corrected(tempi);
    
    temp_med_distance = sqrt(sum((tempB_corrected_med-med_A).^2,1));
    [~,I_med] = sort(temp_med_distance,'ascend');
    temp_candidateIndex = I_med(1:candidateNum);
    temp_med_distance_top = temp_med_distance(temp_candidateIndex) + 0.1;
    
    temptempBoolIndex = temp_med_distance_top > min(threshold_medDisMax,threshold_med_sqrtRatio*sqrt(tempB_corrected_npix));
    %temptempBoolIndex = temp_med_distance_top > threshold_medDisMax;       
    temp_candidateNum = sum(~temptempBoolIndex);
    temp_candidateIndex = I_med(1:temp_candidateNum);

    temp_med_distance_top = temp_med_distance(temp_candidateIndex) + 0.1;    
    temp_candidateScore = (max(temp_med_distance_top)./temp_med_distance_top) * weight.med;
    
    temp_npix_distance = abs(tempB_corrected_npix-npix_A(temp_candidateIndex));
    [~,I_npix] = sort(temp_npix_distance,'ascend');
    temp_candidateScore(I_npix) = temp_candidateScore(I_npix) + (temp_candidateNum:-1:1) * weight.npix;
    
    temp_xy_ratio_distance = abs(tempB_corrected_xy_ratio-xy_ratio_A(temp_candidateIndex));
    [~,I_xy_ratio] = sort(temp_xy_ratio_distance,'ascend');
    temp_candidateScore(I_xy_ratio) = temp_candidateScore(I_xy_ratio) + (temp_candidateNum:-1:1) * weight.xy_ratio;
    
    
    temp_covariance = covMatrix_AB(temp_candidateIndex,tempi)';
    %temp_sse = meanSSEMatrix_AB(temp_candidateIndex,tempi)';
    temp_sse = stdSSEMatrix_AB(temp_candidateIndex,tempi)';        
    [temp_sort,I_sse] = sort(temp_sse,'ascend');    
    temp_addScore = (temp_candidateNum:-1:1) * weight.sseMatrix;
    temp_addScore(temp_covariance(I_sse)<0.1) = 0;%0.3-->0.2-->0.1-->0.2
    temp_candidateScore(I_sse) = temp_candidateScore(I_sse) + temp_addScore;
    
    temp_candidateScore(temp_covariance<0) = 0;
    temp_candidateScore(temp_covariance<eps) = temp_candidateScore(temp_covariance<eps)/2;
    if length(temp_candidateScore) >= 3
        %temp_candidateScore(temp_candidateScore < prctile(temp_candidateScore,50)) = 0;
        temp_candidateScore(temp_candidateScore < prctile(temp_candidateScore,50)) = 0.1;        
    end    
    
    if length(temp_candidateScore) > 0 %#ok<*ISMT>
        for tempj=1:length(temp_candidateScore)
            if temp_covariance(tempj) < eps

                temptemp_med_threshold = min(7,0.8*sqrt(tempB_corrected_npix));%7               
                temp2 = temp_med_distance_top(tempj) > temptemp_med_threshold;
                %temp3 = temp_npix_distance(tempj) > tempB_corrected_npix * 0.56;            
                %temp3 = false;                
                if tempB_corrected_npix > 30
                    temp3 = temp_npix_distance(tempj) > tempB_corrected_npix * 0.56;
                    temp4 = temp_xy_ratio_distance(tempj) > tempB_corrected_xy_ratio * 0.35;
                else
                    temp3 = false;
                    temp4 = false;
                end                
                temp234 = temp2 | temp3 | temp4;                
                if temp234 == 1
                    %temp_candidateScore(tempj) = 0;
                    temp_candidateScore(tempj) = 0.1;
                end
            end
        end
    end 
    
    
    [~,I_scoreMax] = max(temp_candidateScore);
    %if temp_med_distance_top(I_scoreMax) > min(threshold_medDisMax,threshold_med_sqrtRatio*sqrt(tempB_corrected_npix))
    %    continue
    %end
    
    %temp_bad_med_index = temp_med_distance_top > min(threshold_medDisMax,threshold_med_sqrtRatio*sqrt(tempB_corrected_npix));
    %if sum(temp_bad_med_index) > 0
    %    a = 1;
    %    temp_candidateScore(temp_bad_med_index) = 0;
    %end
    
    candidate_B_corrected.index(1:temp_candidateNum,tempi) = temp_candidateIndex;
    candidate_B_corrected.score(1:temp_candidateNum,tempi) = temp_candidateScore;
    
    if tempi == 195
        a = 1; %#ok<*NASGU>
    end
    
    temptempIndex_suite2p = tempExactCellIndex_suite2p_currentFOV_B(tempi);
    if temptempIndex_suite2p == 1175
        a = 1;
        temp_candidateIndex_suite2p = tempExactCellIndex_suite2p_currentFOV_A(temp_candidateIndex);
        a = 1;
    end    
    
    a = 1;
end
validCandidateCount_B_corrected = sum(sum(candidate_B_corrected.score,1)>0);
fprintf('validCandidateCount_B_corrected = %d.\n',validCandidateCount_B_corrected);


%% Calculate valid_pairIndex_AtoB
confusionMatrix_AtoB = zeros(roiNum_A,roiNum_B);
for tempi=1:roiNum_A
    temp_A_index = candidate_A.index(:,tempi);
    if sum(temp_A_index) == 0
        continue
    end
    temp_A_score = candidate_A.score(:,tempi);
    
    temp_candidateNum = sum(temp_A_index>0);
    temp_A_index_valid = temp_A_index(1:temp_candidateNum);
    temp_A_score_valid = temp_A_score(1:temp_candidateNum);
    
    [~,I_score] = sort(temp_A_score_valid,'descend');    
    
    %confusionMatrix_AtoB(tempi,temp_A_index_valid(I_score)) = (temp_candidateNum:-1:1);
    %confusionMatrix_AtoB(tempi,temp_A_index_valid(temp_A_score_valid<eps)) = 0;
    
    confusionMatrix_AtoB(tempi,temp_A_index_valid) = temp_A_score_valid;
        
end
confusionMatrix_BtoA = zeros(roiNum_B,roiNum_A);
for tempi=1:roiNum_B
    temp_B_index = candidate_B_corrected.index(:,tempi);
    if sum(temp_B_index) == 0
        continue
    end
    temp_B_score = candidate_B_corrected.score(:,tempi);
    
    temp_candidateNum = sum(temp_B_index>0);
    temp_B_index_valid = temp_B_index(1:temp_candidateNum);
    temp_B_score_valid = temp_B_score(1:temp_candidateNum);
    
    [~,I_score] = sort(temp_B_score_valid,'descend');    
    
    %confusionMatrix_BtoA(tempi,temp_B_index_valid(I_score)) = (temp_candidateNum:-1:1);
    %confusionMatrix_BtoA(tempi,temp_B_index_valid(temp_B_score_valid<eps)) = 0;
    
    confusionMatrix_BtoA(tempi,temp_B_index_valid) = temp_B_score_valid;
    
end
confusionMatrix_BtoA_trans = confusionMatrix_BtoA';

% boolMatrix_ABpair = (confusionMatrix_AtoB > 0) & (confusionMatrix_BtoA_trans>0);

temp_confusionSum = confusionMatrix_AtoB + confusionMatrix_BtoA_trans;
boolMatrix_ABpair = (temp_confusionSum > 0.2) & (confusionMatrix_AtoB > 0) & (confusionMatrix_BtoA_trans > 0);

a1 = sum(boolMatrix_ABpair,2);
a2 = sum(boolMatrix_ABpair,1);

confusionMatrix_AtoB_shrink = confusionMatrix_AtoB .* boolMatrix_ABpair;
confusionMatrix_BtoA_trans_shrink = confusionMatrix_BtoA_trans .* boolMatrix_ABpair;


confusionMatrix_AtoB_shrink_raw = confusionMatrix_AtoB_shrink;
confusionMatrix_BtoA_trans_shrink_raw = confusionMatrix_BtoA_trans_shrink;

%% valid_pairIndex_AandB
if if_AB_uniquePair == 1
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
    valid_pairIndex_AtoB = -1*ones(1,roiNum_A);
    for tempi=1:roiNum_A
        temp_findIndex = find(confusionMatrix_AtoB_shrink_maxFilter(tempi,:)>0);
        if isempty(temp_findIndex) == 0
            valid_pairIndex_AtoB(tempi) = temp_findIndex;
        end
    end
    valid_pairIndex_AandB = zeros(valid_pairCount,2);
    valid_pairIndex_AandB(:,1) = find(valid_pairIndex_AtoB>-1);
    valid_pairIndex_AandB(:,2) = valid_pairIndex_AtoB(valid_pairIndex_AandB(:,1));
elseif if_AB_uniquePair == 0
    confusionMatrix_AB_shrink_sum = confusionMatrix_AtoB_shrink + confusionMatrix_BtoA_trans_shrink;
    
    %% valid_pairIndex_AandB_AtoB
    %temp_max = max(confusionMatrix_AB_shrink_sum,[],2);
    %confusionMatrix_AtoB_shrink_maxFilter = confusionMatrix_AB_shrink_sum;
    %confusionMatrix_AtoB_shrink_maxFilter(confusionMatrix_AB_shrink_sum < temp_max) = 0;
    
    [temp_max,I] = max(confusionMatrix_AB_shrink_sum,[],2);
    confusionMatrix_AtoB_shrink_maxFilter = zeros(size(confusionMatrix_AB_shrink_sum));
    for tempi=1:roiNum_A
        confusionMatrix_AtoB_shrink_maxFilter(tempi,I(tempi)) = temp_max(tempi);
    end    
    
    valid_pairCount_AtoB = sum(confusionMatrix_AtoB_shrink_maxFilter>0,'all');
    valid_pairIndex_AtoB = -1*ones(1,roiNum_A);
    for tempi=1:roiNum_A
        temp_findIndex = find(confusionMatrix_AtoB_shrink_maxFilter(tempi,:)>0);
        if isempty(temp_findIndex) == 0
            valid_pairIndex_AtoB(tempi) = temp_findIndex;
        end
    end
    valid_pairIndex_AandB_AtoB = zeros(valid_pairCount_AtoB,2);
    valid_pairIndex_AandB_AtoB(:,1) = find(valid_pairIndex_AtoB>-1);
    valid_pairIndex_AandB_AtoB(:,2) = valid_pairIndex_AtoB(valid_pairIndex_AandB_AtoB(:,1));
    
    
    %% valid_pairIndex_AandB_BtoA
    %temp_max = max(confusionMatrix_AB_shrink_sum,[],1);
    %confusionMatrix_BtoA_trans_shrink_maxFilter = confusionMatrix_AB_shrink_sum;
    %confusionMatrix_BtoA_trans_shrink_maxFilter(confusionMatrix_AB_shrink_sum < temp_max) = 0;
    
    [temp_max,I] = max(confusionMatrix_AB_shrink_sum,[],1);
    confusionMatrix_BtoA_trans_shrink_maxFilter = zeros(size(confusionMatrix_AB_shrink_sum));
    for tempi=1:roiNum_B
        confusionMatrix_BtoA_trans_shrink_maxFilter(I(tempi),tempi) = temp_max(tempi);
    end
    
    valid_pairCount_BtoA = sum(confusionMatrix_BtoA_trans_shrink_maxFilter>0,'all');
    valid_pairIndex_BtoA = -1*ones(1,roiNum_B);
    for tempi=1:roiNum_B
        temp_findIndex = find(confusionMatrix_BtoA_trans_shrink_maxFilter(:,tempi)>0);
        if isempty(temp_findIndex) == 0
            valid_pairIndex_BtoA(tempi) = temp_findIndex;
        end
    end
    valid_pairIndex_AandB_BtoA = zeros(valid_pairCount_BtoA,2);
    valid_pairIndex_AandB_BtoA(:,2) = find(valid_pairIndex_BtoA>-1);
    valid_pairIndex_AandB_BtoA(:,1) = valid_pairIndex_BtoA(valid_pairIndex_AandB_BtoA(:,2));
    
    
    %% valid_pairIndex_AandB
    a1 = unique(valid_pairIndex_AandB_AtoB(:,1));
    a2 = unique(valid_pairIndex_AandB_AtoB(:,2));
    a3 = unique(valid_pairIndex_AandB_BtoA(:,1));
    a4 = unique(valid_pairIndex_AandB_BtoA(:,2));
    fprintf('num_A_AtoB=%d,num_B_AtoB=%d,num_A_BtoA=%d,num_B_BtoA=%d.\n',length(a1),length(a2),length(a3),length(a4));
    
    valid_pairIndex_AandB_raw = [valid_pairIndex_AandB_AtoB;valid_pairIndex_AandB_BtoA];
    [M,I] = sort(valid_pairIndex_AandB_raw(:,1));
    valid_pairIndex_AandB_raw_Asorted = valid_pairIndex_AandB_raw(I,:);
    
    valid_pairIndex_AandB_raw_Asorted_unique = unique(valid_pairIndex_AandB_raw_Asorted,'rows');
    aa1 = unique(valid_pairIndex_AandB_raw_Asorted_unique(:,1));
    aa2 = unique(valid_pairIndex_AandB_raw_Asorted_unique(:,2));
    
    [M,I] = sort(valid_pairIndex_AandB_raw_Asorted_unique(:,2));
    valid_pairIndex_AandB_raw_Bsorted_unique = valid_pairIndex_AandB_raw_Asorted_unique(I,:);
    bb1 = unique(valid_pairIndex_AandB_raw_Bsorted_unique(:,1));
    bb2 = unique(valid_pairIndex_AandB_raw_Bsorted_unique(:,2));
    
    valid_pairIndex_AandB_raw_Asorted_unique;
    valid_pairIndex_AandB_raw_Bsorted_unique;
    
    valid_pairCount = size(valid_pairIndex_AandB_raw_Asorted_unique,1);

    
    confusionMatrix_AB_shrink_sum;
    valid_pairIndex_AandB_raw_Asorted_unique_withScore = zeros(valid_pairCount,4);
    valid_pairIndex_AandB_raw_Asorted_unique_withScore(:,1:2) = valid_pairIndex_AandB_raw_Asorted_unique;
    for tempi=1:valid_pairCount        
        temptempx = valid_pairIndex_AandB_raw_Asorted_unique(tempi,1);
        temptempy = valid_pairIndex_AandB_raw_Asorted_unique(tempi,2);        
        temp_score = confusionMatrix_AB_shrink_sum(temptempx,temptempy);
        temp_dis = medDisMatrix_AB(temptempx,temptempy);        
        valid_pairIndex_AandB_raw_Asorted_unique_withScore(tempi,3) = temp_score;
        valid_pairIndex_AandB_raw_Asorted_unique_withScore(tempi,4) = temp_dis;        
    end
    
    [M,I] = sort(valid_pairIndex_AandB_raw_Asorted_unique_withScore(:,2));
    valid_pairIndex_AandB_raw_Bsorted_unique_withScore = valid_pairIndex_AandB_raw_Asorted_unique_withScore(I,:);
    
    valid_pairIndex_AandB_raw_Asorted_unique_withScore;
    
    temp_array = valid_pairIndex_AandB_raw_Asorted_unique_withScore;
    %for tempi
    %[temp_min,temp_max] = bounds(temp_array(:,2));
    
    tempInValidBoolIndex = false(valid_pairCount,1);
    
    [temp_min,temp_max] = bounds(temp_array(:,1));    
    for tempi=temp_min:temp_max
        temptempBoolIndex1 = temp_array(:,1) == tempi;
        temptempIndex1 = find(temp_array(:,1) == tempi);        
        if sum(temptempBoolIndex1) <= 1
            continue
        end
        %[M,I] = max(temp_array(temptempBoolIndex1,3)); 
        %[M,I] = min(temp_array(temptempBoolIndex1,4));
        
        if sum(temptempBoolIndex1) <= 3
            [M,I] = min(temp_array(temptempBoolIndex1,4));            
        elseif sum(temptempBoolIndex1) >= 4
            temptempIndex1Z = temptempIndex1(temp_array(temptempBoolIndex1,4) < median(temp_array(temptempBoolIndex1,4)));
            if length(temptempIndex1Z) > 0
                [M,IZ] = max(temp_array(temptempIndex1Z,3));
                I = find(temptempIndex1 == temptempIndex1Z(IZ));
            elseif length(temptempIndex1Z) == 0
                [M,I] = min(temp_array(temptempBoolIndex1,4));
            end
        end
        
        temptempIndex2 = temptempIndex1;
        temptempIndex2(I) = [];
        
        tempInValidBoolIndex(temptempIndex2) = true;        
    end
    
    
    [temp_min,temp_max] = bounds(temp_array(:,2));    
    for tempi=temp_min:temp_max
        temptempBoolIndex1 = temp_array(:,2) == tempi;
        temptempIndex1 = find(temp_array(:,2) == tempi);        
        if sum(temptempBoolIndex1) <= 1
            continue
        end
        %[M,I] = max(temp_array(temptempBoolIndex1,3)); 
        %[M,I] = min(temp_array(temptempBoolIndex1,4));  
        
        if sum(temptempBoolIndex1) <= 3
            [M,I] = min(temp_array(temptempBoolIndex1,4));             %#ok<*ASGLU>
        elseif sum(temptempBoolIndex1) >= 4
            temptempIndex1Z = temptempIndex1(temp_array(temptempBoolIndex1,4) < median(temp_array(temptempBoolIndex1,4)));
            if length(temptempIndex1Z) > 0
                [M,IZ] = max(temp_array(temptempIndex1Z,3));
                I = find(temptempIndex1 == temptempIndex1Z(IZ));
            elseif length(temptempIndex1Z) == 0
                [M,I] = min(temp_array(temptempBoolIndex1,4));
            end
        end
        
        temptempIndex2 = temptempIndex1;
        temptempIndex2(I) = [];
        
        tempInValidBoolIndex(temptempIndex2) = true;        
    end
    
    valid_pairIndex_AandB_raw_Asorted_unique_withScore_valid = valid_pairIndex_AandB_raw_Asorted_unique_withScore(~tempInValidBoolIndex,:);
    
    [M,I] = sort(valid_pairIndex_AandB_raw_Asorted_unique_withScore_valid(:,2));
    valid_pairIndex_AandB_raw_Bsorted_unique_withScore_valid = valid_pairIndex_AandB_raw_Asorted_unique_withScore_valid(I,:);
    
    aaa1 = unique(valid_pairIndex_AandB_raw_Asorted_unique_withScore_valid(:,1));
    aaa2 = unique(valid_pairIndex_AandB_raw_Asorted_unique_withScore_valid(:,2));
    
    
    %valid_pairIndex_AandB = valid_pairIndex_AandB_raw_Asorted_unique;    
    valid_pairIndex_AandB = valid_pairIndex_AandB_raw_Asorted_unique_withScore_valid(:,1:2);
    valid_pairCount = size(valid_pairIndex_AandB,1);
    
    
    
    valid_pairIndex_AandB_Asorted_withScore_suite2p = valid_pairIndex_AandB_raw_Asorted_unique_withScore;
    valid_pairIndex_AandB_Asorted_withScore_suite2p(:,1) = ...
        tempExactCellIndex_suite2p_currentFOV_A(valid_pairIndex_AandB_raw_Asorted_unique_withScore(:,1));
    valid_pairIndex_AandB_Asorted_withScore_suite2p(:,2) = ...
        tempExactCellIndex_suite2p_currentFOV_B(valid_pairIndex_AandB_raw_Asorted_unique_withScore(:,2));

    valid_pairIndex_AandB_Bsorted_withScore_suite2p = valid_pairIndex_AandB_raw_Bsorted_unique_withScore;
    valid_pairIndex_AandB_Bsorted_withScore_suite2p(:,1) = ...
        tempExactCellIndex_suite2p_currentFOV_A(valid_pairIndex_AandB_raw_Bsorted_unique_withScore(:,1));
    valid_pairIndex_AandB_Bsorted_withScore_suite2p(:,2) = ...
        tempExactCellIndex_suite2p_currentFOV_B(valid_pairIndex_AandB_raw_Bsorted_unique_withScore(:,2));
    
    
    valid_pairIndex_AandB_Asorted_withScore_valid_suite2p = valid_pairIndex_AandB_raw_Asorted_unique_withScore_valid;
    valid_pairIndex_AandB_Asorted_withScore_valid_suite2p(:,1) = ...
        tempExactCellIndex_suite2p_currentFOV_A(valid_pairIndex_AandB_raw_Asorted_unique_withScore_valid(:,1));
    valid_pairIndex_AandB_Asorted_withScore_valid_suite2p(:,2) = ...
        tempExactCellIndex_suite2p_currentFOV_B(valid_pairIndex_AandB_raw_Asorted_unique_withScore_valid(:,2));


end



fprintf('valid_pairCount = %d.\n',valid_pairCount);
valid_pairIndex_AandB;% the index is cellIndex, not suite2p index


%% Plot
if if_plot == 1
    close all
    
    %% Plot FOV A
    fig = figure('Name','template_A','NumberTitle','off'); %#ok<*NASGU>
    set(gcf,'Position',[400 50 800 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    template_plot = template_A;
    
    temp_template_plot = template_plot / max(template_plot,[],'all');
    
    temp_template2_plot = (temp_template_plot.^0.55)*255;
    
    image(temp_template2_plot);
    hold on
    colormap(gray);
    
    for tempi=1:plotRoiNum
        temp_roi_stat = roi_stats_simplified_A{tempi};
        temp_med = temp_roi_stat.med;
        scatter(temp_med(2),temp_med(1),35,[0.8500 0.3250 0.0980],'filled');
        hold on
    end
    
    axis equal;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    set(gca,'Visible','off');
    
    
    %% Plot FOV B
    fig = figure('Name','template_B','NumberTitle','off');
    set(gcf,'Position',[400 50 800 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    template_plot = template_B;
    
    temp_template_plot = template_plot / max(template_plot,[],'all');
    
    temp_template2_plot = (temp_template_plot.^0.55)*255;
    
    image(temp_template2_plot);
    hold on
    colormap(gray);
    
    for tempi=1:plotRoiNum
        temp_roi_stat = roi_stats_simplified_B{tempi};
        temp_med = temp_roi_stat.med;
        scatter(temp_med(2),temp_med(1),35,[0.8500 0.3250 0.0980],'filled');
        hold on
    end
    
    
    axis equal;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    set(gca,'Visible','off');
    
    %% Plot FOV B_corrected
    fig = figure('Name','template_B_corrected','NumberTitle','off');
    set(gcf,'Position',[400 50 800 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    template_plot = template_B_corrected;
    
    temp_template_plot = template_plot / max(template_plot,[],'all');
    
    temp_template2_plot = (temp_template_plot.^0.55)*255;
    
    image(temp_template2_plot);
    hold on
    colormap(gray);
    
    for tempi=1:plotRoiNum
        temp_roi_stat = roi_stats_simplified_B_corrected{tempi};
        temp_med = temp_roi_stat.med;
        scatter(temp_med(2),temp_med(1),35,[0.8500 0.3250 0.0980],'filled');
        hold on
    end
    
    axis equal;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    set(gca,'Visible','off');
end