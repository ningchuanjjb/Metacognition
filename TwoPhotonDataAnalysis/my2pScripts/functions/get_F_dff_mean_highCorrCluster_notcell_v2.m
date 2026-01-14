function [F_dff_median_highCorrCluster_notcell,F_dff_mean_highCorrCluster_notcell,highCorrClusterIndex_notcell_suite2p] = ...
    get_F_dff_mean_highCorrCluster_notcell_v2(F_dff_notcell,temp_binSize,corr_threshold,iscell)

% F_dff_notcell;

% temp_binSize = 4;
% corr_threshold = 0.5;

fprintf('Now compute highCorrCluster...');


F_dff_notcell_2 = F_dff_notcell(:,1:floor(size(F_dff_notcell,2)/temp_binSize)*temp_binSize);
F_dff_notcell_2_3d = reshape(F_dff_notcell_2,[size(F_dff_notcell_2,1) temp_binSize size(F_dff_notcell_2,2)/temp_binSize]);
F_dff_notcell_2_bin = squeeze(mean(F_dff_notcell_2_3d,2));

F_dff_notcell_2_bin_T = F_dff_notcell_2_bin';

clear F_dff_notcell_2 F_dff_notcell_2_3d F_dff_notcell_2_bin

corr_F_dff_notcell = zeros(size(F_dff_notcell_2_bin_T,2),size(F_dff_notcell_2_bin_T,2));
t0 = tic;
parfor tempi=1:size(F_dff_notcell_2_bin_T,2)
    corr_F_dff_notcell(tempi,:) = corr(F_dff_notcell_2_bin_T(:,tempi),F_dff_notcell_2_bin_T); %#ok<*PFBNS>
    a = 1;
end
% toc(t0);
fprintf('Get highCorrCluster! Time is %.1f secs.\n',toc(t0));
a = 1;

% corr_F_dff_notcell;



corrBool_F_dff_notcell = corr_F_dff_notcell >= corr_threshold;

corrBoolSum_F_dff_notcell = sum(corrBool_F_dff_notcell,2);

[M,I] = max(corrBoolSum_F_dff_notcell);

highCorrClusterIndex_notcell = find(corrBool_F_dff_notcell(I,:)==true)';

F_dff_highCorrCluster_notcell = F_dff_notcell(highCorrClusterIndex_notcell,:);
F_dff_mean_highCorrCluster_notcell = mean(F_dff_highCorrCluster_notcell,1);
F_dff_median_highCorrCluster_notcell = median(F_dff_highCorrCluster_notcell,1);

index_notcell = find(iscell(:,1)==0);
index_notcell_suite2p = index_notcell-1;
highCorrClusterIndex_notcell_suite2p = index_notcell_suite2p(highCorrClusterIndex_notcell);

a = 1;

