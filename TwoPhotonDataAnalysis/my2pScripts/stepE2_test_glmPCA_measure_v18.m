close all

glm_beta_flag = 3;


% 1 B_P;
% 2 B_oneMinusP;
% 3 B_P0;
% 4 B_oneMinusP0;

x1_flag = 1;
x1_proj_flag = 2;
x2_flag = 1;
x2_proj_flag = 3;

% if_plot_x_B0_Fmean1 = 1;

if_vectorDis0_vectorProjRatio1 = 1;

norm_coeff = 2;
dimNum = 3;
if_usePCA = 1;

if_loadCoeff_model0_behavior1 = 1;


if glm_beta_flag == 1
    glm_beta = glm_beta_lengthx_baselineBin_multiFOV;
    glm_beta_shuffled = glm_beta_lengthx_baselineBin_multiFOV_shuffled;
elseif glm_beta_flag == 2
    glm_beta = glm_beta_lengthx_delay1Bin_multiFOV;
    glm_beta_shuffled = glm_beta_lengthx_delay1Bin_multiFOV_shuffled;
elseif glm_beta_flag == 3
    glm_beta = glm_beta_lengthx_delay2Bin_multiFOV;
    glm_beta_shuffled = glm_beta_lengthx_delay2Bin_multiFOV_shuffled;
end


if if_plot_proj_B0_Fmean1_trial2 == 1
    F_designMatrixMean_baselineBin_multiFOV;
    F_designMatrixMean_delay1Bin_multiFOV;
    F_designMatrixMean_delay2Bin_multiFOV;
    
    a = 1;
    
    if glm_beta_flag == 1
        F_mean = F_designMatrixMean_baselineBin_multiFOV;
    elseif glm_beta_flag == 2
        F_mean = F_designMatrixMean_delay1Bin_multiFOV;
    elseif glm_beta_flag == 3
        F_mean = F_designMatrixMean_delay2Bin_multiFOV;
    end
    a = 1;
    
    F_P = F_mean(:,(1:numFrames)+numFrames*0);
    F_oneMinusP = F_mean(:,(1:numFrames)+numFrames*1);
    F_P0 = F_mean(:,(1:numFrames)+numFrames*2);
    F_oneMinusP0 = F_mean(:,(1:numFrames)+numFrames*3);
end


if exist('glm_r2_lengthx_delay1Bin_multiFOV','var') == 1
    if_plot_histogram = 0;
    if glm_beta_flag == 1
        glm_r2 = glm_r2_lengthx_baselineBin_multiFOV;
    elseif glm_beta_flag == 2
        glm_r2 = glm_r2_lengthx_delay1Bin_multiFOV;
    elseif glm_beta_flag == 3
        glm_r2 = glm_r2_lengthx_delay2Bin_multiFOV;
    end
    glm_r2_hit = glm_r2(:,1);
    glm_r2_miss = glm_r2(:,2);
    glm_r2_falseAlarm = glm_r2(:,3);
    glm_r2_correctRejection = glm_r2(:,4);
    
    glm_r2_hit_prctile = sum(glm_r2_hit<0)/length(glm_r2_hit);
    glm_r2_miss_prctile = sum(glm_r2_miss<0)/length(glm_r2_miss);
    glm_r2_falseAlarm_prctile = sum(glm_r2_falseAlarm<0)/length(glm_r2_falseAlarm);
    glm_r2_correctRejection_prctile = sum(glm_r2_correctRejection<0)/length(glm_r2_correctRejection);
    
    if if_plot_histogram == 1
        fig1 = figure('Name','Fig1','NumberTitle','off');
        set(gcf,'Position',[35+000 35+000 1430 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
        
        nexttile
        h1 = histogram(glm_r2_hit);
        xlabel(sprintf('r2(hit)'),'fontsize',16);
        title(sprintf('%%%.1f less than 0',glm_r2_hit_prctile*100),'fontsize',16);
        nexttile
        h2 = histogram(glm_r2_miss);
        xlabel(sprintf('r2(miss)'),'fontsize',16);
        title(sprintf('%%%.1f less than 0',glm_r2_miss_prctile*100),'fontsize',16);
        nexttile
        h3 = histogram(glm_r2_falseAlarm);
        xlabel(sprintf('r2(falseAlarm)'),'fontsize',16);
        title(sprintf('%%%.1f less than 0',glm_r2_falseAlarm_prctile*100),'fontsize',16);
        nexttile
        h4 = histogram(glm_r2_correctRejection);
        xlabel(sprintf('r2(correctRejection)'),'fontsize',16);
        title(sprintf('%%%.1f less than 0',glm_r2_correctRejection_prctile*100),'fontsize',16);
    end
end
a = 1;

glm_beta_raw = glm_beta;
glm_beta = glm_beta_raw(r2Valid_boolIndex_inOne,:);
glm_beta_shuffled_raw = glm_beta_shuffled;
glm_beta_shuffled = glm_beta_shuffled_raw(r2Valid_boolIndex_inOne,:);



% [M1,I1] = sort(glm_r2_lengthx_delay2Bin_multiFOV(:,1),'descend');
% [M2,I2] = sort(glm_r2_lengthx_delay2Bin_multiFOV(:,2),'descend');
% [M3,I3] = sort(glm_r2_lengthx_delay2Bin_multiFOV(:,3),'descend');
% [M4,I4] = sort(glm_r2_lengthx_delay2Bin_multiFOV(:,4),'descend');
%
% % temptmep_beta = glm_beta(I1(1:30),:);
% temptmep_beta = glm_beta(I2(1:30),:);

B_P = glm_beta(:,(1:numFrames)+numFrames*0);
B_oneMinusP = glm_beta(:,(1:numFrames)+numFrames*1);
B_P0 = glm_beta(:,(1:numFrames)+numFrames*2);
B_oneMinusP0 = glm_beta(:,(1:numFrames)+numFrames*3);

B_P_shuffled = glm_beta_shuffled(:,(1:numFrames)+numFrames*0);
B_oneMinusP_shuffled = glm_beta_shuffled(:,(1:numFrames)+numFrames*1);
B_P0_shuffled = glm_beta_shuffled(:,(1:numFrames)+numFrames*2);
B_oneMinusP0_shuffled = glm_beta_shuffled(:,(1:numFrames)+numFrames*3);

B_multi = cell(1,4);
B_shuffled_multi = cell(1,4);
F_multi = cell(1,4);

B_multi{1} = B_P;
B_multi{2} = B_oneMinusP;
B_multi{3} = B_P0;
B_multi{4} = B_oneMinusP0;

B_shuffled_multi{1} = B_P_shuffled;
B_shuffled_multi{2} = B_oneMinusP_shuffled;
B_shuffled_multi{3} = B_P0_shuffled;
B_shuffled_multi{4} = B_oneMinusP0_shuffled;

F_multi{1} = F_P;
F_multi{2} = F_oneMinusP;
F_multi{3} = F_P0;
F_multi{4} = F_oneMinusP0;

x1 = B_multi{x1_flag}';
x1_shuffled = B_shuffled_multi{x1_flag}';
x1_proj = B_multi{x1_proj_flag}';
x2 = B_multi{x2_flag}';
x2_shuffled = B_shuffled_multi{x2_flag}';
x2_proj = B_multi{x2_proj_flag}';

% if if_plot_x_B0_Fmean1 == 1 && if_plot_proj_B0_Fmean1_trial2 == 1
if if_plot_proj_B0_Fmean1_trial2 == 1
    %x1 = F_multi{x1_flag}';
    x1_proj = F_multi{x1_proj_flag}';
    %x2 = F_multi{x2_flag}';
    x2_proj = F_multi{x2_proj_flag}';
end

a = 1;


%% x1
x = x1;
x_shuffled = x1_shuffled;
x_proj = x1_proj;
[temp_coeff,temp_score,temp_latent,temp_tsquared,temp_explained,temp_mu] = pca(x); %#ok<*ASGLU>

% temp_score_shuffled = x_shuffled(1,:) * temp_coeff;
temp_score_shuffled = (x_shuffled(1,:)-temp_mu) * temp_coeff;
% temp_score_proj = x_proj * temp_coeff;
temp_score_proj = (x_proj-mean(x_proj,1)) * temp_coeff;


if if_plot_x_B0_Fmean1 == 1 && if_plot_proj_B0_Fmean1_trial2 == 1
    x_new = F_multi{x1_flag}';
    temp_score_new = (x_new-mean(x_new,1)) * temp_coeff;
    temp_score = temp_score_new;
end

dis_B_x1 = zeros(numFrames,1);
dis_B_x1_proj = zeros(numFrames,1);
for tempi=1:numFrames
    if if_vectorDis0_vectorProjRatio1 == 0
        if if_usePCA == 1
            %temptemp_score = temp_score(tempi,1:dimNum);
            %temptemp_score = temp_score(tempi,1:dimNum) - temp_score_shuffled(1:dimNum);
            temptemp_score = temp_score(tempi,1:dimNum) - temp_score_proj(tempi,1:dimNum);
            
            %temptemp_score_proj = temp_score_proj(tempi,1:dimNum);
            %temptemp_score_proj = temp_score_proj(tempi,1:dimNum) - temp_score_shuffled(1:dimNum);
            %temptemp_score_proj = 0 - temp_score_shuffled(1:dimNum);
        elseif if_usePCA == 0
            %temptemp_score = x1(tempi,:);
            %temptemp_score = x(tempi,:) - x_shuffled(tempi,:);
            %temptemp_score = x(tempi,:) - x2(tempi,:);
            
            %temptemp_score_proj = x2(tempi,:) - x_shuffled(tempi,:);
        end
        dis_B_x1(tempi) = norm(temptemp_score,norm_coeff);
        %dis_B_x1_proj(tempi) = norm(temptemp_score_proj,norm_coeff);
        
    elseif if_vectorDis0_vectorProjRatio1 == 1
        if if_usePCA == 1
            a = 1;
            v1 = temp_score(tempi,1:dimNum); % vector 1
            v2 = temp_score_proj(tempi,1:dimNum); % vector 2
            
            theta = (acos(dot(v1,v2) / (norm(v1)*norm(v2))) / pi)*180;
            v2_to_v1_projRatio = norm(v2)*cos(pi*theta/180)/norm(v1);
            
            dis_B_x1(tempi) = v2_to_v1_projRatio;
            %dis_B_x1(tempi) = theta;
        elseif if_usePCA == 0
            
            
        end
        
    end
end
dis_B_x1_raw = dis_B_x1;
% dis_B_x1 = dis_B_x1 ./ dis_B_x1_proj;
% dis_B_x1 = dis_B_x1 - dis_B_x1_proj;

a = 1;


%% x2
x = x2;
x_shuffled = x2_shuffled;
x_proj = x2_proj;
[temp_coeff,temp_score,temp_latent,temp_tsquared,temp_explained,temp_mu] = pca(x);

% temp_score_shuffled = x_shuffled(1,:) * temp_coeff;
temp_score_shuffled = (x_shuffled(1,:)-temp_mu) * temp_coeff;
% temp_score_proj = x_proj * temp_coeff;
temp_score_proj = (x_proj-mean(x_proj,1)) * temp_coeff;

if if_plot_x_B0_Fmean1 == 1 && if_plot_proj_B0_Fmean1_trial2 == 1
    x_new = F_multi{x2_flag}';
    temp_score_new = (x_new-mean(x_new,1)) * temp_coeff;
    temp_score = temp_score_new;
end

dis_B_x2 = zeros(numFrames,1);
dis_B_x2_proj = zeros(numFrames,1);
for tempi=1:numFrames
    if if_vectorDis0_vectorProjRatio1 == 0
        if if_usePCA == 1
            %temptemp_score = temp_score(tempi,1:dimNum);
            %temptemp_score = temp_score(tempi,1:dimNum) - temp_score_shuffled(1:dimNum);
            temptemp_score = temp_score(tempi,1:dimNum) - temp_score_proj(tempi,1:dimNum);
            
            %temptemp_score_proj = temp_score_proj(tempi,1:dimNum);
            %temptemp_score_proj = temp_score_proj(tempi,1:dimNum) - temp_score_shuffled(1:dimNum);
        elseif if_usePCA == 0
            %temptemp_score = x2(tempi,:);
            %temptemp_score = x(tempi,:) - x_shuffled(tempi,:);
            
            %temptemp_score_proj = x1(tempi,:) - x_shuffled(tempi,:);
        end
        dis_B_x2(tempi) = norm(temptemp_score,norm_coeff);
        %dis_B_x2_proj(tempi) = norm(temptemp_score_proj,norm_coeff);
        
    elseif if_vectorDis0_vectorProjRatio1 == 1
        if if_usePCA == 1
            a = 1;
            v1 = temp_score(tempi,1:dimNum); % vector 1
            v2 = temp_score_proj(tempi,1:dimNum); % vector 2
            
            theta = (acos(dot(v1,v2) / (norm(v1)*norm(v2))) / pi)*180;
            v2_to_v1_projRatio = norm(v2)*cos(pi*theta/180)/norm(v1);
            
            dis_B_x2(tempi) = v2_to_v1_projRatio;
            %dis_B_x2(tempi) = theta;            
        elseif if_usePCA == 0
            
            
        end        
        
    end
end
dis_B_x2_raw = dis_B_x2;
% dis_B_x2 = dis_B_x2 ./ dis_B_x2_proj;
% dis_B_x2 = dis_B_x2 - dis_B_x2_proj;
% dis_B_x2 = dis_B_x2_proj;

%%
a = 1;

dis_B_x1 = dis_B_x1';
dis_B_x2 = dis_B_x2';

a = 1;
if if_loadCoeff_model0_behavior1 == 0
    %     tempName = 'FittingResults_GLMlikeLocProdModel_20230705A.mat'; %#ok<*UNRCH>
    tempName = 'FittingResults_GLMlikeLocProdModel_20230706A.mat'; %#ok<*UNRCH>
    S = load(['C:\ASDROOT\STUDY\temp','\data\',tempName]);
    FittingResults = S.FittingResults;
    ModelParams = FittingResults.ModelParams;
    coeff_p = ModelParams.coeff_p;
    coeff_p0 = ModelParams.coeff_p0;
elseif if_loadCoeff_model0_behavior1 == 1
    tempName = 'errorTransitionMatrix_eachLoc_20230808A.mat'; %#ok<*UNRCH>
    S = load(['C:\ASDROOT\STUDY\temp','\data\',tempName]);
    errorTransitionMatrix_eachLoc = S.errorTransitionMatrix_eachLoc;
    coeff_p = 1-errorTransitionMatrix_eachLoc.deleteOne;
    coeff_p0 = errorTransitionMatrix_eachLoc.addOne;
end

dis_B_x1;
dis_B_x2;

a = 1;

dis_B_x1_norm = norm(dis_B_x1,2);
dis_B_x2_norm = norm(dis_B_x2,2);



dis_B_x1_raw = dis_B_x1;
% dis_B_x1 = (dis_B_x1_raw - dis_B_x2);
% dis_B_x1 = dis_B_x2 + dis_B_x1_raw;
% dis_B_x1 = dis_B_x1_raw ./ dis_B_x2;
% dis_B_x1 = dis_B_x2 ./ dis_B_x1_raw;
% dis_B_x1 = dis_B_x1_raw .* dis_B_x2;
% dis_B_x1 = dis_B_x1_raw./norm(dis_B_x1_raw,2) - dis_B_x2./norm(dis_B_x2,2);
% dis_B_x1 = dis_B_x1_raw - dis_B_x2;



dis_B_x2_raw = dis_B_x2;
% dis_B_x2 = dis_B_x1 - dis_B_x2_raw;
% dis_B_x2 = dis_B_x2_raw + dis_B_x1;
% dis_B_x2 = dis_B_x1 ./ dis_B_x2_raw;
% dis_B_x2 = dis_B_x2_raw ./ dis_B_x1;
% dis_B_x2 = dis_B_x1 .* dis_B_x2_raw;

x = dis_B_x1';
y = coeff_p';
[r_P,p_P] = corr(x,y);
% cftool(x,y);

x = dis_B_x2';
y = coeff_p0';
[r_P0,p_P0] = corr(x,y);
% cftool(x,y);

fprintf('r_P=%.4f, p_P=%.4f.\n',r_P,p_P);
fprintf('r_P0=%.4f, p_P0=%.4f.\n',r_P0,p_P0);

fig2 = figure('Name','Fig2','NumberTitle','off');
set(gcf,'Position',[35+300 35+200 900 450]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>

multi_rgbColor = ...
    [228,26,28;
    55,126,184;
    77,175,74;
    152,78,163;
    255,127,0;
    255,255,51]/255;

backgounrdColor = [1 1 1]*0.825;%0.875


nexttile
for tempi=1:numFrames
    scatter(dis_B_x1(tempi),coeff_p(tempi),48,multi_rgbColor(tempi,:),'filled', ...
        'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'MarkerEdgeColor',[0 0 0]);%36
    hold on
end
grid on
set(gca,'color',backgounrdColor);
xlabel(sprintf('Miss distance'),'fontsize',16)
ylabel(sprintf('P (Behavior)'),'fontsize',16)


nexttile

[temp_x_min,temp_x_max] = bounds(dis_B_x2);

for tempi=1:numFrames
    temp_dis_delay2_B_P0 = dis_B_x2(tempi);
    if tempi > 1
        if abs(temp_dis_delay2_B_P0-dis_B_x2(tempi-1)) < (temp_x_max-temp_x_min) * 0.02
            temp_dis_delay2_B_P0 = temp_dis_delay2_B_P0 + (temp_x_max-temp_x_min) * 0.02;
        end
    end
    scatter(temp_dis_delay2_B_P0,coeff_p0(tempi),48,multi_rgbColor(tempi,:),'filled', ...
        'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'MarkerEdgeColor',[0 0 0]);%36
    hold on
end
grid on
set(gca,'color',backgounrdColor);
xlabel(sprintf('False alarm distance'),'fontsize',16)
ylabel(sprintf('P0 (Behavior)'),'fontsize',16)


a = 1;







